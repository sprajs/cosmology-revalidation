"""New broad, exact-CAMB numerical training support; never posterior draws."""
import os
os.environ['OMP_NUM_THREADS']='1'
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
from datetime import datetime,timezone
import hashlib
import json
import multiprocessing
from pathlib import Path
import sys
import time
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
EXTERNAL=HERE.parent/'external_probes'
sys.path.insert(0,str(EXTERNAL))
import camb
from adapter import reference_point
from modern_adapter import modern_info
from fast_lensing import use_fast_lensing

DIRECTORY=ROOT/'.work/unified-cosmology/inference/broad-spectral'
RESULT=ROOT/'studies/unified_cosmology/results/inference'
DESIGN=HERE/'broad-spectral-design.json'
COORDS=['100thetaMC','ombh2','omch2','logA','ns','tau','w','wa']
PHYSICAL=['H0','ombh2','omch2','logA','ns','tau','w','wa']
SPECTRA=['tt','ee','bb','te','pp']
_MODEL=None


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_identities():
    paths=[Path(__file__),EXTERNAL/'modern_adapter.py',EXTERNAL/'adapter.py',EXTERNAL/'fast_lensing.py']
    return {str(p.relative_to(ROOT)):sha(p) for p in paths}


def initialize():
    global _MODEL
    from cobaya.model import get_model
    info=use_fast_lensing(modern_info())
    _MODEL=get_model(info)
    _MODEL.add_requirements({'CAMBdata':None})


def exact(job):
    split,index,coords=job
    file=DIRECTORY/split/f'{index:04d}.npz';recfile=file.with_suffix('.json')
    identities=source_identities()
    if recfile.exists():
        rec=json.loads(recfile.read_text())
        assert rec['requested_coordinates']==coords.tolist()
        assert rec['source_sha256']==identities and rec['design_sha256']==sha(DESIGN)
        if rec['status']=='finite_exact':assert sha(ROOT/rec['file'])==rec['sha256']
        return rec
    started=time.monotonic()
    rec=dict(split=split,index=index,requested_coordinates=coords.tolist(),
        source_sha256=identities,design_sha256=sha(DESIGN))
    try:
        vals=dict(zip(COORDS,coords));extra=modern_info()['theory']['camb']['extra_args']
        cosm={k:vals[k] for k in ['ombh2','omch2','ns','tau','w','wa']}
        pars=camb.set_params(cosmomc_theta=vals['100thetaMC']/100.,
            As=1e-10*np.exp(vals['logA']),**cosm,**extra)
        point=reference_point(_MODEL);point.update(cosm,H0=pars.H0,logA=vals['logA'])
        result=_MODEL.logposterior(point)
        if not np.isfinite(result.logpost):
            rec.update(status='nonfinite_prior_theory_or_likelihood',
                physical_point={k:float(point[k]) for k in PHYSICAL},
                loglikes=[float(v) if np.isfinite(v) else None for v in result.loglikes])
        else:
            dls=_MODEL.provider.get_Cl(ell_factor=True)
            full=_MODEL.provider.get_CAMBdata();p=full.Params
            actual=np.array([100*full.cosmomc_theta(),p.ombh2,p.omch2,np.log(1e10*p.InitPower.As),
                p.InitPower.ns,p.Reion.optical_depth,p.DarkEnergy.w,p.DarkEnergy.wa])
            arrays=np.array([dls[k] for k in SPECTRA]);assert np.isfinite(arrays).all()
            np.savez_compressed(file,coordinates=actual,requested_coordinates=coords,
                physical_point=np.array([point[k] for k in PHYSICAL]),spectra=arrays,
                ell=dls['ell'],loglikes=result.loglikes,derived_rdrag=np.array(full.get_derived_params()['rdrag']))
            rec.update(status='finite_exact',file=str(file.relative_to(ROOT)),sha256=sha(file),
                physical_point={k:float(point[k]) for k in PHYSICAL},actual_coordinates=actual.tolist(),
                likelihood_names=list(_MODEL.likelihood),loglikes=list(map(float,result.loglikes)),
                provider_Dl_length=len(dls['ell']),CAMB_Params_max_l=int(p.max_l),
                CAMB_Params_max_eta_k=float(p.max_eta_k),
                finalized_theory_extra_args=dict(_MODEL.theory['camb'].extra_args))
    except Exception as error:
        rec.update(status='error',exception=type(error).__name__,message=str(error))
    rec['seconds']=time.monotonic()-started
    temp=recfile.with_suffix('.json.tmp');temp.write_text(json.dumps(rec,indent=2,allow_nan=False)+'\n');temp.replace(recfile)
    return rec


def make_design():
    original=EXTERNAL/'spectral-training-design.json';old=json.loads(original.read_text())
    centre=np.array(old['centre']);scale=np.diag([1.,1.,1.,1.,1.,1.,3.,3.])
    chol=np.array(old['coordinate_cholesky'])@scale
    design=dict(frozen_utc=datetime.now(timezone.utc).isoformat(),
        timing='New numerical support declared after slow luminosity-shift pilot fallbacks, before new exact spectra or holdout errors.',
        coordinate_names=COORDS,physical_names=PHYSICAL,spectra_order=SPECTRA,
        centre=centre.tolist(),coordinate_cholesky=chol.tolist(),coordinate_covariance=(chol@chol.T).tolist(),
        transformation='original_coordinate_cholesky @ diag(1,1,1,1,1,1,3,3); original centre unchanged',
        original_design_sha256=sha(original),counts={'train':768,'holdout':128},
        seeds={'train':2727701,'holdout':2727702},
        draw_rule='Independent whitened standard Gaussian vectors; independently90percent width1.5 and10percent width2.5. Every requested point retained; no rejection/redraw of invalid priors or theory domains.',
        fit={'degree':3,'features':165,'minimum_finite_train':330,
             'method':'Unregularized least squares by full-rank SVD, separate RMS output scaling. All finite training rows, no holdout fitting or hyperparameter selection.',
             'feature_singular_value_ratio_min':1e-10,'model_file':str((DIRECTORY/'surrogate-cubic-v1.npz').relative_to(ROOT))},
        theory='Unchanged modern_adapter native CAMB settings and scientific priors; exact validated fast-lensing matrix reassociation only.',
        scientific_scope='A wider computational proposal, not a wider or narrower scientific prior. Heldout errors at fixed nuisance values and exact posterior correction are required; no posterior inference from training counts.',
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in [EXTERNAL/'modern_adapter.py',EXTERNAL/'adapter.py',EXTERNAL/'fast_lensing.py',HERE/'spectral_surrogate.py']},
        modern_acquisition_sha256=sha(ROOT/'studies/unified_cosmology/results/external_probes/modern-acquisition.json'))
    if DESIGN.exists():
        frozen=json.loads(DESIGN.read_text());design['frozen_utc']=frozen['frozen_utc'];assert frozen==design
    else:DESIGN.write_text(json.dumps(design,indent=2)+'\n')
    DIRECTORY.mkdir(parents=True,exist_ok=True);jobs=[]
    for split,count in design['counts'].items():
        (DIRECTORY/split).mkdir(exist_ok=True)
        rng=np.random.default_rng(design['seeds'][split]);draw=rng.normal(size=(count,8))
        widths=np.where(rng.random(count)<.1,2.5,1.5)
        coordinates=centre+(draw*widths[:,None])@chol.T
        jobs.extend((split,j,coords) for j,coords in enumerate(coordinates))
    arrays=DIRECTORY/'design-arrays.npz'
    if not arrays.exists():np.savez(arrays,centre=centre,coordinate_cholesky=chol,
        train=np.array([j[2] for j in jobs if j[0]=='train']),holdout=np.array([j[2] for j in jobs if j[0]=='holdout']))
    else:
        cache=np.load(arrays)
        for split in design['counts']:assert np.array_equal(cache[split],np.array([j[2] for j in jobs if j[0]==split]))
    return design,jobs


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--freeze',action='store_true')
    parser.add_argument('--workers',type=int,default=8);parser.add_argument('--limit',type=int)
    args=parser.parse_args();design,jobs=make_design()
    if args.freeze:
        print(json.dumps({'design':str(DESIGN),'sha256':sha(DESIGN),'planned':len(jobs),'sources':source_identities()},indent=2));return
    started=time.monotonic();initial=source_identities()
    if args.limit:jobs=jobs[:args.limit]
    rows=[]
    with ProcessPoolExecutor(max_workers=args.workers,initializer=initialize,
                            mp_context=multiprocessing.get_context('spawn')) as pool:
        futures=[pool.submit(exact,job) for job in jobs]
        for future in as_completed(futures):
            rows.append(future.result())
            if len(rows)%16==0 or len(rows)==len(jobs):
                print(json.dumps(dict(completed=len(rows),planned=len(jobs),seconds=time.monotonic()-started,
                    finite=sum(r['status']=='finite_exact' for r in rows))),flush=True)
    assert initial==source_identities()
    counts={split:{status:sum(r['split']==split and r['status']==status for r in rows)
                  for status in sorted({r['status'] for r in rows})} for split in ['train','holdout']}
    first=next((r for r in rows if r['status']=='finite_exact'),{})
    result=dict(status='complete_exact_acquisition' if len(rows)==896 else 'partial_acquisition',
        counts=counts,seconds=time.monotonic()-started,workers=args.workers,OMP_threads_per_worker=1,
        manifest_directory=str(DIRECTORY.relative_to(ROOT)),
        theory_metadata={k:first.get(k) for k in ['provider_Dl_length','CAMB_Params_max_l','CAMB_Params_max_eta_k','finalized_theory_extra_args','likelihood_names']},
        source_sha256=initial,design_sha256=sha(DESIGN),design_arrays_sha256=sha(DIRECTORY/'design-arrays.npz'),
        rows=[dict(split=r['split'],index=r['index'],status=r['status'],manifest_sha256=sha(DIRECTORY/r['split']/f"{r['index']:04d}.json"))
              for r in sorted(rows,key=lambda v:(v['split'],v['index']))],
        scope='Native exact spectra/likelihood acquisition, not emulator accuracy or posterior qualification.')
    (RESULT/'broad-spectral-training.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
