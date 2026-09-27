"""Bounded exact-spectrum training/held-out acquisition, never posterior draws.

Four independent CAMB workers, one OpenMP/BLAS thread each. Native likelihood
and theory are evaluated without approximation. All points/rejections persist.
"""
import os
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
from datetime import datetime,timezone
import json
import multiprocessing
from pathlib import Path
import time
import numpy as np
import camb
from scipy.linalg import solve_triangular
from adapter import WORK,ROOT,reference_point
from modern_adapter import modern_info
from acquire import HERE,RESULTS,sha

COORDS=['100thetaMC','ombh2','omch2','logA','ns','tau','w','wa']
PHYSICAL=['H0','ombh2','omch2','logA','ns','tau','w','wa']
SPECTRA=['tt','ee','bb','te','pp']
DIRECTORY=WORK/'spectral-training'
_MODEL=None


def initialize():
    global _MODEL
    from cobaya.model import get_model
    _MODEL=get_model(modern_info())
    _MODEL.add_requirements({'CAMBdata':None})


def exact(job):
    split,index,coords=job
    file=DIRECTORY/split/f'{index:04d}.npz';recfile=file.with_suffix('.json')
    if recfile.exists():
        rec=json.loads(recfile.read_text())
        assert rec['requested_coordinates']==coords.tolist()
        assert rec['training_code_sha256']==sha(__file__)
        assert rec['adapter_sha256']==sha(HERE/'modern_adapter.py')
        return rec
    start=time.monotonic()
    rec={'split':split,'index':index,'requested_coordinates':coords.tolist(),
         'training_code_sha256':sha(__file__),'adapter_sha256':sha(HERE/'modern_adapter.py')}
    try:
        # CAMB's acoustic-angle root is solved before spectral computation.
        vals=dict(zip(COORDS,coords));args=modern_info()['theory']['camb']['extra_args']
        cosm={k:vals[k] for k in ['ombh2','omch2','ns','tau','w','wa']}
        pars=camb.set_params(cosmomc_theta=vals['100thetaMC']/100,
                             As=1e-10*np.exp(vals['logA']),**cosm,**args)
        point=reference_point(_MODEL);point.update(cosm,H0=pars.H0,logA=vals['logA'])
        result=_MODEL.logposterior(point)
        if not np.isfinite(result.logpost):
            rec.update(status='nonfinite_prior_theory_or_likelihood',physical_point={k:float(point[k]) for k in PHYSICAL},
                       loglikes=[float(v) if np.isfinite(v) else None for v in result.loglikes])
        else:
            dls=_MODEL.provider.get_Cl(ell_factor=True)
            full=_MODEL.provider.get_CAMBdata();p=full.Params
            actual=np.array([100*full.cosmomc_theta(),p.ombh2,p.omch2,np.log(1e10*p.InitPower.As),
                             p.InitPower.ns,p.Reion.optical_depth,p.DarkEnergy.w,p.DarkEnergy.wa])
            arrays=np.array([dls[k] for k in SPECTRA])
            assert np.isfinite(arrays).all()
            np.savez_compressed(file,coordinates=actual,requested_coordinates=coords,
                physical_point=np.array([point[k] for k in PHYSICAL]),
                spectra=arrays,ell=dls['ell'],loglikes=result.loglikes,
                derived_rdrag=np.array(full.get_derived_params()['rdrag']))
            rec.update(status='finite_exact',file=str(file.relative_to(ROOT)),sha256=sha(file),
                physical_point={k:float(point[k]) for k in PHYSICAL},actual_coordinates=actual.tolist(),
                likelihood_names=list(_MODEL.likelihood),loglikes=list(map(float,result.loglikes)),
                provider_Dl_length=len(dls['ell']),CAMB_Params_max_l=int(p.max_l),
                CAMB_Params_max_eta_k=float(p.max_eta_k),
                finalized_theory_extra_args=dict(_MODEL.theory['camb'].extra_args))
    except Exception as e:
        rec.update(status='error',exception=type(e).__name__,message=str(e))
    rec['seconds']=time.monotonic()-start
    recfile.write_text(json.dumps(rec,indent=2)+'\n')
    return rec


def make_design():
    p=WORK/'proposal-cpl-dovekie.covmat'
    names=p.read_text().splitlines()[0].lstrip('# ').split()
    c=np.loadtxt(p);ind=[names.index(k) for k in PHYSICAL];c=c[np.ix_(ind,ind)]
    source=json.loads((RESULTS/'proposal.json').read_text())
    J=np.eye(8);J[0,:]=0
    for key,value in source['derivatives'].items():J[0,PHYSICAL.index(key)]=value
    cov=J@c@J.T;chol=np.linalg.cholesky(cov)
    centre=np.array([source['reference_100thetaMC'],.02237,.12,3.044,.9649,.0544,-.85,-.55])
    design={'frozen_utc':datetime.now(timezone.utc).isoformat(),'coordinate_names':COORDS,
        'physical_names':PHYSICAL,'spectra_order':SPECTRA,'centre':centre.tolist(),
        'coordinate_covariance':cov.tolist(),'coordinate_cholesky':chol.tolist(),
        'counts':{'train':512,'holdout':96},'seeds':{'train':2727601,'holdout':2727602},
        'draw_rule':'Independent zero-mean Gaussian whitened draws:90percent sd1.5 and10percent sd2.5; no posterior label, no rejection/redraw to hide invalid points.',
        'proposal_sha256':sha(p),'proposal_transform_sha256':sha(RESULTS/'proposal.json'),
        'adapter_sha256':sha(HERE/'modern_adapter.py'),'modern_acquisition_sha256':sha(RESULTS/'modern-acquisition.json'),
        'scope':'Training support is a local numerical acceleration design, not a scientific prior. Exact CAMB posterior correction and held-out validation are required before using an emulator scientifically.'}
    dest=HERE/'spectral-training-design.json'
    if dest.exists():
        old=json.loads(dest.read_text());design['frozen_utc']=old['frozen_utc'];assert old==design
    else:dest.write_text(json.dumps(design,indent=2)+'\n')
    jobs=[]
    for split,count in design['counts'].items():
        (DIRECTORY/split).mkdir(parents=True,exist_ok=True)
        rng=np.random.default_rng(design['seeds'][split]);x=rng.normal(size=(count,8))
        widths=np.where(rng.random(count)<.1,2.5,1.5)
        coords=centre+(x*widths[:,None])@chol.T
        jobs += [(split,i,v) for i,v in enumerate(coords)]
    np.savez(DIRECTORY/'design-arrays.npz',centre=centre,coordinate_covariance=cov,coordinate_cholesky=chol,
             train=np.array([j[2] for j in jobs if j[0]=='train']),holdout=np.array([j[2] for j in jobs if j[0]=='holdout']))
    return design,jobs


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--limit',type=int);parser.add_argument('--workers',type=int,default=4)
    args=parser.parse_args();DIRECTORY.mkdir(parents=True,exist_ok=True)
    design,jobs=make_design();start=time.monotonic()
    if args.limit:jobs=jobs[:args.limit]
    rows=[]
    with ProcessPoolExecutor(max_workers=args.workers,initializer=initialize,
                             mp_context=multiprocessing.get_context('spawn')) as pool:
        futures=[pool.submit(exact,j) for j in jobs]
        for f in as_completed(futures):
            r=f.result();rows.append(r)
            if len(rows)%16==0 or len(rows)==len(jobs):
                print(json.dumps({'completed':len(rows),'planned':len(jobs),'seconds':time.monotonic()-start,
                                  'finite':sum(r['status']=='finite_exact' for r in rows)}),flush=True)
    counts={s:{status:sum(r['split']==s and r['status']==status for r in rows)
               for status in sorted(set(r['status'] for r in rows))} for s in ['train','holdout']}
    first=next((r for r in rows if r['status']=='finite_exact'),{})
    out={'status':'complete_exact_acquisition' if len(rows)==608 else 'partial_acquisition',
         'counts':counts,'seconds':time.monotonic()-start,'workers':args.workers,'OMP_threads_per_worker':os.environ.get('OMP_NUM_THREADS'),
         'manifest_directory':str(DIRECTORY.relative_to(ROOT)),
         'theory_metadata':{k:first.get(k) for k in ['provider_Dl_length','CAMB_Params_max_l','CAMB_Params_max_eta_k','finalized_theory_extra_args','likelihood_names']},
         'design_sha256':sha(HERE/'spectral-training-design.json'),'code_sha256':sha(__file__),
         'rows':[{'split':r['split'],'index':r['index'],'status':r['status'],
                  'manifest_sha256':sha(DIRECTORY/r['split']/f"{r['index']:04d}.json")} for r in sorted(rows,key=lambda q:(q['split'],q['index']))],
         'scope':'Exact likelihood/spectrum acquisition only, not emulator validation or posterior inference.'}
    (RESULTS/'spectral-training.json').write_text(json.dumps(out,indent=2)+'\n')

if __name__=='__main__':main()
