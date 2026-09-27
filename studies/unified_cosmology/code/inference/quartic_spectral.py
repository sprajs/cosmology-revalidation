"""Declared quartic repair with merged training and a new untouched holdout.

The earlier broad holdout informed this repair; it is not reused as validation.
Scientific priors and native theory stay unchanged. This does not launch chains.
"""
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
import time
import numpy as np
import broad_spectral_training as native
from spectral_surrogate import polynomial,polynomial_exponents,replace_theory

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
WORK=ROOT/'.work/unified-cosmology/inference/quartic-spectral'
RESULT=ROOT/'studies/unified_cosmology/results/inference'
DESIGN=HERE/'quartic-spectral-design.json'
ORIGINAL=ROOT/'.work/unified-cosmology/external-probes/spectral-training'
BROAD=ROOT/'.work/unified-cosmology/inference/broad-spectral'
NATIVE_IDENTITIES=native.source_identities


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def identities():
    return dict(NATIVE_IDENTITIES(),**{str(Path(__file__).relative_to(ROOT)):sha(__file__)})


def initialize():
    native.DIRECTORY=WORK;native.DESIGN=DESIGN;native.source_identities=identities
    native.initialize()


def evaluate(job):return native.exact(job)


def design_and_jobs():
    old_path=HERE.parent/'external_probes/spectral-training-design.json'
    broad_path=HERE/'broad-spectral-design.json'
    old=json.loads(old_path.read_text());broad=json.loads(broad_path.read_text())
    centre=np.array(old['centre']);chol=np.array(broad['coordinate_cholesky'])
    plan=dict(frozen_utc=datetime.now(timezone.utc).isoformat(),
        timing='Declared after the broad cubic holdout exposed large errors near the native high-redshift domain boundary; earlier holdouts are observed diagnostics only.',
        question='Can a higher-order numerical proposal retain wider support and recover likelihood accuracy without changing the scientific target?',
        centre=centre.tolist(),coordinate_cholesky=chol.tolist(),coordinate_names=native.COORDS,
        training={'original_finite':509,'broad_finite':689,'total_expected':1198,
                  'sources':['original/train','broad/train'],'holdout_rows_used':0},
        fit={'degree':4,'features':495,'minimum_finite_train':990,
             'method':'One unregularized full-rank SVD fit, with per-output RMS scaling, all1198 finite training rows. No data-dependent weighting or scientific-prior change.',
             'feature_singular_value_ratio_min':1e-10,'model_file':str((WORK/'surrogate-quartic-v1.npz').relative_to(ROOT))},
        counts={'holdout_original':64,'holdout_broad':64},seeds={'holdout_original':2727801,'holdout_broad':2727802},
        draw_rule='Independent whitened Gaussian draws,90percent sd1.5 and10percent sd2.5, with original or broad Cholesky according to split. No redraw or replacement of invalid prior/theory points.',
        holdout_cholesky={'holdout_original':old['coordinate_cholesky'],'holdout_broad':broad['coordinate_cholesky']},
        frozen_inputs_sha256={str(p.relative_to(ROOT)):sha(p) for p in [old_path,broad_path,
            RESULT/'broad-spectral-training.json',RESULT/'broad-spectral-fit.json',RESULT/'broad-spectral-holdout.json',
            HERE/'spectral_surrogate.py',HERE.parent/'external_probes/modern_adapter.py',HERE.parent/'external_probes/fast_lensing.py']},
        scope='Numerical proposal only. The original and broad cubic models, all attempted spectra and exposed holdouts remain intact. Native scientific settings/domain remain unchanged; exact posterior correction and convergence are mandatory. No chains launched by this module.')
    if DESIGN.exists():
        saved=json.loads(DESIGN.read_text());plan['frozen_utc']=saved['frozen_utc'];assert saved==plan
    else:DESIGN.write_text(json.dumps(plan,indent=2)+'\n')
    WORK.mkdir(exist_ok=True,parents=True);jobs=[]
    for split,count in plan['counts'].items():
        (WORK/split).mkdir(exist_ok=True)
        rng=np.random.default_rng(plan['seeds'][split]);normal=rng.normal(size=(count,8))
        widths=np.where(rng.random(count)<.1,2.5,1.5)
        requested=centre+(normal*widths[:,None])@np.array(plan['holdout_cholesky'][split]).T
        jobs.extend((split,j,point) for j,point in enumerate(requested))
    return plan,jobs


def acquire(workers):
    plan,jobs=design_and_jobs();started=time.monotonic();rows=[];source=identities()
    with ProcessPoolExecutor(max_workers=workers,initializer=initialize,mp_context=multiprocessing.get_context('spawn')) as pool:
        for future in as_completed([pool.submit(evaluate,job) for job in jobs]):
            rows.append(future.result())
            if len(rows)%16==0:print(json.dumps(dict(completed=len(rows),planned=128,seconds=time.monotonic()-started)),flush=True)
    assert identities()==source
    statuses=sorted({r['status'] for r in rows})
    result=dict(status='complete_exact_acquisition',counts={split:{s:sum(r['split']==split and r['status']==s for r in rows) for s in statuses} for split in plan['counts']},
        seconds=time.monotonic()-started,workers=workers,OMP_threads_per_worker=1,
        source_sha256=source,design_sha256=sha(DESIGN),
        rows=[dict(split=r['split'],index=r['index'],status=r['status'],manifest_sha256=sha(WORK/r['split']/f"{r['index']:04d}.json")) for r in sorted(rows,key=lambda r:(r['split'],r['index']))],
        scope='New native holdout acquisition only; no accuracy values used to fit or select the quartic model.')
    assert len(rows)==128
    (RESULT/'quartic-spectral-acquisition.json').write_text(json.dumps(result,indent=2)+'\n')


def fit():
    plan,_=design_and_jobs();dependencies={str(DESIGN.relative_to(ROOT)):sha(DESIGN)}
    points,spectra,counts=[],[],{}
    for label,folder,expected in [('original',ORIGINAL,509),('broad',BROAD,689)]:
        count=0
        for manifest in sorted((folder/'train').glob('*.json')):
            record=json.loads(manifest.read_text());dependencies[str(manifest.relative_to(ROOT))]=sha(manifest)
            assert record['split']=='train'
            if record['status']!='finite_exact':continue
            path=ROOT/record['file'];assert sha(path)==record['sha256']
            dependencies[str(path.relative_to(ROOT))]=sha(path)
            data=np.load(path,allow_pickle=False)
            assert record['CAMB_Params_max_l']==10251 and record['provider_Dl_length']==10152
            points.append(data['coordinates']);spectra.append(data['spectra']);count+=1
        assert count==expected,(label,count,expected)
        counts[label]=count
    assert len(points)==1198 and len(points)>=plan['fit']['minimum_finite_train']
    centre=np.array(plan['centre']);chol=np.array(plan['coordinate_cholesky'])
    white=np.linalg.solve(chol,(np.array(points)-centre).T).T
    exponents=polynomial_exponents(degree=4);assert len(exponents)==495
    features=polynomial(white,exponents);spectra=np.array(spectra)
    length=spectra.shape[-1];target=spectra.reshape(len(points),-1)
    scale=np.maximum(np.sqrt(np.mean(target**2,axis=0)),1e-30)
    u,s,vh=np.linalg.svd(features,full_matrices=False)
    assert s[-1]/s[0]>plan['fit']['feature_singular_value_ratio_min']
    coefficients=((vh.T/s)@u.T)@(target/scale)
    residual=features@coefficients-target/scale
    output=ROOT/plan['fit']['model_file'];assert not output.exists()
    np.savez_compressed(output,centre=centre,coordinate_cholesky=chol,exponents=exponents,
        coefficients=coefficients,output_scale=scale,length=length)
    dependencies.update({str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'spectral_surrogate.py']})
    record=dict(status='fitted_not_posterior_qualified',created_utc=datetime.now(timezone.utc).isoformat(),degree=4,
        training_counts=counts,training_rows=len(points),holdout_rows_used=0,coefficients_per_output=495,
        feature_condition_number=float(s[0]/s[-1]),training_scaled_rms=float(np.sqrt(np.mean(residual**2))),
        training_scaled_max=float(np.max(abs(residual))),model_file=str(output.relative_to(ROOT)),model_sha256=sha(output),
        source_sha256=dependencies,scope='Declared quartic numerical repair; never a new scientific likelihood. Requires untouched mixed-width holdout checks and native exact posterior correction.')
    output.with_suffix('.json').write_text(json.dumps(record,indent=2)+'\n')
    (RESULT/'quartic-spectral-fit.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({k:v for k,v in record.items() if k!='source_sha256'},indent=2))


def check():
    from cobaya.model import get_model
    plan,_=design_and_jobs();fit_record=json.loads((RESULT/'quartic-spectral-fit.json').read_text())
    model_file=ROOT/plan['fit']['model_file'];assert sha(model_file)==fit_record['model_sha256']
    acquisition_path=RESULT/'quartic-spectral-acquisition.json';acquisition=json.loads(acquisition_path.read_text())
    assert acquisition['status']=='complete_exact_acquisition' and acquisition['design_sha256']==sha(DESIGN)
    dependencies={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),DESIGN,acquisition_path,
        RESULT/'quartic-spectral-fit.json',model_file,HERE/'spectral_surrogate.py']}
    rows=[]
    with get_model(native.use_fast_lensing(replace_theory(native.modern_info(),model_file))) as model:
        ref=native.reference_point(model);theory=model.theory['spectral_surrogate']
        for record in acquisition['rows']:
            split,index=record['split'],record['index'];manifest=WORK/split/f'{index:04d}.json'
            assert sha(manifest)==record['manifest_sha256'];meta=json.loads(manifest.read_text())
            dependencies[str(manifest.relative_to(ROOT))]=sha(manifest)
            if meta['status']!='finite_exact':
                rows.append(dict(split=split,index=index,status='native_'+meta['status']));continue
            path=ROOT/meta['file'];assert sha(path)==meta['sha256'];dependencies[str(path.relative_to(ROOT))]=sha(path)
            before=theory.exact_calls;started=time.monotonic()
            result=model.logposterior(dict(ref,**meta['physical_point']))
            if not np.isfinite(result.logpost):
                rows.append(dict(split=split,index=index,status='nonfinite_proposal'));continue
            native_like=dict(zip(meta['likelihood_names'],meta['loglikes']))
            candidate=dict(zip(model.likelihood,result.loglikes));assert set(native_like)==set(candidate)
            delta={k:float(candidate[k]-v) for k,v in native_like.items()}
            rows.append(dict(split=split,index=index,status='compared',delta_loglike_by_component=delta,
                delta_loglike_total=sum(delta.values()),native_fixed_reference_CMB_BAO_loglike=sum(meta['loglikes']),
                w=meta['physical_point']['w'],wa=meta['physical_point']['wa'],
                exact_fallback=theory.exact_calls>before,seconds=time.monotonic()-started))
    summaries={}
    for split in list(plan['counts'])+['combined']:
        group=[r for r in rows if r['status']=='compared' and (split=='combined' or r['split']==split)]
        d=np.array([r['delta_loglike_total'] for r in group])
        summaries[split]=dict(finite_comparisons=len(group),delta_loglike_rms=float(np.sqrt(np.mean(d*d))) if len(d) else None,
            delta_loglike_quantiles=np.quantile(d,[0,.025,.16,.5,.84,.975,1]).tolist() if len(d) else None,
            exact_fallbacks=sum(r['exact_fallback'] for r in group))
    result=dict(status='checked_requires_exact_posterior_correction' if not any(r['status']=='nonfinite_proposal' for r in rows) else 'nonfinite_proposal_failures',
        created_utc=datetime.now(timezone.utc).isoformat(),summaries=summaries,rows=rows,
        model_sha256=sha(model_file),dependencies_sha256=dependencies,
        scope='First untouched validation of the fixed quartic repair, split into original and broader numerical coordinate distributions. No posterior qualification follows without exact correction and convergence.')
    destination=RESULT/'quartic-spectral-holdout.json';assert not destination.exists()
    destination.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(summaries,indent=2))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['freeze','acquire','fit','check'])
    parser.add_argument('--workers',type=int,default=6);args=parser.parse_args()
    if args.action=='freeze':
        design_and_jobs();print(json.dumps(dict(design_sha256=sha(DESIGN),source_sha256=identities()),indent=2))
    elif args.action=='acquire':acquire(args.workers)
    elif args.action=='fit':fit()
    else:check()


if __name__=='__main__':main()
