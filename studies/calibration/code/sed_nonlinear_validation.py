"""Full frozen constructed nonlinear gate, max two worker processes."""
from pathlib import Path
import argparse
import concurrent.futures
import datetime
import hashlib
import importlib.util
import json
import multiprocessing
import os
import platform
import sys
import time

for variable in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[variable]='1'
import numpy as np
import pandas as pd
from scipy.linalg import solve_triangular

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'runs/research_2026_09_26/sed_nonlinear_validation'
BASE=ROOT/'runs/research_2026_09_26/astra_design/validation1020'
PROTOCOL=ROOT/'docs/research-2026-09-26/sed-nonlinear-validation.md'
spec=importlib.util.spec_from_file_location('gate',ROOT/'scripts/research_2026_09_26/sed_nonlinear_refit.py')
gate=importlib.util.module_from_spec(spec);spec.loader.exec_module(gate)
WORKER={}


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def atomic_json(path,value):
    temp=Path(str(path)+'.tmp')
    temp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    temp.replace(path)


def initialize_worker():
    model,bands,paths,zp=gate.sed.fr.build_model()
    with np.load(OUT/'frozen-coefficients.npz') as d:
        coefficients=d['sed_coefficients'];observer=d['observer_griz']
    WORKER.update(model=model,bands=bands,offsets={str(x['Filter Name'])[-1]:float(x['Primary Mag']) for x in zp},coefficients=coefficients,observer=observer)


def support_check(engine,theta,mode):
    phase=(engine.t-engine.t0-theta[3])/(1+engine.z)
    source=engine.model.source
    phase_support=bool(phase.min()>=source.minphase() and phase.max()<=source.maxphase())
    wavelength_support=True;negative_share=0.;negative_epochs=0;maxdm=0.
    # Reassert parameters at this solution before inspecting spectral support.
    flux=engine.flux(theta,mode)
    for idx,wave,weights,shapes,observer in engine.grids:
        rest=wave/(1+engine.z)
        wavelength_support &= bool(rest.min()>=source.minwave() and rest.max()<=source.maxwave())
        spectral=engine.model.flux(engine.t[idx],wave)
        dm=np.polynomial.legendre.legvander(np.tanh(phase[idx]/20),2)@shapes
        maxdm=max(maxdm,float(np.max(abs(dm[:,weights>0]))))
        if mode=='sed':spectral*=np.exp(-gate.K*dm)
        elif mode=='observer':spectral*=np.exp(-gate.K*observer)
        fw=spectral*weights
        neg=np.sum(abs(np.minimum(fw,0)),axis=1)/np.maximum(np.sum(abs(fw),axis=1),1e-300)
        negative_share=max(negative_share,float(neg.max()))
        negative_epochs+=int(np.sum(neg>1e-6))
    return dict(gate=bool(phase_support and wavelength_support and np.all(np.isfinite(flux))),
                phase_support=phase_support,wavelength_support=wavelength_support,
                phase_min=float(phase.min()),phase_max=float(phase.max()),
                source_phase_domain=[float(source.minphase()),float(source.maxphase())],
                sed_max_abs_mag_at_fitted_phase=maxdm,
                negative_absolute_spectral_photon_share_max=negative_share,
                epochs_with_negative_spectral_share_gt_1e_minus6=negative_epochs,
                nonpositive_integrated_model_epochs=int(np.sum(flux<=0)))


def one_object(task):
    cid,official_sensitivity=task
    started=time.perf_counter()
    objpath=BASE/'objectives'/f'objective_{cid}.npz'
    cachepath=BASE/'analysis/objects'/f'{cid}.npz'
    try:
        with np.load(objpath) as data:
            keys=['MJD','band','model_flux','data_flux','data_fluxerr','zHEL','MWEBV','parameters_x0_x1_c_t0','frozen_flux_covariance']
            nominal={k:data[k] for k in keys}
        e=gate.Engine(nominal,WORKER['model'],WORKER['bands'],WORKER['offsets'],WORKER['coefficients'],WORKER['observer'])
        native=e.flux(np.zeros(4));cov=nominal['frozen_flux_covariance'];official=nominal['model_flux']
        chol=np.linalg.cholesky(cov);white=lambda v:solve_triangular(chol,v,lower=True)
        full=np.arange(len(native));order=pd.DataFrame({'t':nominal['MJD'],'b':nominal['band']}).sort_values(['t','b'],kind='stable').index.to_numpy()
        with np.load(cachepath) as d:cached=d['native_flux_model']
        mean_error=float(np.max(abs(native[order]-cached)/np.maximum(abs(cached),1e-30)))
        assert mean_error<1e-9,(cid,mean_error)
        targets={'native_noiseless':native,'observed_flux':nominal['data_flux']}
        if official_sensitivity:targets['official_mean_sensitivity']=official
        fit_details=[];responses=[];predictions={};target_comparisons={};support=[]
        for target_name,target in targets.items():
            if target_name=='native_noiseless':
                base=np.zeros(4);r=white(e.flux(base)-target)
                jw=white(e.jac(base));singular=np.linalg.svd(jw,compute_uv=False)
                closure=bool(float(r@r)<1e-10 and singular[-1]>singular[0]*1e-10)
                detail=dict(gate=closure,known_exact_zero_verified=True,best_chi2=float(r@r),starts=[])
            else:
                base,detail=gate.fit_pair(e,target,cov,full,'nominal')
            base_support=support_check(e,base,'nominal')
            detail['gate']=bool(detail['gate'] and base_support['gate'])
            fit_details.append(dict(target=target_name,mode='nominal',**detail))
            support.append(dict(target=target_name,mode='nominal',**base_support))
            baseflux=e.flux(base);predictions[target_name+'__nominal']=baseflux
            jw=white(e.jac(base));inverse=np.linalg.pinv(jw,rcond=1e-10)
            changes={}
            for mode in ['observer','sed']:
                alt,adetail=gate.fit_pair(e,target,cov,full,mode)
                asupport=support_check(e,alt,mode)
                adetail['gate']=bool(adetail['gate'] and asupport['gate'])
                fit_details.append(dict(target=target_name,mode=mode,**adetail))
                support.append(dict(target=target_name,mode=mode,**asupport))
                shift=alt-base
                displacement=e.flux(base,mode)-baseflux
                finite_tangent=-inverse@white(displacement)
                infinitesimal=-inverse@white(e.flux(base,mode,derivative=True))
                fitted=e.flux(alt,mode);predictions[target_name+'__'+mode]=fitted
                changes[mode]=white(fitted-baseflux)
                responses.append(dict(CID=cid,zHEL=e.z,target=target_name,mode=mode,
                                      paired_gate=bool(detail['gate'] and adetail['gate']),
                                      baseline_theta=base.tolist(),alternative_theta=alt.tolist(),actual_shift=shift.tolist(),
                                      finite_mean_tangent_shift=finite_tangent.tolist(),infinitesimal_tangent_shift=infinitesimal.tolist(),
                                      nonlinear_standardized_mag=float(gate.D@shift),
                                      finite_mean_tangent_standardized_mag=float(gate.D@finite_tangent),
                                      infinitesimal_tangent_standardized_mag=float(gate.D@infinitesimal),
                                      nonlinear_minus_finite_tangent_mag=float(gate.D@(shift-finite_tangent)),
                                      baseline_chi2=float(detail['best_chi2']),alternative_chi2=float(adetail['best_chi2']),
                                      descriptive_profile_log_score_change=float(.5*(detail['best_chi2']-adetail['best_chi2']))))
            a,b=changes['observer'],changes['sed']
            target_comparisons[target_name]=dict(observer_information=float(a@a),sed_information=float(b@b),cross=float(a@b),difference_information=float((a-b)@(a-b)))
        result=dict(CID=cid,zHEL=e.z,epochs=len(native),status='complete',official_mean_sensitivity=official_sensitivity,
                    all_gates_pass=all(f['gate'] for f in fit_details),fit_details=fit_details,responses=responses,support=support,
                    target_comparisons=target_comparisons,
                    mean_check=dict(native_archive_relative_error=mean_error,native_minus_official_exactC_norm=float(np.linalg.norm(white(native-official))),
                                    native_minus_official_max_quoted_sigma=float(np.max(abs(native-official)/nominal['data_fluxerr'])),
                                    native_nonpositive_epochs=int(np.sum(native<=0)),official_nonpositive_epochs=int(np.sum(official<=0))),
                    input_sha256={str(objpath.relative_to(ROOT)):sha(objpath),str(cachepath.relative_to(ROOT)):sha(cachepath)})
        np.savez_compressed(OUT/'objects'/f'{cid}-predictions.npz',**predictions)
        result['prediction_sha256']=sha(OUT/'objects'/f'{cid}-predictions.npz')
    except Exception as exc:
        import traceback
        result=dict(CID=cid,status='error',all_gates_pass=False,error_type=type(exc).__name__,error=str(exc),traceback=traceback.format_exc())
    result['elapsed_seconds']=time.perf_counter()-started
    atomic_json(OUT/'objects'/f'{cid}.json',result)
    return {k:result[k] for k in ['CID','status','all_gates_pass','elapsed_seconds']}


def freeze():
    OUT.mkdir(parents=True,exist_ok=False);(OUT/'objects').mkdir()
    model,bands,paths,zp=gate.sed.fr.build_model()
    coefpath=ROOT/'runs/research_2026_09_26/sed_identification/amplitude-constrained-coefficients.npz'
    sixpath=ROOT/'runs/research_2026_09_26/sed_nonlinear_refit/cohort.csv'
    edgepath=ROOT/'runs/research_2026_09_26/sed_identification/nonpositive-model-epochs.csv'
    memberpath=ROOT/'runs/research_2026_09_26/shared_distance_response_12/contrast-membership.csv'
    cohort=pd.read_csv(BASE/'cohort.csv',dtype={'CID':str})[['CID','zHEL']]
    cohort.to_csv(OUT/'cohort.csv',index=False)
    six=pd.read_csv(sixpath,dtype={'CID':str});edge=pd.read_csv(edgepath,dtype={'CID':str})
    sensitivity=sorted(set(six.CID)|set(edge.CID));assert len(sensitivity)==14
    membership=pd.read_csv(memberpath,dtype={'CID':str})
    assert set(membership.CID)==set(cohort.CID) and membership.low_quartile.sum()==membership.high_quartile.sum()==255
    membership.to_csv(OUT/'contrast-membership.csv',index=False)
    with np.load(coefpath) as d:coeff=d['broad_drift_discovery_rms0.05'][:,3]
    with np.load(BASE/'frozen-discovery-coefficients.npz') as d:observer=d['gauge_griz']@d['basis_mean']
    np.savez_compressed(OUT/'frozen-coefficients.npz',sed_coefficients=coeff,observer_griz=observer)
    (OUT/'protocol.md').write_bytes(PROTOCOL.read_bytes());(OUT/'executed_source.py').write_bytes(Path(__file__).read_bytes())
    inputs=[Path(__file__),ROOT/'scripts/research_2026_09_26/sed_nonlinear_refit.py',ROOT/'scripts/research_2026_09_26/sed_identification.py',ROOT/'scripts/salt_dust_audit/flux_response.py',BASE/'cohort.csv',BASE/'frozen-discovery-coefficients.npz',coefpath,sixpath,edgepath,memberpath]+paths
    result=dict(frozen_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),cohort_count=len(cohort),official_mean_sensitivity_ids=sensitivity,
                source_sha256=sha(__file__),protocol_sha256=sha(OUT/'protocol.md'),
                coefficient_sha256=sha(OUT/'frozen-coefficients.npz'),
                inputs_sha256={str(p.relative_to(ROOT)):sha(p) for p in inputs},workers_max=2)
    atomic_json(OUT/'freeze.json',result)


def aggregate():
    cohort=pd.read_csv(OUT/'cohort.csv',dtype={'CID':str});membership=pd.read_csv(OUT/'contrast-membership.csv',dtype={'CID':str})
    records=[json.loads((OUT/'objects'/f'{cid}.json').read_text()) for cid in cohort.CID]
    errors=[r for r in records if r['status']!='complete']
    failed=[r['CID'] for r in records if not r['all_gates_pass']]
    responses=[row for r in records if r['status']=='complete' for row in r['responses']]
    frame=pd.DataFrame(responses);frame.to_csv(OUT/'paired-responses.csv',index=False)
    support=[dict(CID=r['CID'],**row) for r in records if r['status']=='complete' for row in r['support']]
    pd.DataFrame(support).to_csv(OUT/'support-ledger.csv',index=False)
    pd.DataFrame([dict(CID=r['CID'],zHEL=r['zHEL'],**r['mean_check']) for r in records if r['status']=='complete']).to_csv(OUT/'native-mean-checks.csv',index=False)
    pd.DataFrame([{k:r[k] for k in ['CID','status','all_gates_pass','elapsed_seconds']} for r in records]).to_csv(OUT/'timing-and-gates.csv',index=False)
    low=set(membership.loc[membership.low_quartile,'CID']);high=set(membership.loc[membership.high_quartile,'CID'])
    contrasts={};comparisons={};scores={}
    all_pass=not failed and len(records)==1020
    # No success-only aggregation: withhold headline full-cohort results on failure.
    if all_pass:
        for target in ['native_noiseless','observed_flux']:
            o=sum(r['target_comparisons'][target]['observer_information'] for r in records)
            s=sum(r['target_comparisons'][target]['sed_information'] for r in records)
            cross=sum(r['target_comparisons'][target]['cross'] for r in records)
            diff=sum(r['target_comparisons'][target]['difference_information'] for r in records)
            comparisons[target]=dict(observer_norm=float(np.sqrt(o)),sed_norm=float(np.sqrt(s)),difference_norm=float(np.sqrt(diff)),relative_squared_mismatch=float(diff/o),cosine=float(cross/np.sqrt(o*s)))
            for mode in ['observer','sed']:
                data=frame[(frame.target==target)&(frame['mode']==mode)]
                assert len(data)==1020
                contrasts[target+'__'+mode]={k:float(data.loc[data.CID.isin(high),k].mean()-data.loc[data.CID.isin(low),k].mean()) for k in ['nonlinear_standardized_mag','finite_mean_tangent_standardized_mag','infinitesimal_tangent_standardized_mag']}
                contrasts[target+'__'+mode]['nonlinear_minus_finite_tangent_abs_quantiles_50_95_99_100']=np.quantile(abs(data.nonlinear_minus_finite_tangent_mag),[.5,.95,.99,1]).tolist()
                if target=='observed_flux':scores[mode]=dict(descriptive_profile_log_score_change=float(data.descriptive_profile_log_score_change.sum()),baseline_chi2=float(data.baseline_chi2.sum()),alternative_chi2=float(data.alternative_chi2.sum()))
    details=[f for r in records if r['status']=='complete' for f in r['fit_details']]
    starts=[s for f in details for s in f.get('starts',[])]
    checks=dict(optimizer_runs=len(starts),max_fisher_gradient=max((s['fisher_gradient_norm'] for s in starts),default=None),
                max_jacobian_halfstep_relative_error=max((s['jacobian_halfstep_relative_error'] for s in starts),default=None),
                max_start_mag_disagreement=max((f.get('start_standardized_difference_mag',0) for f in details),default=None),
                max_start_chi2_disagreement=max((f.get('start_chi2_difference',0) for f in details),default=None))
    result=dict(scope='Frozen constructed nonlinear ambiguity, fixed C/mask; no physical correction or model evidence.',objects=len(records),epochs=sum(r.get('epochs',0) for r in records),all_gates_pass=all_pass,failed_ids=failed,errors=errors,
                exact_high255_low255_contrasts=contrasts,nonlinear_mean_change_comparisons=comparisons,descriptive_observed_local_fit_scores=scores,
                verification=checks,total_object_cpu_wall_seconds=float(sum(r['elapsed_seconds'] for r in records)),
                timing_quantiles_seconds=np.quantile([r['elapsed_seconds'] for r in records],[.5,.9,.99,1]).tolist(),
                official_mean_sensitivity_objects=sum(r.get('official_mean_sensitivity',False) for r in records),
                limitations=['Same independent sncosmo generator/fitter; exact SNANA covariance does not make the mean exact SNANA.',
                             'No population prior probability for the constructed SED or observer modes.',
                             'Fixed covariance and retrospective accepted masks; no covariance feedback, training, classifier/selection or BBC.',
                             'Observed objective changes are descriptive local fitted-flux scores, not integrated evidence.'])
    atomic_json(OUT/'result.json',result)
    manifest=dict(source_sha256=sha(__file__),freeze_sha256=sha(OUT/'freeze.json'),result_sha256=sha(OUT/'result.json'),
                  object_results_sha256={r['CID']:sha(OUT/'objects'/f"{r['CID']}.json") for r in records},
                  outputs_sha256={str(p.relative_to(ROOT)):sha(p) for p in OUT.iterdir() if p.is_file() and p.name!='manifest.json'})
    atomic_json(OUT/'manifest.json',manifest)
    print(json.dumps({k:result[k] for k in ['objects','epochs','all_gates_pass','failed_ids','exact_high255_low255_contrasts','nonlinear_mean_change_comparisons']},indent=2),flush=True)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--freeze',action='store_true');parser.add_argument('--run',action='store_true');parser.add_argument('--workers',type=int,default=2);args=parser.parse_args()
    assert 1<=args.workers<=2
    if args.freeze:freeze()
    if not args.run:return
    config=json.loads((OUT/'freeze.json').read_text())
    assert sha(__file__)==config['source_sha256'] and sha(OUT/'protocol.md')==config['protocol_sha256'] and sha(OUT/'frozen-coefficients.npz')==config['coefficient_sha256']
    for name,digest in config['inputs_sha256'].items():assert sha(ROOT/name)==digest,name
    cohort=pd.read_csv(OUT/'cohort.csv',dtype={'CID':str});sensitivity=set(config['official_mean_sensitivity_ids'])
    tasks=[];prior=0
    for cid in cohort.CID:
        p=OUT/'objects'/f'{cid}.json'
        if p.exists():
            saved=json.loads(p.read_text())
            for name,digest in saved.get('input_sha256',{}).items():assert sha(ROOT/name)==digest,name
            prior+=1
        else:tasks.append((cid,cid in sensitivity))
    start=time.perf_counter();done=0
    print('pending',len(tasks),'already saved',prior,'workers',args.workers,flush=True)
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers,mp_context=multiprocessing.get_context('spawn'),initializer=initialize_worker) as executor:
        futures=[executor.submit(one_object,task) for task in tasks]
        for future in concurrent.futures.as_completed(futures):
            info=future.result();done+=1
            elapsed=time.perf_counter()-start
            progress=dict(completed=prior+done,pending=len(tasks)-done,elapsed_this_invocation_seconds=elapsed,
                          estimated_remaining_seconds=elapsed/done*(len(tasks)-done),last=info,workers=args.workers)
            atomic_json(OUT/'progress.json',progress)
            if done%25==0 or not info['all_gates_pass']:
                print(json.dumps(progress),flush=True)
    aggregate()


if __name__=='__main__':main()
