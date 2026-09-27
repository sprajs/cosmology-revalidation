"""Independent Gaussian accounting plus sealed fake-cache HF consumer controls.

No observational cosmology, CAMB background, spectra, or Cobaya models execute.
Consumer tests explicitly mock qualification, target identity and distance calls;
actual files, native payload seals, cache seals and replay are exercised.
"""
import argparse
import ast
import copy
from contextlib import ExitStack
import importlib.metadata
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
from scipy.special import logsumexp, log_ndtr
from scipy.stats import norm

import calibration_holdout_bridge as bridge
from exact_correction import record_digest

ROOT=bridge.ROOT
B=bridge.base


def refuse(call):
    try:call()
    except (AssertionError,ValueError,KeyError,FileNotFoundError):return
    raise AssertionError('Malformed or unsupported fixture accepted.')


def conditional(x, hf, chi):
    # A synthetic77-dimensional proper conditional factor: one scalar normal
    # and76 independent standard-normal ordinates equal to zero.
    z=(x-.1)/.7
    ll=float(norm.logpdf(x,.1,.7)+76*norm.logpdf(0))
    return dict(noncalibrator_loglike=float(hf),noncalibrator_chi2=float(chi),
        conditional_calibrator_loglike=ll,conditional_calibrator_chi2=float(z*z),
        full_loglike=float(hf+ll),full_loglike_closure=0.,released_full_loglike_closure=0.,
        calibration_contrast_mag=float(x-.1),contrast_sigma_mag=.7,contrast_z=float(z),
        conditional_logcdf=float(log_ndtr(z)),conditional_logsf=float(log_ndtr(-z)),
        rows_noncalibrator=1580,rows_calibrator=77)


def synthetic_population(design):
    rng=np.random.default_rng(273291);n=8000
    x=rng.normal(size=n);aux=rng.normal(size=n)
    records=[];backgrounds=[];rows=[]
    for a,u in zip(x,aux):
        # Target N(0,1) prior times HF N(.2,2), all other terms cancel.
        prior=float(norm.logpdf(a)+norm.logpdf(u));hf=float(norm.logpdf(a,.2,2));chi=float(((a-.2)/2)**2)
        record=dict(status='finite',point={'w':float(a),'A_fg':float(u)},
            derived={'q0':float(a-.25),'q05':float(.8*a),'q1':float(.5*a),'j0':1.,'sn_chi2':float(a*a)},
            exact_loglikes={'released_sn':0.,'synthetic_CMB_BAO':0.},
            proposal_loglikes={'released_sn':0.,'synthetic_CMB_BAO':0.},
            exact_logpost=prior,proposal_logpost=prior,log_weight=0.)
        record['payload_sha256']=record_digest(record)
        bg={'old_SN_reconstructed_loglike':0.,'HF_loglike':hf,'HF_chi2':chi,
            'holdout':conditional(a,hf,chi),'background_calls':1}
        records.append(record);backgrounds.append(bg)
        rows.append(bridge.replace_SN(record,bg,record['exact_loglikes'],design))
    return records,backgrounds,rows,np.repeat(np.arange(4),n//4),x


def pure_controls(design,gates):
    records,backgrounds,rows,groups,x=synthetic_population(design)
    result=bridge.summarize(records,rows,groups,gates)
    assert result['qualified_under_declared_numerical_gates'],result['failed_gates']
    mu=.04;sd=np.sqrt(.8);posterior=result['posterior']['w']
    assert abs(posterior['mean']-mu)<.025 and abs(posterior['sd']-sd)<.025
    independent=norm.logpdf(x,.2,2)
    assert np.array_equal(independent,np.array([r['target_logweight']for r in rows]))
    weights=np.exp(independent-logsumexp(independent))
    assert abs(result['weight_diagnostics']['weighted_summaries_for_diagnostics']['sn_chi2']['mean']-weights@(((x-.2)/2)**2))<1e-13
    pred=result['withheld_calibration_predictive'];assert pred['status']=='qualified_predictive_precision'
    truecdf=norm.cdf((mu-.1)/np.sqrt(.7**2+sd**2))
    actualcdf=np.exp(pred['calibration_contrast_predictive_tail']['log_CDF'])
    assert abs(actualcdf-truecdf)<.015
    truell=float(norm.logpdf(.1,mu,np.sqrt(.7**2+sd**2))+76*norm.logpdf(0))
    assert abs(pred['conditional_calibrator_log_predictive_density']-truell)<.03
    altered=[]
    for record,bg in zip(records,backgrounds):
        change=copy.deepcopy(bg);c=change['holdout']
        # Keep HF density exactly fixed, while changing withheld observations
        # radically. Coherent CDF/SF and conditional Gaussian likelihood remain.
        z=c['contrast_z']+35
        c.update(conditional_calibrator_loglike=float(-.5*z*z+76*norm.logpdf(0)),
                 conditional_calibrator_chi2=float(z*z),calibration_contrast_mag=float(z*.7),
                 contrast_z=float(z),conditional_logcdf=float(log_ndtr(z)),conditional_logsf=float(log_ndtr(-z)))
        c['full_loglike']=c['noncalibrator_loglike']+c['conditional_calibrator_loglike']
        altered.append(bridge.replace_SN(record,change,record['exact_loglikes'],design))
    assert [r['target_logweight']for r in altered]==[r['target_logweight']for r in rows]
    result2=bridge.summarize(records,altered,groups,gates)
    for key in ['posterior','conditional_sign_fractions','weighted_covariance','paired_mean_changes_from_parent','weight_diagnostics']:
        assert result2[key]==result[key],key
    assert result2['withheld_calibration_predictive']!=pred
    assert result2['qualified_under_declared_numerical_gates']
    # Independent-chain disagreement can be hidden by increasingly fine
    # within-chain batches. This affects prediction precision, not HF weights.
    chain_rows=copy.deepcopy(rows)
    for row,group in zip(chain_rows,groups):
        probability=.5 if group==3 else .2
        value=row['holdout'];value['conditional_logcdf']=float(np.log(probability))
        value['conditional_logsf']=float(np.log1p(-probability))
        value['conditional_calibrator_loglike']=float(np.log(probability)-70.)
    chain_result=bridge.summarize(records,chain_rows,groups,gates)
    assert chain_result['qualified_under_declared_numerical_gates']
    assert chain_result['posterior']==result['posterior'] and chain_result['weight_diagnostics']==result['weight_diagnostics']
    chain_prediction=chain_result['withheld_calibration_predictive']
    assert chain_prediction['calibration_contrast_predictive_tail'] is None
    assert chain_prediction['conditional_calibrator_log_predictive_density'] is None
    # Independently derive the chain influence standard error from raw weights.
    integrand=np.where(groups==3,.5,.2)
    contribution=weights*integrand/(weights@integrand)
    chain_sums=np.array([(contribution[groups==g]-weights[groups==g]).sum()for g in range(4)])
    independent_chain_error=float(np.sqrt(4*np.var(chain_sums,ddof=1)))
    assert independent_chain_error>.2
    assert all(value<.2 for value in chain_prediction['predictive_precision_diagnostics']['CDF']['relative_batch_MCSE'].values())
    mismatch=copy.deepcopy(backgrounds[0]);mismatch['old_SN_reconstructed_loglike']=.001
    bad=bridge.replace_SN(records[0],mismatch,records[0]['exact_loglikes'],design)
    assert bad['status']=='source_SN_density_mismatch' and bad['target_logweight']==rows[0]['target_logweight']
    failed=bridge.summarize(records,[bad]+rows[1:],groups,gates)
    assert failed['posterior'] is None and failed['withheld_calibration_predictive'] is None
    collapsed=copy.deepcopy(rows);collapsed[0]['target_logweight']+=1000
    failed=bridge.summarize(records,collapsed,groups,gates)
    assert 'raw_weight_ESS'in failed['failed_gates'] and failed['posterior']is None and failed['withheld_calibration_predictive']is None
    # The complete production interface rejects either closure failure and
    # every nonfinite measured/derived conditional quantity before caching.
    model=SimpleNamespace(release=SimpleNamespace(z_hd_noncalibrator=np.array([.1])),evaluate=lambda d:backgrounds[0]['holdout'])
    interface=bridge.NoncalibratorInterface(model,design)
    assert interface.evaluate([1])['loglike']==rows[0]['target_SN_loglike']
    refused=0
    keys=list(backgrounds[0]['holdout'])
    for key in keys:
        if key.startswith('rows_'):continue
        broken=copy.deepcopy(backgrounds[0]['holdout']);broken[key]=float('nan')
        with patch.object(model,'evaluate',return_value=broken):refuse(lambda:interface.evaluate([1]))
        refused+=1
    for key in ['full_loglike_closure','released_full_loglike_closure']:
        broken=copy.deepcopy(backgrounds[0]['holdout']);broken[key]=1.01e-7
        with patch.object(model,'evaluate',return_value=broken):refuse(lambda:interface.evaluate([1]))
        refused+=1
    return (records,backgrounds,rows,groups),{
        'synthetic_points':len(rows),'fixed_seed':273291,
        'Gaussian_recovery':{'expected_mean':mu,'mean':posterior['mean'],'expected_sd':float(sd),'sd':posterior['sd'],
          'expected_predictive_CDF':float(truecdf),'predictive_CDF':float(actualcdf),
          'expected_log_predictive_density':truell,'log_predictive_density':pred['conditional_calibrator_log_predictive_density']},
        'independent_density_accounting_bitwise':True,'withheld_values_cannot_change_target_weights_or_HF_posterior':True,
        'predictive_precision_failure_does_not_unqualify_HF_posterior':True,
        'chain_only_predictive_failure_preserves_qualified_HF':True,
        'independent_chain_only_relative_MCSE':independent_chain_error,
        'source_closure_failure_and_collapsed_weights_withhold_all':True,'finite_and_two_closure_refusals':refused}


def cache_controls(population):
    records,backgrounds,rows,groups=population
    selected_indices=np.r_[np.arange(100),np.arange(2000,2100),np.arange(4000,4100),np.arange(6000,6100)]
    with tempfile.TemporaryDirectory(dir=ROOT/'.work')as directory:
        d=Path(directory);folder=d/'chain';native=folder/'exact';native.mkdir(parents=True)
        with patch.object(bridge,'summarize_run',side_effect=ValueError('Unqualified parent')),\
             patch.object(bridge,'validation_guard',side_effect=AssertionError('Qualifier must be first')):
            refuse(lambda:bridge.actual(folder,d/'absent',d/'refused'))
        assert not(d/'refused').exists()
        data=d/'synthetic-SN.npz'
        np.savez(data,zHD=[.1,.2,.3],zHEL=[.101,.199,.301],MU=[38.,39.,40.],covariance=np.eye(3)*.03)
        selected=[records[i]for i in selected_indices]
        chosen={records[i]['point']['w']:backgrounds[i]for i in selected_indices}
        settings={'model':'lcdm','evolution':'none','sample':'dovekie','calibration':'official_planck'}
        selection=native/'selection.json';selection.write_text(json.dumps({'settings':settings,'points':[r['point']for r in selected],'groups':np.repeat(np.arange(4),100).tolist()}))
        summary=d/'summary.json';summary.write_text(json.dumps({'selection_path':B.relative(selection)}))
        frozen={'identity':'synthetic-only','source_sha256':{},'versions':{'numpy':importlib.metadata.version('numpy')},'assets':{},'sample_sha256':B.digest(data)}
        manifest=folder/'run-0.json';manifest.write_text(json.dumps({'target_identity':frozen}))
        paths=[manifest,selection,summary]
        for i,r in enumerate(selected):
            path=native/f'{i:05d}.json';path.write_text(json.dumps(r));paths.append(path)
        parent={'qualified_under_declared_numerical_gates':True,'target_identity':frozen['identity'],
                'input_sha256':{B.relative(p):B.digest(p)for p in paths},'settings':settings}
        config={'likelihood':{'released_sn':{'data_file':str(data)},'synthetic_CMB_BAO':{}},'theory':{'camb':{'extra_args':{}}}}
        fake=SimpleNamespace(release=SimpleNamespace(data=np.zeros(1657),calibrator=np.arange(1657)<77))
        target={'identity':'synthetic-HF-only-target','assets':{},'versions':frozen['versions'],'source_sha256':{},'source_calibration_release':{'input_and_audit_sha256':{}}}
        with patch.object(bridge,'summarize_run',return_value=parent),patch.object(bridge,'validation_guard',return_value={}),\
             patch.object(B,'configurations',return_value=(config,config)),patch.object(bridge,'identify_target',return_value=target),\
             patch.object(bridge.holdout,'ReleasedHoldout',return_value=fake),\
             patch.object(bridge,'background_densities',side_effect=lambda p,*a:chosen[p['w']]):
            cache=d/'cache';one=bridge.actual(folder,summary,cache)
            assert one['posterior']is None and one['withheld_calibration_predictive']is None and 'minimum_exact_points'in one['failed_gates']
            with patch.object(bridge,'background_densities',side_effect=AssertionError('Cached replay must not compute')):
                two=bridge.actual(folder,summary,cache)
            assert one==two
            # Replay a second fresh synthetic cache with arbitrarily changed
            # withheld diagnostics. Every native/target weight stays identical.
            chosen2=copy.deepcopy(chosen)
            for value in chosen2.values():
                c=value['holdout'];c['conditional_calibrator_loglike']-=10000.;c['calibration_contrast_mag']+=20.;c['contrast_z']+=20/.7
                c['conditional_logcdf']=float(log_ndtr(c['contrast_z']));c['conditional_logsf']=float(log_ndtr(-c['contrast_z']))
                c['full_loglike']=c['noncalibrator_loglike']+c['conditional_calibrator_loglike']
            with patch.object(bridge,'background_densities',side_effect=lambda p,*a:chosen2[p['w']]):
                alternative=bridge.actual(folder,summary,d/'changed-withheld')
            assert one['weight_diagnostics']==alternative['weight_diagnostics']
            for i in range(400):
                a=json.loads((cache/f'{i:05d}.json').read_text());b=json.loads((d/'changed-withheld'/f'{i:05d}.json').read_text())
                assert a['target_logweight']==b['target_logweight'] and a['target_loglikes']==b['target_loglikes']
            point=cache/'00000.json';raw=point.read_bytes();value=json.loads(raw);value['holdout']['conditional_logcdf']-=1.;point.write_text(json.dumps(value))
            refuse(lambda:bridge.actual(folder,summary,cache));point.write_bytes(raw)
            ledger=cache/'record-hashes.json';ledgerraw=ledger.read_bytes();ledger.unlink();point.write_text(json.dumps(value))
            refuse(lambda:bridge.actual(folder,summary,cache));point.write_bytes(raw);ledger.write_bytes(ledgerraw)
            old=manifest.read_bytes();manifest.write_bytes(old+b' ');refuse(lambda:bridge.actual(folder,summary,d/'bad-parent'));manifest.write_bytes(old)
            old=data.read_bytes();data.write_bytes(b'changed');refuse(lambda:bridge.actual(folder,summary,d/'bad-data'));data.write_bytes(old)
            old=(native/'00000.json').read_bytes();bad=json.loads(old);bad['exact_logpost']+=.01;(native/'00000.json').write_text(json.dumps(bad))
            # Even a freshly updated outer mocked qualifier cannot hide a native
            # payload seal violation at the inner consumer boundary.
            pcopy=copy.deepcopy(parent);pcopy['input_sha256'][B.relative(native/'00000.json')]=B.digest(native/'00000.json')
            with patch.object(bridge,'summarize_run',return_value=pcopy):refuse(lambda:bridge.actual(folder,summary,d/'bad-native'))
            (native/'00000.json').write_bytes(old)
    return {'points':400,'qualification_target_distance_interfaces_explicitly_mocked':True,
            'native_seals_and_both_cache_seals_exercised':True,'unchanged_replay_bitwise_no_new_backgrounds':True,
            'changed_withheld_fields_same_target_weights_and_diagnostics':True,
            'changed_parent_data_native_payload_and_cached_conditional_fields_refused':True,
            'unqualified_parent_refused_before_prerequisites_and_work':True}


def configuration_controls(design):
    model=bridge.holdout.ReleasedHoldout() # Data decomposition only, no distances evaluated.
    expected=model.release.original_indices[~model.release.calibrator].tolist()
    surrogate=ROOT/'.work/unified-cosmology/inference/quartic-spectral/surrogate-quartic-v1.npz'
    outputs=[]
    for name in ['lcdm','cpl']:
        settings=dict(model=name,evolution='none',sample='dovekie',calibration='official_planck',fast_lensing=True,gpu=False,surrogate=str(surrogate))
        frozen={'configuration':B.canonical(B.parent_configuration(model=name,evolution='none',sample='dovekie',calibration='official_planck',surrogate=surrogate))}
        old,new=B.configurations(settings,frozen,design)
        left,right=map(copy.deepcopy,[B.canonical(old),B.canonical(new)])
        left['likelihood'].pop('released_sn');right['likelihood'].pop('released_sn');assert left==right
        target=bridge.identify_target(new,model,design,{})
        sn=target['configuration']['likelihood']['released_sn']
        assert sn['selected_original_row_indices']==expected and sn['rows']==1580
        assert not sn['calibrator_likelihood_in_target'] and 'external'not in sn
        assert target['sample_semantics']['name']=='pantheon_shoes_hubble_flow_holdout'
        assert not target['sample_semantics']['absolute_Cepheid_calibrator_means_in_target']
        for mapping in [target['source_sha256'],target['source_calibration_release']['input_and_audit_sha256']]:B.verify_hashes(mapping)
        for key,value in [('evolution','linear'),('sample','pantheon'),('gpu',True),('calibration','independent')]:
            refuse(lambda key=key,value=value:B.configurations(dict(settings,**{key:value}),frozen,design))
        outputs.append({'model':name,'non_SN_priors_and_physics_identical':True,'target_identity':target['identity'],'rows':1580,'withheld_rows':77})
    return outputs


def background_controls(design):
    import camb
    zold=np.array([.32,.11,.32]);zhel=np.array([.321,.109,.322]);znew=np.resize([.27,.07,.27,.41],1580)
    observed=np.array([40.,38.,40.01]);cov=np.array([[.04,.006,.002],[.006,.06,.004],[.002,.004,.05]])
    old=B.IntegratedLuminosity(zold,zhel,observed,cov);requested=[];received=[];kw=[]
    def distance(z):return 4317*np.asarray(z)/(1+np.asarray(z))
    def get_distance(z):requested.append(np.array(z));return distance(z)
    def evaluate(da):received.append(np.array(da));return conditional(.2,-20.,3.)
    def parameters(**kwargs):kw.append(kwargs);return SimpleNamespace(WantTransfer=False)
    def background(p):assert p.WantTransfer;return SimpleNamespace(angular_diameter_distance=get_distance)
    model=SimpleNamespace(release=SimpleNamespace(z_hd_noncalibrator=znew),evaluate=evaluate)
    point=dict(H0=69.,ombh2=.0224,omch2=.12,ns=.97,tau=.055,logA=3.04,w=-.92,wa=.17)
    with patch.object(camb,'set_params',side_effect=parameters),patch.object(camb,'get_background',side_effect=background),\
         patch.object(B,'native_thermal_parameters',return_value=SimpleNamespace(WantTransfer=True)):
        got=bridge.background_densities(point,{'AccuracyBoost':1.},old,model,design)
    assert len(requested)==len(received)==len(kw)==1
    assert np.array_equal(requested[0],np.unique(np.r_[zold,znew])) and np.array_equal(received[0],distance(znew))
    r=5*np.log10(distance(zold)*(1+zold)*(1+zhel))+25-observed;P=np.linalg.inv(cov);one=np.ones(3)
    chi=r@P@r-(r@P@one)**2/(one@P@one)
    ll=-.5*(chi+np.linalg.slogdet(cov)[1]+np.log(one@P@one)+2*np.log(2*np.pi))
    error=abs(ll-got['old_SN_reconstructed_loglike']);assert error<1e-11
    assert got['HF_loglike']==-20. and got['HF_chi2']==3. and 'anchored_SN'not in got
    assert got['distance_unit']=='Mpc' and got['H0_unit']=='km/s/Mpc' and got['background_calls']==1
    assert abs(kw[0]['As']-1e-10*np.exp(point['logA']))<1e-25
    return {'fake_background_calls':1,'real_background_calls':0,'ordered1580_rows_restored':True,'unique_union_redshifts':requested[0].tolist(),
            'old_SN_independent_GLS_loglike_error':float(error),'only_HF_factor_reaches_weight_kernel':True}


def prerequisite_controls():
    # Exercise the real guard with tiny synthetic report files and live source
    # hashes; no production report is rewritten or accepted as observational.
    with tempfile.TemporaryDirectory(dir=ROOT/'.work')as td:
        d=Path(td);kr=d/'kernel.json';vr=d/'bridge.json'
        src=[Path(bridge.__file__),Path(bridge.holdout.__file__),bridge.DESIGN,Path(__file__),bridge.HERE/'calibration_holdout_review.py']
        hashes={B.relative(p):B.digest(p)for p in src}
        kernel={'status':'passed_independent_calibration_holdout_review','source_sha256':hashes,'input_sha256':{}}
        kr.write_text(json.dumps(kernel))
        own={'status':'passed_synthetic_calibration_holdout_bridge_checks_no_observational_evaluation',
             'source_sha256':hashes,'input_sha256':{B.relative(kr):B.digest(kr)},'CMB_spectrum_calls':0,'CAMB_background_calls':0,'observational_bridge_points':0}
        vr.write_text(json.dumps(own))
        with patch.object(bridge,'KERNEL_REVIEW',kr),patch.object(bridge,'VALIDATION',vr):
            assert len(bridge.validation_guard())==2
            for kind in ['status','source','kernel','calls','missing']:
                changed=copy.deepcopy(own)
                if kind=='status':changed['status']='pending'
                if kind=='source':changed['source_sha256'][B.relative(Path(bridge.__file__))]='0'*64
                if kind=='kernel':changed['input_sha256'][B.relative(kr)]='0'*64
                if kind=='calls':changed['CMB_spectrum_calls']=1
                if kind=='missing':changed['source_sha256'].pop(B.relative(Path(__file__)))
                vr.write_text(json.dumps(changed));refuse(bridge.validation_guard)
            vr.write_text(json.dumps(own));kr.unlink();refuse(bridge.validation_guard)
    return {'live_source_bindings_required':True,'pending_stale_missing_report_and_changed_kernel_refused':True,'refusals':6}


def validate():
    import camb,cobaya.model
    attempted=[]
    def forbidden(*a,**k):attempted.append(True);raise AssertionError('Real physical call forbidden in validation.')
    source_paths=[Path(__file__),Path(bridge.__file__),Path(bridge.holdout.__file__),bridge.DESIGN,bridge.GATES,
        bridge.HERE/'calibration_holdout_review.py',Path(B.__file__),B.DESIGN,
        bridge.HERE/'anchored_adapter.py',bridge.HERE/'anchored-design.json',bridge.HERE/'expansion_history.py',
        bridge.HERE/'luminosity_sensitivity.py',bridge.HERE/'measurement_summary.py',bridge.HERE/'probe_omission.py',
        bridge.HERE/'exact_correction.py',bridge.HERE/'target_identity.py',bridge.HERE/'modern_fast.py',bridge.HERE/'modern_run.py',
        bridge.HERE.parent/'distance_ladder/calibration_interface.py']
    sources={B.relative(p):B.digest(p)for p in source_paths}
    for p in source_paths:
        if p.suffix=='.py':ast.parse(p.read_text())
    review=json.loads(bridge.KERNEL_REVIEW.read_text())
    assert review['status']=='passed_independent_calibration_holdout_review'
    B.verify_hashes(review['source_sha256']);B.verify_hashes(review['input_sha256'])
    with ExitStack()as stack:
        for name in ['set_params','get_background','get_results','get_transfer_functions']:stack.enter_context(patch.object(camb,name,forbidden))
        stack.enter_context(patch.object(cobaya.model,'get_model',forbidden));stack.enter_context(patch.object(cobaya.model.Model,'__init__',forbidden))
        design=json.loads(bridge.DESIGN.read_text());gates=json.loads(bridge.GATES.read_text())['overlap_gates']
        population,pure=pure_controls(design,gates)
        cache=cache_controls(population)
        config=configuration_controls(design)
        background=background_controls(design)
        prerequisite=prerequisite_controls()
    assert not attempted
    B.verify_hashes(sources);B.verify_hashes(review['source_sha256']);B.verify_hashes(review['input_sha256'])
    return {'status':'passed_synthetic_calibration_holdout_bridge_checks_no_observational_evaluation',
        'Gaussian_and_isolation_controls':pure,'synthetic_consumer_cache':cache,'live_configuration_identity_controls':config,
        'fake_background_operator':background,'prerequisite_guard':prerequisite,
        'source_sha256':sources,'input_sha256':{B.relative(bridge.KERNEL_REVIEW):B.digest(bridge.KERNEL_REVIEW),**review['input_sha256']},
        'observational_bridge_points':0,'CAMB_background_calls':0,'CMB_spectrum_calls':0,
        'physical_entrypoints_forbidden':True,
        'scope':'Source-bound synthetic accounting and cache/qualification integrity only. Real source-SN closure, overlap and predictive contribution gates remain required; no observed calibration tail or cosmological measurement.'}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    result=validate();a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':result['status'],'native_calls':0}))


if __name__=='__main__':main()
