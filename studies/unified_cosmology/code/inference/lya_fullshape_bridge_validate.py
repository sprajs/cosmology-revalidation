"""Synthetic-only Lyalpha replacement, no-background consumer integrity checks."""
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
from scipy.special import logsumexp
from scipy.stats import norm

import lya_fullshape_bridge as bridge
from exact_correction import record_digest

ROOT=bridge.ROOT


def refuse(call):
    try:call()
    except(AssertionError,ValueError,KeyError,FileNotFoundError):return
    raise AssertionError('Malformed fixture was accepted.')


def pure_controls(design,gates):
    rng=np.random.default_rng(273293);n=8000;groups=np.repeat(np.arange(4),n//4)
    x=rng.normal(size=n);aux=rng.normal(size=n)
    records=[];rows=[]
    for a,u in zip(x,aux):
        prior=float(norm.logpdf(a)+norm.logpdf(u));old=-.5*a*a;new=-.5*((a-.2)/2)**2
        # Proposal is N(0,1); source BAO is cancelled in the proposal residual
        # component so component+prior accounting remains exact in this toy.
        record={'status':'finite','point':{'w':float(a),'A_fg':float(u)},
            'derived':{'q0':float(a-.25),'q05':float(.8*a),'q1':float(.5*a),'j0':1.,'sn_chi2':float(a*a+.1)},
            'exact_loglikes':{bridge.COMPONENT:float(old),'synthetic_CMB_SN':0.},
            'proposal_loglikes':{bridge.COMPONENT:0.,'synthetic_CMB_SN':0.},
            'exact_logpost':float(prior+old),'proposal_logpost':prior,'log_weight':float(old)}
        record['payload_sha256']=record_digest(record)
        row=bridge.replace_BAO(record,old,new,record['exact_loglikes'],design)
        records.append(record);rows.append(row)
    result=bridge.summarize_variant(records,rows,groups,gates)
    assert result['qualified_under_declared_numerical_gates'],result['failed_gates']
    lw=-.5*((x-.2)/2)**2
    error=float(np.max(abs(lw-[r['target_logweight']for r in rows])));assert error<1e-15
    p=result['posterior']['w'];assert abs(p['mean']-.04)<.025 and abs(p['sd']-np.sqrt(.8))<.025
    weights=np.exp(lw-logsumexp(lw))
    assert abs(result['weight_diagnostics']['weighted_summaries_for_diagnostics']['sn_chi2']['mean']-weights@(x*x+.1))<1e-13
    for record,row in zip(records[:100],rows[:100]):
        old=record['exact_loglikes'][bridge.COMPONENT]
        null=bridge.replace_BAO(record,old,old,record['exact_loglikes'],design)
        assert null['target_logweight']==record['log_weight']
        bad=bridge.replace_BAO(record,old+.001,row['target_BAO_loglike'],record['exact_loglikes'],design)
        assert bad['status']=='source_BAO_density_mismatch'and bad['target_logweight']==row['target_logweight']
    shifted=np.array([bridge.replace_BAO(r,r['exact_loglikes'][bridge.COMPONENT],row['target_BAO_loglike']+700,r['exact_loglikes'],design)['target_logweight']for r,row in zip(records,rows)])
    constant=float(np.max(abs(np.exp(shifted-logsumexp(shifted))-weights)));assert constant<1e-13
    bad=copy.deepcopy(rows);bad[0]['status']='source_BAO_density_mismatch'
    fail=bridge.summarize_variant(records,bad,groups,gates);assert fail['posterior']is None and not fail['qualified_under_declared_numerical_gates']
    bad=copy.deepcopy(rows);bad[0]['target_logweight']+=1000
    fail=bridge.summarize_variant(records,bad,groups,gates);assert fail['posterior']is None and 'raw_weight_ESS'in fail['failed_gates']
    r=copy.deepcopy(records[0]);r['proposal_logpost']+=.01;r['log_weight']-=.01
    refuse(lambda:bridge.replace_BAO(r,0.,0.,r['exact_loglikes'],design))
    for value in [float('nan'),float('inf')]:refuse(lambda:bridge.replace_BAO(records[0],value,0.,records[0]['exact_loglikes'],design))
    return records,rows,groups,{'points':n,'seed':273293,'maximum_direct_density_error':error,
        'Gaussian_mean':p['mean'],'expected_mean':.04,'Gaussian_sd':p['sd'],'expected_sd':float(np.sqrt(.8)),
        'normalized_constant_shift_error':constant,'SN_chi2_preserved':True,
        'source_closure_failure_not_substituted_into_weights':True,'all_failed_support_withholds_posterior':True}


def configuration_controls(design):
    surrogate=ROOT/'.work/unified-cosmology/inference/quartic-spectral/surrogate-quartic-v1.npz'
    cases=[]
    for model,evolution,gpu in [('lcdm','none',False),('cpl','none',False),('cpl','linear',True),('cpl','smooth01',True)]:
        settings=dict(model=model,evolution=evolution,gpu=gpu,fast_lensing=True,sample='dovekie',calibration='official_planck',surrogate=str(surrogate))
        factory=bridge.configuration_factory(settings);options={k:settings[k]for k in ['model','evolution','sample','calibration']}
        frozen={'configuration':bridge.canonical(factory(**options,surrogate=surrogate))}
        native,target=bridge.configurations(settings,frozen,design)
        left,right=copy.deepcopy(native),copy.deepcopy(target);left['likelihood'].pop(bridge.COMPONENT);right['likelihood'].pop(bridge.COMPONENT)
        assert left==right and target['likelihood'][bridge.COMPONENT]['replaces_existing_z2p33_pair']
        bad=copy.deepcopy(frozen);bad['configuration']['params']['H0']['prior']['max']+=1
        refuse(lambda:bridge.configurations(settings,bad,design))
        for key,value in [('sample','pantheon'),('calibration','paper_literal'),('evolution','smooth03')]:
            refuse(lambda key=key,value=value:bridge.configurations(dict(settings,**{key:value}),frozen,design))
        changed=copy.deepcopy(native);changed['likelihood'][bridge.COMPONENT]={'cov_file':'other'}
        with patch.object(bridge,'configuration_factory',return_value=lambda **kw:copy.deepcopy(frozen['configuration'])if kw.get('surrogate')else changed):
            refuse(lambda:bridge.configurations(settings,frozen,design))
        cases.append({'model':model,'evolution':evolution,'gpu':gpu,'non_BAO_physics_priors_and_likelihoods_identical':True})
    return cases


def background_controls(design):
    import camb
    model=bridge.kernel.GaussianBAOReplacement.from_release()
    req=[];params=[]
    def da(z):return 4100*np.asarray(z)/(1+np.asarray(z))
    def hub(z):return 70*(1+np.asarray(z))**1.45
    def background(p):
        assert p.WantTransfer
        def angular(z):req.append(np.array(z));return da(z)
        return SimpleNamespace(angular_diameter_distance=angular,hubble_parameter=hub,get_derived_params=lambda:{'rdrag':147.1})
    def set_params(**kw):params.append(kw);return SimpleNamespace(WantTransfer=False)
    point=dict(H0=69.,ombh2=.0224,omch2=.12,ns=.97,tau=.055,logA=3.04,w=-.92,wa=.17)
    with patch.object(camb,'set_params',side_effect=set_params),patch.object(camb,'get_background',side_effect=background),\
         patch.object(bridge,'native_thermal_parameters',return_value=SimpleNamespace(WantTransfer=True)):
        p=bridge.background_prediction(point,{'AccuracyBoost':1.},model.z)
    assert len(req)==len(params)==1 and np.array_equal(req[0],np.unique(model.z))
    assert np.array_equal(p['DM_Mpc'],da(model.z)*(1+model.z))and np.array_equal(p['H_km_s_Mpc'],hub(model.z))
    got=model.evaluate(p['DM_Mpc'],p['H_km_s_Mpc'],p['rdrag_Mpc']);bridge.evaluation_checks(got,model.new,design)
    prediction=[]
    for z,dm,h,typ in zip(model.z,p['DM_Mpc'],p['H_km_s_Mpc'],model.row_types):
        dh=299792.458/h
        prediction.append({'DM_over_rs':dm/147.1,'DH_over_rs':dh/147.1,'DV_over_rs':(z*dm*dm*dh)**(1/3)/147.1}[typ])
    difference=np.asarray(prediction)-model.mean
    old=-.5*difference@np.linalg.inv(model.C)@difference
    assert abs(got['old_loglike']-old)<1e-10
    assert len(got['variants'])==33
    maxima=0.
    for name,numbers in bridge.kernel.variants().items():
        dm,dh,sm,sh,rho=numbers;C=model.C.copy();C[np.ix_(model.order,model.order)]=np.outer([sm,sh],[sm,sh])*[[1,rho],[rho,1]]
        mean=model.mean.copy();mean[model.order]=[dm,dh];r=np.asarray(prediction)-mean
        direct=-.5*r@np.linalg.inv(C)@r
        maxima=max(maxima,abs(got['variants'][name]['new_loglike']-direct))
    assert maxima<1e-10
    for field in ['old_loglike','old_block_closure']:
        bad=copy.deepcopy(got);bad[field]=float('nan');refuse(lambda:bridge.evaluation_checks(bad,model.new,design))
    bad=copy.deepcopy(got);bad['variants']['nominal']['new_block_closure']=2e-9
    refuse(lambda:bridge.evaluation_checks(bad,model.new,design))
    return model,{'fake_background_calls':1,'real_background_calls':0,'ordered13_rows_and_redshift_union_pass':True,
        'old_native_Cobaya_quadratic_independent_error':float(abs(got['old_loglike']-old)),
        'all33_full_covariance_independent_error':float(maxima),'actual_DM_DH_row_indices':model.order.tolist(),
        'all_variants_reuse_identical_predictions':True}


def cache_controls(records,model):
    index=np.r_[np.arange(100),np.arange(2000,2100),np.arange(4000,4100),np.arange(6000,6100)]
    selected=[records[i]for i in index]
    with tempfile.TemporaryDirectory(dir=ROOT/'.work')as td:
        d=Path(td);folder=d/'chain';native=folder/'exact';native.mkdir(parents=True)
        with patch.object(bridge,'summarize_run',side_effect=ValueError('unqualified')),patch.object(bridge,'validation_guard',side_effect=AssertionError('too early')):
            refuse(lambda:bridge.actual(folder,d/'absent',d/'refused'))
        assert not(d/'refused').exists()
        settings={'model':'lcdm','evolution':'none','sample':'dovekie','calibration':'official_planck','fast_lensing':True}
        selection=native/'selection.json';selection.write_text(json.dumps({'settings':settings,'points':[r['point']for r in selected],'groups':np.repeat(np.arange(4),100).tolist()}))
        summary=d/'summary.json';summary.write_text(json.dumps({'selection_path':bridge.relative(selection)}))
        frozen={'identity':'synthetic-parent','source_sha256':{},'versions':{'numpy':importlib.metadata.version('numpy')}}
        manifest=folder/'run-0.json';manifest.write_text(json.dumps({'target_identity':frozen}));paths=[manifest,selection,summary]
        for i,r in enumerate(selected):
            path=native/f'{i:05d}.json';path.write_text(json.dumps(r));paths.append(path)
        parent={'qualified_under_declared_numerical_gates':True,'target_identity':frozen['identity'],'settings':settings,
            'input_sha256':{bridge.relative(p):bridge.digest(p)for p in paths}}
        config={'likelihood':{bridge.COMPONENT:{},'synthetic_CMB_SN':{}},'theory':{'camb':{'extra_args':{}}}}
        target={'identity':'synthetic-child','source_sha256':{},'versions':frozen['versions']}
        calls=[]
        def prediction(point,*a):calls.append(point);return {'DM_Mpc':[point['w']],'H_km_s_Mpc':[1.],'rdrag_Mpc':1.,'background_calls':1}
        def evaluate(dm,*a):
            x=dm[0];old=-.5*x*x
            return {'old_loglike':old,'old_chi2':-2*old,'retained_chi2':0.,'old_lya_chi2':-2*old,'old_block_closure':0.,'prediction':[x],
                'variants':{name:{'new_loglike':-.5*((x-.2)/2)**2-i*1e-5,'new_chi2':((x-.2)/2)**2+2*i*1e-5,
                                  'new_lya_chi2':((x-.2)/2)**2+2*i*1e-5,'new_block_closure':0.,'delta_loglike':-.5*((x-.2)/2)**2-i*1e-5-old}
                            for i,name in enumerate(model.new)}}
        fake=SimpleNamespace(z=model.z,new=model.new,evaluate=evaluate)
        with patch.object(bridge,'summarize_run',return_value=parent),patch.object(bridge,'validation_guard',return_value={}),\
             patch.object(bridge,'configurations',return_value=(config,config)),patch.object(bridge,'load_released',return_value=(fake,{})),\
             patch.object(bridge,'identify_target',return_value=target),patch.object(bridge,'background_prediction',side_effect=prediction):
            cache=d/'cache';one=bridge.actual(folder,summary,cache)
            assert len(calls)==400 and len(one['variants'])==33
            assert all(not value['qualified_under_declared_numerical_gates']and value['posterior']is None for value in one['variants'].values())
            with patch.object(bridge,'background_prediction',side_effect=AssertionError('No new background for replay')):two=bridge.actual(folder,summary,cache)
            assert one==two
            row=cache/'00000.json';raw=row.read_bytes();altered=json.loads(raw);altered['variants']['nominal']['target_logweight']+=.1;row.write_text(json.dumps(altered))
            refuse(lambda:bridge.actual(folder,summary,cache));row.write_bytes(raw)
            ledger=cache/'record-hashes.json';rawledger=ledger.read_bytes();ledger.unlink();row.write_text(json.dumps(altered))
            refuse(lambda:bridge.actual(folder,summary,cache));row.write_bytes(raw);ledger.write_bytes(rawledger)
            before=manifest.read_bytes();manifest.write_bytes(before+b' ');refuse(lambda:bridge.actual(folder,summary,d/'changed-parent'));manifest.write_bytes(before)
            before=(native/'00000.json').read_bytes();r=json.loads(before);r['exact_logpost']+=.01;(native/'00000.json').write_text(json.dumps(r))
            fresh=copy.deepcopy(parent);fresh['input_sha256'][bridge.relative(native/'00000.json')]=bridge.digest(native/'00000.json')
            with patch.object(bridge,'summarize_run',return_value=fresh):refuse(lambda:bridge.actual(folder,summary,d/'changed-native'))
    return {'synthetic_records':400,'variant_targets':33,'background_calls_for_all_variants':400,
        'cache_replay_bitwise_no_extra_calls':True,'changed_parent_native_cache_ledger_and_payload_refused':True,
        'qualification_target_identity_distance_interfaces_explicitly_mocked':True,'all_underpowered_variants_withhold_posterior':True,
        'parent_refused_before_prerequisites_or_work':True}


def prerequisite_controls():
    paths=[Path(bridge.__file__),Path(__file__),Path(bridge.kernel.__file__),bridge.DESIGN]
    sources={bridge.relative(p):bridge.digest(p)for p in paths}
    with tempfile.TemporaryDirectory(dir=ROOT/'.work')as td:
        path=Path(td)/'validation.json'
        report={'status':'passed_synthetic_lya_fullshape_bridge_checks_no_observational_evaluation',
            'source_sha256':sources,'input_sha256':{bridge.relative(p):bridge.digest(p)for p in [bridge.KERNEL_VALIDATION,bridge.SOURCES]},
            'CMB_spectrum_calls':0,'CAMB_background_calls':0,'observational_bridge_points':0}
        path.write_text(json.dumps(report))
        with patch.object(bridge,'VALIDATION',path):
            assert len(bridge.validation_guard())==3
            for kind in ['status','missing_source','stale_source','stale_kernel','physical_calls']:
                bad=copy.deepcopy(report)
                if kind=='status':bad['status']='pending'
                if kind=='missing_source':bad['source_sha256'].pop(bridge.relative(Path(__file__)))
                if kind=='stale_source':bad['source_sha256'][bridge.relative(Path(bridge.__file__))]='0'*64
                if kind=='stale_kernel':bad['input_sha256'][bridge.relative(bridge.KERNEL_VALIDATION)]='0'*64
                if kind=='physical_calls':bad['observational_bridge_points']=1
                path.write_text(json.dumps(bad));refuse(bridge.validation_guard)
            path.unlink();refuse(bridge.validation_guard)
    return {'source_report_and_kernel_bindings_required':True,'stale_incomplete_pending_or_missing_prerequisites_refused':True,'refusals':6}


def validate():
    import camb,cobaya.model
    source_paths=[Path(__file__),Path(bridge.__file__),Path(bridge.kernel.__file__),bridge.DESIGN,bridge.GATES,
        Path(bridge.summary_kernel.__file__),bridge.HERE/'probe_omission.py',bridge.HERE/'measurement_summary.py',
        bridge.HERE/'exact_correction.py',bridge.HERE/'target_identity.py',bridge.HERE/'expansion_history.py',
        bridge.HERE/'modern_run.py',bridge.HERE/'modern_fast.py',bridge.HERE/'modern_gpu.py',bridge.HERE/'luminosity_sensitivity.py']
    sources={bridge.relative(p):bridge.digest(p)for p in source_paths}
    _,inputs=bridge.load_released()
    for p in source_paths:
        if p.suffix=='.py':ast.parse(p.read_text())
    kernel_report=json.loads(bridge.KERNEL_VALIDATION.read_text())
    # Exact required status is fixed after the independently authored kernel
    # validator has frozen; no observational success is inferred from this.
    assert kernel_report['status']=='passed'
    assert all(kernel_report['source_sha256'].get(bridge.relative(p))==bridge.digest(p)for p in [Path(bridge.kernel.__file__),bridge.DESIGN])
    bridge.verify_hashes(kernel_report['source_sha256'])
    bridge.verify_hashes(kernel_report['input_sha256'])
    attempted=[]
    def forbidden(*a,**k):attempted.append(True);raise AssertionError('No physical calls in synthetic validation.')
    with ExitStack()as stack:
        for name in ['set_params','get_background','get_results','get_transfer_functions']:stack.enter_context(patch.object(camb,name,forbidden))
        stack.enter_context(patch.object(cobaya.model,'get_model',forbidden));stack.enter_context(patch.object(cobaya.model.Model,'__init__',forbidden))
        design=json.loads(bridge.DESIGN.read_text());gates=json.loads(bridge.GATES.read_text())['overlap_gates']
        records,rows,groups,pure=pure_controls(design,gates)
        config=configuration_controls(design)
        model,background=background_controls(design)
        cache=cache_controls(records,model)
        prerequisite=prerequisite_controls()
    assert not attempted;bridge.verify_hashes(sources);bridge.verify_hashes(inputs)
    return {'status':'passed_synthetic_lya_fullshape_bridge_checks_no_observational_evaluation',
        'pure_accounting_and_Gaussian_controls':pure,'configuration_controls':config,'fake_background_controls':background,'sealed_consumer_controls':cache,
        'prerequisite_guard_controls':prerequisite,
        'source_sha256':sources,'input_sha256':inputs,'CMB_spectrum_calls':0,'CAMB_background_calls':0,'observational_bridge_points':0,
        'physical_entrypoints_forbidden':True,
        'scope':'Synthetic accounting/configuration/cache controls only; all real parent, source BAO closure and33 separate overlap qualifications remain required.'}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    result=validate();args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':result['status'],'CMB_spectrum_calls':0}))


if __name__=='__main__':main()
