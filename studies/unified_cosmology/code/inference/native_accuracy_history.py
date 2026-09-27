"""Expansion and grey luminosity from a fresh typed native2 posterior."""
import argparse
from contextlib import contextmanager
import copy
import importlib.metadata
import json
from pathlib import Path
import time
from unittest.mock import patch
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.linalg import cho_solve
from scipy.special import logsumexp
import expansion_history as expansion
import luminosity_history as luminosity
from luminosity_sensitivity import IntegratedLuminosity, weight_diagnostics
from likelihood import ReleasedDistances
from measurement_summary import weighted_fraction
from native_accuracy_qualify import QualifiedNativeAccuracy2, qualify_run
import native_accuracy_runtime as rt

ROOT=rt.ROOT;HERE=Path(__file__).resolve().parent
DESIGN=HERE/'native-accuracy-history-design.json'
VALIDATION=ROOT/'studies/unified_cosmology/results/inference/native-accuracy-history-validation.json'


def sources():
    paths=[Path(__file__),DESIGN,HERE/'native_accuracy_history_validate.py',HERE/'native_accuracy_qualify.py',
           HERE/'native_accuracy_runtime.py',HERE/'expansion_history.py',expansion.DESIGN,
           HERE/'luminosity_history.py',HERE/'luminosity_sensitivity.py',HERE/'luminosity-sensitivity-design.json',
           HERE/'likelihood.py',HERE/'measurement_summary.py',HERE/'exact_correction.py']
    return {rt.relative(p):rt.digest(p) for p in paths}


def validation_guard():
    report=json.loads(VALIDATION.read_text())
    assert report['status']=='passed_synthetic_native_accuracy2_history_validation'
    assert report['source_sha256']==sources() and report['physical_calls']==0
    return {rt.relative(VALIDATION):rt.digest(VALIDATION)}


def typed_parent(work,summary):
    target=qualify_run(work,summary)
    assert type(target) is QualifiedNativeAccuracy2 and target.numerical_accuracy==2
    assert target.posterior_summary['status']=='qualified_native_accuracy2_importance_posterior'
    assert len(target.records)==len(target.selected_points)==len(target.locations)==len(target.logweights)==2000
    assert np.array_equal(target.groups,np.repeat(np.arange(4),500))
    weights=np.exp(target.logweights-logsumexp(target.logweights))
    assert np.array_equal(weights,target.normalized_weights)
    for i,row in enumerate(target.records):
        assert row['index']==i and row['point']==target.selected_points[i]
        assert row['log_weight']==target.logweights[i] and row['proposal_logpost']==target.proposal_logposts[i]
        assert row['numerical_target_identity']==target.numerical_target_identity
        assert row['status']=='finite_native_accuracy2'
    return target


class LuminosityKernel:
    def __init__(self,sn,evolution,redshifts):
        assert evolution in {'none','linear','smooth01','smooth03'}
        self.sn=sn;self.evolution=evolution;self.z=np.asarray(redshifts,dtype=float)
        self.basis=CubicSpline(sn.knots,np.eye(5)[:,1:],bc_type='natural')(self.z)
        self.coefficient_covariance=None;self.covariance=np.zeros((len(self.z),len(self.z)))
        if evolution.startswith('smooth'):
            self.coefficient_covariance=cho_solve(sn.gaussian[evolution][0],np.eye(4))
            self.covariance=self.basis@self.coefficient_covariance@self.basis.T
        assert np.max(abs(self.basis[self.z==0]))==0

    def from_prediction(self,prediction,point):
        sn=self.sn;residual=sn.project(np.asarray(prediction)-sn.observed)
        coefficient_mean=None
        if self.evolution.startswith('smooth'):
            assert point.get('epsilon',0.)==0.
            score=sn.spline.T@residual
            coefficient_mean=-self.coefficient_covariance@score
            mean=self.basis@coefficient_mean
            chi2=float(residual@residual-score@self.coefficient_covariance@score)
            lognorm=sn.lognorm+sn.gaussian[self.evolution][1]
        else:
            epsilon=point.get('epsilon',0.)
            if self.evolution=='none':assert epsilon==0.
            assert np.isfinite(epsilon) and -.5<=epsilon<=.5
            residual=sn.project(np.asarray(prediction)+epsilon*np.log1p(sn.z)/np.log(2)-sn.observed)
            mean=epsilon*np.log1p(self.z)/np.log(2)
            chi2=float(residual@residual);lognorm=sn.lognorm
        return {'mean_mag':mean.tolist(),'coefficient_mean_mag':None if coefficient_mean is None else coefficient_mean.tolist(),
                'SN_loglike':float(-.5*(chi2+lognorm)),'SN_chi2':chi2}


def load_sn(target,design):
    config=target.native_configuration['likelihood']['released_sn']
    assert config['external'] is ReleasedDistances,'Anchored or alternate SN likelihood is outside this consumer.'
    evolution=target.parent_settings['evolution'];sigma={'none':0.,'linear':0.,'smooth01':.1,'smooth03':.3}[evolution]
    assert config['smooth_sigma']==sigma
    epsilon=target.native_configuration['params']['epsilon']
    if evolution=='linear':assert epsilon['prior']=={'min':-.5,'max':.5}
    else:assert epsilon==0.
    path=Path(config['data_file']);path=path if path.is_absolute() else ROOT/path
    assert rt.relative(path) in target.input_sha256 and rt.digest(path)==target.input_sha256[rt.relative(path)]
    with np.load(path,allow_pickle=False) as data:
        sn=IntegratedLuminosity(*(data[k] for k in ['zHD','zHEL','MU','covariance']))
    assert np.all(sn.z>0) and np.isfinite(sn.zhel).all() and np.all(1+sn.zhel>0)
    kernel=LuminosityKernel(sn,evolution,design['luminosity_redshifts'])
    kernel.declared_extra=copy.deepcopy(target.native_configuration['theory']['camb']['extra_args'])
    return kernel,path


def background_options(native,declared):
    actual=native['value']['metadata'];extra=copy.deepcopy(actual['finalized_theory_extra_args'])
    # Cobaya consumes halofit_version separately when building CAMBParams. It
    # also increases lmax to satisfy the finalized likelihood requests. Restore
    # the consumed declared option; never infer max_l from padded Cl length.
    assert set(declared)-set(extra) <= {'halofit_version'}
    for k,v in declared.items():
        if k=='lmax':assert extra[k]>=v
        elif k=='halofit_version' and k not in extra:extra[k]=v
        else:assert rt.canonical(extra[k])==rt.canonical(v),(k,'changed native setting')
    assert extra['omk']==0 and extra['dark_energy_model']=='ppf'
    for name in ['AccuracyBoost','lAccuracyBoost','lSampleBoost']:assert extra[name]==2
    return extra


def configuration_checks(target):
    assert set(target.native_configuration['theory'])=={'camb'}
    declared=target.native_configuration['theory']['camb']['extra_args']
    for row in target.records:
        background_options(row,declared)
        actual=row['value']['metadata']
        assert actual['ell_factor'] is False and actual['units']=='FIRASmuK2'
        assert actual['pp']=='dimensionless C_phi_phi' and actual['CAMB_Params_max_l']>0


def numerical_history(background,native,design):
    """Same kernels and numerical thresholds as the frozen original history."""
    result=expansion.background_history(background,design['redshift_grid'],design['finite_difference_step'])
    z=np.array(design['redshift_grid']);checks=result['numerical_checks']
    for name in ['q','j']:
        checks[name+'_native_scaled_error']=float(max(abs(result['history'][name][int(np.flatnonzero(z==zz)[0])]-native[name+label])/(1+abs(native[name+label]))
                    for zz,label in [(0.,'0'),(.5,'05'),(1.,'1')]))
    checks['rdrag_relative_error']=abs(result['scalar']['rdrag_Mpc']/native['rdrag']-1)
    checks['omegam_absolute_error']=abs(result['scalar']['omegam']-native['omegam'])
    fields=[('q_step_scaled_error','q_step_and_native_scaled_error_max'),('q_native_scaled_error','q_step_and_native_scaled_error_max'),
       ('j_step_scaled_error','j_step_and_native_scaled_error_max'),('j_native_scaled_error','j_step_and_native_scaled_error_max'),
       ('rdrag_relative_error','rdrag_relative_error_max'),('omegam_absolute_error','omegam_absolute_error_max'),
       ('friedmann_relative_error','friedmann_relative_error_max'),('age_internal_absolute_Gyr','age_internal_absolute_Gyr_max'),
       ('age_integral_absolute_Gyr_error','age_integral_absolute_Gyr_error_max')]
    result['failed_numerical_gates']=[k for k,g in fields if not np.isfinite(checks[k]) or checks[k]>design['numerical_gates'][g]]
    return result


def consume_background(background,point,native,kernel,design,expansion_design):
    result=numerical_history(background,native['native2_derived'],expansion_design)
    metadata=native['value']['metadata'];assert background.Params.max_l==metadata['CAMB_Params_max_l']
    unique,inverse=np.unique(kernel.sn.z,return_inverse=True)
    da=np.asarray(background.angular_diameter_distance(unique))[inverse]
    assert np.isfinite(da).all() and np.all(da>0)
    prediction=5*np.log10(da*(1+kernel.sn.z)*(1+kernel.sn.zhel))+25
    result['luminosity']=kernel.from_prediction(prediction,point)
    result['SN_loglike_absolute_error']=abs(result['luminosity']['SN_loglike']-native['native2_loglikes']['released_sn'])
    result['SN_chi2_absolute_error']=abs(result['luminosity']['SN_chi2']-native['native2_derived']['sn_chi2'])
    result['native_H_relative_error']=float(np.max(abs(background.hubble_parameter(rt.GRID)/np.array(native['value']['H'])-1)))
    for field,key in [('SN_loglike_absolute_error','SN_loglike_absolute_max'),('SN_chi2_absolute_error','SN_chi2_absolute_max'),
                      ('native_H_relative_error','background_native_H_relative_max')]:
        if not np.isfinite(result[field]) or result[field]>design[key]:result['failed_numerical_gates'].append(field)
    rt.consumer.finite_tree(result)
    return result


@contextmanager
def background_only():
    import camb
    import cobaya.model
    def forbidden(*a,**kw):raise RuntimeError('Accuracy2 histories permit background calculations only.')
    from contextlib import ExitStack
    with ExitStack() as stack:
        for owner,name in [(camb,'get_results'),(camb,'get_transfer_functions'),(cobaya.model,'get_model'),
                           (cobaya.model.Model,'__init__'),(camb.CAMBdata,'calc_power_spectra')]:
            stack.enter_context(patch.object(owner,name,forbidden))
        yield


def calculate(point,native,kernel,design,expansion_design,counter):
    import camb
    extra=background_options(native,kernel.declared_extra)
    cosmology={k:point[k] for k in ['H0','ombh2','omch2','ns','tau']}
    cosmology.update(As=1e-10*np.exp(point['logA']),w=point.get('w',-1.),wa=point.get('wa',0.))
    with background_only():
        parameters=camb.set_params(**cosmology,**extra)
        adjusted=expansion.native_thermal_parameters(parameters)
        counter['background_invocations']+=1
        background=camb.get_background(adjusted)
        result=consume_background(background,point,native,kernel,design,expansion_design)
    result['background_extra_args']=extra
    result['thermal_adapter']={'input_WantTransfer':bool(parameters.WantTransfer),'background_WantTransfer':bool(adjusted.WantTransfer)}
    return result


def bindings(target,design,expansion_design):
    configuration_checks(target)
    kernel,path=load_sn(target,design)
    value={'schema':'typed-native2-history-lineage-v1','numerical_target_identity':target.numerical_target_identity,
       'parent_work':rt.relative(target.plan_path.parent),'parent_summary':rt.relative(target.summary_path),
       'parent_input_sha256':target.input_sha256,'source_sha256':sources(),'validation_sha256':validation_guard(),
       'design':design,'expansion_design':expansion_design,'native_configuration':rt.canonical(target.native_configuration),
       'parent_settings':target.parent_settings,'scientific_versions':target.plan['evidence']['versions'],
       'sample_path':rt.relative(path),'sample_sha256':rt.digest(path)}
    return value,kernel


def prepare(work,summary,cache):
    target=typed_parent(work,summary)
    design=json.loads(DESIGN.read_text());old=json.loads(expansion.DESIGN.read_text())
    value,_=bindings(target,design,old);value['identity']=rt.identity(value)
    cache=Path(cache).resolve();assert cache.is_relative_to(ROOT/'.work') and not cache.exists()
    cache.mkdir(parents=True);rt.write_new(cache/'lineage.json',value)
    return value


def verify(cache):
    lineage=rt.read(Path(cache)/'lineage.json')
    target=typed_parent(ROOT/lineage['parent_work'],ROOT/lineage['parent_summary'])
    value,kernel=bindings(target,json.loads(DESIGN.read_text()),json.loads(expansion.DESIGN.read_text()))
    assert rt.plain(lineage)==dict(value,identity=rt.identity(value))
    return lineage,target,kernel


def point_binding(lineage,index,native):
    return {'identity':lineage['identity'],'index':index,'native2_record_path':native['record_path'],
            'native2_record_sha256':native['record_sha256'],'point':native['point']}


def read_records(cache,lineage,target):
    cache=Path(cache);manifest=cache/'background-records.json'
    if manifest.exists():rt.verify_hashes(rt.read(manifest)['files'])
    rows=[];files={}
    for i,native in enumerate(target.records):
        path=cache/f'{i:05d}.json';claim=cache/f'{i:05d}.claim.json';binding=point_binding(lineage,i,native)
        if path.exists():
            row=rt.read(path);assert row['binding']==binding and rt.read(claim)['binding']==binding
            assert row['claim_sha256']==rt.digest(claim)
            assert row['background_invocations'] in (0,1)
            if row['status']=='finite_native_accuracy2_history':assert row['background_invocations']==1 and not row['failed_numerical_gates']
        else:
            row={'binding':binding,'status':'incomplete_attempt_no_retry','unknown_background_count':claim.exists(),'background_invocations':0}
        for p in [path,claim]:
            if p.exists():files[rt.relative(p)]=rt.digest(p)
        rows.append(row)
    if (cache/'execution.json').exists():files[rt.relative(cache/'execution.json')]=rt.digest(cache/'execution.json')
    return rows,files


def summaries(rows,target,kernel,design,expansion_design):
    failures=[{'index':i,'status':r['status'],'gates':r.get('failed_numerical_gates'),'error':r.get('error')}
              for i,r in enumerate(rows) if r['status']!='finite_native_accuracy2_history']
    result={'status':'failed_native_accuracy2_history' if failures else 'qualified_native_accuracy2_histories',
            'failures':failures,'points':len(rows),'numerical_target_identity':target.numerical_target_identity,
            'parent_proposal_target_identity':target.parent_proposal_target_identity,'settings':target.parent_settings,
            'expansion':None,'luminosity':None,'CMB_spectrum_calls':0,
            'known_background_invocations':sum(r['background_invocations'] for r in rows),
            'unknown_count_attempts':sum(bool(r.get('unknown_background_count')) for r in rows)}
    if failures:return result
    weights=target.normalized_weights;lw=target.logweights;groups=target.groups
    arrays={k:np.array([r['history'][k] for r in rows]) for k in ['H_km_s_Mpc','q','j']}
    scalars={k:np.array([r['scalar'][k] for r in rows]) for k in rows[0]['scalar']}
    means=np.array([r['luminosity']['mean_mag'] for r in rows]);average=weights@means;deviations=means-average
    values=dict(scalars)
    for name,array in dict(arrays,B_mean_mag=means).items():values.update({f'{name}_grid{i}':x for i,x in enumerate(array.T)})
    gates=target.posterior_summary['importance_diagnostics']['weighted_stability_gates']
    stability=weight_diagnostics(lw,values,groups,gates)
    result['history_weighted_stability']=stability
    if not stability['qualified_overlap']:
        result['status']='failed_native_accuracy2_history_stability';return result
    result['expansion']={'redshift':expansion_design['redshift_grid'],
        'history':{k:expansion.pointwise_summary(x,weights,expansion_design['quantile_probabilities']) for k,x in arrays.items()},
        'scalar':{k:expansion.pointwise_summary(x,weights,expansion_design['quantile_probabilities']) for k,x in scalars.items()},
        'acceleration_fraction_by_redshift':[weighted_fraction(x<0,weights,groups) for x in arrays['q'].T],
        'any_acceleration_on_reported_past_grid':weighted_fraction(np.any(arrays['q']<0,axis=1),weights,groups)}
    q0=np.array([r['native2_derived']['q0'] for r in target.records]);qdev=q0-weights@q0
    result['luminosity']={'redshift':design['luminosity_redshifts'],
        'pointwise_magnitude_posterior':[luminosity.gaussian_mixture_summary(means[:,j],float(kernel.covariance[j,j]),weights) for j in range(means.shape[1])],
        'covariance_mag2':(kernel.covariance+(deviations*weights[:,None]).T@deviations).tolist(),
        'conditional_covariance_mag2':kernel.covariance.tolist(),
        'conditional_coefficient_covariance_mag2':None if kernel.coefficient_covariance is None else kernel.coefficient_covariance.tolist(),
        'covariance_with_native2_q0_mag':((deviations*weights[:,None]).T@qdev).tolist(),
        'B_at_zero_imposed':True,'baseline_zero_imposed':kernel.evolution=='none',
        'data_redshift_range':[float(kernel.sn.z.min()),float(kernel.sn.z.max())],
        'outside_observed_redshift_range':((kernel.z<kernel.sn.z.min())|(kernel.z>kernel.sn.z.max())).tolist()}
    result['numerical_check_maxima']={k:max(r['numerical_checks'][k] for r in rows) for k in rows[0]['numerical_checks'] if k!='age_CAMB_default_minus_integral_Gyr'}
    result['age_CAMB_default_minus_integral_Gyr_range']=[min(r['numerical_checks']['age_CAMB_default_minus_integral_Gyr'] for r in rows),max(r['numerical_checks']['age_CAMB_default_minus_integral_Gyr'] for r in rows)]
    result['SN_loglike_max_absolute_error']=max(r['SN_loglike_absolute_error'] for r in rows)
    result['SN_chi2_max_absolute_error']=max(r['SN_chi2_absolute_error'] for r in rows)
    result['native_H_max_relative_error']=max(r['native_H_relative_error'] for r in rows)
    return result


def execute(cache,output):
    cache=Path(cache);lineage,target,kernel=verify(cache)
    assert not (cache/'execution.json').exists() and not any(cache.glob('0*.json'))
    rt.write_new(cache/'execution.json',{'identity':lineage['identity'],'retry_permitted':False})
    for i,native in enumerate(target.records):
        path=cache/f'{i:05d}.json';claim=cache/f'{i:05d}.claim.json';binding=point_binding(lineage,i,native)
        rt.write_new(claim,{'binding':binding});counter={'background_invocations':0};started=time.monotonic()
        try:
            row=calculate(native['point'],native,kernel,lineage['design'],lineage['expansion_design'],counter)
            row['status']='failed_native_accuracy2_history_checks' if row['failed_numerical_gates'] else 'finite_native_accuracy2_history'
        except Exception as exc:row={'status':'failed_native_accuracy2_history_evaluation','error':repr(exc)}
        row.update(binding=binding,claim_sha256=rt.digest(claim),background_invocations=counter['background_invocations'],seconds=time.monotonic()-started)
        rt.write_new(path,row)
        if (i+1)%100==0:print(json.dumps({'background_slots':i+1}),flush=True)
    # Fresh immutable byte/version checks at completion; no second background.
    rt.verify_hashes(target.input_sha256);rt.verify_hashes(lineage['source_sha256'])
    rt.verify_hashes(lineage['validation_sha256'])
    assert rt.read(cache/'lineage.json')==lineage,'History lineage changed during execution.'
    assert all(importlib.metadata.version(k)==v for k,v in lineage['scientific_versions'].items())
    rows,files=read_records(cache,lineage,target)
    rt.write_new(cache/'background-records.json',{'identity':lineage['identity'],'files':files})
    result=summaries(rows,target,kernel,lineage['design'],lineage['expansion_design'])
    result.update(lineage_path=rt.relative(cache/'lineage.json'),lineage_sha256=rt.digest(cache/'lineage.json'),
                  record_manifest_path=rt.relative(cache/'background-records.json'),record_manifest_sha256=rt.digest(cache/'background-records.json'),
                  source_sha256=lineage['source_sha256'],parent_input_sha256=target.input_sha256,
                  interpretation=lineage['design']['qualification_limits'])
    rt.write_new(cache/'summary.json',result);rt.write_new(output,result)
    return result


def report(cache,output=None):
    cache=Path(cache);lineage,target,kernel=verify(cache);rows,files=read_records(cache,lineage,target)
    manifest=rt.read(cache/'background-records.json')
    assert manifest['identity']==lineage['identity'] and manifest['files']==files
    result=summaries(rows,target,kernel,lineage['design'],lineage['expansion_design'])
    saved=rt.plain(rt.read(cache/'summary.json'))
    for k,v in result.items():assert saved[k]==v
    assert saved['lineage_sha256']==rt.digest(cache/'lineage.json') and saved['record_manifest_sha256']==rt.digest(cache/'background-records.json')
    assert saved['source_sha256']==sources() and saved['parent_input_sha256']==target.input_sha256
    if output is not None:assert rt.plain(rt.read(output))==saved
    return saved


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','execute','report'])
    p.add_argument('--work',type=Path);p.add_argument('--summary',type=Path);p.add_argument('--cache',type=Path,required=True);p.add_argument('--output',type=Path)
    a=p.parse_args()
    if a.action=='prepare':assert a.work and a.summary;result=prepare(a.work,a.summary,a.cache)
    elif a.action=='execute':assert a.output;result=execute(a.cache,a.output)
    else:result=report(a.cache,a.output)
    print(json.dumps({'status':result.get('status','prepared_without_physical_calls'),'identity':result.get('identity')}))
