"""Synthetic persistent-pool/cache/typed-qualifier validation; no physics calls.

The original observational parent and numerical prerequisite consumers are
explicitly mocked in the integration fixture. Their production paths are never
relaxed. Files, seals, slot replay, accounting and final new qualifier are real.
"""
import argparse
import copy
from contextlib import ExitStack
import json
from pathlib import Path
from types import SimpleNamespace
import tempfile
from unittest.mock import patch
import numpy as np
from scipy.stats import norm
import native_accuracy_runtime as rt
import native_accuracy_pilot as pilot
import native_accuracy_execute as execute
import native_accuracy_qualify as qualifier
from exact_correction import record_digest


def refused(fn):
    try:fn()
    except (AssertionError,ValueError,FileExistsError,RuntimeError,KeyError):return True
    raise AssertionError('Invalid synthetic state was admitted.')


def spectra(point):
    ell=np.arange(20);x=point['x']
    return {'ell':ell,**{k:np.r_[0.,(np.arange(1,20)+1.)**-2*(1+.001*x)*(i+1)]
                         for i,k in enumerate(rt.contracts.SPECTRA)}}


EXTRA={'AccuracyBoost':2,'lAccuracyBoost':2,'lSampleBoost':2,'lens_potential_accuracy':4,'lmax':19}
META={'finalized_theory_extra_args':EXTRA,'CAMB_Params_max_l':19,'provider_length':20,
      'units':'FIRASmuK2','ell_factor':False,'pp':'dimensionless C_phi_phi'}


def value_at(point):
    x=point['x'];like=-x*x/2-.02*x
    return {'loglikes':{'synthetic':like},'logpriors':[0.],'logpost':like,
            'derived':{'test_derived':2*x+3},'H':(70*(1+.001*x)*np.sqrt(1+rt.GRID)).tolist(),
            'metadata':copy.deepcopy(META)}


class FakeModel:
    """An analytic scalar Gaussian stand-in, never a cosmological model."""
    def __init__(self,counter,fail_at=None):
        self.counter=counter;self.fail_at=fail_at;self.local_calls=0;self.current=None
        counter['constructions']+=1
        self.likelihood={'synthetic':None};self.theory={'camb':SimpleNamespace(extra_args=EXTRA.copy()),'camb.transfers':None}
        self.parameterization=SimpleNamespace(derived_params=lambda:['test_derived'])
        self.provider=SimpleNamespace(get_CAMBdata=self.background,get_Cl=self.get_cl)
    def logposterior(self,point):
        self.counter['invocations']+=1;self.local_calls+=1
        if self.local_calls==self.fail_at:raise ValueError('Injected synthetic density failure.')
        self.current=point;v=value_at(point)
        return SimpleNamespace(logpost=v['logpost'],loglikes=list(v['loglikes'].values()),logpriors=v['logpriors'],derived=list(v['derived'].values()))
    def background(self):
        return SimpleNamespace(Params=SimpleNamespace(max_l=19),hubble_parameter=lambda z:70*(1+.001*self.current['x'])*np.sqrt(1+z))
    def get_cl(self,ell_factor,units):
        assert units=='FIRASmuK2'
        raw=spectra(self.current)
        return rt.contracts.dl_from_raw(raw) if ell_factor else raw
    def close(self):pass


def fixture(root):
    source=root/'source';source.mkdir();native=source/'original';native.mkdir()
    rng=np.random.default_rng(274021);x=norm.ppf((np.arange(500)+.5)/500);rng.shuffle(x)
    allx=np.tile(x,4);allx[167]=allx[166];allx[1229]=allx[1228]
    points=[{'x':float(v)} for v in allx]
    locations=[{'chain':f'chain.{i//500+1}.txt','row':i%500,'expanded_index':3*(i%500)} for i in range(2000)]
    locations[167]=copy.deepcopy(locations[166]);locations[1229]=copy.deepcopy(locations[1228])
    hashes={}
    for i,point in enumerate(points):
        q=-point['x']**2/2
        row={'status':'finite','index':i,'point':point,'exact_logpost':q,'proposal_logpost':q,
             'log_weight':0.,'exact_loglikes':{'synthetic':q},'proposal_loglikes':{'synthetic':q},
             'derived':{'test_derived':2*point['x']-7},'target_identity':'synthetic_original_target'}
        row['payload_sha256']=record_digest(row);path=native/f'{i:05d}.json';path.write_text(json.dumps(row))
        hashes[rt.relative(path)]=rt.digest(path)
    selection={'points':points,'groups':np.repeat(np.arange(4),500).tolist(),'locations':locations}
    selection_path=native/'selection.json';selection_path.write_text(json.dumps(selection))
    summary_path=source/'original-summary.json'
    summary_path.write_text(json.dumps({'selection_path':rt.relative(selection_path),'native_record_sha256':hashes}))
    original_plan=source/'original-screen-plan.json';rt.write_new(original_plan,{'correction_summary_path':rt.relative(summary_path)})
    refs=[]
    for audit,index in enumerate(np.linspace(0,1999,32,dtype=int)):
        index=int(index);point=points[index];v=value_at(point);record=source/f'old-high2-{audit:02}.json';spec=source/f'old-high2-{audit:02}.npz'
        rt.save_spectrum(spec,rt.contracts.dl_from_raw(spectra(point)))
        row={'status':'finite_native_precision','failed_checks':[],'point':point,'audit_index':audit,
             'high_loglikes':v['loglikes'],'high_logpriors':v['logpriors'],'high_logpost':v['logpost'],'high_derived':v['derived'],
             'background':{'high_H_km_s_Mpc':v['H']},'finalized_theory_extra_args':EXTRA,'CAMB_Params_max_l':19,'provider_Dl_length':20}
        rt.write_new(record,row)
        refs.append({'audit_index':audit,'parent_index':index,'point':point,'native_record_path':rt.relative(native/f'{index:05d}.json'),
                     'native_record_sha256':rt.digest(native/f'{index:05d}.json'),'original_location':locations[index],
                     'record_path':rt.relative(record),'record_sha256':rt.digest(record),
                     'spectrum_path':rt.relative(spec),'spectrum_sha256':rt.digest(spec),
                     'original_status':row['status'],'original_failed_checks':[]})
    screen=source/'screen.json';review=source/'review.json'
    screen.write_text('{"synthetic_boundary":true}');review.write_text('{"synthetic_boundary":true}')
    validation=source/'validation.json';validation.write_text('{"synthetic_validation_boundary":true}')
    info={'theory':{'camb':{'extra_args':EXTRA.copy()}},'likelihood':{'synthetic':{}},'params':{'x':{'prior':{'min':-10,'max':10}}}}
    evidence={'receipt_identity':'explicitly_mocked_parent_receipt','parent_proposal_identity':'synthetic_original_target',
              'screen_path':rt.relative(screen),'review_path':rt.relative(review),'original_plan_path':rt.relative(original_plan),
              'settings':{'synthetic':True},'references':refs,'native_configuration':rt.canonical(info),'versions':{},
              'input_sha256':rt.merge(hashes,{rt.relative(p):rt.digest(p) for p in source.iterdir() if p.is_file()},
                                    {rt.relative(selection_path):rt.digest(selection_path)}),
              'original_screen_status':'synthetic','thermal_reviewed_status':'numerical_sensitivity_requires_followup'}
    return evidence,info,screen,review,validation


def run_worker_fixture(work,plan,verify,builder,counter,fail=None):
    rt.write_new(work/'execution.json',{'plan_identity':plan['identity'],'coordinator_pid':9000,'workers':plan['workers'],'retry_permitted':False})
    for request in plan['requests']:
        if request['origin']=='persistent_native':rt.write_new(rt.request_paths(work,request['index'])['ticket'],{'binding':rt.request_identity(plan,request)})
    for worker in range(plan['workers']):
        pid=10000+worker;log=work/f'worker-{worker}.log';log.write_text('Synthetic analytic model; no external process or physics.\n')
        with patch.object(rt,'construct_model',lambda info:FakeModel(counter,fail)),patch.object(rt.os,'getpid',lambda:pid):
            ok=rt.persistent_worker(work,worker,verify,builder)
        rt.write_new(work/f'worker-{worker}-process.json',{'plan_identity':plan['identity'],'worker':worker,'pid':pid,
            'returncode':0 if ok else 1,'log_path':rt.relative(log),'log_sha256':rt.digest(log)})


def main_run():
    base=rt.ROOT/'.work/unified-cosmology/native-accuracy-execution-validation';base.mkdir(parents=True,exist_ok=True)
    root=Path(tempfile.mkdtemp(prefix='case-',dir=base))
    physical_calls=[]
    import camb
    import cobaya.model
    def forbid(*a,**kw):physical_calls.append('forbidden');raise RuntimeError('Physical call prohibited in validator.')
    report={'scope':'Synthetic analytic model and explicitly mocked observational prerequisite boundary; no observational numerical qualification.'}
    with ExitStack() as guards:
        for owner,name in [(camb,'get_results'),(camb,'get_background'),(camb,'get_transfer_functions'),(cobaya.model,'get_model'),(cobaya.model.Model,'__init__')]:
            guards.enter_context(patch.object(owner,name,forbid))
        evidence,info,screen,review,validation=fixture(root)
        guards.enter_context(patch.object(rt,'evidence',lambda s,r:(copy.deepcopy(evidence),copy.deepcopy(info))))
        guards.enter_context(patch.object(rt,'validation_guard',lambda:{rt.relative(validation):rt.digest(validation)}))
        for name,value in json.loads(rt.DESIGN.read_text())['environment'].items():guards.enter_context(patch.dict(rt.os.environ,{name:value}))
        work=root/'pilot';plan=pilot.prepare(screen,review,work)
        assert [[r['audit_index'] for r in plan['requests'] if r['worker']==w] for w in range(2)]==[[0,31,0,31],[31,0,31,0]]
        counter={'constructions':0,'invocations':0};run_worker_fixture(work,plan,pilot.verify_plan,pilot.build_record,counter)
        result,files=pilot.replay(work);rt.seal_summary(work,plan,result,files)
        assert result['status']=='passed_persistent_accuracy2_replay_pilot' and counter=={'constructions':2,'invocations':8}
        checked=pilot.verify_completed(work)
        report['persistent_pilot']={'constructs':2,'density_requests':8,'replayed_records':len(checked['summary']['requests'])}
        assert refused(lambda:rt.execute_workers(work,'native_accuracy_pilot.py'))
        # Replayed source values, multipoles and each numerical gate fail independently.
        value,raw=rt.old_reference(evidence['references'][0]);mutation_checks=[]
        for key in ['loglikes','logpriors','logpost','derived','H','metadata','spectrum','ell']:
            changed=copy.deepcopy(value);changed_raw={k:v.copy() for k,v in raw.items()}
            if key=='loglikes':changed[key]['synthetic']+=1e-3
            elif key=='logpriors':changed[key][0]+=1e-3
            elif key=='logpost':changed[key]+=1e-3
            elif key=='derived':changed[key]['test_derived']+=1e-3
            elif key=='H':changed[key][0]*=1.001
            elif key=='metadata':changed[key]['CAMB_Params_max_l']+=1
            elif key=='spectrum':changed_raw['tt'][2]*=1.001
            else:changed_raw['ell'][2]+=1
            def check():
                result=rt.compare(changed,changed_raw,value,raw,plan['design']['pilot_tolerances']);assert not result['failed_checks']
            assert refused(check);mutation_checks.append(key)
        report['replay_mutations_rejected']=mutation_checks
        # Claim without result is never a known zero or a permitted rerun.
        interrupted=root/'interrupted';ip=pilot.prepare(screen,review,interrupted)
        paths=rt.request_paths(interrupted,0);rt.write_new(paths['ticket'],{'binding':rt.request_identity(ip,ip['requests'][0])})
        rt.write_new(paths['claim'],{'binding':rt.request_identity(ip,ip['requests'][0]),'pid':1})
        assert rt.incomplete(interrupted,ip,ip['requests'][0])['unknown_native_invocation_count']
        assert refused(lambda:rt.execute_workers(interrupted,'native_accuracy_pilot.py'))
        failure=root/'failure';fp=pilot.prepare(screen,review,failure);fc={'constructions':0,'invocations':0}
        run_worker_fixture(failure,fp,pilot.verify_plan,pilot.build_record,fc,fail=2)
        fr,ff=pilot.replay(failure);rt.seal_summary(failure,fp,fr,ff)
        assert fc=={'constructions':2,'invocations':4} and fr['status']=='failed_or_incomplete_persistent_accuracy2_pilot'
        assert sum(r.get('known_unstarted',False) for r in fr['requests'])==4
        assert refused(lambda:pilot.verify_completed(failure))
        report['failure_accounting']={'injected_failures':2,'known_requests':4,'unstarted':4,'interrupted_claim_unknown':True,'retry_refused':True}
        # Full prerequisite boundary is mocked, but the successful real pilot is bound.
        refinement=root/'synthetic-refinement';refinement.mkdir();(refinement/'evidence.json').write_text('{"synthetic":true}')
        def prerequisite_fixture(ev,r,p):
            assert ev==evidence and Path(p)==work and Path(r)==refinement
            prior=pilot.verify_completed(p)
            return {'refinement_work':rt.relative(r),'pilot_work':rt.relative(p),'refinement':{'explicitly_mocked':True},
                    'pilot_plan_identity':prior['plan']['identity'],'pilot_summary_sha256':rt.digest(Path(p)/'summary.json'),
                    'input_sha256':rt.merge(prior['input_sha256'],{rt.relative(refinement/'evidence.json'):rt.digest(refinement/'evidence.json')})}
        guards.enter_context(patch.object(execute,'prerequisites',prerequisite_fixture))
        full=root/'full';full_plan=execute.prepare(screen,review,refinement,work,full,4)
        execute.reuse_fixed32(full,full_plan)
        count={'constructions':0,'invocations':0};run_worker_fixture(full,full_plan,execute.verify_plan,execute.build_record,count)
        result,files,rows=execute.replay(full);rt.seal_summary(full,full_plan,result,files)
        assert count=={'constructions':4,'invocations':1968}
        assert result['status']=='qualified_native_accuracy2_importance_posterior'
        typed=qualifier.qualify_run(full,full/'summary.json')
        assert type(typed) is qualifier.QualifiedNativeAccuracy2 and typed.numerical_accuracy==2
        assert typed.native_configuration==info and typed.numerical_target_identity!=typed.parent_proposal_target_identity
        assert len(typed.records)==len(typed.locations)==2000 and result['reused_fixed32']==32
        derived_error=max(abs(r['native2_derived']['test_derived']-r['native1_derived']['test_derived']-10) for r in rows)
        assert derived_error<2e-14
        assert typed.locations[166]==typed.locations[167] and typed.selected_points[166]==typed.selected_points[167]
        assert typed.locations[1228]==typed.locations[1229] and typed.selected_points[1228]==typed.selected_points[1229]
        x=np.array([p['x'] for p in typed.selected_points]);w=typed.normalized_weights
        mean=float(w@x);variance=float(w@(x-mean)**2)
        assert abs(mean+.02)<.01 and abs(variance-1)<.02
        assert np.max(abs(typed.logweights+.02*x))<1e-14
        report['full_synthetic_cohort']={'slots':2000,'reused':32,'new_requests':1968,'persistent_constructions':4,
           'preserved_repeated_slots':2,'weighted_x_mean':mean,'analytic_gaussian_x_mean':-.02,
           'weighted_x_variance':variance,'analytic_gaussian_x_variance':1.,
           'raw_weight_ess':result['importance_diagnostics']['raw_weight_ess'],
           'fresh_native2_derived_verified':True,'fresh_derived_arithmetic_max_error':derived_error,'typed_identity_distinct':True}
        # Statistical gates stay active even when all densities are finite.
        concentration=copy.deepcopy(rows)
        for i,row in enumerate(concentration):row['log_weight']=0. if i<500 else -1000.
        failed=execute.summarize_rows(concentration,typed.groups)
        assert failed['status']=='failed_native_accuracy2_importance_gates' and failed['posterior'] is None
        assert failed['importance_diagnostics']['failed_gates']
        assert execute.summarize_rows(rows[:1999],typed.groups[:1999])['posterior'] is None
        report['statistical_refusals']=['chain_weight_concentration','fewer_than_2000']
        # A resealed record cannot bypass the previously frozen file ledger.
        target=full/'00000.json';before=target.read_bytes();data=rt.read(target);data['native2_logpost']+=.1
        data['payload_sha256']=rt.identity(data);target.write_text(json.dumps(data))
        assert refused(lambda:qualifier.qualify_run(full,full/'summary.json'));target.write_bytes(before)
        spectrum_path=rt.ROOT/rows[0]['spectrum_path'];before=spectrum_path.read_bytes();spectrum_path.write_bytes(before+b'changed')
        assert refused(lambda:qualifier.qualify_run(full,full/'summary.json'));spectrum_path.write_bytes(before)
        summary=full/'summary.json';before=summary.read_bytes();data=rt.read(summary);data['numerical_target_identity']='false'
        data['payload_sha256']=rt.identity(data);summary.write_text(json.dumps(data))
        assert refused(lambda:qualifier.qualify_run(full,summary));summary.write_bytes(before)
        src=validation;before=src.read_bytes();src.write_text('{"changed":true}')
        assert refused(lambda:qualifier.qualify_run(full,summary));src.write_bytes(before)
        qualifier.qualify_run(full,summary)
        report['tamper_refusals']=['resealed_record','spectrum_bytes','resealed_summary_identity','bound_validation_bytes']
    assert physical_calls==[]
    report.update(status='passed_synthetic_native_accuracy_execution_validation',source_sha256=rt.sources(),
                  physical_calls=0,synthetic_work=rt.relative(root),
                  boundary_limits='Observational parent/32-screen/refinement receipts mocked explicitly; this validates execution/typed consumer mechanics and Gaussian algebra, not real numerical adequacy or actual OS/GPU efficiency.')
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=rt.VALIDATION);a=p.parse_args()
    result=main_run();a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':result['status'],'physical_calls':result['physical_calls'],'full':result['full_synthetic_cohort']}))
