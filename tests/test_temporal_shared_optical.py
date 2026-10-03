"""Hand-written transport fixtures; these tests do not qualify the optical physics.

The parent runs this source under its compute lease. No old generated receipt is
needed. Process tests execute tiny Python helpers only, never compiler/native/reference.
"""
import copy
from decimal import Decimal as D
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import stat
import sys
import tarfile
import tempfile
import types
import unittest
from unittest import mock

PACKET=Path(__file__).resolve().parent
if not (PACKET/'controller.py').exists():
    PACKET=Path(__file__).resolve().parents[1]/'experiments'/'temporal-shared-optical-control'
def load_source(name,path):
    spec=importlib.util.spec_from_file_location(name,path); module=importlib.util.module_from_spec(spec)
    sys.modules[name]=module; spec.loader.exec_module(module); return module
V=load_source('temporal_optical_verifier_under_test',PACKET/'verify_sdk.py')
prior_verifier=sys.modules.get('verify_sdk')
sys.modules['verify_sdk']=V
try: C=load_source('temporal_optical_controller_under_test',PACKET/'controller.py')
finally:
    if prior_verifier is None: del sys.modules['verify_sdk']
    else: sys.modules['verify_sdk']=prior_verifier


def request_fixture():
    return V.parse(V.REQUEST_PROTOTYPE_TEXT.encode())


def native_fixture():
    """Synthetic closed emitter-shaped data; values are transport examples only."""
    def outcome(value=None,availability=0,status='ok'):
        return {'availability':availability,'status':status,'value':value}
    def temporal(control=None):
        rows=[]
        for i in range(8):
            e=(i%4)//2; lo,hi=((9,11),(12,15))[e]; photons=1; admission='ok'; coverage='full' if e==0 else 'partial'; duration=hi-lo; covered=2; fraction=D(1) if e==0 else D(2)/3
            band=2*(i//4)+i%2
            if control=='no-time-overlap' and i==0: lo,hi=20,22; photons=0; covered=0; fraction=0; coverage='no_overlap'
            if control=='zero-duration' and i==0: hi=lo; duration=0; admission='outside_domain'; coverage='unassessed'; covered=0; fraction=0
            if control=='invalid-redshift' and i==0: admission='outside_domain'; coverage='unassessed'; covered=0; fraction=0
            if control=='required-state-row-refusal' and i==7: admission='invalid_input'; band=4; coverage='unassessed'; covered=0; fraction=0
            if control=='temporal-work-refusal': admission='work_limit'; coverage='unassessed'; covered=0; fraction=0
            rows.append({'index':i,'grid_index':0,'band_index':band,'source_epoch_second':10,'observer_lower_second':lo,'observer_upper_second':hi,
                         'admission_status':admission,'observer_duration_second':duration,'covered_observer_second':covered,'covered_fraction':fraction,'coverage':coverage,'segment_work':0,
                         'photons':outcome(photons,1) if admission=='ok' else outcome(None,2,admission),'mean_flux':outcome(),'energy':outcome()})
        if control=='temporal-row-admission-refusal': rows=[]
        okay=control in (None,'no-time-overlap')
        return {'status':'work_limit' if control=='temporal-row-admission-refusal' else 'ok','segment_work':0,'required_rows_admitted':okay,'rows':rows}
    def batch(channel=0,noise=0,photons=1,single=False,badqe=False,failed=None,zero=False):
        rows=[]
        for i in range(1 if single else 6):
            detected=(i<3 and not(i==2 and channel==3)) if not single else not zero
            status=failed or 'ok'
            rows.append({'status':status,'detected':detected,'electron_count':(2,3,3,2)[channel] if detected and noise==0 else None,
                         'measured_adu':D(3 if channel in (0,3) else '3.5') if detected and noise else None,
                         'threshold_adu':D('1.5'),'zero_probability':zero,'log_value':None if failed or zero else D(-1),'log_error':0,'detection_probability':None if failed else D('0.5')})
        return {'status':'outside_domain' if badqe else 'ok','photons':photons,'QE':D('1.010000000000000009') if badqe else D('0.5' if channel%2 else '0.75'),
                'full_exposure_second':2 if channel<2 else 3,'sigma_electrons':0 if noise==0 else D('0.75'),'poisson_terms':0,'omitted_tail':0,'rows':[] if badqe else rows}
    def reducer(index=0,status='ok'):
        item={'record_index':index,'selected_only':index==1,'status':status,'conditional':[{'state_id':'S'+str(i),'mass':D('0.25' if i==0 else '0.75'),'status':status,'log_record':-1,'record_error':0,'log_event':-1,'event_error':0} for i in range(2)],
              'joint_log_error':0,'event_log_error':0,'value_log_error':0}
        for key in ('log_joint','log_event','log_value','product_per_exposure_log_record','product_per_exposure_log_event','record_relative_difference','event_relative_difference'):
            item[key]=None if status!='ok' else (D('0.02') if 'relative_difference' in key else (0 if key=='log_value' and index==1 else -1))
        return item
    groups=[]
    for noise in range(2):
        groups.append({'sigma_electrons':0 if noise==0 else D('0.75'),
                       'attempts':[{'state_index':i//4,'channel_index':i%4,'lambda_electrons':1,'batches':[batch(i%4,noise,p) for p in (1,D('0.9999999999968'),D('1.0000000000032'))]} for i in range(8)],
                       'records':[reducer(i) for i in range(3)]})
    request=request_fixture(); axes=request['shape_contract']['axes']
    singles=[batch(noise=1,single=True,badqe=True),batch(noise=1,single=True,failed='conditioning_budget_exceeded'),batch(noise=1,single=True,failed='invalid_input'),batch(single=True,zero=True)]
    external=[]
    for i,status in enumerate(('outside_domain','conditioning_budget_exceeded','invalid_input')):
        external.append({'id':axes['reducer_control_ids'][i],'retained_state_masses':[D('0.25'),D('0.75')],
                         'replacement_batch':batch(3,1,badqe=i==0,failed=status) if i<2 else None,'reducer':reducer(0 if i<2 else 2,status)})
    return {'schema':'external-temporal-optical-sdk-control/v1','model_id':'finite_bilinear_rest_spectral_time_zero_outside_full_observer_exposure_mean',
            'detector_model_id':'DETECTOR/Poisson-arrivals-fixed-QE-Gaussian-read/v1','selection_id':'DETECTOR/threshold-joint-or-selected-censoring/v1',
            'sampled_empirical_sensitivity':D('5e-12'),'temporal_photon_endpoint_sensitivity':D('3.2e-12'),
            'prepared_status':'ok','grid_count':1,'band_count':4,'prepare_seconds':0,'temporal_evaluate_seconds':0,'total_native_seconds':0,
            'retained_temporal_payload_bytes':1,'detector_six_row_payload_bound_bytes':1,'primary':temporal(),'detector_controls':groups,
            'primary_poisson_terms':0,'all_poisson_terms':0,'all_temporal_segment_work':0,
            'actual_temporal_controls':[{'id':name,'native':temporal(name),'joint_withheld':name!='no-time-overlap'} for name in axes['temporal_control_ids']],
            'actual_detector_controls':[{'id':name,'native':singles[i]} for i,name in enumerate(axes['detector_control_ids'])],
            'actual_external_reducer_controls':external}


def reference_payload_fixture():
    request=request_fixture(); types=request['shape_contract']['reference_types']
    frequency={key:'1e-60' if 'error' in key else '1' for key in types['FrequencyRow']}
    discrete={'count_probability':'0.2','count_probability_error':'1e-60','below':'0.5','below_error':'1e-60'}
    continuous={key:'1e-60' if 'error' in key or 'tail' in key else '0.5' for key in types['ContinuousDetectorReferenceRow']}
    records=[]
    for i in range(6):
        row={key:'1e-60' if 'error' in key else ('0.02' if 'difference' in key else '-1') for key in types['ReferenceRecord']}
        row.update(sigma_index=i//3,record_index=i%3,conditional_record=['0.2','0.3'],conditional_event=['0.5','0.5'],positive_joint='0.25',positive_event='0.5')
        if i%3==1: row['log_value']='0'
        records.append(row)
    return {'schema':'original-high-precision-temporal-optical-reference/v1','mpmath_version':'1.3.0','decimal_digits':60,
            'shared_ancestry':['declared equations and exact SI constants','mpmath arithmetic/GL node generator; no native kernels'],
            'frequency_seconds':0,'elapsed_seconds':0,'integrand_evaluations':{'time_frequency':30720,'characteristic_function':73728},
            'photons':['1']*8,'photon_errors':['1e-60']*8,'frequency_rows':[copy.deepcopy(frequency) for _ in range(8)],'lambda_electrons':['1']*8,
            'detector_rows':[{'sigma_electrons':0,'rows':[copy.deepcopy(discrete) for _ in range(8)]},{'sigma_electrons':D('0.75'),'rows':[copy.deepcopy(continuous) for _ in range(8)]}], 'records':records}


def reference_fixture(payload=None,status='completed'):
    request=request_fixture(); payload=payload or reference_payload_fixture(); hashes={name:'a'*64 for name in ('consumer.cpp','reference.py','controller.py','verify_sdk.py','contract.snapshot.json')}
    runtime={'placeholder':'mocked getter; no numerical runtime test'}
    identities={'request_sha256':'b'*64,'original_contract_sha256':V.CONTRACT_SHA,'admission_sha256':request['origin']['admission_sha256'],'original_reference_sha256':request['origin']['original_sources'][2]['sha256'],'reference_script_sha256':hashes['reference.py'],'current_source_sha256':hashes}
    partial={'completed_groups':['frequency','photons','lambda','detector-sigma0','detector-sigma075','records'],'photon_prefix':payload['photons'],'photon_error_prefix':payload['photon_errors'],
             'frequency_prefix':payload['frequency_rows'],'lambda_prefix':payload['lambda_electrons'],'detector_prefix':payload['detector_rows'],'active_sigma_index':None,'active_sigma_rows':[],'record_prefix':payload['records'],'control_failure':None}
    error=None
    if status=='refused':
        message='late source identity failure'; error={'kind':'admission','stage':'source-terminal-drift','message':message,'message_truncated':False,'message_sha256':V.digest(message.encode()),'details':None}
    return {'schema_version':1,'interface_id':'temporal-shared-optical-reference-transport/v1','status':status,'error':error,'identities':identities,
            'axes':{key:value for key,value in request['shape_contract']['axes'].items() if key!='unit_policy'},
            'settings':{'scientific':request['shape_contract']['unchanged_scientific_settings'],'resources':request['resources']},
            'runtime_before':runtime,'runtime_after':copy.deepcopy(runtime),'payload':payload,'partial':partial,
            'work':{'elapsed_seconds':0,'time_frequency_evaluations':30720,'characteristic_function_evaluations':73728,'cpu_user_seconds':0,'cpu_system_seconds':0,'maximum_rss_kib':0},
            'qualification':{'scope':'empirical-synthetic-fixed-temporal-shared-optical-control/v1','native_outputs_consumed':False,'physical_qualification':False,'observational_qualification':False,'original_input_certificate':False,'limits':['Synthetic fixture only.']}},hashes,runtime


def numeric_wire(value):
    """Fixture JSON tokens retain Decimal numbers, without binary64 conversion."""
    if type(value) is D: return str(value)
    if type(value) is dict: return '{'+','.join(json.dumps(key)+':'+numeric_wire(item) for key,item in value.items())+'}'
    if type(value) is list: return '['+','.join(numeric_wire(item) for item in value)+']'
    return json.dumps(value,allow_nan=False)


def reference_refusal_prefix(detector_count,records,active=None,active_rows=0):
    value,hashes,runtime=reference_fixture(status='refused'); payload=value['payload']; value['payload']=None
    partial=value['partial']; partial['completed_groups']=partial['completed_groups'][:3+detector_count]
    partial['detector_prefix']=copy.deepcopy(payload['detector_rows'][:detector_count])
    partial['record_prefix']=copy.deepcopy(payload['records'][:records]); partial['active_sigma_index']=active
    partial['active_sigma_rows']=copy.deepcopy(payload['detector_rows'][active]['rows'][:active_rows]) if active is not None else []
    return value,hashes,runtime


class TypedAdmissionTests(unittest.TestCase):
    def test_failed_conditional_state_cannot_earn_successful_reducer(self):
        n=native_fixture(); record=n['detector_controls'][0]['records'][0]
        record['conditional'][0]['status']='outside_domain'
        with self.assertRaises(ValueError): V.native_admission(n,request_fixture())
        record['status']='outside_domain'
        for key in ('log_joint','log_event','log_value','product_per_exposure_log_record','product_per_exposure_log_event','record_relative_difference','event_relative_difference'): record[key]=None
        self.assertFalse(V.native_admission(n,request_fixture()))
        self.assertEqual(record['conditional'][0]['log_record'],-1)
        self.assertEqual(record['conditional'][1]['status'],'ok')

    def test_reference_interleaved_earned_partial_phases(self):
        cases=[reference_refusal_prefix(0,0,0,2),reference_refusal_prefix(1,3),reference_refusal_prefix(1,3,1,4)]
        for value,hashes,runtime in cases:
            with self.subTest(partial=value['partial']): self.assertFalse(self.admit_reference(value,hashes,runtime))
        value,hashes,runtime=cases[2]
        value['partial']['control_failure']={'stage':'conditional-positive-reference-gate','sigma_index':1,'channel_index':4,
            'lambda':'1','density':'0','density_error':'1','nondetection':'0.5','nondetection_error':'1e-60',
            'coarse_refined_attempt':reference_payload_fixture()['detector_rows'][1]['rows'][4]}
        self.assertFalse(self.admit_reference(value,hashes,runtime))

    def test_reference_completion_rejects_active_failure_and_false_prefixes(self):
        value,hashes,runtime=reference_fixture()
        failure={'stage':'positive-mixture-error-gate','value':'0.1','error':'1','conditional_values':['0.1','0.1'],'conditional_errors':['1','1']}
        mutations=[lambda r:r['partial'].update(control_failure=failure),lambda r:r['partial'].update(active_sigma_index=1,active_sigma_rows=[]),
                   lambda r:r['partial']['frequency_prefix'].pop(),lambda r:r['partial']['photon_error_prefix'].pop()]
        for mutate in mutations:
            candidate=copy.deepcopy(value); mutate(candidate)
            with self.subTest(mutate=mutate),self.assertRaises(ValueError): self.admit_reference(candidate,hashes,runtime)
        value,hashes,runtime=reference_refusal_prefix(1,3,1,2)
        for mutate in (lambda r:r['partial'].update(active_sigma_index=0),lambda r:r['partial']['record_prefix'].pop(),
                       lambda r:r['partial']['frequency_prefix'].pop(),lambda r:r['partial'].update(detector_prefix=[])):
            candidate=copy.deepcopy(value); mutate(candidate)
            with self.subTest(mutate=mutate),self.assertRaises(ValueError): self.admit_reference(candidate,hashes,runtime)
    def test_exact_longdouble_and_small_refinement_survive(self):
        value=V.parse(b'{"native":0.100000000000000005551,"reference":"0.100000000000000005552"}')
        self.assertEqual(D(value['reference'])-value['native'],D('1e-21'))
        self.assertEqual(V.decimal_string('1.00000000000000000000000000000000000000000000000000000000001')-D(1),D('1e-59'))
    def test_duplicate_nonfinite_bool_and_decimal_grammar(self):
        for blob in (b'{"x":1,"x":2}',b'{"x":NaN}',b'{"x":Infinity}',b'{"x":1e1001}'):
            with self.subTest(blob=blob),self.assertRaises(ValueError): V.parse(blob)
        for value in (True,1.0,D('NaN'),D('Infinity')):
            with self.subTest(value=value),self.assertRaises(ValueError): V.number(value)
        for text in ('+1','01','1.','NaN','Infinity',' 1','1e1001'):
            with self.subTest(text=text),self.assertRaises(ValueError): V.decimal_string(text)
        with self.assertRaises(ValueError): V.exact({'row':True},{'row':1})
    def test_bounded_depth_nodes_and_bytes(self):
        for blob in (b'['*34+b'0'+b']'*34,b'['+b'0,'*65536+b'0]',b'"'+b'x'*(V.MIB)+b'"'):
            with self.assertRaises(ValueError): V.parse(blob)
    def test_frozen_request_types_policies_unknown_keys(self):
        request=request_fixture()
        request.update(request_id='temporal-shared-optical-reproducible-request/v1',interface_id='temporal-shared-optical-bounded-controller/v1',status='reviewed-native-synthetic-control',execution=None)
        for name in ('reference','controller','sdk_verifier'): request['source_ports'][name]['sha256']='a'*64
        def identity(path,limit=V.MIB,expected=None): return b'',{'sha256':expected}
        with mock.patch.object(V,'file_identity',side_effect=identity):
            self.assertIs(V.request_admission(request,Path('/synthetic')),request)
            mutations=[lambda r:r['resources'].update(threads=True),lambda r:r['resources'].update(native_cpu_seconds=181),
                       lambda r:r.update(command='injected'),lambda r:r['shape_contract']['axes'].update(extra_axis=True)]
            for mutate in mutations:
                changed=copy.deepcopy(request); mutate(changed)
                with self.subTest(mutate=mutate),self.assertRaises(ValueError): V.request_admission(changed,Path('/synthetic'))
    def test_native_fixture_and_changed_axes_status_work_domains(self):
        baseline=native_fixture(); request=request_fixture(); self.assertTrue(V.native_admission(baseline,request))
        mutations=[lambda n:n['primary']['rows'].reverse(),lambda n:n['detector_controls'][0]['attempts'][0].update(channel_index=1),
                   lambda n:n['primary']['rows'][0]['photons'].update(value=D('NaN')),lambda n:n['primary'].update(required_rows_admitted=True,rows=[]),
                   lambda n:n.update(all_poisson_terms=1),lambda n:n['primary']['rows'][0]['mean_flux'].update(availability=1,value=1),
                   lambda n:n['actual_temporal_controls'].reverse(),lambda n:n['actual_external_reducer_controls'][0]['reducer'].update(log_joint=0),
                   lambda n:n['detector_controls'][0]['attempts'][0]['batches'][1].update(photons=D('1.01')),lambda n:n.update(grid_count=True)]
        for mutate in mutations:
            candidate=copy.deepcopy(baseline); mutate(candidate)
            with self.subTest(mutate=mutate),self.assertRaises(ValueError): V.native_admission(candidate,request)
    def test_failed_required_attempt_is_preserved_unavailable(self):
        n=native_fixture(); batch=n['detector_controls'][0]['attempts'][0]['batches'][0]
        batch['rows'][0].update(status='work_limit',log_value=None,detection_probability=None)
        self.assertFalse(V.native_admission(n,request_fixture()))
        self.assertEqual(batch['rows'][0]['status'],'work_limit'); self.assertIsNone(batch['rows'][0]['log_value'])
    def admit_reference(self,value,hashes,runtime):
        with mock.patch.object(V,'runtime_admission',side_effect=lambda r,q:r):
            return V.reference_admission(value,request_fixture(),'b'*64,hashes,runtime)
    def test_reference_earned_complete_payload_on_refusal(self):
        value,hashes,runtime=reference_fixture(status='refused'); original=value['payload']
        self.assertFalse(self.admit_reference(value,hashes,runtime)); self.assertIs(value['payload'],original)
    def test_reference_earned_detector_complement_zero_equality_exceeded_rejected(self):
        for sigma in (0,1):
            field='below' if sigma==0 else 'below64'
            for below,error in (('1','0'),('0.75','0.25'),('0.75','0.26')):
                value,hashes,runtime=reference_fixture()
                # Independent synchronized copies isolate the missing positivity gate,
                # rather than failing the already-reviewed payload/prefix identity.
                value['partial']['detector_prefix']=copy.deepcopy(value['payload']['detector_rows'])
                for groups in (value['payload']['detector_rows'],value['partial']['detector_prefix']):
                    groups[sigma]['rows'][0].update({field:below,'below_error':error})
                self.assertEqual(value['payload']['detector_rows'],value['partial']['detector_prefix'])
                gate='discrete' if sigma==0 else 'continuous'
                with self.subTest(sigma=sigma,below=below,error=error),self.assertRaisesRegex(ValueError,gate+' positive reference gate'):
                    self.admit_reference(value,hashes,runtime)
    def test_reference_earned_detector_complement_positive_margin_admitted(self):
        cases=(('0.75','0.24'),
               ('0.75','0.24999999999999999999999999999999999999999999999999999999999'),
               ('0.99999999999999999999999999999999999999999999999999999999999','0'))
        for sigma in (0,1):
            field='below' if sigma==0 else 'below64'
            for below,error in cases:
                value,hashes,runtime=reference_fixture()
                value['partial']['detector_prefix']=copy.deepcopy(value['payload']['detector_rows'])
                for groups in (value['payload']['detector_rows'],value['partial']['detector_prefix']):
                    groups[sigma]['rows'][0].update({field:below,'below_error':error})
                with self.subTest(sigma=sigma,below=below,error=error):
                    self.assertTrue(self.admit_reference(value,hashes,runtime))
                    self.assertEqual(value['payload']['detector_rows'],value['partial']['detector_prefix'])
    def test_reference_failed_detector_complement_diagnostics_retained(self):
        for sigma in (0,1):
            field='below' if sigma==0 else 'below64'
            for below,error in (('1','0'),('0.75','0.25'),('0.75','0.26'),('1.01','0')):
                value,hashes,runtime=reference_refusal_prefix(sigma,0 if sigma==0 else 3,sigma,2)
                attempt=copy.deepcopy(reference_payload_fixture()['detector_rows'][sigma]['rows'][2])
                attempt.update({field:below,'below_error':error})
                density='count_probability' if sigma==0 else 'density64_per_adu'
                density_error='count_probability_error' if sigma==0 else 'density_error_per_adu'
                failure={'stage':'conditional-positive-reference-gate','sigma_index':sigma,'channel_index':2,
                         'lambda':'1','density':attempt[density],'density_error':attempt[density_error],
                         'nondetection':below,'nondetection_error':error,'coarse_refined_attempt':attempt}
                value['partial']['control_failure']=failure
                message=json.dumps(failure,sort_keys=True,separators=(',',':'))
                value['error']={'kind':'reference','stage':failure['stage'],'message':message,
                                'message_truncated':False,'message_sha256':V.digest(message.encode()),'details':copy.deepcopy(failure)}
                saved=copy.deepcopy(value)
                with self.subTest(sigma=sigma,below=below,error=error):
                    self.assertFalse(self.admit_reference(value,hashes,runtime))
                    self.assertIsNone(value['payload'])
                    self.assertEqual(value,saved)
                    self.assertIs(value['partial']['control_failure']['coarse_refined_attempt'],attempt)
    def test_reference_closed_order_positivity_boolean_work(self):
        value,hashes,runtime=reference_fixture(); self.assertTrue(self.admit_reference(value,hashes,runtime))
        mutations=[lambda r:r['payload']['records'].reverse(),lambda r:r['payload']['photon_errors'].__setitem__(0,'-1'),
                   lambda r:r['work'].update(time_frequency_evaluations=30721),lambda r:r['payload'].update(decimal_digits=True),
                   lambda r:r['qualification'].update(original_input_certificate=True),lambda r:r['payload']['records'][0].update(positive_joint_error='1'),
                   lambda r:r['partial']['completed_groups'].reverse(),lambda r:r['payload']['photons'].__setitem__(0,'+1'),
                   lambda r:r['axes'].update(extra=True),lambda r:r['identities']['current_source_sha256'].update(**{'consumer.cpp':'c'*64})]
        for mutate in mutations:
            candidate=copy.deepcopy(value); mutate(candidate)
            with self.subTest(mutate=mutate),self.assertRaises(ValueError): self.admit_reference(candidate,hashes,runtime)
    def test_original_comparison_169_and_real_error_charge(self):
        n=native_fixture(); ref=reference_payload_fixture(); result=C.compare(n,ref)
        self.assertEqual(result['check_count'],169); self.assertTrue(result['accepted'])
        ref['records'][0]['joint_log_error']='1e-10'
        result=C.compare(n,ref); self.assertFalse(result['accepted'])
        gate=next(c for c in result['checks'] if c['name']=='log_joint-0-0')
        self.assertEqual(gate['earned_reference_log_error'],D('1e-10'))


class CaptureAndReceiptTests(unittest.TestCase):
    def test_shared_consumed_and_recovery_quota_preserves_fallback_reserve(self):
        for recovery in (False,True):
            with self.subTest(recovery=recovery),tempfile.TemporaryDirectory() as tmp:
                path=Path(tmp); request=request_fixture(); request['resources']['attempt_store_bytes']=65536+8192
                artifacts={}; capture=C.Capture(path,artifacts,{}); runner=C.Runner(path,request,artifacts,capture)
                blob=b'{"earned":1}'; live=path/'native.stdout.log'; V.write_new(live,blob)
                if recovery:
                    capture.retain('native-stdout',blob,live,parsed=True)
                    consumed=path/artifacts['native-stdout']['path']; os.chmod(consumed,0o600); consumed.write_bytes(b'{"earned":9}')
                remaining=request['resources']['attempt_store_bytes']-65536-runner.store_bytes()
                V.write_new(path/'filler',b'x'*remaining)
                if recovery:
                    capture.terminal(); self.assertTrue(any(f['stage']=='consumed-buffer-recovery' for f in capture.failures))
                    self.assertFalse((path/'native-stdout.recovered.consumed.json').exists())
                else:
                    with self.assertRaises(ValueError): capture.retain('native-stdout',blob,live,parsed=True)
                    self.assertFalse((path/'native-stdout.consumed.json').exists())
                result=C.terminal_record(path,{'request_sha256':'a'*64,'too_large':'x'*(V.MIB+1)},artifacts,runner.room,runner.total_room)
                self.assertFalse(result['numerical_accepted']); self.assertLessEqual((path/'record.json').stat().st_size,65536)
                self.assertLessEqual(runner.store_bytes(),request['resources']['attempt_store_bytes'])

    def test_fallback_checks_actual_total_even_if_reserve_was_externally_consumed(self):
        request=request_fixture(); request['resources']['attempt_store_bytes']=65536+128
        runner=C.Runner(self.path,request,self.artifacts,self.capture)
        V.write_new(self.path/'external',b'x'*request['resources']['attempt_store_bytes'])
        with self.assertRaises(ValueError): C.terminal_record(self.path,{'large':'x'*(V.MIB+1)},self.artifacts,runner.room,runner.total_room)
        self.assertFalse((self.path/'record.json').exists()); self.assertEqual(runner.store_bytes(),request['resources']['attempt_store_bytes'])
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(); self.path=Path(self.tmp.name); self.artifacts={}; self.capture=C.Capture(self.path,self.artifacts,{})
    def tearDown(self):
        for p in self.path.rglob('*'):
            if p.is_dir(): os.chmod(p,0o700)
            else: os.chmod(p,0o600)
        os.chmod(self.path,0o700); self.tmp.cleanup()
    def consume(self,label='native-stdout',blob=b'{"earned":1}'):
        live=self.path/(label+'.log'); V.write_new(live,blob)
        value=self.capture.retain(label,blob,live,parsed=True); return live,value
    def overwrite(self,path,blob): os.chmod(path,0o600); path.write_bytes(blob)
    def test_pass_then_changed_malformed_deleted_stdout_or_stderr(self):
        for label,change in ((label,change) for label in ('native-stdout','reference-output') for change in ('changed','malformed','deleted','stderr')):
            with self.subTest(label=label,change=change),tempfile.TemporaryDirectory() as tmp:
                path=Path(tmp); artifacts={}; capture=C.Capture(path,artifacts,{})
                live=path/('reference.json' if label=='reference-output' else 'native.stdout.log'); blob=b'{"earned":1}'; V.write_new(live,blob)
                parsed=capture.retain(label,blob,live,parsed=True)
                err=path/'native.stderr.log'; V.write_new(err,b'original warning\n'); capture.retain('native-stderr',b'original warning\n',err)
                comparison={'accepted':True,'checks':[{'name':'synthetic-earned-gate','passed':True}]}
                record={'numerical_accepted':True,'failures':[],'comparison':comparison}
                target=err if change=='stderr' else live; os.chmod(target,0o600)
                if change=='deleted': target.unlink()
                else: target.write_bytes(b'{"earned":2}' if change=='changed' else (b'{' if change=='malformed' else b'changed warning\n'))
                C.apply_capture_gate(record,capture)
                self.assertFalse(record['numerical_accepted']); self.assertEqual(record['status'],'failed')
                self.assertIs(record['comparison'],comparison); self.assertTrue(comparison['accepted']); self.assertEqual(parsed,{'earned':1})
                self.assertTrue(any(f['stage']=='raw-output-drift' for f in record['failures']))
                self.assertEqual(V.read(path/artifacts[label]['path']),blob)
    def test_consumed_drift_recovers_once_original_memory(self):
        live,parsed=self.consume(); consumed=self.path/self.artifacts['native-stdout']['path']; self.overwrite(consumed,b'{"earned":9}')
        record={'numerical_accepted':True,'failures':[]}; C.apply_capture_gate(record,self.capture)
        recovered=self.path/self.artifacts['native-stdout']['path']; self.assertEqual(recovered.name,'native-stdout.recovered.consumed.json')
        self.assertEqual(V.read(recovered),b'{"earned":1}'); self.assertEqual(parsed,{'earned':1}); self.assertEqual(V.read(consumed),b'{"earned":9}')
        self.assertFalse(record['numerical_accepted']); self.assertEqual(stat.S_IMODE(recovered.stat().st_mode),0o444)
    def test_malformed_partial_bytes_and_recovery_quota(self):
        live,parsed=self.consume(blob=b'{"partial":'); self.assertIsNone(parsed); self.assertEqual(V.read(self.path/self.artifacts['native-stdout']['path']),b'{"partial":')
        self.overwrite(self.path/self.artifacts['native-stdout']['path'],b'changed')
        self.capture.quota=mock.Mock(side_effect=ValueError('quota'))
        self.capture.terminal(); self.assertTrue(any(f['stage']=='consumed-buffer-recovery' for f in self.capture.failures))
    def test_record_fallback_and_fresh_immutability(self):
        record={'status':'completed','numerical_accepted':True,'request_sha256':'a'*64,'huge':'x'*(V.MIB+1)}
        result=C.terminal_record(self.path,record,{},lambda n:None)
        self.assertEqual(result['interface_id'],'temporal-shared-optical-terminal-refusal/v1'); self.assertFalse(result['numerical_accepted'])
        self.assertLessEqual((self.path/'record.json').stat().st_size,65536); self.assertEqual(V.parse(V.read(self.path/'record.json'))['status'],'failed')
        receipt=self.path/'record.json'; original_bytes=V.read(receipt,65536)
        original_sha=hashlib.sha256(original_bytes).hexdigest(); original_mode=stat.S_IMODE(receipt.stat().st_mode)
        self.assertEqual(original_mode,0o444)
        # Oversize admission fails before exclusive creation; the prior fallback survives.
        with self.assertRaisesRegex(ValueError,'^terminal receipt exceeds1MiB$'):
            C.terminal_record(self.path,record,{},lambda n:None)
        after_oversize=V.read(receipt,65536)
        self.assertEqual(after_oversize,original_bytes)
        self.assertEqual(hashlib.sha256(after_oversize).hexdigest(),original_sha)
        self.assertEqual(stat.S_IMODE(receipt.stat().st_mode),original_mode)
        # A separately valid small record reaches O_EXCL and cannot replace the fallback.
        small={'status':'completed','numerical_accepted':True,'request_sha256':'b'*64}
        self.assertLessEqual(len(V.encode(small)),V.MIB)
        with self.assertRaises(FileExistsError): C.terminal_record(self.path,small,{},lambda n:None)
        after_small=V.read(receipt,65536)
        self.assertEqual(after_small,original_bytes)
        self.assertEqual(hashlib.sha256(after_small).hexdigest(),original_sha)
        self.assertEqual(stat.S_IMODE(receipt.stat().st_mode),original_mode)
    def test_fallback_reference_list_is_bounded_by_actual_encoding(self):
        artifacts={('label'+str(i)+'x'*120):{'path':'p'+str(i)+'x'*250,'bytes':1,'sha256':'a'*64} for i in range(200)}
        record={'request_sha256':'b'*64,'oversize':'x'*(V.MIB+1)}
        result=C.terminal_record(self.path,record,artifacts,lambda n:None)
        self.assertLessEqual(len(V.read(self.path/'record.json')),65536)
        self.assertLessEqual(len(result['artifacts']),128); self.assertEqual(result['unavailable_artifact_count']+len(result['artifacts']),200)
    def test_symlink_path_and_single_slug_refusal(self):
        (self.path/'link').symlink_to(self.path)
        with self.assertRaises(ValueError): V.walk_directory(self.path/'link')
        with self.assertRaises(OSError): V.write_new(self.path/'link',b'{}')
        for slug in ('../x','a/b','','a'*65,'é','-x'):
            with self.assertRaises(ValueError): C.execute(slug)
    def test_sealing_and_exact_json_input_bytes(self):
        V.write_new(self.path/'reference-input.json',V.encode({'row':1})); V.write_new(self.path/'consumer',b'synthetic',mode=0o555)
        self.assertEqual(V.parse(V.read(self.path/'reference-input.json')),{'row':1}); C.seal(self.path)
        self.assertEqual(stat.S_IMODE((self.path/'reference-input.json').stat().st_mode),0o444)
        self.assertEqual(stat.S_IMODE((self.path/'consumer').stat().st_mode),0o555)
    def test_closed_terminal_gates_retain_169_after_drift(self):
        live,parsed=self.consume(); self.overwrite(live,b'{"earned":2}')
        internal={'failures':[],'numerical_accepted':True,'execution_complete':True,'comparison_accepted_before_terminal':True,
                  'comparison_artifact':{'path':'comparisons.json','bytes':10,'sha256':'a'*64},'qualification':'empirical-synthetic-fixed-temporal-shared-optical-control/v1',
                  'inference_status':'blocked','request_sha256':'b'*64,'started_utc':'synthetic','ended_utc':'synthetic'}
        C.apply_capture_gate(internal,self.capture)
        record=C.final_receipt(internal,'synthetic',self.artifacts,self.capture,{},True,True)
        expected=['schema_version','interface_id','created_utc','finished_utc','attempt','request','source_before','source_after','sdk_before','sdk_precompile','sdk_prenative','sdk_prereference','sdk_after','controller_runtime_before','controller_runtime_after','reference_runtime_before','reference_runtime_after','subprocesses','native','reference','comparisons','output_inventory','gates','errors','qualification','seal']
        self.assertEqual(set(record),set(expected)); self.assertEqual(record['gates']['original169_comparisons']['status'],'passed')
        self.assertEqual(record['gates']['raw_output_unchanged']['status'],'failed'); self.assertEqual(record['gates']['numerical']['status'],'failed')
        self.assertEqual(record['comparisons'],internal['comparison_artifact']); self.assertEqual(record['native']['consumed_bytes_sha256'],V.digest(b'{"earned":1}'))


class BoundedProcessTests(unittest.TestCase):
    def setUp(self):
        # Actual helper executable bytes are consumed; no production SDK scanning.
        self.artifact_patch=mock.patch.object(C.Runner,'consumption_map',side_effect=lambda label,args,git:{'files':{'helper':V.file_identity(Path(args[0]).resolve(),67108864)[1]},'errors':[]})
        self.artifact_patch.start()
    def tearDown(self): self.artifact_patch.stop()
    def runner(self,tmp):
        path=Path(tmp); (path/'tmp').mkdir(); request=request_fixture(); artifacts={}; capture=C.Capture(path,artifacts,{})
        return C.Runner(path,request,artifacts,capture),capture
    def test_nonzero_partial_json_and_stderr_retained(self):
        with tempfile.TemporaryDirectory() as tmp:
            runner,capture=self.runner(tmp)
            code='import sys; print("{\\"partial\\":1}"); print("failure diagnostic",file=sys.stderr); sys.exit(7)'
            with self.assertRaises(C.ChildFailure): runner.launch('native',[sys.executable,'-c',code],json_output=True)
            child=runner.children[0]; self.assertEqual(child['returncode'],7); self.assertEqual(child['status'],'nonzero-exit')
            self.assertEqual(capture.parsed['native-stdout'],{'partial':1}); self.assertIn(b'failure diagnostic',capture.buffers['native-stderr'])
            self.assertEqual(child['streams']['stdout']['sha256'],V.digest(capture.buffers['native-stdout']))
    def test_timeout_and_output_overflow_are_refusals(self):
        with tempfile.TemporaryDirectory() as tmp:
            runner,capture=self.runner(tmp)
            # Mock only the watchdog clock; the helper is actually killed and reaped.
            real=C.time.monotonic; calls=[0]
            def clock(): calls[0]+=1; return real()+200 if calls[0]>2 else real()
            with mock.patch.object(C.time,'monotonic',side_effect=clock),self.assertRaises(C.ChildFailure):
                runner.launch('native',[sys.executable,'-c','import time; print("partial",flush=True); time.sleep(10)'])
            self.assertEqual(runner.children[0]['status'],'timeout'); self.assertIsNotNone(runner.children[0]['wait4'])
        with tempfile.TemporaryDirectory() as tmp:
            runner,capture=self.runner(tmp)
            with self.assertRaises(C.ChildFailure): runner.launch('native',[sys.executable,'-c','import sys; sys.stdout.write("x"*1048577); sys.stdout.flush()'])
            child=runner.children[0]; self.assertEqual(child['status'],'output-limit'); self.assertFalse(child['streams']['stdout']['complete'])
            self.assertEqual(len(capture.buffers['native-stdout']),V.MIB)
            closed=C.child_records(runner,capture,False,False)
            slot=next(row for row in closed if row['name']=='native')
            self.assertEqual(slot['stdout']['observed_bytes'],V.MIB+1); self.assertEqual(slot['stdout']['bytes'],V.MIB)
            self.assertEqual(slot['termination_cause'],'stdout_limit')
    def test_before_popen_refusal_has_slot_and_no_fabricated_streams(self):
        with tempfile.TemporaryDirectory() as tmp:
            runner,capture=self.runner(tmp)
            refused={'files':{'helper':None},'errors':[{'role':'helper','message':'changed before launch'}]}
            with mock.patch.object(runner,'consumption_map',return_value=refused),mock.patch.object(C.subprocess,'Popen') as popen,self.assertRaises(C.ChildFailure):
                runner.launch('native',[sys.executable,'-c','pass'])
            popen.assert_not_called(); rows=C.child_records(runner,capture,False,False)
            self.assertEqual([row['name'] for row in rows],runner.request['resources']['schedule']); self.assertEqual(len(rows),28)
            native=next(row for row in rows if row['name']=='native'); self.assertEqual(native['launch_status'],'launch_failed'); self.assertIsNone(native['stdout'])
            self.assertTrue(all(row['parse_error']['message'] for row in rows if row['launch_status']=='not_started'))

    def test_immediate_after_check_detects_consumed_change_before_restoration(self):
        with tempfile.TemporaryDirectory() as tmp:
            runner,capture=self.runner(tmp); helper=Path(tmp)/'mutable-helper'
            original=('#!'+sys.executable+'\nfrom pathlib import Path\nPath(__file__).write_text("changed helper bytes\\n")\nprint("earned output")\n').encode()
            helper.write_bytes(original); helper.chmod(0o700)
            with self.assertRaises(C.ChildFailure): runner.launch('native',[str(helper)])
            child=runner.children[0]; self.assertEqual(child['returncode'],0); self.assertEqual(child['status'],'identity-failed')
            self.assertNotEqual(child['artifacts_before'],child['artifacts_after']); self.assertIn('identity_failure',child)
            helper.write_bytes(original)
            self.assertEqual(child['artifacts_before']['files']['helper']['sha256'],V.digest(original))
            self.assertNotEqual(child['artifacts_after']['files']['helper']['sha256'],V.digest(original))

    def test_interrupt_after_parent_reaped_kills_descendant_and_closes_pipes(self):
        with tempfile.TemporaryDirectory() as tmp:
            runner,capture=self.runner(tmp); reaped=[False]; injected=[False]; killed=[]
            real_wait=C.os.wait4; real_kill=C.os.killpg; real_selector=C.selectors.DefaultSelector
            def wait(pid,flags):
                result=real_wait(pid,flags)
                if result[0]: reaped[0]=True
                return result
            def kill(pid,sig): killed.append((pid,sig,reaped[0])); return real_kill(pid,sig)
            class InterruptingSelector:
                def __init__(self): self.inner=real_selector()
                def register(self,*args): return self.inner.register(*args)
                def unregister(self,*args): return self.inner.unregister(*args)
                def get_map(self): return self.inner.get_map()
                def close(self): return self.inner.close()
                def select(self,timeout):
                    if reaped[0] and not injected[0]: injected[0]=True; raise KeyboardInterrupt('fixture after direct-child exit')
                    return self.inner.select(timeout)
            helper='import os,time\npid=os.fork()\nif pid:\n print(pid,flush=True)\n os._exit(0)\ntime.sleep(10)\n'
            with mock.patch.object(C.os,'wait4',side_effect=wait),mock.patch.object(C.os,'killpg',side_effect=kill),mock.patch.object(C.selectors,'DefaultSelector',InterruptingSelector),self.assertRaises(C.ChildFailure):
                runner.launch('native',[sys.executable,'-c',helper])
            self.assertTrue(injected[0]); self.assertTrue(any(was_reaped for _,_,was_reaped in killed))
            slot=next(row for row in C.child_records(runner,capture,False,False) if row['name']=='native')
            self.assertEqual(slot['launch_status'],'started'); self.assertEqual(slot['termination_cause'],'interrupted')
            self.assertEqual(slot['exit_code'],0); self.assertIsNotNone(slot['user_cpu_seconds'])
            pid=int(capture.buffers['native-stdout'].strip()); status=Path('/proc')/str(pid)/'stat'
            if status.exists(): self.assertEqual(status.read_text().split(') ',1)[1].split()[0],'Z')
    def test_no_undeclared_or_repeated_launch(self):
        with tempfile.TemporaryDirectory() as tmp:
            runner,capture=self.runner(tmp)
            with self.assertRaises(ValueError): runner.launch('undeclared',[sys.executable,'-c','pass'])
            runner.used.add('native')
            with self.assertRaises(ValueError): runner.launch('native',[sys.executable,'-c','pass'])


class MockedOrchestrationTests(unittest.TestCase):
    """Exercise execute() with handwritten bytes and a bounded mocked launcher.

    SDK/discovery/runtime/build operations are expressly mocked transport fixtures,
    not evidence for their real gates. Actual parsing, source tar admission,
    raw capture, typed scientific-shape admission,169comparison and sealing run.
    """
    def run_case(self,failure=None):
        tmp=tempfile.TemporaryDirectory(); root=Path(tmp.name); packet=root/'experiments'/'temporal-shared-optical-control'; packet.mkdir(parents=True)
        request=request_fixture()
        request.update(request_id='temporal-shared-optical-reproducible-request/v1',interface_id='temporal-shared-optical-bounded-controller/v1',status='reviewed-native-synthetic-control',execution=None)
        for name in ('controller.py','verify_sdk.py','consumer.cpp','reference.py'): (packet/name).write_bytes(('synthetic source fixture '+name+'\n').encode())
        (packet/'contract.snapshot.json').write_bytes(b'{"synthetic_transport_fixture":true}\n')
        (packet/'README.md').write_bytes(b'Synthetic fixture; no physics execution.\n')
        raw=numeric_wire(request).encode(); (packet/'request.json').write_bytes(raw)
        (packet/'experiment.json').write_bytes(V.encode({'request_sha256':V.digest(raw)}))
        if failure==('metadata','request-pin'): (packet/'experiment.json').write_bytes(V.encode({'request_sha256':'0'*64}))
        files={str(path.relative_to(root)):path.read_bytes() for path in packet.iterdir()}
        tree=b''.join(('100644 blob '+hashlib.sha1(b'blob '+str(len(blob)).encode()+b'\0'+blob).hexdigest()+'\t'+name+'\0').encode() for name,blob in sorted(files.items()))
        head='a'*40; payload=reference_payload_fixture(); getter={'scope':'mock metadata-only getter fixture'}; base_runner=C.Runner; runners=[]
        class FixtureRunner(base_runner):
            def __init__(self,*args): super().__init__(*args); runners.append(self)
            def launch(self,label,args,git=False,allow=(0,),json_output=False):
                if label not in self.request['resources']['schedule'] or label in self.used: raise AssertionError('fixture undeclared/repeated launch')
                self.used.add(label)
                child={'label':label,'argv':args,'cwd':str(self.attempt),'environment':self.env,'limits':{'fixture':'mocked launch'},'launched':True,'status':'completed',
                       'returncode':0,'wall_seconds':0,'wait4':{'cpu_user_seconds':0,'cpu_system_seconds':0,'maximum_rss_kib':0},'streams':{},'started_utc':'fixture','finished_utc':'fixture'}
                self.children.append(child)
                if failure==('before-popen','native') and label=='native':
                    child.update(launched=False,status='launch-failed',returncode=None,wait4=None,exception={'kind':'IdentityError','message':'fixture before-Popen consumed-file refusal'})
                    raise C.ChildFailure('native before-Popen fixture refusal')
                blob=b''
                if label.startswith('repro-'):
                    if label.endswith('-head'): blob=(head+'\n').encode()
                    elif label.endswith('-branch'): blob=b'codex/synthetic-fixture\n'
                    elif label.endswith('-tree'): blob=tree
                    elif label=='repro-source-archive':
                        with tarfile.open(self.attempt/'source.tar','w') as archive:
                            for name,content in sorted(files.items()):
                                member=tarfile.TarInfo(name); member.size=len(content); member.mode=0o644; archive.addfile(member,io.BytesIO(content))
                elif label.startswith('sdk-'): blob=(self.request['sdk']['revision']+'\n').encode() if label.endswith('-head') else b''
                elif label in ('runtime-before','runtime-after-finally'): blob=V.encode(getter)
                elif label=='discovery': blob=V.encode({'synthetic_discovery_fixture':True})
                elif label=='compiler-version': blob=b'synthetic compiler fixture; no compiler executed\n'
                elif label=='compile': V.write_new(self.attempt/'consumer',b'synthetic ELF identity fixture, never executed',mode=0o555)
                elif label=='native': blob=numeric_wire(native_fixture()).encode()
                elif label=='reference':
                    reference,_,_=reference_fixture(payload=copy.deepcopy(payload)); hashes=C.source_port_hashes(self.attempt/'source'/'experiments'/'temporal-shared-optical-control',self.attempt/'contract.snapshot.json')
                    reference['runtime_before']=copy.deepcopy(getter); reference['runtime_after']=copy.deepcopy(getter)
                    reference['identities'].update(request_sha256=V.digest(raw),reference_script_sha256=hashes['reference.py'],current_source_sha256=hashes)
                    output=numeric_wire(reference).encode(); V.write_new(self.attempt/'reference.json',output)
                    blob=V.encode({'interface_id':'temporal-shared-optical-reference-summary/v1','status':'completed','request_sha256':V.digest(raw),'reference_json_sha256':V.digest(output),'reference_json_bytes':len(output)})
                failed=failure is not None and failure==('nonzero',label)
                if failed: child.update(status='nonzero-exit',returncode=7)
                for stream,content in [('stdout',blob),('stderr',b'fixture retained failure\n' if failed else b'')]:
                    path=self.attempt/(label+'.'+stream+'.log'); item=V.write_new(path,content)
                    item.update(complete=True,observed_bytes=len(content)); child['streams'][stream]=item; self.artifacts[label+'-'+stream]=item
                    self.capture.retain(label+'-'+stream,content,path,parsed=json_output and stream=='stdout')
                if failed: raise C.ChildFailure(label+': literal synthetic nonzero failure')
                return blob
        def sdk_map(request,phase,run):
            # Fixed independent head/status calls are retained even on refusal.
            errors=[]
            for suffix,args in [('head',['rev-parse','--verify','HEAD']),('status',['status','--porcelain=v1','-z','--untracked-files=no'])]:
                try: run('sdk-'+phase+'-'+suffix,['git',*V.GIT_OPTIONS,'-C',request['sdk']['source_root'],*args],git=True)
                except Exception as exc: errors.append({'stage':suffix,'message':str(exc)})
            return {'scope':'mock SDK map fixture, no real SDK qualification','errors':errors}
        patches=[mock.patch.object(C,'ROOT',root),mock.patch.object(C,'PACKET',packet),mock.patch.object(C,'RESULTS',root/'results'/'temporal-shared-optical-control'),
                 mock.patch.object(C,'Runner',FixtureRunner),mock.patch.object(C.sys,'dont_write_bytecode',True),
                 mock.patch.object(V,'request_admission',side_effect=lambda r,p:r),mock.patch.object(V,'fingerprint',side_effect=sdk_map),
                 mock.patch.object(V,'admit_fingerprint',side_effect=lambda value:value),mock.patch.object(V,'artifact_checks',return_value={'scope':'mock artifact fixture'}),
                 mock.patch.object(V,'discovery_admission',side_effect=lambda value,r:value),mock.patch.object(V,'runtime_admission',side_effect=lambda value,r:value),
                 mock.patch.object(C,'controller_runtime',return_value={'scope':'mock controller map fixture','decimal_precision':80}),mock.patch.object(C,'loader_map',return_value={'scope':'mock loader fixture'})]
        try:
            for patch in patches: patch.start()
            result=C.execute('fresh')
            attempt=root/'results'/'temporal-shared-optical-control'/'fresh'
            children=V.parse(V.read(attempt/result['subprocesses']['path']))
            comparisons=V.parse(V.read(attempt/result['comparisons']['path'])) if result['comparisons'] is not None else None
            source_after=V.parse(V.read(attempt/result['source_after']['path'])) if result['source_after'] is not None else None
            outputs={'result':result,'children':children,'comparisons':comparisons,'source_after':source_after,'used':set(runners[0].used) if runners else set(),
                     'launched':[item['label'] for item in runners[0].children] if runners else []}
            self.assertEqual(V.parse(V.read(attempt/'record.json'))['interface_id'],'temporal-shared-optical-attempt-record/v1')
            self.assertEqual(stat.S_IMODE((attempt/'record.json').stat().st_mode),0o444)
            return outputs
        finally:
            for patch in reversed(patches): patch.stop()
            for path in sorted(root.rglob('*'),key=lambda p:len(p.parts)):
                path.chmod(0o700 if path.is_dir() else 0o600)
            root.chmod(0o700); tmp.cleanup()

    def test_complete_synthetic_orchestration_has_exact28_slots_and169checks(self):
        output=self.run_case(); request=request_fixture()
        self.assertEqual(output['result']['gates']['numerical']['status'],'passed')
        self.assertEqual([row['name'] for row in output['children']],request['resources']['schedule']); self.assertEqual(len(output['used']),28)
        self.assertTrue(all(row['launch_status']=='started' for row in output['children']))
        self.assertEqual(output['comparisons']['check_count'],169); self.assertTrue(output['comparisons']['accepted'])

    def test_early_refusal_and_before_popen_refusal_keep_final_and_skipped_slots(self):
        for failure in [('nonzero','repro-initial-head'),('before-popen','native')]:
            with self.subTest(failure=failure):
                output=self.run_case(failure); rows=output['children']
                self.assertEqual(len(rows),28); self.assertEqual(len({row['name'] for row in rows}),28)
                self.assertEqual(output['result']['gates']['numerical']['status'],'failed')
                for label in ('repro-final-head','repro-final-branch','repro-final-status','repro-final-tree','sdk-final-head','sdk-final-status'):
                    self.assertIn(label,output['used'])
                self.assertTrue(any(row['launch_status']=='not_started' and row['parse_error']['message'] for row in rows))
                self.assertIsNone(output['comparisons'])

    def test_final_head_failure_retains_independent_metadata_and_original169(self):
        output=self.run_case(('nonzero','repro-final-head'))
        for label in ('repro-final-branch','repro-final-status','repro-final-tree'): self.assertIn(label,output['used'])
        self.assertEqual(output['source_after']['branch'],'codex/synthetic-fixture'); self.assertIsNone(output['source_after']['head'])
        self.assertEqual(output['result']['gates']['source_unchanged']['status'],'failed')
        self.assertEqual(output['result']['gates']['numerical']['status'],'failed')
        self.assertEqual(output['comparisons']['check_count'],169); self.assertTrue(output['comparisons']['accepted'])
        self.assertEqual(output['result']['gates']['original169_comparisons']['status'],'passed')

    def test_failed_experiment_request_pin_keeps28_unstarted_slots(self):
        output=self.run_case(('metadata','request-pin'))
        self.assertEqual(len(output['children']),28); self.assertEqual(output['used'],set())
        self.assertTrue(all(row['launch_status']=='not_started' and row['stdout'] is None and row['parse_error']['message'] for row in output['children']))
        self.assertEqual(output['result']['gates']['numerical']['status'],'failed')


class SourceAdmissionTests(unittest.TestCase):
    def test_independent_final_metadata_survives_failed_head_child(self):
        calls=[]; head='a'*40
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); blob=b'reviewed source\n'; (root/'file.py').write_bytes(blob)
            oid=hashlib.sha1(b'blob '+str(len(blob)).encode()+b'\0'+blob).hexdigest(); tree=('100644 blob '+oid+'\tfile.py\0').encode()
            def run(label,args,git=False,allow=(0,)):
                calls.append((label,args))
                if label.endswith('-head'): raise C.ChildFailure('literal failed final HEAD')
                if label.endswith('-branch'): return b'codex/fixture\n'
                if label.endswith('-status'): return b''
                return tree
            value=C.repro_map(root,'final',run,head)
            self.assertEqual([label for label,_ in calls],['repro-final-head','repro-final-branch','repro-final-status','repro-final-tree'])
            self.assertIsNone(value['head']); self.assertEqual(value['tree_revision'],head)
            self.assertEqual(value['sources']['file.py']['sha256'],V.digest(blob))
            with self.assertRaises(V.IdentityError): C.admit_source(value)
            calls.clear(); value=C.repro_map(root,'initial',run)
            self.assertEqual(len(calls),3); self.assertIsNone(value['tree_sha256']); self.assertIsNone(value['sources'])

    def test_controller_runtime_mapping_sets_and_missing_origins(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'already-in-inventory.so'; path.write_bytes(b'mapped module fixture')
            module=types.SimpleNamespace(__spec__=types.SimpleNamespace(origin=str(path),submodule_search_locations=None),__cached__=None)
            no_map=C.controller_runtime_collect([('fixture',module)],'','Threads:\t1\n')
            with_map=C.controller_runtime_collect([('fixture',module)],'00-01 r-xp 0 00:00 0 '+str(path),'Threads:\t1\n')
            self.assertEqual(no_map['inventory'],with_map['inventory']); self.assertNotEqual(no_map['mapped_libraries'],with_map['mapped_libraries'])
            self.assertEqual(with_map['decimal_precision'],80)
            path.unlink()
            with self.assertRaises(FileNotFoundError): C.controller_runtime_collect([], '00-01 r-xp 0 00:00 0 '+str(path),'Threads:\t1\n')
            with self.assertRaises(FileNotFoundError): C.controller_runtime_collect([('fixture',module)],'','Threads:\t1\n')
            with self.assertRaises(ValueError): C.controller_runtime_collect([], '00-01 r-xp 0 00:00 0 '+str(path)+' (deleted)','Threads:\t1\n')

    def test_controller_runtime_aggregate_guard_precedes_file_consumption(self):
        with tempfile.TemporaryDirectory() as tmp:
            oversized=Path(tmp)/'oversized.so'
            # Sparse metadata fixture: no65MiB buffer or scientific array is read.
            with oversized.open('wb') as stream: stream.truncate(67108865)
            module=types.SimpleNamespace(__spec__=types.SimpleNamespace(origin=str(oversized),submodule_search_locations=None),__cached__=None)
            real_identity=V.file_identity; consumed=[]
            def identity(path,*args): consumed.append(Path(path)); return real_identity(path,*args)
            with mock.patch.object(V,'file_identity',side_effect=identity),self.assertRaises(ValueError):
                C.controller_runtime_collect([('fixture',module)],'','Threads:\t1\n')
            self.assertNotIn(oversized,consumed)

    def test_reference_declared_aggregate_refuses_before_any_inventory_read(self):
        request=request_fixture(); expected=request['reference_runtime']; types=request['shape_contract']['runtime_types']
        value={key:None for key in types['RuntimeFingerprint']}
        value.update(schema_version=1,interface_id='temporal-shared-optical-reference-runtime/v1',actual_thread_count=1,thread_environment=request['resources']['environment'])
        value['python']={key:None for key in types['PythonIdentity']}
        for key in ('version','implementation','requested_executable','resolved_executable','prefix','base_prefix'): value['python'][key]=expected[key]
        value['python'].update(executable_sha256=expected['executable_sha256'],dont_write_bytecode=True,executable_bytes=1)
        value['mpmath']={key:None for key in types['MpmathIdentity']}
        value['mpmath'].update(version='1.3.0',backend='python',decimal_digits=60,gmpy2_loaded=False,original_subset_map_sha256=expected['legacy_package_inventory_sha256'],binary_precision_bits=203)
        value['inventory']=[{'path':'/absolute-fixture-'+str(i),'bytes':16777217,'sha256':'a'*64} for i in range(4)]
        with mock.patch.object(V,'file_identity') as reader,self.assertRaises(ValueError): V.runtime_admission(value,request)
        reader.assert_not_called()
    def test_meaningful_changed_committed_bytes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); content=b'reviewed-source\n'; (root/'file.py').write_bytes(content)
            oid=hashlib.sha1(b'blob '+str(len(content)).encode()+b'\0'+content).hexdigest()
            tree=('100644 blob '+oid+'\tfile.py\0').encode()
            actual=C.source_map(root,tree); self.assertEqual(actual['file.py']['sha256'],V.digest(content))
            (root/'file.py').write_bytes(b'reviewed-sourcf\n')
            with self.assertRaises(ValueError): C.source_map(root,tree)
    def test_archive_traversal_and_links_reject(self):
        for name,kind in (('../escape','file'),('link','link')):
            with tempfile.TemporaryDirectory() as tmp:
                attempt=Path(tmp)
                with tarfile.open(attempt/'source.tar','w') as tar:
                    member=tarfile.TarInfo(name)
                    if kind=='link': member.type=tarfile.SYMTYPE; member.linkname='../escape'
                    tar.addfile(member)
                with self.assertRaises(ValueError): C.unpack_snapshot(attempt,{'sources':{}})
    def test_actual311_inventory_changed_bytes_and_missed_cmake(self):
        # Entirely local scanner fixture; no compiler/discovery/physics execution.
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); request=request_fixture(); sdk=request['sdk']; sdk['source_root']=sdk['artifacts_root']=str(root)
            sources={}
            names=['src/f'+str(i)+'.rs' for i in range(306)]+['cpp/test.cmake',*sdk['source_explicit']]
            for name in names:
                path=root/name; path.parent.mkdir(parents=True,exist_ok=True); content=('reviewed '+name+'\n').encode(); path.write_bytes(content); sources[name]=V.digest(content)
            (root/'cpp/include').mkdir(parents=True); sdk['headers']=[]; sdk['guides']=[]
            for name in ('archive','cli','compiler','standard_library'):
                path=root/'artifacts'/name; path.parent.mkdir(exist_ok=True); content=('synthetic '+name).encode(); path.write_bytes(content)
                sdk[name]['path']=str(path); sdk[name]['sha256']=V.digest(content)
            loader=root/'artifacts'/'loader'; loader.write_bytes(b'synthetic loader'); request['native_dependencies_proposed']['loader_path']=str(loader); request['native_dependencies_proposed']['loader_sha256']=V.digest(b'synthetic loader')
            manifest={'git_head':sdk['revision'],'git_status':'','sources':sources,'profile':'synthetic-scanner-fixture'}
            build={k:v for k,v in manifest.items() if k not in ('build_id','git_head','git_status')}
            sdk['build_id']=manifest['build_id']=V.digest(json.dumps(build,sort_keys=True,separators=(',',':')).encode())
            path=root/sdk['manifest']['path']; path.parent.mkdir(parents=True); raw=json.dumps(manifest).encode(); path.write_bytes(raw); sdk['manifest']['sha256']=V.digest(raw)
            def run(label,args,git=False): return (sdk['revision']+'\n').encode() if label.endswith('-head') else b''
            value=V.fingerprint(request,'initial',run); V.admit_fingerprint(value); self.assertEqual(len(value['sources']),311); self.assertIn('cpp/test.cmake',value['sources'])
            (root/'cpp/test.cmake').write_bytes(b'meaningfully changed cmake\n')
            refused=V.fingerprint(request,'final',run)
            with self.assertRaises(V.IdentityError): V.admit_fingerprint(refused)
            self.assertEqual(refused['sources']['cpp/test.cmake']['sha256'],V.digest(b'meaningfully changed cmake\n')); self.assertEqual(len(refused['sources']),311)
            self.assertIsNotNone(refused['artifacts']['archive']); self.assertIsNotNone(refused['artifacts']['compiler'])
            (root/'cpp/test.cmake').unlink(); (root/'cpp/other.cmake').write_bytes(b'new')
            with self.assertRaises(ValueError): V.admit_fingerprint(V.fingerprint(request,'final',run))
    def test_fifo_read_refuses_without_waiting(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'fifo'; os.mkfifo(path)
            with self.assertRaises(ValueError): V.read(path)
    def test_runtime_complete_sources_modules_maps_and_actual_bytes(self):
        # Getter-shaped local source identities; no mpmath import/scientific call.
        with tempfile.TemporaryDirectory() as tmp:
            base=Path(tmp); package=base/'mpmath'; package.mkdir(); info=base/'mpmath-1.3.0.dist-info'; info.mkdir()
            contents={'mpmath/__init__.py':b'original init\n','mpmath/ctx.py':b'original context\n','mpmath-1.3.0.dist-info/LICENSE':b'license\n',
                      'mpmath-1.3.0.dist-info/METADATA':b'metadata\n','mpmath-1.3.0.dist-info/RECORD':b'record\n','mpmath-1.3.0.dist-info/WHEEL':b'wheel\n',
                      'mpmath-1.3.0.dist-info/top_level.txt':b'mpmath\n','python':b'python fixture\n','pyvenv.cfg':b'venv fixture\n','libfixture.so':b'library fixture\n'}
            for name,blob in contents.items(): (base/name).write_bytes(blob)
            inventory=[{'path':str(base/name),'bytes':len(blob),'sha256':V.digest(blob)} for name,blob in sorted(contents.items())]
            lookup={item['path']:item for item in inventory}; request=request_fixture(); expected=request['reference_runtime']
            expected['resolved_executable']=expected['requested_executable']=str(base/'python'); expected['prefix']=expected['base_prefix']=str(base)
            expected['executable_sha256']=lookup[str(base/'python')]['sha256']; expected['pyvenv_cfg_sha256']=lookup[str(base/'pyvenv.cfg')]['sha256']; expected['mpmath_init_sha256']=lookup[str(package/'__init__.py')]['sha256']
            selected={item['path']:item['sha256'] for item in inventory if '/mpmath/' in item['path'] or '.dist-info/' in item['path']}
            expected['legacy_package_inventory_sha256']=V.digest(json.dumps(selected,sort_keys=True,separators=(',',':')).encode())
            python={key:expected[key] for key in ('version','implementation','requested_executable','resolved_executable','prefix','base_prefix')}
            python.update(executable_bytes=lookup[str(base/'python')]['bytes'],executable_sha256=expected['executable_sha256'],pyvenv_cfg=lookup[str(base/'pyvenv.cfg')],dont_write_bytecode=True)
            runtime={'schema_version':1,'interface_id':'temporal-shared-optical-reference-runtime/v1','python':python,
                     'mpmath':{'version':'1.3.0','module_origin':str(package/'__init__.py'),'module_root':str(package),'backend':'python','decimal_digits':60,'binary_precision_bits':203,'gmpy2_loaded':False,'original_subset_map_sha256':expected['legacy_package_inventory_sha256']},
                     'thread_environment':request['resources']['environment'],'actual_thread_count':1,
                     'imported_modules':[{'name':'mpmath','kind':'source','origin':str(package/'__init__.py'),'namespace_paths':[]},{'name':'sys','kind':'built-in','origin':None,'namespace_paths':[]}],
                     'mapped_libraries':[lookup[str(base/'libfixture.so')]],'inventory':inventory,
                     'inventory_sha256':V.digest(json.dumps(inventory,sort_keys=True,separators=(',',':')).encode()),
                     'numeric_configuration':{'arithmetic_backend':'mpmath-libmp-python','decimal_digits':60,'binary_precision_bits':203,'active_native_numeric_backends':[],'gmpy2_loaded':False,'byteorder':'little'}}
            self.assertIs(V.runtime_admission(runtime,request),runtime)
            changes=[lambda r:r['python'].update(executable_bytes=1),lambda r:r['imported_modules'].reverse(),
                     lambda r:r['mapped_libraries'][0].update(sha256='c'*64),lambda r:r['numeric_configuration'].update(active_native_numeric_backends=['foreign']),
                     lambda r:r.update(actual_thread_count=True),lambda r:r['thread_environment'].update(OMP_NUM_THREADS='2')]
            for change in changes:
                changed=copy.deepcopy(runtime); change(changed)
                with self.subTest(change=change),self.assertRaises(ValueError): V.runtime_admission(changed,request)
            discovered=[0]; actual_scandir=V.os.scandir
            class Entry:
                path=str(package/'unconsumed-extra.txt')
                def is_symlink(self): return False
                def is_dir(self,follow_symlinks=False): return False
                def is_file(self,follow_symlinks=False): return True
            class Entries:
                def __enter__(self): return self
                def __exit__(self,*args): return False
                def __iter__(self):
                    for _ in range(8192): discovered[0]+=1; yield Entry()
            def scandir(path): return Entries() if Path(path)==package else actual_scandir(path)
            with mock.patch.object(V.os,'scandir',side_effect=scandir),self.assertRaises(ValueError): V.runtime_admission(runtime,request)
            self.assertEqual(discovered[0],4097)  # Guard during discovery, before exhausting8192.
            (package/'ctx.py').write_bytes(b'meaningfully changed context\n')
            with self.assertRaises(V.IdentityError): V.runtime_admission(runtime,request)


if __name__=='__main__': unittest.main()
