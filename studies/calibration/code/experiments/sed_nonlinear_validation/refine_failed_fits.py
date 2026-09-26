"""Separate finite-difference refinement; original objective/gates unchanged."""
from pathlib import Path
import hashlib,importlib.util,json,shutil
import numpy as np
from scipy.linalg import solve_triangular

ROOT=Path(__file__).resolve().parents[3]
BASE=Path(__file__).resolve().parent
source=ROOT/'scripts/research_2026_09_26/sed_nonlinear_validation.py'
spec=importlib.util.spec_from_file_location('v',source);v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
primary=json.loads((BASE/'result.json').read_text());failed=primary['failed_ids']
assert failed and not primary['errors']
RESCUE=BASE/'refined_steps';RESCUE.mkdir(exist_ok=False);(RESCUE/'objects').mkdir()
shutil.copy2(BASE/'frozen-coefficients.npz',RESCUE/'frozen-coefficients.npz')
shutil.copy2(BASE/'rescue_protocol.md',RESCUE/'protocol.md')
shutil.copy2(Path(__file__),RESCUE/'executed_source.py')
original_engine=v.gate.Engine
class FineEngine(original_engine):
    def jac(self,theta,mode='nominal',factor=1.):
        return super().jac(theta,mode,factor=.1*factor)
v.gate.Engine=FineEngine
v.OUT=RESCUE;v.initialize_worker()
config=json.loads((BASE/'freeze.json').read_text());sensitivity=set(config['official_mean_sensitivity_ids'])
reports=[]
for cid in failed:
    info=v.one_object((cid,cid in sensitivity))
    result=json.loads((RESCUE/'objects'/f'{cid}.json').read_text())
    assert result['status']=='complete'
    with np.load(v.BASE/'objectives'/f'objective_{cid}.npz') as d:
        nominal={k:d[k] for k in ['MJD','band','model_flux','data_flux','data_fluxerr','zHEL','MWEBV','parameters_x0_x1_c_t0','frozen_flux_covariance']}
    engine=FineEngine(nominal,v.WORKER['model'],v.WORKER['bands'],v.WORKER['offsets'],v.WORKER['coefficients'],v.WORKER['observer'])
    chol=np.linalg.cholesky(nominal['frozen_flux_covariance']);white=lambda x:solve_triangular(chol,x,lower=True)
    targets={'native_noiseless':engine.flux(np.zeros(4)),'observed_flux':nominal['data_flux'],'official_mean_sensitivity':nominal['model_flux']}
    checks=[]
    for fit in result['fit_details']:
        for start_index,start in enumerate(fit.get('starts',[])):
            theta=np.array(start['theta']);mode=fit['mode'];target=targets[fit['target']]
            residual=white(engine.flux(theta,mode)-target)
            gradients=[]
            for factor in [1.,.5,.25]:
                jw=white(engine.jac(theta,mode,factor));u,s,vt=np.linalg.svd(jw,full_matrices=False)
                gradients.append(float(np.linalg.norm(u.T@residual)))
            directions=vt.T/s
            objective=lambda x: .5*float(white(engine.flux(x,mode)-target)@white(engine.flux(x,mode)-target))
            cost_gradients=[]
            for epsilon in [1e-3,3e-4,1e-4,3e-5]:
                g=np.array([(objective(theta+epsilon*directions[:,j])-objective(theta-epsilon*directions[:,j]))/(2*epsilon) for j in range(4)])
                cost_gradients.append(float(np.linalg.norm(g)))
            passed=all(x<1e-4 for x in gradients) and all(x<1e-4 for x in cost_gradients[-2:])
            checks.append(dict(target=fit['target'],mode=mode,start=start_index,jacobian_gradient_norms=gradients,objective_directional_gradient_norms=cost_gradients,gate=passed))
    original=json.loads((BASE/'objects'/f'{cid}.json').read_text())
    old={(x['target'],x['mode']):x for x in original['responses']}
    shifts=[abs(x['nonlinear_standardized_mag']-old[(x['target'],x['mode'])]['nonlinear_standardized_mag']) for x in result['responses']]
    report=dict(CID=cid,original_gate=original['all_gates_pass'],refined_original_gates_pass=result['all_gates_pass'],extra_gates_pass=all(x['gate'] for x in checks),max_standardized_response_change_mag=max(shifts),checks=checks)
    reports.append(report)
    print(json.dumps({k:value for k,value in report.items() if k!='checks'}),flush=True)
resolution=dict(initial_failure_count=len(failed),refined_count=len(reports),remaining_failed_ids=[r['CID'] for r in reports if not(r['refined_original_gates_pass'] and r['extra_gates_pass'])],records=reports,
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),protocol_sha256=hashlib.sha256((RESCUE/'protocol.md').read_bytes()).hexdigest())
v.atomic_json(RESCUE/'refinement-verification.json',resolution)
if resolution['remaining_failed_ids']:raise SystemExit('Refinement did not pass; preserve failure.')
VIEW=BASE/'resolved';VIEW.mkdir(exist_ok=False);(VIEW/'objects').mkdir()
for name in ['cohort.csv','contrast-membership.csv','freeze.json','protocol.md','frozen-coefficients.npz']:
    (VIEW/name).symlink_to(Path('..')/name)
for cid in v.pd.read_csv(BASE/'cohort.csv',dtype={'CID':str}).CID:
    origin=RESCUE if cid in failed else BASE
    for suffix in ['.json','-predictions.npz']:
        (VIEW/'objects'/f'{cid}{suffix}').symlink_to(origin/'objects'/f'{cid}{suffix}')
v.OUT=VIEW;v.aggregate()
resolved=json.loads((VIEW/'result.json').read_text());resolved['numerical_refinement']=resolution
v.atomic_json(VIEW/'result.json',resolved)
manifest=json.loads((VIEW/'manifest.json').read_text());manifest['result_sha256']=v.sha(VIEW/'result.json');manifest['numerical_refinement_source_sha256']=v.sha(__file__);manifest['numerical_refinement_protocol_sha256']=v.sha(RESCUE/'protocol.md');manifest['numerical_refinement_verification_sha256']=v.sha(RESCUE/'refinement-verification.json')
v.atomic_json(VIEW/'manifest.json',manifest)
