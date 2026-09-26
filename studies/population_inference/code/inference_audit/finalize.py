from pathlib import Path
import hashlib,json,datetime
R=Path(__file__).resolve().parents[3];O=R/'phase2/inference_audit'
def record(p):return {'path':str(p.relative_to(R)),'bytes':p.stat().st_size,'sha256':hashlib.file_digest(p.open('rb'),'sha256').hexdigest()}
inputs=[R/('scripts/phase2/hierarchy/'+n+'.py') for n in ['core','conditional','prepare_data','validate','validate_conditional','coverage','inference','forward_summary']]+[R/('docs/phase2/'+n+'.md') for n in ['conditional-prediction-plan','generative-model','independent-engine-validation']]+[R/'phase2/hierarchy/conditional-cohort.json',R/'phase2/official/portable_pilot/double_parameters_and_hessian.json',R/'phase2/independent_flux/independent_refits.json',R/'phase2/hierarchy/coverage/summary.json',R/'phase2/hierarchy/coverage/replicates.json',R/'phase2/hierarchy/validation/results.json',R/'phase2/hierarchy/validation/manifest.json',R/'phase2/hierarchy/uv.lock']
for name in ['conditioned-multistart-best','recovered-mask']:
 inputs.extend([R/f'phase2/hierarchy/data/{name}/data.npz',R/f'phase2/hierarchy/data/{name}/rows.csv'])
entries=[record(p) for p in inputs]
for p in inputs:
 if p.suffix in ['.py','.md']:
  h=hashlib.file_digest(p.open('rb'),'sha256').hexdigest();q=O/'provenance'/(h+p.suffix)
  if not q.exists():q.write_bytes(p.read_bytes())
(O/'final_input_manifest.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'files':entries,'heldout_prediction_outputs_read':False},indent=2)+'\n')
tail=json.load(open(O/'stable_tail_validation.json'));cov=json.load(open(O/'covariance_audit.json'));assert tail['passed'];assert tail['core_sha256']==record(R/'scripts/phase2/hierarchy/core.py')['sha256'];assert all(r['frozen_membership_identical'] for r in cov['cohorts'])
old=json.load(open(R/'phase2/assumptions/outputs_manifest.json'));bad=[]
for f in old['files']:
 if record(R/f['path'])['sha256']!=f['sha256']:bad.append(f['path'])
assert not bad,bad
status={'passed':True,'stable_tail_core_hash_matches':True,'previous_measurement_audit_outputs_unchanged':len(old['files']),'fixed_cohort_size':1063,'fixed_test_size':213,'real_benchmark_scores_read':False,'core_edit_authorized_by_parent':True,'core_edit_scope':'scaled negative-tail normal CDF helper and algebraically equivalent exGaussian density/CDF uses','resume':'See docs/phase2/inference-audit.md','status':'bounded audit complete; paused work requested by user, parent owns goal pause'}
(O/'verification.json').write_text(json.dumps(status,indent=2)+'\n')
outputs=sorted([p for p in O.rglob('*') if p.is_file() and p.name!='outputs_manifest.json' and p.suffix!='.log']+list((R/'scripts/phase2/inference_audit').glob('*.py'))+[R/'docs/phase2/inference-audit.md'])
(O/'outputs_manifest.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'files':[record(p) for p in outputs],'other_changed_artifacts':[record(R/'scripts/phase2/hierarchy/core.py'),record(R/'phase2/hierarchy/validation/results.json'),record(R/'phase2/hierarchy/validation/manifest.json')]},indent=2)+'\n')
print(json.dumps(status,indent=2));print('Manifested',len(outputs),'files')
