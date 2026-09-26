from pathlib import Path
import json,hashlib,shutil,ast,difflib,importlib.util
R=Path('/home/szymon/Documents/ChatGPT/supernova');P=R/'runs/research_2026_09_26/raisin_profile_solver_review/csp-filter-interpretation/native-design/timing-followon-design/prospective-engineering';O=P/'noiseless-checker-recovery';(O/'gates').mkdir(exist_ok=False)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
s=(P/'run_engineering_v5_fit_readme.py').read_text();n=s.replace("'freeze-fit-readme.json'","'freeze-checker-resume.json'").replace("P/'fits-readme'","P/'fits-resume'").replace("P/'fits-readme/NIR_true12'","P/'fits-resume/NIR_true12'").replace("'-failure-fit-readme.json'","'-failure-checker-resume.json'")
n=n.replace("int(x['CUTFLAG_SNANA'])==0","int(x['CUTFLAG_SNANA'])==3")
n=n.replace("float(np.float32(x['MJD']))","float(x['MJD'])")
n=n.replace('fitter native R4 measurement mismatch','fitter native FITS measurement mismatch')
n=n.replace('same native data/order in every callback; native fitter may R4-round exposure MJD.','same native data/order in every callback; FITS exposure MJD remains R8.')
n=n.replace("save(w/'gate.json',", "save(P/'noiseless-checker-recovery/gates'/(w.name+'-gate.json'),")
a=" original=run_fit('noiseless_joint_original',data,'grizJH',peak_step=2,original=True);instrumented=run_fit('noiseless_joint',data,'grizJH',peak_step=2)"
b=" original=P/'fits-readme/noiseless_joint_original';instrumented=P/'fits-readme/noiseless_joint' # reuse exact completed joint jobs; no rerun"
assert a in n;n=n.replace(a,b);ast.parse(n);(P/'run_engineering_v5_checker_resume.py').write_text(n)
(O/'executor.patch').write_text(''.join(difflib.unified_diff(s.splitlines(True),n.splitlines(True),fromfile='run_engineering_v5_fit_readme.py',tofile='run_engineering_v5_checker_resume.py')))
# Only existing-log checking: importing this script does not invoke native helpers.
spec=importlib.util.spec_from_file_location('resume_check',P/'run_engineering_v5_checker_resume.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
a,b=m.compare(P/'fits-readme/noiseless_joint_original',P/'fits-readme/noiseless_joint',117,['1'],0,0)
for c in a:
 assert len(a[c])==len(b[c])
 for x,y in zip(a[c],b[c]):assert m.np.array_equal(x['array'],y['array']) and m.np.array_equal(x['W'],y['W']) and x['objective']==y['objective']
_,detail=m.block_review(P/'fits-readme/noiseless_joint',117,['1'])
truth=json.loads((P/'generation-gate.json').read_text())['branches']['noiseless'][0]['DLMU_true']
for d in detail:
 assert abs(d['D']-truth)<=.001 and abs(d['peak']-57707.80078125)<=.01 and d['max_standardized_model_residual']<=.02
summary={'pass_existing_joint_checks':True,'not_full_noiseless_gate_NIR_pending':True,'completed_joint_jobs_reused_not_rerun':True,'truth_D':truth,'detail':detail,'D_error':detail[0]['D']-truth,'peak_error':detail[0]['peak']-57707.80078125,'strict_support_record_count':1,'support_duplicate_guard_unchanged':True,'native_calls':0}
save(O/'existing-joint-check.json',summary)
shutil.copyfile('/tmp/prepare_checker_resume.py',O/'prepare.py');activity=json.loads((P/'native-activity-v5.json').read_text());shutil.copyfile(P/'native-activity-v5.json',O/'activity-at-amendment.json')
proof={'status':'Additive checker correction and resume freeze; no new native execution','scientific_protocol_sha256':sha(P/'protocol-v5.json'),'semantic_corrections':[{'old':'CUTFLAG_SNANA==0','new':'CUTFLAG_SNANA==3 and unchanged ERRFLAG_FIT==0','source':'snana.car:1679 identifies passbits0,1;21625 sets bit0 for successful SN cuts;9337 sets bit1 on last-iteration ERRFLAG0. Zero meant neither condition passed, so the original success checker was inverted.','scope':'Correct interpretation of frozen success intent; no event/epoch cut or failure threshold relaxed.'},{'old':'Compare observed MJD to float32(PHOT.MJD)','new':'Compare observed MJD exactly to double(PHOT.MJD)','source':'sndata.h208 MJD double; sntools_dataformat_fits.c614 writes1D and3398 reads RD_SNFITSIO_DBL; snana.car808 defines SNLC8_MJD. Existing117 row logs independently match native FITS double exactly.','scope':'FITS route differs from previously audited rounded text/plot export; no epoch changed or tolerance introduced. Header PEAKMJD, flux and error remain native R4.'}],'support':'Current authoritative joint log has one summary; emitter registers exactly one atexit callback. Strict missing/duplicate/overflow/nonfinite/bad-call checks remain unchanged. Earlier apparent duplication was from printing the same line twice during inspection.','reuse':'Completed original/support joint logs, FITRES, NML and execution records are hash-bound and rechecked. Resume invokes only the pending noiseless NIR job, then exactly the already declared noisy jobs if their stage is released. It never regenerates or reruns the successful joint jobs.','outputs':'New native work uses fits-resume/ at same relative depth. Checker gate reports go to this additive review folder, not old fit folders. Original failure record persists.','thresholds':'All distance, peak, covariance, stationarity, mask, support and model-residual tolerances unchanged. Existing joint checks pass; pending NIR prevents a full noiseless pass.','full_null':'full_null_check.py independently compares every saved model/entry/measurement/weight/objective record, all FITRES science rows, exact HEAD arrays and same PHOT file for true-peak versus null arms. It is prepared before noisy outcomes; it is not yet an executed null result.','budget':{'cap_seconds':120,'carried_native_seconds':sum(x['wall_seconds'] for x in activity)},'release':'Root release must bind original protocol-v5.json and freeze-checker-resume.json. No native execution in preparation.'}
save(O/'amendment.json',proof)
frozen=json.loads((P/'freeze-fit-readme.json').read_text());frozen['status']='Corrected source semantics and reuse-only noiseless resume';frozen['authoritative_executor']='run_engineering_v5_checker_resume.py';frozen['external_execution_release_required']['freeze_sha256']='SHA256 of freeze-checker-resume.json'
paths=[P/'run_engineering_v5_checker_resume.py']+[x for x in O.iterdir() if x.is_file()]
# Frozen completed native artifacts, but not writable postprocessing gate files.
for branch in ['noiseless_joint_original','noiseless_joint']:
 paths += [x for x in (P/'fits-readme'/branch).iterdir() if x.is_file()]
for q in paths:frozen['files'][str(q.relative_to(R))]=sha(q)
save(P/'freeze-checker-resume.json',frozen)
for name,h in frozen['files'].items():assert sha(R/name)==h,name
print(json.dumps({'protocol_sha256':sha(P/'protocol-v5.json'),'amendment_sha256':sha(O/'amendment.json'),'executor_sha256':sha(P/'run_engineering_v5_checker_resume.py'),'freeze_sha256':sha(P/'freeze-checker-resume.json'),'files':len(frozen['files']),'null_checker_sha256':sha(O/'full_null_check.py'),'existing_joint':summary},indent=2))
