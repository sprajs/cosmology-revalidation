from pathlib import Path
import json,hashlib,shutil
H=Path(__file__).resolve().parent;Q=H.parent;P=Q.parent;R=Path('/home/szymon/Documents/ChatGPT/supernova');sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
shutil.copyfile(Q/'native-activity.json',H/'carried-native-activity.json');carry=json.loads((H/'carried-native-activity.json').read_text());spent=sum(x['wall_seconds'] for x in carry);assert spent==11.42092445801245
proto=json.loads((Q/'execution-protocol.json').read_text())
proto.update(status='Additive start-only stage protocol; each stage requires root release; no execution yet',identity='Exact restricted-domain estimator plus opt-in first-iteration MINUIT-local start test. Shared native INIVAL/prior center unchanged.',parent_execution_protocol_sha256=sha(Q/'execution-protocol.json'),patch_sha256=sha(H/'start-only.patch'),binary_sha256=sha(H/'build/bin/snlc_fit.exe'),
 stages={'identity':'New binary with PROSP_MINUIT_PEAK_SHIFT absent and0 for nominal joint12 and joint9; all4 runs must be exactly identical to each corresponding already completed saved scientific state, FITRES and domain/support records before baseline reuse.',
 'starts':'New minus/plus localMINUIT±2 jobs, both exact nominal joint12 NML initialization. Stored MINUIT readback, unchanged first common prior/entry state, all original numerical/support/convergence gates. No NML±2 substitution.',
 'nir':'Only after identity+starts gates: reuse nominal joint12 peak; run true/estimated NIR12/9 and postinit D±.2, then full null12. Offset environment absent for every NIR/default job. Full null before effect-table creation.'},
 resources={'workers':1,'native_total_seconds':120,'per_process_seconds':40,'prior_used_seconds':spent,'remaining_seconds_at_freeze':120-spent,'planned_additional_processes':15,'identity_jobs':4,'true_start_jobs':2,'NIR_jobs':9,'build_separate_seconds':json.loads((H/'build-result.json').read_text())['wall_seconds'],'generation_calls':0,'activity':'carried-native-activity.json is immutable prefix; native-activity.json appends every subsequent attempted process including failures'},
 source_state='MNPARM gets a separate local scalar. Common INIVAL, prior center, native bounds, masks, covariance initialization and other parameters remain untouched; only iteration1 peak entry is offset. Enabled-only MNPOUT readback verifies actual storage.',
 weight_interpretation='ITER1 diagonal model variance is parameter-dependent. W at different parameter points need not equal. Same prehook deterministic state and objective function at common coordinates are preserved by identical NML/source before the hook, with first CSP_ENTRY and pre-grid bounds exactly compared. Later C/prior may follow different histories; this is iterative-estimator path stability, not global likelihood uniqueness.',
 actual_start_gates={'first_shift_days':[-2,2],'first_shift_precision_tolerance_days':0.004,'later_shift':'exact0 relative to each iteration native source','first_shared_CSP_ENTRY':'all tokens exact against saved nominaljoint12','first_pregrid_bounds_initializer':'all tokens exact','prior_center':'every callback common prior center equals source entry; first is identical to nominal','stored_parameter':'PKMJD, positive internalindex, fixed lower/upper bounds unchanged; MNPOUT stored and proposed agree to original.004day tolerance','coverage':'one proposed and one MNPOUT record for everyCID and everyiteration; missing/duplicate/abort records fail'},
 mutation_scope='New root-level P/fits-start-only/ and P/derived-start-only/ plus this owned directory. All prior frozen inputs/results untouched. No reseeds, new photons, sample changes, scientificgate relaxation or result selection.',
 old_initializer_run='Preserve fits-restricted/joint_minus and its failure. It changed the pre-grid initializer but actual optimizer entries differed0 for all8; it is not used as true-start evidence.',
 reuse='Nominaljoint12/9 and existing noiseless gates can be reused only after exact disabledidentity. The old nominal peak remains the prescribed input to NIR adapters; never select among localstart fits by Q.',
 environment='Always remove inherited PROSP_MINUIT_PEAK_SHIFT. Set only−2/+2 on true-start branches,0 on explicitzero identity controls, absent elsewhere. Existing Dshift cleared independently.',
 outputs='identity-gate.json; starts-gate.json; per-case gates with domain/actualstart records; full-null-identity.json; engineering-result.json only after all checks; failures and cumulative activity preserved.',
 root_release_schema={'protocol_sha256':'SHA256 execution-protocol.json','freeze_sha256':'SHA256 execution-freeze.json','stages':['identity OR starts OR nir']})
proto['fits']['starts']='True localMNPARM peak±2 around each native post-grid entry in iteration1 only, with nominal NML and shared prior unchanged. Later entries native. NIR D±.2 unchanged.'
proto['fits']['gate_order']=['previous generation/noiseless gates','disabledjoint12/9 fullscientificidentity','true local±2 actualstoredentry and numericalgates','NIR numerical/support/fullnull','engineeringeffecttable']
(H/'execution-protocol.json').write_text(json.dumps(proto,indent=2)+'\n')
files=json.loads((Q/'execution-freeze.json').read_text())['files'].copy()
for p,h in files.items():assert sha(R/p)==h,p
new=[H/x for x in ['execution-protocol.json','run_start_only.py','prepare_executor.py','freeze_execution.py','full_identity.py','schema.json','build-protocol.json','build-freeze.json','build-result.json','start-only.patch','source/snana.car','source/genmag_snoopy.c','source/prosp_minuit_start.h','build/bin/snlc_fit.exe','build/src/snana.car','build/src/genmag_snoopy.c','build/src/prosp_minuit_start.h','test-result.json','executor-unit-result.json','test_executor.py','carried-native-activity.json']]
new += [Q/x for x in ['noiseless-gate.json','identity-gate.json'] if (Q/x).exists()]
new += [Q/'native-activity.json',Q/'noisy-failure.json',Q/'run_restricted.py',P/'noiseless-checker-recovery/full_null_check.py']
for label in ['joint12','joint9']:
 new += [P/'fits-restricted'/label/x for x in ['native.log','fit.FITRES.TEXT','fit.nml','execution.json']]
for p in new:files[str(p.relative_to(R))]=sha(p)
(H/'execution-freeze.json').write_text(json.dumps({'files':files},indent=2)+'\n')
assert not (P/'fits-start-only').exists() and not (P/'derived-start-only').exists()
for v in [R/'phase2/pte/restricted-peak-engineering/start-only-hook/build',R/'phase2/pte/fit-private-lookup/SNDATA_ROOT']:assert len(str(v))<120
f=json.loads((H/'build-freeze.json').read_text())
for p,h in f['source_files'].items():assert sha(R/p)==h,p
report={'pass':True,'native_calls':0,'frozen_files':len(files),'carried_native_seconds':spent,'source_preserved_after_compile':True,'old_binary_sha256':sha(Q/'build/bin/snlc_fit.exe'),'new_binary_sha256':sha(H/'build/bin/snlc_fit.exe'),'protocol_sha256':sha(H/'execution-protocol.json'),'freeze_sha256':sha(H/'execution-freeze.json'),'runner_sha256':sha(H/'run_start_only.py')}
(H/'execution-preflight.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
