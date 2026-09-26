from pathlib import Path
import json,hashlib,ast,re,collections
from astropy.io import fits
Q=Path(__file__).resolve().parent;P=Q.parent;R=Path('/home/szymon/Documents/ChatGPT/supernova')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
old=json.loads((P/'protocol-v5.json').read_text())
proto={k:old[k] for k in ['question','scope','truth','cadence','fits','gates','estimands']}
proto.update(status='Frozen new-estimator stages; each awaits an explicit root release; no new photons',
 identity='Private restricted peak estimator conditional on all117 exposure times and fixed redshift; original unrestricted8 support failure remains failed',
 parent_protocol_sha256=sha(P/'protocol-v5.json'),corrected_native_semantics='CUTFLAG_SNANA3 are two pass bits with ERRFLAG_FIT0; FITS MJD remains R8. Corrections already source-verified in noiseless-checker-recovery, not new scientific tolerance changes.',
 domain_build_protocol_sha256=sha(Q/'build-protocol.json'),binary_sha256=sha(Q/'build/bin/snlc_fit.exe'),
 stages={'disabled':'New binary absent and explicit0 on unchanged saved noiseless joint input; exact full CSP/model/W/objectives/FITRES/support/HEAD/PHOT identity versus fits-readme/noiseless_joint. Stop if either fails.',
 'noiseless':'Active1 same noiseless joint12 then NIR12; original recovery, complete mask, covariance/stationarity and source-support gates.',
 'noisy':'Active1 exact existing8 photons, joint12/9 and peak starts±2; then NIRtruth/estimated each12/9 and postinit Dstarts±.2; NIR null12. Full saved-model/weight/objective/FITRES null identity and domain-schema checks precede engineering effect table.'},
 mutation_scope='Only new fits-restricted/, derived-restricted/ and this owned directory; all earlier runs, native inputs and binaries unchanged. Estimated and null HEAD copies alter only PEAKMJD; physical PHOT stays same file.',
 domain='All117 SNLC8_MJD, loaded SNooPy[-20,70] and KCOR[-20,85], fixed native zHEL; actual peak bounds[57671.34299936581,57715.067000181196]. They are physical table limits, not a truth-centered timing prior or empirical validation. All phases before any physical mean must pass. Final ±4days remains an independent engineering recovery guard.',
 domain_gates='Exact logged raw/safe bounds at each iteration, configured count and callback iteration coverage, finite support, same bound coordinates, no abort/reject/unknown record, equal hard-guard and independent support-audit mean-call counts. Every final CSP peak inside bounds; report final margins separately from trial minimum.',
 failed_history='Original first8 noisy joint evaluated unsupported means for CIDs4/7 and remains closed; same8 reused only to validate a separately specified estimator, not as new independent MonteCarlo draws. Historical LCPLOT timing gate also remains failed.',
 resources={'workers':1,'new_estimator_cumulative_fit_seconds':120,'per_process_seconds':40,'planned_new_fit_processes':17,'disabled':2,'noiseless':2,'noisy':13,'previous_failed_estimator_seconds':'Recorded separately in native-activity-v5.json; not erased and not charged/reused in this explicitly new120s allowance.','build_seconds_separate':json.loads((Q/'build-result.json').read_text())['wall_seconds'],'generation_calls':0},
 failure='Preserve all attempted outputs and new activity ledger; no overwritten cases, refill, reseed, CID/epoch deletion, relaxed gate, automatic bound extension or lowest-Q branch selection. A failed stage stops later stages. No64 draw protocol is frozen or authorized.',
 outputs='Per-case gates with final/search margins; disabled exact identity reports; active noiseless gate; exact full null report; engineering-result.json only after every noisy gate passes.',
 root_release_schema={'protocol_sha256':'SHA256 execution-protocol.json','freeze_sha256':'SHA256 execution-freeze.json','stages':['disabled OR noiseless OR noisy']})
proto['cadence']['fixedmask']='Exact117 joint and6NIR at every callback. One metadata-derived fixed peak domain for all117 applies in every arm; no flux/SNR cuts or epoch expansion.'
proto['fits']['peak_estimator']='Same native free D+peak with hard metadata-domain limits, fixed shape/AV/RV, same-photon grizJH. Native5day prior remains iteratively recentered. No prior-free/external timing claim.'
proto['fits']['gate_order']=['prior generation gate already passed','disabled binary exact full-state identity','active noiseless recovery','active eight joint numerical/support gates','active NIR numerical/support/full-null gates','engineering effect table']
proto['gates']['native_flags']='Complete FITRES membership, ERRFLAG_FIT0 and CUTFLAG_SNANA3 pass bits.'
proto['gates']['null']='All CSP entry/mean/measurement/error/phase/W/objective tokens, FITRES science rows, support and domain records exact; all HEAD table fields equal and identical PHOT file.'
(Q/'execution-protocol.json').write_text(json.dumps(proto,indent=2)+'\n')
# Freeze only existing authoritative source/data. Future products cannot be included yet.
files=[Q/x for x in ['execution-protocol.json','run_restricted.py','prepare_executor.py','freeze_execution.py','build-protocol.json','build-freeze.json','build-result.json','restricted-peak.patch','build/bin/snlc_fit.exe','build/src/snlc_fit.car','build/src/genmag_snoopy.c','build/src/prosp_peak_domain.h']]
files += [P/x for x in ['generation-gate.json','protocol-v5.json','run_engineering_v5_checker_resume.py','noiseless-checker-recovery/full_null_check.py','inputs-v2/fit-template-v3.nml','inputs-v2/cadence.simlib','inputs-v2/cadence-ledger.csv','inputs-v2/manifest.json','fit-private-lookup/SNDATA_ROOT/SIM/PATH_SNDATA_SIM.LIST','fits-readme/noiseless_joint/fit.nml','fits-readme/noiseless_joint/native.log','fits-readme/noiseless_joint/fit.FITRES.TEXT']]
files += [R/'runs/research_2026_09_26/astra_design/raisin_timing_assets/instrumentation_2021/parse_native.py']
for ds in ['ledger','noiseless']:
 d=P/'readme-adapted'/ds/'PTE';files += [p for p in d.iterdir() if p.is_file()]
 files.append(P/'generation-v5'/ds/'native.log')
for x in [P/'inputs-v2/kcor.fits',P/'inputs-v2/snoopy.B18']:
 files += [x] if x.is_file() else [p for p in x.rglob('*') if p.is_file()]
# Full preserved SNDATA reference inventory used by previous private-root adapter.
root=P/'fit-private-lookup/SNDATA_ROOT'
for folder in ['filters','standards']:
 d=root/folder
 if d.exists():files += [p for p in d.rglob('*') if p.is_file()]
freeze={'files':{str(p.relative_to(R)):sha(p) for p in sorted(set(files))}}
(Q/'execution-freeze.json').write_text(json.dumps(freeze,indent=2)+'\n')
# Read-only preflight: no native calls.
checks={}
for ds in ['ledger','noiseless']:
 d=P/'readme-adapted'/ds/'PTE';h=list(d.glob('*HEAD.FITS*'));p=list(d.glob('*PHOT.FITS*'));assert len(h)==len(p)==1
 with fits.open(h[0]) as hf,fits.open(p[0]) as pf:
  assert 'REDSHIFT_HELIO' in hf[1].columns.names
  count=8 if ds=='ledger' else 1;assert len(hf[1].data)==count
  for row in hf[1].data:
   qq=pf[1].data[int(row['PTROBS_MIN'])-1:int(row['PTROBS_MAX'])];assert len(qq)==117
  checks[ds]={'events':count,'FITS_MJD_format':pf[1].columns['MJD'].format,'PHOT_sha256':sha(p[0]),'HEAD_sha256':sha(h[0])}
assert (root/'SIM/PATH_SNDATA_SIM.LIST').read_text().strip()==''
for x in [R/'phase2/pte/restricted-peak-engineering/build',R/'phase2/pte/fit-private-lookup/SNDATA_ROOT']:assert len(str(x))<120
# No path adaptation needed: new workdirs retain the old root-level depth.
for data in ['../../readme-adapted/noiseless','../../readme-adapted/ledger']:
 hypothetical=P/'fits-restricted'/'not-created';assert (hypothetical/data).resolve().exists()
assert not (P/'fits-restricted').exists() and not (P/'derived-restricted').exists(),'fresh outputs required'
s=ast.parse((Q/'run_restricted.py').read_text());assert len([n for n in s.body if isinstance(n,ast.FunctionDef) and n.name=='native'])==1
report={'pass':True,'native_calls':0,'preflight':checks,'frozen_files':len(freeze['files']),'protocol_sha256':sha(Q/'execution-protocol.json'),'freeze_sha256':sha(Q/'execution-freeze.json'),'runner_sha256':sha(Q/'run_restricted.py')}
(Q/'execution-preflight.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
