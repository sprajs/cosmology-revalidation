from pathlib import Path
import hashlib,json,subprocess
P=Path(__file__).resolve().parent;R=Path('/home/szymon/Documents/ChatGPT/supernova');B=P.parent/'fit-support/build'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
sources={str(p.relative_to(R)):sha(p) for p in sorted((P/'build/src').iterdir()) if p.is_file() and p.suffix not in ['.o','.f']}
diffs=[]
for p in sorted((B/'src').iterdir()):
 if p.is_file() and p.suffix not in ['.o','.f']:
  q=P/'build/src'/p.name
  if sha(p)!=sha(q):diffs.append(p.name)
assert diffs==['genmag_snoopy.c','snlc_fit.car'],diffs
assert sha(B/'bin/snlc_fit.exe')=='d3ed187ec43968894ec651a6f94a8a9e1258058e2155dc9f523e03cf9fd56ade'
protocol={
 'status':'Build and static checks only authorized; no native fit or generation release',
 'identity':'Private separate restricted-peak estimator, not a repair of the failed unrestricted8 outcome',
 'source':'v11_04d 10ec91297e4482d593cb5d3d055b10d4aa915071 plus previously validated output/start hooks and this explicit estimator domain',
 'base_binary_sha256':sha(B/'bin/snlc_fit.exe'),
 'changed_existing_source_files':diffs,'added_source':'prosp_peak_domain.h',
 'mode':{'variable':'PROSP_HARD_PEAK_DOMAIN','absent_or_0':'Return before any native state writes, physical checks or output','1':'Metadata hard bounds and pre-physical-mean aborts','other':'Fail explicitly'},
 'domain':{'rows':117,'inputs':'All SNLC8_MJD R8 physical rows, actual fixed INIVAL zHEL, actual loaded SNooPy and KCOR phase tables; no brightness/fit/truth-driven bound',
 'formula':'max_i(MJD_i-(1+z)*tmax) <= absolute_peak <= min_i(MJD_i-(1+z)*tmin)',
 'current_absolute_bounds':[57671.34299936581,57715.067000181196],
 'parameter_coordinates':'Subtract MJDOFF; check R8 FCN and R4 mask arithmetic at both endpoints; bounded nextafter inward only if required for roundoff. Current input requires zero adjustment.',
 'invariants':'Same full row metadata/order, fixed z, MJDOFF, model/KCOR limits and native INIBND must remain exact across all iterations. Unsupported model/free-z/count aborts.',
 'initialization':'Hook after first-iteration user/SIM coordinate overrides, before FITINI_COV and FITINI_ADJUST; assertion before subsequent FITINI_PARVAL early return',
 'native_minuit':'INIBND sent unchanged to MNPARM; first-iteration initial grid remains native and must satisfy domain guards',
 'guard':'After native prior-only/out-of-bound returns, FCN checks fixed peak/redshift before its physical means. Every USRFUN additionally checks finite phase/z/shape/nuisances against loaded tables before nearest-filter/template/KCOR/extinction operations, irrespective of PROSP_FIT_SUPPORT.'},
 'unchanged':'All117 signed photons, same native noise/error pairs and HEAD files; fixed shape/AV/RV, iteratively recentered prior, native covariance loop/objective, no clipping/dropping/refilling/reseeding. Same six NIR rows per draw. This build generates no photons.',
 'gates':'All prior numerical/measurement/support/noiseless/null gates including final |peak-truth|<=4 remain unchanged; source-corrected CUTFLAG3 and FITS R8 MJD semantics remain authoritative.',
 'comparison_plan_before_noisy_results':[
 'New private binary with mode absent and then0 on exact saved noiseless joint input: require all CSP entry/row/model/W/objective tokens and FITRES science rows exactly equal saved support build; no tolerance. Existing original-support identity retained.',
 'Enable1 for same noiseless joint and fixed-peak NIR: require original mean/D/peak recovery, state/arithmetic/SPD/iteration/mask/support gates; record changed estimator coordinates rather than require enabled fit identical.',
 'Only after a separate root release: re-evaluate same8 failed-branch photon realizations under enabled estimator, preserving unrestricted outputs. All draws retained; joint9/12 and all declared starts, NIR true/fittedpeak and all D starts, full true/null model/W/objective/FITRES equality. Stop on any failure; no lowest-Q branch choice.',
 'No64 design or new draws follow automatically. No conditional response is a historical/survey correction.'
 ],
 'resources':{'build_workers':1,'build_seconds_cap':120,'new_native_fits_authorized':0,'new_photons_authorized':0},
 'schema':{
 'PROSP_PEAK_BOUNDS':['CID','iteration','n_all_metadata','zHEL','MJDOFF','phase_lower','phase_upper','KCOR_lower','KCOR_upper','raw_absolute_lower','raw_absolute_upper','safe_absolute_lower','safe_absolute_upper','parameter_lower','parameter_upper','initial_absolute_peak','initial_distance_lower_days','initial_distance_upper_days'],
 'PROSP_DOMAIN_SUMMARY':['CID','configured','bounds_assertion_count','physical_mean_calls','FCN_physical_calls','min_phase','max_phase','minimum_phase_edge_distance_rest_days','min_FCN_absolute_peak','max_FCN_absolute_peak','minimum_peak_edge_distance_observer_days'],
 'PROSP_DOMAIN_FCN_REJECT':['CID','iteration','parameter_peak','MJDOFF','zHEL'],
 'PROSP_DOMAIN_MEAN_REJECT':['CID','iteration','filter','fitdata_index','phase','Tobs','z','D','shape','AV','RV','phase_lower','phase_upper','shape_lower','shape_upper'],
 'PROSP_DOMAIN_ABORT':['CID','reason']},
 'source_and_command_hashes':sources,
 'build_command_sha256':sha(P/'build.sh'),
 'patch_sha256':sha(P/'restricted-peak.patch'),
 'static_check_sha256':sha(P/'static_check.py'),'static_result_sha256':sha(P/'static-result.json'),
 'tools':{a:subprocess.check_output([str(R/'phase2/official/build/sysroot/usr/bin/gfortran') if a=='gfortran' else a,'--version'],text=True).splitlines()[0] for a in ['gcc','gfortran','make']},
 'boundary_reporting':'Search minimum endpoint distance is not final active-bound status. Compute final margins independently from CSP_OBJECTIVE and logged bounds; do not reject an interior solution merely because an earlier trial touched an endpoint.',
 'preserved_static_failure':'First standalone -Werror build caught a logging-format field count and harness indentation; preserved under first-static-format-failure, corrected before native compile. No native fits.'
}
(P/'build-protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
files=[P/'build-protocol.json',P/'build.sh',P/'restricted-peak.patch',P/'static_check.py',P/'static-result.json',P/'prepare.py',P/'freeze_build.py']
freeze={'files':{str(p.relative_to(R)):sha(p) for p in files},'source_files':sources}
(P/'build-freeze.json').write_text(json.dumps(freeze,indent=2)+'\n')
print('protocol',sha(P/'build-protocol.json'));print('freeze',sha(P/'build-freeze.json'));print('source_count',len(sources))
