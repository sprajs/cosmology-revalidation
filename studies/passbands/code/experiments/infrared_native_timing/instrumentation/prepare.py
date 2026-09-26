from pathlib import Path
import hashlib,json,shutil,difflib
R=Path('/home/szymon/Documents/ChatGPT/supernova')
P=R/'runs/research_2026_09_26/raisin_profile_solver_review/csp-filter-interpretation/native-design/timing-followon-design/prospective-engineering'
I=P/'instrumentation';I.mkdir(exist_ok=False)
B=P/'build-v11_04d';S=B/'src/snlc_sim.c';old=S.read_text()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
shutil.copy2(S,I/'snlc_sim.original.c');shutil.copy2(B/'bin/snlc_sim.exe',I/'snlc_sim.original.exe')
schema={}
def printf(tag, cols, expr, fmts=None):
    schema[tag]=cols
    if fmts is None:fmts=['%.17g']*len(cols)
    return 'printf("'+tag+' '+' '.join(fmts)+'\\n", '+', '.join(expr)+');'
noise_cols=['attempt','cid','epoch','band','field','mjd','trest','flux_true_pe','flux_observed_pe','src_variance_pe2','sky_variance_pe2','zp_variance_pe2','template_source_pe','template_sky_variance_pe2','host_variance_pe2','reported_variance_before_realization_pe2','reported_variance_adjustment_pe2','reported_variance_final_pe2','calc_reported_variance_pe2','true_S_variance_pe2','true_SZ_variance_pe2','true_T_variance_pe2','true_F_variance_pe2','true_SUM_variance_pe2','gauss_S','gauss_T','gauss_F','shift_SZ_pe','shift_T_pe','shift_F_pe','Npe_per_fluxcal','ADU_per_Npe','flux_ADU','error_ADU','NEA','saturation_excess_pe','random_template_option','smearflag_flux','redcov_index']
noise_expr=['NGENLC_TOT','GENLC.CID','epoch','BAND','FIELD','SIMLIB_OBS_GEN.MJD[epoch]','GENLC.epoch_rest[epoch]','fluxTrue','fluxObs','FLUXNOISE->SQSIG_SRC','FLUXNOISE->SQSIG_SKY','FLUXNOISE->SQSIG_ZP','FLUXNOISE->SQSIG_TSRC','FLUXNOISE->SQSIG_TSKY','FLUXNOISE->SQSIG_HOST_PHOT','SIG_FINAL_TRUE*SIG_FINAL_TRUE','FLUXNOISE->SQSIG_RAN','FLUXNOISE->SQSIG_FINAL_DATA','FLUXNOISE->SQSIG_CALC_DATA']+['FLUXNOISE->SQSIG_FINAL_TRUE[TYPE_FLUXNOISE_'+x+']' for x in ['S','SZ','T','F','SUM']]+['GAURAN_SZ','GAURAN_T','GAURAN_F','SHIFT_SZ','SHIFT_T','SHIFT_F','FLUXNOISE->Npe_over_FLUXCAL','FLUXNOISE->NADU_over_Npe','GENLC.flux[epoch]','GENLC.fluxerr_data[epoch]','FLUXNOISE->NEA','npe_above_sat','OVP','INPUTS.SMEARFLAG_FLUX','FLUXNOISE->INDEX_REDCOV']
fmt=['%d','%d','%d','%s','%s']+['%.17g']*(len(noise_cols)-9)+['%d']*4
assert len(noise_cols)==len(noise_expr)==len(fmt)
noise=printf('PROSP_NOISE',noise_cols,noise_expr,fmt)
anchor='  // check for really crazy flux values\n  check_crazyFlux(epoch, FLUXNOISE);'
assert old.count(anchor)==1
new=old.replace(anchor,anchor+'\n\n  /* Output-only prospective noise ledger: no model/RNG calls. */\n  if (getenv("PROSP_LEDGER") != NULL) {\n    '+noise+'\n  }')
anchor='  SNDATA.SEARCH_PEAKMJD = GENLC.PEAKMJD_SMEAR;'
eventcols=['attempt','cid','libid','peak_true','peak_header','zhel_true','zcmb_true','DLMU_true','shape','AV','RV','MWEBV_map','MWEBV_true','MWEBV_error','nepoch','forced_accept']
eventexpr=['NGENLC_TOT','GENLC.CID','GENLC.SIMLIB_ID','GENLC.PEAKMJD','GENLC.PEAKMJD_SMEAR','GENLC.REDSHIFT_HELIO','GENLC.REDSHIFT_CMB','GENLC.DLMU','GENLC.SHAPEPAR','GENLC.AV','GENLC.RV','GENLC.MWEBV','GENLC.MWEBV_SMEAR','GENLC.MWEBV_ERR','GENLC.NEPOCH','GENLC.ACCEPTFLAG_FORCE']
event=printf('PROSP_EVENT',eventcols,eventexpr,['%d']*3+['%.17g']*11+['%d']*2)
rowcols=['attempt','cid','epoch','mjd','filter_index','fluxcal_native_R4','fluxcal_error_native_R4','photflag','obsflag_gen','obsflag_write','obsflag_peak','obsflag_template','source_simlib_index','source_exposure_id']
rowexpr=['NGENLC_TOT','GENLC.CID','prosp_ep','SNDATA.MJD[prosp_ep]','GENLC.IFILT_OBS[prosp_ep]','(double)SNDATA.FLUXCAL[prosp_ep]','(double)SNDATA.FLUXCAL_ERRTOT[prosp_ep]','SNDATA.PHOTFLAG[prosp_ep]','GENLC.OBSFLAG_GEN[prosp_ep]','GENLC.OBSFLAG_WRITE[prosp_ep]','GENLC.OBSFLAG_PEAK[prosp_ep]','GENLC.OBSFLAG_TEMPLATE[prosp_ep]','SIMLIB_OBS_GEN.ISTORE_RAW[prosp_ep]','SIMLIB_OBS_GEN.IDEXPT[prosp_ep]']
row=printf('PROSP_ROW',rowcols,rowexpr,['%d']*3+['%.17g','%d','%.17g','%.17g']+['%d']*7)
assert new.count(anchor)==1
new=new.replace(anchor,anchor+'\n\n  if (getenv("PROSP_LEDGER") != NULL) {\n    '+event+'\n    int prosp_ep;\n    for(prosp_ep=1; prosp_ep<=GENLC.NEPOCH; prosp_ep++) {\n      '+row+'\n    }\n  }')
anchor='    NGENLC_TOT++;'
assert new.count(anchor)==1
new=new.replace(anchor,anchor+'\n    if(getenv("PROSP_LEDGER") != NULL) { printf("PROSP_ATTEMPT %d %d\\n", NGENLC_TOT, ilc); }')
schema['PROSP_ATTEMPT']=['actual_attempt','loop_counter']
anchor='  ilc = ilc_orig = *ILC;'
# source actual spacing
if anchor not in new:anchor='  ilc = ilc_orig = *ILC ;'
assert new.count(anchor)==1
new=new.replace(anchor,anchor+'\n  if(getenv("PROSP_LEDGER") != NULL) { printf("PROSP_REJECT %d %d %d %s\\n", NGENLC_TOT, *ILC, GENLC.CID, REJECT_STAGE); }')
schema['PROSP_REJECT']=['actual_attempt','loop_counter','cid','stage']
S.write_text(new);(I/'snlc_sim.instrumented.c').write_text(new)
(I/'output-only.patch').write_text(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='original/snlc_sim.c',tofile='instrumented/snlc_sim.c')))
(I/'schema.json').write_text(json.dumps(schema,indent=2)+'\n')
# Short, project-local entry point avoids both C and Fortran path buffers.
alias=R/'phase2/pte';assert not alias.exists();alias.symlink_to(P, target_is_directory=True)
ledger={'status':'frozen before instrumented compile or generation','purpose':'output-only binary64 internal truth/noise and native R4 flux/error ledger; no model/RNG/noise/fit change','original_generator_sha256':sha(I/'snlc_sim.original.exe'),'original_source_sha256':sha(I/'snlc_sim.original.c'),'instrumented_source_sha256':sha(S),'patch_sha256':sha(I/'output-only.patch'),'schema_sha256':sha(I/'schema.json'),'source_commit':'10ec91297e4482d593cb5d3d055b10d4aa915071','short_alias':str(alias),'canonical_target':str(P),'SNANA_DIR':str(alias/'build-v11_04d'),'generation_authorized':False,'build_cap_seconds':180,'prior_build_seconds':18.751008435006952,'remaining_build_seconds':161.24899156499305,'equivalence_gate':'same fixed 8 attempts/seed with original and instrumented binaries, compare all native FITS columns exactly excluding run/provenance filename header cards; ledger-enabled must not change source outputs','guard':'PROSP_LEDGER environment presence, absent default emits no added output'}
(I/'protocol.json').write_text(json.dumps(ledger,indent=2)+'\n')
print(json.dumps({'instrumentation_protocol_sha256':sha(I/'protocol.json'),'source_sha256':sha(S),'alias':str(alias),'SNANA_DIR_chars':len(str(alias/'build-v11_04d'))},indent=2))
