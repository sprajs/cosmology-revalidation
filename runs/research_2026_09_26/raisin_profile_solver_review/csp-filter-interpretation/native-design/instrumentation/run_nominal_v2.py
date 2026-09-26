from pathlib import Path
import json,hashlib,shutil,os,subprocess,time
ROOT=Path.cwd();O=Path(__file__).resolve().parent;B=O/'SNANA-v11_04k-output';P=O/'pilot-v2';S=ROOT/'runs/research_2026_09_26/csp_native_filter_response/cohort-preparation/pilot';R=ROOT/'runs/research_2026_09_26/csp_native_filter_response/fits/nominal';sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert not P.exists();d=P/'fits/nominal';d.mkdir(parents=True);(P/'data').symlink_to(S/'data',target_is_directory=True)
for name in ['fit.nml','vpec.list']:shutil.copyfile(S/'fits/nominal'/name,d/name)
inputs=[B/'bin/snlc_fit.exe',B/'src/snlc_fit.car',O/'build-amendment-v2.json',O/'output-only-v2.patch',d/'fit.nml',d/'vpec.list',Path(__file__),R/'fit.FITRES.TEXT',R/'fit.LCPLOT.TEXT']+list((S/'data/nominal/CSPDR3_RAISIN').glob('*'))
proto={'scope':'Output-only nominal equivalence, no changed-filter arm. Exact prepared NML and vpec copied, data symlink to same prepared pilot nominal data.','hashes':{str(p.relative_to(ROOT)):sha(p) for p in inputs if p.is_file()},'timeout_seconds':30,'schemas':{'CSP_ENTRY':['CID','ITER','NFITDATA','USE_FITCOV','LREPEAT_ITER','D_entry','shape_entry','AV_entry','peak_entry_absolute','D_step','shape_step','AV_step','peak_step','D_previous','shape_previous','AV_previous','peak_previous_absolute'],'CSP_ROW':['CID','ITER','accepted_index','source_epoch','band','MJD','Trest','modelF','model_mag_error','dataF','data_error','z','MWEBV','rest_filter_mean_wavelength','fudge_flux_error','model_flux_error','rest_lambda_fit_min','rest_lambda_fit_max'],'CSP_OBJECTIVE':['CID','ITER','NFITDATA','USE_FITCOV','totalQ','priorQ','sigmaQ','D','shape','AV','peak_absolute','peak_prior_center_absolute','search_peak','MJDOFF'],'CSP_WROW':['CID','ITER','row','inverse_covariance_row...'],'CSP_WDIAG':['CID','ITER','row','inverse_variance']},'blocks':'Rows followed by OBJECTIVE and W. Preserve every block and all entry markers.'}
(P/'protocol.json').write_text(json.dumps(proto,indent=2)+'\n');print('frozen',sha(P/'protocol.json'),flush=True)
env=os.environ.copy();env.update(SNANA_DIR='/tmp/snana-csp-audit-8c5f0d',SNDATA_ROOT=str(ROOT/'phase2/official/inputs/SNDATA_ROOT'),LD_LIBRARY_PATH=str(ROOT/'phase2/official/build/sysroot/usr/lib'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
t=time.monotonic()
with (d/'fit.log').open('x') as f:p=subprocess.run([str(B/'bin/snlc_fit.exe'),'fit.nml'],cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=30)
r={'returncode':p.returncode,'elapsed_seconds':time.monotonic()-t,'protocol_sha256':sha(P/'protocol.json'),'log_sha256':sha(d/'fit.log'),'science_rows_equal':False,'LCPLOT_bytes_equal':False}
if p.returncode==0:
 rows=lambda p:[x for x in p.read_text().splitlines() if x.startswith(('SN:','VARNAMES:'))]
 r['science_rows_equal']=rows(d/'fit.FITRES.TEXT')==rows(R/'fit.FITRES.TEXT');r['LCPLOT_bytes_equal']=(d/'fit.LCPLOT.TEXT').read_bytes()==(R/'fit.LCPLOT.TEXT').read_bytes()
r['pass']=r['returncode']==0 and r['science_rows_equal'] and r['LCPLOT_bytes_equal'];(P/'execution-result.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
