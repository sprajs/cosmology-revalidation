from pathlib import Path
import shutil,json,hashlib,datetime
P=Path(__file__).resolve().parent;R=Path('/home/szymon/Documents/ChatGPT/supernova');f=P/'runner.py';freeze=P/'fits/endpoints/freeze.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
assert not (P/'runner-v1.py').exists();assert not (P/'fits/endpoints/fit.log').exists();shutil.copy2(f,P/'runner-v1.py');shutil.copy2(freeze,P/'fits/endpoints/freeze-v1.json')
s=f.read_text();assert '(gradient<0)' in s;s=s.replace('(gradient<0)','(gradient<=0)')
old="  gate=bool(actual==expected and np.all(np.isfinite(a)) and np.all(np.isfinite(W)) and np.linalg.eigvalsh(b['C']).min()>0 and abs(cl)<=1e-8 and ob['MJDOFF']==0)"
new="""  headers={line.split(':',1)[0]:line.split(':',1)[1].strip().split()[0] for line in (R/case['file']).read_text().splitlines() if ':' in line and line.split(':',1)[1].strip() and not line.startswith('OBS:')}
  expected_peak=float(np.float32(float(headers['PEAKMJD'])))
  fixed_native=ob['shape']==1.0 and ob['AV']==0.0 and ob['peak_absolute']==expected_peak and b['entry']['shape_step']==0 and b['entry']['AV_step']==0 and b['entry']['peak_step']==0
  gate=bool(actual==expected and b['NFITDATA']==len(raw) and fixed_native and np.all(np.isfinite(a)) and np.all(np.isfinite(W)) and np.linalg.eigvalsh(b['C']).min()>0 and abs(cl)<=1e-8 and ob['MJDOFF']==0)"""
assert old in s;s=s.replace(old,new)
s=s.replace("'full_state_gate':gate,'objective':ob", "'full_state_gate':gate,'native_NFITDATA':b['NFITDATA'],'expected_input_NUMOBS':len(raw),'fixed_native_coordinates_gate':fixed_native,'expected_native_fixed_peak':expected_peak,'objective':ob")
old="  clone=True\n  if case['kind']=='clone_control':"
new="""  # FITRES does not necessarily expose NUMOBS; never invent that source column.
  raw_text=(R/case['file']).read_text().splitlines();n_input=sum(x.startswith('OBS:') for x in raw_text);header_peak=next(x.split()[1] for x in raw_text if x.startswith('PEAKMJD:'));peak_token=format(float(np.float32(float(header_peak))),'.4f')
  fixed_expected={'STRETCH':'1','AV':'0','RV':'1.51800','PKMJD':peak_token,'PKMJDINI':peak_token,'STRETCHERR':'0','AVERR':'0','RVERR':'0','PKMJDERR':'0.0000'}
  fixed_gate=all(float(fit[c][k])==float(v) for k,v in fixed_expected.items());ndof_gate=float(fit[c]['NDOF'])==n_input-1 and float(fit[c]['NDOF'])==float(published[cid]['NDOF']);numobs_present=fit[c].get('NUMOBS');numobs_gate=numobs_present is None or float(numobs_present)==n_input
  numeric=numeric and fixed_gate and ndof_gate and numobs_gate
  clone=True
  if case['kind']=='clone_control':"""
assert old in s;s=s.replace(old,new)
s=s.replace("'clone_gate':bool(clone)})", "'clone_gate':bool(clone),'FITRES_NUMOBS':numobs_present,'NUMOBS_source_note':'Absent column retained as null; native NFITDATA and reconstructed input count independently recorded.','native_NFITDATA':end['native_NFITDATA'],'expected_input_NUMOBS':n_input,'FITRES_NDOF':fit[c]['NDOF'],'source_FITRES_NDOF':published[cid]['NDOF'],'NDOF_gate':ndof_gate,'fixed_expected':fixed_expected,'fixed_actual':{k:fit[c][k] for k in fixed_expected},'fixed_gate':fixed_gate})")
f.write_text(s)
amend={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'Before any rounding fit, fix zero-gradient corner tie to high for both directions as frozen protocol requires; explicitly record and gate native NFITDATA/input row count, FITRES NDOF and fixed shape/AV/RV/peak/error fields. No domains, cohort, cap or scientific inputs changed.','original_runner_sha256':sha(P/'runner-v1.py'),'new_runner_sha256':sha(f),'original_endpoint_freeze_sha256':sha(P/'fits/endpoints/freeze-v1.json'),'input_case_ledger_sha256_unchanged':sha(P/'fits/endpoints/case-ledger.json'),'NUMOBS':'Not present in current FITRES schema: report null and exact native/input counts, never substitute a fabricated source field.'}
(P/'preexecution-amendment.json').write_text(json.dumps(amend,indent=2)+'\n');d=json.loads(freeze.read_text());d['files'][str(f.relative_to(R))]=sha(f);d['files'][str((P/'preexecution-amendment.json').relative_to(R))]=sha(P/'preexecution-amendment.json');d['preexecution_amendment_sha256']=sha(P/'preexecution-amendment.json');freeze.write_text(json.dumps(d,indent=2)+'\n');print(json.dumps({'runner_sha256':sha(f),'freeze_sha256':sha(freeze),'amendment_sha256':sha(P/'preexecution-amendment.json')},indent=2))
