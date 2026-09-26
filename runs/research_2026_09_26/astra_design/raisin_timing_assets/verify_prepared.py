from pathlib import Path
import json,gzip,hashlib,csv,struct,collections
R=Path.cwd();P=Path(__file__).parent;O=P/'baseline_2021';d=json.loads((O/'baseline-freeze.json').read_text())
for n,h in d['files'].items():assert hashlib.sha256((R/n).read_bytes()).hexdigest()==h
old=json.loads((P/'public_metadata/author-20211111-tree.json').read_text());by={r['path']:r for r in old['tree']}
assets={'kcor/kcor_DES_NIR.fits':R/'sources/repos/djones1040__RAISIN_DataRelease@a383c4b/kcor/kcor_DES_NIR.fits','snoopy.B18/SNooPy_test.fits':R/'sources/repos/djones1040__RAISIN_DataRelease@a383c4b/model/snoopy.B18/SNooPy_B18.fits','snoopy.B18/snoopy.info':R/'phase2/official/inputs/SNDATA_ROOT/models/snoopy/snoopy.B18/snoopy.info','fit/sim/REFAC_DES_RAISIN_NIR.nml':R/'runs/research_2026_09_26/raisin_simulation_assets/timing/REFAC_DES_RAISIN_NIR.nml'}
for n,p in assets.items():b=p.read_bytes();assert len(b)==by[n]['size'];assert hashlib.sha1(f'blob {len(b)}\0'.encode()+b).hexdigest()==by[n]['sha']
with gzip.open(R/'runs/research_2026_09_26/raisin_simulation_execution_review/nominal_fitres/lcplot-feasibility/FITOPT000.LCPLOT.gz','rt') as f:
 rows=collections.defaultdict(list)
 for l in f:
  a=l.split()
  if not a:continue
  if a[0]=='VARNAMES:':names=a[1:]
  if a[0]=='OBS:':
   r=dict(zip(names,a[1:]))
   if r['DATAFLAG']=='1' and int(r['CID'])<=8:rows[r['CID']].append(r)
ledger=[]
for cid in map(str,range(1,9)):
 a=O/'fits/baseline/data/DES_RAISIN_SIM'/f'{cid}.DAT';b=O/'fits/baseline_copy/data/DES_RAISIN_SIM'/f'{cid}.DAT';assert a.read_bytes()==b.read_bytes();lines=a.read_text().splitlines();obs=[l.split()[1:] for l in lines if l.startswith('OBS:')];assert len(obs)==len(rows[cid]);assert obs==[[q['MJD'],q['BAND'],'NULL',q['FLUXCAL'],q['FLUXCAL_ERR']] for q in rows[cid]]
 pk=next(l.split()[1] for l in lines if l.startswith('PEAKMJD:'));pk32=struct.unpack('f',struct.pack('f',float(pk)))[0];ledger.append({'CID':cid,'rows':len(obs),'peak_token':pk,'native_header_float32':pk32,'float32_minus_token_day':pk32-float(pk)})
with (P/'prepared-float32-ledger.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(ledger[0]));w.writeheader();w.writerows(ledger)
assert not any((O/'fits').rglob('fit.log'))
res={'pass':True,'frozen_file_hashes_verified':len(d['files']),'contemporaneous_model_KCOR_info_NML_Git_blobs_verified':4,'byteidentical_pairs':8,'source_epoch_tokens_unchanged':sum(r['rows'] for r in ledger),'no_native_fit_logs':True,'protocol_sha256':hashlib.sha256((O/'protocol.json').read_bytes()).hexdigest(),'freeze_sha256':hashlib.sha256((O/'baseline-freeze.json').read_bytes()).hexdigest()};(P/'prepared-verification.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
