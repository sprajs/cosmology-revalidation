"""Metadata/source audit only: never inspect variant observed-flux residuals."""
from pathlib import Path
import json,hashlib,difflib
import numpy as np,pandas as pd
from astropy.io import fits
ROOT=Path(__file__).resolve().parents[3];OUT=Path(__file__).resolve().parent/'noise_variant_design';OUT.mkdir(exist_ok=True);BASE=ROOT/'phase2/literature/simulations';FWD=ROOT/'phase2/hierarchy/forward'
ARMS=['P21','P21_rho000','P21_rho090','P21_noisetrue120'];sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();hashes={};tables={};attempts={};fitsmeta={};heads={};phot={};records=[]
cols=['generated_attempt_index','CID','GENZ','GALZTRUE','LIBID','RA','DEC','MWEBV','MU','PEAKMJD','GALID','LOGMASS_TRUE','LOGSFR_TRUE','r_obs_auto','obs_gr_auto','SALT2mB','SALT2x1','SALT2c','SALT2alpha','SALT2beta','AV','RV','MAGSMEAR_COH','PEAKMAG_g','PEAKMAG_r','PEAKMAG_i','PEAKMAG_z','NOBS','TRESTMIN','TRESTMAX','TGAPMAX']
for arm in ARMS:
 version='PH2_pilot02_'+arm;mp=BASE/'manifests'/f'{version}.json';m=json.loads(mp.read_text());assert m['exit_code']==0 and m['success_marker'];hnum=0
 for item in [m['input'],m['executable'],m['log'],*m['outputs']]:
  p=ROOT/item['path'];h=sha(p);assert h==item['sha256'],p;hashes[str(p.relative_to(ROOT))]=h;hnum+=1
 for p in [BASE/f'{version}-generated.csv.gz',FWD/f'{arm}-attempts.csv.gz',FWD/f'{arm}-fitted.csv.gz',mp]:hashes[str(p.relative_to(ROOT))]=sha(p)
 g=pd.read_csv(BASE/f'{version}-generated.csv.gz');assert len(g)==26518 and g.generated_attempt_index.is_unique;tables[arm]=g
 a=pd.read_csv(FWD/f'{arm}-attempts.csv.gz');assert np.array_equal(a.generated_attempt_index,g.generated_attempt_index);assert np.array_equal(a.CID,g.CID);attempts[arm]=a
 fm=pd.read_csv(FWD/f'{arm}-fitted.csv.gz');assert fm.generated_attempt_index.is_unique;fitsmeta[arm]=fm
 stem=BASE/'outputs'/version/version;head=fits.getdata(str(stem)+'_HEAD.FITS',1);ph=fits.getdata(str(stem)+'_PHOT.FITS',1);heads[arm]={int(row['SNID']):row for row in head};phot[arm]=ph
 assert set(heads[arm])==set(g.loc[(g.SIM_EFFMASK==5)&(g.CUTMASK==4095),'CID'])
 assert set(fm.generated_attempt_index)==set(a.loc[a.fit_pass,'generated_attempt_index'])
 eq={c:bool(np.array_equal(tables['P21'][c],g[c])) for c in cols};assert all(eq.values())
 records.append(dict(arm=arm,verified_manifest_entries=hnum,generated=len(g),written=int(a.written.sum()),fitted=int(a.fit_pass.sum()),quality=int(a.basic_quality_pass.sum()),all_declared_latent_design_columns_equal=all(eq.values())))
# Compare every common-written epoch in recorded order. True SIM_MAGOBS only;
# do not read flux, errors, PHOTFLAG or fit residuals to select the next cohort.
columns=['MJD','BAND','CCDNUM','FIELD','PSF_SIG1','PSF_SIG2','PSF_RATIO','SKY_SIG','SKY_SIG_T','RDNOISE','ZEROPT','ZEROPT_ERR','GAIN','XPIX','YPIX','SIM_MAGOBS'];phot_checks=[]
for arm in ARMS[1:]:
 common=sorted(set(heads['P21'])&set(heads[arm]));n=0;fail=[]
 for cid in common:
  hb=heads['P21'][cid];ha=heads[arm][cid];b=phot['P21'][int(hb['PTROBS_MIN'])-1:int(hb['PTROBS_MAX'])];a=phot[arm][int(ha['PTROBS_MIN'])-1:int(ha['PTROBS_MAX'])];n+=len(b)
  if len(a)!=len(b):fail.append((cid,'NOBS'));continue
  for c in columns:
   if not np.array_equal(a[c],b[c]):fail.append((cid,c))
 assert not fail,(arm,fail[:10]);phot_checks.append(dict(arm=arm,common_written=len(common),epochs_checked=n,columns_checked=columns,failures=fail))
# Preserve all-attempt paired selection changes independently of pilot eligibility.
ledger=attempts['P21'][['generated_attempt_index','CID','GENZ','LIBID']].copy()
for arm in ARMS:
 for stage in ['written','fit_pass','basic_quality_pass']:ledger[arm+'_'+stage]=attempts[arm][stage].to_numpy()
ledger.to_csv(OUT/'full-selection-ledger.csv',index=False)
flips=[]
for arm in ARMS[1:]:
 for stage in ['written','fit_pass','basic_quality_pass']:
  b=attempts['P21'][stage].to_numpy(bool);v=attempts[arm][stage].to_numpy(bool);flips.append(dict(arm=arm,stage=stage,both=int((b&v).sum()),newly_selected=int((~b&v).sum()),lost=int((b&~v).sum()),net=int(v.sum()-b.sum())))
pd.DataFrame(flips).to_csv(OUT/'selection-flips.csv',index=False)
# Freeze pilot by common archived-fit metadata, not by classifier or residual.
common=set.intersection(*[set(fitsmeta[a].generated_attempt_index) for a in ARMS]);eligible=sorted(common)
for arm in ARMS:
 q=fitsmeta[arm].set_index('generated_attempt_index').loc[eligible];good=np.isfinite(q[['x0','x1_double','c_double','t0_double','zHEL']]).all(axis=1)&(q.x0>0)&q.zHEL.between(.05,1.2);eligible=sorted(set(eligible)&set(q.index[good]))
q=fitsmeta['P21'].set_index('generated_attempt_index').loc[eligible].copy();q['field']=[str(phot['P21'][int(heads['P21'][int(cid)]['PTROBS_MIN'])-1]['FIELD']).strip() for cid in q.CID];q['zbin']=pd.cut(q.zHEL,[.05,.2,.35,.5,.65,.8,1.2],include_lowest=True,labels=False).astype(int);q['cell']=q.field+'_'+q.zbin.astype(str);q['hash_rank']=[hashlib.sha256(f'20260926-noise-variant-pilot:{int(i)}'.encode()).hexdigest() for i in q.index]
first=q.sort_values('hash_rank').groupby('cell',sort=True).head(4);chosen=set(first.index);assert len(chosen)<=256
for i in q.sort_values('hash_rank').index:
 if len(chosen)==256:break
 chosen.add(i)
assert len(chosen)==256;chosen=sorted(chosen);cohort=q.loc[chosen].reset_index();cohort.to_csv(OUT/'cohort.csv',index=False,float_format='%.17g')
for arm in ARMS:
 out=fitsmeta[arm].set_index('generated_attempt_index').loc[chosen].reset_index();out=out.merge(cohort[['generated_attempt_index','field','zbin','cell','hash_rank']],on='generated_attempt_index',validate='one_to_one');out.to_csv(OUT/f'{arm}-cohort.csv',index=False,float_format='%.17g');(OUT/f'{arm}-cids.txt').write_text('\n'.join(out.CID.astype(str))+'\n')
inputs={arm:(BASE/'inputs'/f'PH2_pilot02_{arm}.input').read_text() for arm in ARMS};diff=''
for arm in ARMS[1:]:diff+=''.join(difflib.unified_diff(inputs['P21'].splitlines(True),inputs[arm].splitlines(True),fromfile='P21',tofile=arm))
(OUT/'input-differences.patch').write_text(diff)
for p in [Path(__file__),ROOT/'phase2/official/build/SNANA-2fe0f56/src/snlc_sim.c',ROOT/'phase2/official/build/SNANA-2fe0f56/src/sntools_fluxErrModels.c',BASE/'inputs/DES_FLUXERRMODEL_TRUE120_DATA100.DAT',ROOT/'phase2/official/inputs/SNDATA_ROOT/simlib/DES/DES-SN5YR_DES_FLUXERRMODEL_SIM.DAT']:hashes[str(p.relative_to(ROOT))]=sha(p)
report=dict(arms=records,phot_checks=phot_checks,common_archived_fits=len(common),eligible_common_fits=len(eligible),cohort_count=len(cohort),cohort_cells=int(cohort.cell.nunique()),fields=sorted(cohort.field.unique()),rng_boundary='Same verified binary/seed/attempt design. Source draws flux-noise variates before deterministic rho mixing, which consumes no RNG. Stored latent/cadence/noiseless epoch equality verified; raw random arrays were not saved, so no direct bitwise variate-array comparison is claimed.',cohort_selection='Four smallest fixed hashes perfield/z cell then fill to256 by same hash. Common archived fitted finite-coordinate support in all4 arms, before classifier/native residual scores. This conditions on all4 historical fit outcomes, not a survey causal sample.')
(OUT/'inventory.json').write_text(json.dumps(report,indent=2)+'\n');(OUT/'manifest.json').write_text(json.dumps({'inputs_sha256':hashes,'outputs_sha256':{p.name:sha(p) for p in OUT.iterdir() if p.is_file() and p.name!='manifest.json'}},indent=2)+'\n');print(json.dumps(report,indent=2))
