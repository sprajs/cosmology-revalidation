from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
from astropy.io import fits
P=Path(__file__).resolve().parent;ROOT=P.parents[2];O=P/'simulation_design';assert not O.exists();O.mkdir()
inputs=[];summary={};edges=[.05,.2,.35,.5,.65,.8,1.2]
for arm in ['P21','G10']:
 fpath=ROOT/f'phase2/hierarchy/forward/{arm}-fitted.csv.gz';apath=ROOT/f'phase2/hierarchy/forward/{arm}-attempts.csv.gz';f=pd.read_csv(fpath,dtype={'CID':str});att=pd.read_csv(apath)
 version='PH2_pilot02_'+arm;d=ROOT/'phase2/literature/simulations/outputs'/version;hp=d/(version+'_HEAD.FITS');pp=d/(version+'_PHOT.FITS');h=fits.getdata(hp,1);phot=fits.getdata(pp,1)
 field={str(q['SNID']).strip():str(phot[int(q['PTROBS_MIN'])-1]['FIELD']).strip() for q in h};f['field_first_phot']=f.CID.map(field);assert f.field_first_phot.notna().all();f['field']=f.field_first_phot
 eligible=f[np.isfinite(f[['x0','x1','c','PKMJD','zHEL']]).all(axis=1)&(f.x0>0)&f.zHEL.between(.05,1.2)].copy();eligible['zbin']=pd.cut(eligible.zHEL,edges,include_lowest=True,labels=False)
 assert eligible.zbin.notna().all();eligible['cell']=eligible.field+'_'+eligible.zbin.astype(int).astype(str);eligible['selection_hash']=eligible.CID.map(lambda cid:hashlib.sha256(f'20260926-sim-residual-v1|{arm}|{cid}'.encode()).hexdigest())
 ranked=eligible.sort_values('selection_hash');first=ranked.groupby('cell',sort=False).head(4);assert len(first)<=256
 chosen=pd.concat([first,ranked[~ranked.CID.isin(first.CID)].head(256-len(first))]).sort_values(['field','zHEL','CID']).reset_index(drop=True);assert len(chosen)==256 and chosen.CID.is_unique
 chosen.to_csv(O/f'{arm}-cohort.csv',index=False);(O/f'{arm}-cids.txt').write_text('\n'.join(chosen.CID)+'\n')
 summary[arm]={'generated':len(att),'written':int(att.written.sum()),'archived_fit':int(att.fit_pass.sum()),'archived_basic_quality':int(att.basic_quality_pass.sum()),'eligible_archived_fits':len(eligible),'pilot':len(chosen),'pilot_archived_basic_quality':int(chosen.basic_quality_pass.sum()),'pilot_field_counts':chosen.field.value_counts().to_dict(),'pilot_cell_counts':chosen.cell.value_counts().to_dict(),'full_fit_firstPHOT_vs_FITRES_field_disagreement':int((f.field_first_phot!=f.FIELD).sum())}
 inputs += [fpath,apath,hp,pp]
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();inputs += [P/'simulation-residual-protocol.md',Path(__file__),P/'validation1020/frozen-discovery-coefficients.npz']
report={'metadata_only_freeze':True,'no_simulated_epoch_residual_scores_inspected':True,'redshift_edges':edges,'arms':summary,'inputs_sha256':{str(p.relative_to(ROOT)):sha(p) for p in inputs},'cohort_sha256':{p.name:sha(p) for p in O.iterdir() if p.is_file()}}
(O/'manifest.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(summary,indent=2))
