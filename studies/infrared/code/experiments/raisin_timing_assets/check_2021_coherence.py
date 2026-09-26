from pathlib import Path
from decimal import Decimal as D
import gzip,csv,json,collections,hashlib
R=Path.cwd();P=Path(__file__).parent;O=P/'author_20211111'
def read(p,tag):
 out=[]
 with gzip.open(p,'rt') as f:
  for l in f:
   a=l.split()
   if not a:continue
   if a[0]=='VARNAMES:':keys=a[1:]
   if a[0]==tag:assert len(a)==len(keys)+1;out.append(dict(zip(keys,a[1:])))
 return out
f=read(O/'nir.FITRES.gz','SN:');g=read(O/'optnir.FITRES.gz','SN:');fd={r['CID']:r for r in f};gd={r['CID']:r for r in g};assert len(fd)==len(f) and len(gd)==len(g)
pl=read(R/'runs/research_2026_09_26/raisin_simulation_execution_review/nominal_fitres/lcplot-feasibility/FITOPT000.LCPLOT.gz','OBS:');pd=collections.defaultdict(list)
for r in pl:pd[r['CID']].append(r)
meta={r['plot_CID']:r for r in csv.DictReader((R/'runs/research_2026_09_26/astra_design/raisin_simulation_replay_review/metadata-candidates-500.csv').open())}
rows=[]
for cid,rr in sorted(pd.items(),key=lambda x:int(x[0])):
 clock={D(r['MJD'])-D(r['Tobs']) for r in rr if r['DATAFLAG']=='0'};assert len(clock)==1;clock=clock.pop();data=[r for r in rr if r['DATAFLAG']=='1'];fit=fd.get(cid)
 rows.append({'CID':cid,'NIR_present':fit is not None,'accepted_rows':len(data),'NDOF_plus1':int(fit['NDOF'])+1,'clock':str(clock),'PKMJD':fit['PKMJD'],'clock_gap':str(clock-D(fit['PKMJD'])),'clock_matches':abs(clock-D(fit['PKMJD']))<=D('.0056'),'identified_cadence':meta[cid]['identified_SIMLIB'],'FITRES_SIM_LIBID':fit['SIM_LIBID'],'cadence_matches':meta[cid]['identified_SIMLIB']==fit['SIM_LIBID']})
with (P/'historical-coherence-500.csv').open('w',newline='') as h:
 w=csv.DictWriter(h,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
common=sorted(fd.keys()&gd.keys(),key=int);truth=[k for k in f[0] if k.startswith('SIM_')]+['PKMJDINI','zHEL','zCMB','zHD','MWEBV'];mismatch={k:sum(fd[c][k]!=gd[c][k] for c in common) for k in truth if k in g[0]}
res={'source_commit':'aaa709ead7a7339d56a2a4604d78991329d9368f','plot_git_blob_unchanged':'d81d64d3f1f395ad9dbbab71837d89f4811320af','NIR_rows':len(f),'joint_rows':len(g),'common_FITRES_CIDs':len(common),'sameCID_count_matches':sum(r['accepted_rows']==r['NDOF_plus1'] for r in rows),'sameCID_clock_matches':sum(r['clock_matches'] for r in rows),'sameCID_cadence_matches':sum(r['cadence_matches'] for r in rows),'max_abs_clock_gap_day':str(max(abs(D(r['clock_gap'])) for r in rows)),'NIR_joint_metadata_mismatch_counts':mismatch,'scope':'Metadata coherence, not baseline native reproduction or recovery of original HEAD/PHOT. No flux/Q/D matching used.'}
(P/'historical-coherence-result.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
