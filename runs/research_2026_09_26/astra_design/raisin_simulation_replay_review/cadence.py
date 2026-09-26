from pathlib import Path
import gzip,csv,json,collections
import numpy as np
R=Path.cwd();P=Path(__file__).parent;S=R/'runs/research_2026_09_26/raisin_simulation_execution_review/nominal_fitres'
def table(p,tag):
 with gzip.open(p,'rt') as f:
  out=[]
  for l in f:
   a=l.split()
   if not a:continue
   if a[0]=='VARNAMES:':cols=a[1:]
   if a[0]==tag:out.append(dict(zip(cols,a[1:])))
 return out
fit=table(S/'nir.FITRES.gz','SN:');fd={int(r['CID']):r for r in fit};pl=table(S/'lcplot-feasibility/FITOPT000.LCPLOT.gz','OBS:');pd=collections.defaultdict(list)
for r in pl:pd[int(r['CID'])].append(r)
lib={};active=None
lp=R/'runs/research_2026_09_26/raisin_sign_source/sim/simlibs/DES_RAISIN.simlib'
for line in lp.read_text().splitlines():
 a=line.split()
 if not a:continue
 if a[0]=='LIBID:':active=int(a[1]);lib[active]={'epochs':[]}
 if a[0]=='REDSHIFT:':lib[active].update(z=float(a[1]),peak=float(a[3]))
 if a[0]=='S:' and a[3] in ['J','H']:lib[active]['epochs'].append((a[3],float(a[1])))
rows=[];fc=np.array([float(r['PKMJD']) for r in fit]);fq=np.array([float(r['FITCHI2']) for r in fit]);fn=np.array([int(r['NDOF'])+1 for r in fit]);fl=np.array([int(r['SIM_LIBID']) for r in fit]);fids=np.array([int(r['CID']) for r in fit])
for cid,rr in sorted(pd.items()):
 d=[r for r in rr if r['DATAFLAG']=='1'];m=[r for r in rr if r['DATAFLAG']=='0'];clock=np.median([float(r['MJD'])-float(r['Tobs']) for r in m]);q=sum(float(r['CHI2']) for r in d)
 # All data rows must have distinct same-band nearest cadence epochs. Subset only: no assertion about missing/accepted epochs.
 gaps={}
 for lid,v in lib.items():
  used=set();ds=[]
  for e in d:
   choices=[(abs(float(e['MJD'])-t),i) for i,(b,t) in enumerate(v['epochs']) if b==e['BAND'] and i not in used]
   if not choices:ds.append(float('inf'));break
   gap,i=min(choices);used.add(i);ds.append(gap)
  gaps[lid]=max(ds)
 candidates={str(t):[k for k,v in gaps.items() if v<=t] for t in [.003,.0056,.01]}
 c=candidates['0.0056'];source=fd[cid];n=len(d)
 mask=(abs(fc-clock)<=.0056)&(fn==n)&np.isin(fl,c)
 # Q matching is only descriptive; marginal per-row CHI2 need not sum a correlated fit objective.
 qmask=mask&(abs(fq-q)<=n*.005+1e-4)
 row={'CID':cid,'data_rows':n,'model_clock':clock,'sameCID_LIBID':source['SIM_LIBID'],'cadence_candidates_0p003':','.join(map(str,candidates['0.003'])),'cadence_candidates_0p0056':','.join(map(str,c)),'cadence_candidates_0p01':','.join(map(str,candidates['0.01'])),'sameCID_cadence_max_gap':gaps[int(source['SIM_LIBID'])],'clock_count_cadence_candidate_count':int(mask.sum()),'clock_count_cadence_first20':','.join(map(str,fids[mask][:20])),'extra_sumCHI2_candidate_count_descriptive':int(qmask.sum()),'extra_sumCHI2_CIDs_descriptive':','.join(map(str,fids[qmask])),'sameCID_sumCHI2_difference_descriptive':q-float(source['FITCHI2']),'sameCID_clock_gap':float(clock-float(source['PKMJD']))}
 rows.append(row)
with (P/'cadence-500.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
res={'source_SIMLIBs':len(lib),'epochs_JH':{str(k):len(v['epochs']) for k,v in lib.items()},'plot_objects':len(rows),'sameCID_cadence_matches':{str(t):sum(r['sameCID_cadence_max_gap']<=t for r in rows) for t in [.003,.0056,.01]},'candidate_cardinalities':{t:dict(collections.Counter(0 if not r[k] else len(r[k].split(',')) for r in rows)) for t,k in [('0.003','cadence_candidates_0p003'),('0.0056','cadence_candidates_0p0056'),('0.01','cadence_candidates_0p01')]},'sameCID_peak_gap_abs_quantiles':{str(t):float(np.quantile(abs(np.array([r['sameCID_clock_gap'] for r in rows])),t)) for t in [0,.5,.9,1]},'clock_count_cadence_ambiguous_or_no_candidate':sum(r['clock_count_cadence_candidate_count']!=1 for r in rows),'descriptive_extraQ_matching_cardinality':dict(collections.Counter(r['extra_sumCHI2_candidate_count_descriptive'] for r in rows)),'first20':rows[:20]}
(P/'cadence-result.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps({k:v for k,v in res.items() if k!='first20'},indent=2))
