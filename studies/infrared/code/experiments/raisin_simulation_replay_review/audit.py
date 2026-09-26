from pathlib import Path
import gzip,csv,json,collections,hashlib
import numpy as np
R=Path.cwd();P=Path(__file__).parent;S=R/'runs/research_2026_09_26/raisin_simulation_execution_review/nominal_fitres'
FIT=S/'nir.FITRES.gz';LC=S/'lcplot-feasibility/FITOPT000.LCPLOT.gz'
def parse(p,tag):
 op=gzip.open if p.suffix=='.gz' else open
 with op(p,'rt') as f:
  rows=[];cols=None
  for line in f:
   a=line.split()
   if not a:continue
   if a[0]=='VARNAMES:':cols=a[1:]
   if a[0]==tag:
    assert cols and len(a)-1==len(cols);rows.append(dict(zip(cols,a[1:])))
 return rows
fit=parse(FIT,'SN:');plots=parse(LC,'OBS:');fd={int(x['CID']):x for x in fit};assert len(fd)==len(fit)
pd=collections.defaultdict(list)
for r in plots:pd[int(r['CID'])].append(r)
ledger=[]
for cid,rr in sorted(pd.items()):
 data=[r for r in rr if int(r['DATAFLAG'])==1];model=[r for r in rr if int(r['DATAFLAG'])==0];rej=[r for r in rr if int(r['DATAFLAG'])==-1]
 pkdata=np.array([float(r['MJD'])-float(r['Tobs']) for r in data]);pkmodel=np.array([float(r['MJD'])-float(r['Tobs']) for r in model]);pk=float(np.median(pkmodel))
 row={'CID':cid,'data':len(data),'rejected':len(rej),'model':len(model),'bands':','.join(sorted({r['BAND'] for r in data})),'IFITs':','.join(sorted({r['IFIT'] for r in rr})),'model_peak_median':pk,'model_peak_span':float(np.ptp(pkmodel)),'data_peak_median':float(np.median(pkdata)),'data_peak_span':float(np.ptp(pkdata)),'data_model_peak_max_difference':float(np.max(abs(pkdata-pk))),'same_CID_NDOF_plus1':int(fd[cid]['NDOF'])+1,'same_CID_peak':float(fd[cid]['PKMJD']),'same_CID_peak_difference':pk-float(fd[cid]['PKMJD']),'same_CID_truth_peak_difference':pk-float(fd[cid]['SIM_PKMJD']),'sum_plot_CHI2':sum(float(r['CHI2']) for r in data),'same_CID_FITCHI2':float(fd[cid]['FITCHI2'])}
 for off in range(-3,4):
  if cid+off in fd:
   row[f'offset{off}_peak_gap']=pk-float(fd[cid+off]['PKMJD']);row[f'offset{off}_count_gap']=len(data)-(int(fd[cid+off]['NDOF'])+1)
  else:row[f'offset{off}_peak_gap']=None;row[f'offset{off}_count_gap']=None
 # Conditional candidate audit only: find all table rows with the rounded inferred clock.
 candidates=[k for k,v in fd.items() if abs(float(v['PKMJD'])-pk)<=.0056]
 row['clock_candidate_count']=len(candidates);row['first20_clock_CIDs']='|'.join(map(str,candidates[:20]))
 ledger.append(row)
with (P/'coherence-500.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(ledger[0]));w.writeheader();w.writerows(ledger)
res={'FITRES_rows':len(fit),'LCPLOT_rows':len(plots),'plot_objects':len(pd),'plot_CID_range':[min(pd),max(pd)],'DATAFLAG_counts':dict(collections.Counter(r['DATAFLAG'] for r in plots)),'IFIT_counts':dict(collections.Counter(r['IFIT'] for r in plots)),'sameCID_count_matches':sum(r['data']==r['same_CID_NDOF_plus1'] for r in ledger),'sameCID_peak_matches_0p0056d':sum(abs(r['same_CID_peak_difference'])<=.0056 for r in ledger),'max_internal_data_model_clock_gap':max(r['data_model_peak_max_difference'] for r in ledger),'max_model_clock_span':max(r['model_peak_span'] for r in ledger),'offsets':{str(o):{'available':sum(r[f'offset{o}_peak_gap'] is not None for r in ledger),'peak_matches':sum(r[f'offset{o}_peak_gap'] is not None and abs(r[f'offset{o}_peak_gap'])<=.0056 for r in ledger),'count_matches':sum(r[f'offset{o}_count_gap']==0 for r in ledger)} for o in range(-3,4)},'first20':ledger[:20]}
(P/'coherence-result.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps({k:v for k,v in res.items() if k!='first20'},indent=2));print('first20:',[(r['CID'],r['data'],r['model_peak_median'],r['first20_clock_CIDs']) for r in ledger[:20]])
