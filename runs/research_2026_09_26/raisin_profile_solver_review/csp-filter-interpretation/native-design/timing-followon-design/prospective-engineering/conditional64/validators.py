from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import json,hashlib,re,importlib.util
import numpy as np
from astropy.io import fits
P=Path(__file__).resolve().parent;Q=P;H=P;R=Path('/home/szymon/Documents/ChatGPT/supernova')
source=R/'runs/research_2026_09_26/astra_design/raisin_timing_assets/instrumentation_2021/parse_native.py'
s=importlib.util.spec_from_file_location('original_native_parser',source);parser=importlib.util.module_from_spec(s);s.loader.exec_module(parser)
def check(ok,msg):
 if not ok:raise RuntimeError(msg)
def save(p,d):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def read(branch):
 with fits.open(P/'datasets/ledger/PTE/PTE_HEAD.FITS') as f:h=f[1].data.copy();hc=f[1].columns
 with fits.open(P/'datasets/ledger/PTE/PTE_PHOT.FITS') as f:p=f[1].data.copy();pc=f[1].columns
 return h,p,hc,pc
def base_block_review(w,expected_n,expected_cids,support=True):
 blocks=parser.parse(w/'native.log');by={};detail=[]
 table=[];names=None
 for line in (w/'fit.FITRES.TEXT').read_text().splitlines():
  v=line.split()
  if v and v[0]=='VARNAMES:':names=v[1:]
  elif v and v[0]=='SN:':check(names is not None and len(v)-1==len(names),'FITRES schema');table.append(dict(zip(names,v[1:])))
 check({x['CID'] for x in table}==set(expected_cids) and len(table)==len(expected_cids),'FITRES complete membership');check(all(int(x['ERRFLAG_FIT'])==0 and int(x['CUTFLAG_SNANA'])==3 for x in table),'native error/cut flag')
 nml=(w/'fit.nml').read_text();fixed='INISTP_PEAKMJD = 0' in nml;arm_peaks={}
 if fixed:
  datapath=__import__('re').search(r"PRIVATE_DATA_PATH = '([^']+)'",nml).group(1);headers=list((w/datapath/'PTE').glob('*HEAD.FITS*'));check(len(headers)==1,'arm HEAD')
  with fits.open(headers[0]) as hh:arm_peaks={str(x['SNID']).strip():float(x['PEAKMJD']) for x in hh[1].data}
 for b in blocks:
  by.setdefault(b['CID'],[]).append(b);check(b['NFITDATA']==expected_n,'accepted cadence size');ob=b['objective'];
  if fixed:check(ob['peak_absolute']==arm_peaks[b['CID']],'fixed peak differs from assigned native HEAD')
  x=b['array'];W=b['W'];C=b['C'];r=x[:,4]-x[:,2];q=float(r@W@r)+ob['priorQ']+ob['sigmaQ'];check(abs(q-ob['totalQ'])<1e-7,'objective arithmetic');check(np.linalg.eigvalsh(C).min()>0,'C not SPD');check(np.isfinite(x).all() and np.isfinite(C).all(),'nonfinite state');check(x[:,1].min()>=-20 and x[:,1].max()<=70,'finalcallback phase');check(ob['shape']==1 and ob['AV']==0,'fixed nuisance changed')
 check(set(by)==set(expected_cids),'complete membership')
 branch='noiseless' if w.name.startswith('noiseless') else 'ledger';head,phot,_,_=read(branch)
 for hdr in head:
  c=str(hdr['SNID']).strip();native_rows=phot[int(hdr['PTROBS_MIN'])-1:int(hdr['PTROBS_MAX'])];wanted=[x for x in native_rows if expected_n==117 or str(x['BAND']).strip() in 'JH'];observed=by[c][-1]['rows'];check(len(wanted)==len(observed),'input row count')
  remaining=[(str(x['BAND']).strip(),float(x['MJD']),float(x['FLUXCAL']),float(x['FLUXCALERR'])) for x in wanted]
  for x in observed:
   key=(x['band'],x['MJD'],x['dataF'],x['data_error']);check(key in remaining,'fitter native FITS measurement mismatch '+repr(key));remaining.remove(key)
  check(not remaining,'missing fitter input exposure')
 sup={}
 if support:
  for line in (w/'native.log').read_text().splitlines():
   v=line.split()
   if not v:continue
   check(v[0] not in ['PROSP_SUPPORT_FIRST_BAD','PROSP_SUPPORT_OVERFLOW'],'all-call support violation')
   if v[0]=='PROSP_SUPPORT':check(len(v)==11,'supportschema');check(v[1] not in sup,'duplicate support summary');sup[v[1]]=v[2:];check(int(v[2])>0 and int(v[3])==int(v[4])==0,'support counters')
  check(set(sup)==set(expected_cids),'missing allcall summary')
 for cid,bb in by.items():
  b=bb[-1];prev=next(x for x in reversed(bb[:-1]) if x['ITER']<b['ITER']);ob=b['objective'];last=prev['objective'];check(abs(ob['D']-last['D'])<=.001,'last iterations D');check(abs(ob['peak_absolute']-last['peak_absolute'])<=.01,'last iterations peak');check(abs(ob['peak_absolute']-57707.80078125)<=4,'final peak guard')
  # same native data/order in every callback; FITS exposure MJD remains R8.
  identity=lambda z:[(r['source_epoch'],r['band'],r['MJD'],r['dataF'],r['data_error']) for r in z['rows']]
  check(all(identity(x)==identity(b) for x in bb),'state accepted rows changed');L=np.linalg.cholesky(prev['C']);d=np.linalg.solve(L,b['C']-prev['C']);d=np.linalg.solve(L,d.T).T;metric=float(np.linalg.norm(d,2));check(metric<=.001,'last iterations whitened C')
  x=b['array'];f=x[:,2];y=x[:,4];a=float(f@b['W']@y/(f@b['W']@f));check(a>0,'nonpositive amplitude');da=-2.5*np.log10(a);check(abs(da)<=.001,'fixedC amplitude stationarity');detail.append({'CID':cid,'D':ob['D'],'peak':ob['peak_absolute'],'dataQ':ob['totalQ']-ob['priorQ']-ob['sigmaQ'],'priorQ':ob['priorQ'],'sigmaQ':ob['sigmaQ'],'fixedC_delta_D':da,'last_C_whitened_norm':metric,'max_standardized_model_residual':float(np.max(np.abs((x[:,4]-x[:,2])/x[:,5]))),'mean_calls':int(sup[cid][0]) if support else None})
 return by,detail

def domain_review(w,by,cids):
 mode=json.loads((w/'execution.json').read_text())['hard_domain_mode']
 lines=[x.split() for x in (w/'native.log').read_text().splitlines() if x.startswith(('PROSP_DOMAIN','PROSP_PEAK_BOUNDS'))]
 if mode!=1:
  check(not lines,'disabled mode emitted domain records');return {'enabled':False}
 schema=json.loads((Q/'build-protocol.json').read_text())['schema'];bounds={};summaries={}
 for v in lines:
  check(v[0] in ['PROSP_PEAK_BOUNDS','PROSP_DOMAIN_SUMMARY'],'domain abort/reject or unknown record')
  check(len(v)-1==len(schema[v[0]]),'domain schema mismatch');d=dict(zip(schema[v[0]],v[1:]));c=d['CID']
  check(c in cids,'unexpected domain CID')
  if v[0]=='PROSP_PEAK_BOUNDS':bounds.setdefault(c,[]).append(d)
  else:check(c not in summaries,'duplicate domain summary');summaries[c]=d
 check(set(bounds)==set(summaries)==set(cids),'missing domain records')
 head,phot,_,_=read('noiseless' if w.name.startswith('noiseless') else 'ledger')
 nit=int(re.search(r'NFIT_ITERATION\s*=\s*(\d+)',(w/'fit.nml').read_text()).group(1));out=[]
 for c in cids:
  hh=next(x for x in head if str(x['SNID']).strip()==c);tt=phot[int(hh['PTROBS_MIN'])-1:int(hh['PTROBS_MAX'])]['MJD'];z=float(hh['REDSHIFT_HELIO'])
  lo=float(np.max(tt-(1+z)*70));hi=float(np.min(tt+(1+z)*20));check(len(tt)==117,'domain metadata count')
  bb=bounds[c];ss=summaries[c]
  check(len(bb)==nit and {int(x['iteration']) for x in bb}==set(range(1,nit+1)),'bounds iteration coverage')
  check({b['ITER'] for b in by[c]}==set(range(1,nit+1)),'callback iteration coverage')
  for x in bb:
   check(int(x['n_all_metadata'])==117 and float(x['zHEL'])==z,'bound metadata changed')
   check([float(x['phase_lower']),float(x['phase_upper']),float(x['KCOR_lower']),float(x['KCOR_upper'])]==[-20,70,-20,85],'table phase bounds changed')
   check(float(x['raw_absolute_lower'])==float(x['safe_absolute_lower'])==lo and float(x['raw_absolute_upper'])==float(x['safe_absolute_upper'])==hi,'bound arithmetic/inward adjustment mismatch')
   off=float(x['MJDOFF']);check(off==0,'frozen engineering MJDOFF changed');check(float(x['parameter_lower'])==lo-off and float(x['parameter_upper'])==hi-off,'MJDOFF bound coordinates')
   pk=float(x['initial_absolute_peak']);check(float(x['initial_distance_lower_days'])==pk-lo and float(x['initial_distance_upper_days'])==hi-pk,'initial endpoint distances')
   check(lo<=pk<=hi,'initial peak outside domain')
  check(int(ss['configured'])==1 and int(ss['bounds_assertion_count'])==nit,'domain setup incomplete')
  check(int(ss['physical_mean_calls'])>0 and int(ss['FCN_physical_calls'])>0,'missing guarded physical calls')
  check(float(ss['min_phase'])>=-20 and float(ss['max_phase'])<=70 and float(ss['minimum_phase_edge_distance_rest_days'])>=0,'guarded phase failure')
  check(float(ss['min_FCN_absolute_peak'])>=lo and float(ss['max_FCN_absolute_peak'])<=hi and float(ss['minimum_peak_edge_distance_observer_days'])>=0,'guarded peak failure')
  support=next(x.split() for x in (w/'native.log').read_text().splitlines() if x.startswith('PROSP_SUPPORT '+c+' '));check(int(ss['physical_mean_calls'])==int(support[2]),'guard versus audit mean-call coverage')
  final=by[c][-1]['objective']['peak_absolute'];check(lo<=final<=hi,'final peak outside domain')
  out.append(dict(CID=c,absolute_bounds=[lo,hi],final_peak=final,final_lower_margin_days=final-lo,final_upper_margin_days=hi-final,final_exactly_on_boundary=final in [lo,hi],search_minimum_edge_distance_days=float(ss['minimum_peak_edge_distance_observer_days']),mean_calls=int(ss['physical_mean_calls']),FCN_physical_calls=int(ss['FCN_physical_calls'])))
 return {'enabled':True,'details':out,'interpretation':'Search endpoint contact is separate from final active-bound status.'}
def records(path,tag):return [x.split() for x in path.read_text().splitlines() if x.startswith(tag+' ')]

def shift_review(w,by,cids,domain):
 execution=json.loads((w/'execution.json').read_text());shift=execution.get('minuit_peak_shift')
 a=records(w/'native.log','PROSP_MNPARM_PEAK');b=records(w/'native.log','PROSP_MNPOUT_PEAK')
 aborts=records(w/'native.log','PROSP_MNPARM_ABORT');check(not aborts,'MINUIT hook abort')
 if shift in [None,0]:check(not a and not b,'default hook emitted start records');return {'active':False}
 check(shift in [-2,2],'undeclared peak shift')
 nit=int(re.search(r'NFIT_ITERATION\s*=\s*(\d+)',(w/'fit.nml').read_text()).group(1));check(nit==12,'start branch iteration count')
 schema=json.loads((H/'schema.json').read_text());parsed={}
 for tag,rows in [('PROSP_MNPARM_PEAK',a),('PROSP_MNPOUT_PEAK',b)]:
  out={}
  for row in rows:
   check(len(row)-1==len(schema[tag]),'start record schema');d=dict(zip(schema[tag],row[1:]));key=(d['CID'],int(d['iteration']));check(key not in out,'duplicate start record');out[key]=d
  check(set(out)=={(c,it) for c in cids for it in range(1,nit+1)},'start record coverage');parsed[tag]=out
 original=P/'fits/joint12'
 # Full source entry identity at iteration1, before the start-only hook.
 own_entries={x[1]:x for x in records(w/'native.log','CSP_ENTRY:') if int(x[2])==1}
 base_entries={x[1]:x for x in records(original/'native.log','CSP_ENTRY:') if int(x[2])==1}
 check(own_entries==base_entries,'first common initialization/prior state changed')
 own_bounds=[x for x in records(w/'native.log','PROSP_PEAK_BOUNDS') if int(x[2])==1]
 base_bounds=[x for x in records(original/'native.log','PROSP_PEAK_BOUNDS') if int(x[2])==1]
 check(own_bounds==base_bounds,'first pre-grid metadata/initialization changed')
 details=[]
 for c in cids:
  dd=next(x for x in domain['details'] if x['CID']==c);lo,hi=dd['absolute_bounds']
  for it in range(1,nit+1):
   aa=parsed['PROSP_MNPARM_PEAK'][c,it];bb=parsed['PROSP_MNPOUT_PEAK'][c,it]
   cb=next(x for x in by[c] if x['ITER']==it);entry=cb['entry']['peak_entry_absolute']
   # Objective schema position12 is the actual common INIVAL prior center.
   obj=next(x for x in records(w/'native.log','CSP_OBJECTIVE:') if x[1]==c and int(x[2])==it)
   prior=float(obj[12]);source=float(aa['native_source_parameter']);stored=float(bb['stored_MINUIT_parameter']);proposed=float(aa['proposed_MNPARM_parameter']);expected=shift if it==1 else 0
   check(source==entry==prior==float(aa['unchanged_prior_center_parameter'])==float(bb['unchanged_prior_center_parameter'])==float(bb['native_source_parameter']),'shared prior source changed')
   check(float(aa['requested_shift_days'])==shift and float(aa['applied_shift_days'])==expected,'applied shift metadata')
   check(proposed==source+expected and proposed==float(bb['proposed_MNPARM_parameter']),'proposed shift arithmetic')
   check(abs(stored-source-expected)<.004 and abs(stored-proposed)<.004,'actual stored start separation')
   if it>1:check(stored==source,'later optimizer entry changed')
   check(bb['parameter_name']=='PKMJD' and int(bb['internal_index'])>0,'MINUIT stored parameter identity')
   check(float(aa['lower_bound'])==float(bb['input_lower_bound'])==float(bb['stored_lower_bound'])==lo,'lower bound changed')
   check(float(aa['upper_bound'])==float(bb['input_upper_bound'])==float(bb['stored_upper_bound'])==hi,'upper bound changed')
   check(float(aa['distance_lower_days'])==proposed-lo and float(aa['distance_upper_days'])==hi-proposed,'entry endpoint margins')
   if it==1:details.append({'CID':c,'source_and_prior':source,'actual_MINUIT_start':stored,'actual_shift_days':stored-source,'requested_shift_days':shift})
 return {'active':True,'first_source_CSP_ENTRY_and_pregrid_bounds_exact':True,'details':details,'scope':'First objective function/state preserved; W may depend on coordinates. Later C/prior histories may differ.'}
