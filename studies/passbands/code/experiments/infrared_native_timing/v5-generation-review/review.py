from pathlib import Path
import json,hashlib,collections,math
import numpy as np
from astropy.io import fits
O=Path(__file__).resolve().parent;P=O.parent;R=Path('/home/szymon/Documents/ChatGPT/supernova')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,x):(O/name).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
schema=json.loads((P/'instrumentation/schema-v2.json').read_text())
def readlog(branch):
 d={k:[] for k in schema}
 for line in (P/'generation-v5'/branch/'native.log').read_text().splitlines():
  v=line.split()
  if v and v[0] in schema:
   assert len(v)==len(schema[v[0]])+1
   x=dict(zip(schema[v[0]],v[1:]))
   for k in x:
    if k not in ['band','field','stage']:x[k]=float(x[k])
   d[v[0]].append(x)
 return d
def load(branch):
 d=P/'generation-v5'/branch/'output/PTE';h=next(d.glob('*HEAD.FITS'));q=next(d.glob('*PHOT.FITS'))
 with fits.open(h) as f:head=f[1].data.copy()
 with fits.open(q) as f:phot=f[1].data.copy()
 return head,phot,[h,q]
inputs=[];tables={}
for b in ['original','ledger','noiseless']:
 h,q,pp=load(b);tables[b]=(h,q);inputs+=pp+[P/'generation-v5'/b/'native.log',P/'generation-v5'/b/'execution.json']
for x,y in zip(tables['original'],tables['ledger']):
 assert x.dtype==y.dtype and x.shape==y.shape
 for key in x.dtype.names:
  assert np.array_equal(x[key],y[key],equal_nan=True) if x[key].dtype.kind in 'fc' else np.array_equal(x[key],y[key])
cadence={}
for line in (P/'inputs-v2/cadence.simlib').read_text().splitlines():
 v=line.split()
 if v and v[0]=='S:':
  key=(float(v[1]),v[3]);assert key not in cadence;cadence[key]={'gain':float(v[4]),'readnoise':float(v[5]),'skysig':float(v[6]),'zp':float(v[10])}
assert len(cadence)==117
maxima=collections.defaultdict(float);serialization=[];objects=[];reference=None;reported=[];all_counts={}
def eq(name,a,b):
 delta=abs(a-b)/max(1.,abs(a),abs(b));maxima[name]=max(maxima[name],delta);assert delta<2e-12,(name,a,b,delta)
for branch,n in [('ledger',8),('noiseless',1)]:
 d=readlog(branch);head,phot=tables[branch]
 assert len(head)==n and len(d['PROSP_EVENT'])==n and not d['PROSP_REJECT']
 assert [x['actual_attempt'] for x in d['PROSP_ATTEMPT']]==list(range(1,n+1))
 assert [x['loop_counter'] for x in d['PROSP_ATTEMPT']]==list(range(1,n+1))
 assert {int(x['cid']) for x in d['PROSP_EVENT']}==set(range(1,n+1))
 all_counts[branch]={k:len(v) for k,v in d.items()}
 for cid in range(1,n+1):
  ev=next(x for x in d['PROSP_EVENT'] if x['cid']==cid);h=next(x for x in head if int(x['SNID'])==cid)
  assert ev['attempt']==cid and ev['libid']==11 and ev['nepoch']==123 and ev['forced_accept']==0
  assert ev['peak_true']==ev['peak_header']==float(h['PEAKMJD'])==57707.80078125
  assert ev['zhel_true']==ev['zcmb_true']==float(h['REDSHIFT_HELIO'])==float(h['REDSHIFT_FINAL'])==0.453000009059906
  assert ev['shape']==1 and ev['AV']==0 and ev['RV']==1.5180000066757202
  assert ev['MWEBV_map']==ev['MWEBV_true']==.006 and ev['MWEBV_error']==0
  assert float(h['MWEBV'])==float(np.float32(.006)) and h['MWEBV_ERR']==0
  rr=[x for x in d['PROSP_ROW'] if x['cid']==cid];noise=[x for x in d['PROSP_NOISE'] if x['cid']==cid]
  assert [x['epoch'] for x in rr]==list(range(1,124))
  assert [x['epoch'] for x in noise]==list(range(1,118))
  wr={int(x['epoch']):x for x in rr if x['obsflag_write']==1};assert len(wr)==117
  assert all(x['obsflag_gen']==x['obsflag_write']==0 and x['obsflag_peak']==1 and x['mjd']==-9 for x in rr[117:])
  rows=phot[int(h['PTROBS_MIN'])-1:int(h['PTROBS_MAX'])];assert len(rows)==int(h['NOBS'])==117
  bykey={(float(x['MJD']),x['BAND'].strip()):x for x in rows};assert set(bykey)==set(cadence) and len(bykey)==117
  means={}
  for x in noise:
   ep=int(x['epoch']);key=(x['mjd'],x['band']);y=bykey[key];w=wr[ep]
   assert x['attempt']==cid and w['attempt']==cid and x['mjd']==w['mjd']
   assert x['field']=='DES16E1dcx' and y['FIELD'].strip()=='DES16E1dcx'
   assert w['obsflag_gen']==1 and w['obsflag_peak']==w['obsflag_template']==0 and w['photflag']==int(y['PHOTFLAG'])==0
   assert float(y['FLUXCAL'])==w['fluxcal_native_R4'] and float(y['FLUXCALERR'])==w['fluxcal_error_native_R4']
   assert x['saturation_excess_pe']<0 and x['random_template_option']==0 and x['redcov_index']==-9
   assert x['smearflag_flux']==(1 if branch=='ledger' else 0)
   assert all(x[k]==0 for k in ['zp_variance_pe2','template_source_pe','template_sky_variance_pe2','host_variance_pe2','true_T_variance_pe2','true_F_variance_pe2','shift_T_pe','shift_F_pe'])
   eq('source_flux_is_source_variance',x['flux_true_pe'],x['src_variance_pe2'])
   q=x['src_variance_pe2']+x['sky_variance_pe2']
   for k in ['true_S_variance_pe2','true_SZ_variance_pe2','true_SUM_variance_pe2','calc_reported_variance_pe2','reported_variance_before_realization_pe2']:eq(k,q,x[k])
   eq('SZ_shift',math.sqrt(q)*x['gauss_SZ'],x['shift_SZ_pe'])
   eq('observed_noise_equation',x['flux_true_pe']+x['shift_SZ_pe'],x['flux_observed_pe'])
   adjustment=max(x['flux_observed_pe'],0)-x['flux_true_pe']
   eq('reported_adjustment',adjustment,x['reported_variance_adjustment_pe2'])
   eq('reported_final',q+adjustment,x['reported_variance_final_pe2'])
   eq('ADU_flux',x['flux_observed_pe']*x['ADU_per_Npe'],x['flux_ADU'])
   eq('ADU_error',math.sqrt(x['reported_variance_final_pe2'])*x['ADU_per_Npe'],x['error_ADU'])
   eq('phase', (x['mjd']-ev['peak_true'])/(1+ev['zhel_true']),x['trest'])
   scale=10**(-.4*(float(y['ZEROPT'])-27.5))
   pred_flux=np.float32(float(np.float32(x['flux_ADU']))*scale)
   pred_error=np.float32(float(np.float32(x['error_ADU']))*scale)
   assert pred_flux==y['FLUXCAL'] and pred_error==y['FLUXCALERR'],(branch,cid,ep,'R4 serialization')
   # These are equality checks of generator expectations, not an inferred empirical variance.
   means[key]=[x[k] for k in ['flux_true_pe','src_variance_pe2','sky_variance_pe2','true_SUM_variance_pe2','reported_variance_before_realization_pe2','Npe_per_fluxcal','NEA']]
   reported.append(x['reported_variance_final_pe2']/q)
   if branch=='noiseless':assert x['gauss_SZ']==x['gauss_T']==x['gauss_F']==x['shift_SZ_pe']==x['reported_variance_adjustment_pe2']==0
  if reference is None:reference=means
  assert means==reference,'mean or true variance changed between draws/noiseless'
  objects.append({'branch':branch,'CID':cid,'written_rows':len(rows),'noise_rows':len(noise),'peak_placeholders':6,'negative_rows':int(np.sum(rows['FLUXCAL']<0)),'zero_rows':int(np.sum(rows['FLUXCAL']==0)),'truth_D':ev['DLMU_true']})
# The repair does not alter any prior source expression beyond two extents.
source_identity=json.loads((P/'repair-v5/source-identity-check.json').read_text());assert all(x['undo_exact_two_bounds_repairs_recovers_previous_source'] for x in source_identity)
for name in ['original','ledger']:
 old=(P/'instrumentation'/('snlc_sim.original.c' if name=='original' else 'snlc_sim.instrumented.c')).read_text()
 new=(P/'repair-v5'/(name+'.c')).read_text().replace(',SNRMAX_FILT[MXCUTWIN_SNRMAX+1][MXFILTINDX]',',SNRMAX_FILT[MXCUTWIN_SNRMAX][MXFILTINDX]').replace('MEM=NEP; if (NEP < MXFILTINDX) { MEM=MXFILTINDX; }\n  MEM++; // one-based TLIST/INDEX_SORT require indices 1..MEM-1','MEM=NEP; if (NEP < MXFILTINDX) { MEM=MXFILTINDX+1; }')
 assert old==new
result={'pass':True,'no_native_calls':True,'independent_parser_no_executor_import':True,'original_vs_ledger_all_native_columns_exact':True,'records':all_counts,'objects':objects,'equation_max_scaled_errors':dict(maxima),'all_1053_written_flux_and_error_values_reproduce_source_R4_conversion_exactly':True,'fixed_true_mean_and_variance_identical_across_all_eight_draws_and_noiseless':True,'reported_to_true_variance_ratio_range_all_rows': [min(reported),max(reported)],'repair_value_preservation':'Undoing only the two extent changes recovers prior source bytes exactly. Actual original-vs-ledger table identity tests output-hook noninterference on identical repaired sources. No runtime equivalence to the crashing unpatched generator is claimed.','scope':'Engineering generator/count/noise/serialization checks only. No noiseless fitter or survey bias gate is certified here.'}
save('result.json',result)
inputs.extend([P/'inputs-v2/cadence.simlib',P/'instrumentation/schema-v2.json',P/'repair-v5/original.c',P/'repair-v5/ledger.c',O/'review.py',O/'result.json'])
save('manifest.json',{'files':{str(p.relative_to(R)):sha(p) for p in inputs}})
print(json.dumps(result,indent=2))
