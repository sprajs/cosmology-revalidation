"""Parse output-only CSP audit blocks, preserving callback multiplicity."""
from pathlib import Path
import json,hashlib,sys
import numpy as np
O=Path(__file__).resolve().parent
ROW=['MJD','Trest','modelF','model_mag_error','dataF','data_error','z','MWEBV','rest_filter_mean_wavelength','fudge_flux_error','model_flux_error','rest_lambda_fit_min','rest_lambda_fit_max']
OBJ=['totalQ','priorQ','sigmaQ','D','shape','AV','peak_absolute','peak_prior_center_absolute','search_peak','MJDOFF']
ENTRY=['D_entry','shape_entry','AV_entry','peak_entry_absolute','D_step','shape_step','AV_step','peak_step','D_previous','shape_previous','AV_previous','peak_previous_absolute']
def parse(path):
 blocks=[];b=None;entry=None
 def finish():
  if b is None:return
  n=b['NFITDATA'];assert len(b['rows'])==n and len(b['w'])==n
  assert [r['accepted_index'] for r in b['rows']]==list(range(1,n+1))
  a=np.array([[r[k] for k in ROW] for r in b['rows']]);W=np.array(b['w']) if b['USE_FITCOV'] else np.diag(np.array(b['w'])[:,0]);assert W.shape==(n,n)
  b['array']=a;b['W']=W;b['C']=np.linalg.inv(W);blocks.append(b)
 for line in path.read_text().splitlines():
  t=line.split()
  if not t or not t[0].startswith('CSP_'):continue
  tag=t[0]
  if tag=='CSP_ENTRY:':entry={'CID':t[1],'ITER':int(t[2]),'NFITDATA':int(t[3]),'USE_FITCOV':t[4]=='T','LREPEAT_ITER':t[5]=='T',**dict(zip(ENTRY,map(float,t[6:])))};assert len(t)==18
  elif tag=='CSP_START_SHIFT:':
   assert len(t)==6
  elif tag=='CSP_ROW:':
   if b is None or 'objective' in b:
    finish();b={'CID':t[1],'ITER':int(t[2]),'entry':entry.copy(),'rows':[],'w':[]}
   assert (b['CID'],b['ITER'])==(t[1],int(t[2]));assert len(t)==19
   b['rows'].append({'accepted_index':int(t[3]),'source_epoch':int(t[4]),'band':t[5],**dict(zip(ROW,map(float,t[6:])))})
  elif tag=='CSP_OBJECTIVE:':
   assert len(t)==15 and b is not None
   assert (b['CID'],b['ITER'])==(t[1],int(t[2]));b['NFITDATA']=int(t[3]);b['USE_FITCOV']=t[4]=='T';b['objective']=dict(zip(OBJ,map(float,t[5:])))
  elif tag in ['CSP_WROW:','CSP_WDIAG:']:
   assert (b['CID'],b['ITER'])==(t[1],int(t[2]));assert int(t[3])==len(b['w'])+1;b['w'].append(list(map(float,t[4:])))
  else:raise ValueError(tag)
 finish();return blocks

def review(path,out):
 out.mkdir(exist_ok=False);blocks=parse(path);ledger=[];arrays={}
 for i,b in enumerate(blocks):
  x=b['array'];W=b['W'];C=b['C'];ob=b['objective'];y=x[:,ROW.index('dataF')];f=x[:,ROW.index('modelF')];r=y-f;Q=float(r@W@r);amp=float(f@W@y/(f@W@f));delta=float(-2.5*np.log10(amp)) if amp>0 else None
  rows=dict(block=i,CID=b['CID'],ITER=b['ITER'],n=b['NFITDATA'],entry=b['entry'],objective=ob,Q_data_reconstructed=Q,Q_total_reconstructed=Q+ob['priorQ']+ob['sigmaQ'],Q_closure=Q+ob['priorQ']+ob['sigmaQ']-ob['totalQ'],W_asymmetry=float(np.max(np.abs(W-W.T))),C_min_eigenvalue=float(np.linalg.eigvalsh(C).min()),amplitude=amp,frozen_C_delta_D=delta,phase_min=float(x[:,1].min()),phase_max=float(x[:,1].max()),nonpositive_model_count=int(np.sum(f<=0)))
  ledger.append(rows)
  pre=f'b{i:03d}';arrays.update({pre+'_rows':x,pre+'_W':W,pre+'_C':C,pre+'_bands':np.array([r['band'] for r in b['rows']]),pre+'_source_epoch':np.array([r['source_epoch'] for r in b['rows']])})
 pairs=[]
 for cid in dict.fromkeys(b['CID'] for b in blocks):
  ids=[i for i,b in enumerate(blocks) if b['CID']==cid];a,b=blocks[ids[-2]],blocks[ids[-1]]
  identity=all(np.array_equal(a[k],b[k]) for k in ['array']) # flux changes; use data/identity columns explicitly below
  idx=[ROW.index(x) for x in ['MJD','dataF','data_error']]
  identity=np.array_equal(a['array'][:,idx],b['array'][:,idx]) and [r['band'] for r in a['rows']]==[r['band'] for r in b['rows']]
  gap=b['objective']['D']-a['objective']['D'];Cdelta=float(np.linalg.norm(b['C']-a['C'])/np.linalg.norm(a['C'])) if identity else None
  pairs.append({'CID':cid,'final_two_blocks':ids[-2:],'identity':identity,'delta_D':gap,'C_relative_Frobenius_change':Cdelta,'last_two_D_gate':abs(gap)<=.001,'final_frozen_C_D_gate':abs(ledger[ids[-1]]['frozen_C_delta_D'])<=.001})
 result={'source_log':str(path),'log_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'rows_schema':ROW,'blocks':ledger,'last_two_iteration_comparisons':pairs,'max_abs_objective_closure':max(abs(r['Q_closure']) for r in ledger),'all_objectives_close':all(abs(r['Q_closure'])<1e-8 for r in ledger),'all_last_two_and_amplitude_gates':all(r['last_two_D_gate'] and r['final_frozen_C_D_gate'] for r in pairs)}
 np.savez_compressed(out/'native.npz',**arrays);(out/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['blocks','rows_schema']},indent=2))
if __name__=='__main__':review(Path(sys.argv[1]),Path(sys.argv[2]))
