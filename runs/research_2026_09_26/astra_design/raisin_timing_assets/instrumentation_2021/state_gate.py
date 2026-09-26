"""State-only equivalence and initialization gates. No paired timing scores."""
from pathlib import Path
import json,hashlib,collections,sys,numpy as np
from parse_native import parse,ROW
P=Path(__file__).resolve().parent;B=P.parent/'baseline_2021';R=Path('/home/szymon/Documents/ChatGPT/supernova')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read_table(p):
 names=None;out={}
 for line in p.read_text().splitlines():
  if line.startswith('VARNAMES:'):names=line.split()[1:]
  if line.startswith('SN:'):
   r=dict(zip(names,line.split()[1:]));out[r['CID']]=r
 return out

def main():
 assert not (P/'state-result.json').exists()
 baseline=B/'fits/baseline';source=read_table(baseline/'fit.FITRES.TEXT');reports={};arrays={};comparisons=[]
 for cond in ['absent','zero','plus','minus']:
  w=P/'fits'/cond;assert json.loads((w/'execution.json').read_text())['returncode']==0
  fit=read_table(w/'fit.FITRES.TEXT');assert set(fit)==set(source)
  for cid,r in fit.items():assert r['ERRFLAG_FIT']=='0'
  science=lambda p:'\n'.join(x for x in p.read_text().splitlines() if x.startswith(('VARNAMES:','SN:')))
  exact=science(w/'fit.FITRES.TEXT')==science(baseline/'fit.FITRES.TEXT') and (w/'fit.LCPLOT.TEXT').read_bytes()==(baseline/'fit.LCPLOT.TEXT').read_bytes()
  if cond in ['absent','zero']:assert exact,(cond,'output equivalence')
  blocks=parse(w/'fit.log');ledger=[];bycid=collections.defaultdict(list)
  for i,b in enumerate(blocks):
   a=b['array'];W=b['W'];C=b['C'];ob=b['objective'];r=a[:,4]-a[:,2];Q=float(r@W@r);closure=Q+ob['priorQ']+ob['sigmaQ']-ob['totalQ'];f=a[:,2];amp=float(f@W@a[:,4]/(f@W@f));delta=float(-2.5*np.log10(amp)) if amp>0 else None
   assert np.all(np.isfinite(a)) and np.all(np.isfinite(W));assert ob['MJDOFF']==0
   asym=float(np.max(np.abs(W-W.T)));eig=float(np.linalg.eigvalsh(C).min());assert eig>0 and asym<1e-10 and abs(closure)<1e-8
   # Native full MJD stays in R8 through objective. Data/error are exactly R4 ingested values.
   dat=baseline/'data/DES_RAISIN_SIM'/f"{b['CID']}.DAT";raw=[x.split() for x in dat.read_text().splitlines() if x.startswith('OBS:')]
   expected=collections.Counter((x[2],float(x[1]),float(np.float32(float(x[4]))),float(np.float32(float(x[5])))) for x in raw)
   got=collections.Counter((x['band'],x['MJD'],x['dataF'],x['data_error']) for x in b['rows']);assert got==expected,(cond,b['CID'],b['ITER'])
   assert ob['shape']==1 and ob['AV']==0 and ob['peak_absolute']==float(source[b['CID']]['PKMJD']) or (ob['shape']==1 and ob['AV']==0 and format(ob['peak_absolute'],'.4f')==source[b['CID']]['PKMJD'])
   q={'CID':b['CID'],'ITER':b['ITER'],'block':i,'entry':b['entry'],'objective':ob,'Q_data_reconstructed':Q,'objective_closure':closure,'C_min_eigenvalue':eig,'W_asymmetry':asym,'frozen_C_delta_D':delta,'nonpositive_model_rows':int(np.sum(f<=0))};ledger.append(q);bycid[b['CID']].append(q)
   pre=f'{cond}_{i:03d}';arrays.update({pre+'_rows':a,pre+'_W':W,pre+'_C':C,pre+'_bands':np.array([x['band'] for x in b['rows']])})
  reports[cond]={'exact_default_science_outputs':exact,'blocks':ledger,'FITRES':fit}
  for cid in source:
   unique={b['ITER']:b for b in bycid[cid]};assert set(unique)=={1,2,3};last=unique[3];prev=unique[2]
   comparisons.append({'condition':cond,'CID':cid,'last_iteration_D_change':last['objective']['D']-prev['objective']['D'],'frozen_C_delta_D':last['frozen_C_delta_D'],'iteration_gate':abs(last['objective']['D']-prev['objective']['D'])<=.001,'frozen_C_gate':last['frozen_C_delta_D'] is not None and abs(last['frozen_C_delta_D'])<=.001})
 # Check actual hook and fixed coordinates against absent first entries, independent of final results.
 starts=[]
 for cid in source:
  first=lambda c:next(b['entry'] for b in reports[c]['blocks'] if b['CID']==cid and b['ITER']==1)
  ref=first('absent');Ds=[]
  for c,s in [('zero',0.),('plus',.2),('minus',-.2)]:
   en=first(c);assert abs((en['D_entry']-ref['D_entry'])-s)<1e-12
   for k in ['shape_entry','AV_entry','peak_entry_absolute','D_step','shape_step','AV_step','peak_step']:assert en[k]==ref[k]
   Ds.append(float(reports[c]['FITRES'][cid]['DLMAG']))
  Ds.append(float(reports['absent']['FITRES'][cid]['DLMAG']));starts.append({'CID':cid,'DLMAG_range':max(Ds)-min(Ds),'start_gate':max(Ds)-min(Ds)<=.001})
 result={'status':'computationally_complete','protocol_sha256':sha(P/'state-protocol.json'),'default_equivalence_pass':True,'full_state_objective_and_mask_pass':True,'actual_start_shift_gate':True,'all_numerical_gates_pass':all(x['iteration_gate'] and x['frozen_C_gate'] for x in comparisons) and all(x['start_gate'] for x in starts),'iteration_comparisons':comparisons,'start_spreads':starts,'conditions':reports,'scope':'Native algorithm state audit, no timing intervention or population inference.'}
 np.savez_compressed(P/'state-arrays.npz',**arrays);(P/'state-result.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({k:v for k,v in result.items() if k not in ['conditions','iteration_comparisons']},indent=2))
if __name__=='__main__':main()
