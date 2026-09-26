"""Frozen source-precision sensitivity, fixed eight original objects; no timing intervention."""
from pathlib import Path
import json,hashlib,shutil,subprocess,os,time,sys,collections,importlib.util,argparse
import numpy as np
P=Path(__file__).resolve().parent;R=Path('/home/szymon/Documents/ChatGPT/supernova');B=P.parent/'baseline_2021';I=P.parent/'instrumentation_2021'
spec=importlib.util.spec_from_file_location('native_parser',I/'parse_native.py');native=importlib.util.module_from_spec(spec);spec.loader.exec_module(native)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def load(p):return json.loads(p.read_text())
def table(p):
 out={};names=None
 for line in p.read_text().splitlines():
  if line.startswith('VARNAMES:'):names=line.split()[1:]
  if line.startswith('SN:'):
   d=dict(zip(names,line.split()[1:]));out[d['CID']]=d
 return out

def source_gate():
 for f,h in load(P/'protocol.json')['inputs_sha256'].items():assert sha(R/f)==h,f
 assert load(I/'state-result.json')['all_numerical_gates_pass']
 return load(P/'precision-ledger.json')['coordinates']

def make_work(stage,cases):
 w=P/'fits'/stage;w.mkdir(parents=True,exist_ok=False);src=B/'fits/baseline';v=w/'data/DES_RAISIN_SIM';v.mkdir(parents=True)
 shutil.copy2(src/'fit.nml',w/'fit.nml');(w/'kcor.fits').symlink_to((src/'kcor.fits').resolve());shutil.copytree(src/'snoopy.B18',w/'snoopy.B18',symlinks=True)
 shutil.copy2(src/'data/DES_RAISIN_SIM/DES_RAISIN_SIM.README',v/'DES_RAISIN_SIM.README')
 for case in cases:
  original=src/'data/DES_RAISIN_SIM'/f"{case['original_CID']}.DAT";lines=original.read_text().splitlines();epoch=0
  for i,line in enumerate(lines):
   if line.startswith('SNID:'):lines[i]='SNID: '+case['new_CID']
   if line.startswith('OBS:'):
    epoch+=1;t=line.split()
    for ch in case['changes']:
     if ch['epoch']==epoch:t[{'MJD':1,'FLUXCAL':4,'FLUXCALERR':5}[ch['field']]]=format(ch['value'],'.17g')
    lines[i]=' '.join(t)
   for ch in case['changes']:
    if ch['epoch'] is None and line.startswith(ch['field']+':'):
     t=line.split();t[1]=format(ch['value'],'.17g');lines[i]=' '.join(t)
  f=v/f"{case['new_CID']}.DAT";f.write_text('\n'.join(lines)+'\n');case['file']=str(f.relative_to(R));case['sha256']=sha(f)
 (v/'DES_RAISIN_SIM.LIST').write_text(''.join(f"{x['new_CID']}.DAT\n" for x in cases));save(w/'case-ledger.json',cases)
 files=[P/'runner.py',P/'protocol.json',P/'precision-ledger.json',I/'parse_native.py',I/'state-result.json',I/'build/bin/snlc_fit.exe',w/'fit.nml',w/'case-ledger.json']+list(v.iterdir())
 save(w/'freeze.json',{'stage':stage,'files':{str(f.relative_to(R)):sha(f) for f in files},'protocol_sha256':sha(P/'protocol.json')});print(stage,len(cases),sha(w/'freeze.json'))

def prepare_endpoints():
 coords=source_gate();cases=[]
 def add(cid,kind,changes,**extra):cases.append({'new_CID':str(100001+len(cases)),'original_CID':cid,'kind':kind,'changes':changes,**extra})
 for cid in map(str,range(1,9)):add(cid,'clone_control',[])
 for j,c in enumerate(coords):
  if c['input_low']==c['input_high']:continue
  for side in ['low','high']:add(c['CID'],'coordinate_endpoint',[{'field':c['field'],'epoch':c['epoch'],'value':c['input_'+side]}],coordinate_index=j,side=side)
 assert len(cases)==248;make_work('endpoints',cases)


def run(stage):
 source_gate();w=P/'fits'/stage;f=load(w/'freeze.json');release=load(P/'execution-release.json');assert release['protocol_sha256']==sha(P/'protocol.json') and release['runner_sha256']==sha(P/'runner.py')
 if stage=='endpoints':assert release['endpoint_freeze_sha256']==sha(w/'freeze.json')
 else:
  assert load(P/'endpoint-result.json')['all_gates_pass'];assert f['protocol_sha256']==release['protocol_sha256']
 for name,h in f['files'].items():assert sha(R/name)==h,name
 assert not (w/'fit.log').exists()
 before=load(I/'state-protocol.json')['resources']['earlier_native_seconds']+sum(load(x)['wall_seconds'] for x in (I/'fits').glob('*/execution.json'))
 stage_spent=sum(load(x)['wall_seconds'] for x in (P/'fits').glob('*/execution.json'));cfg=load(P/'protocol.json')['resources'];allowed=min(cfg['per_process_seconds'],cfg['stage_seconds']-stage_spent,cfg['original_total_pilot_seconds']-before-stage_spent);assert allowed>0
 env=os.environ.copy();env.pop('CSP_POSTINIT_DSHIFT',None);env.update(SNANA_DIR=str(I/'build'),SNDATA_ROOT=str(R/'phase2/official/inputs/SNDATA_ROOT'),LD_LIBRARY_PATH=str(R/'phase2/official/build/sysroot/usr/lib'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1');start=time.monotonic()
 with (w/'fit.log').open('x') as log:
  try:code=subprocess.run([str(I/'build/bin/snlc_fit.exe'),'fit.nml'],cwd=w,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=allowed).returncode
  except subprocess.TimeoutExpired:code='timeout'
 save(w/'execution.json',{'returncode':code,'wall_seconds':time.monotonic()-start,'timeout_seconds':allowed,'native_seconds_before':before+stage_spent,'release_sha256':sha(P/'execution-release.json'),'freeze_sha256':sha(w/'freeze.json'),'runner_sha256':sha(P/'runner.py')});assert code==0


def analyze(stage):
 w=P/'fits'/stage;cases=load(w/'case-ledger.json');bynew={c['new_CID']:c for c in cases};fit=table(w/'fit.FITRES.TEXT');assert set(fit)==set(bynew)
 baseline=load(I/'state-result.json')['conditions']['absent'];refs={}
 for b in baseline['blocks']:refs[(b['CID'],b['ITER'])]=b
 blocks=native.parse(w/'fit.log');groups=collections.defaultdict(dict);arrays={};state=[]
 for j,b in enumerate(blocks):
  case=bynew[b['CID']];a=b['array'];W=b['W'];ob=b['objective'];y=a[:,4];model=a[:,2];res=y-model;qd=float(res@W@res);cl=qd+ob['priorQ']+ob['sigmaQ']-ob['totalQ'];amp=float(model@W@y/(model@W@model));ad=float(-2.5*np.log10(amp)) if amp>0 else None
  actual=collections.Counter((r['band'],r['MJD'],r['dataF'],r['data_error']) for r in b['rows']);raw=[line.split() for line in (R/case['file']).read_text().splitlines() if line.startswith('OBS:')];expected=collections.Counter((x[2],float(x[1]),float(np.float32(float(x[4]))),float(np.float32(float(x[5])))) for x in raw)
  gate=bool(actual==expected and np.all(np.isfinite(a)) and np.all(np.isfinite(W)) and np.linalg.eigvalsh(b['C']).min()>0 and abs(cl)<=1e-8 and ob['MJDOFF']==0)
  row={'new_CID':b['CID'],'original_CID':case['original_CID'],'ITER':b['ITER'],'native_D':ob['D'],'native_data_Q':qd,'objective_closure':cl,'frozen_C_delta_D':ad,'full_state_gate':gate,'objective':ob,'entry':b['entry']};state.append(row);groups[b['CID']][b['ITER']]=row
  arrays[f'b{j:04d}_rows']=a;arrays[f'b{j:04d}_W']=W;arrays[f'b{j:04d}_band']=np.array([r['band'] for r in b['rows']])
 results=[];published=table(B/'fits/baseline/fit.FITRES.TEXT')
 # Original archivedQ (cap anchor) was explicitly preserved in baseline gate perCID.
 import gzip
 arch={};names=None
 with gzip.open(P.parent/'author_20211111/nir.FITRES.gz','rt') as f:
  for line in f:
   if line.startswith('VARNAMES:'):names=line.split()[1:]
   if line.startswith('SN:'):
    r=dict(zip(names,line.split()[1:]));
    if r['CID'] in published:arch[r['CID']]=r
 for case in cases:
  c=case['new_CID'];cid=case['original_CID'];g=groups[c];assert set(g)=={1,2,3};end=g[3];ref=refs[(cid,3)];deltaD=end['native_D']-ref['objective']['D'];deltaQ=end['native_data_Q']-ref['Q_data_reconstructed'];capQ=max(.01,.001*float(arch[cid]['FITCHI2']));iteration=abs(end['native_D']-g[2]['native_D']);fg=end['frozen_C_delta_D'];numeric=all(x['full_state_gate'] for x in g.values()) and fit[c]['ERRFLAG_FIT']=='0' and iteration<=.001 and fg is not None and abs(fg)<=.001
  clone=True
  if case['kind']=='clone_control':
   clone=all(fit[c][k]==v for k,v in published[cid].items() if k!='CID') and deltaD==0 and deltaQ==0
   # All native model/data/W arrays exact by coordinate against absent branch.
   refsblocks=[b for b in native.parse(I/'fits/absent/fit.log') if b['CID']==cid]
   current=[b for b in blocks if b['CID']==c]
   clone=clone and len(refsblocks)==len(current) and all(np.array_equal(x['array'],y['array']) and np.array_equal(x['W'],y['W']) for x,y in zip(refsblocks,current))
  results.append({**case,'native_D':end['native_D'],'native_data_Q':end['native_data_Q'],'delta_D':deltaD,'delta_data_Q':deltaQ,'data_Q_cap':capQ,'distance_cap_pass':abs(deltaD)<=.001,'data_Q_cap_pass':abs(deltaQ)<=capQ,'iteration_D_change':iteration,'frozen_C_delta_D':fg,'native_gate':bool(numeric),'clone_gate':bool(clone)})
 result={'stage':stage,'all_gates_pass':all(r['distance_cap_pass'] and r['data_Q_cap_pass'] and r['native_gate'] and r['clone_gate'] for r in results),'cases':results,'native_state':state,'protocol_sha256':sha(P/'protocol.json'),'runner_sha256':sha(P/'runner.py'),'freeze_sha256':sha(w/'freeze.json'),'max_abs_delta_D':max(abs(x['delta_D']) for x in results),'max_abs_delta_data_Q':max(abs(x['delta_data_Q']) for x in results)}
 out=P/('endpoint-result.json' if stage=='endpoints' else 'corner-result.json');assert not out.exists();save(out,result);np.savez_compressed(w/'states.npz',**arrays);print(json.dumps({k:v for k,v in result.items() if k not in ['cases','native_state']},indent=2))


def prepare_corners():
 coords=source_gate();r=load(P/'endpoint-result.json');assert r['all_gates_pass'];ends={(c['coordinate_index'],c['side']):c for c in r['cases'] if c['kind']=='coordinate_endpoint'};cases=[];envelopes=[]
 for cid in map(str,range(1,9)):
  for stat,key in [('D','native_D'),('Q','native_data_Q')]:
   for direction in ['max','min']:
    changes=[];prediction=0.;terms=[]
    for j,c in enumerate(coords):
     if c['CID']!=cid or c['input_low']==c['input_high']:continue
     gradient=(ends[(j,'high')][key]-ends[(j,'low')][key])/(c['input_high']-c['input_low']);want_high=(gradient>=0) if direction=='max' else (gradient<0);value=c['input_high' if want_high else 'input_low'];shift=gradient*(value-c['baseline_input']);prediction+=shift;changes.append({'field':c['field'],'epoch':c['epoch'],'value':value});terms.append({'coordinate_index':j,'gradient':gradient,'linear_contribution':shift})
    cases.append({'new_CID':str(200001+len(cases)),'original_CID':cid,'kind':'directed_corner','statistic':stat,'direction':direction,'linear_prediction':prediction,'changes':changes});envelopes.append({'CID':cid,'statistic':stat,'direction':direction,'linear_prediction':prediction,'terms':terms})
 for cid in map(str,range(1,9)):
  cc=[c for c in coords if c['CID']==cid and c['field']=='MJD'];assert all(len(c['SIMLIB_candidates'])==1 for c in cc);cases.append({'new_CID':str(200001+len(cases)),'original_CID':cid,'kind':'candidate_SIMLIB_times','changes':[{'field':'MJD','epoch':c['epoch'],'value':float(c['SIMLIB_candidates'][0]['MJD'])} for c in cc]})
 assert len(cases)==40;save(P/'linear-envelopes.json',envelopes);make_work('corners',cases)

if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('action',choices=['prepare_endpoints','prepare_corners','run_endpoints','run_corners','analyze_endpoints','analyze_corners']);x=a.parse_args().action
 if x=='prepare_endpoints':prepare_endpoints()
 elif x=='prepare_corners':prepare_corners()
 elif x.startswith('run_'):run(x.split('_',1)[1])
 else:analyze(x.split('_',1)[1])
