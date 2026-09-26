"""Root-released five-job processing test; immutable source, no branch selection."""
from pathlib import Path
import sys,json,hashlib,os,subprocess,time,argparse
from collections import Counter,defaultdict
import numpy as np
ROOT=Path.cwd();O=Path(__file__).resolve().parent;I=O.parent/'instrumentation';V=O.parent/'convergence-diagnostic';sys.path.insert(0,str(I));from parse_native import parse
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def groups(bs):
 out={}
 for b in bs:out.setdefault(b['CID'],[]).append(b)
 return out

def exact(a,b):return a['objective']==b['objective'] and a['entry']==b['entry'] and np.array_equal(a['array'],b['array']) and np.array_equal(a['W'],b['W']) and [(r['band'],r['source_epoch']) for r in a['rows']]==[(r['band'],r['source_epoch']) for r in b['rows']]
def physical(a,b,rawmap):
 old=Counter((r['MJD'],r['dataF'],r['data_error'],r['band']) for r in a['rows']);expect=Counter()
 for key,n in old.items():
  targets=rawmap[(a['CID'],key)];assert len(targets)>=n,(a['CID'],key,n)
  if len(set(targets))==1:counts=Counter({targets[0]:n})
  elif len(targets)==n:counts=Counter(targets)
  else:raise RuntimeError('Partial ambiguousduplicate mapping; no arbitrary assignment')
  for band,c in counts.items():expect[(*key[:3],band)]+=c
 new=Counter((r['MJD'],r['dataF'],r['data_error'],r['band']) for r in b['rows']);return expect==new

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--execute',action='store_true');a=ap.parse_args();assert a.execute,'Requires root release and explicit --execute'
 p=json.loads((O/'protocol.json').read_text())
 for n,h in p['inputs_sha256'].items():assert sha(ROOT/n)==h,n
 nominal=groups(parse(V/'fits/iter12_default/fit.log'));assert set(nominal)==set(p['membership'])
 nominal_shifted={n:groups(parse(V/'fits'/n/'fit.log')) for n in ['iter12_minus','iter12_plus']}
 rawmap=defaultdict(list)
 for r in json.loads((O/'input-map.json').read_text()):rawmap[(r['CID'],tuple(r['key_MJD_dataF_dataE_band']))].append(r['new_band'])
 binary=I/'SNANA-v11_04k-output/bin/snlc_fit.exe';assert Path('/tmp/snana-csp-audit-8c5f0d').resolve()==binary.parents[1]
 env=os.environ.copy();env.update(SNANA_DIR='/tmp/snana-csp-audit-8c5f0d',SNDATA_ROOT=str(ROOT/'phase2/official/inputs/SNDATA_ROOT'),LD_LIBRARY_PATH=str(ROOT/'phase2/official/build/sysroot/usr/lib'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1');env.pop('CSP_POSTINIT_DSHIFT',None)
 used=0;executions=[]
 for j in p['jobs']:
  d=O/'fits'/j['name'];assert not(d/'fit.log').exists();e=env.copy()
  if j['CSP_POSTINIT_DSHIFT'] is not None:e['CSP_POSTINIT_DSHIFT']=str(j['CSP_POSTINIT_DSHIFT'])
  t=time.monotonic()
  with (d/'fit.log').open('x') as f:r=subprocess.run([str(binary),'fit.nml'],cwd=d,env=e,stdout=f,stderr=subprocess.STDOUT,timeout=min(60,120-used))
  dt=time.monotonic()-t;used+=dt;x={'job':j['name'],'returncode':r.returncode,'seconds':dt,'cumulative_new_seconds':used,'overall_native_seconds':p['resource']['carried_native_seconds']+used,'protocol_sha256':sha(O/'protocol.json'),'log_sha256':sha(d/'fit.log')};executions.append(x);dump(d/'execution.json',x);dump(O/'execution-ledger.json',{'jobs':executions,'additional_seconds':used,'overall_seconds':p['resource']['carried_native_seconds']+used});print(x,flush=True)
  if r.returncode:raise RuntimeError('Nativefailure retained; stop')
  current=groups(parse(d/'fit.log'));assert set(current)==set(p['membership']);assert all(len(bs)==j['iterations'] for bs in current.values())
  columns=None;fr=[]
  for line in (d/'fit.FITRES.TEXT').read_text().splitlines():
   t=line.split()
   if t and t[0]=='VARNAMES:':columns=t[1:]
   elif t and t[0]=='SN:':fr.append(dict(zip(columns,t[1:],strict=True)))
  assert len(fr)==42 and {r['CID'] for r in fr}==set(p['membership'])
  assert all(int(r['ERRFLAG_FIT'])==0 and float(r['STRETCH'])==1 and float(r['AV'])==0 and all(float(r[k])==0 for k in ['STRETCHERR','AVERR','PKMJDERR']) for r in fr)
  ledger=[]
  for cid,bs in current.items():
   for index,b in enumerate(bs):
    n=nominal[cid][index];rows=physical(n,b,rawmap) if j['arm']=='known_Jdw_to_j' else [(r['MJD'],r['dataF'],r['data_error'],r['band']) for r in n['rows']]==[(r['MJD'],r['dataF'],r['data_error'],r['band']) for r in b['rows']]
    # Shifted controls compare with correspondingly shifted nominal run.
    control_ref=n
    if j['CSP_POSTINIT_DSHIFT'] is not None:
     name='iter12_minus' if j['CSP_POSTINIT_DSHIFT']<0 else 'iter12_plus'
     if cid in p['controls_10']:control_ref=nominal_shifted[name][cid][index]
    control=exact(control_ref,b) if cid in p['controls_10'] or j['arm']=='nominal_copy' else True
    ledger.append({'CID':cid,'ITER':b['ITER'],'physical_rows':rows,'control_exact':control})
  dump(d/'structural-gates.json',{'pass':all(v['physical_rows'] and v['control_exact'] for v in ledger),'checks':ledger});assert all(v['physical_rows'] and v['control_exact'] for v in ledger)
  if j['name']=='changed12_default':
   nine=groups(parse(O/'fits/changed09_default/fit.log'));ok=all(exact(a,b) for cid in current for a,b in zip(nine[cid],current[cid][:9],strict=True));dump(d/'prefix-gate.json',{'changed09_to_changed12_exact':ok});assert ok
  if j['name']=='nominal12_copy':
   science=lambda q:[l for l in q.read_text().splitlines() if l.startswith(('SN:','VARNAMES:'))]
   ok=science(d/'fit.FITRES.TEXT')==science(V/'fits/iter12_default/fit.FITRES.TEXT') and (d/'fit.LCPLOT.TEXT').read_bytes()==(V/'fits/iter12_default/fit.LCPLOT.TEXT').read_bytes();dump(d/'copy-gate.json',{'science_and_LCPLOT_exact':ok});assert ok
 print('Five jobs complete. Numericalgates and certifiedresponse stillrequire frozenpostprocessing.',flush=True)
if __name__=='__main__':main()
