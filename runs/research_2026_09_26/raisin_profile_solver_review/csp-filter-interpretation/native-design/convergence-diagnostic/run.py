"""Run only after root releases frozen protocol. No scientific option switches."""
from pathlib import Path
import argparse,json,hashlib,os,subprocess,time,sys
import numpy as np
ROOT=Path.cwd();O=Path(__file__).resolve().parent;I=O.parent/'instrumentation';sys.path.insert(0,str(I));from parse_native import parse
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def groups(bs):
 d={}
 for b in bs:d.setdefault(b['CID'],[]).append(b)
 return d
def prefix(a,b,n):
 aa,bb=groups(a),groups(b);assert aa.keys()==bb.keys();checks=[]
 for cid in aa:
  for x,y in zip(aa[cid][:n],bb[cid][:n],strict=True):
   checks.append({'CID':cid,'ITER':x['ITER'],'equal':x['objective']==y['objective'] and x['entry']==y['entry'] and np.array_equal(x['array'],y['array']) and np.array_equal(x['W'],y['W']) and [(r['band'],r['source_epoch']) for r in x['rows']]==[(r['band'],r['source_epoch']) for r in y['rows']]})
 return {'pass':all(x['equal'] for x in checks),'checks':checks}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--execute',action='store_true');args=ap.parse_args();assert args.execute,'Execution requires explicit root release and --execute.'
 p=json.loads((O/'protocol.json').read_text());a=json.loads((O/'pre-execution-amendment.json').read_text())
 for n,h in {**p['inputs_sha256'],**a['additional_sha256']}.items():assert sha(ROOT/n)==h,n
 binary=I/'SNANA-v11_04k-output/bin/snlc_fit.exe';assert Path('/tmp/snana-csp-audit-8c5f0d').resolve()==binary.parents[1]
 env=os.environ.copy();env.update(SNANA_DIR='/tmp/snana-csp-audit-8c5f0d',SNDATA_ROOT=str(ROOT/'phase2/official/inputs/SNDATA_ROOT'),LD_LIBRARY_PATH=str(ROOT/'phase2/official/build/sysroot/usr/lib'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1');env.pop('CSP_POSTINIT_DSHIFT',None)
 used=0;execs=[];old=parse(ROOT/'runs/research_2026_09_26/csp_native_filter_response/cohort-execution/full/fits/nominal/fit.log');expected=set(groups(old));assert len(expected)==42
 for j in p['jobs']:
  d=O/'fits'/j['name'];assert not(d/'fit.log').exists();e=env.copy()
  if j['CSP_POSTINIT_DSHIFT'] is not None:e['CSP_POSTINIT_DSHIFT']=str(j['CSP_POSTINIT_DSHIFT'])
  t=time.monotonic()
  with (d/'fit.log').open('x') as f:r=subprocess.run([str(binary),'fit.nml'],cwd=d,env=e,stdout=f,stderr=subprocess.STDOUT,timeout=min(60,120-used))
  dt=time.monotonic()-t;used+=dt;x={'job':j['name'],'returncode':r.returncode,'seconds':dt,'cumulative_new_seconds':used,'overall_carried_seconds':a['carried_native_seconds']+used,'log_sha256':sha(d/'fit.log'),'protocol_sha256':sha(O/'protocol.json'),'amendment_sha256':sha(O/'pre-execution-amendment.json')};execs.append(x);dump(d/'execution.json',x);dump(O/'execution-ledger.json',{'jobs':execs,'additional_seconds':used,'overall_seconds':a['carried_native_seconds']+used});print(x,flush=True)
  if r.returncode:raise RuntimeError('Native failed; preserve output and stop.')
  bs=parse(d/'fit.log');assert set(groups(bs))==expected
  cols=None;fr=[]
  for l in (d/'fit.FITRES.TEXT').read_text().splitlines():
   t=l.split()
   if t and t[0]=='VARNAMES:':cols=t[1:]
   elif t and t[0]=='SN:':fr.append(dict(zip(cols,t[1:],strict=True)))
  assert {v['CID'] for v in fr}==expected and len(fr)==42
  assert all(int(v['ERRFLAG_FIT'])==0 and float(v['STRETCH'])==1 and float(v['AV'])==0 and float(v['STRETCHERR'])==0 and float(v['AVERR'])==0 and float(v['PKMJDERR'])==0 for v in fr)
  if j['name']=='iter09_default':z=prefix(old,bs,3);dump(d/'prefix-gate.json',z);assert z['pass'],'Original3 prefix differs; stop before interpretation.'
  elif j['name']=='iter12_default':z=prefix(parse(O/'fits/iter09_default/fit.log'),bs,9);dump(d/'prefix-gate.json',z);assert z['pass'],'Nine-step prefix differs; stop before interpretation.'
 print('Native stage done; all outputs retained. Numerical and support certification requires separate frozen audit.',flush=True)
if __name__=='__main__':main()
