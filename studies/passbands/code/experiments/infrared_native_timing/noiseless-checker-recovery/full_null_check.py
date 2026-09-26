"""Read-only exact saved-state identity. No executor/parser import, native calls or fitting."""
from pathlib import Path
import argparse,collections,hashlib,json,re
import numpy as np
from astropy.io import fits

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def records(p):
 out=[];blocks=[];b=None;entry=None
 def done():
  if b is None:return
  n=b['n'];assert len(b['rows'])==len(b['w'])==n
  assert [int(x[3]) for x in b['rows']]==list(range(1,n+1))
  assert [int(x[3]) for x in b['w']]==list(range(1,n+1))
  rows=np.array([[float(v) for v in x[6:]] for x in b['rows']]);assert rows.shape==(n,13)
  w=np.array([[float(v) for v in x[4:]] for x in b['w']]);w=w if b['cov'] else np.diag(w[:,0]);assert w.shape==(n,n)
  assert np.isfinite(rows).all() and np.isfinite(w).all();assert np.linalg.eigvalsh(w).min()>0
  r=rows[:,4]-rows[:,2];o=b['objective'];q=float(r@w@r)+float(o[6])+float(o[7]);b['qerr']=abs(q-float(o[5]));assert b['qerr']<1e-7
  blocks.append(b)
 for line in p.read_text().splitlines():
  t=line.split()
  if not t:continue
  if t[0].startswith('CSP_'):
   out.append(t)
   if t[0]=='CSP_ENTRY:':assert len(t)==18;entry=t
   elif t[0]=='CSP_START_SHIFT:':assert len(t)==6
   elif t[0]=='CSP_ROW:':
    assert len(t)==19
    if b is None or 'objective' in b:
     done();b={'cid':t[1],'iter':int(t[2]),'rows':[],'w':[],'entry':entry}
    assert t[1]==b['cid'] and int(t[2])==b['iter'];b['rows'].append(t)
   elif t[0]=='CSP_OBJECTIVE:':assert len(t)==15;b.update(n=int(t[3]),cov=t[4]=='T',objective=t)
   elif t[0] in ['CSP_WROW:','CSP_WDIAG:']:assert t[1]==b['cid'] and int(t[2])==b['iter'];b['w'].append(t)
   else:raise AssertionError(t[0])
 done();assert blocks
 return out,blocks

def compare(left,right,cids,n,support_both=True):
 l,lb=records(left/'native.log');r,rb=records(right/'native.log');assert l==r,'CSP model/measurement/entry/objective/weight records differ'
 assert {b['cid'] for b in lb}==set(cids)
 for cid in cids:
  bb=[b for b in lb if b['cid']==cid];assert {b['iter'] for b in bb}==set(range(1,13));assert all(b['n']==n for b in bb)
 def science(p):
  lines=[x for x in p.read_text().splitlines() if x.startswith(('VARNAMES:','SN:'))];assert lines and lines[0].startswith('VARNAMES:')
  names=lines[0].split()[1:];rows=[dict(zip(names,x.split()[1:])) for x in lines[1:]];assert {x['CID'] for x in rows}==set(cids) and len(rows)==len(cids)
  assert all(int(x['ERRFLAG_FIT'])==0 and int(x['CUTFLAG_SNANA'])==3 for x in rows)
  return lines
 assert science(left/'fit.FITRES.TEXT')==science(right/'fit.FITRES.TEXT'),'full FITRES science rows differ'
 if support_both:
  ss=[]
  for p in [left,right]:
   lines=[x.split() for x in (p/'native.log').read_text().splitlines() if x.startswith('PROSP_SUPPORT')]
   assert len(lines)==len(cids) and {x[1] for x in lines}==set(cids)
   assert all(x[0]=='PROSP_SUPPORT' and len(x)==11 and int(x[2])>0 and int(x[3])==int(x[4])==0 for x in lines)
   ss.append(lines)
  assert ss[0]==ss[1],'all-call support summary differs'
 def tables(p):
  m=re.search(r"PRIVATE_DATA_PATH\s*=\s*'([^']+)'",(p/'fit.nml').read_text());assert m
  d=p/m.group(1)/'PTE';h=list(d.glob('*HEAD.FITS*'));q=list(d.glob('*PHOT.FITS*'));assert len(h)==len(q)==1
  with fits.open(h[0]) as f:head=f[1].data.copy()
  with fits.open(q[0]) as f:phot=f[1].data.copy()
  return h[0],q[0],head,phot
 lh,lp,hh,pp=tables(left);rh,rp,h2,p2=tables(right)
 for a,b in [(hh,h2),(pp,p2)]:
  assert a.dtype==b.dtype and a.shape==b.shape
  for k in a.dtype.names:assert np.array_equal(a[k],b[k],equal_nan=True) if a[k].dtype.kind in 'fc' else np.array_equal(a[k],b[k])
 assert lp.resolve()==rp.resolve() and sha(lp)==sha(rp),'null must use the same PHOT file'
 return {'pass':True,'left':str(left),'right':str(right),'expected_cids':cids,'rows_per_callback':n,'callbacks':len(lb),'exact_CSP_token_records':len(l),'full_model_mean_error_phase_nuisance_entry_and_W_equal':True,'C_equal_by_identical_invertible_W':True,'full_total_prior_sigma_objectives_equal':True,'all_FITRES_science_rows_equal':True,'HEAD_all_columns_equal':True,'PHOT_same_file_and_hash':True,'allcall_support_equal_checked':support_both,'max_abs_Q_reconstruction':max(b['qerr'] for b in lb),'hashes':{str(p):sha(p) for p in [left/'native.log',right/'native.log',left/'fit.FITRES.TEXT',right/'fit.FITRES.TEXT',lh,rh,lp]}}
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('left',type=Path);a.add_argument('right',type=Path);a.add_argument('--cids',required=True);a.add_argument('--rows',type=int,required=True);a.add_argument('--output',type=Path,required=True);a.add_argument('--instrumentation-control',action='store_true');x=a.parse_args()
 z=compare(x.left,x.right,x.cids.split(','),x.rows,not x.instrumentation_control);x.output.write_text(json.dumps(z,indent=2)+'\n');print(json.dumps({k:v for k,v in z.items() if k!='hashes'},indent=2))
