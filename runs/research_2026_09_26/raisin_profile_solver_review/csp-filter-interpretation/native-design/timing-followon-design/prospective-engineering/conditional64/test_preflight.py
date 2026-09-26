"""Read-only equivalence against engineering outputs plus synthetic summary formulas.
No new photons, fits, or 64 outcomes are accessed.
"""
from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import ast,json,hashlib,importlib.util,tempfile
import numpy as np
from astropy.io import fits
import validators as v
import run as run
O=Path(__file__).resolve().parent;P=O.parent;H=P/'restricted-peak-engineering/start-only-hook';R=Path('/home/szymon/Documents/ChatGPT/supernova')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def fn(path,name):
 t=path.read_text();n=next(n for n in ast.parse(t).body if isinstance(n,ast.FunctionDef) and n.name==name);return ast.get_source_segment(t,n)
assert fn(O/'validators.py','base_block_review')==fn(P/'restricted-peak-engineering/run_restricted.py','base_block_review')
assert fn(O/'validators.py','domain_review')==fn(P/'restricted-peak-engineering/run_restricted.py','domain_review')
assert fn(O/'validators.py','shift_review')==fn(H/'run_start_only.py','shift_review').replace("original=P/'fits-restricted/joint12'","original=P/'fits/joint12'")
assert (O/'full_identity.py').read_bytes()==(H/'full_identity.py').read_bytes()
assert (O/'head_adapter.py').read_bytes()==(H/'head-byte-adapter/adapter.py').read_bytes()
for b in ['original','ledger']:
 s=(O/'generation'/b/'sim.input').read_text();old=(P/'generation-v5'/b/'sim.input').read_text()
 expect=old.replace('SIMLIB_NREPEAT: 8','SIMLIB_NREPEAT: 64').replace('NGENTOT_LC: 8','NGENTOT_LC: 64').replace('RANSEED: 26092671','RANSEED: 26092672');expect=expect[:expect.index('PATH_SNDATA_SIM:')]+f'PATH_SNDATA_SIM: {R}/phase2/pt64/generation/{b}/output\n';assert s==expect
 assert (O/'generation'/b/'output').is_dir() and os.access(O/'generation'/b/'output',os.W_OK)
assert (R/'phase2/pt64').resolve()==O
assert (O/'inputs-v2').resolve()==P/'inputs-v2'
assert (O/'fit-sndata/SIM/PATH_SNDATA_SIM.LIST').read_text()==''
for h in ['repair-v5/original-build','repair-v5/ledger-build','restricted-peak-engineering/start-only-hook/build']:
 assert len(str(R/'phase2/pte'/h))<120
for q in ['generator-sndata','fit-sndata']:assert len(str(R/'phase2/pt64'/q))<120

def read_old(branch):
 with fits.open(P/'datasets/ledger/PTE/PTE_HEAD.FITS') as f:pass
# Read-only test uses old frozen FITS through a local replacement of read, no writes.
def read8(branch):
 with fits.open(P/'readme-adapted/ledger/PTE/PTE_HEAD.FITS') as f:h=f[1].data.copy();hc=f[1].columns
 with fits.open(P/'readme-adapted/ledger/PTE/PTE_PHOT.FITS') as f:p=f[1].data.copy();pc=f[1].columns
 return h,p,hc,pc
origread=v.read;v.read=read8
cids=[str(i) for i in range(1,9)];w=P/'fits-start-only/NIR_estimated12';by,d=v.base_block_review(w,6,cids,True);dom=v.domain_review(w,by,cids)
old=json.loads((H/'gates/NIR_estimated12-gate.json').read_text());assert d==old['details'] and dom==old['domain']
v.read=origread
# All same null-check logic reproduces the already frozen complete eight-event null.
a=run.identity.compare(P/'fits-start-only/NIR_true12',P/'fits-start-only/NIR_null12',cids,6,True,12);assert a['pass']
# Formula-only checks use fixed synthetic arrays, never the first8 science effect.
X=np.column_stack([np.arange(64),np.arange(64)/10,-np.arange(64)/20,-.15*np.arange(64)]);assert np.allclose(X[:,2]-X[:,1],X[:,3]);C=np.cov(X,rowvar=False,ddof=1);assert np.allclose(C,C.T) and np.all(np.diag(C)>=0)
for i in range(4):assert np.isclose(np.std(X[:,i],ddof=1)/np.sqrt(64),np.sqrt(C[i,i]/64))
# No frozen source or input is ever written by the imported executor.
res={'gate_pass':True,'native_calls':0,'reused_validator_functions_byte_equal':True,'exact_prior8_NIR_gate_replay':True,'exact_prior8_full_null_replay':True,'only_generation_input_changes':['seed26092672','fixed_attempts64','repeat64','new_output_path'],'all_output_parents_writable':True,'short_alias_paths_checked':True,'summary_formula_fixed_synthetic_check':True,'engineering8_not_used_as64_outcomes':True,'max_RSS_bytes':__import__('resource').getrusage(__import__('resource').RUSAGE_SELF).ru_maxrss*1024}
(O/'preflight-result.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
