from pathlib import Path
import json,hashlib,ast,shutil,difflib
R=Path('/home/szymon/Documents/ChatGPT/supernova');P=(R/'phase2/pte').resolve();H=P/'restricted-peak-engineering/start-only-hook';O=P/'conditional64';O.mkdir(exist_ok=False);S=R/'phase2/pt64';S.symlink_to(O)
def save(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def extract(path,names):
 t=path.read_text();a=ast.parse(t);return '\n\n'.join(ast.get_source_segment(t,n) for n in a.body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef)) and n.name in names)+'\n'
(O/'inputs-v2').symlink_to(P/'inputs-v2')
for mode in ['generator-sndata','fit-sndata']:
 d=O/mode;d.mkdir()
 for f in (P/'fit-private-lookup/SNDATA_ROOT').iterdir():
  if f.name!='SIM':(d/f.name).symlink_to(f.resolve())
 (d/'SIM').mkdir();(d/'SIM/PATH_SNDATA_SIM.LIST').write_text('')
for branch in ['original','ledger']:
 w=O/'generation'/branch;w.mkdir(parents=True);(w/'output').mkdir()
 s=(P/'generation-v5'/branch/'sim.input').read_text();s=s.replace('SIMLIB_NREPEAT: 8','SIMLIB_NREPEAT: 64').replace('NGENTOT_LC: 8','NGENTOT_LC: 64').replace('RANSEED: 26092671','RANSEED: 26092672')
 s=s[:s.index('PATH_SNDATA_SIM:')]+f'PATH_SNDATA_SIM: {S}/generation/{branch}/output\n';(w/'sim.input').write_text(s)
 (O/(branch+'-input.patch')).write_text(''.join(difflib.unified_diff((P/'generation-v5'/branch/'sim.input').read_text().splitlines(True),s.splitlines(True),fromfile='engineering8/'+branch,tofile='conditional64/'+branch)))
shutil.copy2(H/'schema.json',O/'schema.json');shutil.copy2(P/'restricted-peak-engineering/build-protocol.json',O/'build-protocol.json')
shutil.copy2(H/'full_identity.py',O/'full_identity.py');shutil.copy2(H/'head-byte-adapter/adapter.py',O/'head_adapter.py')
old=P/'restricted-peak-engineering/run_restricted.py'
header='''from pathlib import Path
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
def save(p,d):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\\n')
def read(branch):
 with fits.open(P/'datasets/ledger/PTE/PTE_HEAD.FITS') as f:h=f[1].data.copy();hc=f[1].columns
 with fits.open(P/'datasets/ledger/PTE/PTE_PHOT.FITS') as f:p=f[1].data.copy();pc=f[1].columns
 return h,p,hc,pc
'''
s=header+extract(old,['base_block_review','domain_review'])+extract(H/'run_start_only.py',['records','shift_review'])
s=s.replace("original=P/'fits-restricted/joint12'","original=P/'fits/joint12'")
(O/'validators.py').write_text(s)
save(O/'reused-source-functions.json',{'base_block_review':str(old.relative_to(R)),'domain_review':str(old.relative_to(R)),'records':str((H/'run_start_only.py').relative_to(R)),'shift_review':'same source except reference directory fits/joint12; every threshold unchanged','parser':'exact original imported parser, no compression or arithmetic changes','full_identity':'byte-identical copied start-only full_identity.py','head_adapter':'byte-identical copied passing exact-byte adapter.py'})
shutil.copy2('/tmp/prepare_conditional64.py',O/'prepare.py')
print(O)
