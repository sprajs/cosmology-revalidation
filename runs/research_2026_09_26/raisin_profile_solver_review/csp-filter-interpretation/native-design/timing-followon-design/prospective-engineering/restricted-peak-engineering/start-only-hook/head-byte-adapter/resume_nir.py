"""Additive administrative resume: identical frozen NIR routine with raw HEAD adapter.
No native work on import. Execution requires a new root release.
"""
from pathlib import Path
import argparse,importlib.util,json,hashlib
from adapter import patch_head
O=Path(__file__).resolve().parent;H=O.parent;R=Path('/home/szymon/Documents/ChatGPT/supernova')
s=importlib.util.spec_from_file_location('frozen_start_executor',H/'run_start_only.py');e=importlib.util.module_from_spec(s);s.loader.exec_module(e)
sha=e.sha;save=e.save;check=e.check

def adapt_head(kind,peaks):
 check(kind in ['estimated','null'],'unknown adapter kind')
 out=e.P/'derived-start-only-bytes'/kind;out.mkdir(parents=True,exist_ok=False)
 src,h,p=e.files('ledger');dst=out/'PTE';dst.mkdir()
 for f in src.iterdir():
  if f==h:continue
  if f.is_file():(dst/f.name).symlink_to(f.resolve())
 result=patch_head(h,dst/h.name,peaks,require_null=kind=='null')
 check((dst/p.name).resolve()==p.resolve(),'PHOT link changed')
 save(out/'adapter-ledger.json',{'source_PHOT_sha256':sha(p),'source_PHOT_same_file':True,'changed_column_only':'PEAKMJD','byte_patch':result})
 return '../../derived-start-only-bytes/'+kind

def main(release):
 d=json.loads(release.read_text())
 check(d['protocol_sha256']==sha(O/'protocol.json') and d['freeze_sha256']==sha(O/'freeze.json'),'additive release hash')
 check(d['stages']==['nir'],'only NIR resume released')
 for f,h in json.loads((O/'freeze.json').read_text())['files'].items():check(sha(R/f)==h,'additive hash changed '+f)
 e.release('nir',H/'root-release-nir.json')
 carried=json.loads((O/'carried-native-activity.json').read_text());now=json.loads((H/'native-activity.json').read_text())
 check(now==carried,'unexpected intervening native activity')
 check(not (e.P/'derived-start-only-bytes').exists(),'byte adapter already attempted')
 check(not list((e.P/'fits-start-only').glob('NIR_*')),'NIR fit already attempted')
 e.adapt_head=adapt_head
 try:
  e.run_nir()
  save(O/'resume-result.json',{'gate_pass':True,'new_root_release_sha256':sha(release),'existing_engineering_result_sha256':sha(H/'engineering-result.json'),'original_NIR_routine_unchanged':True})
 except Exception as err:
  save(O/'resume-failure.json',{'exception':type(err).__name__,'message':str(err),'release_sha256':sha(release),'protocol_sha256':sha(O/'protocol.json')});raise

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--release',type=Path,required=True);a=ap.parse_args();main(a.release)
