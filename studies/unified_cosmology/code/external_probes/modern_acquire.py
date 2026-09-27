"""Install and hash-pin the contemporary public CMB assets.

First create .modern-venv and install modern-requirements-lock.txt. Repository
dependencies are pinned Git revisions. This script downloads real ACT data via
its public installer and verifies the exact data trees used by the adapter.
"""
import argparse
from datetime import datetime,timezone
import importlib
import importlib.metadata
import json
from pathlib import Path
from adapter import WORK,PACKAGES,ROOT
from acquire import sha,tree,HERE,RESULTS

REPOS={
 'act_dr6_cmbonly':('ACTCollaboration/DR6-ACT-lite','627aeafb88ae5ad1aa66b406bea2d65cfa66a27d'),
 'act_dr6_lenslike':('ACTCollaboration/act_dr6_lenslike','5d54308ec60acbf545fd4e46a10091cc6f8a2237'),
 'candl':('Lbalkenhol/candl','650db0a6a0a2febed1e57350f0b991b537a4d2f7'),
 'candl_data':('Lbalkenhol/candl_data','c4eb8fbc5b0eb259ba487494b3cff918eac5c47a'),
 'spt_candl_data':('SouthPoleTelescope/spt_candl_data','efe35b7bfd815a115d5123e210a6849a4175bc81')}

def main():
    p=argparse.ArgumentParser();p.add_argument('--verify-only',action='store_true');a=p.parse_args()
    if not a.verify_only:
        from cobaya.install import install
        install('act_dr6_cmbonly','act_dr6_cmbonly.PlanckActCut','act_dr6_lenslike.ACTDR6LensLike',
                path=str(PACKAGES),no_progress_bars=True,no_set_global=True)
    assets={};inventories={}
    for name,folder,url in [
        ('ACTDR6_CMB',PACKAGES/'data/ACTDR6CMBonly','https://lambda.gsfc.nasa.gov/data/act/pspipe/sacc_files/dr6_data_cmbonly.tar.gz'),
        ('ACT_Planck_lensing',PACKAGES/'data/ACT_dr6_likelihood/v1.2','https://lambda.gsfc.nasa.gov/data/suborbital/ACT/ACT_dr6/likelihood/data/ACT_dr6_likelihood_v1.2.tgz')]:
        rec,inv=tree(folder);assert rec['files']>0
        assets[name]={'url':url,**rec};inventories[name]=inv
    for module,(repo,revision) in REPOS.items():
        folder=Path(importlib.import_module(module).__file__).parent
        rec,inv=tree(folder)
        assets[module]={'url':f'https://github.com/{repo}/tree/{revision}',
                       'git_revision':revision,**rec}
        inventories[module]=inv
    outfile=RESULTS/'modern-acquisition.json'
    if outfile.exists():
        old=json.loads(outfile.read_text())
        assert old['assets']==assets,'Modern dependency/data identities changed; do not silently refresh.'
    invfile=WORK/'modern-file-inventory.json';invfile.write_text(json.dumps(inventories,indent=2)+'\n')
    rec={'status':'verified','created_utc':datetime.now(timezone.utc).isoformat(),'assets':assets,
         'versions':{n:importlib.metadata.version(n) for n in ['camb','cobaya','candl-like','candl-data','spt-candl-data','sacc','numpy','scipy']},
         'requirements_lock_sha256':sha(HERE/'modern-requirements-lock.txt'),
         'inventory_path':str(invfile.relative_to(ROOT)),'inventory_sha256':sha(invfile),
         'code_sha256':sha(__file__),
         'scope':'Installed libraries include additional unused data; adapter selects only SPT2025 TnE and SPT2023 lensing, with exact config recorded separately. ACT data are NASA real-data assets, never the bundled simulation fallback. Shared Planck/BAO assets are pinned in acquisition.json.'}
    outfile.write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))

if __name__=='__main__':main()
