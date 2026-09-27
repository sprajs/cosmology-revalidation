"""Restore/verify frozen modern assets without rewriting their design identity.

modern_acquire.py records initial acquisition. This repeat-use entrypoint
preserves that immutable record, including the pre-training dependency lock.
"""
import argparse,importlib,json
from pathlib import Path
from modern_acquire import REPOS
from acquire import ROOT,RESULTS,WORK,PACKAGES,sha,tree


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--verify-only',action='store_true');args=parser.parse_args()
    manifest=RESULTS/'modern-acquisition.json';expected=json.loads(manifest.read_text())
    if not args.verify_only:
        from cobaya.install import install
        install('act_dr6_cmbonly','act_dr6_cmbonly.PlanckActCut','act_dr6_lenslike.ACTDR6LensLike',
                path=str(PACKAGES),no_progress_bars=True,no_set_global=True)
    paths={'ACTDR6_CMB':PACKAGES/'data/ACTDR6CMBonly',
           'ACT_Planck_lensing':PACKAGES/'data/ACT_dr6_likelihood/v1.2'}
    paths.update({name:Path(importlib.import_module(name).__file__).parent for name in REPOS})
    verified={}
    for name,path in paths.items():
        summary,_=tree(path);original={k:expected['assets'][name][k] for k in summary}
        assert summary==original,(name,summary,original)
        verified[name]=summary
    # The original complete selected-runtime lock is also immutable.
    assert sha(Path(__file__).with_name('modern-requirements-lock.txt'))==expected['requirements_lock_sha256']
    out={'status':'passed','frozen_acquisition_sha256':sha(manifest),'verified_assets':verified,
         'code_sha256':sha(__file__),'manifest_rewritten':False,
         'scope':'Restore public data and verify installed selected scientific libraries; supplemental sampling packages do not alter these trees.'}
    (RESULTS/'modern-restoration.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'status':'passed','verified_assets':len(verified),'manifest_rewritten':False}))

if __name__=='__main__':main()
