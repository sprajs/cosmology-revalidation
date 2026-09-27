#!/usr/bin/env python3
"""Pin the released flux/calibration/model relationship without double correction."""
import json,hashlib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from astropy.io import fits
from common import ROOT,WORK,RESULTS,sha
from acquire import get

def main():
    tree=json.loads((WORK/'release-tree.json').read_text());commit=json.loads((WORK/'release-commit.json').read_text())['sha'];data=WORK/'SNDATA_ROOT_2026-04-10'
    selected=[r for r in tree['tree'] if r['type']=='blob' and r['path'].startswith('2_LCFIT_MODEL/SALT3.DOVEKIE/') and '/plots/' not in r['path']]
    def one(item):
        path=item['path'];r=get(f'https://raw.githubusercontent.com/des-science/DES-SN5YR/{commit}/{path}',WORK/'release'/path);b=Path(r['path']).read_bytes();assert hashlib.sha1(f'blob {len(b)}\0'.encode()+b).hexdigest()==item['sha']
        bundle=data/'models/SALT3/SALT3.DOVEKIE'/Path(path).name
        r['official_bundle_matches_git_bytes']=bundle.exists() and sha(bundle)==r['sha256'];r['bundle_path']=str(bundle.relative_to(ROOT));return r
    with ThreadPoolExecutor(max_workers=4) as pool:models=list(pool.map(one,selected))
    passbands=[]
    for p in sorted((data/'kcor/Dovekie').glob('calib_DES-SN5YR*.fits*')):
        with fits.open(p) as f:
            names=[h.name for h in f];record={'path':str(p.relative_to(ROOT)),'sha256':sha(p),'HDUs':names}
            record['extensions']=[{'name':h.name,'rows':len(h.data) if h.data is not None else 0,'columns':list(h.columns.names) if hasattr(h,'columns') else []} for h in f]
        passbands.append(record)
    # A source-level field registry: author calibration systematics enter one total
    # released-distance covariance, not repeated independent noise at each stage.
    result={'code_sha256':sha(__file__),'status':'passed_asset_identity_audit','release_commit':commit,'SALT3_nominal_files':models,'calibration_files':passbands,'flux_semantics':{'SMP':'FLUXCAL=10^[-.4(m−27.5)], with chromatic/DCR and hostsurfacebrightness error corrections already applied as documented in the author README.','calibration_offsets':'AB+Fragilistic offsets are not applied in the released SMP flux; use the matched calibration/model inputs rather than add the offsets twice.','Dovekie_released_distances':'Total covariance includes the released calibration, lightcurve, selection and contamination variations; do not append MUERR variance or reuse them as independent calibration likelihoods.','cross_epoch_covariance':'SMP release inspected here supplies epoch errors and flags but no dense cross-epoch covariance product.'},'scope':'Assetidentity and interpretation, not rawdetector or calibration retraining closure.'}
    assert all(r['official_bundle_matches_git_bytes'] for r in models if not r['path'].endswith('SALT3.INFO'))
    (RESULTS/'calibration-registry.json').write_text(json.dumps(result,indent=2)+'\n');print(len(models),len(passbands))
if __name__=='__main__':main()
