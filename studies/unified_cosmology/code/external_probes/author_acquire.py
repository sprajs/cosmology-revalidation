"""Acquire the source/header reconstruction assets without redistributing them."""
import argparse,json,tarfile
from pathlib import Path
from author_config import BASE,CSL,RELEASE,fetch
from acquire import RESULTS,ROOT,sha,tree


def main():
    p=argparse.ArgumentParser();p.add_argument('--verify-only',action='store_true');a=p.parse_args()
    urls={}
    for name,path in {
        'act_dr6_lite_interface.py':'likelihood/act-dr6-lite/act_dr6_lite_interface.py',
        'candl_cosmosis_interface.py':'likelihood/candl/candl_cosmosis_interface.py',
        'planck_py_interface.py':'likelihood/planck_py/planck_py_interface.py',
        'planck_lite_py.py':'likelihood/planck_py/planck_lite_py.py'}.items():
        urls[f'modules/{name}']=f'https://raw.githubusercontent.com/joezuntz/cosmosis-standard-library/{CSL}/{path}'
    cslpaths=['boltzmann/camb/camb_interface.py','utility/exclude_w0_wa/w0wa_sum_prior.py',
        'likelihood/act-dr6-lens/act_dr6_lenslike_interface.py','likelihood/act-dr6-lens/get-act-data.sh',
        'utility/consistency/consistency.py','utility/consistency/consistency_interface.py','likelihood/bao/desi-dr2/desi_dr2.py']
    for folder in ['2018','2018_low_ell']:
        for suffix in ['.dat','_cov.npz','_weights.dat']:
            cslpaths.append(f'likelihood/planck_py/data/{folder}/planck_lite_2018_v22{suffix}')
    for path in cslpaths:urls[f'CSL/{path}']=f'https://raw.githubusercontent.com/joezuntz/cosmosis-standard-library/{CSL}/{path}'
    urls.update({
        'Dovekie_cosmosis_likelihood.py':f'https://raw.githubusercontent.com/des-science/DES-SN5YR/{RELEASE}/5_COSMOLOGY/Dovekie_cosmosis_likelihood.py',
        'SPT3G_D1_TnE_lite-historical.yaml':'https://raw.githubusercontent.com/SouthPoleTelescope/spt_candl_data/68d690697344cd533087af821f025f23dc6c5e52/spt_candl_data/SPT3G_D1_TnE_v0/lite/SPT3G_D1_TnE_lite.yaml',
        'gaussian_likelihood.py':'https://raw.githubusercontent.com/joezuntz/cosmosis/238df945450ceb56241a89569f1092956d3df5d6/cosmosis/gaussian_likelihood.py',
        'ACT_dr6_likelihood_v1.1.tgz':'https://lambda.gsfc.nasa.gov/data/suborbital/ACT/ACT_dr6/likelihood/data/ACT_dr6_likelihood_v1.1.tgz'})
    # Keep the exact public core revision separate from the unrecorded historical installation.
    dest=RESULTS/'author-acquisition.json';old=json.loads(dest.read_text()) if dest.exists() else None
    assets={}
    for rel,url in urls.items():
        path=BASE/rel
        if a.verify_only:assert path.exists(),rel
        rec=fetch(url,path)
        if old:assert old['assets'][rel]==rec,rel
        assets[rel]=rec
    if not (BASE/'act-lensing/v1.1').exists():
        with tarfile.open(BASE/'ACT_dr6_likelihood_v1.1.tgz') as tar:tar.extractall(BASE/'act-lensing',filter='data')
    summary,inventory=tree(BASE/'act-lensing/v1.1');invfile=BASE/'lensing-v1.1-inventory.json'
    invfile.write_text(json.dumps(inventory,indent=2)+'\n')
    if old:assert old['ACT_lensing_extracted_tree']==summary
    import spt_candl_data,yaml
    historical=yaml.safe_load((BASE/'SPT3G_D1_TnE_lite-historical.yaml').read_text())
    current=yaml.safe_load(Path(spt_candl_data.SPT3G_D1_TnE_lite).read_text())
    out={'status':'verified','assets':assets,'ACT_lensing_extracted_tree':summary,
        'lensing_inventory_path':str(invfile.relative_to(ROOT)),'lensing_inventory_sha256':sha(invfile),
        'historical_SPT_lite_priors':historical['priors'],'historical_current_SPT_yaml_equal':historical==current,
        'missing_runtime_inputs':['Uncommitted local modifications, if any, at author cwd Git HEAD',
            'Original likelihood/dovekie-sn/dovekie_sn_likelihood.py absent at that public revision',
            'Author installed CAMB, ACT and candl versions and complete dependency lock'],
        'code_sha256':sha(__file__)}
    dest.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'assets':len(assets),'historical_current_SPT_yaml_equal':historical==current}))

if __name__=='__main__':main()
