"""Install versioned public likelihood assets locally and verify immutable trees."""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
from pathlib import Path
from adapter import ROOT, WORK, PACKAGES

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / 'results/external_probes'
COMPONENTS = ['planck_2018_highl_plik.TTTEEE', 'planck_2018_lowl.TT',
              'planck_2018_lowl.EE', 'planck_2018_lensing.native',
              'planck_2018_highl_plik.TTTEEE_lite_native', 'bao.desi_dr2']
ORIGINS = {
 'planck_2018': 'https://pla.esac.esa.int/pla-sl/data-action?COSMOLOGY.COSMOLOGY_OID=151902',
 'planck_2018_lowT_native': 'https://github.com/CobayaSampler/planck_native_data/releases/download/v1/planck_2018_lowT.zip',
 'planck_2018_lowE_native': 'https://github.com/CobayaSampler/planck_native_data/releases/download/v1/planck_2018_lowE.zip',
 'planck_2018_pliklite_native': 'https://github.com/CobayaSampler/planck_native_data/releases/download/v1/plik_lite_2018_AL.zip',
 'planck_supp_data_and_covmats': 'https://github.com/CobayaSampler/planck_supp_data_and_covmats/archive/refs/tags/v2.1.tar.gz',
 'bao_data': 'https://github.com/CobayaSampler/bao_data/archive/refs/tags/v2.6.tar.gz',
}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def tree(path):
    inventory = {str(p.relative_to(path)): {'sha256': sha(p), 'bytes': p.stat().st_size}
                 for p in sorted(Path(path).rglob('*'))
                 if p.is_file() and '__pycache__' not in p.parts and p.suffix != '.pyc'}
    encoded = json.dumps(inventory, sort_keys=True, separators=(',', ':')).encode()
    return {'files': len(inventory), 'bytes': sum(v['bytes'] for v in inventory.values()),
            'tree_sha256': hashlib.sha256(encoded).hexdigest()}, inventory


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--verify-only', action='store_true')
    a = p.parse_args()
    for name, version in [('cobaya','3.6.2'), ('camb','1.6.6'), ('clipy-like','0.15')]:
        assert importlib.metadata.version(name) == version, name
    if not a.verify_only:
        from cobaya.install import install
        install(*COMPONENTS, path=str(PACKAGES), no_progress_bars=True,
                no_set_global=True, skip_global=True)
    inventories = {}
    assets = {}
    for name, url in ORIGINS.items():
        rec, inv = tree(PACKAGES / 'data' / name)
        assert rec['files'] > 0, name
        assets[name] = {'url': url, **rec}
        inventories[name] = inv
    code, inv = tree(PACKAGES/'code/planck/clipy/clipy')
    assets['clipy_source'] = {
        'url': 'https://github.com/benabed/clipy/archive/refs/tags/clipy_0.15.tar.gz', **code}
    inventories['clipy_source'] = inv
    RESULTS.mkdir(parents=True, exist_ok=True)
    dest = RESULTS/'acquisition.json'
    old = None
    if dest.exists():
        old = json.loads(dest.read_text())
        assert old['assets'] == assets, 'Acquired likelihood bytes changed; do not silently refresh their scientific identity.'
    inventory_file = WORK/'asset-file-inventory.json'
    inventory_file.write_text(json.dumps(inventories, indent=2)+'\n')
    versions = {p: importlib.metadata.version(p) for p in ['cobaya','camb','clipy-like','numpy','scipy','astropy']}
    out = {'status':'verified', 'created_utc':old['created_utc'] if old else datetime.now(timezone.utc).isoformat(),
           'assets':assets, 'versions':versions, 'code_sha256':sha(__file__),
           'requirements_lock_sha256':sha(HERE/'requirements-lock.txt'),
           'file_inventory':str(inventory_file.relative_to(ROOT)),
           'file_inventory_sha256':sha(inventory_file),
           'scope':'Public likelihood data and third-party code are downloaded into ignored .work; tree identity includes every acquired data file.'}
    dest.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'assets':len(assets),'data_files':sum(x['files'] for x in assets.values()),'versions':versions}))


if __name__ == '__main__':
    main()
