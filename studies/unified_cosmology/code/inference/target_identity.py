"""Bind a cosmological target to actual source, likelihood bytes and configuration."""
import hashlib
import importlib
import importlib.metadata
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
WORK=ROOT/'.work/unified-cosmology/external-probes'


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def canonical(value):
    if isinstance(value,dict):return {k:canonical(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [canonical(v) for v in value]
    if isinstance(value,type):return value.__module__+'.'+value.__qualname__
    if isinstance(value,Path):return str(value.resolve())
    return value


def identify(info,sample_file,model_file=None):
    packages=WORK/'packages'
    recorded={}
    for name,invfile in [('primary','asset-file-inventory.json'),('modern','modern-file-inventory.json')]:
        inventory=json.loads((WORK/invfile).read_text())
        for group,files in inventory.items():
            if name=='primary':
                folder=packages/('code/planck/clipy/clipy' if group=='clipy_source' else 'data/'+group)
            elif group=='ACTDR6_CMB':folder=packages/'data/ACTDR6CMBonly'
            elif group=='ACT_Planck_lensing':folder=packages/'data/ACT_dr6_likelihood/v1.2'
            else:folder=Path(importlib.import_module(group).__file__).parent
            for relative,expected in files.items():
                actual=digest(folder/relative)
                assert actual==expected['sha256'],f'Likelihood asset changed: {group}/{relative}'
            recorded[name+'/'+group]=hashlib.sha256(json.dumps(files,sort_keys=True).encode()).hexdigest()
    sources=[HERE/p for p in ['modern_run.py','likelihood.py','late_geometry.py','spectral_surrogate.py','target_identity.py']]
    sources += [HERE.parent/'external_probes'/p for p in ['adapter.py','modern_adapter.py']]
    record={'configuration':canonical(info),'assets':recorded,
        'source_sha256':{str(p.relative_to(ROOT)):digest(p) for p in sources},
        'sample_sha256':digest(sample_file),
        'versions':{p:importlib.metadata.version(p) for p in ['camb','cobaya','numpy','scipy','candl-like','sacc']}}
    if model_file:record['surrogate_sha256']=digest(model_file)
    record['identity']=hashlib.sha256(json.dumps(record,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    return record
