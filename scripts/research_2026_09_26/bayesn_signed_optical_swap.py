"""Run the frozen forward gate with a synthetic sealed NIR payload, without opening the real sealed NIR file."""
from __future__ import annotations
import hashlib
import json
import shutil
from pathlib import Path
import numpy as np
import bayesn_signed_optical_gate as gate

def main():
    root=gate.OUT
    scratch=root/'synthetic-nir-swap'
    scratch.mkdir(exist_ok=True)
    for name in ('optical-payload.json','release-filter-config.yaml','execution-protocol.json','RAISIN_DES_J.dat'):
        shutil.copy2(root/name,scratch/name)
    synthetic=dict(scope='Synthetic sentinel only; no observed NIR data read.',rows=[dict(MJD=1.,band='J',flux=-1e9,error=1e-3)])
    (scratch/'nir-sealed.json').write_text(json.dumps(synthetic,indent=2)+'\n')
    manifest=json.loads((root/'adapter-manifest.json').read_text())
    manifest['outputs']['nir-sealed.json']=gate.sha(scratch/'nir-sealed.json')
    (scratch/'adapter-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    gate.OUT=scratch
    gate.main()
    equal={}
    for bins in (300,600,1200,2400):
        with np.load(root/f'forward-{bins}.npz') as a,np.load(scratch/f'forward-{bins}.npz') as b:
            equal[str(bins)]=bool(np.array_equal(a['prediction'],b['prediction']) and
                                   np.array_equal(a['errors'],b['errors']) and np.array_equal(a['mask'],b['mask']))
    assert all(equal.values())
    result=dict(scope='Structural numerical invariance to synthetic NIR payload swap; no observed NIR opened by this runner.',
                optical_input_sha256=gate.sha(scratch/'optical-payload.json'),synthetic_NIR_sha256=gate.sha(scratch/'nir-sealed.json'),
                predictions_and_error_mask_exact_equal=equal,forward_gate_pass=bool(json.loads((scratch/'forward-algebra-result.json').read_text())['gates_pass']))
    (root/'synthetic-nir-swap-result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))
if __name__=='__main__':main()
