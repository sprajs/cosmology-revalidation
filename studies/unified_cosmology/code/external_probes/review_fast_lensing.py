"""Independent synthetic native-vs-compressed ACT/Planck response check."""
import sys, json, hashlib
from pathlib import Path
import numpy as np
from acquire import ROOT, RESULTS
import act_dr6_lenslike.act_dr6_lenslike as native_module
from fast_lensing import BinnedResponse
from act_dr6_lenslike.act_dr6_lenslike import get_corrected_clkk
rng = np.random.default_rng(927615)
errors = []
for _ in range(20):
    n = 9
    d = {'include_planck': True, 'likelihood_corrections': True,
         'fiducial_cl_kk': rng.normal(size=n)}
    for s in ['tt', 'ee', 'bb', 'te']:
        d['fiducial_cl_' + s] = rng.normal(size=n)
    for suffix, key, m in [('', 'binmat_act', 4), ('_planck', 'binmat_planck', 3)]:
        d[key] = rng.normal(size=(m,n))
        d['fAL' + suffix] = rng.uniform(.5, 2, n)
        d['fAL' + suffix][:2] = 0
        d['dAL_dC' + suffix] = rng.normal(size=(4,n,n))
        for s in ['kk', 'tt', 'ee', 'bb', 'te']:
            d['dN1_' + s + suffix] = rng.normal(size=(n,n))
    kk = rng.normal(size=n)
    cmb = {s:rng.normal(size=n) for s in ['tt', 'ee', 'bb', 'te']}
    native = np.concatenate([
        d[key] @ get_corrected_clkk(d, kk, cmb['tt'], cmb['te'], cmb['ee'], cmb['bb'], suffix)
        for suffix,key in [('', 'binmat_act'), ('_planck', 'binmat_planck')]])
    compressed = BinnedResponse(d).predict(kk, cmb)
    errors.append(float(abs(compressed - native).max()))
assert max(errors) < 1e-12
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
source = Path(__file__).resolve().with_name('fast_lensing.py')
native_source = Path(native_module.__file__).resolve()
result = {'status':'passed', 'seed':927615, 'synthetic_matrix_cases':20,
          'nonzero_fiducial_low_ell_and_zero_fAL_low_ell':True,
          'maximum_absolute_band_difference':max(errors), 'case_errors':errors,
          'source_sha256':{str(p.resolve().relative_to(ROOT)):sha(p) for p in [Path(__file__),source,native_source]},
          'authorship':'Separately authored independent review; only import/output paths adapted for repository reproduction.',
          'scope':'Independent algebra/orientation/lowell stress check; curves are not physical cosmologies. Does not test posterior convergence or omitted measurement covariance.'}
(RESULTS/'fast-lensing-independent.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
