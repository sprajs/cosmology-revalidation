#!/usr/bin/env python3
"""Compile unmodified SNANA extinction routines, compare historical/current curves.

Only diagnostics are stubbed; source function bodies are copied byte-for-byte.
No SNANA original input/source is changed. All derivatives use E(B-V) mag.
"""
from __future__ import annotations
import ctypes
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
CURRENT = ROOT / 'sources/repos/RickKessler__SNANA/src'
HISTORICAL = ROOT / 'phase2/official/build/SNANA-2fe0f56/src'
OUT = ROOT / 'runs/salt_dust_audit/snana_extinction'


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def build_library(out=OUT):
    out.mkdir(parents=True, exist_ok=True)
    # Extract contiguous extinction-routine section, leaving the original math intact.
    current = (CURRENT / 'MWgaldust.c').read_text()
    historical = (HISTORICAL / 'MWgaldust.c').read_text()
    section = current[current.index('double GALextinct(double'):current.index('// ========== FUNCTION TO RETURN EBV(SFD)')]
    old = historical[historical.index('double GALextinct(double'):historical.index('// ========== FUNCTION TO RETURN EBV(SFD)')]
    old = old.replace('GALextinct(', 'GALextinct_historical(', 1)
    # Historical option constant 99 must remain 99, not today's -99 symbol.
    old = old.replace('OPT_MWCOLORLAW_FITZ99', '99')
    prelude = '''#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "MWgaldust.h"
#define SEV_FATAL 4
char c1err[4096], c2err[4096];
void errmsg(int severity, int ignored, char *fnam, char *one, char *two) {
  fprintf(stderr, "%s: %s; %s\\n", fnam, one, two); abort();
}
void concat_callfun_plus_fnam(char *one, char *two, char *out) {
  snprintf(out, 60, "%s", two);
}
'''
    src = out / 'snana_extinction_extracted.c'
    src.write_text(prelude + section + '\n' + old)
    so = out / 'libsnana_extinction.so'
    cmd = ['cc', '-O2', '-fPIC', '-shared', '-I', str(CURRENT), str(src), '-lm', '-o', str(so)]
    run = subprocess.run(cmd, text=True, capture_output=True)
    (out / 'build.log').write_text('COMMAND ' + ' '.join(cmd) + '\n' + run.stdout + run.stderr)
    run.check_returncode()
    return so


class SnanaExtinction:
    def __init__(self, library=None):
        self.library_path = Path(library) if library else build_library()
        self.lib = ctypes.CDLL(str(self.library_path))
        self.lib.GALextinct.argtypes = [ctypes.c_double] * 3 + [ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.c_char_p]
        self.lib.GALextinct.restype = ctypes.c_double
        self.lib.GALextinct_historical.argtypes = [ctypes.c_double] * 3 + [ctypes.c_int]
        self.lib.GALextinct_historical.restype = ctypes.c_double
        self.pars = (ctypes.c_double * 4)(-99, -99, -99, -99)

    def __call__(self, wave, rv=3.1, ebv=1.0, option=99, historical=False):
        arr = np.asarray(wave, dtype=float)
        if not (np.all(arr >= 1000) and np.all(arr <= (25000 if historical or option == -99 else 35000))):
            raise ValueError('Requested wavelength outside this audit wrapper validated range')
        if not 0.01 <= rv <= 8.0:
            raise ValueError('RV outside released-grid audit range; RV<2 or RV>6 is extrapolation')
        call = self.lib.GALextinct_historical if historical else self.lib.GALextinct
        if historical:
            values = [call(rv, rv * ebv, float(w), option) for w in arr.flat]
        else:
            values = [call(rv, rv * ebv, float(w), option, self.pars, b'audit') for w in arr.flat]
        return np.asarray(values).reshape(arr.shape)


def main():
    ext = SnanaExtinction()
    waves = np.arange(2000., 10001., 10.)
    rvs = [1.0, 1.5, 2.0, 2.5, 3.1, 4.0, 5.0, 6.0]
    records = []
    historic_max = 0.0
    independent_max = 0.0
    try:
        import extinction
    except ImportError as exc:
        raise RuntimeError('Run with phase2/env-official/bin/python (extinction is pinned there)') from exc
    for rv in rvs:
        exact = ext(waves, rv)
        approx = ext(waves, rv, option=-99)
        historical = ext(waves, rv, historical=True)
        reference = extinction.fitzpatrick99(waves, rv, rv)
        historic_max = max(historic_max, float(np.max(np.abs(approx-historical))))
        independent_max = max(independent_max, float(np.max(np.abs(exact-reference))))
        for w, e, a, h, ref in zip(waves, exact, approx, historical, reference):
            records.append((w, rv, e, a, h, ref, e-a))
    table = np.asarray(records)
    np.savetxt(OUT / 'wavelength_response.csv', table, delimiter=',', comments='', header='wavelength_A,RV,dA_dEBV_exact,dA_dEBV_approx,dA_dEBV_historical,dA_dEBV_extinction_py,exact_minus_approx_mag_per_EBV')
    # Natural derivative matrix for monochromatic magnitudes; not fitted distances.
    pivot = np.array([3500., 4000., 4400., 5500., 6500., 8000.])
    rows = []
    for rv in rvs:
        delta = ext(pivot, rv) - ext(pivot, rv, option=-99)
        exact = ext(pivot, rv)
        # dA/dRv at fixed EBV=0.1, not fixed AV.
        if 1 < rv < 6:
            drv = (ext(pivot, rv+1e-4, .1)-ext(pivot, rv-1e-4, .1))/2e-4
        else:
            drv = np.full_like(pivot, np.nan)
        for w,e,d,dr in zip(pivot,exact,delta,drv):
            rows.append((rv,w,e,dr,.1*d))
    np.savetxt(OUT / 'monochromatic_response_matrix.csv', rows, delimiter=',', comments='', header='RV,wavelength_A,dA_dEBV_mag_per_mag,dA_dRV_at_EBV0p1,exact_minus_approx_at_EBV0p1_mag')
    stats=[]
    for rv in rvs:
        m=(table[:,1]==rv)&(table[:,0]>=3500)&(table[:,0]<=8000)
        d=table[m,6]
        stats.append({'RV':rv,'max_abs_delta_A_at_EBV0p1_mag':float(.1*np.max(np.abs(d))), 'lambda_max_A':float(table[m,0][np.argmax(np.abs(d))]), 'rms_delta_A_at_EBV0p1_mag':float(.1*np.sqrt(np.mean(d*d)))})
    result = {'meaning':'Monochromatic extinction implementation response, not a Hubble-distance correction or empirical dust uncertainty.', 'current_commit':'886408a4e171896db5eaa97e735a655f50cec2db', 'historical_build':'SNANA-2fe0f56', 'option_semantics':{'historical_99':'polynomial approximation', 'current_99':'F99 natural cubic spline', 'current_minus99':'historical polynomial approximation'},'checks':{'current_approx_vs_historical_max_abs_mag_per_EBV':historic_max,'current_exact_vs_extinction_py_max_abs_mag_per_EBV':independent_max}, 'grid':{'wavelength_A':[2000,10000,10], 'RV':rvs, 'low_RV_caution':'RV<2 is outside current header recommended RV minimum; source explicitly does not enforce it.'}, 'optical_3500_8000':stats,'versions':{'python':sys.version,'numpy':np.__version__,'extinction':extinction.__version__}}
    assert historic_max < 1e-12, historic_max
    assert independent_max < 1e-7, independent_max
    (OUT/'results.json').write_text(json.dumps(result, indent=2)+'\n')
    files=[Path(__file__).resolve(), CURRENT/'MWgaldust.c', CURRENT/'MWgaldust.h', HISTORICAL/'MWgaldust.c', HISTORICAL/'MWgaldust.h', ext.library_path, OUT/'snana_extinction_extracted.c', OUT/'wavelength_response.csv', OUT/'monochromatic_response_matrix.csv', OUT/'results.json']
    (OUT/'manifest.json').write_text(json.dumps({'files':[{'path':str(p.relative_to(ROOT)), 'sha256':sha256(p),'bytes':p.stat().st_size} for p in files]},indent=2)+'\n')
    print(json.dumps(result, indent=2))

if __name__ == '__main__':
    main()
