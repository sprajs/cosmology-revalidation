#!/usr/bin/env python3
"""Regenerate a historical real design, compare solver residuals.

--prepare needs the physical-age environment and recovered SDSS data.
Without it, only NumPy/SciPy are required, enabling version comparisons.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import sys,json,hashlib
from pathlib import Path
import numpy as np
import scipy
from scipy.optimize import nnls,lsq_linear
ROOT=Path(__file__).resolve().parents[4];WORK=ROOT/'.work/physical-ages'

def main():
    path=WORK/'nnls-regression-fixture.npz'
    if '--prepare' in sys.argv:
        from fit import inputs
        from model import Grid, Library
        item=next(x for x in inputs('sdss')[0] if x['id']=='1003262072500283392')
        # Preserve the exact physical representation that exposed the bug.
        # Production photometry now uses better Gauss quadrature; changing this
        # regression fixture would conflate a solver test with that correction.
        lib=Library(); pars=Grid(lib).parameters[14]
        metal=list(lib.data['metallicity']).index(pars['metallicity'])
        _,spectra,_=lib.at_redshift(item['z'])
        a=lib.photometry(spectra[metal],item['z'],item['ebv'],pars['tau'],pars['power'],integration='native')
        error=np.array(item['error']);a=a/error[:,None];a/=np.linalg.norm(a,axis=0)
        y=np.array(item['flux'])/error
        np.savez(path,A=a,y=y)
    fixture=np.load(path);a=fixture['A'];y=fixture['y'];acopy=a.copy();ycopy=y.copy()
    x,reported=nnls(acopy,ycopy,maxiter=2000)
    direct=float(np.linalg.norm(a@x-y)**2)
    fit=lsq_linear(a,y,bounds=(0,np.inf),method='bvls',tol=1e-11,max_iter=1000)
    result=dict(scipy=scipy.__version__,numpy=np.__version__,shape=list(a.shape),object='1003262072500283392',nuisance_cell=14,
        fixture_definition='Historical observer-grid trapezoid quadrature, deliberately retained only for this solver regression; final physical fits use Gauss quadrature.',
        nnls_reported_squared_residual=float(reported**2),nnls_recomputed_squared_residual=direct,
        nnls_residual_identity_passed=bool(abs(reported**2-direct)<1e-8),nnls_input_mutation_max=float(max(abs(a-acopy).max(),abs(y-ycopy).max())),
        bounded_variable_squared_residual=float(np.linalg.norm(a@fit.x-y)**2),bounded_variable_KKT=float(fit.optimality),bounded_variable_success=bool(fit.success),
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        boundary='A documented numerical failure of the tested version/design; it is not a new astrophysical result. New age fitting uses column-scaled BVLS and recomputes the objective. Older BAO/RAISIN uses are separately audited on their own designs.')
    output=ROOT/f'studies/host_ages/results/physical_ages/solver-audit-scipy-{scipy.__version__}-numpy-{np.__version__}.json'
    output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
