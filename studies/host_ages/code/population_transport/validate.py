#!/usr/bin/env python3
"""Independent algebra checks and identity verification for transport analyses."""
from pathlib import Path
import hashlib,json,datetime
import numpy as np
from scipy.optimize import brentq
from scipy.linalg import cho_factor,cho_solve
from run import weighted_mean_range, ROOT, OUT


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    checked=[]
    for name in ['run.json','titan-run.json']:
        rec=json.loads((OUT/name).read_text())
        for section in ['input_sha256','additional_input_sha256','code_sha256','output_sha256']:
            for file,expected in rec.get(section,{}).items():
                assert sha(ROOT/file)==expected,file
                checked.append(file)
    rng=np.random.default_rng(81829)
    gaps=[]
    for n in [5,29,196]:
        x=rng.uniform(.04,12,n)
        for ratio in [1,2,5]:
            observed=weighted_mean_range(x,ratio)
            # Independent fractional-optimization root: the optimal weights are
            # maximal above the unknown optimum for a maximum, below for minimum.
            roots=[]
            for maximum in [False,True]:
                def f(m):
                    w=np.where((x>m) if maximum else (x<m),ratio,1.)
                    return np.sum(w*(x-m))
                roots.append(brentq(f,min(x),max(x)))
            gaps.append(float(np.max(np.abs(np.array(roots)-observed))))
    assert max(gaps)<1e-10
    # Known-age slope information agrees with full augmented GLS inversion.
    n=53;r=rng.normal(size=(n,n));c=r@r.T/n+np.eye(n)*.02
    x=np.column_stack([np.ones(n),rng.normal(size=(n,3))]);a=rng.normal(size=n)
    cf=cho_factor(c);ix=cho_solve(cf,x)
    ar=a-x@np.linalg.solve(x.T@ix,ix.T@a)
    inverse_information=1/(ar@cho_solve(cf,ar))
    d=np.column_stack([x,a]);direct=np.linalg.inv(d.T@cho_solve(cf,d))[-1,-1]
    assert abs(inverse_information-direct)<1e-12
    summary=json.loads((OUT/'summary.json').read_text())
    identities=[]
    for policy,cases in summary['mapping_scenarios'].items():
        for case in cases:
            if case['mapping'] in ['affine','mean_preserving_scatter']:
                identities.append(abs(case['departure_from_simple_compensation_mag']))
    assert max(identities)<1e-12
    titan=json.loads((OUT/'titan-summary.json').read_text())
    assert titan['source_column_max_abs_difference']<1e-13
    assert titan['rows']==8610 and not titan['final_sample_reproduced']
    out=dict(checked_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        verified_file_identities=len(checked),unique_verified_files=len(set(checked)),
        selection_extrema_independent_fractional_optimization_max_error_gyr=max(gaps),
        GLS_partial_information_variance_identity_abs_error=float(abs(inverse_information-direct)),
        affine_and_mean_preserving_scatter_max_compensation_error_mag=max(identities),
        author_snapshot_original_property_max_difference=titan['source_column_max_abs_difference'],
        passed=True,
        scope='Numerical/algebra and provenance checks; not observational identification or coverage validation of the physical age model.',
        validator_sha256=sha(Path(__file__)))
    (OUT/'validation.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
