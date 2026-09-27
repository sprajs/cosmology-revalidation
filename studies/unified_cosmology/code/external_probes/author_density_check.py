"""Independent normalization and calibration-diagnostic summary of author rows.

No fitted constant. This is a post-inspection diagnostic, not a claim that the
public source was the exact runtime or that every author row was re-evaluated.
"""
import json
import numpy as np
from scipy.linalg import cho_factor,cho_solve
from acquire import ROOT,RESULTS,sha


def main():
    source=RESULTS/'author-density-reconstruction.json';r=json.loads(source.read_text())
    assert len(r['rows'])==5
    file=ROOT/'.work/unified-cosmology/survey-selection/normalized/dovekie-total.npz'
    with np.load(file) as f:C=f['covariance']
    ones=np.ones(len(C));u=cho_solve(cho_factor(C,lower=True),ones);A=float(ones@u)
    normalization=float(.5*np.log(A/(2*np.pi)))
    alternatives={}
    for field in ['density_difference','alternative_Aplanck_in_Planck_only_difference','alternative_Aplanck_in_Planck_and_ACT_difference']:
        residual=np.array([row[field] for row in r['rows']])-normalization
        alternatives[field]={'absolute_loglike_residuals':residual.tolist(),'max_abs_loglike_residual':float(np.max(abs(residual))),
            'max_minus_min_loglike_residual':float(np.ptp(residual))}
    out={'status':'numerical_reconstruction_checked_at_five_rows','selection':r['selection'],
        'rows':r['rows_selected'],'code_sha256':sha(__file__),'source_result_sha256':sha(source),'normalized_SN_data_sha256':sha(file),
        'SN_absolute_magnitude_normalization':{'A':A,'half_log_A_over_2pi':normalization,'fitted_additive_constant':False},
        'calibration_alternatives':alternatives,'interpretation':[
            'Literal public header/source misses the recorded likelihood shape.',
            'Shared A_planck calibration for both Planck and ACT is a source-informed post-inspection diagnostic, evaluated on the same exact spectra.',
            'Absolute agreement after the independently derived SN normalization is reported numerically; no tolerance was retroactively called a predeclared scientific gate.',
            'Agreement on five retained rows supports a useful numerical reconstruction, not identity of unrecorded author source or dependency versions.',
            'The original declared modern target and training data remain unchanged.'],
        'dependency_records':{name:sha(RESULTS/name) for name in ['author-acquisition.json','author-configuration.json','modern-acquisition.json']},
        'exact_author_runtime_recovered':False,'can_claim_independent_cosmological_measurement':False}
    (RESULTS/'author-density-check.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

if __name__=='__main__':main()
