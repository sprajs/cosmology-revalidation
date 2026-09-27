"""Count distinct released probability fields on the exact distance cohort.

No probability is re-fitted, and no covariance or Gaussian reference degrees
of freedom are changed. The count comparison does not reconstruct BBC.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
from lib.records import fitres

REGISTRY = ROOT/'studies/unified_cosmology/results/survey_selection/dovekie-interface.json'
WORK = ROOT/'.work/unified-cosmology/survey-selection'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def relative(path):
    return str(Path(path).relative_to(ROOT))


def run():
    paths = [WORK/'release/4_DISTANCES_COVMAT'/name for name in
             ['DES-Dovekie_HD.csv', 'DES-Dovekie_Metadata.csv']]
    ledger_path = WORK/'normalized/dovekie-ledger.csv'
    vector_path = WORK/'normalized/dovekie-total.npz'
    registry = json.loads(REGISTRY.read_text())
    expected = registry['inputs'] | registry['outputs']
    for p in paths+[ledger_path, vector_path]:
        assert digest(p) == expected[relative(p)], 'Changed released cohort input: '+str(p)
    hd, all_metadata = [fitres(p) for p in paths]
    assert hd.index.is_unique and all_metadata.index.is_unique
    assert len(hd) == 1820 == registry['rows']
    assert set(hd.index) == set(all_metadata.index), 'Different metadata physical-event cohort.'
    meta = all_metadata.loc[hd.index]
    ledger = pd.read_csv(ledger_path, dtype={'CID': str})
    assert np.array_equal(ledger.row, np.arange(1820))
    assert np.array_equal(ledger.CID, hd.index.astype(str))
    for key in ['IDSURVEY', 'zHD', 'zHEL', 'MU']:
        assert np.array_equal(hd[key], meta[key]), 'HD/metadata mismatch: '+key
    with np.load(vector_path, allow_pickle=False) as data:
        assert np.array_equal(data['CID'], hd.index.astype(str))
        assert data['covariance'].shape == data['precision'].shape == (1820, 1820)
        for key in ['zHD', 'zHEL', 'MU']:
            assert np.array_equal(data[key], hd[key]), 'Normalized vector mismatch: '+key
    beams = hd.PROBIA_BEAMS.to_numpy(float)
    cc = meta.PROBCC_BEAMS.to_numpy(float)
    assert np.isfinite(beams).all() and np.isfinite(cc).all()
    assert ((beams >= 0)&(beams <= 1)).all() and ((cc >= 0)&(cc <= 1)).all()
    closure = float(np.max(abs(beams-(1-cc))))
    assert closure <= 5.1e-6, 'Released probability complements differ beyond printed rounding.'
    des = meta.IDSURVEY.to_numpy() == 10
    assert des.sum() == 1623 and (~des).sum() == 197
    assert np.all(beams[~des] == 1) and np.all(cc[~des] == 0)
    assert (meta.loc[~des, 'TYPE'] == 1).all(), 'Non-DES sample is not the expected spectroscopic type.'
    classifier = meta.loc[des, 'PROB_SNNV19'].to_numpy(float)
    assert np.isfinite(classifier).all() and ((classifier >= 0)&(classifier <= 1)).all()
    classifier_sum = float(classifier.sum()+(~des).sum())
    inputs = {relative(p): digest(p) for p in paths+[ledger_path, vector_path, REGISTRY]}
    result = {'status': 'passed_released_probability_count_accounting',
        'cohort': {'rows': 1820, 'DES': 1623, 'spectroscopic_non_DES': 197,
                   'CID_unique': True, 'metadata_identifier_survey_redshift_distance_join': 'exact',
                   'normalized_vector_order': 'exact', 'covariance_shape': [1820, 1820]},
        'released_BEAMS_posterior_diagnostic': {
            'HD_column': 'PROBIA_BEAMS', 'HD_header_comment_alias': 'PROB1A_BEAMS',
            'sum': float(beams.sum()), 'DES_sum': float(beams[des].sum()),
            'metadata_complement_column': '1-PROBCC_BEAMS',
            'metadata_complement_sum': float((1-cc).sum()),
            'maximum_per_row_complement_difference': closure},
        'classifier_input_probability_plus_spectroscopic_sample': {
            'DES_column': 'PROB_SNNV19', 'DES_sum': float(classifier.sum()),
            'spectroscopic_contribution': 197, 'sum': classifier_sum,
            'rounds_to_paper_1684': int(round(classifier_sum)) == 1684},
        'paper': {'url': 'https://arxiv.org/html/2511.07517v3',
                 'section': '10.5 and Appendix C', 'printed_effective_count': 1684},
        'interpretation': 'The classifier-input sum rounds to the printed count; the released BBC posterior diagnostic does not. This identifies a plausible field distinction, not the unreleased author counting implementation or an author error. Neither probability sum replaces the 1820-row Gaussian rank in a replication model.',
        'covariance_modified': False, 'Gaussian_reference_df_modified': False,
        'physical_model_or_background_calls': 0, 'input_sha256': inputs,
        'source_sha256': {relative(p): digest(p) for p in [Path(__file__).resolve(), ROOT/'lib/records.py']}}
    assert all(digest(ROOT/k) == value for k, value in inputs.items())
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists(), 'Preserve previous results.'
    value = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'status': value['status'], 'output': str(args.output)}))
