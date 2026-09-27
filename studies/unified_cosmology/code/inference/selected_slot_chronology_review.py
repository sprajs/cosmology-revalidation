"""Reconstruct stratified correction slots without changing their multiplicity."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
RESULTS = ROOT/'studies/unified_cosmology/results/inference'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rel(path):
    return str(Path(path).resolve().relative_to(ROOT))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--correction-summary', type=Path,
                        default=RESULTS/'modern-lcdm-independence-exact-correction.json')
    parser.add_argument('--output', type=Path, default=RESULTS/'selected-slot-chronology-lcdm-review.json')
    args = parser.parse_args()
    summary = json.loads(args.correction_summary.read_text())
    assert summary['status'] == 'passed_importance_weight_gates'
    selection_path = ROOT/summary['selection_path']
    assert sha(selection_path) == summary['selection_sha256']
    selection = json.loads(selection_path.read_text())
    diagnostic_path = ROOT/selection['diagnostics_path']
    assert sha(diagnostic_path) == selection['diagnostics_sha256']
    diagnostic = json.loads(diagnostic_path.read_text())
    assert diagnostic['status'] == 'passed'
    source = Path(__file__).with_name('exact_correction.py')
    assert sha(source) == summary['code_sha256']
    bound = {rel(p):sha(p) for p in [args.correction_summary, selection_path, diagnostic_path, source]}
    bound.update(selection['chain_sha256'])
    assert all(sha(ROOT/p) == expected for p,expected in bound.items())
    names = list(selection['points'][0])
    length = diagnostic['equal_chain_length_for_diagnostics']
    groups = np.asarray(selection['groups'])
    assert set(groups) == {0,1,2,3}
    random = np.random.default_rng(272639)
    records = []
    for group in range(4):
        slots = np.flatnonzero(groups == group)
        assert np.array_equal(slots,np.arange(slots[0],slots[-1]+1))
        n = len(slots)
        locations = [selection['locations'][int(i)] for i in slots]
        chain_names = {v['chain'] for v in locations}
        assert len(chain_names) == 1
        paths = [ROOT/p for p in selection['chain_sha256'] if Path(p).name in chain_names]
        assert len(paths) == 1
        path = paths[0]
        header = path.open().readline().lstrip('#').split()
        raw = np.loadtxt(path,ndmin=2)
        assert np.isfinite(raw[:,0]).all() and np.array_equal(raw[:,0],raw[:,0].astype(int))
        assert (raw[:,0] > 0).all()
        origin = np.repeat(np.arange(len(raw)),raw[:,0].astype(int))
        assert length >= n and len(origin)-int(.3*len(origin)) >= length
        discarded = len(origin)-length
        indices = np.floor((np.arange(n)+random.random(n))*length/n).astype(int)
        chronology = np.array([p['expanded_index'] for p in locations])
        assert np.array_equal(chronology,discarded+indices)
        assert np.all(np.diff(chronology)>=0)
        for slot, expanded in zip(slots,chronology):
            row = int(origin[expanded])
            assert selection['locations'][int(slot)]['row'] == row
            point = {name:float(raw[row,header.index(name)]) for name in names}
            assert point == selection['points'][int(slot)]
            assert -float(raw[row,header.index('minuslogpost')]) == selection['stored_logposts'][int(slot)]
        # A given integer bin may be split by a boundary between real strata.
        # Sum exact geometric inclusion probabilities over all strata.
        width = length/n
        multiplicity = np.zeros(length)
        for j in range(n):
            lo,hi = j*width,(j+1)*width
            for k in range(int(np.floor(lo)),min(length,int(np.ceil(hi)))):
                multiplicity[k] += max(0,min(hi,k+1)-max(lo,k))/width
        error = float(np.max(abs(multiplicity-n/length)))
        assert error < 1e-12
        ties = np.flatnonzero(np.diff(chronology)==0)
        split = int(.8*n)
        records.append({'chain':group,'selected_slots':n,'unique_expanded_indices':len(set(chronology)),
                        'discarded_expanded_states':discarded,'stratum_width':width,
                        'reconstructed_indices_points_RLE_rows_and_stored_density':True,
                        'expected_index_multiplicity':n/length,'maximum_multiplicity_error':error,
                        'adjacent_tied_slot_pairs':[[int(i),int(i+1)] for i in ties],
                        'crosses_80percent_split':bool(split-1 in ties),
                        'crosses_10block_boundary':bool(any((i+1)%(n//10)==0 for i in ties)),
                        'crosses_20block_boundary':bool(any((i+1)%(n//20)==0 for i in ties))})
    assert all(sha(ROOT/p) == expected for p,expected in bound.items())
    output = {'status':'passed_independent_stratified_slot_chronology_review',
              'seed':272639,'expanded_length_per_chain':length,'selected_slots':len(groups),
              'chains':records,'input_sha256':bound,'source_sha256':{rel(__file__):sha(__file__)},
              'numpy_version':np.__version__,'model_background_spectrum_calls':0,
              'observed_parameters_reestimated':False,
              'meaning':[
                  'Independent uniform positions in disjoint real strata can floor to the same integer state in adjacent strata.',
                  'All selected slots and their raw importance weights must remain; deduplicating changes the estimator.',
                  'Nondecreasing expanded indices preserve chronology. They do not turn repeated or correlated states into independent samples.',
                  'Conditional on the frozen expanded chain, expected index multiplicity is uniform. Self-normalized importance ratios retain usual finite-sample ratio error.',
                  'A generic chronological train/test split must disclose or prevent the same selected state straddling its boundary. No such boundary tie occurs in this cohort.',
                  'Chronological block bootstrap is an approximate dependence diagnostic, not a proof of independent data or global posterior coverage.'
              ]}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':output['status'],'selected_slots':len(groups),'chains':records}))


if __name__ == '__main__':
    main()
