"""Supplemental quantile precision with repeated chronological selection slots.

The frozen block bootstrap and midpoint weighted-quantile estimator are imported
unchanged. Repeated expanded indices are legitimate stratified selection slots;
they are neither removed nor merged. This is numerical precision, not a new
posterior qualification or extra astrophysical uncertainty.
"""
import argparse
import json
from pathlib import Path

import numpy as np

import quantile_precision as original

ROOT = original.ROOT
HERE = Path(__file__).resolve().parent


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


def integer(value):
    return isinstance(value, int) and not isinstance(value, bool) and value >= 0


def chronology(selection):
    """Validate slot identity/order without sorting, combining, or removing rows."""
    groups = np.asarray(selection['groups'])
    locations, points = selection['locations'], selection['points']
    assert groups.ndim == 1 and len(groups) == len(locations) == len(points)
    assert all(integer(x) for x in selection['groups'])
    assert set(groups) == {0, 1, 2, 3}
    chain_paths = sorted(selection['chain_sha256'])
    assert len(chain_paths) == 4
    chain_names = [Path(path).name for path in chain_paths]
    assert len(set(chain_names)) == 4
    counts = [int(np.sum(groups == group)) for group in range(4)]
    assert counts[0] >= 500 and len(set(counts)) == 1
    assert np.array_equal(groups, np.repeat(np.arange(4), counts)), 'Chain groups are reordered or interleaved.'
    stored = selection['stored_logposts']
    assert len(stored) == len(groups) and np.isfinite(stored).all()
    summaries = {}
    for group in range(4):
        indices = np.flatnonzero(groups == group)
        local = [locations[i] for i in indices]
        assert all(item['chain'] == chain_names[group] for item in local), 'Location chain does not match group.'
        for item in local:
            assert integer(item['expanded_index']) and integer(item['row'])
        positions = np.array([item['expanded_index'] for item in local], dtype=np.int64)
        rows = np.array([item['row'] for item in local], dtype=np.int64)
        assert np.all(np.diff(positions) >= 0), 'Expanded selection positions decrease.'
        assert np.all(np.diff(rows) >= 0), 'Original RLE row positions decrease.'
        repeated = np.flatnonzero(np.diff(positions) == 0)
        for index in repeated:
            left, right = int(indices[index]), int(indices[index+1])
            assert locations[left]['row'] == locations[right]['row'], 'Repeated expanded slot has a different RLE row.'
            assert points[left] == points[right], 'Repeated expanded slot has a different parameter point.'
            assert stored[left] == stored[right], 'Repeated expanded slot has a different stored proposal density.'
        summaries[str(group)] = {'chain': chain_names[group], 'selection_slots': len(indices),
            'unique_expanded_indices': len(np.unique(positions)),
            'repeated_expanded_index_slots': int(len(indices)-len(np.unique(positions))),
            'unique_original_RLE_rows': len(np.unique(rows)),
            'first_expanded_index': int(positions[0]), 'last_expanded_index': int(positions[-1])}
    return summaries


def verify_inputs(inputs):
    for path, expected in inputs.items():
        assert original.digest(ROOT/path) == expected, 'Qualified input changed: '+path


def actual(folder, summary_path):
    # Must precede cache/selection access: fresh full original parent qualification.
    from measurement_summary import summarize_run
    folder, summary_path = Path(folder).resolve(), Path(summary_path).resolve()
    qualified = summarize_run(folder, summary_path)
    assert qualified['qualified_under_declared_numerical_gates'] is True
    verify_inputs(qualified['input_sha256'])
    parent_summary = json.loads(summary_path.read_text())
    selection_path = (ROOT/parent_summary['selection_path']).resolve()
    assert selection_path.name == 'selection.json' and selection_path.parent.parent == folder
    assert relative(summary_path) in qualified['input_sha256']
    assert relative(selection_path) in qualified['input_sha256']
    assert original.digest(selection_path) == parent_summary['selection_sha256']
    selection = json.loads(selection_path.read_text())
    chronology_report = chronology(selection)
    groups = np.asarray(selection['groups'])
    records = []
    for index in range(len(groups)):
        path = selection_path.parent/f'{index:05d}.json'
        assert relative(path) in qualified['input_sha256']
        assert original.digest(path) == qualified['input_sha256'][relative(path)]
        record = json.loads(path.read_text())
        assert record['index'] == index and record['point'] == selection['points'][index]
        records.append(record)
    names = qualified['weighted_covariance']['parameter_order']
    values = {name: np.array([dict(r['point'], **r['derived'])[name] for r in records]) for name in names}
    # No deduplication, sorting, weight changes, alternate seeds, or estimator changes.
    report = original.block_quantile_precision(values, [r['log_weight'] for r in records], groups)
    baseline_error = 0.
    for name in names:
        expected = qualified['posterior'][name]['quantiles_025_16_50_84_975']
        error = float(np.max(abs(np.array(report['baseline_quantiles'][name])-expected)))
        baseline_error = max(baseline_error, error)
        assert error <= 1e-12, 'Original qualified measurement and diagnostic quantiles differ.'
    verify_inputs(qualified['input_sha256'])
    report.update(target_identity=qualified['target_identity'], settings=qualified['settings'],
        parent_measurement_gates_passed=True, design=original.DESIGN,
        chronology={'rule': 'Nondecreasing expanded indices within unchanged contiguous chain groups; repeated slots remain explicit.',
                    'chains': chronology_report,
                    'total_repeated_expanded_index_slots': sum(x['repeated_expanded_index_slots'] for x in chronology_report.values()),
                    'rows_dropped_or_combined': 0},
        qualified_measurement_baseline_max_abs_error=baseline_error,
        input_sha256=qualified['input_sha256'],
        source_sha256={relative(p): original.digest(p) for p in
                       [Path(__file__), Path(original.__file__), HERE/'measurement_summary.py']},
        original_estimator_source_sha256=original.digest(original.__file__))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--chain-folder', type=Path, required=True)
    parser.add_argument('--correction-summary', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = actual(args.chain_folder, args.correction_summary)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'status': report['status'], 'output': str(args.output),
                      'repeated_selection_slots': report['chronology']['total_repeated_expanded_index_slots']}))


if __name__ == '__main__':
    main()
