"""Synthetic-only chronology/provenance and unchanged-bootstrap controls."""
import argparse
from contextlib import ExitStack
import copy
import json
from pathlib import Path
import tempfile
from unittest.mock import patch

import numpy as np

import quantile_precision_ordered as consumer

original = consumer.original
ROOT = consumer.ROOT
OUTPUT = ROOT/'studies/unified_cosmology/results/inference/quantile-precision-ordered-validation.json'


def write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, allow_nan=False)+'\n')


def quantiles(x, raw_logweight):
    # Independent normalized midpoint CDF, including every repeated slot.
    weights = np.exp(raw_logweight-np.max(raw_logweight))
    weights /= np.sum(weights)
    order = np.argsort(x)
    xx, ww = x[order], weights[order]
    cdf = np.concatenate(([0.], np.cumsum(ww)))
    centers = .5*(cdf[:-1]+cdf[1:])
    return np.interp([.025, .16, .5, .84, .975], centers, xx)


def fixture(work):
    folder = work/'chains'; exact = folder/'selected'
    paths = []
    rng = np.random.default_rng(273709)
    x = rng.normal(size=2000); y = .2+.7*x+rng.normal(size=2000)*.3
    positions = np.tile(np.arange(500)*2, 4)
    rows = np.tile(np.arange(500), 4)
    # Two adjacent selections share an expanded index, with unchanged source identity.
    for right in [500+67, 1000+104]:
        positions[right] = positions[right-1]; rows[right] = rows[right-1]
        x[right] = x[right-1]; y[right] = y[right-1]
    lw = .5*np.sin(x)+.2*y
    groups = np.repeat(np.arange(4), 500)
    chain_sha = {}
    for i in range(4):
        p = folder/f'chain.{i+1}.txt'; p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(f'Explicit synthetic parent chain {i}\n')
        paths.append(p); chain_sha[consumer.relative(p)] = original.digest(p)
    points = [{'x': float(v)} for v in x]
    records = []
    for i in range(2000):
        row = {'index': i, 'point': points[i], 'derived': {'y': float(y[i])}, 'log_weight': float(lw[i])}
        p = exact/f'{i:05d}.json'; write(p, row); paths.append(p); records.append(row)
    selection = {'groups': groups.tolist(), 'points': points,
        'locations': [{'chain': f'chain.{i//500+1}.txt', 'row': int(rows[i]), 'expanded_index': int(positions[i])} for i in range(2000)],
        'chain_sha256': chain_sha, 'stored_logposts': (-x*x).tolist()}
    selection_path = exact/'selection.json'; write(selection_path, selection); paths.append(selection_path)
    summary_path = work/'correction.json'
    write(summary_path, {'selection_path': consumer.relative(selection_path), 'selection_sha256': original.digest(selection_path)})
    paths.append(summary_path)
    qualified = {'qualified_under_declared_numerical_gates': True, 'target_identity': 'explicit-synthetic-fixture',
        'settings': {'synthetic_only': True}, 'weighted_covariance': {'parameter_order': ['x', 'y']},
        'posterior': {name: {'quantiles_025_16_50_84_975': quantiles(v, lw).tolist()} for name, v in [('x', x), ('y', y)]},
        'input_sha256': {consumer.relative(p): original.digest(p) for p in paths}}
    return folder, summary_path, selection, qualified, {'x': x, 'y': y}, lw, groups


def refused(call):
    try:
        call()
    except (AssertionError, ValueError):
        return True
    raise AssertionError('Required refusal did not occur')


def independent_bootstrap(values, lw, groups):
    answer = {}
    for count, seed in [(10, 928271), (20, 928272)]:
        rng = np.random.default_rng(seed)
        draws = {name: [] for name in values}
        for _ in range(500):
            selected = []
            for group in range(4):
                positions = np.flatnonzero(groups == group)
                width = len(positions)//count
                for block in rng.choice(count, size=count, replace=True):
                    selected.extend(positions[block*width:(block+1)*width])
            selected = np.asarray(selected)
            assert len(selected) == len(groups)
            for name, values_one in values.items():
                draws[name].append(quantiles(values_one[selected], lw[selected]))
        answer[str(count)] = {name: np.asarray(rows) for name, rows in draws.items()}
    return answer


def run():
    import measurement_summary
    work = ROOT/'.work/unified-cosmology/inference/quantile-ordered-validation'
    work.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=work) as temporary:
        folder, summary, selection, parent, values, lw, groups = fixture(Path(temporary))
        calls = []
        def qualifier(*args):
            calls.append(args)
            return copy.deepcopy(parent)
        with patch.object(measurement_summary, 'summarize_run', qualifier):
            report = consumer.actual(folder, summary)
            assert len(calls) == 1 and calls[0] == (folder.resolve(), summary.resolve())
            # Preserve explicit evidence that the original strictly increasing guard refuses this fixture.
            assert refused(lambda: original.actual(folder, summary))
        assert report['points'] == 2000
        assert report['chronology']['rows_dropped_or_combined'] == 0
        assert report['chronology']['total_repeated_expanded_index_slots'] == 2
        assert [report['chronology']['chains'][str(i)]['unique_expanded_indices'] for i in range(4)] == [500, 499, 499, 500]
        direct, saved_draws = original.block_quantile_precision(values, lw, groups, return_draws=True)
        for key, value in direct.items():
            assert report[key] == value
        independent = independent_bootstrap(values, lw, groups)
        maxima = {'weighted_quantile': 0., 'bootstrap_draw': 0., 'bootstrap_SD': 0., 'bootstrap_interval': 0.}
        for name, v in values.items():
            maxima['weighted_quantile'] = max(maxima['weighted_quantile'], float(np.max(abs(quantiles(v, lw)-report['baseline_quantiles'][name]))))
            for count in ['10', '20']:
                draw = independent[count][name]
                maxima['bootstrap_draw'] = max(maxima['bootstrap_draw'], float(np.max(abs(draw-saved_draws[count][name]))))
                result = report['partitions'][count]['marginals'][name]
                maxima['bootstrap_SD'] = max(maxima['bootstrap_SD'], float(np.max(abs(draw.std(axis=0, ddof=1)-result['bootstrap_SD']))))
                maxima['bootstrap_interval'] = max(maxima['bootstrap_interval'], float(np.max(abs(np.quantile(draw, [.025, .975], axis=0).T-result['bootstrap_quantiles_025_975']))))
        assert max(maxima.values()) < 1e-12, maxima
        controls = {}
        changes = {
            'decreasing_expanded_position': lambda s: s['locations'][100].update(expanded_index=0),
            'decreasing_RLE_row': lambda s: s['locations'][100].update(row=0),
            'wrong_group_chain': lambda s: s['locations'][600].update(chain='chain.1.txt'),
            'interleaved_groups': lambda s: s['groups'].__setitem__(100, 1),
            'fractional_position': lambda s: s['locations'][100].update(expanded_index=200.5),
            'repeated_slot_wrong_point': lambda s: s['points'][567].update(x=999.),
            'repeated_slot_wrong_RLE_row': lambda s: s['locations'][567].update(row=67),
            'repeated_slot_wrong_stored_density': lambda s: s['stored_logposts'].__setitem__(567, 999.),
            'missing_location': lambda s: s['locations'].pop(),
        }
        for name, change in changes.items():
            modified = copy.deepcopy(selection); change(modified)
            controls[name] = refused(lambda: consumer.chronology(modified))
        strict = copy.deepcopy(selection)
        for i in range(2000):
            strict['locations'][i]['expanded_index'] = i%500
            strict['locations'][i]['row'] = i%500
        strict_report = consumer.chronology(strict)
        assert sum(x['repeated_expanded_index_slots'] for x in strict_report.values()) == 0
        controls['strictly_increasing_positive_control'] = True
        with patch.object(measurement_summary, 'summarize_run', side_effect=AssertionError('Unqualified parent')):
            controls['failed_parent_before_selection_read'] = refused(lambda: consumer.actual(Path(temporary)/'missing', Path(temporary)/'missing-summary'))
        with patch.object(measurement_summary, 'summarize_run', qualifier):
            native_path = folder/'selected/00000.json'; original_bytes = native_path.read_bytes()
            native_path.write_text('{}\n')
            controls['changed_native_record'] = refused(lambda: consumer.actual(folder, summary))
            native_path.write_bytes(original_bytes)
            parent['posterior']['x']['quantiles_025_16_50_84_975'][0] += .01
            controls['qualified_baseline_mismatch'] = refused(lambda: consumer.actual(folder, summary))
            parent['posterior']['x']['quantiles_025_16_50_84_975'][0] -= .01
        # Merging equal observations alters the declared midpoint interpolation estimator.
        xx = np.array([0., 0., 1., 4.]); ww = np.array([.1, .1, .3, .5])
        explicit = original.weighted_quantiles(xx, ww)
        merged = original.weighted_quantiles(np.array([0., 1., 4.]), np.array([.2, .3, .5]))
        assert np.max(abs(explicit-merged)) > .01
        merge_difference = float(np.max(abs(explicit-merged)))
    sources = [Path(__file__), Path(consumer.__file__), Path(original.__file__), consumer.HERE/'measurement_summary.py']
    return {'schema': 'ordered-quantile-consumer-synthetic-validation-v1', 'status': 'passed',
        'scope': 'Synthetic input files and explicit mocked parent qualifier; actual cosmological source qualification is delegated unchanged to measurement_summary.summarize_run.',
        'physical_model_calls': 0, 'background_calls': 0, 'observational_consumer_runs': 0,
        'selected_slots': 2000, 'chain_count': 4, 'repeated_expanded_index_slots': 2,
        'unique_expanded_indices_per_chain': [500, 499, 499, 500],
        'bootstrap_partitions': [10, 20], 'replicates_per_partition': 500,
        'all_original_bootstrap_fields_bitwise_equal': True,
        'independent_arithmetic_max_abs_errors': maxima,
        'merging_repeated_observations_changes_midpoint_quantiles_max_abs': merge_difference,
        'controls': controls,
        'source_sha256': {consumer.relative(p): original.digest(p) for p in sources}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUTPUT)
    args = parser.parse_args()
    import camb
    import cobaya
    import cobaya.model
    def forbidden(*args, **kwargs):
        raise AssertionError('Physical calculations forbidden')
    with ExitStack() as stack:
        for module, name in [(camb, 'get_background'), (camb, 'get_results'), (camb, 'get_transfer_functions'),
                             (cobaya, 'get_model'), (cobaya.model, 'get_model')]:
            stack.enter_context(patch.object(module, name, forbidden))
        report = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'status': report['status'], 'output': str(args.output), 'sha256': original.digest(args.output)}))


if __name__ == '__main__':
    main()
