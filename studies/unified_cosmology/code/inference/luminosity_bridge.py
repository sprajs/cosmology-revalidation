"""Exact SN-density bridge between declared Gaussian luminosity-prior targets."""
import argparse
import copy
import hashlib
import importlib.metadata
import json
from pathlib import Path
import time

import numpy as np
from scipy.special import logsumexp

from luminosity_sensitivity import IntegratedLuminosity, exact_background_record, weight_diagnostics
from measurement_summary import summarize_run, weighted_fraction
from late_geometry import sample_path
from target_identity import canonical

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
DESIGN = HERE/'luminosity-bridge-design.json'
GATES = HERE/'luminosity-sensitivity-design.json'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


def identity(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def sn_loglike(background, evolution):
    value = background['baseline_SN_loglike']
    if evolution != 'none':
        value += background['log_likelihood_ratios'][evolution]
    return float(value)


def bridge_record(parent, background, source, target, tolerance):
    source_loglike = sn_loglike(background, source)
    target_loglike = sn_loglike(background, target)
    closure = source_loglike-parent['exact_loglikes']['released_sn']
    delta = target_loglike-source_loglike
    total = float(parent['log_weight']+delta)
    result = {'source_SN_loglike': source_loglike, 'target_SN_loglike': target_loglike,
              'source_SN_loglike_closure': float(closure), 'SN_log_likelihood_ratio': float(delta),
              'parent_exact_proposal_logweight': parent['log_weight'], 'target_logweight': total}
    if not all(np.isfinite(value) for value in result.values()):
        raise ArithmeticError('Nonfinite bridge density or weight.')
    result['status'] = 'finite_bridge' if abs(closure) <= tolerance else 'source_SN_density_mismatch'
    return result


def summarize_bridge(parent_records, rows, groups, gates):
    failures = [{'index': index, 'status': row['status'], 'error': row.get('error')}
                for index, row in enumerate(rows) if row['status'] != 'finite_bridge']
    if failures:
        return {'status': 'failed_bridge_evaluation', 'failures': failures, 'posterior': None}
    lw = np.array([row['target_logweight'] for row in rows])
    parent_lw = np.array([row['log_weight'] for row in parent_records])
    names = sorted(set(parent_records[0]['point']) | set(parent_records[0]['derived']))
    values = {name: np.array([dict(row['point'], **row['derived'])[name] for row in parent_records]) for name in names}
    diagnostic = weight_diagnostics(lw, values, np.asarray(groups), gates)
    failures = list(diagnostic['failed_gates'])
    if len(rows) < gates['minimum_exact_points']:
        failures.append('minimum_exact_points')
    weights = np.exp(lw-logsumexp(lw))
    parent_weights = np.exp(parent_lw-logsumexp(parent_lw))
    matrix = np.column_stack([values[name] for name in names])
    centered = matrix-weights@matrix
    covariance = (centered*weights[:, None]).T@centered
    signs = None
    if not failures:
        signs = {name+'_below_zero': weighted_fraction(values[name] < 0., weights, np.asarray(groups))
                 for name in ['q0', 'q05', 'q1', 'j0'] if name in values}
    return {'status': 'qualified_conditional_bridge' if not failures else 'insufficient_bridge_overlap_or_stability',
            'qualified_under_declared_numerical_gates': not failures,
            'points': len(rows), 'failed_gates': failures, 'weight_diagnostics': diagnostic,
            'posterior': diagnostic['weighted_summaries_for_diagnostics'] if not failures else None,
            'diagnostic_weighted_summaries': diagnostic['weighted_summaries_for_diagnostics'],
            'conditional_sign_fractions': signs,
            'weighted_covariance': {'parameter_order': names, 'matrix': covariance.tolist(),
                                    'normalization': 'Untrimmed normalized posterior weights, no Bessel correction.',
                                    'qualified': not failures},
            'maximum_normalized_weight_change_from_parent': float(np.max(abs(weights-parent_weights))),
            'SN_log_likelihood_ratio_quantiles': np.quantile([r['SN_log_likelihood_ratio'] for r in rows], [0, .025, .5, .975, 1]).tolist(),
            'source_SN_loglike_max_absolute_closure': max(abs(r['source_SN_loglike_closure']) for r in rows)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--chain-folder', type=Path, required=True)
    parser.add_argument('--correction-summary', type=Path, required=True)
    parser.add_argument('--target', choices=['smooth01', 'smooth03'], default='smooth03')
    parser.add_argument('--cache', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    parent = summarize_run(args.chain_folder, args.correction_summary)
    summary = json.loads(args.correction_summary.read_text())
    selection_path = ROOT/summary['selection_path']
    assert selection_path.parent.parent == args.chain_folder.resolve()
    selection = json.loads(selection_path.read_text())
    settings = selection['settings']
    design = json.loads(DESIGN.read_text()); gates = json.loads(GATES.read_text())['overlap_gates']
    assert settings['evolution'] in design['allowed_sources'], 'An explicit linear-coefficient parent cannot use this marginal bridge.'
    manifest = json.loads((args.chain_folder/'run-0.json').read_text())
    frozen = manifest['target_identity']
    for path, expected in frozen['source_sha256'].items():
        assert digest(ROOT/path) == expected, 'Scientific source changed.'
    for name, version in frozen['versions'].items():
        assert importlib.metadata.version(name) == version, 'Scientific environment changed.'
    if settings.get('fast_lensing'):
        from modern_fast import configuration
    else:
        from modern_run import configuration
    options = {key: settings[key] for key in ['model', 'evolution', 'sample', 'calibration']}
    source_configuration = configuration(**options)
    # A bridge may change only the analytically integrated luminosity covariance.
    target_options = dict(options, evolution=args.target)
    target_configuration = configuration(**target_options)
    changed = copy.deepcopy(canonical(source_configuration))
    changed['likelihood']['released_sn']['smooth_sigma'] = design['Gaussian_coefficient_priors_mag'][args.target]
    assert changed == canonical(target_configuration), 'A supposedly cancelling non-SN factor changed.'
    proposal_configuration = configuration(**options, surrogate=Path(settings['surrogate']))
    assert canonical(proposal_configuration) == frozen['configuration'], 'Parent target configuration changed.'
    data_path = sample_path(settings['sample'])
    assert digest(data_path) == frozen['sample_sha256']
    with np.load(data_path, allow_pickle=False) as data:
        sn = IntegratedLuminosity(*(data[key] for key in ['zHD', 'zHEL', 'MU', 'covariance']))
    target_description = {'native_configuration': canonical(target_configuration),
                          'source_sha256': frozen['source_sha256'], 'versions': frozen['versions'],
                          'likelihood_assets': frozen['assets'], 'sample_sha256': frozen['sample_sha256']}
    target_description['identity'] = identity(target_description)
    sources = [Path(__file__), DESIGN, GATES, HERE/'luminosity_sensitivity.py',
               HERE/'measurement_summary.py', HERE/'exact_correction.py', HERE/'late_geometry.py',
               HERE/'modern_run.py', HERE/'modern_fast.py', HERE/'target_identity.py']
    hashes = {relative(path): digest(path) for path in sources}
    lineage = {'qualified_parent_inputs': parent['input_sha256'], 'bridge_source_sha256': hashes,
               'parent_proposal_target_identity': frozen['identity'], 'native_target': target_description,
               'target_evolution': args.target, 'source_evolution': settings['evolution'],
               'SN_data_path': relative(data_path), 'SN_data_sha256': digest(data_path)}
    cache_identity = identity(lineage)
    cache = args.cache.resolve(); assert cache.is_relative_to(ROOT/'.work')
    cache.mkdir(parents=True, exist_ok=True)
    lineage_path = cache/'lineage.json'
    if lineage_path.exists():
        assert json.loads(lineage_path.read_text()) == lineage, 'Different bridge cache identity.'
    else:
        lineage_path.write_text(json.dumps(lineage, indent=2)+'\n')
    cache_manifest = cache/'bridge-records.json'
    if cache_manifest.exists():
        for path, expected in json.loads(cache_manifest.read_text()).items():
            assert digest(ROOT/path) == expected, 'Previously recorded bridge cache changed.'
    parent_records = []; rows = []; cache_hashes = {}; started = time.monotonic()
    for index, point in enumerate(selection['points']):
        native_path = selection_path.parent/f'{index:05d}.json'
        native = json.loads(native_path.read_text()); parent_records.append(native)
        path = cache/f'{index:05d}.json'
        if path.exists():
            row = json.loads(path.read_text()); checksum = row.pop('payload_sha256')
            assert identity(row) == checksum, 'Bridge cache payload changed.'
            row['payload_sha256'] = checksum
            assert row['identity'] == cache_identity and row['index'] == index
            assert row['native_record_sha256'] == digest(native_path)
        else:
            try:
                background = exact_background_record(point, settings, sn)
                json.dumps(background, allow_nan=False)
                row = bridge_record(native, background, settings['evolution'], args.target,
                                    design['source_SN_loglike_closure_absolute_tolerance'])
                row['background'] = background
            except Exception as error:
                row = {'status': 'exception', 'error': repr(error)}
            row.update(identity=cache_identity, index=index, native_record_sha256=digest(native_path), point=point)
            row['payload_sha256'] = identity(row)
            temp = path.with_suffix('.part'); temp.write_text(json.dumps(row, indent=2, allow_nan=False)+'\n'); temp.replace(path)
        rows.append(row); cache_hashes[relative(path)] = digest(path)
        if (index+1) % 100 == 0:
            print(json.dumps({'bridge_points': index+1, 'seconds': time.monotonic()-started}), flush=True)
    cache_manifest.write_text(json.dumps(cache_hashes, indent=2)+'\n')
    result = summarize_bridge(parent_records, rows, selection['groups'], gates)
    result.update(source_settings=parent['settings'], target_settings=target_options,
                  parent_proposal_target_identity=frozen['identity'], native_target=target_description,
                  bridge_identity=cache_identity, source_sha256=hashes,
                  lineage_path=relative(lineage_path), lineage_sha256=digest(lineage_path),
                  cache_manifest_path=relative(cache_manifest), cache_manifest_sha256=digest(cache_manifest),
                  correction_summary=relative(args.correction_summary), correction_summary_sha256=digest(args.correction_summary),
                  qualification='A separately identified conditional target via exact SN density ratios, not a target chain or new native-correction campaign. Finite overlap and chain checks cannot establish unseen mode coverage. Gaussian luminosity priors are not measured age corrections.',
                  CMB_spectrum_calls=0)
    for mapping in [hashes, frozen['source_sha256'], parent['input_sha256'], {relative(data_path): frozen['sample_sha256']}]:
        assert all(digest(ROOT/path) == expected for path, expected in mapping.items()), 'Bridge input/source changed during evaluation.'
    assert all(importlib.metadata.version(name) == version for name, version in frozen['versions'].items())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'status': result['status'], 'source': settings['evolution'], 'target': args.target}))


if __name__ == '__main__':
    main()
