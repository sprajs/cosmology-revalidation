"""Synthetic algebra, chronology and immutable-cache tests; no model calls."""
import argparse
from contextlib import ExitStack
import copy
import importlib.metadata
import json
from pathlib import Path
import tempfile
from unittest.mock import patch

import numpy as np
from scipy.integrate import quad
from scipy.stats import multivariate_normal, multivariate_t, norm

import anchored_proposal as proposal

SEED = 273731


def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def refuses(call, message=None):
    try:
        call()
    except (AssertionError, ValueError, KeyError, FileNotFoundError) as error:
        if message is not None:
            assert message in str(error), repr(error)
        return type(error).__name__
    raise AssertionError('Invalid synthetic fixture was accepted.')


def algebra():
    rng = np.random.default_rng(SEED)
    cases = []
    for case, dimension in enumerate([1, 2, 5, 13, 15]*2):
        nchain = 500+25*(case % 3)
        n = 4*nchain
        a = rng.normal(size=(dimension, dimension))
        covariance = a@a.T+.3*np.eye(dimension)
        reference = proposal.FrozenMixture(rng.normal(size=dimension), covariance,
                                          [f'x{k}' for k in range(dimension)])
        points = rng.normal(size=(n, dimension))@reference.L.T+reference.mean
        groups = np.repeat(np.arange(4), nchain)
        chronology = np.tile(np.arange(nchain)*3+7, 4)
        logweights = rng.normal(0, .7, n)
        shuffle = rng.permutation(n)
        points, groups, chronology, logweights = [v[shuffle] for v in [points, groups, chronology, logweights]]
        fitted, mask, stats = proposal.fit(points, logweights, groups, chronology, reference)
        independent_mask = np.zeros(n, bool)
        for group in range(4):
            ordered = sorted([i for i in range(n) if groups[i] == group], key=lambda i: chronology[i])
            independent_mask[ordered[:(4*len(ordered))//5]] = True
        assert np.array_equal(mask, independent_mask)
        assert min(x['withheld'] for x in stats['partition']) >= 100
        raw = np.exp(logweights[mask]-max(logweights[mask]))
        weights = raw/raw.sum()
        mean = np.sum(points[mask]*weights[:, None], axis=0)
        centered = points[mask]-mean
        # Independent physical-coordinate identity, avoiding the producer's
        # whiten/recolour operations and triangular solve.
        cov = .95*np.einsum('i,ij,ik->jk', weights, centered, centered)+.05*covariance
        mean_error = float(np.max(abs(fitted.mean-mean)))
        covariance_error = float(np.max(abs(fitted.covariance-cov)))
        white = np.linalg.solve(reference.L, (points[mask]-reference.mean).T).T
        wc = white-np.sum(weights[:, None]*white, axis=0)
        eigen = np.linalg.eigvalsh(np.einsum('i,ij,ik->jk', weights, wc, wc))
        eigen_error = float(np.max(abs(eigen-stats['empirical_white_covariance_eigenvalues'])))
        assert max(mean_error, covariance_error, eigen_error) < 2e-10
        assert abs(stats['raw_training_weight_ESS']-1/(weights@weights)) < 2e-10
        altered = points.copy(); altered[~mask] = rng.normal(size=((~mask).sum(), dimension))*1e40
        altered_weights = logweights.copy(); altered_weights[~mask] = rng.normal(size=(~mask).sum())*1e4
        heldout, new_mask, new_stats = proposal.fit(altered, altered_weights, groups, chronology, reference)
        assert np.array_equal(heldout.mean, fitted.mean) and np.array_equal(heldout.covariance, fitted.covariance)
        assert np.array_equal(mask, new_mask) and new_stats == stats
        shifted, _, shifted_stats = proposal.fit(points, logweights+1000., groups, chronology, reference)
        offset_error = max(float(np.max(abs(shifted.mean-fitted.mean))),
                           float(np.max(abs(shifted.covariance-fitted.covariance))))
        assert offset_error < 2e-10
        test = rng.normal(size=(256, dimension))@reference.L.T*3+reference.mean
        independent_pdf = np.logaddexp(
            np.log(.9)+multivariate_normal.logpdf(test, fitted.mean, fitted.covariance*1.05**2),
            np.log(.1)+multivariate_t.logpdf(test, fitted.mean, fitted.covariance*2.4, df=5))
        density_error = float(np.max(abs(fitted.logpdf(test)-independent_pdf)))
        assert density_error < 2e-10
        mass_error = None
        if dimension == 1:
            mass = quad(lambda x: np.exp(fitted.logpdf([x])), -np.inf, np.inf, epsabs=1e-10)[0]
            mass_error = abs(mass-1.)
            assert mass_error < 1e-9
        cases.append(dict(dimension=dimension, points=n, mean_error=mean_error,
                          covariance_error=covariance_error, white_eigenvalue_error=eigen_error,
                          common_logweight_offset_error=offset_error, mixture_logdensity_error=density_error,
                          scalar_integrated_mass_error=mass_error, heldout_changes_leave_fit_bitwise_equal=True))
    failures = {}
    duplicate = chronology.copy()
    ids = np.flatnonzero(groups == 0); duplicate[ids[1]] = duplicate[ids[0]]
    failures['duplicate_chronology'] = refuses(lambda: proposal.fit(points, logweights, groups, duplicate, reference))
    failures['omitted_chain'] = refuses(lambda: proposal.partition(np.minimum(groups, 2), chronology))
    failures['fewer_than_2000'] = refuses(lambda: proposal.fit(points[:1999], logweights[:1999], groups[:1999], chronology[:1999], reference))
    invalid = chronology.astype(float); invalid[0] += .2
    failures['noninteger_chronology'] = refuses(lambda: proposal.partition(groups, invalid))
    concentrated = np.full(len(logweights), -1000.); concentrated[np.flatnonzero(mask)[0]] = 0.
    failures['insufficient_training_ESS'] = refuses(lambda: proposal.fit(points, concentrated, groups, chronology, reference), 'Insufficient')
    invalid = points.copy(); invalid[0, 0] = np.nan
    failures['nonfinite_coordinates'] = refuses(lambda: proposal.fit(invalid, logweights, groups, chronology, reference))
    invalid = logweights.copy(); invalid[0] = -np.inf
    failures['nonfinite_logweights'] = refuses(lambda: proposal.fit(points, invalid, groups, chronology, reference))
    return {'seed': SEED, 'cases': cases, 'invalid_inputs_refused': failures}


def consumer_fixture(directory):
    """Run actual cache replay and freeze; mock only qualification/physics I/O.

    Synthetic native/bridge records are created without invoking any background.
    The real parent diagnostics are deliberately not certified by this fixture.
    """
    import anchored_bridge as bridge
    from exact_correction import record_digest
    rng = np.random.default_rng(SEED+1)
    n = 2000; groups = np.repeat(np.arange(4), 500)
    x = rng.normal(size=n); x[groups == 3] += 1.5
    auxiliary = rng.normal(size=n)
    settings = dict(model='lcdm', evolution='none', sample='dovekie', calibration='official_planck')
    params = {k: {'prior': {'min': -10., 'max': 10.}} for k in ['w', 'A_fg']}
    folder = directory/'chain'; native_dir = folder/'exact'; native_dir.mkdir(parents=True)
    cache = directory/'bridge'; cache.mkdir()
    fake_data = directory/'data-placeholder'; fake_data.write_text('Synthetic fixture only, no observational values.')
    frozen = {'identity': 'synthetic-parent', 'source_sha256': {}, 'assets': {'synthetic': True},
              'versions': {'numpy': importlib.metadata.version('numpy')}, 'sample_sha256': proposal.digest(fake_data)}
    target = {'identity': 'synthetic-anchored', 'assets': frozen['assets'], 'versions': frozen['versions'],
              'source_sha256': {}, 'anchored_calibration': {'input_and_audit_sha256': {}},
              'configuration': {'params': params}}
    native_records = []
    for i in range(n):
        point = dict(w=float(x[i]), A_fg=float(auxiliary[i]))
        ll = dict(released_sn=float(norm.logpdf(x[i], .2, 1.)), synthetic_CMB_BAO=float(norm.logpdf(x[i])))
        prior = float(norm.logpdf(auxiliary[i]))
        record = dict(status='finite', point=point, derived=dict(q0=float(x[i]-.1), q05=float(x[i]*.8),
                     q1=float(x[i]*.5), j0=1., sn_chi2=float((x[i]-.2)**2)), exact_loglikes=ll,
                     proposal_loglikes=ll, exact_logpost=prior+sum(ll.values()),
                     proposal_logpost=prior+sum(ll.values()), log_weight=0.)
        record['payload_sha256'] = record_digest(record)
        dump(native_dir/f'{i:05d}.json', record); native_records.append(record)
    selection_path = native_dir/'selection.json'; correction_path = directory/'correction.json'
    dump(selection_path, dict(settings=settings, points=[r['point'] for r in native_records], groups=groups.tolist(),
                             locations=[{'expanded_index': int(i % 500)*2} for i in range(n)]))
    dump(correction_path, {'selection_path': proposal.relative(selection_path)})
    dump(folder/'run-0.json', {'target_identity': frozen})
    parent_paths = [folder/'run-0.json', selection_path, correction_path, *sorted(native_dir.glob('[0-9]*.json'))]
    parent = {'input_sha256': {proposal.relative(p): proposal.digest(p) for p in parent_paths}, 'settings': settings}
    paths = [Path(bridge.__file__), bridge.DESIGN, bridge.GATES,
             *[proposal.HERE/f for f in ['anchored_adapter.py','anchored-design.json','expansion_history.py',
               'luminosity_sensitivity.py','measurement_summary.py','probe_omission.py','exact_correction.py',
               'target_identity.py','modern_fast.py','modern_run.py']],
             proposal.HERE.parent/'distance_ladder/calibration_interface.py']
    hashes = {proposal.relative(p): proposal.digest(p) for p in paths}
    lineage = dict(qualified_parent_inputs=parent['input_sha256'], source_sha256=hashes,
                   parent_proposal_target_identity=frozen['identity'], native_target=target,
                   parent_settings=settings, parent_SN_sha256=proposal.digest(fake_data),
                   all_non_SN_factors_and_priors_identical=True)
    cache_identity = bridge.identity(lineage)
    bridge.payload(cache/'lineage.json', lineage)
    design = json.loads(bridge.DESIGN.read_text())
    manifest = {proposal.relative(cache/'lineage.json'): proposal.digest(cache/'lineage.json')}
    for i, record in enumerate(native_records):
        background = dict(old_SN_reconstructed_loglike=record['exact_loglikes']['released_sn'],
                          anchored_SN=dict(loglike=float(norm.logpdf(x[i], .35, .08)), chi2=float(((x[i]-.35)/.08)**2)),
                          background_calls=0, explicitly_synthetic=True)
        row = bridge.replace_SN(record, background, record['exact_loglikes'], design)
        row.update(background=background, identity=cache_identity, index=i, point=record['point'],
                   native_record_sha256=proposal.digest(native_dir/f'{i:05d}.json'))
        row['payload_sha256'] = bridge.identity(row)
        path = cache/f'{i:05d}.json'; bridge.payload(path, row)
        manifest[proposal.relative(path)] = proposal.digest(path)
    bridge.payload(cache/'record-hashes.json', manifest)
    reference_dir = directory/'reference'; reference_dir.mkdir()
    np.savez(reference_dir/'proposal.npz', mean=[0., 0.], cov=[[1., .2], [.2, 1.5]], names=list(params))
    snapshot = reference_dir/'synthetic-source'; snapshot.write_text('not a chain; synthetic reference')
    dump(reference_dir/'proposal.json', dict(source_sha256=proposal.digest(proposal.HERE/'independence_proposal.py'),
         design_sha256=proposal.digest(proposal.HERE/'independence-sampling-design.json'),
         snapshot_files={'synthetic-source': {'sha256': proposal.digest(snapshot)}},
         proposal_sha256=proposal.digest(reference_dir/'proposal.npz'), parent_target_identity=frozen['identity']))
    validation_path = directory/'explicitly-synthetic-validator-prerequisite.json'
    dump(validation_path, dict(status='passed_synthetic_anchored_proposal_validation',
         proposal_source_sha256=proposal.sources(), source_sha256={proposal.relative(Path(__file__)): proposal.digest(__file__)}))
    config = {'likelihood': {'released_sn': {'data_file': str(fake_data)}, 'synthetic_CMB_BAO': {}},
              'theory': {'camb': {'extra_args': {}}}}
    report_path = directory/'bridge-report.json'; output = directory/'trained'
    with ExitStack() as stack:
        stack.enter_context(patch.object(bridge, 'summarize_run', return_value=parent))
        stack.enter_context(patch.object(bridge, 'configurations', return_value=(config, config)))
        stack.enter_context(patch.object(bridge.anchored, 'identify', return_value=target))
        stack.enter_context(patch.object(bridge, 'load_SN_interfaces', return_value=(None, None)))
        stack.enter_context(patch.object(proposal, 'VALIDATION', validation_path))
        report = bridge.actual(folder, correction_path, cache)
        assert report['status'] == 'insufficient_anchored_overlap_or_stability'
        assert 'raw_weight_ESS' in report['failed_gates']
        assert report['posterior'] is None and report['source_SN_max_absolute_loglike_closure'] == 0.
        dump(report_path, report)
        frozen_record = proposal.freeze(report_path, reference_dir, output)
        assert frozen_record['bridge_status'] == report['status']
        assert frozen_record['statistics']['posterior_qualification'] is False
        assert frozen_record['training_settings']['sample'] == 'pantheon_shoes_anchored'
        loaded, metadata = proposal.load(output, {'params': params})
        assert metadata['training_target_identity'] == target['identity']
        refuses(lambda: proposal.freeze(report_path, reference_dir, output))
        # The real replay must call the parent qualifier, and may not promote an
        # unqualified parent merely because its cached bridge rows are complete.
        with patch.object(bridge, 'summarize_run', side_effect=AssertionError('unqualified-parent')):
            refuses(lambda: proposal.freeze(report_path, reference_dir, directory/'bad-parent'), 'unqualified-parent')
        assert not (directory/'bad-parent').exists()
        missing = cache/'00000.json'; original = missing.read_bytes(); missing.unlink()
        refuses(lambda: proposal.freeze(report_path, reference_dir, directory/'missing-record'))
        missing.write_bytes(original)
        assert not (directory/'missing-record').exists()
        saved = report_path.read_bytes(); changed = copy.deepcopy(report)
        changed['source_SN_max_absolute_loglike_closure'] = 1.
        dump(report_path, changed)
        refuses(lambda: proposal.freeze(report_path, reference_dir, directory/'nonreplaying-report'), 'does not reproduce')
        report_path.write_bytes(saved)
        changed = copy.deepcopy(report); changed['status'] = 'failed_anchored_bridge_evaluation'; dump(report_path, changed)
        refuses(lambda: proposal.freeze(report_path, reference_dir, directory/'failed-child'))
        report_path.write_bytes(saved)
        # A sealed row and refreshed file ledger cannot hide a density-closure
        # failure: actual replay must reject the otherwise allowed report.
        ledger_path = cache/'record-hashes.json'; saved_ledger = ledger_path.read_bytes()
        row = json.loads(original); row.pop('payload_sha256')
        row['status'] = 'source_SN_density_mismatch'; row['source_SN_loglike_closure'] = .01
        row['payload_sha256'] = bridge.identity(row); dump(missing, row)
        bad_ledger = json.loads(saved_ledger); bad_ledger[proposal.relative(missing)] = proposal.digest(missing)
        dump(ledger_path, bad_ledger)
        changed = copy.deepcopy(report); changed['cache_manifest_sha256'] = proposal.digest(ledger_path)
        dump(report_path, changed)
        refuses(lambda: proposal.freeze(report_path, reference_dir, directory/'bad-density-row'), 'does not reproduce')
        missing.write_bytes(original); ledger_path.write_bytes(saved_ledger); report_path.write_bytes(saved)
        assert not (directory/'bad-density-row').exists()

    record_path = output/'proposal.json'; initial_record = json.loads(record_path.read_text())
    tamper_tests = []
    for key, value in [('source_sha256', {}), ('training_target_identity', 'wrong'), ('bridge_status', 'qualified_conditional_anchored_bridge'),
                       ('training_settings', {}), ('training_parent_settings', {}), ('reference_training_target_identity', 'wrong'),
                       ('scope', 'qualified posterior'), ('bridge_report_path', 'unbound.json')]:
        changed = copy.deepcopy(initial_record); changed[key] = value; dump(record_path, changed)
        refuses(lambda: proposal.load(output, {'params': params})); tamper_tests.append(key)
        dump(record_path, initial_record)
    for path in [output/'proposal.npz', output/'training.npz', output/'reference-proposal.npz', cache/'00001.json']:
        saved = path.read_bytes(); path.write_bytes(saved+b'tampered')
        refuses(lambda: proposal.load(output, {'params': params})); tamper_tests.append(path.name)
        path.write_bytes(saved)
    # Reseal a changed training mask at the file-checksum layer: deterministic
    # reconstruction must still reject the omitted/extra training assignment.
    training = output/'training.npz'; original = training.read_bytes()
    with np.load(training, allow_pickle=False) as z:
        arrays = {key: np.array(z[key]) for key in z}
    arrays['training_mask'][0] = ~arrays['training_mask'][0]
    np.savez(training, **arrays)
    changed = copy.deepcopy(initial_record); changed['snapshot_sha256']['training.npz'] = proposal.digest(training)
    dump(record_path, changed)
    refuses(lambda: proposal.load(output, {'params': params})); tamper_tests.append('resealed_changed_training_partition')
    training.write_bytes(original); dump(record_path, initial_record)
    refuses(lambda: proposal.load(output, {'params': dict(reversed(list(params.items())))}))
    altered = copy.deepcopy(params); altered['w']['periodic'] = True
    refuses(lambda: proposal.load(output, {'params': altered}))
    proposal.load(output, {'params': params})
    return dict(points=n, training_points=frozen_record['statistics']['training_points'],
                withheld_points=frozen_record['statistics']['withheld_points'], child_status=report['status'],
                child_scientific_failed_gates=report['failed_gates'], training_weight_ESS=frozen_record['statistics']['raw_training_weight_ESS'],
                child_weight_ESS=report['weight_diagnostics']['raw_weight_ESS'],
                failed_child_remains_unqualified=True, complete_actual_bridge_cache_replayed=True,
                parent_qualification_and_physics_interfaces_explicitly_mocked=True,
                synthetic_validation_prerequisite_explicitly_mocked_by_temporary_record=True,
                unqualified_parent_missing_record_failed_child_nonreplaying_report_refused=True,
                resealed_source_density_failure_refused=True,
                fresh_output_only=True, tamper_tests_refused=tamper_tests,
                parameter_order_and_periodic_coordinates_refused=True)


def validate():
    import camb
    import cobaya.model
    import anchored_bridge as bridge
    before = proposal.sources()
    attempted = []
    def forbidden(*args, **kwargs):
        attempted.append(True)
        raise AssertionError('This validator forbids scientific models and backgrounds.')
    with ExitStack() as stack:
        for owner, name in [(camb, 'set_params'), (camb, 'get_background'), (camb, 'get_results'),
                            (camb, 'get_transfer_functions'), (cobaya.model, 'get_model'),
                            (cobaya.model.Model, '__init__'), (bridge, 'background_densities')]:
            stack.enter_context(patch.object(owner, name, forbidden))
        mathematics = algebra()
        with tempfile.TemporaryDirectory(prefix='anchored-proposal-synthetic-', dir=proposal.ROOT/'.work') as directory:
            consumer = consumer_fixture(Path(directory))
    assert attempted == [] and proposal.sources() == before
    dependencies = [Path(__file__), Path(bridge.__file__), bridge.DESIGN, bridge.GATES,
                    proposal.HERE/'exact_correction.py', proposal.HERE/'probe_omission.py',
                    proposal.HERE/'luminosity_sensitivity.py', proposal.HERE/'measurement_summary.py']
    hashes = dict(before)
    hashes.update({proposal.relative(p): proposal.digest(p) for p in dependencies})
    return dict(status='passed_synthetic_anchored_proposal_validation', algebra=mathematics, consumer=consumer,
                proposal_source_sha256=before, source_sha256=hashes,
                native_models_and_backgrounds_globally_forbidden=True, actual_background_calls=0,
                CMB_spectrum_calls=0, observational_bridge_points=0, sampling_steps=0,
                limitations='Synthetic proposal mathematics and cache-routing tests only. Physics inputs and parent qualification are mocked in the consumer fixture; real freeze must freshly replay a complete qualified-parent bridge. A trained proposal is never an admitted posterior.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=proposal.VALIDATION)
    args = parser.parse_args()
    result = validate(); dump(args.output, result)
    print(json.dumps({'status': result['status'], 'source_sha256': result['source_sha256'],
                      'consumer': result['consumer']}))


if __name__ == '__main__':
    main()
