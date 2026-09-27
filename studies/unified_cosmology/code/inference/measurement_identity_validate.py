"""Verify reporting-time scientific identity checks without evaluating cosmology."""
import argparse
import copy
import importlib
import json
from pathlib import Path
from unittest.mock import patch

import measurement_summary as measurement


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, action='append', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    import camb
    import target_identity
    tested = []
    def no_computation(*args, **kwargs):
        raise AssertionError('Identity validation attempted a cosmological calculation.')
    with patch.object(camb, 'get_results', no_computation), \
         patch.object(camb, 'get_background', no_computation), \
         patch.object(camb, 'get_transfer_functions', no_computation):
        for path in args.manifest:
            manifest = json.loads(path.read_text())
            intact = measurement.verify_current_target(manifest)
            frozen = manifest['target_identity']
            settings = manifest['arguments']
            backend = importlib.import_module('modern_gpu' if settings.get('gpu') else
                                             'modern_fast' if settings.get('fast_lensing') else
                                             'target_identity')
            failures = []
            def rejected(name, changed=manifest):
                try:
                    measurement.verify_current_target(changed)
                except AssertionError:
                    failures.append(name)
                else:
                    raise AssertionError('Identity mutation was accepted: '+name)
            source = measurement.ROOT / next(iter(frozen['source_sha256']))
            sample = Path(frozen['configuration']['likelihood']['released_sn']['data_file'])
            model = Path(settings['surrogate'])
            if not model.is_absolute():
                model = measurement.ROOT / model
            original_digest = measurement.digest
            for name, target in [('source_bytes', source), ('sample_bytes', sample),
                                 ('surrogate_bytes', model)]:
                def altered_digest(path, target=target):
                    return '0'*64 if Path(path).resolve() == target.resolve() else original_digest(path)
                with patch.object(measurement, 'digest', altered_digest):
                    rejected(name)
            changed = copy.deepcopy(manifest)
            package = next(iter(changed['target_identity']['versions']))
            changed['target_identity']['versions'][package] = 'invalid-version-for-test'
            rejected('package_version', changed)
            # The backend must actually rehash released assets: changing a
            # reported asset digest must fail even when manifests are intact.
            original_asset_digest = target_identity.digest
            primary_inputs = {str((measurement.ROOT/name).resolve()) for name in intact}
            def altered_asset_digest(path):
                if str(Path(path).resolve()) not in primary_inputs:
                    return '0'*64
                return original_asset_digest(path)
            with patch.object(target_identity, 'digest', altered_asset_digest):
                rejected('third_party_likelihood_bytes')
            # Exercise the final complete-record comparison independently of
            # the earlier source/sample/package checks.
            altered_identity = copy.deepcopy(frozen)
            altered_identity['identity'] = '0'*64
            with patch.object(backend, 'identify', return_value=altered_identity):
                rejected('recomputed_target_identity')
            assert len(failures) == 6
            tested.append({'manifest': measurement.relative(path),
                           'manifest_sha256': measurement.digest(path),
                           'target_identity': frozen['identity'],
                           'backend': backend.__name__, 'intact_identity_matches': True,
                           'verified_first_party_and_manifest_inputs': len(intact),
                           'rejected_changes': failures})
    result = {'status': 'passed', 'cases': tested,
              'native_spectrum_calls': 0, 'background_calls': 0,
              'mutations': 'Digest/version/backend substitutions isolated to the validator; no source, observations, model or active chain files changed.',
              'scope': 'Current scientific source, input, dependency and runtime identity at reporting. This does not qualify unfinished chains or validate physical assumptions.',
              'source_sha256': {measurement.relative(p): measurement.digest(p) for p in
                                [Path(__file__), Path(measurement.__file__)]}}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': result['status'], 'targets': len(tested),
                      'rejected_mutations': sum(len(row['rejected_changes']) for row in tested)}))


if __name__ == '__main__':
    main()
