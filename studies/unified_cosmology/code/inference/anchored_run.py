"""Prepare immutable anchored target configurations; no sampling/evaluation CLI."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path

from anchored_adapter import (ROOT, HERE, DATA, DESIGN, SAMPLE, CALIBRATION_RECORD,
                              configuration, identify, canonical, digest, relative)


def prepare(output, model, surrogate=None, validation=None):
    output = Path(output).resolve()
    assert output.is_relative_to(ROOT/'.work') and not output.exists(), 'Fresh ignored preparation directory required.'
    environment = {key: os.environ.get(key) for key in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'CLIPY_NOJAX']}
    assert all(value == '1' for value in environment.values()), 'Use the unchanged one-thread NumPy scientific runtime.'
    if validation is None:
        validation = ROOT/'studies/unified_cosmology/results/inference/anchored-validation.json'
    validation = Path(validation).resolve()
    checked = json.loads(validation.read_text())
    assert checked['status'] == 'passed_configuration_identity_and_provider_checks_no_inference'
    for path, expected in checked['source_sha256'].items():
        assert digest(ROOT/path) == expected, 'Plumbing validation is stale: '+path
    assert checked['calibration_audit_sha256'] == digest(CALIBRATION_RECORD)
    info = configuration(model=model, surrogate=surrogate)
    identity = identify(info, DATA, surrogate)
    manifest = {'schema': 'prepared-anchored-target-v1', 'created_utc': datetime.now(timezone.utc).isoformat(),
        'arguments': {'model': model, 'evolution': 'none', 'sample': SAMPLE, 'calibration': 'official_planck',
                      'surrogate': str(Path(surrogate).resolve()) if surrogate else None, 'fast_lensing': True, 'gpu': False},
        'target_identity': identity, 'runtime_environment': environment,
        'preparation_source_sha256': {relative(p): digest(p) for p in [Path(__file__), HERE/'anchored_adapter.py', HERE/'anchored_validate.py', DESIGN]},
        'validation_path': relative(validation), 'validation_sha256': digest(validation),
        'status': 'prepared_only_no_model_evaluation_or_sampling',
        'native_or_background_calls': 0, 'chain_qualification': None,
        'routing_limit': 'Existing original-sample drivers/correction/measurement consumers must not consume this manifest. An explicitly anchored factory and qualifier are required before inference.'}
    output.mkdir(parents=True)
    for name, value in [('target.json', manifest), ('configuration.json', canonical(info))]:
        with (output/name).open('x') as stream:
            json.dump(value, stream, indent=2, allow_nan=False); stream.write('\n')
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['prepare'])
    parser.add_argument('--model', choices=['lcdm', 'cpl'], required=True)
    parser.add_argument('--surrogate', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args.output, args.model, args.surrogate)
    print(json.dumps({'status': result['status'], 'target_identity': result['target_identity']['identity']}))


if __name__ == '__main__':
    main()
