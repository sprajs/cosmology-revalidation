"""Run unchanged native density correction while retaining its computed spectra.

This optional companion changes neither the correction selection nor its density,
weights or qualification. Sidecars permit later lensing comparisons to reuse the
native spectra. An absent/failed sidecar never triggers an extra native call here.
"""
from __future__ import annotations
import hashlib
import json
import os
from pathlib import Path

import numpy as np

import exact_correction as core

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
_evaluate = core.evaluate
_initialize = core.initialize
_SOURCE_ENV = 'COSMOLOGY_SPECTRAL_CAPTURE_SOURCES'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


def dependencies():
    return {relative(path): digest(path) for path in [Path(__file__),
            HERE/'exact_correction.py', HERE/'joint_lensing_bridge.py',
            HERE/'joint-lensing-bridge-design.json']}


_START_SOURCES = json.loads(os.environ[_SOURCE_ENV]) if _SOURCE_ENV in os.environ else dependencies()


def check_startup_sources():
    assert dependencies() == _START_SOURCES, 'Capture implementation changed after driver startup.'


def initialize_with_capture(settings):
    check_startup_sources()
    _initialize(settings)
    check_startup_sources()


def sidecar_path(native_path):
    native_path = Path(native_path)
    return native_path.parent/'spectra'/native_path.name


def read_sidecar(native_path, require_current_sources=True):
    native_path = Path(native_path).resolve()
    parent = json.loads(native_path.read_text()); core.verify_record(parent)
    path = sidecar_path(native_path)
    row = json.loads(path.read_text()); core.verify_record(row)
    assert row['parent_native_record_sha256'] == digest(native_path)
    assert row['parent_payload_sha256'] == parent['payload_sha256']
    assert row['index'] == parent['index'] and row['point'] == parent['point']
    assert row['correction_identity'] == parent['target_identity']
    if require_current_sources:
        assert row['source_sha256'] == dependencies(), 'Spectral capture source changed.'
    if row['status'] == 'captured_from_parent_native_evaluation':
        assert parent['status'] == 'finite'
        assert row['CLIPY_NOJAX'] == '1'
        assert digest(ROOT/row['spectra_path']) == row['spectra_sha256']
        with np.load(ROOT/row['spectra_path'], allow_pickle=False) as spectra:
            validate_spectra(spectra)
    return row


def validate_spectra(spectra):
    assert set(spectra) == {'ell', 'tt', 'ee', 'te', 'bb', 'pp'}
    ell = np.asarray(spectra['ell'])
    assert ell.ndim == 1 and len(ell) >= 3102
    assert np.array_equal(ell, np.arange(len(ell)))
    assert all(np.asarray(spectra[k]).shape == ell.shape and np.isfinite(spectra[k]).all()
               for k in spectra)


def evaluate_with_spectra(task):
    check_startup_sources()
    index, point, destination, correction_identity, _ = task
    native = Path(destination)/f'{index:05d}.json'
    existed = native.exists()
    row = _evaluate(task)
    check_startup_sources()
    # The original function remains the sole producer of numerical corrections.
    # In particular, the optional cache cannot alter a failure or its weights.
    stored = json.loads(native.read_text()); core.verify_record(stored)
    assert core.record_digest(row) == stored['payload_sha256']
    sidecar = sidecar_path(native)
    if sidecar.exists():
        read_sidecar(native)
        return row
    sidecar.parent.mkdir(exist_ok=True)
    capture = {'index': index, 'point': point, 'correction_identity': correction_identity,
               'parent_native_record_sha256': digest(native),
               'parent_payload_sha256': row['payload_sha256'],
               'source_sha256': dict(_START_SOURCES), 'CLIPY_NOJAX': os.environ.get('CLIPY_NOJAX'),
               'additional_native_evaluations': 0}
    if existed:
        capture['status'] = 'unavailable_from_preexisting_native_record'
    elif row['status'] != 'finite':
        capture['status'] = 'unavailable_parent_not_finite'
    else:
        try:
            assert capture['CLIPY_NOJAX'] == '1', 'Use the baseline NumPy likelihood backend.'
            # The exact and proposal Cobaya models own separate provider states.
            # Independent validation checks this after proposal evaluation too.
            raw = core._exact.provider.get_Cl(ell_factor=False, units='FIRASmuK2')
            spectra = {k: np.array(raw[k], copy=True) for k in ['ell', 'tt', 'ee', 'te', 'bb', 'pp']}
            validate_spectra(spectra)
            from joint_lensing_bridge import source_normalization, DESIGN
            constants = source_normalization(core._exact, json.loads(DESIGN.read_text()))
            assert np.isfinite(list(constants.values())).all()
            archive = sidecar.with_suffix('.npz')
            assert not archive.exists(), 'Preserve an orphaned spectral archive; no overwrite.'
            temporary = archive.with_suffix('.part')
            with temporary.open('xb') as stream:
                np.savez_compressed(stream, **spectra)
            temporary.replace(archive)
            capture.update(status='captured_from_parent_native_evaluation',
                           spectra_path=relative(archive), spectra_sha256=digest(archive),
                           source_normalized_gaussian_constants=constants,
                           spectral_units='raw C_l: TT/EE/TE/BB FIRASmuK2, pp dimensionless potential')
        except Exception as error:
            capture.update(status='capture_failed_parent_density_retained', error=repr(error))
    capture['payload_sha256'] = core.record_digest(capture)
    temporary = sidecar.with_suffix('.part')
    temporary.write_text(json.dumps(capture, indent=2, allow_nan=False)+'\n')
    temporary.replace(sidecar)
    check_startup_sources()
    read_sidecar(native)
    return row


def main():
    assert os.environ.get('CLIPY_NOJAX') == '1'
    check_startup_sources()
    old = os.environ.get(_SOURCE_ENV)
    # Spawned workers receive the parent's snapshot, not a new post-evaluation
    # reading of possibly edited source files.
    os.environ[_SOURCE_ENV] = json.dumps(_START_SOURCES, sort_keys=True)
    core.initialize = initialize_with_capture
    core.evaluate = evaluate_with_spectra
    try:
        core.main()
        check_startup_sources()
    finally:
        if old is None:
            os.environ.pop(_SOURCE_ENV, None)
        else:
            os.environ[_SOURCE_ENV] = old


if __name__ == '__main__':
    main()
