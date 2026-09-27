#!/usr/bin/env python3
"""Read-only evaluation audit of the released ACT/Planck/SPT joint lensing product.

Downloads only the small official joint/MUSE releases to ignored work storage.
Reuses ACT v1.2 response files without changing them. No CAMB calls or inference.
"""
from __future__ import annotations
import argparse
import gc
import hashlib
import io
import json
from pathlib import Path
import sys
import tarfile
import urllib.request
import warnings
import zipfile

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
WORK = ROOT / '.work/unified-cosmology/joint-lensing-availability'
RELEASES = {
    'spt_act_likelihood-1.0.tar.gz': (
        'https://lambda.gsfc.nasa.gov/data/suborbital/act_spt_joint/spt_act_likelihood-1.0.tar.gz',
        'c6948403fec7ff8c3a7d2ad17ea12ae53977623eaedf2369db6481256f014d56'),
    'muse_3g_like_march_2025.zip': (
        'https://lambda.gsfc.nasa.gov/data/suborbital/SPT/muse_3g_like_march_2025.zip',
        '4933ceb993ceb92e48b2191b8ccb1595d182bc4bc3862c1fbad61d5f82136e9b'),
}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(2**20), b''):
            h.update(b)
    return h.hexdigest()


def acquire(work):
    work.mkdir(parents=True, exist_ok=True)
    for name, (url, digest) in RELEASES.items():
        path = work / name
        if not path.exists():
            with urllib.request.urlopen(url, timeout=60) as r:
                payload = r.read()
            if hashlib.sha256(payload).hexdigest() != digest:
                raise ValueError(f'Official release changed: {url}')
            path.write_bytes(payload)
        if sha(path) != digest:
            raise ValueError(f'Input checksum mismatch: {path}')
    release = work / 'release'
    with tarfile.open(work / 'spt_act_likelihood-1.0.tar.gz') as tar:
        for member in tar.getmembers():
            if not member.isfile():
                continue
            relative = Path(*Path(member.name).parts[1:])
            dest = release / relative
            if not dest.resolve().is_relative_to(release.resolve()):
                raise ValueError('Unsafe archive path')
            payload = tar.extractfile(member).read()
            if dest.exists() and dest.read_bytes() != payload:
                raise ValueError(f'Existing release file changed: {dest}')
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(payload)
    return release


def audit(work, act_data):
    release = acquire(work)
    sys.path.insert(0, str(release))
    import act_dr6_spt_lenslike as native
    directory = release / 'act_dr6_spt_lenslike/data/v1.2'
    response = act_data / 'like_corrs'
    if not response.is_dir():
        raise FileNotFoundError('Acquire the official ACT DR6 v1.2 response assets first')
    link = directory / 'like_corrs'
    if link.exists() and link.resolve() != response.resolve():
        raise ValueError('Response directory does not match requested input')
    if not link.exists():
        link.symlink_to(response.resolve(), target_is_directory=True)
    ell, tt, ee, bb, te = np.loadtxt(response / 'cosmo2017_10K_acc3_lensedCls.dat', unpack=True)
    ellp, pp = np.loadtxt(response / 'cosmo2017_10K_acc3_lenspotentialCls.dat', usecols=[0, 5], unpack=True)
    factor = 2 * np.pi / ell / (ell + 1)
    args = (ellp, pp * 2 * np.pi / 4, ell, tt * factor, ee * factor, te * factor, bb * factor)
    expected = {('baseline', True): 38.17, ('baseline', False): 38.67,
                ('extended', True): 41.27, ('extended', False): 41.73}
    cases = []
    for suffix in ['baseline', 'extended']:
        for lens_only in [True, False]:
            variant = 'actplanckspt3g_' + suffix
            with warnings.catch_warnings(record=True) as caught:
                data = native.load_data(variant, ddir=str(directory), lens_only=lens_only,
                                        like_corrections=not lens_only)
            ll, prediction = native.generic_lnlike(data, *args, return_theory=True)
            cov = data['cov']
            nact = data['binmat_act'].shape[0]
            corr = cov / np.sqrt(np.outer(np.diag(cov), np.diag(cov)))
            hartlap = (400 - len(cov) - 2) / 399
            precision = np.linalg.solve(cov, np.eye(len(cov))) * hartlap
            residual = data['data_binned_clkk'] - prediction
            independently_evaluated = float(residual @ precision @ residual)
            covariance_names = ['act_planck', 'act_spt', 'planck_spt']
            blocks = [(slice(0, nact), slice(nact, nact + 9)),
                      (slice(0, nact), slice(nact + 9, None)),
                      (slice(nact, nact + 9), slice(nact + 9, None))]
            independent_cov = np.zeros_like(cov)
            for s in [slice(0, nact), slice(nact, nact + 9), slice(nact + 9, None)]:
                independent_cov[s, s] = cov[s, s]
            independent_chi2 = float(residual @ np.linalg.solve(independent_cov, residual) * hartlap)
            npz = np.load(directory / 'muse_likelihood.npz')
            assert np.array_equal(cov[-16:, -16:], npz['cov_kk'])
            assert np.allclose(-2 * ll, independently_evaluated, rtol=0, atol=1e-12)
            assert abs(-2 * ll - expected[(suffix, lens_only)]) < .05
            symmetric_corr = (corr + corr.T) / 2
            asymmetry = float(abs(corr - corr.T).max())
            assert asymmetry < 1e-12
            assert np.linalg.eigvalsh(symmetric_corr).min() > 0
            supports = {}
            for label in ['act', 'planck', 'spt']:
                matrix = data['binmat_' + label]
                support = np.flatnonzero(np.any(matrix != 0, axis=0))
                supports[label] = [int(support.min()), int(support.max())]
            cases.append(dict(variant=variant, lens_only=lens_only,
                block_order=['ACT', 'Planck', 'SPT_MUSE'], block_lengths=[nact, 9, 16],
                chi2=float(-2 * ll), author_test_chi2=expected[(suffix, lens_only)],
                independent_quadratic_agreement=float(abs(-2 * ll - independently_evaluated)),
                hypothetical_block_diagonal_chi2_same_hartlap=independent_chi2,
                hartlap_precision_factor=hartlap,
                min_correlation_eigenvalue=float(np.linalg.eigvalsh(symmetric_corr).min()),
                covariance_max_asymmetry_in_correlation_units=asymmetry,
                covariance_symmetric_to_normalized_tolerance=bool(asymmetry < 1e-12),
                max_abs_cross_correlation={k: float(abs(corr[b]).max()) for k, b in zip(covariance_names, blocks)},
                bin_window_nonzero_multipoles=supports,
                spt_covariance_exactly_npz=True, warnings=[str(w.message) for w in caught]))
            del data
            gc.collect()
    # Data-only execution isolates a public API pitfall; this is not a science fit.
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        full0 = native.load_data('actplanckspt3g_baseline', ddir=str(directory), lens_only=False, like_corrections=False, indep=False)
        full1 = native.load_data('actplanckspt3g_baseline', ddir=str(directory), lens_only=False, like_corrections=False, indep=True)
    obj = object.__new__(native.ACTDR6LensLike)
    requested = obj.get_requirements()
    muse_details = {}
    import h5py
    with zipfile.ZipFile(work / 'muse_3g_like_march_2025.zip') as archive:
        name = 'muse_3g_like_march_2025/src/muse3glike/dat/90_150_220.h5'
        with h5py.File(io.BytesIO(archive.read(name)), 'r') as h5:
            muse_details = dict(original_data_dimensions=int(h5['data'].shape[0]),
                pp_bandpowers_match_original=bool(np.array_equal(npz['d_pp'], h5['data'][:16])),
                band_windows_match_original=bool(np.array_equal(npz['bpwf'], h5['BPWF']['ϕϕ'][:].T)),
                original_raw_covariance_equals_joint_gaussian_covariance=bool(np.array_equal(npz['cov_pp'], h5['covariance'][:16, :16])),
                joint_to_original_raw_covariance_diagonal_ratio=(np.diag(npz['cov_pp']) / np.diag(h5['covariance'][:16, :16])).tolist(),
                original_axis=h5['axis'].asstr()[()],
                original_release_scope='Transformed-space Gaussian including delensed EE and systematics; not interchangeable with Qu joint Gaussian lensing compression')
    inputs = {str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p): sha(p)
              for p in sorted(release.rglob('*')) if p.is_file() and '__pycache__' not in p.parts}
    # All available response inputs are pinned, including additional unused fiducial metadata.
    for p in sorted(response.iterdir()):
        if p.is_file():
            inputs[str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p)] = sha(p)
    return dict(status='released_likelihood_evaluation_reproduced_no_inference', native_camb_calls=0,
        scientific_target_changed=False, official_downloads={k: dict(url=v[0], sha256=v[1]) for k, v in RELEASES.items()},
        source_commit_crosscheck='cfd88b8ebe24e969d407637fd2ac20a31220ac29',
        source_sha256={str(Path(__file__).relative_to(ROOT)): sha(__file__)}, input_sha256=inputs,
        cases=cases, muse_original_release_crosscheck=muse_details,
        interface_pitfalls=dict(full_mode_indep_flag_leaves_covariance_unchanged=bool(np.array_equal(full0['cov'], full1['cov'])),
            default_cobaya_requirements=requested,
            required_explicit_spectra=['pp', 'tt', 'ee', 'te', 'bb'],
            generic_interface_validated=True,
            message='Ensure all response spectra are provided explicitly; do not infer independent full-mode covariance from indep=True'),
        inference_boundary='Released joint Gaussian lensing covariance includes analytic SPT cross blocks; no primary-lensing, lensing-BAO or lensing-SN cross block is supplied. No reconstruction-level covariance regeneration attempted.')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--work', type=Path, default=WORK)
    ap.add_argument('--act-data', type=Path, default=ROOT / '.work/unified-cosmology/external-probes/packages/data/ACT_dr6_likelihood/v1.2')
    ap.add_argument('--output', type=Path, default=ROOT / 'studies/unified_cosmology/results/inference/joint-lensing-availability.json')
    a = ap.parse_args()
    result = audit(a.work.resolve(), a.act_data.resolve())
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({'status': result['status'], 'cases': len(result['cases']), 'native_camb_calls': 0}))


if __name__ == '__main__':
    main()
