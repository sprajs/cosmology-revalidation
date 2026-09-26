"""Verify acquired author simulation assets; no simulation or outcome fitting."""
from pathlib import Path
import hashlib
import json
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'runs/research_2026_09_26/raisin_simulation_assets'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def active_rows(path, key):
    return [line.split()[1:] for line in path.read_text().splitlines()
            if line.strip().startswith(key + ':')]


def efficiency(path):
    result = {}
    band = None
    for line in path.read_text().splitlines():
        bits = line.split()
        if not bits:
            continue
        if bits[0] == 'FILTER:':
            band = bits[1]
            result[band] = []
        if bits[0] == 'SNR:':
            result[band].append([float(x) for x in bits[1:3]])
    return result


def main():
    acquisition = json.loads((OUT / 'acquisition.json').read_text())
    verified = []
    for row in acquisition['files']:
        path = ROOT / row['path']
        data = path.read_bytes()
        assert sha(path) == row['sha256']
        assert hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest() == row['git_blob']
        verified.append(row['source_path'])
    author = OUT / 'author'
    des = author / 'sim/inputs/DES'
    oir = author / 'OIR.J19/OIR.INFO'
    bundled = ROOT / 'phase2/official/inputs/SNDATA_ROOT/models/OIR/OIR.J19/OIR.INFO'
    sigma = np.array(active_rows(oir, 'COLOR_SIGMA'), dtype=float).ravel()
    corr = np.array(active_rows(oir, 'COLOR_CORMAT'), dtype=float)
    scale = float(active_rows(oir, 'COLOR_SIGMA_SCALE')[0][0])
    assert sigma.shape == (7,) and corr.shape == (7, 7)
    node_cov = 1.3 * (corr * np.outer(scale * sigma, scale * sigma) + 1e-9 * np.diag(np.diag(corr)))
    pipeline = efficiency(des / 'SEARCHEFF_PIPELINE_RAISIN.DAT')
    comparison = efficiency(des / 'SEARCHEFF_PIPELINE_DES.DAT')
    shared = {}
    for band in 'griz':
        ra = np.array(pipeline[band]); de = np.array(comparison[band])
        assert np.array_equal(ra[:-1], de)
        assert np.array_equal(ra[-1], [1e6, 1.])
        shared[band] = {'shared_rows': len(de), 'additional_row': ra[-1].tolist()}
    nominal_spec = np.array(active_rows(des / 'SEARCHEFF_SPEC_DES_Moller_G10_v7.DAT', 'SPECEFF'), float)
    alternative_spec = np.array(active_rows(author / 'sim/flatdist/DES/SEARCHEFF_SPEC_DES_Moller_G10_v7.DAT', 'SPECEFF'), float)
    common_spec = {float(x): float(y) for x, y in alternative_spec}
    matches = [(x, y, common_spec[x]) for x, y in nominal_spec if x in common_spec]
    spec_delta = max(abs(y - other) for x, y, other in matches)
    error_rows = np.array(active_rows(des / 'DES3YR_SIM_ERRORFUDGES.DAT', 'ROW'), float)
    source = ROOT / 'sources/repos/RickKessler__SNANA@v11_04k/src/sntools_genSmear.c'
    logic = ROOT / 'phase2/official/inputs/SNDATA_ROOT/models/searcheff/SEARCHEFF_PIPELINE_LOGIC.DAT'
    result = {
        'scope': 'Source/configuration availability and numerical content only; historical execution remains unlinked.',
        'verified_files': verified,
        'OIR_author_identical_to_bundled_2024': oir.read_bytes() == bundled.read_bytes(),
        'OIR_sigma': sigma.tolist(), 'OIR_correlation_eigenvalues': np.linalg.eigvalsh(corr).tolist(),
        'OIR_native_node_covariance_eigenvalues': np.linalg.eigvalsh(node_cov).tolist(),
        'OIR_node_covariance_rule': 'v11_04k source applies COLOR_SIGMA_SCALE, sigma outer product times correlation, diagonal 1e-9, then covariance factor 1.3. These are spectral-node values before wavelength/filter integration.',
        'pipeline_shared_content': shared,
        'spectroscopic_nominal_rows': len(nominal_spec),
        'spectroscopic_flatdist_rows': len(alternative_spec),
        'spectroscopic_shared_magnitude_rows': len(matches),
        'spectroscopic_max_shared_probability_difference': spec_delta,
        'error_map_rows': len(error_rows),
        'error_scale_tabulated_range': [float(error_rows[:, 1].min()), float(error_rows[:, 1].max())],
        'DES_bundled_logic': [line for line in logic.read_text().splitlines() if line.startswith('DES:')],
        'not_established': [
            'Historical effective input/command manifest and original SNDATA_ROOT version',
            'Joint distribution of cadence, extraction error, flags, population, HST follow-up and event retention',
            'Original full-DIFFIMG exposure metadata for omitted rows',
            'Regenerated distance-bias correction and cosmology'
        ],
        'inputs_sha256': {str(p.relative_to(ROOT)): sha(p) for p in [
            OUT / 'acquisition.json', oir, bundled, source, logic, Path(__file__)
        ]},
    }
    (OUT / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
    (OUT / 'executed_source.py').write_bytes(Path(__file__).read_bytes())
    files = {str(p.relative_to(OUT)): sha(p) for p in OUT.rglob('*') if p.is_file() and p.name != 'manifest.json'}
    (OUT / 'manifest.json').write_text(json.dumps({'files_sha256': files}, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
