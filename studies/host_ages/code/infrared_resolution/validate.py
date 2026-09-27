#!/usr/bin/env python3
"""Independent linear-algebra checks and predeclared beam-model sensitivities.

These validate a conditional information calculation, not a dust luminosity
measurement. The empirical radial profile is the shared physical input.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
from scipy.linalg import qr, svd, lstsq
from scipy.ndimage import map_coordinates, gaussian_filter1d
import analyze as reference

ROOT = reference.ROOT
WORK = reference.WORK
OUT = reference.OUT
INPUT = reference.INPUT


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def quantiles(values):
    finite = [v for v in values if v is not None and np.isfinite(v)]
    summary = dict(zip(['p05', 'median', 'p95'], map(float, np.quantile(finite, [.05, .5, .95])))) if finite else dict(p05=None, median=None, p95=None)
    return dict(summary, identified=len(finite), numerically_unresolved=len(values)-len(finite))


def templates(host, cat, band, profile, neighbour_radius=2., fit_radius=4., pixel_average=False):
    """Reconstruct positions from the full source catalogue, without the KD tree."""
    coord = SkyCoord(float(host.deep_RA)*u.deg, float(host.deep_DEC)*u.deg)
    allcoord = SkyCoord(cat.RA.to_numpy()*u.deg, cat.DEC.to_numpy()*u.deg)
    sep = coord.separation(allcoord).arcsec
    fw = reference.FWHM[band]
    use = (sep <= neighbour_radius*fw) & (cat.ID.to_numpy() != int(host.deep_ID))
    idx = np.flatnonzero(use)
    idx = idx[np.argsort(sep[idx])]
    prefix = {'C3': 'CDFS-SWIRE-NEST', 'X3': 'XMM-LSS-NEST'}[host.field]
    with fits.open(INPUT/'help'/f'{prefix}_image_{band}_SMAP_v6.0.fits') as hdus:
        wcs = WCS(hdus['image'].header)
        scale = abs(wcs.pixel_scale_matrix[0, 0])*3600
        hx, hy = wcs.world_to_pixel(coord)
        nx, ny = wcs.world_to_pixel(allcoord[idx])
        radius = int(np.ceil(fit_radius*fw/scale))
        y, x = np.meshgrid(np.arange(round(float(hy))-radius, round(float(hy))+radius+1),
                          np.arange(round(float(hx))-radius, round(float(hx))+radius+1), indexing='ij')
        inside = (x >= 0) & (x < hdus['image'].data.shape[1]) & (y >= 0) & (y < hdus['image'].data.shape[0])
        x, y = x[inside], y[inside]
        good = ((hdus['flag'].data[y, x] == 0) & np.isfinite(hdus['image'].data[y, x]) &
                np.isfinite(hdus['error'].data[y, x]) & (hdus['error'].data[y, x] > 0))
        x, y = (x[good]-hx)*scale, (y[good]-hy)*scale
        circle = x*x+y*y <= (fit_radius*fw)**2
        x, y = x[circle], y[circle]
    dx, dy = np.r_[0., (nx-hx)*scale], np.r_[0., (ny-hy)*scale]
    a = np.zeros((len(x), len(dx)))
    offsets = [-1/3, 0., 1/3] if pixel_average else [0.]
    for ox in offsets:
        for oy in offsets:
            a += np.interp(np.hypot(x[:, None]+ox*scale-dx, y[:, None]+oy*scale-dy),
                           profile['radius'], profile['profile'], right=0.) / len(offsets)**2
    sky = np.column_stack([np.ones(len(x)), x/fw, y/fw])
    # Independent least-squares elimination, rather than the primary QR projector.
    a -= sky @ lstsq(sky, a, lapack_driver='gelsy')[0]
    return a, sep[idx], sky


def solve(a, tolerance=1e-10):
    nuisance = a[:, 1:]
    u, s, _ = svd(nuisance, full_matrices=False, lapack_driver='gesvd')
    rank = int(np.sum(s > tolerance*s[0])) if len(s) else 0
    h = a[:, 0]
    residual = h-u[:, :rank]@(u[:, :rank].T@h)
    relative = np.linalg.norm(residual)/np.linalg.norm(h)
    return (float(1/relative) if relative > 1e-10 else np.nan), residual, rank


def main():
    result = json.loads((OUT/'information.json').read_text())
    for name, expected in result['dependencies_sha256'].items():
        assert digest(ROOT/name) == expected, name
    frame = pd.read_csv(ROOT/result['output'], dtype={'sn_id': str})
    assert digest(ROOT/result['output']) == result['output_sha256']
    assert len(frame) == 795 and frame.sn_id.nunique() == 265
    assert not frame.duplicated(['sn_id', 'band_um']).any()
    table = pd.read_csv(INPUT/'des-deep-crosswalk.csv', dtype={'SNID': str, 'deep_ID': 'Int64'})
    table = table[table.primary_match & table.selected_Dovekie & table.deep_quality & table.host_dlr_lt4]
    maps = pd.read_csv(INPUT/'des-spire-photometry.csv', dtype={'sn_id': str})
    maps = maps[maps.eightband_331 & maps.valid_map_measurement & (maps.common_sky_offsets >= 64)]
    counts = maps.groupby('sn_id').band_um.nunique()
    assert set(frame.sn_id) == set(counts[counts == 3].index)
    # Selection is independent of the computed inflation and of infrared flux.
    selected = sorted(set(frame.sn_id), key=lambda x: hashlib.sha256(x.encode()).hexdigest())[:12]
    beams, beam_records = reference.beams()
    extended_beams = {}
    for record, name in zip(beam_records, ['PSW', 'PMW', 'PLW']):
        band = record['band_um']
        with fits.open(WORK/f'0x5000241aL_{name}_bgmod10_1arcsec.fits.gz') as hdus:
            cx, cy = record['fitted_center_xy']
            radius = np.arange(0., 5*reference.FWHM[band]+.25, .25)
            theta = np.arange(720)*2*np.pi/720
            values = map_coordinates(hdus['image'].data,
                                     [cy+radius[:, None]*np.sin(theta), cx+radius[:, None]*np.cos(theta)], order=1)
            profile = values.mean(axis=1)
            profile /= profile[0]
            extended_beams[band] = {'radius': radius, 'profile': profile}
            assert np.allclose(profile[:len(beams[band]['profile'])], beams[band]['profile'], atol=1e-12, rtol=0.)
    cats = {}
    for field in ['C3', 'X3']:
        cat = pq.read_table(INPUT/'des'/f'Y3_DEEP_FIELDS_PHOTOM-SN-{field}-0000.parquet',
                            columns=['ID', 'RA', 'DEC', 'FLAGS', 'MASK_FLAGS', 'KNN_CLASS']).to_pandas()
        cats[field] = cat[(cat.FLAGS == 0) & (cat.MASK_FLAGS == 0) & (cat.KNN_CLASS == 1)]
    checks = []
    rng = np.random.default_rng(927042)
    for sn in selected:
        host = table[table.SNID == sn].iloc[0]
        for band in reference.BANDS:
            a, sep, sky = templates(host, cats[host.field], band, beams[band])
            v, residual, rank = solve(a)
            primary = frame[(frame.sn_id == sn) & (frame.band_um == band)].iloc[0]
            error = abs(v/primary.all_candidates_noise_inflation-1)
            assert error < 2e-7
            assert len(sep) == primary.neighbours_within_2FWHM
            assert abs(sep[0]-primary.nearest_separation_arcsec) < 1e-8
            q, r = qr(a[:, 1:], mode='economic', pivoting=False)
            qr_residual = a[:, 0]-q@(q.T@a[:, 0])
            qv = np.linalg.norm(a[:, 0])/np.linalg.norm(qr_residual)
            qr_error = abs(qv/v-1)
            assert qr_error < 2e-7
            # Original unprojected constraints: sky and every neighbour must be
            # orthogonal to the host residual, and its dot with h equals r.r.
            norm = np.linalg.norm(a[:, 0])
            orthogonality = float(np.max(np.abs(a[:, 1:].T@residual)) /
                                  (np.linalg.norm(a[:, 1:], axis=0).max()*np.linalg.norm(residual)))
            assert orthogonality < 2e-6
            influence = residual/(residual@a[:, 0])
            se = float(np.linalg.norm(influence))
            assert abs(se/(v/norm)-1) < 2e-6
            # Noise injections sample the exact linear estimator. They are a
            # calculation check under assumed white noise, not coverage on sky.
            noise = rng.normal(size=(len(a), 2000))
            estimates = influence@noise
            zmean = float(estimates.mean()/(se/np.sqrt(len(estimates))))
            ratio = float(estimates.std(ddof=1)/se)
            assert abs(zmean) < 5 and .91 < ratio < 1.09
            variations = {}
            for tag, nr, fr in [('neighbours_1p5', 1.5, 4.), ('neighbours_3', 3., 4.),
                                ('fit_radius_3', 2., 3.), ('fit_radius_5', 2., 5.)]:
                alt, _, _ = templates(host, cats[host.field], band, beams[band], nr, fr)
                variations[tag] = solve(alt)[0]/v
            variations['rank_tolerance_1e8'] = solve(a, 1e-8)[0]/v
            variations['rank_tolerance_1e12'] = solve(a, 1e-12)[0]/v
            variations['rank_tolerance_1e6'] = solve(a, 1e-6)[0]/v
            for tag, profile, average in [('pixel_average_3x3', beams[band], True),
                                           ('beam_support_5', extended_beams[band], False)]:
                alt, _, _ = templates(host, cats[host.field], band, profile, pixel_average=average)
                variations[tag] = solve(alt)[0]/v
            shortened = {key: value[beams[band]['radius'] <= 3*reference.FWHM[band]]
                         for key, value in beams[band].items() if key in ['radius', 'profile']}
            alt, _, _ = templates(host, cats[host.field], band, shortened)
            variations['beam_support_3'] = solve(alt)[0]/v
            smoothed = {'radius': beams[band]['radius'],
                        'profile': gaussian_filter1d(beams[band]['profile'], 1./.25, mode='reflect')}
            alt, _, _ = templates(host, cats[host.field], band, smoothed)
            variations['radial_smoothing_1arcsec'] = solve(alt)[0]/v
            # A host template lying in the nuisance span has no identifiable
            # amplitude. Do not report inverse floating-point residuals as a
            # huge but measured finite uncertainty.
            variations = {key: (float(value) if np.isfinite(value) else None)
                          for key, value in variations.items()}
            checks.append(dict(sn_id=sn, band_um=band, primary_relative_error=error,
                               QR_relative_error=qr_error, neighbour_rank=rank,
                               original_unit_orthogonality=orthogonality,
                               noise_mean_standard_errors=zmean, noise_sd_to_expected=ratio,
                               sensitivity_ratio=variations))
    # Independent continuous Gaussian overlap from numerical 2-D integration.
    axis = np.arange(-80., 80.01, .2)
    y, x = np.meshgrid(axis, axis, indexing='ij')
    analytic_checks = []
    for separation in [2., 5., 10., 20.]:
        sigma = 8.
        h = np.exp(-(x*x+y*y)/(2*sigma*sigma)).ravel()
        n = np.exp(-((x-separation)**2+y*y)/(2*sigma*sigma)).ravel()
        rho = float(h@n/np.linalg.norm(h)/np.linalg.norm(n))
        exact = np.exp(-separation**2/(4*sigma*sigma))
        assert abs(rho-exact) < 1e-12
        fwhm = separation*np.sqrt(4*np.log(2)/(-np.log(.75)))
        assert abs(1/np.sqrt(1-np.exp(-4*np.log(2)*(separation/fwhm)**2))-2) < 1e-12
        analytic_checks.append(dict(separation=separation, overlap_error=abs(rho-exact)))
    sensitivities = {}
    for band in reference.BANDS:
        group = [c for c in checks if c['band_um'] == band]
        sensitivities[str(band)] = {key: quantiles([c['sensitivity_ratio'][key] for c in group])
                                    for key in group[0]['sensitivity_ratio']}
    output = dict(status='passed', completed_utc=datetime.now(timezone.utc).isoformat(),
                  selected_hosts=selected, checked_templates=len(checks),
                  white_noise_injections_per_template=2000, checks=checks,
                  sensitivity_ratios=sensitivities, gaussian_checks=analytic_checks,
                  scope='Conditional linear algebra and model-boundary sensitivity. White-noise recovery does not establish a correct image covariance or astrophysical source model.',
                  dependencies_sha256={str(p.relative_to(ROOT)): digest(p) for p in
                    [OUT/'information.json', Path(__file__), Path(reference.__file__), Path(__file__).with_name('validation-design.json')]})
    (OUT/'validation.json').write_text(json.dumps(output, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k: output[k] for k in ['status', 'checked_templates', 'sensitivity_ratios']}, indent=2))


if __name__ == '__main__':
    main()
