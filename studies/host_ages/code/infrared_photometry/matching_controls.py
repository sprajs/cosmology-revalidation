#!/usr/bin/env python3
"""Independent spherical matches and local shifted-position controls.

The controls measure local catalogue coincidence rates, not posterior host
association probabilities or infrared image deblending. No flux is fitted.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u

ROOT = Path(__file__).resolve().parents[4]
WORK = ROOT/'.work/infrared-photometry'
OUT = ROOT/'studies/host_ages/results/infrared_photometry'
HERE = Path(__file__).parent
CATALOGUES = {'chandra_cat_f05': 'C3', 'chandra_24_cat_f05': 'C3',
              'servscdfsi12': 'C3', 'xmm_cat_s05': 'X3',
              'xmm_24_cat_s05': 'X3', 'servsxmmi12': 'X3'}


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def separation(ra, dec, source_ra, source_dec):
    r1, d1, r2, d2 = np.deg2rad(ra), np.deg2rad(dec), np.deg2rad(source_ra), np.deg2rad(source_dec)
    hav = np.sin((d2-d1)/2)**2+np.cos(d1)*np.cos(d2)*np.sin((r2-r1)/2)**2
    return np.rad2deg(2*np.arctan2(np.sqrt(np.clip(hav, 0, 1)), np.sqrt(np.clip(1-hav, 0, 1))))*3600


def shifted(ra, dec, offset_arcsec, bearing_degrees):
    r, d = np.deg2rad(ra), np.deg2rad(dec)
    delta = np.deg2rad(offset_arcsec/3600)
    bearing = np.deg2rad(bearing_degrees)
    d2 = np.arcsin(np.sin(d)*np.cos(delta)+np.cos(d)*np.sin(delta)*np.cos(bearing))
    r2 = r+np.arctan2(np.sin(bearing)*np.sin(delta)*np.cos(d), np.cos(delta)-np.sin(d)*np.sin(d2))
    return np.rad2deg(r2) % 360, np.rad2deg(d2)


def main():
    design = json.loads((HERE/'matching-control-design.json').read_text())
    hosts = pd.read_csv(WORK/'hosts.csv', dtype={'SNID': str, 'deep_ID': str})
    candidates = pd.read_csv(WORK/'candidates.csv', dtype={'source_key': str, 'cntr': str})
    assert len(hosts) == 265 and hosts.SNID.is_unique and hosts.deep_ID.is_unique
    assert not candidates.duplicated(['catalogue', 'source_key']).any()
    assert set(candidates.catalogue) == set(CATALOGUES)
    radii = np.asarray(design['radii_arcsec'])
    controls = [(r, b) for r in design['position_controls']['radii_arcsec']
                for b in design['position_controls']['bearings_degrees']]
    assert max(r for r, _ in controls)+max(radii) < design['position_controls']['candidate_query_radius_arcsec']
    rows, summaries, comparisons = [], [], []
    rng = np.random.default_rng(927061)
    for catalogue, field in CATALOGUES.items():
        data = candidates[candidates.catalogue == catalogue]
        assert set(data.field) == {field}
        parent = hosts[hosts.field == field]
        actual = np.zeros((len(parent), len(radii)), dtype=int)
        control = np.zeros((len(parent), len(controls), len(radii)), dtype=int)
        for i, host in enumerate(parent.itertuples(index=False)):
            distance = separation(host.deep_RA, host.deep_DEC, data.ra.to_numpy(), data.dec.to_numpy())
            local = distance <= design['position_controls']['candidate_query_radius_arcsec']+1e-6
            source = data[local]
            distance = distance[local]
            actual[i] = (distance[:, None] <= radii).sum(axis=0)
            nearest = int(np.argmin(distance)) if len(distance) else None
            for j, (offset, bearing) in enumerate(controls):
                ra, dec = shifted(host.deep_RA, host.deep_DEC, offset, bearing)
                assert abs(float(separation(host.deep_RA, host.deep_DEC, ra, dec))-offset) < 1e-7
                delta = separation(ra, dec, source.ra.to_numpy(), source.dec.to_numpy())
                control[i, j] = (delta[:, None] <= radii).sum(axis=0)
            # Compare an independent library implementation for every primary
            # nearest source, retaining all empty candidate neighborhoods.
            if nearest is not None:
                match = source.iloc[nearest]
                a = SkyCoord(host.deep_RA*u.deg, host.deep_DEC*u.deg)
                b = SkyCoord(match.ra*u.deg, match.dec*u.deg)
                error = abs(a.separation(b).arcsec-distance[nearest])
                assert error < 1e-7
                comparisons.append(float(error))
            for k, radius in enumerate(radii):
                rows.append(dict(SNID=host.SNID, deep_ID=host.deep_ID, field=field,
                                 catalogue=catalogue, radius_arcsec=float(radius),
                                 candidates_in_radius=int(actual[i, k]),
                                 nearest_source_key=source.iloc[nearest].source_key if nearest is not None else None,
                                 nearest_separation_arcsec=float(distance[nearest]) if nearest is not None else None,
                                 shifted_positions=len(controls),
                                 shifted_with_match=int(np.sum(control[i, :, k] > 0)),
                                 shifted_with_multiple=int(np.sum(control[i, :, k] > 1)),
                                 candidates_within75=len(source)))
        resamples = rng.integers(0, len(parent), size=(2000, len(parent)))
        for k, radius in enumerate(radii):
            observed = (actual[:, k] > 0).astype(float)
            shifted_fraction = (control[:, :, k] > 0).mean(axis=1)
            excess = observed-shifted_fraction
            interval = np.quantile(excess[resamples].mean(axis=1), [.025, .975])
            summaries.append(dict(catalogue=catalogue, field=field, hosts=len(parent), radius_arcsec=float(radius),
                                  actual_with_match=int(observed.sum()),
                                  actual_unique=int(np.sum(actual[:, k] == 1)),
                                  actual_multiple=int(np.sum(actual[:, k] > 1)),
                                  shifted_match_fraction=float(shifted_fraction.mean()),
                                  actual_minus_shifted_fraction=float(excess.mean()),
                                  conditional_host_bootstrap_95=list(map(float, interval))))
    frame = pd.DataFrame(rows)
    output = WORK/'matching-controls.csv'
    frame.to_csv(output, index=False)
    dependencies = [WORK/'hosts.csv', WORK/'candidates.csv', HERE/'matching-control-design.json', Path(__file__), output]
    record = dict(status='passed', created_utc=datetime.now(timezone.utc).isoformat(),
                  host_catalogue_configurations=len(rows)//len(radii), rows=len(rows),
                  shifted_positions_per_configuration=len(controls), radii_arcsec=list(radii),
                  independent_nearest_angle_checks=len(comparisons),
                  max_spherical_separation_difference_arcsec=max(comparisons),
                  summaries=summaries,
                  scope='Raw catalogue positional coincidence and local shifted-position controls. No flux-quality filtering, calibrated association probability or deblending. Host resampling conditions on the source list and offsets and omits shared-source or field correlation. Unknown footprint and correlated galaxy environments remain.',
                  dependencies_sha256={str(p.relative_to(ROOT)): sha(p) for p in dependencies})
    (OUT/'matching-controls.json').write_text(json.dumps(record, indent=2, allow_nan=False)+'\n')
    print(json.dumps([row for row in summaries if row['radius_arcsec'] == 2.], indent=2))


if __name__ == '__main__':
    main()
