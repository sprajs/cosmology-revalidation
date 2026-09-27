#!/usr/bin/env python3
"""Plot empirical beam widths and conditional nearest-neighbour information."""
from datetime import datetime, timezone
import json
import hashlib
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.ndimage import map_coordinates

ROOT = Path(__file__).resolve().parents[4]
OUT = ROOT/'studies/host_ages/results/infrared_resolution'
WORK = ROOT/'.work/infrared-resolution'
BANDS = [250, 350, 500]


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    source = OUT/'information.json'
    result = json.loads(source.read_text())
    table_path = ROOT/result['output']
    assert sha(table_path) == result['output_sha256']
    table = pd.read_csv(table_path)
    profiles = {}
    inputs = []
    for record, name in zip(result['beam_records'], ['PSW', 'PMW', 'PLW']):
        path = WORK/f'0x5000241aL_{name}_bgmod10_1arcsec.fits.gz'
        assert sha(path) == record['source_sha256']
        inputs.append(path)
        with fits.open(path) as hdus:
            cx, cy = record['fitted_center_xy']
            radius = np.arange(0., 4*record['nominal_FWHM_arcsec']+.25, .25)
            angle = np.arange(720)*2*np.pi/720
            values = map_coordinates(hdus['image'].data,
                                     [cy+radius[:, None]*np.sin(angle), cx+radius[:, None]*np.cos(angle)], order=1)
            profile = values.mean(axis=1)
            profiles[record['band_um']] = {'radius': radius, 'profile': profile/profile[0]}
    colors = ['#006d77', '#c67c00', '#8351a8']
    with plt.rc_context({'font.size': 11, 'axes.spines.top': False,
                         'axes.spines.right': False, 'savefig.facecolor': 'white'}):
        figure, axes = plt.subplots(1, 2, figsize=(11.6, 4.7))
        for band, color in zip(BANDS, colors):
            p = profiles[band]
            axes[0].plot(p['radius'], p['profile'], label=f'{band} μm', color=color, lw=2)
            values = np.sort(table.loc[table.band_um == band, 'nearest_pair_radial_noise_inflation'])
            assert len(values) == 265 and np.isfinite(values).all()
            axes[1].step(values, np.arange(1, len(values)+1)/len(values), where='post', color=color, lw=2,
                         label=f'{band} μm: median {np.median(values):.2f}')
        separation = result['summaries'][0]['nearest_separation_arcsec']['median']
        axes[0].axvline(separation, color='#444444', ls='--', lw=1)
        axes[0].text(22, .70, f'Median nearest galaxy\nseparation: {separation:.2f} arcsec', fontsize=10)
        axes[0].set(xlim=(0, 65), ylim=(-.03, 1.04), xlabel='Angular radius (arcsec)',
                    ylabel='Empirical beam / central value', title='A. Galaxies overlap within the measured beam')
        axes[0].legend(frameon=False, loc='upper right')
        axes[1].axvline(2, color='#777777', ls='--', lw=1)
        axes[1].set(xscale='log', xlim=(1, 35), ylim=(0, 1.02), xlabel='Host-amplitude standard-error inflation',
                    ylabel='Fraction of 265 hosts', title='B. Freeing just the nearest neighbour')
        axes[1].set_xticks([1, 2, 5, 10, 20], labels=['1', '2', '5', '10', '20'])
        axes[1].legend(frameon=False, loc='lower right')
        figure.text(.5, .025, 'Conditional geometry with white pixel noise; no host dust flux is inferred. More neighbours and PSF/covariance uncertainty remain.',
                    ha='center', fontsize=9)
        figure.tight_layout(rect=(0, .075, 1, 1))
        output = ROOT/'docs/figures/infrared-source-separation.png'
        figure.savefig(output, dpi=180)
        plt.close(figure)
    paths = [source, table_path, Path(__file__), output]+inputs
    record = dict(created_utc=datetime.now(timezone.utc).isoformat(), figure=str(output.relative_to(ROOT)),
                  dependencies_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths},
                  scope='Empirical calibration-source radial beams and nearest-pair conditional information; population percentiles are not confidence intervals.')
    (OUT/'figure.json').write_text(json.dumps(record, indent=2)+'\n')
    print(output)


if __name__ == '__main__':
    main()
