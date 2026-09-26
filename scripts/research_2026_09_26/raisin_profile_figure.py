"""Plot conditional shape profiles; displayed Q bands are not confidence sets."""
from pathlib import Path
import csv
import hashlib
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
G = ROOT / 'runs/research_2026_09_26/astra_design/raisin_signed_refit/fixed-c-profile/global-profile'
OUT = ROOT / 'runs/research_2026_09_26/raisin_profile_decomposition'
summary = json.loads((G / 'result.json').read_text())
rows = list(csv.DictReader((G / 'shape-profiles.csv').open()))
fig, axes = plt.subplots(2, 1, sharex=True, figsize=(8.3, 7.0), constrained_layout=True)
for arm, label, color in [('A', '65 positive epochs', '#0069a0'), ('B', '74 signed epochs', '#b44820')]:
    r = sorted([r for r in rows if r['metric'] == 'Banchor_' + arm and r['stage'] == 'fine'], key=lambda r: float(r['shape']))
    x = np.array([float(z['shape']) for z in r]); q = np.array([float(z['Q']) for z in r]); d = np.array([float(z['DLMAG']) for z in r])
    best = summary['metrics']['Banchor_' + arm]['finest']
    axes[0].plot(x, q - best['Q'], color=color, label=label, lw=1.7)
    axes[1].plot(x, d, color=color, lw=1.7)
    axes[0].scatter(best['shape'], 0, color=color, s=45, zorder=5)
    axes[1].scatter(best['shape'], best['DLMAG'], color=color, s=45, zorder=5)
    axes[1].annotate(f"{best['DLMAG']:.3f} mag", (best['shape'], best['DLMAG']),
                     xytext=(7, -22), textcoords='offset points', color=color)
axes[0].set_ylim(-0.15, 5); axes[0].set_ylabel('Q minus each arm’s minimum')
axes[0].axhline(1, color='#aaaaaa', ls='--', lw=.8)
axes[0].legend(frameon=False, loc='upper left')
axes[0].set_title('DES16C1cim: the preferred template changes,\nwhile competing distance solutions fit almost equally well', loc='left', fontsize=13)
axes[1].set_ylabel('Profiled distance parameter (mag)'); axes[1].set_xlabel('SNooPy stretch parameter')
axes[1].set_xlim(.7, 1.3)
for ax in axes:
    ax.grid(alpha=.18); ax.spines[['top', 'right']].set_visible(False)
fig.get_layout_engine().set(rect=(0, .08, 1, .92))
fig.text(.08, .025, 'Fixed header peak and common covariance; extinction and amplitude profiled.\nNumerical profile geometry only: the Q=1 guide is not a calibrated confidence interval.', fontsize=9, color='#444444')
for ext in ['png', 'pdf']:
    fig.savefig(OUT / ('conditional-shape-profile.' + ext), dpi=170)
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
(OUT / 'figure-manifest.json').write_text(json.dumps({'files_sha256': {
    str(p.relative_to(ROOT)): sha(p) for p in [Path(__file__), G / 'result.json', G / 'shape-profiles.csv',
    OUT / 'conditional-shape-profile.png', OUT / 'conditional-shape-profile.pdf']}}, indent=2) + '\n')
