#!/usr/bin/env python3
"""Render the manuscript figures directly from verified result records."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'validation/figures'
OUT.mkdir(exist_ok=True)
BASE = 'revalidation-20260926-'


def read(path):
    return json.loads(path.read_text())


plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'savefig.facecolor': 'white', 'figure.facecolor': 'white'})
core = read(ROOT/'validation/reports/core.json')['cosmology']
values = [('Flat ΛCDM · supernovae', core['posterior_q0_mean'],
           [1.5*x-1 for x in core['posterior_95_interval']])]
for name, label in [('joint', 'CPL · supernovae + BAO'),
                    ('qbins', 'Flexible expansion · first bin (z < 0.1)'),
                    ('template', 'CPL + imposed age correction')]:
    q = read(ROOT/'results'/(BASE+name)/'summary.json')['q0']
    values.append((label, q['mean'], [q['q025_median_q975'][i] for i in [0,2]]))
fig, ax = plt.subplots(figsize=(10.5, 4.3), layout='constrained')
for y, (label, mean, ci) in enumerate(values):
    colour = '#ac4c21' if y == 3 else '#176b86'
    ax.plot(ci, [y,y], color=colour, lw=4, solid_capstyle='round')
    ax.scatter([mean], [y], color=colour, edgecolors='white', s=75, zorder=3)
    ax.text(.30, y, f'{mean:+.3f}  [{ci[0]:+.3f}, {ci[1]:+.3f}]', va='center', fontsize=10)
ax.axvline(0, color='#555555', ls='--', lw=1)
ax.set(yticks=range(4), yticklabels=[v[0] for v in values], xlim=(-.68,.87), ylim=(3.65,-.65),
       xlabel='Deceleration parameter: negative = acceleration', title='Different assumptions give different conditional inferences')
ax.set_xticks([-.6,-.4,-.2,0,.2])
ax.grid(axis='x', alpha=.15)
fig.savefig(OUT/'expansion.png', dpi=180)
plt.close(fig)

hst = read(ROOT/'validation/reports/hst.json')
fig, axes = plt.subplots(1,2,figsize=(11,3.9),layout='constrained',gridspec_kw={'width_ratios':[1,1.1]})
for ax, items, title, limits in [
    (axes[0], [(label,hst['science']['visits'][key]['tile_influence']) for key,label in [('search','Science: search'),('template','Science: template')]], 'Two science repeat pairs', (.4,1.08)),
    (axes[1], [('Dark: all eligible pixels',hst['dark_calibrated']['pixel_block_influence']),('Dark: fixed apertures',hst['dark_calibrated']['aperture_16_region_influence'])], 'One native dark pair', (0,4.6))]:
    for y,(label,item) in enumerate(items):
        ax.plot(item['deletion_ratio_range'],[y,y],lw=4,color='#82b5c4',solid_capstyle='round')
        ax.scatter([item['full_ratio']],[y],s=70,color='#176b86',zorder=3)
        ax.text(limits[0]+.03*(limits[1]-limits[0]),y-.24,f"{label}: {item['full_ratio']:.3f}",fontsize=10)
    ax.axvline(1,color='#555555',ls='--',lw=1)
    ax.set(xlim=limits,ylim=(1.5,-.65),yticks=[],xlabel='Squared difference / quoted variance',title=title)
    ax.grid(axis='x',alpha=.15)
fig.suptitle('Spatial influence depends on what the measurement averages',fontsize=13)
fig.supxlabel('Dots: full sample. Lines: fixed-group deletion ranges, not confidence intervals.',fontsize=10)
fig.savefig(OUT/'spatial-influence.png',dpi=180)
plt.close(fig)
