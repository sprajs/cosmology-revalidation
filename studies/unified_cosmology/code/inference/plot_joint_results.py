"""Plot qualified native joint measurements, without mixing physical targets."""
import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[4]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('measurement', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = json.loads(args.measurement.read_text())
    assert report['status'] == 'qualified_conditional_measurements'
    from measurement_summary import summarize_run
    # Revalidate the actual parents, including all numerical gates and files.
    rows = []
    for row in report['measurements']:
        summary_path = ROOT / row['native_correction_summary']
        summary = json.loads(summary_path.read_text())
        folder = (ROOT / summary['selection_path']).parent.parent
        actual = summarize_run(folder, summary_path)
        assert actual == row, 'Reported measurement no longer matches its parents.'
        rows.append(actual)
    assert rows
    assert all(row['settings']['sample'] == 'dovekie' and
               row['settings']['calibration'] == 'official_planck' for row in rows), \
        'This figure caption is defined for the Dovekie/official-calibration targets.'
    labels = []
    for row in rows:
        settings = row['settings']
        model = {'lcdm': 'Flat ΛCDM', 'cpl': 'Evolving dark energy'}[settings['model']]
        evolution = {'none': 'fixed SN luminosity', 'linear': 'free linear brightness drift',
                     'smooth01': 'smooth brightness, σ = 0.1 mag',
                     'smooth03': 'smooth brightness, σ = 0.3 mag'}[settings['evolution']]
        labels.append(model + '\n' + evolution)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, 4, figsize=(13, 2.6 + .8*len(rows)), sharey=True)
    colors = ['#245f9b', '#268073', '#b56c19', '#8b559a', '#be5360']
    keys = ['H0', 'omegam', 'q0', 'j0']
    titles = ['Expansion rate today', 'Matter fraction', 'Acceleration today', 'Change in acceleration']
    xlabels = ['$H_0$ [km s$^{-1}$ Mpc$^{-1}$]', '$\\Omega_m$', '$q_0$', '$j_0$']
    for axis, key, title, xlabel in zip(axes, keys, titles, xlabels):
        if key in ['q0', 'j0']:
            axis.axvline(0., color='#777777', lw=1, ls='--', zorder=0)
        for index, row in enumerate(rows):
            values = np.asarray(row['posterior'][key]['quantiles_025_16_50_84_975'])
            assert np.isfinite(values).all() and np.all(np.diff(values) >= 0)
            y = len(rows)-1-index
            color = colors[index % len(colors)]
            axis.plot(values[[0, 4]], [y, y], color=color, lw=1.5)
            axis.plot(values[[1, 3]], [y, y], color=color, lw=4)
            axis.scatter(values[2], y, s=28, color=color, zorder=3)
        axis.set_title(title, fontsize=11, pad=15)
        axis.set_xlabel(xlabel)
        axis.grid(axis='x', alpha=.15)
        axis.set_ylim(-.6, len(rows)-.4)
        axis.spines['left'].set_visible(False)
        axis.tick_params(axis='y', length=0)
    axes[0].set_yticks(np.arange(len(rows))[::-1], labels)
    fig.suptitle('A shared cosmology from supernovae, BAO and the CMB', fontsize=16, y=.985)
    fig.text(.5, .91, 'Dovekie · DESI DR2 · Planck + ACT + SPT', ha='center', color='#444444')
    fig.text(.5, .07, 'Dots: medians. Thick/thin lines: 68%/95% conditional credible intervals.\n'
             'q < 0: accelerating expansion. j < 0: scale-factor acceleration decreases with time.\n'
             'Brightness alternatives are declared priors, not measured progenitor-age corrections.',
             ha='center', fontsize=9, color='#444444')
    fig.subplots_adjust(left=.245, right=.985, bottom=.25, top=.80, wspace=.31)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180, facecolor='white')
    plt.close(fig)
    provenance = {'status': 'generated_from_qualified_native_measurements',
                  'input': str(args.measurement.resolve().relative_to(ROOT)),
                  'input_sha256': hashlib.sha256(args.measurement.read_bytes()).hexdigest(),
                  'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  'figure': str(args.output.resolve().relative_to(ROOT)),
                  'figure_sha256': hashlib.sha256(args.output.read_bytes()).hexdigest()}
    args.output.with_suffix('.json').write_text(json.dumps(provenance, indent=2)+'\n')


if __name__ == '__main__':
    main()
