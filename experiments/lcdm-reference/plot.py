"""Plot retained external CLASS predictions without rerunning the solver."""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def consumed(path):
    raw = path.read_bytes()
    return raw, {'path': str(path.resolve()), 'bytes': len(raw),
                 'sha256': hashlib.sha256(raw).hexdigest()}


def identity(path):
    return consumed(path)[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--attempt', type=Path, required=True)
    parser.add_argument('--record-sha256', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    record_path = args.attempt / 'record.json'
    record_raw, record_pin = consumed(record_path)
    if record_pin['sha256'] != args.record_sha256:
        raise ValueError('attempt record identity differs')
    record = json.loads(record_raw)
    if record['status'] != 'completed' or record['gates']['execution'] != 'passed':
        raise ValueError('attempt did not complete')
    if [row['id'] for row in record['cases']] != ['anchor', 'precision', 'ns-minus', 'ns-plus']:
        raise ValueError('the four reviewed cases are required in original order')
    curves, input_pins = {}, [record_pin]
    for row in record['cases']:
        p = Path(row['products']['path'])
        raw, pin = consumed(p)
        if any(pin[k] != row['products'][k] for k in ('path', 'bytes', 'sha256')):
            raise ValueError('retained prediction product differs')
        curves[row['id']] = json.loads(raw)
        input_pins.append(pin)
    anchor = curves['anchor']
    ell = np.asarray(anchor['cmb']['ell'])
    for case, product in curves.items():
        if product['cmb']['ell'] != ell.tolist() or product['cmb']['lensed'] is not True:
            raise ValueError('multipole/lensing identity differs')
        if product['cmb']['raw_units'] != 'dimensionless Dl=l(l+1)Cl/(2pi)':
            raise ValueError('unexpected spectrum units')
        for s in ('tt', 'ee', 'te'):
            y = np.asarray(product['cmb']['spectra'][s]['Dl_uK2'])
            if y.shape != ell.shape or not np.isfinite(y).all():
                raise ValueError('invalid CMB plot data')
        for law in ('linear_matter', 'nonlinear_matter'):
            if product[law][0]['z'] != 0:
                raise ValueError('matter reference epoch differs')
    args.output.mkdir(parents=True, exist_ok=False)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.labelcolor': '#252a31', 'text.color': '#252a31',
                         'axes.edgecolor': '#737981', 'grid.color': '#e2e5e9',
                         'grid.linewidth': .6, 'svg.fonttype': 'none',
                         'path.simplify': False})
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.4))
    styles = {'anchor': ('#22639e', '-', 1.6, 'Anchor: n_s = 0.9660499'),
              'precision': ('#737981', ':', 1.2, 'Changed precision; same physics'),
              'ns-minus': ('#c47324', '--', 1.15, 'n_s − 0.005; other inputs fixed'),
              'ns-plus': ('#697f38', '-.', 1.15, 'n_s + 0.005; other inputs fixed')}
    for ax, spectrum in zip(axes.flat[:3], ('tt', 'ee', 'te')):
        for case in ('ns-minus', 'ns-plus', 'anchor', 'precision'):
            color, line, width, label = styles[case]
            ax.plot(ell, curves[case]['cmb']['spectra'][spectrum]['Dl_uK2'],
                    color=color, linestyle=line, linewidth=width, label=label)
        ax.set(title=f'Lensed {spectrum.upper()}', xlabel='Multipole ℓ',
               ylabel='D_ℓ = ℓ(ℓ+1) C_ℓ / (2π) [µK²]', xlim=(2, 2508))
        if spectrum != 'te':
            ax.set_ylim(bottom=0)
        else:
            ax.axhline(0, color='#737981', linewidth=.7)
        ax.grid(axis='y')
    ax = axes.flat[3]
    for case in ('ns-minus', 'ns-plus', 'anchor', 'precision'):
        color, line, width, label = styles[case]
        product = curves[case]['linear_matter'][0]
        ax.loglog(product['k_1_Mpc'], product['P_Mpc3'], color=color,
                  linestyle=line, linewidth=width, label=label)
    nonlinear = anchor['nonlinear_matter'][0]
    ax.loglog(nonlinear['k_1_Mpc'], nonlinear['P_Mpc3'], color='#353b43',
              linestyle=(0, (5, 3)), linewidth=1, label='Anchor Halofit at z = 0')
    ax.set(title='Matter power at z = 0', xlabel='k [Mpc⁻¹]', ylabel='P(k) [Mpc³]')
    ax.grid(which='major', axis='y')
    fig.suptitle('One cosmology: CMB and matter predictions', fontsize=19, y=.977)
    fig.text(.5, .933,
             'CLASS 3.3.0 • H₀ = 67.32117 km s⁻¹ Mpc⁻¹ • ω_b = 0.02238280 • ω_cdm = 0.1201075',
             ha='center', fontsize=10)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    order = [2, 0, 1, 3]
    handles, labels = [handles[i] for i in order], [labels[i] for i in order]
    handles += [ax.lines[-1]]; labels += ['Anchor Halofit at z = 0']
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(.5, .902),
               ncol=3, frameon=False, fontsize=9)
    fig.subplots_adjust(left=.09, right=.975, top=.80, bottom=.18, hspace=.40, wspace=.27)
    fig.text(.09, .048,
             'Illustrative fixed state with one 0.06 eV relic, fixed YHe, HyRec and declared reionization. No observational overlay or fit.\n'
             'Curves retain emitted sample grids. Precision changes are empirical comparisons; solver and interpolation bounds are unavailable.',
             fontsize=9, linespacing=1.5)
    files = []
    for extension in ('png', 'svg'):
        p = args.output / ('class-reference-predictions.' + extension)
        fig.savefig(p, dpi=160, facecolor='white')
        files.append(identity(p))
    plt.close(fig)
    source_pin = identity(Path(__file__))
    manifest = {'kind': 'actual-class-prediction-figure', 'input_identities': input_pins,
                'source': source_pin, 'outputs': files,
                'runtime': {'python': sys.version, 'matplotlib': matplotlib.__version__,
                            'numpy': np.__version__,
                            'packages': {name: importlib.metadata.version(name) for name in
                                         ('matplotlib', 'numpy', 'pillow', 'contourpy', 'fonttools')}},
                'plotted': 'All emitted ell2..2508 lensed TT/EE/TE; original z0 linear grids for all four cases plus anchor Halofit.',
                'selection': 'No downsampling, reordering, interpolation or extrapolation by this plotting source.',
                'error_bounds': None, 'inference': 'not established', 'interpretation': 'not qualified'}
    (args.output / 'figure-manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps({'outputs': files, 'manifest': str(args.output / 'figure-manifest.json')}))


if __name__ == '__main__':
    main()
