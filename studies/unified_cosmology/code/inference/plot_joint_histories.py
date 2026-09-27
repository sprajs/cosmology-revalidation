"""Render verified expansion and luminosity histories as pointwise intervals.

The renderer accepts already audited arrays. The command-line entry point
performs the independent final audit before it allows any scientific figure.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[4]
COLORS = ['#245f9b', '#268073', '#b56c19', '#8b559a']
LABELS = {
    ('lcdm', 'none'): ('Flat ΛCDM', 'Fixed SN brightness'),
    ('cpl', 'none'): ('Evolving dark energy', 'Fixed SN brightness'),
    ('cpl', 'linear'): ('Evolving dark energy', 'Free linear brightness drift'),
    ('cpl', 'smooth01'): ('Evolving dark energy', 'Smooth brightness\nprior σ = 0.1 mag'),
}


def validate_arrays(cases):
    """Check plotting semantics without inferring scientific qualification."""
    assert cases and len(cases) <= 4
    identities = set()
    for case in cases:
        settings = case['settings']
        key = (settings['model'], settings['evolution'])
        assert key in LABELS and key not in identities
        identities.add(key)
        assert settings['sample'] == 'dovekie'
        assert settings['calibration'] == 'official_planck'
        for coordinate, fields in [('z_expansion', ['q', 'j']), ('z_brightness', ['B'])]:
            z = np.asarray(case[coordinate], dtype=float)
            assert z.ndim == 1 and len(z) >= 2 and z[0] == 0
            assert np.isfinite(z).all() and np.all(np.diff(z) > 0)
            for field in fields:
                values = np.asarray(case[field], dtype=float)
                assert values.shape == (len(z), 5) and np.isfinite(values).all()
                assert np.all(np.diff(values, axis=1) >= 0)
        assert np.array_equal(np.asarray(case['B'])[0], np.zeros(5)), 'B(0)=0 is the declared anchor.'
        if settings['evolution'] == 'none':
            assert np.count_nonzero(case['B']) == 0, 'Fixed brightness is an assumption, not an error band.'


def render(cases, output, *, diagnostic_watermark=None):
    validate_arrays(cases)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(len(cases), 3, figsize=(12.6, 2.05*len(cases)+2.2),
                             squeeze=False, sharex='col', sharey='col')
    titles = ['Scale-factor acceleration', 'Evolution of acceleration', 'Supernova brightness drift']
    fields = ['q', 'j', 'B']
    units = ['$q(z)$', '$j(z)$', '$B(z)$ [mag]']
    for i, case in enumerate(cases):
        key = (case['settings']['model'], case['settings']['evolution'])
        color = COLORS[list(LABELS).index(key)]
        for col, (field, unit) in enumerate(zip(fields, units)):
            ax = axes[i, col]
            z = np.asarray(case['z_brightness' if field == 'B' else 'z_expansion'])
            values = np.asarray(case[field])
            ax.axhline(0, color='#777777', linewidth=.8, linestyle='--', zorder=0)
            if field == 'B' and key[1] == 'none':
                ax.plot(z, values[:, 2], color=color, linewidth=1.8)
                ax.text(.97, .13, 'Fixed at zero', ha='right', va='bottom',
                        transform=ax.transAxes, color=color, fontsize=9)
            else:
                ax.fill_between(z, values[:, 0], values[:, 4], color=color, alpha=.13, linewidth=0)
                ax.fill_between(z, values[:, 1], values[:, 3], color=color, alpha=.27, linewidth=0)
                ax.plot(z, values[:, 2], color=color, linewidth=1.7)
            ax.set_ylabel(unit)
            ax.grid(alpha=.13)
            ax.set_xlim(z[0], z[-1])
            if i == 0:
                ax.set_title(titles[col], fontsize=11, pad=13)
            if i == len(cases)-1:
                ax.set_xlabel('Redshift z')
        title, subtitle = LABELS[key]
        axes[i, 0].text(-.39, .66, title+'\n'+subtitle, transform=axes[i, 0].transAxes,
                        ha='right', va='center', fontsize=10, color=color)
    # Common axes within each column make model-dependent broadening visible.
    # Matplotlib's shared autoscaling includes every row, without clipping tails.
    fig.suptitle('Expansion history and its dependence on supernova brightness',
                 fontsize=15, x=.57, y=.98)
    fig.text(.57, .944, 'Dovekie · DESI DR2 · Planck + ACT + SPT', ha='center', color='#444444')
    fig.text(.57, .059,
             'Lines: medians. Dark/light bands: pointwise 68%/95% conditional credible intervals.\n'
             'q < 0: acceleration. j < 0: scale-factor acceleration decreases with time. Positive B: dimmer SNe.\n'
             'Larger redshift = earlier. Brightness drift is a sensitivity model, not an observed age correction.',
             ha='center', va='center', fontsize=9, color='#444444', linespacing=1.55)
    if diagnostic_watermark:
        fig.text(.57, .50, diagnostic_watermark, ha='center', va='center', rotation=18,
                 fontsize=30, color='#9d2525', alpha=.28)
    fig.subplots_adjust(left=.29, right=.99, bottom=.14, top=.87, wspace=.27, hspace=.24)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180, facecolor='white')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    import final_audit
    audit = final_audit.audit_all()
    assert len(audit['cohorts']) == 4, 'The figure is defined for the four declared targets.'
    inputs = dict(audit['source_sha256'])
    cases = []
    for cohort in audit['cohorts']:
        assert not cohort['integrity_errors'], 'Scientific evidence failed identity checks.'
        assert cohort['parent_qualification'] == 'qualified_at_declared_native_accuracy', \
            f"No qualified parent yet for {cohort['label']}: {cohort['parent_qualification']}"
        assert cohort['native_precision_screen'] == 'diagnostic', 'Native precision screen must complete without flags.'
        records = {}
        for child_name in ['expansion_history', 'luminosity_history']:
            child = cohort['verified_children'][child_name]
            assert child['qualified'] is True
            path = ROOT/child['path']
            assert hashlib.sha256(path.read_bytes()).hexdigest() == child['sha256']
            record = json.loads(path.read_text())
            assert record['settings'] == cohort['settings']
            assert record['target_identity'] == cohort['target_identity']
            records[child_name] = record
            inputs[child['path']] = child['sha256']
        expansion, brightness = records['expansion_history'], records['luminosity_history']
        for field in ['q', 'j']:
            assert expansion['history'][field]['quantile_probabilities'] == [.025, .16, .5, .84, .975]
        cases.append({'settings': cohort['settings'], 'z_expansion': expansion['redshift'],
                      'q': expansion['history']['q']['quantiles_by_redshift_or_scalar'],
                      'j': expansion['history']['j']['quantiles_by_redshift_or_scalar'],
                      'z_brightness': brightness['redshift'],
                      'B': [row['quantiles_025_16_50_84_975'] for row in brightness['pointwise_magnitude_posterior']]})
    render(cases, args.output)
    for path, expected in inputs.items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == expected, 'Figure input changed.'
    provenance = {'status': 'generated_from_audited_conditional_histories',
                  'audit_time_utc': audit['audit_time_utc'], 'audit_status': audit['status'],
                  'inputs_sha256': inputs,
                  'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  'figure': str(args.output.resolve().relative_to(ROOT)),
                  'figure_sha256': hashlib.sha256(args.output.read_bytes()).hexdigest(),
                  'settings': [case['settings'] for case in cases],
                  'interpretation': 'Pointwise credible intervals under four separate models. '
                      'Brightness is a sensitivity assumption, not an identified age correction. '
                      'Every parent and both histories revalidated; fixed32 native precision screen has no flags.'}
    args.output.with_suffix('.json').write_text(json.dumps(provenance, indent=2)+'\n')


if __name__ == '__main__':
    main()
