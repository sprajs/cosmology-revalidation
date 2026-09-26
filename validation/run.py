#!/usr/bin/env python3
"""Run the complete research and validation campaign in a new results directory."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--name', required=True, help='New campaign name below results/')
    parser.add_argument('--plan-only', action='store_true')
    args = parser.parse_args()
    base = args.name
    if Path(base).name != base or base in {'.', '..', ''}:
        parser.error('name must be a single directory name')
    variants = [
        ('alternatives/joint-cpl', 'cosmology', 'joint-cpl.json', {}),
        ('alternatives/flexible-expansion', 'cosmology', 'qbins.json', {}),
        ('alternatives/age-catalogue-r19', 'ages', 'ages-r19.json', {}),
        ('alternatives/des-colour', 'des-predictors', 'des-colour.json', {}),
        ('alternatives/age-template', 'cosmology', 'template-sensitivity.json', {}),
        ('robustness/lcdm', 'cosmology', None, {'seed': 41002}),
        ('robustness/joint-cpl', 'cosmology', 'joint-cpl.json', {'seed': 41006}),
        ('robustness/flexible-expansion', 'cosmology', 'qbins.json', {'seed': 41009}),
        ('robustness/age-template', 'cosmology', 'template-sensitivity.json', {'seed': 41011}),
        ('robustness/des-predictors', 'des-predictors', None, {'seed': 2026092602}),
    ]
    plans = []
    for suffix, workflow, filename, changes in variants:
        config = json.loads((ROOT/'configs'/filename).read_text()) if filename else {}
        config.update(changes)
        if suffix.endswith('/age-template'):
            config['correction_csv'] = f'results/{base}/baseline/populations/median_correction.csv'
        plans.append((suffix, workflow, config))
    if args.plan_only:
        print(json.dumps({'baseline_directory': base+'/baseline', 'alternatives_and_seeds': plans}, indent=2))
        return
    destination = ROOT/'results'/base
    if destination.exists():
        parser.error('campaign already exists; choose a new name')
    config_dir = destination/'configurations'
    config_dir.mkdir(parents=True, exist_ok=False)

    def run(*argv):
        subprocess.run([sys.executable, *map(str, argv)], cwd=ROOT, check=True)

    run(ROOT/'research.py', 'verify')
    run(ROOT/'replay.py', '--output', base+'/baseline')
    for suffix, workflow, config in plans:
        path = config_dir/(suffix.replace('/', '-')+'.json')
        path.write_text(json.dumps(config, indent=2)+'\n')
        run(ROOT/'research.py', 'run', workflow, '--name', base+'/'+suffix, '--config', path)
    run(ROOT/'validation/core_checks.py', '--results', base+'/baseline', '--second-seed', base+'/robustness/lcdm')
    run(ROOT/'validation/des_checks.py', '--results', base+'/baseline', '--colour', base+'/alternatives/des-colour', '--second-seed', base+'/robustness/des-predictors')
    run(ROOT/'validation/raisin_checks.py', '--results', base+'/baseline')
    run(ROOT/'validation/hst_checks.py', '--results', base+'/baseline')
    run(ROOT/'validation/record.py', '--name', base)


if __name__ == '__main__':
    main()
