#!/usr/bin/env python3
"""Replay all 19 workflows, five alternatives, five second seeds, then audit.

Use a new prefix. Existing result directories are never overwritten. Reports in
validation/reports are refreshed only after successful calculations; preserve
the dated Git record before starting a new validation campaign.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prefix', required=True)
    parser.add_argument('--plan-only', action='store_true')
    args = parser.parse_args()
    base = args.prefix
    if Path(base).name != base or base in {'.', '..', ''}:
        parser.error('prefix must be a single directory name')
    variants = [
        ('joint', 'cosmology', 'joint-cpl.json', {}),
        ('qbins', 'cosmology', 'qbins.json', {}),
        ('age-r19', 'ages', 'ages-r19.json', {}),
        ('des-colour', 'des-predictors', 'des-colour.json', {}),
        ('template', 'cosmology', 'template-sensitivity.json', {}),
        ('lcdm-second-seed', 'cosmology', None, {'seed': 41002}),
        ('joint-second-seed', 'cosmology', 'joint-cpl.json', {'seed': 41006}),
        ('qbins-second-seed', 'cosmology', 'qbins.json', {'seed': 41009}),
        ('template-second-seed', 'cosmology', 'template-sensitivity.json', {'seed': 41011}),
        ('des-second-seed', 'des-predictors', None, {'seed': 2026092602}),
    ]
    plans = []
    for suffix, workflow, filename, changes in variants:
        config = json.loads((ROOT/'configs'/filename).read_text()) if filename else {}
        config.update(changes)
        if suffix.startswith('template'):
            config['correction_csv'] = f'results/{base}-v1-populations/median_correction.csv'
        plans.append((suffix, workflow, config))
    if args.plan_only:
        print(json.dumps({'default_replay_prefix': base+'-v1', 'alternatives_and_seeds': plans}, indent=2))
        return
    if list((ROOT/'results').glob(base+'-*')):
        parser.error('result prefix already exists; choose a new prefix')
    config_dir = ROOT/'validation/configs'/base
    config_dir.mkdir(parents=True, exist_ok=False)

    def run(*argv):
        subprocess.run([sys.executable, *map(str, argv)], cwd=ROOT, check=True)

    run(ROOT/'research.py', 'verify')
    run(ROOT/'replay.py', '--prefix', base+'-v1')
    for suffix, workflow, config in plans:
        path = config_dir/(suffix+'.json')
        path.write_text(json.dumps(config, indent=2)+'\n')
        run(ROOT/'research.py', 'run', workflow, '--name', base+'-'+suffix, '--config', path)
    run(ROOT/'validation/core_checks.py', '--prefix', base+'-v1', '--second-seed', base+'-lcdm-second-seed')
    run(ROOT/'validation/des_checks.py', '--prefix', base+'-v1', '--colour', base+'-des-colour', '--second-seed', base+'-des-second-seed')
    run(ROOT/'validation/raisin_checks.py', '--prefix', base+'-v1')
    # This audit uses frozen reference outputs; the manifest compares fresh ones.
    run(ROOT/'validation/hst_checks.py')
    run(ROOT/'validation/record.py', '--prefix', base)


if __name__ == '__main__':
    main()
