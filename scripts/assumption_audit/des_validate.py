#!/usr/bin/env python3
"""Validate audit claims and ensure outputs match the inspected source bytes."""
from pathlib import Path
import hashlib
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'runs/assumption_audit/des'
r=json.loads((OUT/'release_checks.json').read_text())
c=json.loads((OUT/'covariance_results.json').read_text())
assert r['calibration_overlap']['relative_stat_whitened_frobenius_remainder'] < 2e-6
assert r['calibration_overlap']['CAL_SALT3_minus_CALSPEC_positive_rank_at_2e_6']==9
assert r['vpec_increment']['min_stat_whitened_eigenvalue'] < -.04
assert r['vpec_increment']['violating_pair_min_eigenvalue_using_unrounded_release'] < -1e-5
assert r['vpec_increment']['total_min_stat_whitened_eigenvalue']>0
assert r['vpec_increment']['max_abs_deviation_from_5_decimal_quantisation_mag2']<1e-8
assert not r['shipped_likelihood_reader']['zHD_available']
assert all(x['n']==10 and x['number_different']==10 for x in r['PS1_Foundation_wavelength_draws'])
rank={x['name']:x['positive_rank'] for x in c['single_systematics']}
assert all(rank[x]==1 for x in ['W22','MWEBV','COLORLAW','P24a','P24b','P24c'])
assert abs(c['native_flat_lcdm_Omega_m']-.3303168)<2e-6
assert c['stat_offdiagonal_max_abs_mag2']<1e-10
for name,script in [('manifest.json','des_covariance.py'),('release_checks_manifest.json','des_release_checks.py')]:
    manifest=json.loads((OUT/name).read_text())
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    assert manifest['script_sha256']==sha(ROOT/'scripts/assumption_audit'/script)
    assert manifest['script_sha256']==sha(ROOT/manifest['code_snapshot'])
    for path,digest in manifest['inputs'].items(): assert sha(ROOT/path)==digest,path
    for path,digest in manifest['outputs'].items(): assert sha(ROOT/path)==digest,path
(OUT/'code'/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
print('PASS: result checks, source/output hashes, and producing-script snapshot hashes')
