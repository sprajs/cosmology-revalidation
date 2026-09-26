#!/usr/bin/env python3
"""Prepared, gate-sequential native CSP 42-object filter response runner.

Preparation performs no native fits. Do not invoke run-baseline until the
separately reviewed pilot instrumentation/stability gate has passed.
"""
import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
PREP = ROOT / 'runs/research_2026_09_26/csp_native_filter_response/cohort-preparation'
DESIGN = ROOT / 'runs/research_2026_09_26/raisin_profile_solver_review/csp-filter-interpretation/native-design'
OUT = ROOT / 'runs/research_2026_09_26/csp_native_filter_response/cohort-execution'
ARCH = DESIGN / 'CSP-nominal.FITRES.gz'
BIN = ROOT / 'phase2/official/build/SNANA-v11_04k/bin/snlc_fit.exe'
PREP_SHA = '5a51cb68baf51755cf52963164e1961f8704784b2abcccec3d659b1487c5fe44'
DESIGN_SHA = '1947050a3ac68f6a995a3d5adc93390266127feb519318e4d4947d9eaf1496a2'
MAX_ACTIVE = 600.0
ALREADY_USED = 1.331

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def dump(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')

def table(text, marker):
    names = None
    for line in text.splitlines():
        if line.startswith('VARNAMES:'):
            names = line.split()[1:]
        elif line.startswith(marker):
            vals = line.split()[1:]
            assert names and len(names) == len(vals), line[:200]
            yield dict(zip(names, vals))

def source_gate():
    assert sha(PREP / 'protocol.json') == PREP_SHA
    assert sha(DESIGN / 'protocol.json') == DESIGN_SHA
    p = json.loads((PREP / 'protocol.json').read_text())
    for path, digest in p['inputs_sha256'].items():
        assert sha(ROOT / path) == digest, path
    assert len(p['jobs']) == 6
    assert p['jobs'][3]['objects'] == 42
    assert p['jobs'][5]['changed_rows'] == 188
    design = json.loads((DESIGN / 'protocol.json').read_text())
    assert len(design['affected_CIDs']) == 32
    assert len(design['unaffected_CIDs']) == 10
    assert set(design['affected_CIDs']).isdisjoint(design['unaffected_CIDs'])
    return p, design

def all_files(path):
    return sorted(p for p in path.rglob('*') if p.is_file())

def prepare():
    prep, design = source_gate()
    assert not OUT.exists(), 'Refuse to replace prepared execution'
    OUT.mkdir()
    shutil.copytree(PREP / 'full', OUT / 'full')
    copied = {}
    for src in all_files(PREP / 'full'):
        rel = src.relative_to(PREP / 'full')
        dst = OUT / 'full' / rel
        assert src.read_bytes() == dst.read_bytes()
        copied[str(rel)] = sha(dst)
    assert len(copied) == sum(1 for _ in all_files(PREP / 'full'))
    shutil.copy2(Path(__file__), OUT / 'executed-source.py')
    protocol = {
        'status': 'Prepared only; no full-cohort native fit run',
        'frozen_design_sha256': DESIGN_SHA,
        'cohort_preparation_sha256': PREP_SHA,
        'runner_sha256': sha(Path(__file__)),
        'copied_file_sha256': copied,
        'archived_raw_fitres_sha256': sha(ARCH),
        'native_binary_sha256': sha(BIN),
        'arms_order': ['nominal', 'nominal_copy', 'known_Jdw_to_j'],
        'active_seconds_cap': MAX_ACTIVE,
        'prior_pilot_active_seconds': ALREADY_USED,
        'baseline_gate': design['baseline_gates']['archive'],
        'mask_gate': design['native']['mask'],
        'controls': design['unaffected_CIDs'],
        'notes': 'External root clearance and independently reviewed pilot stability gate required before run-baseline. Relabel arm requires all42 baseline and copy gates. No hidden CID exclusion.'
    }
    dump(OUT / 'execution-protocol.json', protocol)
    print('prepared_protocol_sha256', sha(OUT / 'execution-protocol.json'))

def execution_gate(gate_path, gate_sha):
    source_gate()
    p = json.loads((OUT / 'execution-protocol.json').read_text())
    assert p['runner_sha256'] == sha(Path(__file__))
    assert p['cohort_preparation_sha256'] == PREP_SHA
    for rel, digest in p['copied_file_sha256'].items():
        assert sha(OUT / 'full' / rel) == digest, rel
    assert gate_path and gate_sha, 'Require independently reviewed pilot gate path and SHA-256'
    assert sha(gate_path) == gate_sha, 'Pilot gate hash mismatch'
    gate = json.loads(Path(gate_path).read_text())
    assert gate.get('pass') is True and gate.get('approved_for_full42') is True, 'Pilot gate not released'
    return p, gate

def env_native():
    build = ROOT / 'phase2/official/build/SNANA-v11_04k'
    env = os.environ.copy()
    env.update(SNANA_DIR=str(build), SNDATA_ROOT=str(ROOT / 'phase2/official/inputs/SNDATA_ROOT'),
               LD_LIBRARY_PATH=str(ROOT / 'phase2/official/build/sysroot/usr/lib'),
               OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    return env

def active_used():
    total = ALREADY_USED
    for arm in ('nominal', 'nominal_copy', 'known_Jdw_to_j'):
        file = OUT / 'full/fits' / arm / 'execution.json'
        if file.exists():
            total += json.loads(file.read_text())['seconds']
    return total

def run_arm(arm, gate_path, gate_sha):
    p, gate = execution_gate(gate_path, gate_sha)
    if arm == 'known_Jdw_to_j':
        result = json.loads((OUT / 'baseline-review.json').read_text())
        assert result['pass'] is True, 'All42 nominal baseline gate required'
        assert result['protocol_sha256'] == sha(OUT / 'execution-protocol.json')
    if arm == 'nominal_copy':
        assert (OUT / 'full/fits/nominal/execution.json').exists()
    if arm == 'known_Jdw_to_j':
        assert (OUT / 'full/fits/nominal_copy/execution.json').exists()
    work = OUT / 'full/fits' / arm
    log = work / 'fit.log'
    assert not log.exists(), 'Never overwrite a native execution'
    remaining = MAX_ACTIVE - active_used()
    assert remaining > 1, 'Active-time cap exhausted'
    start = time.monotonic()
    with log.open('x') as f:
        try:
            proc = subprocess.run([str(BIN), 'fit.nml'], cwd=work, env=env_native(),
                                  stdout=f, stderr=subprocess.STDOUT, timeout=remaining)
            rc = proc.returncode
        except subprocess.TimeoutExpired:
            rc = 124
    elapsed = time.monotonic() - start
    record = {'arm': arm, 'returncode': rc, 'seconds': elapsed,
              'cumulative_active_seconds': active_used() + elapsed,
              'protocol_sha256': sha(OUT / 'execution-protocol.json'),
              'pilot_gate_sha256': gate_sha, 'log_sha256': sha(log),
              'fitres_sha256': sha(work / 'fit.FITRES.TEXT') if (work / 'fit.FITRES.TEXT').exists() else None,
              'lcplot_sha256': sha(work / 'fit.LCPLOT.TEXT') if (work / 'fit.LCPLOT.TEXT').exists() else None}
    dump(work / 'execution.json', record)
    print(json.dumps(record, indent=2))
    assert rc == 0, 'Native execution failure preserved'

def fitres(path):
    return {r['CID']: r for r in table(Path(path).read_text(), 'SN:')}

def archive_fitres():
    with gzip.open(ARCH, 'rt') as f:
        return {r['CID']: r for r in table(f.read(), 'SN:')}

def review_baseline():
    protocol = json.loads((OUT / 'execution-protocol.json').read_text())
    expected = set(json.loads((DESIGN / 'protocol.json').read_text())['affected_CIDs'] +
                   json.loads((DESIGN / 'protocol.json').read_text())['unaffected_CIDs'])
    archived = archive_fitres()
    a = fitres(OUT / 'full/fits/nominal/fit.FITRES.TEXT')
    b = fitres(OUT / 'full/fits/nominal_copy/fit.FITRES.TEXT')
    rows = []
    for cid in sorted(expected):
        old, x, y = archived.get(cid), a.get(cid), b.get(cid)
        row = {'CID': cid, 'archive_present': old is not None,
               'nominal_present': x is not None, 'copy_present': y is not None}
        if old and x and y:
            row.update({'delta_DLMAG': float(x['DLMAG'])-float(old['DLMAG']),
                        'delta_FITCHI2': float(x['FITCHI2'])-float(old['FITCHI2']),
                        'NDOF_equal': x['NDOF'] == old['NDOF'],
                        'ERRFLAG_zero': x['ERRFLAG_FIT'] == y['ERRFLAG_FIT'] == '0',
                        'fixed_parameters': all(float(x[k]) == 0 for k in ('STRETCHERR','AVERR','PKMJDERR','RVERR')) and
                                            float(x['STRETCH']) == 1 and float(x['AV']) == 0 and float(x['RV']) == 1.518,
                        'copy_row_equal': x == y})
            row['archive_gate'] = (abs(row['delta_DLMAG']) <= .001 and
                                   abs(row['delta_FITCHI2']) <= max(.01,.001*abs(float(old['FITCHI2']))) and
                                   row['NDOF_equal'] and row['ERRFLAG_zero'] and row['fixed_parameters'])
        else:
            row['archive_gate'] = False
        rows.append(row)
    lc_a = OUT / 'full/fits/nominal/fit.LCPLOT.TEXT'
    lc_b = OUT / 'full/fits/nominal_copy/fit.LCPLOT.TEXT'
    result = {'pass': len(rows)==42 and set(a)==set(b)==expected and all(r['archive_gate'] and r.get('copy_row_equal') for r in rows)
              and lc_a.read_bytes() == lc_b.read_bytes(),
              'protocol_sha256': sha(OUT / 'execution-protocol.json'),
              'copy_LCPLOT_bytes_equal': lc_a.read_bytes() == lc_b.read_bytes(),
              'rows': rows,
              'still_required_before_relabel': 'Independent reviewed native optimizer/support/final-state gate; root release recorded separately.'}
    dump(OUT / 'baseline-review.json', result)
    print('baseline_pass', result['pass'], 'n', len(rows))

def physical_lc_rows(path):
    rows = {}
    for r in table(Path(path).read_text(), 'OBS:'):
        if r['DATAFLAG'] != '1':
            continue
        cid = r['CID']
        key = (r['MJD'], r['BAND'], r['FLUXCAL'], r['FLUXCAL_ERR'])
        rows.setdefault(cid, []).append(key)
    return rows

def review_relabel():
    base = json.loads((OUT / 'baseline-review.json').read_text())
    assert base['pass'] is True
    design = json.loads((DESIGN / 'protocol.json').read_text())
    a = fitres(OUT / 'full/fits/nominal/fit.FITRES.TEXT')
    b = fitres(OUT / 'full/fits/known_Jdw_to_j/fit.FITRES.TEXT')
    expected = set(design['affected_CIDs'] + design['unaffected_CIDs'])
    mask_a = physical_lc_rows(OUT / 'full/fits/nominal/fit.LCPLOT.TEXT')
    mask_b = physical_lc_rows(OUT / 'full/fits/known_Jdw_to_j/fit.LCPLOT.TEXT')
    allowed = set()
    with (PREP / 'changed-row-ledger.csv').open() as f:
        for row in csv.DictReader(f):
            if row['scope'] == 'full':
                allowed.add((row['CID'], row['MJD']))
    rows = []
    for cid in sorted(expected):
        x,y = a.get(cid), b.get(cid)
        item={'CID':cid,'affected':cid in design['affected_CIDs'],'nominal_present':x is not None,'relabel_present':y is not None}
        if x and y:
            item['delta_DLMAG']=float(y['DLMAG'])-float(x['DLMAG'])
            item['ERRFLAG_zero']=x['ERRFLAG_FIT']==y['ERRFLAG_FIT']=='0'
            item['fixed_parameters']=all(float(y[k])==0 for k in ('STRETCHERR','AVERR','PKMJDERR','RVERR')) and float(y['STRETCH'])==1 and float(y['AV'])==0 and float(y['RV'])==1.518
            # Every native plotted datum must preserve MJD, flux and error; only listed J->j is allowed.
            aa=sorted((m, f, e, ('J' if band=='j' and (cid,m) in allowed else band)) for m,band,f,e in mask_a.get(cid,[]))
            bb=sorted((m, f, e, ('J' if band=='j' and (cid,m) in allowed else band)) for m,band,f,e in mask_b.get(cid,[]))
            item['plotted_physical_mask_equal_after_undo']=aa==bb
            item['unaffected_science_equal']=x==y if not item['affected'] else None
        rows.append(item)
    pass_gate = (set(a)==set(b)==expected and all(r.get('ERRFLAG_zero') and r.get('fixed_parameters') and r.get('plotted_physical_mask_equal_after_undo') for r in rows)
                 and all(r['unaffected_science_equal'] for r in rows if not r['affected']))
    result={'pass':pass_gate,'rows':rows,'protocol_sha256':sha(OUT/'execution-protocol.json'),
            'scope':'Conditional raw fit response at fixed 42-object membership, before BBC, no cosmology fit.'}
    if pass_gate:
        affected=[r['delta_DLMAG'] for r in rows if r['affected']]
        full=[r['delta_DLMAG'] for r in rows]
        result['mean_deltaD_affected32']=sum(affected)/32
        result['mean_deltaD_full42']=sum(full)/42
        result['conditional_high_minus_low_response']=-result['mean_deltaD_full42']
    dump(OUT/'response-review.json',result)
    print('response_gate',pass_gate)

if __name__ == '__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('action', choices=['prepare','run-baseline','review-baseline','run-relabel','review-relabel'])
    ap.add_argument('--arm', choices=['nominal','nominal_copy'])
    ap.add_argument('--gate-path', type=Path)
    ap.add_argument('--gate-sha')
    args=ap.parse_args()
    if args.action=='prepare': prepare()
    elif args.action=='run-baseline':
        assert args.arm in ('nominal','nominal_copy')
        run_arm(args.arm,args.gate_path,args.gate_sha)
    elif args.action=='review-baseline': review_baseline()
    elif args.action=='run-relabel': run_arm('known_Jdw_to_j',args.gate_path,args.gate_sha)
    elif args.action=='review-relabel': review_relabel()
