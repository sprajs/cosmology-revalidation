#!/usr/bin/env python3
"""Frozen pre-outcome source and baseline replay for eight archived DES simulations."""
import argparse
import collections
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time

ROOT = Path('/home/szymon/Documents/ChatGPT/supernova')
OUT = Path('/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/astra_design/raisin_timing_assets/baseline_2021')
DESIGN = OUT / 'protocol.json'
LC = ROOT / 'runs/research_2026_09_26/raisin_simulation_execution_review/nominal_fitres/lcplot-feasibility/FITOPT000.LCPLOT.gz'
ARCH = Path('/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/astra_design/raisin_timing_assets/author_20211111/nir.FITRES.gz')
NML = ROOT / 'runs/research_2026_09_26/raisin_simulation_assets/timing/REFAC_DES_RAISIN_NIR.nml'
SIMLIB = ROOT / 'runs/research_2026_09_26/raisin_sign_source/sim/simlibs/DES_RAISIN.simlib'
RELEASE = ROOT / 'sources/repos/djones1040__RAISIN_DataRelease@a383c4b'
BIN = Path('/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/astra_design/raisin_timing_assets/snana_v11_04d/build/bin/snlc_fit.exe')
PROTOCOL_SHA = '399dc9562d8c19af2e18b153cd887f1a576ec2ace9d04c96aa266c56fcc57185'
CIDS = tuple(str(i) for i in range(1, 9))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2) + '\n')


def source_gate():
    assert sha(DESIGN) == PROTOCOL_SHA
    d = json.loads(DESIGN.read_text())
    assert d['pilot_CIDs'] == list(CIDS)
    for p, digest in d['inputs_sha256'].items():
        assert sha(ROOT / p) == digest, p
    assert not any(re.search(r'(?m)^\s*OPT_SETPKMJD\s*=', NML.read_text()) for _ in [0])
    assert 'FUDGE_DATAERR_SCALE' not in NML.read_text()
    assert 'FLUXERRMODEL_FILE' not in NML.read_text()
    assert 'FUDGE_FLUXERR_SCALE' not in NML.read_text()
    return d


def read_fitres(path):
    d = {}
    with gzip.open(path, 'rt') as f:
        names = None
        for line in f:
            if line.startswith('VARNAMES:'):
                names = line.split()[1:]
            elif line.startswith('SN:'):
                row = dict(zip(names, line.split()[1:]))
                if row['CID'] in CIDS:
                    assert row['CID'] not in d
                    d[row['CID']] = row
    assert set(d) == set(CIDS)
    return d


def read_lc():
    result = collections.defaultdict(list)
    with gzip.open(LC, 'rt') as f:
        names = f.readline().split()[1:]
        for line in f:
            if not line.startswith('OBS:'):
                continue
            row = dict(zip(names, line.split()[1:]))
            if row['CID'] in CIDS and row['DATAFLAG'] == '1':
                assert row['BAND'] in ('J', 'H')
                result[row['CID']].append(row)
    assert set(result) == set(CIDS)
    return result


def simlib_coords():
    d, current = {}, None
    for line in SIMLIB.read_text().splitlines():
        m = re.match(r'LIBID:\s*(\d+)', line)
        if m:
            current = m.group(1)
            assert current not in d
            d[current] = {}
        if current is not None:
            m = re.search(r'\bRA:\s*(\S+)\s+DECL:\s*(\S+)', line)
            if m:
                d[current]['RA'], d[current]['DEC'] = m.group(1), m.group(2)
    assert all('RA' in v and 'DEC' in v for v in d.values())
    return d


def set_key(text, key, val):
    pat = rf'(?m)^(\s*{re.escape(key)}\s*=).*$'
    text, n = re.subn(pat, lambda m: m.group(1) + ' ' + val, text)
    assert n == 1, (key, n)
    return text


def make_phot(cid, fit, rows, coord):
    r = fit[cid]
    lines = [
        'SURVEY: DES', f'SNID: {cid}', 'SNTYPE: 1', 'FAKE: 2',
        'FILTERS: grizJH', f"RA: {coord['RA']}", f"DEC: {coord['DEC']}",
        'PIXSIZE: 0.27', f"MWEBV: {r['MWEBV']}",
        f"REDSHIFT_HELIO: {r['zHEL']} +- {r['zHELERR']}",
        f"REDSHIFT_FINAL: {r['zCMB']} +- {r['zCMBERR']}",
        f"VPEC: {r['VPEC']} +- {r['VPECERR']}",
        f"PEAKMJD: {r['PKMJDINI']}",
        'SIM_MODEL_NAME: snoopy.B18',
        f"SIM_TYPE_INDEX: {r['SIM_TYPE_INDEX']}",
        f"SIM_TEMPLATE_INDEX: {r['SIM_TEMPLATE_INDEX']}",
        f"SIM_LIBID: {r['SIM_LIBID']}",
        f"SIM_NGEN_LIBID: {r['SIM_NGEN_LIBID']}",
        f"SIM_REDSHIFT_CMB: {r['SIM_ZCMB']}",
        f"SIM_VPEC: {r['SIM_VPEC']}",
        f"SIM_PEAKMJD: {r['SIM_PKMJD']}",
        f"SIM_STRETCH: {r['SIM_STRETCH']}",
        f"SIM_AV: {r['SIM_AV']}",
        f"SIM_RV: {r['SIM_RV']}",
        f"SIM_DLMAG: {r['SIM_DLMAG']}",
        'SIMOPT_MWCOLORLAW: 94',
        f'NOBS: {len(rows[cid])}', 'NVAR: 5',
        'VARLIST: MJD FLT FIELD FLUXCAL FLUXCALERR',
    ]
    lines += [f"OBS: {e['MJD']} {e['BAND']} NULL {e['FLUXCAL']} {e['FLUXCAL_ERR']}" for e in rows[cid]]
    lines.append('END:')
    return '\n'.join(lines) + '\n'


def prepare():
    d = source_gate()
    assert not (OUT / 'baseline-freeze.json').exists()
    fit, rows, coords = read_fitres(ARCH), read_lc(), simlib_coords()
    for cid in CIDS:
        assert fit[cid]['SIM_LIBID'] in coords
        assert len(rows[cid]) >= 2
    raw = NML.read_text()
    assert '&SNLCINP' in raw
    raw = '  &SNLCINP' + raw.split('&SNLCINP', 1)[1]
    ledger = []
    for condition in ('baseline', 'baseline_copy'):
        work = OUT / 'fits' / condition
        version = work / 'data' / 'DES_RAISIN_SIM'
        version.mkdir(parents=True, exist_ok=False)
        (version / 'DES_RAISIN_SIM.README').write_text('DOCUMENTATION:\n  PURPOSE: Rounded archived nominal NIR LCPLOT replay, CID 1–8\nDOCUMENTATION_END:\n')
        (version / 'DES_RAISIN_SIM.LIST').write_text(''.join(f'{cid}.DAT\n' for cid in CIDS))
        for cid in CIDS:
            path = version / f'{cid}.DAT'
            path.write_text(make_phot(cid, fit, rows, coords[fit[cid]['SIM_LIBID']]))
            ledger.append({'condition': condition, 'CID': cid, 'file': str(path.relative_to(ROOT)),
                           'sha256': sha(path), 'n_rows': len(rows[cid]),
                           'synthetic_metadata': 'FAKE2, SURVEY/FILTERS, SIM_MODEL_NAME, PIXSIZE, SNTYPE, FIELD NULL, README; SIM_REDSHIFT_HELIO omitted'})
        (work / 'kcor.fits').symlink_to(RELEASE / 'kcor/kcor_DES_NIR.fits')
        (work / 'snoopy.B18').mkdir()
        (work / 'snoopy.B18/snoopy.info').symlink_to(ROOT / 'phase2/official/inputs/SNDATA_ROOT/models/snoopy/snoopy.B18/snoopy.info')
        (work / 'snoopy.B18/SNooPy_test.fits').symlink_to(RELEASE / 'model/snoopy.B18/SNooPy_B18.fits')
        nml = raw.replace('&SNLCINP', "&SNLCINP\n PRIVATE_DATA_PATH = 'data'", 1)
        for key, val in {'VERSION_PHOTOMETRY': "'DES_RAISIN_SIM'", 'KCOR_FILE': "'kcor.fits'",
                         'TEXTFILE_PREFIX': "'fit'", 'SNTABLE_LIST': "'FITRES(text:key) LCPLOT(text:key)'",
                         'FITMODEL_NAME': "'./snoopy.B18'"}.items():
            nml = set_key(nml, key, val)
        (work / 'fit.nml').write_text(nml)
    for cid in CIDS:
        a = OUT / 'fits/baseline/data/DES_RAISIN_SIM' / f'{cid}.DAT'
        b = OUT / 'fits/baseline_copy/data/DES_RAISIN_SIM' / f'{cid}.DAT'
        assert a.read_bytes() == b.read_bytes()
    save(OUT / 'input-ledger.json', ledger)
    paths = [OUT / 'baseline_runner.py', OUT / 'baseline_gate.py', DESIGN, ARCH, LC, NML, SIMLIB, BIN,
             RELEASE / 'kcor/kcor_DES_NIR.fits', RELEASE / 'model/snoopy.B18/SNooPy_B18.fits',
             ROOT / 'phase2/official/inputs/SNDATA_ROOT/models/snoopy/snoopy.B18/snoopy.info', OUT / 'input-ledger.json']
    paths += sorted(p for p in (OUT / 'fits').rglob('*') if p.is_file() and not p.is_symlink())
    save(OUT / 'baseline-freeze.json', {'status': 'frozen before native fits', 'design_sha256': PROTOCOL_SHA,
                                        'files': {str(p.relative_to(ROOT)): sha(p) for p in paths},
                                        'cohort': CIDS, 'metadata_gaps': ['Original SIM_REDSHIFT_HELIO absent from FITRES; not substituted',
                                                                        'Original per-epoch PHOTFLAG absent',
                                                                        'LCPLOT flux/error/MJD rounded'],
                                        'nml_changes': ['PRIVATE_DATA_PATH data location', 'VERSION_PHOTOMETRY same value',
                                                        'KCOR_FILE path only', 'TEXTFILE_PREFIX short output path',
                                                        'SNTABLE_LIST text output format', 'FITMODEL_NAME path only']})
    print(json.dumps({'prepared': len(ledger), 'baseline_freeze_sha256': sha(OUT / 'baseline-freeze.json')}))


def run(condition):
    release = json.loads((OUT / 'execution-release.json').read_text())
    assert release['protocol_sha256'] == PROTOCOL_SHA
    d = source_gate()
    freeze = json.loads((OUT / 'baseline-freeze.json').read_text())
    for p, digest in freeze['files'].items():
        assert sha(ROOT / p) == digest, p
    work = OUT / 'fits' / condition
    assert condition in ('baseline', 'baseline_copy')
    assert not (work / 'fit.log').exists()
    env = os.environ.copy()
    sysroot = ROOT / 'phase2/official/build/sysroot/usr'
    env.update(SNANA_DIR='/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/astra_design/raisin_timing_assets/snana_v11_04d/build',
               SNDATA_ROOT=str(ROOT / 'phase2/official/inputs/SNDATA_ROOT'),
               LD_LIBRARY_PATH=str(sysroot / 'lib'), OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    begin = time.monotonic()
    with (work / 'fit.log').open('x') as stream:
        try:
            proc = subprocess.run([str(BIN), 'fit.nml'], cwd=work, env=env,
                                  stdout=stream, stderr=subprocess.STDOUT, timeout=600)
            code = proc.returncode
        except subprocess.TimeoutExpired:
            code = 'timeout'
    save(work / 'execution.json', {'returncode': code, 'wall_seconds': time.monotonic() - begin,
                                   'design_sha256': PROTOCOL_SHA,
                                   'freeze_sha256': sha(OUT / 'baseline-freeze.json'),
                                   'log_sha256': sha(work / 'fit.log')})
    print(json.dumps({'condition': condition, 'returncode': code,
                      'wall_seconds': json.loads((work / 'execution.json').read_text())['wall_seconds']}))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('phase', choices=['prepare', 'run'])
    parser.add_argument('--condition', choices=['baseline', 'baseline_copy'])
    args = parser.parse_args()
    if args.phase == 'prepare':
        prepare()
    else:
        assert args.condition
        run(args.condition)


if __name__ == '__main__':
    main()
