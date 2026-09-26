"""Separate public training-list membership from existing validation outcomes."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT/'runs/research_2026_09_26/astra_design'
OUT = ROOT/'runs/research_2026_09_26/training_overlap'
ROSTER = ROOT/'sources/updates/2026-09-26-training-roster/roster.csv'
HEAD = ROOT/'sources/repos/des-science__DES-SN5YR@1.3/0_DATA/DES-SN5YR_DES/DES-SN5YR_DES_HEAD.FITS.gz'
PROTOCOL = ROOT/'docs/research-2026-09-26/training-overlap-protocol.md'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, obj):
    path.write_text(json.dumps(obj, indent=2)+'\n')


def design():
    OUT.mkdir(exist_ok=False)
    roster = list(csv.DictReader(ROSTER.open()))
    des = {str(int(r['SNID'])): r for r in roster if r['SURVEY'] == 'DES'}
    assert len(des) == sum(r['SURVEY'] == 'DES' for r in roster)
    head = fits.getdata(HEAD, 1)
    ids = [str(int(str(x).strip())) for x in head['SNID']]
    assert len(ids) == len(set(ids))
    h = dict(zip(ids, head))
    joins = []
    for cid, row in des.items():
        separation = None
        if cid in h:
            # Stable spherical great-circle separation, including RA wrapping.
            ra1, dec1 = np.deg2rad([float(row['RA']), float(row['DEC'])])
            ra2, dec2 = np.deg2rad([float(h[cid]['RA']), float(h[cid]['DEC'])])
            hav = np.sin((dec1-dec2)/2)**2 + np.cos(dec1)*np.cos(dec2)*np.sin((ra1-ra2)/2)**2
            separation = float(2*np.arcsin(np.sqrt(np.clip(hav, 0, 1)))*180/np.pi*3600)
            assert separation < 1., (cid, separation)
        joins.append({'CID': cid, 'found_in_DES_HEAD': cid in h, 'separation_arcsec': separation})
    assert all(x['found_in_DES_HEAD'] for x in joins)
    paths = [ROSTER, HEAD, PROTOCOL, Path(__file__), BASE/'shared43/projected-modes.npz', BASE/'validation1020/cohort.csv']
    with np.load(paths[-2]) as data:
        discovery = data['CID'].astype(str)
    validation = list(csv.DictReader(paths[-1].open()))
    assert len(discovery) == 43 and len(validation) == 1020
    rows = [{'cohort': 'discovery', 'CID': cid, 'field': '', 'zHEL': '', 'in_public_K21_DES_input': cid in des} for cid in discovery]
    rows += [{'cohort': 'validation', 'CID': r['CID'], 'field': r['field'], 'zHEL': r['zHEL'], 'in_public_K21_DES_input': r['CID'] in des} for r in validation]
    assert len({r['CID'] for r in rows}) == 1063
    with (OUT/'membership.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    counts = {cohort: {str(member): sum(r['cohort'] == cohort and r['in_public_K21_DES_input'] == member for r in rows)
                       for member in (True, False)} for cohort in ('discovery', 'validation')}
    save(OUT/'design.json', {'counts': counts, 'training_DES_rows': len(des), 'joins': joins,
                            'max_sky_separation_arcsec': max(x['separation_arcsec'] for x in joins),
                            'input_sha256': {str(p.relative_to(ROOT)): sha(p) for p in paths},
                            'membership_sha256': sha(OUT/'membership.csv'),
                            'scope': 'Public config-listed K21 DES inputs; no exact DES5YR accepted training execution or all-survey physical alias completeness claimed.'})
    print(json.dumps(counts, indent=2))


def score():
    frozen = json.loads((OUT/'design.json').read_text())
    for path, digest in frozen['input_sha256'].items():
        assert sha(ROOT/path) == digest
    assert sha(OUT/'membership.csv') == frozen['membership_sha256']
    members = [r for r in csv.DictReader((OUT/'membership.csv').open()) if r['cohort'] == 'validation']
    scores_path = BASE/'validation1020/analysis/object-scores.csv'
    score_rows = [r for r in csv.DictReader(scores_path.open()) if r['arm'] == 'published_mask']
    scores = {r['CID']: r for r in score_rows}
    assert len(scores) == len(score_rows) == 1020 and set(scores) == {r['CID'] for r in members}
    def summary(rows):
        values = [scores[r['CID']] for r in rows]
        m = sum(float(r['matched_filter']) for r in values)
        info = sum(float(r['information']) for r in values)
        gain = sum(float(r['fixed_prediction_gain']) for r in values)
        assert abs(gain-m+.5*info) < 1e-10
        return {'objects': len(rows), 'M': m, 'I': info, 'G': gain, 'amplitude': m/info if info > 0 else None}
    result = {}
    field_rows = []
    for member in ('True', 'False'):
        subset = [r for r in members if r['in_public_K21_DES_input'] == member]
        result[member] = summary(subset)
        for field in sorted({r['field'] for r in subset}):
            field_rows.append({'in_public_K21_DES_input': member, 'field': field,
                               **summary([r for r in subset if r['field'] == field])})
    all_result = summary(members)
    for key in ('objects', 'M', 'I', 'G'):
        assert abs(result['True'][key]+result['False'][key]-all_result[key]) < 1e-10
    save(OUT/'score.json', {'strata': result, 'all': all_result, 'fields': field_rows,
                           'input_sha256': {str(p.relative_to(ROOT)): sha(p) for p in
                                            [scores_path, OUT/'design.json', OUT/'membership.csv']},
                           'scope': 'Descriptive partition of original fixed-vector result; neither subset is declared independently trained.'})
    print(json.dumps({'strata': result, 'all': all_result}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('phase', choices=['design', 'score'])
    args = parser.parse_args(); (design if args.phase == 'design' else score)()
