"""Independent, read-only verification of frozen RAISIN refit input arms."""
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'runs/research_2026_09_26/astra_design/raisin_signed_refit'
OUT = ROOT / 'runs/research_2026_09_26/raisin_refit_input_review'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_phot(path):
    header, observations = {}, []
    names = None
    for line in path.read_text().splitlines():
        words = line.split()
        if not words:
            continue
        if words[0] == 'OBS:':
            assert names is not None and len(words[1:]) == len(names)
            row = dict(zip(names, words[1:]))
            observations.append((row, line))
        elif words[0] == 'VARLIST:':
            names = words[1:]
            header['VARLIST:'] = tuple(names)
        elif words[0].endswith(':') and words[0] not in ('NOBS:', 'END:'):
            header[words[0]] = tuple(words[1:])
        elif words[0] == 'NOBS:':
            count = int(words[1])
    assert count == len(observations), path
    return header, observations


def key(row):
    return (row.get('FLT', row.get('BAND')), float(row['MJD']),
            float(row['FLUXCAL']), float(row['FLUXCALERR']))


def main():
    OUT.mkdir(exist_ok=True)
    assert not (OUT / 'input-check.json').exists(), 'Preserve completed review'
    manifest_path = SOURCE / 'input-manifest.json'
    manifest = json.loads(manifest_path.read_text())
    verified = {}
    for name, expected in manifest['inputs_sha256'].items():
        actual = digest(ROOT / name)
        assert actual == expected, name
        verified[name] = actual
    cases = list(csv.DictReader((SOURCE / 'cohort.csv').open()))
    assert sorted(c['CID'] for c in cases) == manifest['cohort']
    results = []
    for case in cases:
        cid = case['CID']
        released = ROOT / case['released_path']
        raw = ROOT / case['raw_path']
        header, rrows = read_phot(released)
        _, rawrows = read_phot(raw)
        optical = [r for r, _ in rawrows if key(r)[0] in ('g', 'r', 'i', 'z')]
        positive = [r for r in optical if float(r['FLUXCAL']) > 0]
        negative = [r for r in optical if float(r['FLUXCAL']) <= 0]
        nonoptical = Counter(key(r) for r, _ in rrows if key(r)[0] not in 'griz')
        arms = {}
        for arm in ('R', 'Rcopy', 'A', 'B', 'H'):
            p = SOURCE / 'data' / ('RSR_' + arm) / (cid + '.snana.dat')
            h, rows = read_phot(p)
            assert h == header, (cid, arm, 'header')
            if arm in ('R', 'Rcopy'):
                assert p.read_bytes() == released.read_bytes()
            assert all(float(r['FLUXCALERR']) > 0 for r, _ in rows)
            assert Counter(key(r) for r, _ in rows if key(r)[0] not in 'griz') == nonoptical
            arms[arm] = rows
        assert Counter(key(r) for r, _ in arms['A'] if key(r)[0] in 'griz') == Counter(map(key, positive))
        assert Counter(key(r) for r, _ in arms['B'] if key(r)[0] in 'griz') == Counter(map(key, optical))
        assert [line for r, line in arms['B'] if key(r)[0] not in 'griz' or float(r['FLUXCAL']) > 0] == [line for _, line in arms['A']]
        assert [line for _, line in arms['H'][:len(rrows)]] == [line for _, line in rrows]
        assert Counter(key(r) for r, _ in arms['H'][len(rrows):]) == Counter(map(key, negative))
        for row in negative:
            assert not any(key(row)[0] == key(r)[0] and abs(key(row)[1] - key(r)[1]) <= 0.00055 for r, _ in rrows)
        assert Counter(key(r) for r, _ in arms['H']) == Counter(key(r) for r, _ in rrows) + Counter(map(key, negative))
        phase_neg = [r for r in negative if -7 <= (float(r['MJD']) - float(case['peak_header'])) / (1 + float(case['zHEL'])) <= 45]
        results.append({'CID': cid, 'source_optical': len(optical), 'positive': len(positive),
                        'nonpositive': len(negative), 'header_phase_nonpositive': len(phase_neg),
                        'all_input_gates_pass': True})
    output = {'review': 'Independent parser; no import of arm-generating source; no new data mask',
              'manifest_sha256': digest(manifest_path), 'verified_hashes': verified, 'objects': results,
              'total_source_optical': sum(r['source_optical'] for r in results),
              'total_nonpositive': sum(r['nonpositive'] for r in results),
              'total_header_phase_nonpositive': sum(r['header_phase_nonpositive'] for r in results),
              'limitations': ['Input identity alone does not establish native acceptance, numerical convergence or historical reproduction.',
                              'Author PHOTFLAG provenance stays in sidecars; no new flag cut is applied.',
                              'Appended row ordering requires native ordering or permutation closure.']}
    (OUT / 'input-manifest-snapshot.json').write_bytes(manifest_path.read_bytes())
    (OUT / 'executed-source.py').write_bytes(Path(__file__).read_bytes())
    output['review_source_sha256'] = digest(Path(__file__))
    (OUT / 'input-check.json').write_text(json.dumps(output, indent=2) + '\n')
    print(json.dumps({k: output[k] for k in ('total_source_optical', 'total_nonpositive', 'total_header_phase_nonpositive')}))


if __name__ == '__main__':
    main()
