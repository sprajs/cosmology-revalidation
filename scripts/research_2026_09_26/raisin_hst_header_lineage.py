"""Inspect acquired PRIMARY headers only; never reads image pixels."""
from pathlib import Path
import csv
import hashlib
import json
import re
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[2]
P = ROOT / 'runs/research_2026_09_26/raisin_hst_pixel_feasibility/header_followup'
visits = {r['obs_id']: r for r in csv.DictReader((P.parent / 'visit-groups.csv').open())}
groups, exposures = [], []
for path in sorted(P.rglob('*.primary-header')):
    raw = path.read_bytes()
    h = fits.Header.fromstring(raw.decode('ascii'))
    assert h['NAXIS'] == 0 and len(raw) % 2880 == 0
    root = h['ROOTNAME']
    v = visits[root]
    assert h['FILTER'] == v['filter'] and h['TARGNAME'].upper() == 'DES16E1DCX'
    assert h['ASN_ID'].lower() == root
    inputs = []
    for key in h:
        if re.fullmatch(r'D\d{3}DATA', key):
            prefix = key[:4]
            name = h[key]
            assert name.endswith('_flt.fits[sci,1]')
            inputs.append(name)
            exposures.append(dict(association=root, filter=h['FILTER'], visit=v['visit'],
                constituent=name, exposure_seconds=h[prefix+'DEXP'],
                drizzle_kernel=h[prefix+'KERN'], pixfrac=h[prefix+'PIXF'],
                output_scale=h[prefix+'SCAL'], input_scale=h[prefix+'ISCL'],
                output_units=h[prefix+'OUUN']))
    assert len(inputs) == len(set(inputs)) == h['NDRIZIM'] == 4
    assert abs(sum(x['exposure_seconds'] for x in exposures if x['association'] == root) - h['EXPTIME']) < 1e-6
    groups.append(dict(association=root, filter=h['FILTER'], visit=v['visit'],
        n_constituents=len(inputs), exposure_seconds=h['EXPTIME'],
        processing_date=h['DATE'], pipeline=h['OPUS_VER'], calwf3=h['CAL_VER'],
        crds_context=h['CRDS_CTX'], photflam=h['PHOTFLAM'], photplam=h['PHOTPLAM'],
        photometry_reference=h['IMPHTTAB'], flat_reference=h['PFLTFILE'],
        nonlinearity_reference=h['NLINFILE'], nlincorr=h['NLINCORR'],
        header_sha256=hashlib.sha256(raw).hexdigest(), header_path=str(path.relative_to(ROOT))))
assert len(groups) == 8 and len(exposures) == len({x['constituent'] for x in exposures}) == 32
for name, rows in [('association-processing.csv', groups), ('constituent-exposures.csv', exposures)]:
    with (P / name).open('w') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
summary = dict(associations=len(groups), unique_constituents=len(exposures),
    verified='All32 constituent exposure names,4 per association, and exposure-time sums from current primary headers',
    header_bytes=sum(p.stat().st_size for p in P.rglob('*.primary-header')),
    pixels_read=0, current_processing_dates=sorted({x['processing_date'] for x in groups}),
    calwf3_versions=sorted({x['calwf3'] for x in groups}),
    drizzle_kernels=sorted({x['drizzle_kernel'] for x in exposures}),
    pixfrac=sorted({x['pixfrac'] for x in exposures}),
    output_scale_arcsec=sorted({x['output_scale'] for x in exposures}),
    limitations=['Headers are from2026 reprocessing; original published images not reproduced.',
        'NLINCORR records accumulated-count correction, not proof of count-rate nonlinearity correction.',
        'Constituent names identified, but their actual files/calibrations have not yet been acquired.',
        'Public currentDRZ scale differs from paper0.11arcsec; no equivalence with private subtraction products assumed.'])
(P/'lineage-result.json').write_text(json.dumps(summary,indent=2)+'\n')
manifest = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
    for p in sorted(P.rglob('*')) if p.is_file() and p.name != 'lineage-manifest.json'}
manifest[str(Path(__file__).resolve().relative_to(ROOT))] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
(P/'lineage-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps(summary,indent=2))
