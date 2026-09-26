"""Independent RAW-prefix metadata audit; never accesses an HDU data array."""
from pathlib import Path
import argparse
import hashlib
import json
from astropy.io import fits


def inspect(path):
    with fits.open(path, memmap=True, do_not_scale_image_data=True) as hdus:
        primary = dict(hdus[0].header)
        headers = {(h.header['EXTNAME'], int(h.header['EXTVER'])): h.header.copy()
                   for h in hdus[1:]}
        assert len(headers) == primary['NEXTEND']
        n = primary['NSAMP']
        assert len(headers) == 5*n
        groups = {}
        for i in range(1, n+1):
            sci = headers['SCI', i]
            row = {k: sci.get(k) for k in ['SAMPNUM', 'SAMPTIME', 'DELTATIM',
                'NAXIS1', 'NAXIS2', 'BITPIX', 'BUNIT', 'LTV1', 'LTV2',
                'LTM1_1', 'LTM2_2']}
            assert row['NAXIS1'] == row['NAXIS2'] == 1024
            for ext in ['TIME', 'SAMP']:
                h = headers[ext, i]
                assert h['NAXIS'] == 0 and 'PIXVALUE' in h
                row[ext] = {k: h.get(k) for k in
                            ['PIXVALUE', 'NPIX1', 'NPIX2', 'BITPIX', 'BUNIT']}
            assert row['SAMPNUM'] not in groups
            groups[row['SAMPNUM']] = row
        assert sorted(groups) == list(range(n))
        return primary, groups, {
            ext: [headers[ext, i]['NAXIS'] for i in range(1, n+1)]
            for ext in ['SCI', 'ERR', 'DQ', 'SAMP', 'TIME']}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('pilot', type=Path)
    args = ap.parse_args()
    out = args.pilot/'root-prefix-review'
    out.mkdir(exist_ok=False)
    results = []
    paths = []
    for visit, dark, science in [('search', 'idbx43p7q', 'icxoi1bcq'),
                                  ('template', 'idp247tnq', 'icxoi4hgq')]:
        dp, sp = [args.pilot/'files'/(root+'_raw.fits') for root in [dark, science]]
        paths += [dp, sp]
        d, dg, dl = inspect(dp)
        s, sg, sl = inspect(sp)
        assert d['NSAMP'] == 16 and s['NSAMP'] == 8
        assert {i: dg[i] for i in range(8)} == sg
        for k in ['DETECTOR', 'SUBARRAY', 'SUBTYPE', 'SAMP_SEQ', 'SAMPZERO',
                  'CCDGAIN', 'CCDTAB', 'OSCNTAB', 'NLINFILE', 'BLEVCORR', 'ZOFFCORR']:
            assert d[k] == s[k], k
        assert d['SAMP_SEQ'] == 'SPARS50'
        assert not d['SUBARRAY']
        primary_keys = ['ROOTNAME', 'NSAMP', 'EXPFLAG', 'QUALITY', 'QUALCOM1',
            'QUALCOM2', 'FILTER', 'SAA_EXIT', 'SAA_TIME', 'CCDGAIN', 'CCDTAB',
            'OSCNTAB', 'BPIXTAB', 'NLINFILE', 'ZSIGCORR', 'DARKCORR', 'FLATCORR',
            'BLEVCORR', 'ZOFFCORR', 'UNITCORR', 'CRCORR', 'CAL_VER']
        differences = {k: {'dark': d.get(k), 'science': s.get(k)}
                       for k in primary_keys if d.get(k) != s.get(k)}
        results.append(dict(visit=visit, dark=dark, science=science,
            metadata_prefix_exact=True, common_CCDTAB=d['CCDTAB'],
            common_OSCNTAB=d['OSCNTAB'], common_NLINFILE=d['NLINFILE'],
            prefix=[dg[i] for i in range(8)], primary_differences=differences,
            dark_HDU_layout=dl, science_HDU_layout=sl))
    result = dict(header_gate_pass=True, pixels_scored=False, pairs=results,
        scope='The eight-read timing/geometry prefix matches exactly. '
        'Calibration flags, quality state and bad-pixel references differ. '
        'No gain table values, read covariance, raw-to-calibrated replay, '
        'master-dark independence or illuminated-science transport is certified.')
    (out/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    paths += [Path(__file__), out/'result.json']
    (out/'manifest.json').write_text(json.dumps({str(p.resolve()):
        hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}, indent=2)+'\n')
    print(json.dumps({'header_gate_pass': True, 'pairs': len(results),
                      'pixels_scored': False}, indent=2))


if __name__ == '__main__':
    main()
