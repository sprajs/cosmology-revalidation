"""Copy the first eight chronological RAW reads without reserializing pixels.

Only explicit FITS header cards are replaced. Original inputs remain immutable.
An already-eight-read input is a whole-file byte identity control.
"""
from pathlib import Path
import hashlib
import json
import struct
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[2]
H = ROOT / 'runs/research_2026_09_26/raisin_hst_pixel_feasibility'
OUT = H / 'calwf3_dark_prefix'


def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def replace_card(block, key, value):
    """Change exactly one existing 80-byte card, preserving its comment."""
    indices = [i for i in range(0, len(block), 80)
               if block[i:i+8].decode('ascii').strip() == key]
    assert len(indices) == 1, key
    start = indices[0]
    original = block[start:start+80].decode('ascii')
    card = fits.Card.fromstring(original)
    if card.value == value:
        return block, None
    replacement = fits.Card(key, value, card.comment).image.encode('ascii')
    assert len(replacement) == 80
    result = block[:start] + replacement + block[start+80:]
    return result, dict(keyword=key, original=original, replacement=replacement.decode('ascii'))


def prefix(source, target):
    original = source.read_bytes()
    edits, maps, chunks = [], [], []
    with fits.open(source, memmap=True, do_not_scale_image_data=True) as hdus:
        nsamp = int(hdus[0].header['NSAMP'])
        assert nsamp in (8, 16) and len(hdus) == 1 + 5*nsamp
        assert hdus[0].header['NEXTEND'] == 5*nsamp
        assert hdus[0].header['SAMP_SEQ'] == 'SPARS50'
        assert all('CHECKSUM' not in h.header and 'DATASUM' not in h.header for h in hdus)
        first_version = nsamp-7
        primary_end = hdus[1].fileinfo()['hdrLoc']
        primary = original[:primary_end]
        tlast = float(hdus['SCI', first_version].header['SAMPTIME'])
        assert tlast == 302.934753
        if nsamp != 8:
            newcards = dict(NSAMP=8, NEXTEND=40, EXPTIME=tlast,
                            EXPEND=float(hdus[0].header['EXPSTART']) + tlast/86400.)
            for key, value in newcards.items():
                primary, edit = replace_card(primary, key, value)
                if edit:
                    edits.append(dict(hdu=0, **edit))
        chunks.append(primary)
        for newver, oldver in enumerate(range(first_version, nsamp+1), 1):
            assert hdus['SCI', oldver].header['SAMPNUM'] == 8-newver
            for name in ('SCI', 'ERR', 'DQ', 'SAMP', 'TIME'):
                old = hdus[name, oldver]
                info = old.fileinfo()
                header = original[info['hdrLoc']:info['datLoc']]
                data = original[info['datLoc']:info['datLoc']+info['datSpan']]
                header, edit = replace_card(header, 'EXTVER', newver)
                if edit:
                    edits.append(dict(hdu=f'{name},{oldver}', **edit))
                chunks.extend((header, data))
                maps.append(dict(name=name, oldver=oldver, newver=newver,
                                 data_span=info['datSpan'], data_sha256=hashlib.sha256(data).hexdigest()))
    result = b''.join(chunks)
    with target.open('xb') as f:
        f.write(result)
    with fits.open(source, memmap=True, do_not_scale_image_data=True) as original_fits, \
         fits.open(target, memmap=True, do_not_scale_image_data=True) as new:
        assert len(new) == 41 and new[0].header['NSAMP'] == 8 and new[0].header['NEXTEND'] == 40
        assert new[0].header['EXPFLAG'] == original_fits[0].header['EXPFLAG']
        for item in maps:
            a = original_fits[item['name'], item['oldver']]
            b = new[item['name'], item['newver']]
            assert b.header['EXTVER'] == item['newver']
            for key in a.header:
                if key != 'EXTVER':
                    assert str(a.header[key]) == str(b.header[key]), (item, key)
            info = b.fileinfo()
            data = result[info['datLoc']:info['datLoc']+info['datSpan']]
            assert hashlib.sha256(data).hexdigest() == item['data_sha256']
    if nsamp == 8:
        assert result == original and not edits
    return dict(source=str(source), target=str(target), source_sha256=sha(source),
                target_sha256=sha(target), source_nsamp=nsamp, edits=edits, mapping=maps,
                null_whole_file_exact=(result == original),
                scope='Nominal-time first8 prefix; EXPEND derived bookkeeping, not recovered physical clock')


def main():
    OUT.mkdir(exist_ok=False)
    header = json.loads((H/'dark_ramp_execution/header-gate.json').read_text())
    paths = {x['root']: Path(x['path']) for x in header['records']}
    paths['icxoi1bcq'] = H/'dark_ramp_raw_pilot/files/icxoi1bcq_raw.fits'
    roots = ['icxoi1bcq', 'idbx41onq', 'idbx43p7q']
    plan = dict(source_sha256=sha(__file__), roots=roots,
                input_hashes={str(paths[r]): sha(paths[r]) for r in roots},
                construction='Retain chronological reads0..7, exact original data blocks; oldEXTVER9..16 to1..8 for16-read inputs; four primary bookkeeping edits only',
                no_changes='No calibration/reference/quality flag, SAMPNUM/SAMPTIME/DELTATIM or pixel changes',
                identity_gate='Already8-read source whole-file SHA exact and native SCI/ERR/DQ/SAMP/TIME exact against earlier native replay',
                resources='Three RAW copies under60MB; no native scoring in constructor')
    (OUT/'constructor-protocol.json').write_text(json.dumps(plan, indent=2)+'\n')
    results = []
    for root in roots:
        work = OUT/'work'/root
        work.mkdir(parents=True, exist_ok=False)
        results.append(prefix(paths[root], work/(root+'_raw.fits')))
    (OUT/'constructor-result.json').write_text(json.dumps(results, indent=2)+'\n')
    print(json.dumps([dict(root=r, input_nsamp=x['source_nsamp'], edits=len(x['edits']),
                           null_exact=x['null_whole_file_exact']) for r, x in zip(roots, results)]))


if __name__ == '__main__':
    main()
