"""Execute prefix identity or frozen NORMAL-pair dark calibration, separately."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess
import time
import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[2]
H = ROOT/'runs/research_2026_09_26/raisin_hst_pixel_feasibility'
OUT = H/'calwf3_dark_prefix'
OLD = H/'calwf3_native_replay'


def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def save(name, value):
    with (OUT/name).open('x') as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write('\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', choices=['null', 'dark'])
    args = parser.parse_args()
    binary = OLD/'build-fortran/pkg/wfc3/calwf3.e'
    assert sha(binary) == '230dc5e2c760217e7924f0c11153278dd43e10c22fa73828cad8aff0a9c97d00'
    earlier_protocol = json.loads((OLD/'replay-protocol.json').read_text())
    for path, expected in earlier_protocol['input_hashes'].items():
        assert sha(path) == expected, path
    constructors = json.loads((OUT/'constructor-result.json').read_text())
    for row in constructors:
        assert sha(row['target']) == row['target_sha256']
        assert sha(row['source']) == row['source_sha256']
    if args.mode == 'dark':
        null = json.loads((OUT/'null-result.json').read_text())
        assert null['all_five_components_bit_exact'] and null['whole_RAW_byte_identity']
        release = json.loads((OUT/'dark-release.json').read_text())
        assert release['approved'] is True and release['executor_sha256'] == sha(__file__)
        assert release['constructor_result_sha256'] == sha(OUT/'constructor-result.json')
    roots = ['icxoi1bcq'] if args.mode == 'null' else ['idbx41onq', 'idbx43p7q']
    inputs = {str(Path(__file__)): sha(__file__),
              str(OUT/'constructor-result.json'): sha(OUT/'constructor-result.json'),
              str(OLD/'replay-protocol.json'): sha(OLD/'replay-protocol.json'),
              str(binary): sha(binary)}
    inputs.update({str(OUT/'work'/r/(r+'_raw.fits')): sha(OUT/'work'/r/(r+'_raw.fits')) for r in roots})
    plan = dict(mode=args.mode, roots=roots, input_sha256=inputs,
                refs='Identical pinned iref directory from exact current-science replay; no reference substitutions',
                environment={'iref': str(OLD/'iref')+'/', 'OMP_NUM_THREADS': '1'},
                flags='Exactly input RAW dark flags: ZSIG/DARK/FLAT/PHOT omitted; DQI/ZOFF/BLEV/NLIN/UNIT/CR native',
                clock_scope='Conditional recorded nominal times, not reconstructed physical clocks',
                pair='NORMAL search41onq to43p7q; overlaps original total-RAW cohort',
                next_gate='All native components finite/expected units and shape, fixed reference-mask-only score. Retain native DQ/SAMP/ERR failures and flags, no CR-clipped primary region.',
                per_process_seconds=60, total_native_seconds_cap=120,
                total_work_output_bytes_cap=600000000, workers=1)
    save(args.mode+'-protocol.json', plan)
    runs = []
    for root in roots:
        work = OUT/'work'/root
        assert not (work/(root+'_flt.fits')).exists()
        env = os.environ.copy()
        env.update(plan['environment'])
        t0 = time.monotonic()
        with (work/'native.log').open('x') as log:
            result = subprocess.run([str(binary), root+'_raw.fits'], cwd=work, env=env,
                                    stdout=log, stderr=subprocess.STDOUT, timeout=60)
        row = dict(root=root, returncode=result.returncode, seconds=time.monotonic()-t0)
        row['outputs'] = {str(p): dict(bytes=p.stat().st_size, sha256=sha(p))
                          for p in work.iterdir() if p.is_file()}
        runs.append(row)
        save(root+'-execution.json', row)
        print(json.dumps({k: row[k] for k in ('root','returncode','seconds')}), flush=True)
        assert result.returncode == 0
        assert sum(x['seconds'] for x in runs) <= 120
        assert sum(p.stat().st_size for p in (OUT/'work').rglob('*') if p.is_file()) <= 600000000
    if args.mode == 'null':
        native = OUT/'work/icxoi1bcq/icxoi1bcq_flt.fits'
        reference = OLD/'work/icxoi1bcq/icxoi1bcq_flt.fits'
        compared = []
        with fits.open(native, memmap=False) as a, fits.open(reference, memmap=False) as b:
            for name in ('SCI','ERR','DQ','SAMP','TIME'):
                x, y = a[name,1].data, b[name,1].data
                same = x.shape == y.shape and x.dtype == y.dtype and x.tobytes() == y.tobytes()
                assert a[name,1].header.get('BUNIT') == b[name,1].header.get('BUNIT')
                compared.append(dict(HDU=name, bit_exact=same, pixels=x.size))
        matched = constructors[0]
        record = dict(all_five_components_bit_exact=all(x['bit_exact'] for x in compared),
                      whole_RAW_byte_identity=matched['null_whole_file_exact'], components=compared,
                      native_sha256=sha(native), reference_sha256=sha(reference),
                      reference='Previous native output, not archived post-AstroDrizzle DQ')
        save('null-result.json', record)
        assert record['all_five_components_bit_exact']
        print(json.dumps(record))
    else:
        save('dark-execution-result.json', dict(runs=runs, native_seconds=sum(x['seconds'] for x in runs),
                                               scope='Native execution only, no variance score'))


if __name__ == '__main__':
    main()
