"""Pinned, isolated two-RAW replay; compare all output pixels without selection."""
from pathlib import Path
import hashlib
import json
import os
import shutil
import subprocess
import time
import warnings
import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[2]
H = ROOT/'runs/research_2026_09_26/raisin_hst_pixel_feasibility'
OUT = H/'calwf3_native_replay'


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()


def save(name, value):
    with (OUT/name).open('x') as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write('\n')


def main():
    invpath = H/'calwf3_native_feasibility/inventory.json'
    inv = json.loads(invpath.read_text())
    acquire = json.loads((OUT/'acquisition-result.json').read_text())
    fetched = {r['name']: r for r in acquire['rows']}
    binary = OUT/'build-fortran/pkg/wfc3/calwf3.e'
    assert sha(binary) == '230dc5e2c760217e7924f0c11153278dd43e10c22fa73828cad8aff0a9c97d00'
    version = subprocess.check_output([str(binary), '--version'], text=True).strip()
    assert version == '3.7.3'
    refs = OUT/'iref'
    refs.mkdir(exist_ok=False)
    hashes = {str(invpath): sha(invpath), str(binary): sha(binary),
              str(Path(__file__)): sha(Path(__file__))}
    for path in (OUT/'build-fortran').rglob('*.so'):
        hashes[str(path)] = sha(path)
    for name, r in inv['refs'].items():
        if r['local']:
            source = Path(r['local'][0]['path'])
            expected = r['local'][0]['sha256']
        else:
            source = OUT/'files'/name
            expected = fetched[name]['sha256']
        assert sha(source) == expected
        hashes[str(source)] = expected
        (refs/name).symlink_to(source)
    for obj in inv['objects'].values():
        for key in ['raw', 'flt']:
            path = Path(obj[key]['path'])
            assert sha(path) == obj[key]['sha256']
            hashes[str(path)] = obj[key]['sha256']
    protocol = dict(purpose='Exact current archive RAW-to-FLT numerical closure, no historical reduction or SN photometry claim',
                    roots=list(inv['objects']), binary_version=version,
                    input_hashes=hashes, environment={'iref': str(refs)+'/', 'OMP_NUM_THREADS': '1'},
                    per_native_process_seconds=60, native_total_seconds=120,
                    worker_count=1, new_output_bytes_cap=600000000,
                    primary='All SCI,ERR,DQ,SAMP,TIME pixels bit-identical, dimensions/units/reference names/processing switches equal',
                    diagnostic='Report all float mismatch counts, max absolute difference, RMS and maximum ULP distance. <=8 float32 ULP is separately labelled approximate closure, never exact.',
                    DQ_SAMP_TIME='Exact required; all pixels retained including flagged pixels and nonfinite patterns',
                    failure='Keep both planned roots and all outputs/failures; no threshold change or source patch after scoring')
    save('replay-protocol.json', protocol)
    runs = []
    for root, obj in inv['objects'].items():
        work = OUT/'work'/root
        work.mkdir(parents=True, exist_ok=False)
        raw = Path(obj['raw']['path'])
        shutil.copyfile(raw, work/raw.name)
        assert sha(work/raw.name) == obj['raw']['sha256']
        env = os.environ.copy()
        env.update(protocol['environment'])
        command = [str(binary), raw.name]
        started = time.monotonic()
        with (work/'native.log').open('x') as log:
            try:
                result = subprocess.run(command, cwd=work, env=env, stdout=log,
                                        stderr=subprocess.STDOUT, timeout=60)
                item = dict(root=root, returncode=result.returncode, command=command,
                            seconds=time.monotonic()-started)
            except subprocess.TimeoutExpired:
                item = dict(root=root, returncode=None, timeout=True, command=command,
                            seconds=time.monotonic()-started)
        item['outputs'] = {p.name: dict(bytes=p.stat().st_size, sha256=sha(p))
                           for p in work.iterdir() if p.is_file()}
        runs.append(item)
        save(root+'-execution.json', item)
        assert sum(p.stat().st_size for p in (OUT/'work').rglob('*') if p.is_file()) <= 600000000
        print(json.dumps({k: item[k] for k in ['root','returncode','seconds']}), flush=True)
    comparisons = []
    for item in runs:
        root = item['root']
        produced = OUT/'work'/root/(root+'_flt.fits')
        if item['returncode'] != 0 or not produced.exists():
            comparisons.append(dict(root=root, primary_exact=False, execution_failed=True))
            continue
        archived = Path(inv['objects'][root]['flt']['path'])
        metrics = []
        with warnings.catch_warnings(record=True) as warning_list:
            warnings.simplefilter('always')
            with fits.open(produced, memmap=False) as a, fits.open(archived, memmap=False) as b:
                headers = {}
                for key in ['CAL_VER', *inv['objects'][root]['raw_refs'],
                            *inv['objects'][root]['archive_flt_switches']]:
                    headers[key] = dict(replay=a[0].header.get(key), archive=b[0].header.get(key))
                for name in ['SCI','ERR','DQ','SAMP','TIME']:
                    x = a[name, 1].data
                    y = b[name, 1].data
                    assert x.shape == y.shape == (1014, 1014)
                    assert x.dtype == y.dtype
                    unit_equal = a[name, 1].header.get('BUNIT') == b[name, 1].header.get('BUNIT')
                    samebits = x.tobytes() == y.tobytes()
                    record = dict(HDU=name, pixels=x.size, dtype=str(x.dtype),
                                  unit_equal=unit_equal, bit_exact=samebits,
                                  unequal_values=int(np.sum(x != y)))
                    if x.dtype.kind == 'f':
                        assert x.dtype.itemsize == 4
                        finite = np.isfinite(x) & np.isfinite(y)
                        record['nonfinite_pattern_equal'] = bool(np.array_equal(np.isfinite(x), np.isfinite(y)))
                        gap = x[finite].astype(np.float64)-y[finite].astype(np.float64)
                        record.update(max_absolute_difference=float(np.max(abs(gap))),
                                      RMS_difference=float(np.sqrt(np.mean(gap*gap))))
                        # Ordered float32 bit keys handle both signs and +/-zero.
                        def keys(z):
                            u = z.astype(np.float32).view(np.uint32)
                            return np.where(u & 0x80000000, ~u, u ^ 0x80000000).astype(np.int64)
                        record['max_ULP_difference'] = int(np.max(abs(keys(x[finite])-keys(y[finite]))))
                    metrics.append(record)
            warning_text = sorted({str(w.message) for w in warning_list})
        header_equal = all(v['replay'] == v['archive'] for v in headers.values())
        exact = header_equal and all(m['unit_equal'] and m['bit_exact'] for m in metrics)
        approximate = header_equal and all(m['unit_equal'] and (
            m['bit_exact'] if m['HDU'] in ['DQ','SAMP','TIME'] else
            m.get('nonfinite_pattern_equal', False) and m['max_ULP_difference'] <= 8)
            for m in metrics)
        comparisons.append(dict(root=root, primary_exact=exact,
                                approximate_8ULP=approximate, headers=headers,
                                metrics=metrics, warnings=warning_text))
    save('replay-result.json', dict(primary_exact=all(x['primary_exact'] for x in comparisons),
         native_seconds=sum(x['seconds'] for x in runs), comparisons=comparisons,
         scope=protocol['purpose']))
    print(json.dumps([dict(root=x['root'], exact=x['primary_exact'],
                          approximate=x.get('approximate_8ULP')) for x in comparisons]))


if __name__ == '__main__':
    main()
