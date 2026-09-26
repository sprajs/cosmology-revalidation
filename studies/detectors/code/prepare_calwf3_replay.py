"""Acquire exact public replay dependencies; no calibration or SCI scoring."""
from pathlib import Path
import hashlib
import io
import json
import tarfile
import time
import urllib.request
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[2]
H = ROOT/'runs/research_2026_09_26/raisin_hst_pixel_feasibility'
OUT = H/'calwf3_native_replay'


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()


def save(name, value):
    with (OUT/name).open('x') as f:
        json.dump(value, f, indent=2)
        f.write('\n')


def main():
    OUT.mkdir(exist_ok=False)
    invpath = H/'calwf3_native_feasibility/inventory.json'
    treepath = H/'calwf3_variance_source/tree-3.7.3.json'
    inv = json.loads(invpath.read_text())
    commit = '6a1147d7bdccb7e2a7b73af276f4597c6420fbb8'
    source = dict(name='hstcal-pinned.tar.gz',
                  url='https://api.github.com/repos/spacetelescope/hstcal/tarball/'+commit,
                  size=None, cap=20000000)
    cmake = dict(name='cmake-3.31.6-py3-none-manylinux_2_17_x86_64.manylinux2014_x86_64.whl',
                 url='https://files.pythonhosted.org/packages/59/e8/096984b89133681533650b9078c5ed1c5c9b534e869b5487f22d4de1935c/cmake-3.31.6-py3-none-manylinux_2_17_x86_64.manylinux2014_x86_64.whl',
                 size=27800904, cap=30000000,
                 sha256='1c8b05df0602365da91ee6a3336fe57525b137706c4ab5675498f662ae1dbcec')
    refs = [dict(name=n, url=r['official_head']['url'],
                 size=int(r['official_head']['content_length']),
                 cap=int(r['official_head']['content_length']))
            for n, r in inv['refs'].items() if not r['local']]
    items = [source, cmake]+refs
    protocol = dict(purpose='Private pinned CALWF3 replay dependencies only; no image outcome',
                    source_commit=commit, inventory_sha256=sha(invpath),
                    tree_sha256=sha(treepath), items=items, cap_bytes=800000000,
                    acquisition_seconds=240, request_timeout_seconds=60,
                    one_attempt_per_file=True, existing_files_untouched=True,
                    source_check='Every regular extracted file must match exact Git blob; no links',
                    reference_check='Exact HTTP size, file SHA256, FITS structure and header lineage',
                    build='Separate future private two-worker build capped180seconds; no install into existing env/system',
                    executed_source_sha256=sha(Path(__file__)))
    save('acquisition-protocol.json', protocol)
    (OUT/'files').mkdir()
    started = time.monotonic()
    downloaded = 0
    result = []
    try:
        for item in items:
            assert time.monotonic()-started < protocol['acquisition_seconds']
            target = OUT/'files'/item['name']
            n = 0
            t = time.monotonic()
            record = dict(name=item['name'], url=item['url'])
            result.append(record)
            with urllib.request.urlopen(item['url'], timeout=60) as response:
                record.update(status=response.status, final_url=response.url,
                              headers=dict(response.headers))
                assert response.status == 200
                with target.open('xb') as f:
                    while True:
                        b = response.read(1048576)
                        if not b:
                            break
                        n += len(b)
                        downloaded += len(b)
                        assert n <= item['cap'] and downloaded <= protocol['cap_bytes']
                        assert time.monotonic()-started < protocol['acquisition_seconds']
                        f.write(b)
            record.update(bytes=n, sha256=sha(target), seconds=time.monotonic()-t)
            assert item['size'] is None or n == item['size']
            assert 'sha256' not in item or record['sha256'] == item['sha256']
            if item['name'].endswith('.fits'):
                with fits.open(target, memmap=False) as hd:
                    hd.verify('exception')
                    record['HDU_count'] = len(hd)
                    record['primary'] = {k: hd[0].header.get(k) for k in
                                         ['TELESCOP','INSTRUME','DETECTOR','FILETYPE','USEAFTER']}
            print(json.dumps({k: record[k] for k in ['name','bytes','sha256','seconds']}), flush=True)
        tree = json.loads(treepath.read_text())
        wanted = {x['path']: x for x in tree['tree'] if x['type'] == 'blob'}
        dest = OUT/'source'
        dest.mkdir()
        found = set()
        with tarfile.open(OUT/'files'/source['name']) as tar:
            for member in tar:
                if member.isdir():
                    continue
                assert member.isfile(), member.name
                name = '/'.join(Path(member.name).parts[1:])
                assert name in wanted and name not in found, name
                blob = tar.extractfile(member).read()
                digest = hashlib.sha1(b'blob '+str(len(blob)).encode()+b'\0'+blob).hexdigest()
                assert digest == wanted[name]['sha'], name
                path = dest/name
                path.parent.mkdir(exist_ok=True, parents=True)
                path.write_bytes(blob)
                path.chmod(0o755 if wanted[name]['mode']=='100755' else 0o644)
                found.add(name)
        assert found == set(wanted)
        save('acquisition-result.json', dict(pass_all=True, network_bytes=downloaded,
              seconds=time.monotonic()-started, regular_source_files=len(found), rows=result))
    except Exception as exc:
        save('acquisition-failure.json', dict(exception=type(exc).__name__, message=str(exc),
             bytes=downloaded, seconds=time.monotonic()-started, rows=result))
        raise


if __name__ == '__main__':
    main()
