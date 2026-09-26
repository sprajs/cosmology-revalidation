"""Acquire primary CSP DR3 products and both published errata; no data scores."""
from pathlib import Path
import datetime
import hashlib
import json
import urllib.request
import tarfile

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'runs/research_2026_09_26/csp_dr3_provenance'
FILES = {
    'data-page.html': 'https://csp.obs.carnegiescience.edu/data',
    'CSP_Photometry_DR3.tgz': 'https://csp.obs.carnegiescience.edu/data/CSP_Photometry_DR3.tgz',
    'errata-page.html': 'https://authors.library.caltech.edu/records/0f6hh-4pr18',
    'erratum-2017.pdf': 'https://authors.library.caltech.edu/records/0f6hh-4pr18/files/Krisciunas_2017_AJ_154_278.pdf?download=1&research=20260926',
    'erratum-2020.pdf': 'https://authors.library.caltech.edu/records/0f6hh-4pr18/files/Krisciunas_2020_AJ_160_289.pdf?download=1&research=20260926',
}


def main():
    OUT.mkdir(exist_ok=False)
    ledger = []
    for name, url in FILES.items():
        row = {'name': name, 'url': url, 'retrieved_utc': datetime.datetime.now(datetime.timezone.utc).isoformat()}
        try:
            with urllib.request.urlopen(url, timeout=30) as response:
                raw = response.read(30_000_001)
                assert len(raw) <= 30_000_000
                row.update({'content_type': response.headers.get('Content-Type'),
                            'last_modified': response.headers.get('Last-Modified'),
                            'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest()})
            if name.endswith('.pdf'):
                assert raw.startswith(b'%PDF-')
            (OUT / name).write_bytes(raw)
            row['status'] = 'acquired'
        except Exception as exc:
            row.update(status='failed', error=str(exc))
        ledger.append(row)
        (OUT / 'acquisition.json').write_text(json.dumps({'files': ledger}, indent=2)+'\n')
        print(name, row['status'], row.get('bytes', ''))
    archive = OUT / 'CSP_Photometry_DR3.tgz'
    if archive.exists():
        with tarfile.open(archive) as tar:
            members = []
            for m in tar.getmembers():
                entry = {'name': m.name, 'size': m.size, 'regular_file': m.isfile(), 'mtime': m.mtime}
                if m.isfile():
                    entry['sha256'] = hashlib.sha256(tar.extractfile(m).read()).hexdigest()
                members.append(entry)
        (OUT / 'archive-members.json').write_text(json.dumps(members, indent=2)+'\n')
        print('archive entries', len(members))
    (OUT / 'executed-acquisition.py').write_bytes(Path(__file__).read_bytes())


if __name__ == '__main__':
    main()
