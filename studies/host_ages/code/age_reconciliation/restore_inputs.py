#!/usr/bin/env python3
"""Restore the six frozen comparison inputs without the original local archive."""
import argparse
import hashlib
import json
from pathlib import Path
import tempfile
import urllib.request

ROOT = Path(__file__).resolve().parents[4]
PP = 'sources/repos/PantheonPlusSH0ES__DataRelease@7fc6805/Pantheon+_Data/4_DISTANCES_AND_COVAR/'
PP_URL = 'https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/7fc680548d1ea9fe6ee983a89d8b5635bb5784a8/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/'


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--destination', type=Path, default=ROOT/'.work/age-inputs')
    args = parser.parse_args()
    manifest = json.loads((ROOT/'studies/host_ages/results/age_reconciliation/manifest.json').read_text())
    routes = {
        f'data/host_ages/chung2025/table{i}.dat': (ROOT/f'data/ages/table{i}.dat', None)
        for i in (1, 2)
    }
    for name in ('Pantheon+SH0ES.dat', 'Pantheon+SH0ES_STAT+SYS.cov'):
        routes[PP+name] = (ROOT/'data/distances'/name, PP_URL+name.replace('+', '%2B'))
    for name in ('table2.dat', 'ReadMe'):
        routes['data/host_ages/gupta2011/'+name] = (None, 'https://cdsarc.cds.unistra.fr/ftp/J/ApJ/740/92/'+name)
    assert set(routes) == set(manifest['inputs'])
    restored = []
    for relative, expected in manifest['inputs'].items():
        target = args.destination/relative
        if target.exists():
            if digest(target) != expected:
                raise ValueError(f'Existing input identity differs: {target}')
            restored.append({'path': relative, 'status': 'already verified'})
            continue
        local, url = routes[relative]
        if local is not None and local.exists():
            payload = local.read_bytes()
            source = str(local.relative_to(ROOT))
        elif url:
            request = urllib.request.Request(url, headers={'User-Agent': 'Cosmology-Revalidation/1.0'})
            with urllib.request.urlopen(request, timeout=45) as response:
                payload = response.read()
            source = url
        else:
            raise FileNotFoundError(f'Restore {local} from the frozen main age supplement first; see docs/methods/data.md')
        if hashlib.sha256(payload).hexdigest() != expected:
            raise ValueError(f'Input hash mismatch: {source}')
        target.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(dir=target.parent, delete=False) as stream:
            temporary = Path(stream.name)
            stream.write(payload)
        try:
            target.hardlink_to(temporary)
        finally:
            temporary.unlink(missing_ok=True)
        restored.append({'path': relative, 'source': source, 'sha256': expected})
    print(json.dumps({'archive_argument': str(args.destination), 'verified_inputs': restored}, indent=2))


if __name__ == '__main__':
    main()
