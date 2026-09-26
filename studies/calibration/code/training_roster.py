"""Read the exact config-driven K21 input lists without executing training code."""
import configparser
import csv
import hashlib
import json
import tarfile
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'sources/updates/2026-09-26-training-roster'
PREFIX = 'SALT3TRAIN_K21_PUBLIC/'


def main():
    rows = []
    inputs = {}
    def record(name, data):
        inputs[name] = hashlib.sha256(data).hexdigest()
        return data.decode('utf-8')
    # Stream once; random seeks through gzip otherwise re-decompress the archive
    # for every input file. Only source text/configuration is retained in memory.
    payloads = {}
    with tarfile.open(OUT/'SALT3TRAIN_K21_PUBLIC.tgz', mode='r|gz') as stream:
        for member in stream:
            if member.isfile() and (member.name.endswith(('.LIST', '.conf')) or
                                     '/SALT3TRAIN_K21/SALT3TRAIN_K21_' in member.name):
                assert '..' not in PurePosixPath(member.name).parts
                payloads[member.name] = stream.extractfile(member).read()
    with tarfile.open(OUT/'SALT3TRAIN_K21_PUBLIC.tgz') as archive:
        def get(name):
            member = archive.getmember(name)
            assert member.isfile() and '..' not in PurePosixPath(name).parts
            return record(name, payloads[name])
        conf = get(PREFIX+'Train_SALT3_public.conf')
        config = configparser.ConfigParser()
        config.read_string(conf)
        (OUT/'Train_SALT3_public.conf').write_text(conf)
        (OUT/'training.conf').write_text(get(PREFIX+'training.conf'))
        for snlist in config['iodata']['snlists'].split(','):
            name = PREFIX+snlist.strip()
            listing = get(name)
            for line in listing.splitlines():
                line = line.split('#')[0].strip()
                if not line:
                    continue
                data_path = str(PurePosixPath(name).parent / line)
                text = get(data_path)
                header = {}
                for raw in text.splitlines():
                    if raw.strip().startswith(('OBS:', 'SPECTRUM_ID:', 'SPEC:')):
                        break
                    if ':' in raw and not raw.lstrip().startswith('#'):
                        key, value = raw.split(':', 1)
                        header[key.strip()] = value.split('#', 1)[0].strip()
                rows.append({'archive_path': data_path, 'input_list': snlist.strip(),
                             'SNID': header.get('SNID', ''), 'SURVEY': header.get('SURVEY', ''),
                             'RA': header.get('RA', ''), 'DEC': header.get('DEC', ''),
                             'PEAKMJD': header.get('PEAKMJD', ''),
                             'REDSHIFT_HELIO': header.get('REDSHIFT_HELIO', ''),
                             'IAUC': header.get('IAUC', ''), 'FILTERS': header.get('FILTERS', ''),
                             'NOBS': header.get('NOBS', ''), 'sha256': inputs[data_path]})
        for name in ('SALT3_PARS_INIT.LIST', 'SALT3_PKMJD_INIT.LIST'):
            (OUT/name).write_text(get(PREFIX+name))
    assert all(r['SNID'] for r in rows)
    assert len({r['archive_path'] for r in rows}) == len(rows)
    with (OUT/'roster.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    counts = {survey: sum(r['SURVEY'] == survey for r in rows) for survey in sorted({r['SURVEY'] for r in rows})}
    result = {'files': len(rows), 'unique_string_SNID': len({r['SNID'] for r in rows}),
              'survey_counts': counts, 'input_sha256': inputs,
              'archive_sha256': hashlib.sha256((OUT/'SALT3TRAIN_K21_PUBLIC.tgz').read_bytes()).hexdigest(),
              'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'roster_sha256': hashlib.sha256((OUT/'roster.csv').read_bytes()).hexdigest(),
              'scope': 'Config-listed public K21 input files, not proof of accepted DES5YR training execution or unique physical supernovae.'}
    (OUT/'roster-manifest.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ('files', 'unique_string_SNID', 'survey_counts')}, indent=2))


if __name__ == '__main__':
    main()
