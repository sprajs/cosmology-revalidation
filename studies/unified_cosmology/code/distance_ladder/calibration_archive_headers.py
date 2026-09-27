"""Sequential streaming supplement for the five previously deferred tarballs.

No files are extracted and no inputs are changed. Compressed bytes are hashed
as they are read; complete gzip streams are drained to validate their trailers.
Only candidate text members up to2MiB are inspected in memory. No nested archive
decompression or inference is performed.
"""
import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import tarfile
import time

ROOT = Path(__file__).resolve().parents[4]
PARENT = ROOT/'studies/unified_cosmology/results/distance_ladder/calibration-input-search.json'
PARENT_SHA = '490f75f9407c053a63abd24545e663ca351009de454a118675efee09f49a9ec3'
WORK = ROOT/'.work/unified-cosmology/calibration-archive-headers'
EXACT = 'v6_9_duplicate_cid.cov'
MAX_TEXT = 2*1024**2
TARGET = re.compile(r'v6[_\.]?9|duplicate.?cid|ceph.*cov|cov.*ceph', re.I)
RECIPE = re.compile(r'(ceph|shoes|pantheon|pplus|duplicate|covariance).*\.(py|ipynb|ya?ml|md|txt)$', re.I)
TEXT_NEEDLES = re.compile(r'v6[_\.]?9_duplicate_cid|CEPH_DIST|ceph[^\n]{0,80}cov|cov[^\n]{0,80}ceph', re.I)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


class HashingReader:
    def __init__(self, stream):
        self.stream = stream
        self.hasher = hashlib.sha256()
        self.bytes_read = 0

    def read(self, size=-1):
        value = self.stream.read(size)
        self.hasher.update(value)
        self.bytes_read += len(value)
        return value


def resolve(label, archive):
    if label.startswith('$ARCHIVE/'):
        path = archive/label[len('$ARCHIVE/'):]
        assert path.resolve().is_relative_to(archive.resolve())
        return path
    assert label.startswith('$REPO/')
    path = ROOT/label[len('$REPO/'):]
    assert path.resolve().is_relative_to(ROOT)
    return path


def scan(record, archive, index):
    path = resolve(record['path'], archive)
    started = time.monotonic()
    out = dict(path=record['path'], expected_compressed_bytes=record['bytes'],
               member_count=0, exact_filename_matches=[], candidate_members=[],
               nested_archive_count=0, errors=[], gzip_trailer_validated=False,
               tar_headers_complete=False, extracted_files=0)
    ledger = WORK/f'{index:02d}-members.jsonl'
    try:
        before = path.stat()
        assert before.st_size == record['bytes'], 'Input size changed since the parent filename audit.'
        with path.open('rb') as original, ledger.open('w') as listing:
            reader = HashingReader(original)
            try:
                with gzip.GzipFile(fileobj=reader, mode='rb') as uncompressed:
                    with tarfile.open(fileobj=uncompressed, mode='r|', bufsize=1024*1024) as bundle:
                        last_progress = time.monotonic()
                        for member in bundle:
                            name = PurePosixPath(member.name).name
                            header = dict(name=member.name, size=member.size,
                                          type=member.type.decode('ascii', errors='replace'),
                                          linkname=member.linkname)
                            listing.write(json.dumps(header, separators=(',', ':'))+'\n')
                            out['member_count'] += 1
                            if re.search(r'\.(tar|tgz|tar\.gz|zip)$', name, re.I):
                                out['nested_archive_count'] += 1
                            exact = name == EXACT
                            target = bool(TARGET.search(name))
                            recipe = bool(RECIPE.search(name)) or name.upper() == 'PPLUS.YML'
                            host_matrix = name.lower().startswith('allc_shoes_ceph')
                            if exact:
                                out['exact_filename_matches'].append(header)
                            if exact or target or recipe or host_matrix:
                                candidate = dict(header, exact_name=exact, targeted_name=target,
                                                 recipe_name=recipe, known_host_matrix_name=host_matrix)
                                text_extension = PurePosixPath(name).suffix.lower() in {'.py','.ipynb','.yml','.yaml','.md','.txt','.cov'}
                                if member.isfile() and text_extension and member.size <= MAX_TEXT:
                                    stream = bundle.extractfile(member)
                                    body = stream.read(MAX_TEXT+1)
                                    assert len(body) == member.size
                                    candidate['in_memory_content_sha256'] = hashlib.sha256(body).hexdigest()
                                    try:
                                        text = body.decode('utf-8')
                                    except UnicodeDecodeError:
                                        candidate['content_inspection'] = 'not_utf8_text'
                                    else:
                                        hits = [dict(line=i, text=line[:700]) for i,line in enumerate(text.splitlines(), 1)
                                                if TEXT_NEEDLES.search(line)]
                                        candidate.update(content_inspection='text_inspected_in_memory_no_extraction',
                                                         matching_line_count=len(hits), matching_lines=hits[:30],
                                                         excerpt_limit_reached=len(hits)>30)
                                else:
                                    candidate['content_inspection'] = 'not_a_small_regular_text_member'
                                out['candidate_members'].append(candidate)
                            if time.monotonic()-last_progress > 30:
                                print(json.dumps(dict(archive=index, members=out['member_count'],
                                                      compressed_bytes_read=reader.bytes_read,
                                                      seconds=time.monotonic()-started)), flush=True)
                                last_progress = time.monotonic()
                    out['tar_headers_complete'] = True
                    # Tar ends at its terminator; continue through gzip CRC/ISIZE
                    # and any remaining compressed bytes instead of hashing only
                    # the prefix consumed by tarfile's read-ahead buffer.
                    while uncompressed.read(1024*1024):
                        pass
                    out['gzip_trailer_validated'] = True
            except Exception as error:
                out['errors'].append({'stage': 'tar_or_gzip_stream', 'error': repr(error)})
            finally:
                while reader.read(1024*1024):
                    pass
                out['compressed_sha256'] = reader.hasher.hexdigest()
                out['compressed_bytes_hashed'] = reader.bytes_read
        after = path.stat()
        unchanged = (before.st_size, before.st_mtime_ns, before.st_ino) == (after.st_size, after.st_mtime_ns, after.st_ino)
        out['input_stat_unchanged'] = unchanged
        assert unchanged, 'Archive changed while being inspected.'
        assert out['compressed_bytes_hashed'] == after.st_size
    except Exception as error:
        out['errors'].append({'stage': 'input_or_identity', 'error': repr(error)})
    if ledger.exists():
        out.update(member_ledger_path=str(ledger.relative_to(ROOT)), member_ledger_sha256=sha(ledger))
    out['status'] = 'complete_read_only_header_audit' if not out['errors'] else 'archive_audit_error_preserved'
    out['seconds'] = time.monotonic()-started
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, default=Path.home()/'.local/share/cosmology-revalidation/archive')
    parser.add_argument('--output', type=Path, default=ROOT/'studies/unified_cosmology/results/distance_ladder/calibration-archive-headers.json')
    args = parser.parse_args()
    assert sha(PARENT) == PARENT_SHA
    parent = json.loads(PARENT.read_text())
    deferred = [r for r in parent['local_search']['tar_headers'] if not r['member_names_checked']]
    assert len(deferred) == 5
    WORK.mkdir(parents=True, exist_ok=True)
    source_sha = sha(__file__)
    report = dict(status='running_sequential_archive_header_supplement',
                  created_utc=datetime.now(timezone.utc).isoformat(),
                  parent_path=str(PARENT.relative_to(ROOT)), parent_sha256=PARENT_SHA,
                  source_sha256={str(Path(__file__).relative_to(ROOT)):source_sha},
                  declared_archive_count=5, declared_total_compressed_bytes=sum(r['bytes'] for r in deferred),
                  max_in_memory_text_member_bytes=MAX_TEXT, archives=[],
                  interpretation='Names and limited candidate source content only; no covariance construction is inferred from a filename match.',
                  limitations=['Nested archive payloads are counted but not recursively decompressed.',
                               'Only candidate regular UTF8 text members up to2MiB receive content inspection.',
                               'This closes the five declared archive-header deferrals; it is not an exhaustive search of private production inputs.'],
                  extraction_operations=0, cosmology_fits=0, network_requests=0, parallel_workers=1)
    write(args.output, report)
    for index, record in enumerate(deferred):
        result = scan(record, args.archive, index)
        report['archives'].append(result)
        write(args.output, report)
        print(json.dumps(dict(archive=index, status=result['status'], members=result['member_count'],
                              exact_matches=len(result['exact_filename_matches']),
                              candidates=len(result['candidate_members']), seconds=result['seconds'])), flush=True)
    assert sha(PARENT) == PARENT_SHA and sha(__file__) == source_sha
    failures = [x for x in report['archives'] if x['errors']]
    report['status'] = 'complete_archive_header_supplement' if not failures else 'archive_header_supplement_with_preserved_errors'
    report['all_exact_matches'] = [{'archive': x['path'], 'member': m} for x in report['archives'] for m in x['exact_filename_matches']]
    report['failed_archives'] = len(failures)
    report['total_members'] = sum(x['member_count'] for x in report['archives'])
    report['completed_utc'] = datetime.now(timezone.utc).isoformat()
    write(args.output, report)
    print(json.dumps({k:report[k] for k in ['status','total_members','all_exact_matches','failed_archives']}), flush=True)


if __name__ == '__main__':
    main()
