#!/usr/bin/env python3
"""Inspect three fixed exact-byte source metadata documents.

Recover metadata bytes in memory; scientific input restoration and qualification
are separate operations.
"""
import argparse
import base64
import binascii
import hashlib
import json
import os
from pathlib import Path
import re
import stat
import sys
import zlib

ROOT = Path(__file__).absolute().parents[1]
WIRE_LIMIT = 131072
DECODED_LIMIT = 524288
COMPRESSED_LIMIT = 98304
DEPTH_LIMIT = 32
NODE_LIMIT = 16384
GRAPH_SIZE_LIMIT = 8388608
HEX = re.compile(r'[0-9a-f]{64}\Z')
FORMAT = 'zlib-rfc1950-base64-rfc4648/v1'
WIRE_FIELDS = {'schema_version', 'document_id', 'encoding', 'decoded_bytes',
               'decoded_sha256', 'compressed_bytes', 'compressed_sha256', 'origin', 'payload'}
# Exact byte pins; changing a document requires a separately reviewed source identity.
DOCUMENTS = {
    'legacy-inputs': {
        'path': 'sources/legacy-inputs.encoded.json',
        'wire_bytes': 35386, 'wire_sha256': '64cb6730e48b6d9ae1f7fe5fc35dcff6544de941685db6d0c827bf9cc0be05b8',
        'decoded_bytes': 255902, 'decoded_sha256': '7ed81f37e617ce7bb2638f30596165f46ae7ee81ca31788177e90cf57ad5e56f',
        'compressed_bytes': 26197, 'compressed_sha256': 'aa13165dd6583e54b624948ab841e908938b0c0d06e290bc41bf4fe1f7fa17c3',
        'origin': {'kind': 'committed-source', 'revision': '23ebabd606aaceaca469de59c70ec6d7bed87989',
                   'path': 'sources/legacy-inputs.json'},
    },
    'asset-status': {
        'path': 'sources/actual-data-asset-status.encoded.json',
        'wire_bytes': 50186, 'wire_sha256': '6ec1b15793e5729cb678755e609248f5fdf418cf13de24cfec6a504619b5802c',
        'decoded_bytes': 306466, 'decoded_sha256': 'b13c6c485ede3649f50c260cee60b8b733a76ab0db66ada6746ac21ed9dbfda0',
        'compressed_bytes': 37291, 'compressed_sha256': 'b5022c264e267461f1232fef999fbbac2e3bf6544d24e6f696a4e32444794374',
        'origin': {'kind': 'reviewed-source-proposal', 'revision': None,
                   'path': '.work/actual-data-asset-status-v2-proposal-20261003/manifest.json'},
    },
    'asset-identities': {
        'path': 'sources/actual-data-asset-status-input-identities.encoded.json',
        'wire_bytes': 5664, 'wire_sha256': 'd069001dccc03e247c502d5d900b9751b65983996abbfc6700ed02c42fcad536',
        'decoded_bytes': 17175, 'decoded_sha256': 'bf7d84e2801b56b31de067a44970e73241b6b1248cd765fdaea1135381718556',
        'compressed_bytes': 3893, 'compressed_sha256': '555c1f927b7d30b98c1724c99f2771ab54a84b106fcdbb13feb1dcc9aa209b98',
        'origin': {'kind': 'reviewed-source-proposal', 'revision': None,
                   'path': '.work/actual-data-asset-status-v2-proposal-20261003/input-identities.json'},
    },
}
ROOT_TYPES = {
    'legacy-inputs': {'created_utc': str, 'files': list, 'updated_utc': str, 'format': str,
                      'historical_revision': str, 'scope': str, 'frozen_bundle': dict},
    'asset-status': {'schema': str, 'status': str, 'created_date': str, 'scope': str,
        'successor_lineage': dict, 'scope_relation_semantics': dict, 'rubric': dict,
        'narrow_inventory_basis': dict, 'families': dict, 'contextual_relations': list,
        'scope_map': list, 'joint_probe_gates': list, 'activity': dict, 'assets': list,
        'retained_receipts': list, 'fresh_source_identity_ref': str, 'original_version_policy': str,
        'location_check_summary': dict, 'asset_count': int, 'asset_classification_boundary': str,
        'review_status': str, 'remaining_source_binding_gaps': list},
    'asset-identities': {'schema': str, 'root': str, 'date': str, 'scope': str,
        'source_documents': list, 'source_document_count': int, 'primary_asset_locations': dict,
        'fresh_hash_limit_primary_asset_bytes': int, 'original_inventory_UTC': str,
        'authorization': dict, 'activity': dict, 'changed_document_count': int,
        'canonical_projection_admission_matches': int, 'canonical_projection_admission_mismatches': int,
        'source_hash_phase': str, 'successor_lineage': dict, 'successor_source_document_checks': dict,
        'successor_small_document_semantic_check': dict},
}


class MetadataSourceError(ValueError):
    def __init__(self, document, stage, message, identities=None):
        self.document, self.stage = document, stage
        self.identities = identities
        super().__init__(message)


class _BoundedParser(argparse.ArgumentParser):
    def error(self, message):
        # argparse normally echoes arbitrary supplied arguments without a bound.
        raise MetadataSourceError(None, 'arguments', 'invalid inspection arguments; use --help')


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def _integer(text):
    if len(text.lstrip('-')) > 20:
        raise ValueError('metadata integer token exceeds 20 digits')
    value = int(text)
    if not -(2**63) <= value <= 2**64-1:
        raise ValueError('metadata integer outside signed64/unsigned64 union')
    return value


def _reject_number(text):
    raise ValueError('these frozen metadata documents contain integer JSON numbers only')


def _pairs(items):
    result = {}
    for key, value in items:
        if key in result:
            raise ValueError('duplicate JSON key')
        result[key] = value
    return result


def _json(raw, string_limit):
    # Bound nesting before json.loads, so a bounded byte count cannot bypass depth.
    depth = 0
    quoted = escaped = False
    for value in raw:
        if quoted:
            if escaped:
                escaped = False
            elif value == 92:
                escaped = True
            elif value == 34:
                quoted = False
        elif value == 34:
            quoted = True
        elif value in (91, 123):
            depth += 1
            if depth > DEPTH_LIMIT:
                raise ValueError('metadata nesting limit exceeded')
        elif value in (93, 125):
            depth -= 1
            if depth < 0:
                raise ValueError('metadata nesting is malformed')
    if depth or quoted:
        raise ValueError('metadata structure is incomplete')
    result = json.loads(raw.decode('utf-8', errors='strict'), object_pairs_hook=_pairs,
                        parse_int=_integer, parse_float=_reject_number, parse_constant=_reject_number)
    pending, nodes, graph_bytes, seen = [result], 0, 0, set()
    while pending:
        item = pending.pop()
        nodes += 1
        if nodes > NODE_LIMIT:
            raise ValueError('metadata node limit exceeded')
        if id(item) not in seen:
            seen.add(id(item))
            graph_bytes += sys.getsizeof(item)
            if graph_bytes > GRAPH_SIZE_LIMIT:
                raise ValueError('retained metadata graph shallow-size sum exceeds8MiB')
        if type(item) is str:
            if len(item.encode('utf-8', errors='strict')) > string_limit:
                raise ValueError('metadata string limit exceeded')
        elif type(item) is dict:
            pending.extend(item.keys())
            pending.extend(item.values())
        elif type(item) is list:
            pending.extend(item)
        elif type(item) not in (int, bool, type(None)):
            raise ValueError('unsupported metadata JSON type')
    return result


def _canonical(obj):
    return (json.dumps(obj, sort_keys=True, separators=(',', ':'), ensure_ascii=True,
                       allow_nan=False) + '\n').encode('ascii')


def _read(root, relative, limit):
    # Fixed relative sources, bounded nonblocking regular files, no symlink parent.
    if not root.is_absolute() or len(root.parts) > 64:
        raise ValueError('metadata root is not a bounded absolute path')
    rel = Path(relative)
    if rel.is_absolute() or str(rel) != relative or not 2 <= len(rel.parts) <= 4 or any(
            part in ('', '.', '..') for part in rel.parts):
        raise ValueError('metadata source path is not a canonical confined relative path')
    fd = os.open('/', os.O_RDONLY | os.O_DIRECTORY)
    try:
        for part in root.parts[1:] + rel.parts[:-1]:
            if part in ('', '.', '..'):
                raise ValueError('unsafe metadata path component')
            child = os.open(part, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW, dir_fd=fd)
            os.close(fd)
            fd = child
        leaf = os.open(rel.name, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK, dir_fd=fd)
        try:
            before = os.fstat(leaf)
            if not stat.S_ISREG(before.st_mode) or not 0 < before.st_size <= limit:
                raise ValueError('metadata leaf is nonregular or exceeds its byte bound')
            with os.fdopen(leaf, 'rb', closefd=False) as stream:
                raw = stream.read(limit+1)
            after = os.fstat(leaf)
            if (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns, before.st_ctime_ns) != (
                    after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns, after.st_ctime_ns):
                raise ValueError('metadata leaf changed during its bounded read')
            if len(raw) != before.st_size or len(raw) > limit:
                raise ValueError('metadata consumed byte count differs')
            return raw
        finally:
            os.close(leaf)
    finally:
        os.close(fd)


def _decode(document, raw, spec):
    if len(raw) != spec['wire_bytes'] or len(raw) > WIRE_LIMIT or sha(raw) != spec['wire_sha256']:
        raise ValueError('encoded source identity differs')
    obj = _json(raw, WIRE_LIMIT)
    if type(obj) is not dict or set(obj) != WIRE_FIELDS or _canonical(obj) != raw:
        raise ValueError('encoded source has extra keys, trailing bytes or noncanonical JSON')
    if type(obj['schema_version']) is not int or obj['schema_version'] != 1:
        raise ValueError('encoded schema version is not exact integer1')
    if type(obj['document_id']) is not str or obj['document_id'] != document or obj['encoding'] != FORMAT:
        raise ValueError('encoded source document/encoding differs')
    if type(obj['origin']) is not dict or obj['origin'] != spec['origin']:
        raise ValueError('encoded source origin differs')
    for key, bound in [('decoded_bytes', DECODED_LIMIT), ('compressed_bytes', COMPRESSED_LIMIT)]:
        if type(obj[key]) is not int or not 0 < obj[key] <= bound or obj[key] != spec[key]:
            raise ValueError('encoded source size is not an exact bounded integer')
    for key in ('decoded_sha256', 'compressed_sha256'):
        if type(obj[key]) is not str or HEX.fullmatch(obj[key]) is None or obj[key] != spec[key]:
            raise ValueError('encoded source digest differs from the trusted source pin')
    payload = obj['payload']
    if type(payload) is not str or len(payload) != 4*((obj['compressed_bytes']+2)//3):
        raise ValueError('encoded source payload length/type differs')
    compressed = base64.b64decode(payload.encode('ascii', errors='strict'), validate=True)
    if base64.b64encode(compressed).decode('ascii') != payload:
        raise ValueError('encoded source uses noncanonical base64 padding bits')
    if len(compressed) != obj['compressed_bytes'] or sha(compressed) != obj['compressed_sha256']:
        raise ValueError('compressed source consumed byte identity differs')
    inflater = zlib.decompressobj(wbits=15)
    decoded = inflater.decompress(compressed, obj['decoded_bytes']+1)
    if not inflater.eof or inflater.unused_data or inflater.unconsumed_tail:
        raise ValueError('compressed source is truncated, concatenated, trailing or oversized')
    if len(decoded) != obj['decoded_bytes'] or sha(decoded) != spec['decoded_sha256']:
        raise ValueError('decoded source exact original byte identity differs')
    return decoded


def _schema(document, obj):
    expected = ROOT_TYPES[document]
    if type(obj) is not dict or set(obj) != set(expected) or any(type(obj[key]) is not typ for key, typ in expected.items()):
        raise ValueError('decoded source closed root schema/types differ')
    # Exact decoded digest above binds all deeper fields, statuses, roles and histories.
    if document == 'legacy-inputs':
        if len(obj['files']) != 391:
            raise ValueError('legacy source391 count differs')
        paths = []
        for row in obj['files']:
            if type(row) is not dict or type(row.get('path')) is not str or type(row.get('bytes')) is not int:
                raise ValueError('legacy input record type differs')
            if row['bytes'] < 0 or type(row.get('sha256')) is not str or HEX.fullmatch(row['sha256']) is None:
                raise ValueError('legacy input byte/digest domain differs')
            paths.append(row['path'])
        if len(set(paths)) != 391:
            raise ValueError('legacy source paths are not unique')
    elif document == 'asset-status':
        if obj['asset_count'] != 90 or len(obj['assets']) != 90 or len(obj['retained_receipts']) != 15:
            raise ValueError('asset90/receipt15 identity counts differ')
        if len(obj['scope_map']) != 19 or len(obj['families']) != 11 or len(obj['contextual_relations']) != 4:
            raise ValueError('asset scope/family/context counts differ')
    elif obj['source_document_count'] != 22 or len(obj['source_documents']) != 22:
        raise ValueError('asset companion source22 count differs')


def read_document(document, root=ROOT):
    if type(document) is not str or document not in DOCUMENTS:
        raise MetadataSourceError(document, 'route', 'only the three reviewed source routes are supported')
    spec = DOCUMENTS[document]
    captured = {'expected_transport': {'path': spec['path'], 'bytes': spec['wire_bytes'],
                                       'sha256': spec['wire_sha256']},
                'expected_decoded': {'bytes': spec['decoded_bytes'], 'sha256': spec['decoded_sha256']},
                'consumed_transport': None, 'decoder': None}
    stage = 'decoder-read'
    try:
        decoder = _read(root, 'scripts/metadata_source.py', 32768)
        captured['decoder'] = {'path': 'scripts/metadata_source.py', 'bytes': len(decoder), 'sha256': sha(decoder)}
        stage = 'read'
        raw = _read(root, spec['path'], WIRE_LIMIT)
        captured['consumed_transport'] = {'path': spec['path'], 'bytes': len(raw), 'sha256': sha(raw)}
        stage = 'decode'
        decoded = _decode(document, raw, spec)
        stage = 'schema'
        obj = _json(decoded, 4096)
        _schema(document, obj)
        identity = {'format': 'reproducible-metadata-source/v1', 'document_id': document,
                    'transport': {'path': spec['path'], 'bytes': len(raw), 'sha256': sha(raw)},
                    'decoded': {'bytes': len(decoded), 'sha256': sha(decoded), 'origin': spec['origin']},
                    'decoder': captured['decoder']}
        return obj, decoded, identity
    except (OSError, ValueError, UnicodeError, binascii.Error, zlib.error) as exc:
        raise MetadataSourceError(document, stage, str(exc), captured) from exc


def load_legacy_inputs(root=ROOT):
    return read_document('legacy-inputs', root)[0]


def main():
    document = None
    try:
        if len(sys.argv) > 3 or any(len(arg.encode('utf-8', errors='strict')) > 64 for arg in sys.argv[1:]):
            raise MetadataSourceError(None, 'arguments', 'inspection argument count or byte limit exceeded')
        parser = _BoundedParser(description=__doc__)
        parser.add_argument('document', choices=tuple(DOCUMENTS))
        parser.add_argument('--identity', action='store_true', help='show transport, decoded and decoder identities')
        args = parser.parse_args()
        document = args.document
        _, decoded, identity = read_document(args.document)
        sys.stdout.buffer.write(_canonical(identity) if args.identity else decoded)
    except (MetadataSourceError, UnicodeError) as exc:
        # Bounded refusal, no writes/downloads or fabricated partial source object.
        error = 'inspection arguments must be valid UTF-8' if isinstance(exc, UnicodeError) else str(exc)[:256]
        print(json.dumps({'document_id': document, 'status': 'refused',
                          'stage': getattr(exc, 'stage', 'arguments'),
                          'error': error, 'identities': getattr(exc, 'identities', None)},
                         sort_keys=True), file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
