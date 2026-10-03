"""Metadata-only fixtures; never import native/scientific controllers."""
import base64
import contextlib
import hashlib
import io
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import zlib

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
import metadata_source as M
import restore_inputs


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def canonical(obj):
    return (json.dumps(obj, sort_keys=True, separators=(',', ':'), ensure_ascii=True,
                       allow_nan=False)+'\n').encode('ascii')


def fixture(decoded=b'{"value":"exact\\nmetadata","count":2}\n', compressed=None):
    compressed = zlib.compress(decoded, 9) if compressed is None else compressed
    origin = {'kind':'fixture', 'revision':None, 'path':'fixture.json'}
    obj = {'schema_version':1, 'document_id':'fixture', 'encoding':M.FORMAT,
           'decoded_bytes':len(decoded), 'decoded_sha256':digest(decoded),
           'compressed_bytes':len(compressed), 'compressed_sha256':digest(compressed),
           'origin':origin, 'payload':base64.b64encode(compressed).decode('ascii')}
    spec = {k:obj[k] for k in ('decoded_bytes','decoded_sha256','compressed_bytes','compressed_sha256','origin')}
    raw = canonical(obj)
    spec.update(wire_bytes=len(raw), wire_sha256=digest(raw))
    return obj, spec, raw


def repin_wire(obj, spec):
    raw = canonical(obj)
    spec = dict(spec, wire_bytes=len(raw), wire_sha256=digest(raw))
    return raw, spec


class MetadataSourceTests(unittest.TestCase):
    def test_exact_bytes_include_whitespace_and_final_newline(self):
        original = b'{\n  "count" : 2, "value" : "backslash\\\\end"\n}\n'
        _, spec, raw = fixture(original)
        self.assertEqual(M._decode('fixture', raw, spec), original)
        self.assertNotEqual(original, canonical(json.loads(original)))

    def test_real_named_sources_exact_original_hash_and_counts(self):
        expected = {
            'legacy-inputs':(255902,'7ed81f37e617ce7bb2638f30596165f46ae7ee81ca31788177e90cf57ad5e56f'),
            'asset-status':(306466,'b13c6c485ede3649f50c260cee60b8b733a76ab0db66ada6746ac21ed9dbfda0'),
            'asset-identities':(17175,'bf7d84e2801b56b31de067a44970e73241b6b1248cd765fdaea1135381718556')}
        for name, (size, sha) in expected.items():
            with self.subTest(name=name):
                obj, raw, identity = M.read_document(name)
                self.assertEqual((len(raw),digest(raw)), (size,sha))
                self.assertEqual(obj,json.loads(raw))  # Original integer metadata parser semantics.
                self.assertEqual(identity['decoded']['sha256'],sha)
                self.assertNotEqual(identity['transport']['sha256'],sha)
                self.assertEqual(identity['decoder']['sha256'],digest((ROOT/'scripts/metadata_source.py').read_bytes()))
        self.assertEqual(len(M.read_document('legacy-inputs')[0]['files']),391)
        status = M.read_document('asset-status')[0]
        self.assertEqual((len(status['assets']),len(status['retained_receipts']),len(status['scope_map'])),(90,15,19))
        self.assertEqual(len(M.read_document('asset-identities')[0]['source_documents']),22)

    def test_changed_consumed_wire_fails_trusted_pin(self):
        _, spec, raw = fixture()
        with self.assertRaisesRegex(ValueError,'encoded source identity'):
            M._decode('fixture',raw.replace(b'exact',b'other') if b'exact' in raw else raw[:-2]+b'!\n',spec)

    def test_self_reported_decoded_digest_is_not_authority(self):
        obj,spec,_ = fixture()
        obj['decoded_sha256'] = '0'*64
        raw,spec = repin_wire(obj,spec)  # Make the outer pin pass; retain trusted original raw pin.
        with self.assertRaisesRegex(ValueError,'trusted source pin'):
            M._decode('fixture',raw,spec)

    def test_typed_size_and_version_reject_bool(self):
        for key in ('schema_version','decoded_bytes','compressed_bytes'):
            with self.subTest(key=key):
                obj,spec,_ = fixture()
                obj[key] = True
                raw,spec = repin_wire(obj,spec)
                with self.assertRaises(ValueError):
                    M._decode('fixture',raw,spec)

    def test_closed_document_origin_and_encoding(self):
        for key, value in [('extra',1),('document_id','alias'),('encoding','gzip'),
                           ('origin',{'kind':'fixture','revision':None,'path':'other'})]:
            with self.subTest(key=key):
                obj,spec,_ = fixture();obj[key] = value
                raw,spec = repin_wire(obj,spec)
                with self.assertRaises(ValueError):M._decode('fixture',raw,spec)
        with self.assertRaisesRegex(M.MetadataSourceError,'three reviewed'):
            M.read_document('arbitrary-path')

    def test_noncanonical_json_and_duplicate_keys(self):
        _,spec,raw = fixture()
        for changed in (b' '+raw, raw+b' ', raw[:-2]+b',"schema_version":1}\n'):
            with self.subTest(changed=changed[:40]):
                changed_spec = dict(spec,wire_bytes=len(changed),wire_sha256=digest(changed))
                with self.assertRaises(ValueError):M._decode('fixture',changed,changed_spec)

    def test_canonical_base64_unused_padding_bits_and_alphabet(self):
        for extra in range(8):
            obj,spec,_ = fixture(b'{"value":"'+b'x'*extra+b'"}')
            if obj['payload'].endswith('='):break
        alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/'
        payload = obj['payload']
        self.assertTrue(payload.endswith('='))
        at = -3 if payload.endswith('==') else -2
        modified = payload[:at]+alphabet[alphabet.index(payload[at])|1]+payload[at+1:]
        self.assertEqual(base64.b64decode(modified,validate=True),base64.b64decode(payload,validate=True))
        obj['payload'] = modified;raw,spec = repin_wire(obj,spec)
        with self.assertRaisesRegex(ValueError,'noncanonical base64'):M._decode('fixture',raw,spec)
        for changed in (payload+'\n',payload.replace('=',''),payload+'=='):
            obj['payload'] = changed;raw,altered = repin_wire(obj,spec)
            with self.assertRaises(ValueError):M._decode('fixture',raw,altered)

    def test_complete_stream_refuses_trailing_concatenated_and_truncated(self):
        decoded = b'{"value":"unchanged"}\n';compressed = zlib.compress(decoded,9)
        for changed in (compressed+b'junk',compressed+zlib.compress(b'other'),compressed[:-1]):
            with self.subTest(length=len(changed)):
                _,spec,raw = fixture(decoded,changed)
                with self.assertRaisesRegex(ValueError,'truncated, concatenated, trailing or oversized'):
                    M._decode('fixture',raw,spec)

    def test_bounded_decompression_refuses_highly_compressible_excess(self):
        obj,spec,_ = fixture(b'{"value":"'+b'x'*100000+b'"}')
        obj['decoded_bytes'] = 32;spec['decoded_bytes'] = 32
        raw,spec = repin_wire(obj,spec)
        with self.assertRaisesRegex(ValueError,'oversized'):M._decode('fixture',raw,spec)

    def test_bounded_encoded_and_declared_size(self):
        obj,spec,_ = fixture();obj['decoded_bytes'] = M.DECODED_LIMIT+1
        spec['decoded_bytes'] = obj['decoded_bytes'];raw,spec = repin_wire(obj,spec)
        with self.assertRaisesRegex(ValueError,'bounded integer'):M._decode('fixture',raw,spec)
        with self.assertRaisesRegex(ValueError,'encoded source identity'):
            M._decode('fixture',b'x'*(M.WIRE_LIMIT+1),dict(spec,wire_bytes=M.WIRE_LIMIT+1))

    def test_decoded_duplicate_nonfinite_fraction_and_invalid_utf8(self):
        for raw in (b'{"a":1,"a":2}',b'{"a":NaN}',b'{"a":Infinity}',b'{"a":1.25}',
                    b'{"a":1e1000}',b'{"a":18446744073709551616}',b'"\xff"'):
            with self.subTest(raw=raw),self.assertRaises((ValueError,UnicodeError)):M._json(raw,4096)

    def test_depth_nodes_string_and_graph_object_bounds(self):
        cases = [b'['*33+b'0'+b']'*33, canonical([0]*16384),canonical('x'*4097)]
        for raw in cases:
            with self.subTest(length=len(raw)),self.assertRaises(ValueError):M._json(raw,4096)
        with patch.object(M,'GRAPH_SIZE_LIMIT',1),self.assertRaisesRegex(ValueError,'shallow-size'):
            M._json(b'{}',4096)

    def test_decoded_root_types_reject_boolean_asset_count(self):
        obj,_,_ = M.read_document('asset-status')
        obj['asset_count'] = True
        with self.assertRaisesRegex(ValueError,'schema/types'):M._schema('asset-status',obj)

    def test_path_confinement_rejects_leaf_parent_symlink_and_fifo(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary);sources = root/'sources';sources.mkdir()
            outside = root/'outside';outside.write_bytes(b'private')
            leaf = sources/'value.json';leaf.symlink_to(outside)
            with self.assertRaises(OSError):M._read(root,'sources/value.json',20)
            leaf.unlink();os.mkfifo(leaf)
            with self.assertRaisesRegex(ValueError,'nonregular'):M._read(root,'sources/value.json',20)
            leaf.unlink();sources.rmdir();sources.symlink_to(root)
            with self.assertRaises(OSError):M._read(root,'sources/outside',20)
            with self.assertRaises(ValueError):M._read(root,'../outside',20)
            with self.assertRaises(ValueError):M._read(root,str(outside),20)
            with self.assertRaises(ValueError):M._read(root,'sources//outside',20)

    def test_corrupt_input_preserved_refusal_precedes_restoration(self):
        error = M.MetadataSourceError('legacy-inputs','decode','trusted byte pin differs')
        with patch.object(restore_inputs,'read_document',side_effect=error), \
                patch.object(restore_inputs.urllib.request,'urlopen') as network, \
                patch.object(sys,'argv',['restore_inputs.py','--group','bao']), \
                contextlib.redirect_stdout(io.StringIO()) as out, \
                self.assertRaises(SystemExit) as stopped:
            restore_inputs.main()
        self.assertEqual(stopped.exception.code,1);network.assert_not_called()
        self.assertEqual(json.loads(out.getvalue())['source_admission'],'refused')
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary);p = root/'sources/legacy-inputs.encoded.json';p.parent.mkdir()
            helper = root/'scripts/metadata_source.py';helper.parent.mkdir();helper.write_bytes(b'fixture source')
            p.write_bytes(b'corrupt');p.chmod(0o444)
            original = (p.read_bytes(),p.stat().st_mode)
            with self.assertRaises(M.MetadataSourceError) as refused:M.read_document('legacy-inputs',root)
            self.assertEqual(refused.exception.stage,'decode')
            self.assertEqual(refused.exception.identities['consumed_transport']['sha256'],digest(b'corrupt'))
            self.assertEqual(refused.exception.identities['decoder']['sha256'],digest(b'fixture source'))
            self.assertEqual((p.read_bytes(),p.stat().st_mode),original)

    def test_restore_list_and_dryrun_consume_admitted_metadata_without_downloads(self):
        manifest,_,identity = M.read_document('legacy-inputs')
        for args in (['restore_inputs.py','--list'],['restore_inputs.py','--group','bao','--dry-run']):
            with self.subTest(args=args),patch.object(restore_inputs,'read_document',return_value=(manifest,b'',identity)) as reader, \
                    patch.object(restore_inputs.urllib.request,'urlopen') as network, \
                    patch.object(restore_inputs,'restore') as restoration,patch.object(sys,'argv',args), \
                    contextlib.redirect_stdout(io.StringIO()) as out:
                restore_inputs.main()
            reader.assert_called_once_with('legacy-inputs',restore_inputs.ROOT)
            network.assert_not_called();restoration.assert_not_called()
            if '--list' in args:
                self.assertEqual(sum(v['files'] for v in json.loads(out.getvalue()).values()),391)
            else:
                output = out.getvalue();summary = json.loads(output[output.index('{\n'):])
                expected = [r for r in manifest['files'] if r['group']=='bao']
                self.assertEqual(summary['selected'],len(expected));self.assertEqual(summary['failures'],[])
                rows = [json.loads(line) for line in output.splitlines() if line.startswith('{"path"')]
                self.assertEqual([r['path'] for r in rows],[r['path'] for r in expected])
                self.assertEqual(summary['manifest_source'],identity)

    def test_restoration_preserves_changed_existing_file(self):
        row = {'path':'data/source.txt','bytes':8,'sha256':digest(b'original')}
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary);path = root/row['path'];path.parent.mkdir(parents=True)
            path.write_bytes(b'changed!');path.chmod(0o444)
            before = (path.read_bytes(),path.stat().st_mode)
            with patch.object(restore_inputs.urllib.request,'urlopen') as network, \
                    self.assertRaisesRegex(ValueError,'preserved without replacement'):
                restore_inputs.restore(row,root,root/'unused-archive',{})
            network.assert_not_called();self.assertEqual((path.read_bytes(),path.stat().st_mode),before)

    def test_cli_machine_identity_and_exact_stdout_without_writes(self):
        raw = b'{\n "value": "preserved"\n}\n';identity = {'transport':{'sha256':'1'*64},'decoded':{'sha256':'2'*64}}
        for args, expected in [(['metadata_source.py','asset-status'],raw),
                               (['metadata_source.py','asset-status','--identity'],canonical(identity))]:
            output = io.BytesIO()
            class Stdout:
                buffer = output
            with patch.object(M,'read_document',return_value=({},raw,identity)), \
                    patch.object(sys,'argv',args),patch.object(sys,'stdout',Stdout()):
                self.assertEqual(M.main(),0)
            self.assertEqual(output.getvalue(),expected)

    def test_cli_argument_refusal_is_bounded_json_before_source_reads(self):
        cases = [[], ['unknown'], ['x'*100000], ['legacy-inputs','--'+'x'*100000],
                 ['legacy-inputs','--unknown'], ['legacy-inputs','--identity','extra'],
                 ['\U0001f30c'*10000], ['\udcff']]
        for args in cases:
            with self.subTest(args=[a[:20] for a in args]), \
                    patch.object(sys,'argv',['metadata_source.py',*args]), \
                    patch.object(M,'read_document') as reader, \
                    contextlib.redirect_stderr(io.StringIO()) as err, \
                    contextlib.redirect_stdout(io.StringIO()) as out:
                self.assertEqual(M.main(),1)
            reader.assert_not_called()
            self.assertEqual(out.getvalue(),'')
            self.assertLessEqual(len(err.getvalue().encode('utf-8')),4096)
            self.assertNotIn('\\udcff',err.getvalue())
            refusal = json.loads(err.getvalue())
            self.assertEqual((refusal['status'],refusal['stage']),('refused','arguments'))
            self.assertIsNone(refusal['document_id'])
            self.assertIsNone(refusal['identities'])


if __name__ == '__main__':
    unittest.main()
