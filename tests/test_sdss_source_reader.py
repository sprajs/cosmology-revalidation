"""Light structural/adversarial checks only; no native or scientific execution."""
import copy,gzip,importlib.util,io,json,struct,tempfile,unittest
from unittest.mock import patch
from types import SimpleNamespace
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result)
    return result
fixtures=module('sdss_source_fixtures',Path(__file__).with_name('sdss_source_fixtures.py'))
reader=module('sdss_source_reader',ROOT/'experiments/sdss-released-observer-contract/reader.py')

class StructuralReaderChecks(unittest.TestCase):
    def runtime_fixture(self):
        return {'path':reader.text_identity('/fixture/python'),'bytes':8,'sha256':reader.sha(b'python!!'),'version':reader.text_identity('fixture Python'),'implementation':reader.text_identity('cpython')}
    def fixture(self):
        b,s,r=fixtures.fits_one_row();return reader.GzipReader(b,20000),s,r
    def test_source_hash_and_size_before_parsing(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'source';p.write_bytes(b'original')
            self.assertEqual(reader.admit(p,{'bytes':8,'sha256':reader.sha(b'original')}),b'original')
            with self.assertRaisesRegex(reader.Refusal,'hash differs'):reader.admit(p,{'bytes':8,'sha256':reader.sha(b'changed!')})
            with self.assertRaisesRegex(reader.Refusal,'length differs'):reader.admit(p,{'bytes':7,'sha256':reader.sha(b'original')})
    def test_quoted_slash_and_doubled_quote(self):
        self.assertEqual(reader.value_card("'/source/path' / comment"),'/source/path')
        self.assertEqual(reader.value_card("'O''Brien' / comment"),"O'Brien")
        with self.assertRaises(reader.Refusal):reader.value_card("'unclosed")
    def test_nul_padding_and_raw_negative_flux(self):
        s,expected,raw=self.fixture();meta=reader.table_header(s,expected,64)
        parsed=reader.decode(s.exact(20),meta['columns'])
        self.assertEqual(parsed,{'SNID':'6057','FLUXCAL':-3.25})
        self.assertEqual(reader.raw_strings(raw,meta['columns'])['SNID'],raw[:16].hex())
        self.assertIsNone(meta['columns'][1]['qualified_unit']);self.assertEqual(meta['columns'][1]['source_TUNIT'],'')
        with self.assertRaisesRegex(reader.Refusal,'interior NUL'):reader.ascii_field(b'a\0b')
    def test_exact_closed_schema_refuses_extra_scaling(self):
        b,expected,_=fixtures.fits_one_row({'TSCAL2':'2'})
        original=copy.deepcopy(expected);original[1].pop('TSCAL2')
        with self.assertRaisesRegex(reader.Refusal,'schema differs'):reader.table_header(reader.GzipReader(b,20000),original,64)
    def test_duplicate_columns_and_keys(self):
        b,expected,_=fixtures.fits_one_row({'TTYPE2':'SNID'})
        with self.assertRaisesRegex(reader.Refusal,'duplicate FITS column'):reader.table_header(reader.GzipReader(b,20000),expected,64)
        data=fixtures.head([fixtures.card('NAXIS','0'),fixtures.card('NAXIS','0')])
        with self.assertRaisesRegex(reader.Refusal,'duplicate FITS key'):reader.header(reader.GzipReader(gzip.compress(data),10000),64)
    def test_bounded_gzip_bomb_and_crc(self):
        stream=reader.GzipReader(gzip.compress(b'x'*1024),128)
        with self.assertRaisesRegex(reader.Refusal,'inflated source limit'):stream.exact(129)
        b=bytearray(gzip.compress(b'payload'));b[-8]^=1
        stream=reader.GzipReader(bytes(b),100);stream.exact(7)
        with self.assertRaisesRegex(reader.Refusal,'gzip integrity failure'):stream.end()
    def test_truncation_and_missing_end(self):
        with self.assertRaises(reader.Refusal):reader.GzipReader(gzip.compress(b'x'),20).exact(2)
        data=fixtures.card('NAXIS','0')+fixtures.card('BITPIX','8')+fixtures.card('SIMPLE','T')
        with self.assertRaisesRegex(reader.Refusal,'card limit'):reader.header(reader.GzipReader(gzip.compress(data),1000),2)
    def test_extra_payload_refused(self):
        s=reader.GzipReader(gzip.compress(b'abc'),20);s.exact(2)
        with self.assertRaisesRegex(reader.Refusal,'extra FITS'):s.end()
    def test_pointers_boundaries_and_nobs(self):
        reader.pointer(1,3,3,3,3)
        for values in [(0,2,3,3),(1,4,4,3),(2,1,0,3),(1,3,2,3)]:
            with self.assertRaises(reader.Refusal):reader.pointer(*values)
        with self.assertRaisesRegex(reader.Refusal,'event record limit'):reader.pointer(1,3,3,3,2)
    def test_original_declared_decimal_bounds_not_float_slack(self):
        self.assertTrue(reader.coordinate('1.0005','1.000')['within_printed_bound'])
        self.assertFalse(reader.coordinate('1.000500000001','1.000')['within_printed_bound'])
        self.assertFalse(reader.coordinate('337.672525','337.672')['within_printed_bound'])
        self.assertEqual(reader.coordinate('1.000','1.000')['difference_degree'],'0.000')
    def test_duplicate_join_retained_then_refused(self):
        rows=[{'CID':'0006057'},{'CID':'6057'}];i=reader.index(rows,'CID',True)
        self.assertEqual(len(i['6057']),2)
        with self.assertRaisesRegex(reader.Refusal,'match count 2'):reader.unique(i,'6057','fixture')
        with self.assertRaisesRegex(reader.Refusal,'match count 0'):reader.unique(i,'unknown','fixture')
    def test_text_schema_and_limits(self):
        c={'bounds':{'maximum_text_rows':2,'maximum_text_line_bytes':32},'final_columns':['CID','IDSURVEY'],'selection':{'release_rows':1}}
        self.assertEqual(reader.final_rows(b'CID IDSURVEY\n6057 1\n',c)[0]['CID'],'6057')
        for b in [b'IDSURVEY CID\n1 6057',b'CID IDSURVEY\n6057 1 extra',b'CID IDSURVEY\n1 1\n2 1\n3 1']:
            with self.assertRaises(reader.Refusal):reader.final_rows(b,c)
    def test_nonfinite_source_value_is_refused(self):
        with self.assertRaisesRegex(reader.Refusal,'nonfinite FITS'):reader.decode(struct.pack('>f',float('nan')),[{'name':'FLUXCAL','offset':0,'bytes':4,'format':'1E'}])
    def test_explicit_decimal_context_retains_beyond28_digits(self):
        with reader.decimal.localcontext() as c:
            c.prec=6
            result=reader.exact_decimal_difference('1.00000000000000000000000000000000001','1')
            self.assertEqual(str(result),'1E-35')
            coordinate=reader.coordinate('1.00000000000000000000000000000000001','1.000')
            self.assertEqual(coordinate['difference_degree'],'1E-35')
    def test_output_existing_traversal_and_symlink_refused(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);reader.fresh_output(root/'new-attempt',root)
            (root/'existing').mkdir()
            (root/'link').symlink_to(root/'missing')
            for p in [root/'existing',root/'link',root/'..'/'escape',root/'nested'/'attempt',root/('x'*65)]:
                with self.assertRaises(reader.Refusal):reader.fresh_output(p,root)
    def test_packet_output_root_is_fixed_and_refuses_symlink_parents(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);output=reader.packet_output_root(root)
            self.assertEqual(output,root/'results/sdss-released-observer-contract')
            with self.assertRaises(reader.Refusal):reader.fresh_output(output/'..'/'escape',output)
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);(root/'outside').mkdir();(root/'results').symlink_to(root/'outside',target_is_directory=True)
            with self.assertRaisesRegex(reader.Refusal,'symlink'):reader.packet_output_root(root)
    def test_identity_errors_and_exact_original_pins(self):
        for bad in [{'sha256':None,'error':'missing'},{'bytes':8,'sha256':'not-a-hash'},{'bytes':0,'sha256':reader.sha(b'')}]:
            with self.assertRaisesRegex(reader.Refusal,'identity not admitted'):reader.require_identities({'fixture':bad},'packet')
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);(root/'fixture').write_bytes(b'original')
            pins={'fixture':{'bytes':8,'sha256':reader.sha(b'original')}}
            self.assertEqual(reader.source_identity_errors(reader.source_identities(root,pins),pins),[])
            (root/'fixture').write_bytes(b'changed!')
            self.assertEqual(len(reader.source_identity_errors(reader.source_identities(root,pins),pins)),1)
            (root/'fixture').write_bytes(b'longer-than-bounded')
            self.assertEqual(reader.source_identities(root,pins)['fixture']['error'],'identity byte limit')
    def test_git_source_identity_with_mocked_commands(self):
        root=Path('/fixture-root');paths={n:root/n for n in ['reader','contract','candidate','experiment','README']}
        rows='\n'.join('100644 '+'b'*40+' 0\t'+n for n in paths)+'\n'
        def replies(status=0,listing=rows):
            return [SimpleNamespace(returncode=0,stdout=b'/fixture-root\n'),SimpleNamespace(returncode=0,stdout=('a'*40+'\n').encode()),SimpleNamespace(returncode=status,stdout=b''),SimpleNamespace(returncode=0,stdout=listing.encode())]
        with patch.object(reader.subprocess,'run',side_effect=replies()) as run,patch.object(reader,'bounded_identity',return_value={'git_blob_sha1':'b'*40}):
            value=reader.source_control_identity(root,paths);reader.require_source_control(value)
            self.assertEqual(value['revision'],'a'*40);self.assertEqual(len(value['packet_files']),5)
            self.assertTrue(all(call.kwargs['timeout']==10 for call in run.call_args_list))
        for responses,actual in [(replies(status=1),'b'*40),(replies(listing=rows.splitlines()[0]+'\n'),'b'*40),(replies(),'c'*40)]:
            with patch.object(reader.subprocess,'run',side_effect=responses),patch.object(reader,'bounded_identity',return_value={'git_blob_sha1':actual}):
                value=reader.source_control_identity(root,paths)
                with self.assertRaises(reader.Refusal):reader.require_source_control(value)
    def main_fixture(self,root):
        folder=root/'experiments/sdss-released-observer-contract';folder.mkdir(parents=True)
        candidate=b'candidate';(folder/'candidate.json').write_bytes(candidate)
        (folder/'README.md').write_text('fixture README');(folder/'experiment.json').write_text('{}')
        source=root/'inputs';source.mkdir();(source/'fixture').write_bytes(b'original')
        contract={'schema':1,'candidate_path':'candidate.json','candidate_sha256':reader.sha(candidate),'sources':{'fixture':{'bytes':8,'sha256':reader.sha(b'original')}},'bounds':{'maximum_compressed_total_bytes':8,'maximum_output_bytes':4194304}}
        cb=json.dumps(contract).encode();(folder/'input-contract.json').write_bytes(cb)
        git={'revision':'a'*40,'tracked_checkout_clean':True,'packet_files':{'fixture':{'git_blob_sha1':'b'*40}}}
        result={'events':[{}],'coordinate_failures':[{}],'photometry_slices':{'fixture':{'records':[{}]}},'signed_negative_flux_records':1}
        return folder,source,reader.sha(cb),git,result
    def test_main_retains_lineage_and_failure_on_original_source_drift(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);folder,source,digest,git,result=self.main_fixture(root)
            def regenerate(contract,blobs):
                self.assertEqual(blobs['fixture'],b'original')
                (source/'fixture').write_bytes(b'changed!');return result
            with patch.multiple(reader,ROOT=root,FOLDER=folder,EXPECTED_CONTRACT_SHA256=digest),patch.object(reader,'source_control_identity',return_value=git),patch.object(reader,'runtime_identity',return_value=self.runtime_fixture()),patch.object(reader,'regenerate',side_effect=regenerate),patch('sys.argv',['reader','--input-root',str(source),'--attempt','drift']),patch('builtins.print'):
                with self.assertRaises(SystemExit):reader.main()
            output=root/'results/sdss-released-observer-contract/drift'
            receipt=json.loads((output/'receipt.json').read_bytes())
            self.assertEqual(receipt['status'],'failed_identity_changed')
            self.assertEqual(receipt['original_sources_before']['fixture']['sha256'],reader.sha(b'original'))
            self.assertEqual(receipt['original_sources_after']['fixture']['sha256'],reader.sha(b'changed!'))
            self.assertTrue((output/'lineage.json').is_file())
            self.assertEqual((output/'lineage.json').stat().st_mode & 0o777,0o444)
            self.assertEqual(output.stat().st_mode & 0o777,0o555)
    def test_main_refuses_identical_before_after_packet_identity_errors(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);folder,source,digest,git,result=self.main_fixture(root)
            identity=reader.bounded_identity
            def checked(path,limit):
                if Path(path).name=='README.md':return {'path':str(path),'sha256':None,'error':'fixture missing'}
                return identity(path,limit)
            with patch.multiple(reader,ROOT=root,FOLDER=folder,EXPECTED_CONTRACT_SHA256=digest),patch.object(reader,'source_control_identity',return_value=git),patch.object(reader,'runtime_identity',return_value=self.runtime_fixture()),patch.object(reader,'bounded_identity',side_effect=checked),patch.object(reader,'regenerate') as regenerate,patch('sys.argv',['reader','--input-root',str(source),'--attempt','missing']),patch('builtins.print'):
                with self.assertRaises(SystemExit):reader.main()
                regenerate.assert_not_called()
            output=root/'results/sdss-released-observer-contract/missing'
            receipt=json.loads((output/'receipt.json').read_bytes())
            self.assertTrue(receipt['status'].startswith('failed'));self.assertIn('identity not admitted',receipt['error'])
            self.assertFalse((output/'lineage.json').exists())
    def test_long_input_argv_refuses_before_any_output(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            with patch.object(reader,'ROOT',root),patch('sys.argv',['reader','--input-root','/'+('x'*1024),'--attempt','long-input']):
                with self.assertRaisesRegex(reader.Refusal,'no attempt created'):reader.main()
            self.assertFalse((root/'results').exists())
    def test_long_exception_is_hashed_bounded_and_terminal_readonly(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);folder,source,digest,git,result=self.main_fixture(root)
            diagnostic='\u0001'*100000
            with patch.multiple(reader,ROOT=root,FOLDER=folder,EXPECTED_CONTRACT_SHA256=digest),patch.object(reader,'source_control_identity',return_value=git),patch.object(reader,'runtime_identity',return_value=self.runtime_fixture()),patch.object(reader,'regenerate',side_effect=RuntimeError(diagnostic)),patch('sys.argv',['reader','--input-root',str(source),'--attempt','long-error']),patch('builtins.print'):
                with self.assertRaises(SystemExit):reader.main()
            output=root/'results/sdss-released-observer-contract/long-error';blob=(output/'receipt.json').read_bytes();receipt=json.loads(blob)
            self.assertEqual(receipt['status'],'failed');self.assertLessEqual(len(blob),65536)
            self.assertEqual(receipt['error_identity']['sha256'],reader.sha(diagnostic.encode()))
            self.assertTrue(receipt['error_identity']['truncated']);self.assertLessEqual(len(receipt['error'].encode()),256)
            self.assertEqual(output.stat().st_mode & 0o777,0o555);self.assertEqual((output/'receipt.json').stat().st_mode & 0o777,0o444)
    def test_fallback_receipt_retains_used_hashes_under_unicode_expansion(self):
        identity={'path':reader.text_identity('\u0001'*10000),'bytes':8,'sha256':reader.sha(b'original')}
        packet={name:identity for name in ['reader','contract','candidate','experiment','README']}
        sources={'source-'+str(i):identity for i in range(9)}
        git={'revision':'a'*40,'packet_files':{name:{'git_blob_sha1':'b'*40} for name in packet}}
        receipt={'status':'failed','error':'\u0001'*100000,'contract_sha256':'a'*64,'lineage_sha256':'b'*64,'identities_before':packet,'identities_after':packet,'original_sources_before':sources,'original_sources_after':sources,'admitted_sources':sources,'source_control_before':git,'source_control_after':git,'runtime_before':self.runtime_fixture(),'runtime_after':self.runtime_fixture(),'integrity_errors':['\u0001'*100000]*32}
        fallback,encoded=reader.terminal_receipt(receipt)
        self.assertLessEqual(len(encoded),65536);self.assertEqual(fallback['status'],'failed_receipt_bound')
        self.assertEqual(fallback['lineage_sha256'],'b'*64);self.assertEqual(fallback['admitted_sources'][0]['sha256'],identity['sha256'])
        self.assertEqual(fallback['original_sources_before'][0]['path_text_sha256'],identity['path']['sha256'])
        self.assertEqual(fallback['runtime_before']['executable'][0]['sha256'],self.runtime_fixture()['sha256'])
    def test_main_oversized_receipt_preserves_lineage_and_readonly_failure(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);folder,source,digest,git,result=self.main_fixture(root);git['fixture-extra']='\u0001'*100000
            with patch.multiple(reader,ROOT=root,FOLDER=folder,EXPECTED_CONTRACT_SHA256=digest),patch.object(reader,'source_control_identity',return_value=git),patch.object(reader,'runtime_identity',return_value=self.runtime_fixture()),patch.object(reader,'regenerate',return_value=result),patch('sys.argv',['reader','--input-root',str(source),'--attempt','fallback']),patch('builtins.print'):
                with self.assertRaises(SystemExit):reader.main()
            output=root/'results/sdss-released-observer-contract/fallback';blob=(output/'receipt.json').read_bytes();receipt=json.loads(blob)
            self.assertEqual(receipt['status'],'failed_receipt_bound');self.assertLessEqual(len(blob),65536)
            self.assertEqual(receipt['lineage_sha256'],reader.sha((output/'lineage.json').read_bytes()))
            self.assertEqual((output/'lineage.json').stat().st_mode & 0o777,0o444);self.assertEqual(output.stat().st_mode & 0o777,0o555)
    def test_python_runtime_identity_uses_bounded_exact_executable_bytes(self):
        with tempfile.TemporaryDirectory() as d:
            executable=Path(d)/'python';executable.write_bytes(b'python!!')
            with patch.object(reader.sys,'executable',str(executable)):
                identity=reader.runtime_identity()
            self.assertEqual(identity['bytes'],8);self.assertEqual(identity['sha256'],reader.sha(b'python!!'))
            self.assertEqual(identity['path']['text'],str(executable.resolve()))
            self.assertEqual(identity['version']['sha256'],reader.sha(reader.sys.version.encode('utf-8',errors='surrogatepass')))
    def test_runtime_change_retains_lineage_and_revokes_structural_status(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);folder,source,digest,git,result=self.main_fixture(root)
            before=self.runtime_fixture();after=copy.deepcopy(before);after['sha256']=reader.sha(b'changed!')
            with patch.multiple(reader,ROOT=root,FOLDER=folder,EXPECTED_CONTRACT_SHA256=digest),patch.object(reader,'source_control_identity',return_value=git),patch.object(reader,'runtime_identity',side_effect=[before,after]),patch.object(reader,'regenerate',return_value=result),patch('sys.argv',['reader','--input-root',str(source),'--attempt','runtime-drift']),patch('builtins.print'):
                with self.assertRaises(SystemExit):reader.main()
            output=root/'results/sdss-released-observer-contract/runtime-drift';receipt=json.loads((output/'receipt.json').read_bytes())
            self.assertEqual(receipt['status'],'failed_identity_changed');self.assertIn('Python executable/version identity changed',receipt['integrity_errors'])
            self.assertEqual(receipt['runtime_before']['sha256'],before['sha256']);self.assertEqual(receipt['runtime_after']['sha256'],after['sha256'])
            self.assertEqual(receipt['lineage_sha256'],reader.sha((output/'lineage.json').read_bytes()))

if __name__=='__main__':unittest.main(verbosity=2)
