"""Original bounded reader for the pinned SDSS source-only lineage slice.

No physical flux, exposure, frame, velocity, template or detector calculation.
Only this exact release/schema is admitted; no generic FITS/workflow interface.
"""
from __future__ import annotations
import argparse, csv, decimal, gzip, hashlib, io, json, math, re, struct, subprocess, sys
from pathlib import Path

D = decimal.Decimal
ROOT = Path(__file__).resolve().parents[2]
FOLDER = Path(__file__).resolve().parent
EXPECTED_CONTRACT_SHA256 = '8e7690758e32f5c9535e11cd5f52b365d6d3b915c932b5649636ded044883c16'
class Refusal(ValueError): pass
def need(condition, reason):
    if not condition: raise Refusal(reason)
def sha(blob): return hashlib.sha256(blob).hexdigest()

def text_identity(value,limit=256):
    raw=str(value).encode('utf-8',errors='surrogatepass')
    return {'text':raw[:limit].decode('utf-8',errors='ignore'),'utf8_bytes':len(raw),'sha256':sha(raw),'truncated':len(raw)>limit}

def error_record(error):
    identity=text_identity(error)
    return {'error_type':text_identity(type(error).__name__,64)['text'],'error':identity.pop('text'),'error_identity':identity}

def admit(path, expected):
    """Read bounded bytes once; hashing and parsing refer to this same object."""
    with Path(path).open('rb') as f: blob=f.read(expected['bytes']+1)
    need(len(blob)==expected['bytes'], 'source byte length differs: '+str(path))
    need(sha(blob)==expected['sha256'], 'source hash differs: '+str(path))
    return blob

def value_card(text):
    """FITS scalar lexical value; quoted slash/doubled quotes are retained."""
    s=text.lstrip()
    if not s.startswith("'"): return s.split('/',1)[0].strip()
    result=[];i=1
    while i<len(s):
        if s[i]=="'":
            if i+1<len(s) and s[i+1]=="'": result.append("'");i+=2;continue
            need(not s[i+1:].strip() or s[i+1:].lstrip().startswith('/'), 'junk after FITS string')
            return ''.join(result).strip()
        result.append(s[i]);i+=1
    raise Refusal('unterminated FITS string')

class GzipReader:
    def __init__(self, blob, maximum):
        self.f=gzip.GzipFile(fileobj=io.BytesIO(blob),mode='rb')
        self.maximum=maximum;self.count=0;self.hash=hashlib.sha256()
    def read(self,n):
        need(n>=0 and self.count+n<=self.maximum, 'inflated source limit')
        try: b=self.f.read(n)
        except (OSError,EOFError) as e: raise Refusal('gzip integrity failure: '+str(e)) from e
        self.count+=len(b);self.hash.update(b);return b
    def exact(self,n):
        b=self.read(n);need(len(b)==n, 'truncated FITS payload');return b
    def discard(self,n):
        while n:
            step=min(n,65536);self.exact(step);n-=step
    def end(self):
        need(not self.read(1), 'unexpected extra FITS/gzip payload')
        self.f.close()

def header(stream, maximum_cards):
    fields={};raw=[]
    for _ in range(maximum_cards):
        b=stream.exact(80)
        try:s=b.decode('ascii')
        except UnicodeError as e:raise Refusal('non-ASCII FITS header') from e
        raw.append(s);key=s[:8].strip()
        if key=='END':
            need(not s[8:].strip(), 'malformed FITS END')
            padding=stream.exact((-len(raw)*80)%2880)
            need(all(x==32 for x in padding),'nonblank FITS header padding')
            return fields,raw
        if key in {'','COMMENT','HISTORY','HIERARCH'}:continue
        need(s[8:10]=='= ', 'nonvalue FITS card: '+key)
        need(key not in fields,'duplicate FITS key: '+key)
        fields[key]=value_card(s[10:])
    raise Refusal('FITS header exceeds card limit/no END')

def table_header(stream, expected, maximum_cards):
    primary,raw_primary=header(stream,maximum_cards)
    extension,raw_extension=header(stream,maximum_cards)
    need(primary==expected[0], 'primary FITS schema differs')
    need(extension==expected[1], 'binary table FITS schema differs')
    need(extension.get('XTENSION')=='BINTABLE','not binary table')
    need(extension.get('PCOUNT')=='0' and extension.get('GCOUNT')=='1','unsupported heap/group')
    columns=[];offset=0
    for i in range(1,int(extension['TFIELDS'])+1):
        fmt=extension['TFORM'+str(i)];kind=fmt[-1]
        count=fmt[:-1]
        need(count.isdigit() and 1<=int(count)<=32,'unsupported FITS repeat')
        count=int(count)
        need(kind=='A' or count==1,'numeric vector not admitted')
        need(kind in {'A','E','D','I','J','K'},'unsupported FITS column kind')
        width=count*{'A':1,'E':4,'D':8,'I':2,'J':4,'K':8}[kind]
        name=extension['TTYPE'+str(i)]
        need(name not in {c['name'] for c in columns},'duplicate FITS column')
        columns.append({'name':name,'format':fmt,'offset':offset,'bytes':width,'source_TUNIT':extension.get('TUNIT'+str(i)),'qualified_unit':None})
        offset+=width
    need(offset==int(extension['NAXIS1']),'FITS row width differs')
    return {'primary':primary,'extension':extension,'raw_primary_cards':raw_primary,'raw_extension_cards':raw_extension,'columns':columns,'row_bytes':offset,'rows':int(extension['NAXIS2'])}

def decode(raw, columns):
    row={}
    for c in columns:
        b=raw[c['offset']:c['offset']+c['bytes']];kind=c['format'][-1]
        if kind=='A':
            try:v=ascii_field(b)
            except UnicodeError as e:raise Refusal('non-ASCII FITS row string') from e
        else:
            v=struct.unpack('>'+{'E':'f','D':'d','I':'h','J':'i','K':'q'}[kind],b)[0]
            need(not isinstance(v,float) or math.isfinite(v), 'nonfinite FITS field')
        row[c['name']]=v
    return row

def ascii_field(blob):
    value=blob.decode('ascii').rstrip(' \x00').lstrip(' ')
    need('\x00' not in value,'interior NUL in FITS string')
    return value

def raw_strings(raw,columns):
    return {c['name']:raw[c['offset']:c['offset']+c['bytes']].hex() for c in columns if c['format'].endswith('A')}

def text_lines(blob,bounds):
    try:lines=blob.decode('ascii').splitlines()
    except UnicodeError as e:raise Refusal('non-ASCII source text') from e
    need(len(lines)<=bounds['maximum_text_rows']+1,'text row limit')
    need(all(len(line)<=bounds['maximum_text_line_bytes'] for line in lines),'text line limit')
    return lines
def final_rows(blob,contract):
    lines=text_lines(blob,contract['bounds']);need(bool(lines),'empty final table')
    columns=lines[0].split();need(columns==contract['final_columns'],'final column schema differs')
    rows=[]
    for line in lines[1:]:
        if not line.strip():continue
        values=line.split();need(len(values)==len(columns),'final row width differs')
        rows.append(dict(zip(columns,values)))
    need(len(rows)==contract['selection']['release_rows'],'final row count differs')
    return rows
def update_rows(blob,contract):
    lines=text_lines(blob,contract['bounds']);reader=csv.reader(lines)
    need(next(reader,None)==contract['update_columns'],'update column schema differs')
    rows=[]
    for values in reader:
        if not values:continue
        need(len(values)==len(contract['update_columns']),'update row width differs')
        rows.append(dict(zip(contract['update_columns'],values)))
    return rows
def overrides(blob,key,bounds):
    lines=text_lines(blob,bounds)
    need(lines and lines[0].split()==['VARNAMES:','CID',key],'override schema differs')
    rows=[]
    for line in lines[1:]:
        if not line.strip():continue
        fields=line.split();need(len(fields)==3 and fields[0]=='SN:','override row schema differs')
        rows.append({'CID':fields[1],'value':fields[2]})
    return rows
def join_id(value):
    """Source ledger's literal leading-zero serialization rule, not event dedup."""
    return value.lstrip('0') or '0'
def index(rows,key,normalize=False):
    result={}
    for i,r in enumerate(rows):
        k=join_id(r[key]) if normalize else r[key]
        result.setdefault(k,[]).append((i,r))
    return result
def unique(indexed,key,label):
    matches=indexed.get(key,[]);need(len(matches)==1,label+' source match count '+str(len(matches))+' for '+key)
    return matches[0]
def pointer(lo,hi,nobs,total,maximum=None):
    need(isinstance(lo,int) and isinstance(hi,int) and isinstance(nobs,int),'noninteger PHOT pointers')
    need(1<=lo<=hi<=total and hi-lo+1==nobs,'PHOT pointer/NOBS mismatch')
    if maximum is not None:need(nobs<=maximum,'event record limit')
def exact_decimal_difference(left,right):
    try:
        with decimal.localcontext() as context:
            context.prec=128;context.traps[decimal.Inexact]=True;context.traps[decimal.Rounded]=True
            a,b=D(left),D(right);need(a.is_finite() and b.is_finite(),'nonfinite decimal coordinate')
            return a-b
    except decimal.DecimalException as e:raise Refusal('decimal difference not exactly admitted') from e
def coordinate(source_literal,final_literal):
    s=D(source_literal);v=D(final_literal)
    need(s.is_finite() and v.is_finite(),'nonfinite coordinate')
    with decimal.localcontext() as context:
        context.prec=128;context.traps[decimal.Inexact]=True;context.traps[decimal.Rounded]=True
        half=D(5).scaleb(v.as_tuple().exponent-1)
        difference=exact_decimal_difference(source_literal,final_literal)
        return {'source_literal':source_literal,'final_literal':final_literal,'difference_degree':str(difference),'printed_half_unit_degree':str(half),'within_printed_bound':abs(difference)<=half}

def head_table(blob,contract):
    stream=GzipReader(blob,contract['bounds']['maximum_head_uncompressed_bytes'])
    meta=table_header(stream,contract['fits']['sdss-header-example.dat'],contract['bounds']['maximum_header_cards'])
    need(meta['rows']==contract['selection']['header_rows'],'HEAD row count')
    body=stream.exact(meta['rows']*meta['row_bytes']);stream.discard((-len(body))%2880);stream.end()
    snid=meta['columns'][0];need(snid['name']=='SNID' and snid['format']=='16A','HEAD SNID schema')
    ids={}
    for i in range(meta['rows']):
        raw=body[i*meta['row_bytes']:i*meta['row_bytes']+16]
        try:k=ascii_field(raw)
        except UnicodeError as e:raise Refusal('non-ASCII HEAD SNID') from e
        need(bool(k),'empty HEAD SNID');ids.setdefault(k,[]).append(i)
    meta.update({'inflated_bytes':stream.count,'inflated_sha256':stream.hash.hexdigest()})
    return body,meta,ids

def phot_slice(blob,name,lo,hi,contract):
    stream=GzipReader(blob,contract['bounds']['maximum_phot_uncompressed_bytes'])
    meta=table_header(stream,contract['fits'][name],contract['bounds']['maximum_header_cards'])
    need(meta['rows']==contract['selection']['photometry_rows'],'PHOT row count')
    pointer(lo,hi,hi-lo+1,meta['rows'],contract['bounds']['maximum_event_records'])
    size=meta['row_bytes'];stream.discard((lo-1)*size);records=[];whole=hashlib.sha256()
    for i in range(lo,hi+1):
        raw=stream.exact(size);whole.update(raw)
        records.append({'row_one_based':i,'row_zero_based':i-1,'binary_record_sha256':sha(raw),'raw_string_storage_hex':raw_strings(raw,meta['columns']),'decoded':decode(raw,meta['columns'])})
    stream.discard((meta['rows']-hi)*size)
    stream.discard((-(meta['rows']*size))%2880);stream.end()
    meta.update({'inflated_bytes':stream.count,'inflated_sha256':stream.hash.hexdigest(),'slice_sha256':whole.hexdigest()})
    return {'metadata':meta,'records':records}

def regenerate(contract,blobs):
    final=final_rows(blobs['sn-table.dat'],contract);updates=update_rows(blobs['all_redshifts_PVs.csv'],contract)
    ui=index(updates,'SNID',True);oi={key:index(overrides(blobs[key+'.txt'],key,contract['bounds']),'CID',True) for key in ['REDSHIFT_CMB','REDSHIFT_CMB_ERR','VPEC','VPEC_ERR']}
    body,headmeta,hi=head_table(blobs['sdss-header-example.dat'],contract)
    selected=[(i,r) for i,r in enumerate(final) if r['IDSURVEY']==contract['selection']['IDSURVEY']]
    need(len(selected)==contract['selection']['selected_rows'],'selected row count differs')
    need(len({r['CID'] for _,r in selected})==len(selected),'duplicate selected light-curve CID')
    events=[];failures=[]
    for i,r in selected:
        matches=hi.get(r['CID'],[]);need(len(matches)==1,'HEAD source match count '+str(len(matches))+' for '+r['CID'])
        head_index=matches[0];raw=body[head_index*headmeta['row_bytes']:(head_index+1)*headmeta['row_bytes']];h=decode(raw,headmeta['columns'])
        pointer(h['PTROBS_MIN'],h['PTROBS_MAX'],h['NOBS'],contract['selection']['photometry_rows'])
        update_index,u=unique(ui,join_id(r['CID']),'update')
        source_overrides={}
        for key,field in [('REDSHIFT_CMB','zCMB'),('REDSHIFT_CMB_ERR','zCMBERR'),('VPEC','VPEC'),('VPEC_ERR','VPECERR')]:
            oirow,ov=unique(oi[key],join_id(r['CID']),key)
            source_overrides[key]={'source_row_zero_based':oirow,'raw':ov,'release_field':field,'source_minus_final_exact_decimal':str(exact_decimal_difference(ov['value'],r[field]))}
        coord={axis:coordinate(u[sfield],r[axis]) for axis,sfield in [('RA','RA'),('DEC','Dec')]}
        head_coord={axis:coordinate(str(D.from_float(h[sfield])),r[axis]) for axis,sfield in [('RA','RA'),('DEC','DECL')]}
        for axis,c in coord.items():
            if not c['within_printed_bound']:failures.append({'CID':r['CID'],'release_row_zero_based':i,'axis':axis,**c})
        events.append({'release_row_zero_based':i,'CID':r['CID'],'header_row_zero_based':head_index,'update_row_zero_based':update_index,'header_binary_record_sha256':sha(raw),'head_raw_string_storage_hex':raw_strings(raw,headmeta['columns']),'PHOTFILE':headmeta['primary']['PHOTFILE'],'photometry_one_based_range':[h['PTROBS_MIN'],h['PTROBS_MAX']],'NOBS':h['NOBS'],'IAUC':h['IAUC'],'final':r,'head':h,'update':u,'overrides':source_overrides,'update_final_coordinates':coord,'earlier_head_final_coordinates':head_coord,'physical_event_identity_qualified':False})
    need(sum(e['NOBS'] for e in events)==contract['selection']['linked_photometry_points'],'coupled PHOT point count differs')
    chosen=[e for e in events if e['CID']==contract['selection']['event_CID']];need(len(chosen)==1,'chosen event not unique');event=chosen[0]
    for k in ['release_row_zero_based','header_row_zero_based','update_row_zero_based']:need(event[k]==contract['selection']['event_'+k],'chosen event index differs')
    need(event['photometry_one_based_range']==contract['selection']['event_photometry_one_based'] and event['NOBS']==contract['selection']['event_records'],'chosen event pointers differ')
    lo,up=event['photometry_one_based_range'];named=event['PHOTFILE']+'.gz'
    need(named=='SDSS_allCandidates+BOSS_PHOT.FITS.gz','HEAD PHOTFILE differs')
    slices={n:phot_slice(blobs[n],n,lo,up,contract) for n in [named,'JLA2014_SDSS_DS17_PHOT.FITS.gz']}
    a,b=slices.values();need(a['records']==b['records'],'distinct PHOT sources differ in chosen event records')
    negative_flux_records=sum(r['decoded']['FLUXCAL']<0 for r in a['records'])
    expected=contract['structural_expectations']
    need([(f['CID'],f['axis']) for f in failures]==[(cid,'RA') for cid in expected['retained_RA_failure_CIDs_in_release_order']],'retained coordinate failure identities differ')
    need(negative_flux_records==expected['signed_negative_flux_records_in_CID6057_slice'],'signed source flux record count differs')
    return {'role':'calibrated/fitted source lineage; no physical conversion, fit or inference','source_revision':contract['source_revision'],'candidate_sha256':contract['candidate_sha256'],'source_gates_origin':contract['source_gates_origin'],'HEAD_metadata':headmeta,'events':events,'coordinate_failures':failures,'selected_photometry_points':sum(e['NOBS'] for e in events),'single_event':event,'photometry_slices':slices,'signed_negative_flux_records':negative_flux_records,'physical_units':contract['semantic_nulls'],'unit_blocker':'FITS TUNIT cards present but blank; no EXPTIME column. No clock, duration, area or flux conversion inferred.','selection':'IDSURVEY1 only; all 321 rows and all 115 chosen records retained, including all signed fluxes, source flags and sentinels. No detection threshold inferred.','covariance':'Release row order retained; no covariance subset, independence or unknown cross-covariance substitution.','source_status':'coordinate failures retained; frame/reduction/event qualification remains blocked'}

def fresh_output(path,root):
    target=Path(path).absolute();root=Path(root).resolve()
    need(target.parent==root and re.fullmatch(r'[a-z0-9][a-z0-9-]{0,63}',target.name),'output must be a bounded single-component attempt in the packet result directory')
    need(not target.exists() and not target.is_symlink(),'output already exists; preserve prior attempt')
    return target

def packet_output_root(root):
    root=Path(root).resolve()
    results=root/'results';destination=results/'sdss-released-observer-contract'
    for path in [results,destination]:
        need(not path.is_symlink(),'symlink result directory refused')
        need(not path.exists() or path.is_dir(),'result parent is not a directory')
        need(path.resolve()==path,'result directory escapes repository root')
    results.mkdir(exist_ok=True)
    destination.mkdir(exist_ok=True)
    return destination

def bounded_identity(path,limit,git_blob=False):
    display=text_identity(path)
    try:
        with Path(path).open('rb') as f:b=f.read(limit+1)
        if len(b)>limit:return {'path':display,'sha256':None,'error':'identity byte limit'}
        result={'path':display,'bytes':len(b),'sha256':sha(b)}
        if git_blob:result['git_blob_sha1']=hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()
        return result
    except (OSError,ValueError) as e:return {'path':display,'sha256':None,**error_record(e)}

def require_identities(identities,label):
    for key,value in identities.items():
        need('error' not in value and re.fullmatch(r'[0-9a-f]{64}',value.get('sha256') or '') and isinstance(value.get('bytes'),int) and value['bytes']>0,label+' identity not admitted: '+key)

def source_identities(input_root,sources):
    return {name:bounded_identity(input_root/name,expected['bytes']) for name,expected in sources.items()}

def source_identity_errors(identities,sources):
    return [name+' original source identity differs from pin' for name,expected in sources.items() if 'error' in identities[name] or identities[name].get('bytes')!=expected['bytes'] or identities[name].get('sha256')!=expected['sha256']]

def source_control_identity(root,identity_paths):
    """Read-only Git source gate; ignored results are outside committed inputs."""
    try:
        paths={name:str(path.relative_to(root)) for name,path in identity_paths.items()}
        def query(arguments,maximum):
            process=subprocess.run(['git','-C',str(root),*arguments],stdout=subprocess.PIPE,stderr=subprocess.DEVNULL,timeout=10,check=False)
            need(process.returncode==0,'Git identity command failed: '+arguments[0])
            need(len(process.stdout)<=maximum,'Git identity output limit')
            return process.stdout.decode('ascii').strip()
        top=query(['rev-parse','--show-toplevel'],4096)
        need(Path(top).resolve()==root.resolve(),'Git root differs from packet ROOT')
        revision=query(['rev-parse','--verify','HEAD'],128)
        need(re.fullmatch(r'[0-9a-f]{40}',revision),'invalid committed Git revision')
        query(['diff','--quiet','--no-ext-diff','--no-textconv','HEAD','--'],0)
        entries=query(['ls-files','--stage','--',*paths.values()],8192).splitlines()
        committed={}
        for line in entries:
            metadata,path=line.split('\t',1);mode,blob,stage=metadata.split()
            need(path in paths.values() and path not in committed and mode in {'100644','100755'} and stage=='0' and re.fullmatch(r'[0-9a-f]{40}',blob),'packet Git index identity differs')
            committed[path]={'mode':mode,'git_blob_sha1':blob}
        need(set(committed)==set(paths.values()),'all five packet files must be committed')
        actual={path:bounded_identity(root/path,131072,git_blob=True) for path in paths.values()}
        need(all('error' not in value and value.get('git_blob_sha1')==committed[path]['git_blob_sha1'] for path,value in actual.items()),'packet bytes differ from committed Git blobs')
        return {'root':text_identity(root),'revision':revision,'tracked_checkout_clean':True,'untracked_ignored_storage':'excluded from committed source gate; five packet files must be tracked','packet_files':{name:{'path':path,**committed[path]} for name,path in paths.items()}}
    except (OSError,subprocess.SubprocessError,ValueError) as e:
        return {'root':text_identity(root),'revision':None,**error_record(e)}

def require_source_control(identity):
    need('error' not in identity and identity.get('tracked_checkout_clean') is True and re.fullmatch(r'[0-9a-f]{40}',identity.get('revision') or ''),'clean committed Reproducible source identity not admitted')

def runtime_identity():
    """Bounded identity of the actual Python executable used for this reader."""
    try:
        executable=Path(sys.executable).resolve(strict=True)
        identity=bounded_identity(executable,67108864)
    except (OSError,ValueError,RuntimeError) as e:identity={'path':text_identity(sys.executable),'sha256':None,**error_record(e)}
    identity['version']=text_identity(sys.version)
    identity['implementation']=text_identity(sys.implementation.name,64)
    return identity

def terminal_receipt(receipt):
    """Fixed-cardinality fallback preserves used hashes instead of raising on size."""
    encoded=(json.dumps(receipt,indent=2,allow_nan=False)+'\n').encode()
    if len(encoded)<=65536:return receipt,encoded
    def digest(value,width=64):
        return value if isinstance(value,str) and re.fullmatch('[0-9a-f]{'+str(width)+'}',value) else None
    def compact_records(records,maximum):
        result=[]
        for name,value in list(records.items())[:maximum]:
            path=value.get('path',{})
            count=value.get('bytes')
            result.append({'name':text_identity(name,32),'bytes':count if isinstance(count,int) and 0<=count<=134217728 else None,'sha256':digest(value.get('sha256')),'path_display':text_identity(path.get('text',''),48)['text'] if isinstance(path,dict) else None,'path_text_sha256':digest(path.get('sha256')) if isinstance(path,dict) else None,'had_error':'error' in value})
        return result
    def compact_git(identity):
        return {'revision':digest(identity.get('revision'),40),'had_error':'error' in identity,'packet_files':[{'name':text_identity(name,32),'git_blob_sha1':digest(value.get('git_blob_sha1'),40)} for name,value in list(identity.get('packet_files',{}).items())[:5]]}
    original_error=receipt.get('error_identity',{})
    fallback={'status':'failed_receipt_bound','previous_status':text_identity(receipt.get('status'),64),'scientific_execution':None,'error':'Full terminal receipt exceeded 65536 bytes; bounded identity summary retained.','full_receipt_bytes':len(encoded),'full_receipt_sha256':sha(encoded),'contract_sha256':digest(receipt.get('contract_sha256')),'lineage_sha256':digest(receipt.get('lineage_sha256')),'original_error_text_sha256':digest(original_error.get('sha256')),'output_bytes_before_receipt':receipt.get('output_bytes_before_receipt',0)}
    for key,maximum in [('identities_before',5),('identities_after',5),('original_sources_before',9),('original_sources_after',9),('admitted_sources',9)]:fallback[key]=compact_records(receipt.get(key,{}),maximum)
    for key in ['source_control_before','source_control_after']:fallback[key]=compact_git(receipt.get(key,{}))
    for key in ['runtime_before','runtime_after']:
        identity=receipt.get(key,{})
        fallback[key]={'executable':compact_records({'python':identity},1),'version_display':text_identity(identity.get('version',{}).get('text',''),128)['text'],'version_text_sha256':digest(identity.get('version',{}).get('sha256')),'implementation_text_sha256':digest(identity.get('implementation',{}).get('sha256'))}
    fallback['input_root_text_sha256']=digest(receipt.get('input_root',{}).get('sha256'))
    fallback['input_argument_text_sha256']=digest(receipt.get('input_argument',{}).get('sha256'))
    fallback['integrity_errors']=[text_identity(value,64) for value in receipt.get('integrity_errors',[])[:32]]
    return fallback,(json.dumps(fallback,indent=2,allow_nan=False)+'\n').encode()

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-root',required=True,help='Read-only directory containing the nine exact source basenames; argv at most 1024 UTF-8 bytes')
    parser.add_argument('--attempt',required=True,help='Fresh bounded result component; never reuses an attempt')
    args=parser.parse_args()
    need(re.fullmatch(r'[a-z0-9][a-z0-9-]{0,63}',args.attempt),'invalid attempt component')
    input_argument=text_identity(args.input_root)
    need(input_argument['utf8_bytes']<=1024,'input-root argv exceeds 1024 UTF-8 bytes; no attempt created')
    args.input_root=Path(args.input_root)
    input_root=text_identity(args.input_root.absolute())
    need(input_root['utf8_bytes']<=4096,'absolute input-root exceeds 4096 UTF-8 bytes; no attempt created')
    args.contract=FOLDER/'input-contract.json'
    args.output=fresh_output(packet_output_root(ROOT)/args.attempt,ROOT/'results/sdss-released-observer-contract')
    args.output.mkdir();receipt={'status':'started','scientific_execution':None,'operations':'source-only structural reads, hashes and joins'}
    identity_paths={'reader':Path(__file__).resolve(),'contract':args.contract,'candidate':FOLDER/'candidate.json','experiment':FOLDER/'experiment.json','README':FOLDER/'README.md'}
    receipt['identities_before']={k:bounded_identity(p,131072) for k,p in identity_paths.items()}
    receipt['source_control_before']=source_control_identity(ROOT,identity_paths)
    receipt['runtime_before']=runtime_identity()
    receipt['input_argument']=input_argument;receipt['input_root']=input_root
    output_limit=4194304;sources=None
    try:
        require_identities(receipt['identities_before'],'packet')
        require_source_control(receipt['source_control_before'])
        require_identities({'python':receipt['runtime_before']},'runtime')
        with args.contract.open('rb') as f:cb=f.read(65537)
        need(len(cb)<=65536 and sha(cb)==EXPECTED_CONTRACT_SHA256,'reviewed input contract differs')
        contract=json.loads(cb);receipt['contract_sha256']=sha(cb)
        need(contract['schema']==1,'unknown input contract')
        need(contract['candidate_path']=='candidate.json','candidate snapshot path differs')
        with (FOLDER/contract['candidate_path']).open('rb') as f:candidate=f.read(65537)
        need(len(candidate)<=65536 and sha(candidate)==contract['candidate_sha256'],'pinned candidate differs')
        need(sum(x['bytes'] for x in contract['sources'].values())<=contract['bounds']['maximum_compressed_total_bytes'],'total admitted byte limit')
        need(1<=len(contract['sources'])<=9,'original source identity cardinality exceeds reviewed pins')
        need(contract['bounds']['maximum_output_bytes']==output_limit,'output byte allocation differs')
        sources=contract['sources']
        receipt['original_sources_before']=source_identities(args.input_root,sources)
        require_identities(receipt['original_sources_before'],'original source')
        need(not source_identity_errors(receipt['original_sources_before'],sources),'original source before identity differs from pin')
        blobs={name:admit(args.input_root/name,expected) for name,expected in contract['sources'].items()}
        receipt['admitted_sources']={n:{'bytes':len(b),'sha256':sha(b)} for n,b in blobs.items()}
        result=regenerate(contract,blobs)
        encoded=(json.dumps(result,indent=2,allow_nan=False)+'\n').encode()
        need(len(encoded)<=output_limit-65536,'lineage output byte limit')
        (args.output/'lineage.json').write_bytes(encoded)
        receipt.update({'status':'structural_completed_with_retained_coordinate_failures','lineage_sha256':sha(encoded),'selected_rows':len(result['events']),'coordinate_failures':len(result['coordinate_failures']),'single_event_records':len(next(iter(result['photometry_slices'].values()))['records']),'signed_negative_flux_records':result['signed_negative_flux_records'],'no_physical_qualification':True})
    except Exception as e:
        receipt.update({'status':'failed',**error_record(e)})
    finally:
        receipt['identities_after']={k:bounded_identity(p,131072) for k,p in identity_paths.items()}
        receipt['source_control_after']=source_control_identity(ROOT,identity_paths)
        receipt['runtime_after']=runtime_identity()
        receipt['integrity_errors']=[k+' identity changed' for k in identity_paths if receipt['identities_before'][k]!=receipt['identities_after'][k]]
        receipt['integrity_errors'] += [k+' packet identity error' for k,v in receipt['identities_after'].items() if 'error' in v or not v.get('sha256')]
        if receipt['source_control_before']!=receipt['source_control_after']:receipt['integrity_errors'].append('committed source identity changed')
        if 'error' in receipt['source_control_after']:receipt['integrity_errors'].append('committed source identity error')
        if receipt['runtime_before']!=receipt['runtime_after']:receipt['integrity_errors'].append('Python executable/version identity changed')
        if 'error' in receipt['runtime_after'] or not receipt['runtime_after'].get('sha256'):receipt['integrity_errors'].append('Python executable identity error')
        if sources is not None:
            receipt['original_sources_after']=source_identities(args.input_root,sources)
            receipt['integrity_errors'] += source_identity_errors(receipt['original_sources_after'],sources)
            receipt['integrity_errors'] += [n+' original source identity changed' for n in sources if receipt['original_sources_before'][n]!=receipt['original_sources_after'][n]]
        if receipt['integrity_errors']:receipt['status']='failed_identity_changed'
        receipt['output_bytes_before_receipt']=sum(p.stat().st_size for p in args.output.iterdir() if p.is_file())
        receipt,rb=terminal_receipt(receipt)
        (args.output/'receipt.json').write_bytes(rb)
        for p in args.output.iterdir():p.chmod(0o444)
        args.output.chmod(0o555)
    print(json.dumps(receipt))
    if receipt['status'].startswith('failed'):raise SystemExit(1)
if __name__=='__main__':main()
