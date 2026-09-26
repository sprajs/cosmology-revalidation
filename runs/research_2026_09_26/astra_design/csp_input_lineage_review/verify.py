"""Independent original-archive and byte-level preparation review. No executor imports/fits."""
from pathlib import Path
from decimal import Decimal,localcontext
import csv,tarfile,hashlib,json,re,collections
ROOT=Path.cwd();OUT=Path(__file__).parent
AR=ROOT/'runs/research_2026_09_26/csp_dr3_provenance/CSP_Photometry_DR3.tgz'
SRC=ROOT/'sources/repos/djones1040__RAISIN_DataRelease@a383c4b/photometry/RAISIN/CSPDR3_RAISIN'
PREP=ROOT/'runs/research_2026_09_26/csp_native_filter_response/cohort-preparation'
DESIGN=ROOT/'runs/research_2026_09_26/raisin_profile_solver_review/csp-filter-interpretation/native-design'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def csvread(p):
 with p.open(newline='') as f:return list(csv.DictReader(f))
def dump(n,x):(OUT/n).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def writecsv(n,x):
 with (OUT/n).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(x[0]));w.writeheader();w.writerows(x)
def norm(x):return x.lower().removeprefix('sn')
mapmerge={'Jdw':'J','Hdw':'H'};maprelease={'J':'J','Jrc2':'j','Y':'Y','Ydw':'y','H':'H'}
NIR={'J','Jrc2','Jdw','H','Hdw','Y','Ydw'}
readpaths={AR,Path(__file__),PREP/'protocol.json',DESIGN/'cohort42.csv'}
prep=json.loads((PREP/'protocol.json').read_text());failed_hash=[]
for name,h in prep['inputs_sha256'].items():
 p=ROOT/name;readpaths.add(p)
 if sha(p)!=h:failed_hash.append(name)
assert not failed_hash,failed_hash
assert sha(AR)=='ea337b375a7da6223b11f2dd6b50e8a4546045c3133719cc66cbf1bdf4795026'
raw=[];snpy=[];presence=set()
with tarfile.open(AR) as t:
 for lineno,line in enumerate(t.extractfile('DR3/SN_photo.dat').read().decode().splitlines(),1):
  a=line.split()
  if a[1] in NIR:raw.append((norm(a[0]),a[1],Decimal(a[2]),Decimal(a[3]),Decimal(a[4])))
 for member in t:
  if not(member.isfile() and member.name.endswith('_snpy.txt')):continue
  lines=t.extractfile(member).read().decode().splitlines();name=norm(lines[0].split()[0]);assert name not in presence;presence.add(name);band=None
  for line in lines[1:]:
   a=line.split()
   if not a or a[0].startswith('#'):continue
   if a[0]=='filter':band=a[1];continue
   if band in NIR:snpy.append((name,band,Decimal(a[0]),Decimal(a[1]),Decimal(a[2])))
rawcounter=collections.Counter((n,mapmerge.get(b,b),t,m,e) for n,b,t,m,e in raw);sncounter=collections.Counter(snpy)
assert rawcounter==sncounter
# Only metadata determine physical candidate labels; magnitudes never resolve mixed labels.
phys={}
for n,b,t,m,e in raw:
 rb=maprelease[mapmerge.get(b,b)];phys.setdefault((n,rb,t),set()).add(b)
def released(p):
 readpaths.add(p);rows=[];name=None;vars=None
 for lineidx,line in enumerate(p.read_text().splitlines(),1):
  a=line.split()
  if not a:continue
  if a[0]=='SNID:':name=norm(a[1])
  if a[0]=='VARLIST:':vars=a[1:]
  if a[0]=='OBS:':
   d=dict(zip(vars,a[1:],strict=True))
   if d['FLT'] not in ['J','j','Y','y','H']:continue
   rows.append({'name':name,'band':d['FLT'],'t':Decimal(d['MJD'])-53000,'m':Decimal(d['MAG']),'e':Decimal(d['MAGERR']),'flux':Decimal(d['FLUXCAL']),'fluxerr':Decimal(d['FLUXCALERR']),'line':lineidx})
 assert name is not None
 return name,rows
def expected(n):
 out=[]
 for nn,b,t,m,e in snpy:
  if nn!=n:continue
  with localcontext() as c:
   c.prec=50;f=Decimal(10)**(-Decimal('.4')*(m-Decimal('27.5')));fe=f*(Decimal(10)**(Decimal('.4')*e)-1)
   # Compare published numeric precision, not spelling of exponent signs.
   out.append((maprelease[b],t,Decimal(format(m,'.3f')),Decimal(format(e,'.3f')),Decimal(format(f,'.5e')),Decimal(format(fe,'.5e'))))
 return collections.Counter(out)
listing=SRC/'CSPDR3_RAISIN.LIST';readpaths.add(listing);names=[line.split()[0] for line in listing.read_text().splitlines() if line.strip() and not line.lstrip().startswith('#')]
assert len(names)==len(set(names))==76
objectrows=[];physicalrows=[];sourcepaths={};sourceNIR={}
for fn in names:
 path=SRC/fn;name,rr=released(path);sourcepaths[name]=path;sourceNIR[name]=rr;want=expected(name);actual=collections.Counter((r['band'],r['t'],r['m'],r['e'],r['flux'],r['fluxerr']) for r in rr)
 found=name in presence;hasNIR=any(x[0]==name for x in snpy)
 if found:assert want==actual,name
 objectrows.append({'CID':name,'archive_present':found,'archive_has_NIR':hasNIR,'release_NIR_rows':len(rr),'expected_NIR_rows':sum(want.values()),'full_decimal_magnitude_error_and_flux_identity':found and want==actual})
 for r in rr:
  candidates=phys.get((name,r['band'],r['t']),set());physicalrows.append({'CID':name,'source_line':r['line'],'band':r['band'],'time':str(r['t']),'physical_candidates':'|'.join(sorted(candidates)),'physical_unique':len(candidates)==1,'decimal_converter_identity':found and want==actual})
writecsv('independent-objects.csv',objectrows);writecsv('independent-physical-rows.csv',physicalrows)
unique={b:sum(r['physical_candidates']==b and r['physical_unique'] for r in physicalrows) for b in ['Jdw','Hdw']}
assert unique=={'Jdw':303,'Hdw':288},unique
cohort=csvread(DESIGN/'cohort42.csv');ids={norm(r['CID']) for r in cohort};assert len(ids)==42
selected79=ROOT/'runs/research_2026_09_26/raisin_differential/frozen-membership.csv';readpaths.add(selected79)
published=csvread(selected79);pubcsp={norm(r['CID']) for r in published if r['survey']=='CSP'};assert ids==pubcsp
changed=[];perfile=[];ambiguities=[];controls=[]
for scope,wantids in [('pilot',{'2004ef','2005hc'}),('full',ids)]:
 dirs={arm:PREP/scope/'data'/arm/'CSPDR3_RAISIN' for arm in ['nominal','nominal_copy','known_Jdw_to_j']}
 lists={arm:(di/'CSPDR3_RAISIN.LIST').read_bytes() for arm,di in dirs.items()};assert len(set(lists.values()))==1
 preparednames=[line.split()[0] for line in lists['nominal'].decode().splitlines() if line.strip() and not line.lstrip().startswith('#')]
 assert len(preparednames)==len(set(preparednames))==len(wantids)
 actualids=set()
 for fn in preparednames:
  src=SRC/fn;name,rr=released(src);actualids.add(name);orig=src.read_bytes();lines=orig.splitlines(keepends=True);out=list(lines);changes=0
  for r in rr:
   cand=phys.get((name,r['band'],r['t']),set())
   if r['band']=='J' and cand=={'Jdw'}:
    i=r['line']-1;tokens=list(re.finditer(rb'\S+',lines[i]));b=tokens[2];assert b.group()==b'J';out[i]=lines[i][:b.start()]+b'j'+lines[i][b.end():];changes+=1;changed.append({'scope':scope,'CID':name,'source_line':r['line'],'MJD':str(r['t']+53000),'old_band':'J','new_band':'j'})
   elif scope=='full' and r['band']=='J' and len(cand)>1:
    ambiguities.append({'CID':name,'source_line':r['line'],'MJD':str(r['t']+53000),'labels':'|'.join(sorted(cand)),'unchanged':True})
  expectedbytes=b''.join(out)
  for arm,di in dirs.items():
   f=di/fn;readpaths.add(f);wantbytes=expectedbytes if arm=='known_Jdw_to_j' else orig;assert f.read_bytes()==wantbytes,(scope,arm,name)
  perfile.append({'scope':scope,'CID':name,'changed_J_tokens':changes,'all_other_bytes_identical':True,'nominal_copy_source_byte_identical':True})
  if scope=='full' and changes==0:controls.append(name)
 assert actualids==wantids
 # Prepared fitting configurations differ only in their explicit data-arm path.
 nmls=[]
 for arm in dirs:
  f=PREP/scope/'fits'/arm/'fit.nml';readpaths.add(f);text=f.read_text();normtext=re.sub(r"PRIVATE_DATA_PATH\s*=\s*'[^']*'","PRIVATE_DATA_PATH = '<arm>'",text);nmls.append(normtext)
 assert len(set(nmls))==1
assert len(controls)==10 and len(ambiguities)==4 and all(r['CID']=='2006kf' for r in ambiguities)
assert sum(r['changed_J_tokens'] for r in perfile if r['scope']=='full')==188
assert sum(r['changed_J_tokens'] for r in perfile if r['scope']=='pilot')==6
rootledger=csvread(PREP/'changed-row-ledger.csv')
key=lambda r:(r['scope'],norm(r['CID']),int(r['source_line']),Decimal(r['MJD']),r['old_band'],r['new_band'])
assert collections.Counter(map(key,rootledger))==collections.Counter(map(key,changed))
writecsv('independent-changed-tokens.csv',changed);writecsv('independent-file-gates.csv',perfile);writecsv('independent-ambiguous-unchanged.csv',ambiguities)
rootresult=ROOT/'runs/research_2026_09_26/csp_magnitude_lineage/result.json';readpaths.add(rootresult);v=json.loads(rootresult.read_text())
assert v['raw_NIR_rows']==len(raw) and v['SNpy_NIR_rows']==len(snpy)
assert len(raw)==len(snpy)==5491
found=[r for r in objectrows if r['archive_present']];assert len(found)==73
noNIR=[r['CID'] for r in found if not r['archive_has_NIR']];assert len(noNIR)==8
missing=[r['CID'] for r in objectrows if not r['archive_present']];assert missing==['2012fr','2012ht','2015f']
assert sum(r['release_NIR_rows'] for r in objectrows)==v['release_NIR_rows']==3598
for n in ['v1-before-presence-flag-fix','v2-before-no-NIR-presence-fix']:assert (rootresult.parent/n/'result.json').exists()
res={'pass':True,'scope':'Independent original-archive/magnitude/flux lineage and prepared-byte intervention verification; no native fit or physical calibration/distance inference.','imports_root_executor':False,'verified_preparation_input_hashes':len(prep['inputs_sha256']),'raw_NIR_rows':len(raw),'SNpy_NIR_rows':len(snpy),'full_relabelled_decimal_multiset_identity':True,'listed_RAISIN_objects':len(names),'archive_present_objects':len(found),'archive_present_with_NIR':len(found)-len(noNIR),'archive_present_without_NIR':noNIR,'added_CSP_II_absent_from_archive':missing,'all_archive_present_released_magnitude_error_flux_rows_exact':True,'physical_unique_counts':unique,'cohort_equals_published42CSP':True,'full_changed_tokens':188,'pilot_changed_tokens':6,'full_unchanged_controls':controls,'ambiguous_rows_retained':ambiguities,'all_prepared_nominal_and_copies_byte_identical_to_source':True,'changed_arm_exactly_reconstructed_from_source_metadata':True,'prepared_NMLs_only_arm_data_path_differs':True,'changed_ledger_exact':True,'v1_v2_reporting_revisions_retained':True}
dump('result.json',res);dump('input-manifest.json',{'sha256':{str(p.relative_to(ROOT)):sha(p) for p in sorted(readpaths)}});print(json.dumps(res,indent=2))
