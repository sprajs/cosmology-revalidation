"""Frozen all-row comparison of original CSP tables, SNpy and RAISIN products."""
from pathlib import Path
from decimal import Decimal as D
from collections import Counter, defaultdict
import csv
import hashlib
import json
import sys
import tarfile

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'runs/research_2026_09_26/csp_optical_magnitude_lineage'
ARCHIVE = ROOT/'runs/research_2026_09_26/csp_dr3_provenance/CSP_Photometry_DR3.tgz'
DATA = ROOT/'sources/repos/djones1040__RAISIN_DataRelease@a383c4b/photometry/RAISIN/CSPDR3_RAISIN'
OPTICAL = {'u','g','r','i','B','V','V0','V1'}
MERGE = {'V0':'V'}
MAP = {'u':'u','g':'g','r':'r','i':'i','B':'B','V':'n','V0':'m','V1':'o'}

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,x): p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def norm(x): return str(D(x).normalize())

def prepare():
    OUT.mkdir(exist_ok=False)
    files=[ARCHIVE,Path(__file__),DATA/'CSPDR3_RAISIN.LIST']
    files += [DATA/n for n in (DATA/'CSPDR3_RAISIN.LIST').read_text().splitlines() if n.strip() and not n.startswith('#')]
    dump(OUT/'protocol.json',dict(
        status='Frozen before numerical magnitude/error comparisons; metadataV0→V collapse known; physical transformation history not yet resolved. Separate optical test after NIR outcome.',
        scope='Every OPTICAL row in original archive; all LISTed RAISIN objects; preserve duplicates, unmatched objects/rows and any mismatches. No distance fits or correction inferred from identity alone.',
        raw_to_snpy='Exact Decimal multiset of name, mapped band, original time, magnitude, error. V0→V. No tolerance or matching on brightness.',
        snpy_to_raisin='For each LISTed object match exact name/band/time (+53000 MJD). Compare all MAG/MAGERR and converter-predicted flux/error rounded to five decimal scientific formatting; retain missing and excess multiplicity. Do not assume historical converter execution.',
        flux='10**(-0.4*(mag-27.5)); error=flux*(10**(0.4*magerr)-1), matching tagged converter. Decimal mag/error exact; floating conversion compared after published-format rounding.',
        raw_physical_filter='Disambiguate only where metadata (name,time,mapped band) has one raw physical filter; no photometric matching to choose physical labels.',
        inputs_sha256={str(p.relative_to(ROOT)):sha(p) for p in files} ))
    print('protocol',sha(OUT/'protocol.json'))

def run():
    p=json.loads((OUT/'protocol.json').read_text())
    for n,d in p['inputs_sha256'].items(): assert sha(ROOT/n)==d,n
    raw=[];snpy=[];archive_names=set()
    with tarfile.open(ARCHIVE) as tar:
        for line in tar.extractfile('DR3/SN_photo.dat').read().decode().splitlines():
            a=line.split();assert len(a)==5
            if a[1] in OPTICAL:raw.append((a[0].lower(),a[1],norm(a[2]),norm(a[3]),norm(a[4])))
        for m in tar:
            if not(m.isfile() and m.name.endswith('_snpy.txt')):continue
            lines=tar.extractfile(m).read().decode().splitlines();name=lines[0].split()[0].lower();band=None
            archive_names.add(name)
            for l in lines[1:]:
                a=l.split()
                if not a or a[0].startswith('#'):continue
                if a[0]=='filter':band=a[1];continue
                assert len(a)==3
                if band in OPTICAL:snpy.append((name,band,*map(norm,a)))
    expected=Counter((c,MERGE.get(b,b),t,m,e) for c,b,t,m,e in raw)
    actual=Counter(snpy)
    rawmeta=defaultdict(set)
    for c,b,t,m,e in raw:rawmeta[(c,MERGE.get(b,b),t)].add(b)
    byobj=defaultdict(list)
    for r in snpy:byobj[r[0]].append(r)
    rowledger=[];objects=[];release=[]
    for fname in (DATA/'CSPDR3_RAISIN.LIST').read_text().splitlines():
        if not fname.strip() or fname.startswith('#'):continue
        hdr={};obs=[];cols=None
        for line in (DATA/fname).read_text().splitlines():
            a=line.split()
            if not a:continue
            if a[0]=='VARLIST:':cols=a[1:]
            elif a[0]=='OBS:':obs.append(dict(zip(cols,a[1:])))
            elif ':' in a[0]:hdr[a[0][:-1]]=a[1:]
        cid=hdr['SNID'][0];c='sn'+cid.lower();rr=[]
        for o in obs:
            if o['FLT'] not in MAP.values():continue
            rr.append((o['FLT'],norm(D(o['MJD'])-53000),norm(o['MAG']),norm(o['MAGERR']),norm(o['FLUXCAL']),norm(o['FLUXCALERR'])))
        pred=[]
        for cc,b,t,m,e in byobj.get(c,[]):
            f=10**(-.4*(float(m)-27.5));err=f*(10**(.4*float(e))-1)
            pred.append((MAP[b],t,m,e,norm(f'{f:.5e}'),norm(f'{err:.5e}')))
        pm=Counter((r[0],r[1]) for r in pred);rm=Counter((r[0],r[1]) for r in rr)
        pg=Counter(r[:4] for r in pred);rg=Counter(r[:4] for r in rr)
        pc=Counter(pred);rc=Counter(rr)
        obj=dict(CID=cid,file=fname,archive_found=c in archive_names,archive_has_OPTICAL=c in byobj,expected_OPTICAL_rows=len(pred),released_OPTICAL_rows=len(rr),metadata_equal=pm==rm,mag_error_equal=pg==rg,full_converter_equal=pc==rc,expected_not_released=sum((pc-rc).values()),released_not_expected=sum((rc-pc).values()))
        objects.append(obj)
        # Complete mismatch records, retaining original printed-derived decimal values.
        for label,counts in [('expected_unmatched',pc-rc),('release_unmatched',rc-pc)]:
            for r,n in sorted(counts.items()):rowledger.append(dict(CID=cid,kind=label,multiplicity=n,band=r[0],time=r[1],mag=r[2],magerr=r[3],flux=r[4],fluxerr=r[5]))
        inv={v:k for k,v in MAP.items()}
        for r,n in rc.items():
            physical=rawmeta.get((c,inv[r[0]],r[1]),set())
            release.append(dict(CID=cid,band=r[0],time=r[1],multiplicity=n,raw_physical_filter='|'.join(sorted(physical)),physical_unique=len(physical)==1,converter_identity=pc[r]>=n))
    def save(name,rows,fields=None):
        with (OUT/name).open('w') as f:
            w=csv.DictWriter(f,fieldnames=fields or list(rows[0]));w.writeheader();w.writerows(rows)
    save('objects.csv',objects);save('released-physical-filter-ledger.csv',release)
    save('release-mismatches.csv',rowledger,['CID','kind','multiplicity','band','time','mag','magerr','flux','fluxerr'])
    mismatches=[]
    for label,counts in [('raw_expected_unmatched',expected-actual),('snpy_unmatched',actual-expected)]:
        for r,n in sorted(counts.items()):mismatches.append(dict(kind=label,multiplicity=n,name=r[0],band=r[1],time=r[2],mag=r[3],magerr=r[4]))
    save('archive-mismatches.csv',mismatches,['kind','multiplicity','name','band','time','mag','magerr'])
    changed={}
    for b in ['V0']:
        subset=Counter((c,MERGE[b],t,m,e) for c,bb,t,m,e in raw if bb==b)
        changed[b]=dict(raw_rows=sum(subset.values()),unmatched_after_relabel=sum((subset-actual).values()),release_rows_with_unique_physical_label=sum(r['multiplicity'] for r in release if r['raw_physical_filter']==b),release_rows_unique_and_converter_exact=sum(r['multiplicity'] for r in release if r['raw_physical_filter']==b and r['converter_identity']))
    result=dict(raw_OPTICAL_rows=len(raw),SNpy_OPTICAL_rows=len(snpy),raw_to_SNpy_complete_decimal_identity=expected==actual,raw_unmatched=sum((expected-actual).values()),SNpy_unmatched=sum((actual-expected).values()),changed_filter_groups=changed,listed_objects=len(objects),archive_present_objects=sum(r['archive_found'] for r in objects),metadata_exact_objects=sum(r['metadata_equal'] for r in objects),mag_error_exact_objects=sum(r['mag_error_equal'] for r in objects),converter_exact_objects=sum(r['full_converter_equal'] for r in objects),nonexact_objects=[r for r in objects if not r['full_converter_equal']],release_OPTICAL_rows=sum(r['multiplicity'] for r in release),ambiguous_or_absent_raw_label_rows=sum(r['multiplicity'] for r in release if not r['physical_unique']),protocol_sha256=sha(OUT/'protocol.json'),scope=p['scope'])
    dump(OUT/'result.json',result);(OUT/'executed-source.py').write_bytes(Path(__file__).read_bytes())
    dump(OUT/'manifest.json',{f.name:sha(f) for f in OUT.iterdir() if f.is_file() and f.name!='manifest.json'})
    print(json.dumps(result,indent=2))

if __name__=='__main__':{'prepare':prepare,'run':run}[sys.argv[1]]()
