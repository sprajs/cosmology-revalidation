"""Source-defined date/label check. Magnitudes/errors never parsed."""
from pathlib import Path
from decimal import Decimal
import tarfile,csv,json,hashlib,collections
P=Path(__file__).parent;ROOT=Path.cwd();AR=ROOT/'runs/research_2026_09_26/csp_dr3_provenance/CSP_Photometry_DR3.tgz'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
pr={'scope':'Raw opticalV0 and per-SN SNooPy labels by metadata date only; never read numeric magnitude/error','archive_sha256':sha(AR),'source_reading':'Krisciunas2017 section6.1.1: LC3014 beforeJan14,2006; LC3009 untilJan25; LC9844after. ArchiveREADME agrees. Civildate UTC boundaries MJD53749 and53760; literal integerJD boundary minus2400000.5 sensitivity53748.5 and53759.5.','date_offset':53000,'source_script_sha256':sha(Path(__file__))}
assert not(P/'date-protocol.json').exists();(P/'date-protocol.json').write_text(json.dumps(pr,indent=2)+'\n')
with tarfile.open(AR) as t:
 raw=[];actual=collections.Counter()
 for line in t.extractfile('DR3/SN_photo.dat').read().decode().splitlines():
  a=line.split()
  if a[1]=='V0':raw.append((a[0].lower(),Decimal(a[2])))
 for m in t:
  if not(m.isfile() and m.name.endswith('_snpy.txt')):continue
  lines=t.extractfile(m).read().decode().splitlines();name=lines[0].split()[0].lower();band=None
  for line in lines[1:]:
   a=line.split()
   if not a or a[0].startswith('#'):continue
   if a[0]=='filter':band=a[1];continue
   if band in ['V','V0','V1']:actual[name,Decimal(a[0]),band]+=1
 rows=[]
 for (name,time),mult in sorted(collections.Counter(raw).items()):
  candidates=[b for b in ['V','V0','V1'] if actual[name,time,b]>0];assert len(candidates)==1 and actual[name,time,candidates[0]]==mult
  mjd=time+53000;expected='V0' if mjd<53749 else 'V1' if mjd<53760 else 'V';expectedJD='V0' if mjd<Decimal('53748.5') else 'V1' if mjd<Decimal('53759.5') else 'V'
  rows.append({'CID':name,'time':str(time),'MJD':str(mjd),'multiplicity':mult,'SNpy_label':candidates[0],'civil_date_expected':expected,'literal_JD_expected':expectedJD,'civil_label_agrees':candidates[0]==expected,'literal_JD_label_agrees':candidates[0]==expectedJD})
with (P/'date-label-ledger.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
sumrows=lambda rr:sum(r['multiplicity'] for r in rr)
result={'rawV0_rows':len(raw),'rawV0_metadata_groups':len(rows),'SNpy_label_counts':{b:sumrows([r for r in rows if r['SNpy_label']==b]) for b in ['V0','V1','V']},'civil_date_label_matches':sumrows([r for r in rows if r['civil_label_agrees']]),'literal_JD_label_matches':sumrows([r for r in rows if r['literal_JD_label_agrees']]),'date_ranges_by_label':{b:[str(min(Decimal(r['MJD']) for r in rows if r['SNpy_label']==b)),str(max(Decimal(r['MJD']) for r in rows if r['SNpy_label']==b))] for b in ['V0','V']},'civil_date_mismatches':[r for r in rows if not r['civil_label_agrees']],'literal_JD_mismatches':[r for r in rows if not r['literal_JD_label_agrees']],'magnitudes_or_errors_read':False,'protocol_sha256':sha(P/'date-protocol.json')}
(P/'date-result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
