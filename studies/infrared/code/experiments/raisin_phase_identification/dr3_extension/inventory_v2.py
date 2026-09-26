"""Primary archive metadata linkage; brightness and errors are never parsed."""
from pathlib import Path
import json,csv,hashlib,tarfile,sys,math
O=Path(__file__).resolve().parent;BASE=O.parent;ROOT=BASE.parents[3]
sys.path.insert(0,str(BASE));from inventory import table
# Module-name collision avoided by executing this file as main; imported inventory is parent path first.
sha=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()
def main():
 p=json.loads((O/'protocol-v2.json').read_text());archive=ROOT/p['archive'];assert sha(archive)==p['archive_sha256'] and sha(__file__)==p['source_sha256']
 data=ROOT/'sources/repos/djones1040__RAISIN_DataRelease@a383c4b/photometry/RAISIN/CSPDR3_RAISIN'
 refs={}
 for name in (data/'CSPDR3_RAISIN.LIST').read_text().splitlines():
  h,r=table(data/name);refs[h['SNID'][0].lower()]=(h,r)
 objects=[];matches=[];cadence=[];rawfilters={}
 with tarfile.open(archive) as t:
  for member in t:
   if not(member.isfile() and member.name.endswith('_snpy.txt')):continue
   lines=t.extractfile(member).read().decode().splitlines();head=lines[0].split();cid=head[0][2:];z,ra,dec=map(float,head[1:4]);rows=[];band=None
   for linenum,line in enumerate(lines[1:],2):
    a=line.split()
    if not a or a[0].startswith('#'):continue
    if a[0]=='filter':band=a[1];continue
    # Ignore magnitude and error text completely.
    mjd=float(a[0])+53000.;rows.append(dict(CID=cid,source_line=linenum,raw_filter=band,MJD=mjd));rawfilters[band]=rawfilters.get(band,0)+1
   linked=cid.lower() in refs;rec=dict(CID=cid,archive_member=member.name,epochs=len(rows),matched_listed_name=linked,has_source_peak=False)
   if linked:
    h,rr=refs[cid.lower()];zref=float(h['REDSHIFT_HELIO'][0]);rra=float(h['RA'][0]);rdec=float(h['DEC'][0]);sep=3600*math.hypot((ra-rra)*math.cos(math.radians(dec)),dec-rdec)
    rec.update(coordinate_difference_arcsec=sep,zHEL_difference=z-zref,identity_gate=sep<=5 and abs(z-zref)<=.0002,header_peak=float(h['PEAKMJD'][0]))
    if not rec['identity_gate']:
     objects.append(rec);continue
    for r in rows:
     mapped=p['candidate_filter_map'].get(r['raw_filter']);cand=[x for x in rr if x['band']==mapped and abs(x['MJD']-r['MJD'])<=.00055] if mapped else []
     matches.append(dict(**r,candidate_SNANA_filter=mapped,match_count=len(cand),nearest_abs_delta_MJD=min([abs(x['MJD']-r['MJD']) for x in cand],default=None)))
    for b in ['Y','Ydw','J','Jrc2','Jdw','H','Hdw']:
     rr2=[r for r in rows if r['raw_filter']==b];s=dict(CID=cid,raw_filter=b,epochs=len(rr2),provisional_peak_from_linked_optical_NIR_header=True)
     for label,lo,hi in [('early',-7,7),('late',10,20)]:
      pp=[(r['MJD']-rec['header_peak'])/(1+zref) for r in rr2 if lo<=(r['MJD']-rec['header_peak'])/(1+zref)<=hi];sel=[r for r in rr2 if lo<=(r['MJD']-rec['header_peak'])/(1+zref)<=hi]
      s[label+'_nights']=len(set(int(r['MJD']) for r in sel));s[label+'_span']=max(pp)-min(pp) if pp else 0.
     s['eligible']=all(s[k+'_nights']>=3 and s[k+'_span']>=2 for k in ['early','late']);cadence.append(s)
   objects.append(rec)
 def save(name,rr):
  fields=list(dict.fromkeys(k for r in rr for k in r))
  with(O/name).open('w') as f:w=csv.DictWriter(f,fields);w.writeheader();w.writerows(rr)
 save('objects.csv',objects);save('epoch-time-filter-match.csv',matches);save('exact-filter-cadence.csv',cadence)
 summary=dict(source_objects=len(objects),source_epochs=sum(r['epochs'] for r in objects),raw_filter_counts=rawfilters,
  linked_listed_objects=sum(r['matched_listed_name'] for r in objects),identity_gate_failures=[r for r in objects if r.get('identity_gate') is False],objects_without_linked_peak=[r['CID'] for r in objects if not r['matched_listed_name']],
  mapped_epoch_count=sum(r['candidate_SNANA_filter'] is not None for r in matches),unique_time_filter_matches=sum(r['match_count']==1 for r in matches),ambiguous_matches=sum(r['match_count']>1 for r in matches),unmatched_mapped_epochs=sum(r['candidate_SNANA_filter'] is not None and r['match_count']==0 for r in matches),unmapped_epochs=sum(r['candidate_SNANA_filter'] is None for r in matches),
  eligible_exactY_with_linked_peak=[r['CID'] for r in cadence if r['raw_filter']=='Y' and r['eligible']],
  full134_phase_feasibility='Unevaluated for objects without peak metadata; no brightness-based peak estimated',
  no_observed_flux_error_or_magnitude_parsed=True,source_converter_execution_not_established=True)
 (O/'result.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
