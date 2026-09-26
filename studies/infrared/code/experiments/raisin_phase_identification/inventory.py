"""Metadata-only cadence inventory: never converts FLUXCAL, MAG or their errors."""
from pathlib import Path
import json,csv,hashlib
ROOT=Path(__file__).resolve().parents[4];O=Path(__file__).resolve().parent
REL=ROOT/'sources/repos/djones1040__RAISIN_DataRelease@a383c4b';DATA=REL/'photometry/RAISIN/CSPDR3_RAISIN'
sha=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()
def table(path):
 lines=path.read_text().splitlines();h={};rows=[];cols=None
 for line in lines:
  x=line.split()
  if not x:continue
  if x[0]=='VARLIST:':cols=x[1:]
  elif x[0]=='OBS:':
   # Explicit allowlist prevents observed brightness or uncertainties entering design.
   row=dict(zip(cols,x[1:]));rows.append(dict(MJD=float(row['MJD']),band=row['FLT'],FIELD=row.get('FIELD',''),has_PHOTFLAG='PHOTFLAG' in row))
  elif ':' in x[0]:h[x[0][:-1]]=x[1:]
 return h,rows

def main():
 protocol=json.loads((O/'protocol.json').read_text());assert sha(__file__)==protocol['inventory_source_sha256']
 files=[DATA/x for x in (DATA/'CSPDR3_RAISIN.LIST').read_text().splitlines() if x.strip() and not x.startswith('#')]
 overlap={x['CID'].lower() for x in csv.DictReader((ROOT/'runs/research_2026_09_26/bayesn_distance_identification/M20-training-name-overlap.csv').open())}
 objects=[];bands=[];epochs=[];ids=[]
 for path in files:
  h,rows=table(path);cid=h['SNID'][0];assert cid not in ids;ids.append(cid);z=float(h['REDSHIFT_HELIO'][0]);t0=float(h['PEAKMJD'][0])
  for j,r in enumerate(rows):r.update(CID=cid,source_row=j,phase=(r['MJD']-t0)/(1+z))
  for b in ['Y','y','J','j','H','h']:
   bb=[r for r in rows if r['band']==b];rec=dict(CID=cid,band=b,total_epochs=len(bb),M20_name_overlap=cid.lower() in overlap)
   for label,lo,hi in [('early',-7,7),('late',10,20),('late_extended',15,30)]:
    selected=[r for r in bb if lo<=r['phase']<=hi];nights=sorted(set(int(r['MJD']) for r in selected));ph=[r['phase'] for r in selected]
    rec[label+'_epochs']=len(selected);rec[label+'_nights']=len(nights);rec[label+'_phase_span']=max(ph)-min(ph) if ph else 0.
   rec['cadence_eligible']=all(rec[l+'_nights']>=3 and rec[l+'_phase_span']>=2 for l in ['early','late']);bands.append(rec)
  ob=[r for r in rows if r['band'] in ['B','g','n','m','o']];earlyopt=[r for r in ob if -10<=r['phase']<=10]
  shapeopt=[r for r in ob if 10<r['phase']<=40]
  optical= len(set(int(r['MJD']) for r in earlyopt))>=3 and any(r['phase']<0 for r in earlyopt) and any(r['phase']>0 for r in earlyopt) and len(set(int(r['MJD']) for r in shapeopt))>=3
  elig=[r['band'] for r in bands if r['CID']==cid and r['cadence_eligible']]
  objects.append(dict(CID=cid,file=str(path.relative_to(ROOT)),zHEL=z,header_peak=t0,M20_name_overlap=cid.lower() in overlap,
   primary_Y='Y' in elig,eligible_exact_filters=','.join(elig),has_second_NIR_band=any(b in elig for b in ['J','j','H','h']),optical_cadence_proxy=optical,
   optical_proxy_NOT_precision_or_independence_gate=True,epochs=len(rows),PHOTFLAG_available=all(r['has_PHOTFLAG'] for r in rows)))
  epochs.extend(rows)
 def savecsv(name,rr):
  with(O/name).open('w') as f:w=csv.DictWriter(f,rr[0].keys());w.writeheader();w.writerows(rr)
 savecsv('objects.csv',objects);savecsv('exact-filter-cadence.csv',bands);savecsv('metadata-epochs.csv',epochs)
 eligible=[r for r in objects if r['primary_Y']];paired=[r for r in eligible if r['has_second_NIR_band']]
 summary=dict(listed_objects=len(objects),physical_photometry_files=len(list(DATA.glob('*.DAT'))),listed_epochs=len(epochs),primary_Y_count=len(eligible),primary_Y_ids=[r['CID'] for r in eligible],
   primary_Y_with_second_NIR=len(paired),primary_Y_with_optical_proxy=sum(r['optical_cadence_proxy'] for r in eligible),
   M20_name_overlap_total=sum(r['M20_name_overlap'] for r in objects),M20_name_overlap_Y=sum(r['M20_name_overlap'] for r in eligible),not_in_M20_name_roster_Y=sum(not r['M20_name_overlap'] for r in eligible),
   cadence_gate_pass=len(eligible)>=15 and len(paired)>=10,
   flux_values_or_errors_scored=False,selection_and_independent_training_gates_NOT_closed=True,
   exact_filter_eligible_counts={b:sum(r['band']==b and r['cadence_eligible'] for r in bands) for b in ['Y','y','J','j','H','h']})
 (O/'inventory-result.json').write_text(json.dumps(summary,indent=2)+'\n')
 sources=files+[DATA/'CSPDR3_RAISIN.LIST',DATA/'CSPDR3_RAISIN.README',REL/'kcor/kcor_CSPDR3_BD17.input',REL/'kcor/kcor_CSPDR3_BD17.fits',ROOT/'runs/research_2026_09_26/bayesn_distance_identification/M20-training-name-overlap.csv',Path(__file__),O/'protocol.json']
 (O/'input-manifest.json').write_text(json.dumps({str(f.relative_to(ROOT)):sha(f) for f in sources},indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
