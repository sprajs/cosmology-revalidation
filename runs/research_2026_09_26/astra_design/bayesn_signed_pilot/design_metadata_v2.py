"""Metadata-only planning; never parse observed FLUXCAL values."""
from pathlib import Path
import csv,json,hashlib
import numpy as np
from astropy.io import fits
ROOT=Path.cwd();OUT=ROOT/'runs/research_2026_09_26/astra_design/bayesn_signed_pilot';B=ROOT/'runs/research_2026_09_26/bayesn_distance_identification'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def parse(p):
 h={};rows=[]
 for lineno,l in enumerate(p.read_text().splitlines(),1):
  a=l.split()
  if not a:continue
  if a[0]=='VARLIST:':names=a[1:]
  elif a[0]=='OBS:':
   d=dict(zip(names,a[1:]));rows.append(dict(MJD=float(d['MJD']),band=d.get('FLT',d.get('BAND')),line=lineno))
  elif a[0].endswith(':'):h[a[0][:-1]]=a[1:]
 return h,rows
cp=ROOT/'runs/research_2026_09_26/astra_design/raisin_signed_refit/cohort.csv';co=list(csv.DictReader(cp.open()));paths=[cp,Path(__file__),ROOT/'sources/repos/djones1040__RAISIN_DataRelease@a383c4b/kcor/kcor_DES_NIR.fits']+[ROOT/c[k] for c in co for k in ('raw_path','released_path')]+sorted((B/'release-filters').glob('RAISIN_DES_?.dat'))
protocol={'status':'Metadata only, frozen before these counts; header dates are feasibility anchors only, not valid held-out-NIR inference priors or final masks. No FLUXCAL values parsed.', 'cohort':'all10existingDES16lineageobjects','pilot':'lowest/highest zHEL; tie lexicalCID','reference_safe_phase_intervals':[[-5,35],[0,30]],'final_mask_rule':'To be instantiated only after a new optical-only timing localization under a separately frozen executor; no NIR flux used for clock, initialization or masks.','paths':{str(p.relative_to(ROOT)):sha(p) for p in paths}}
assert not (OUT/'metadata-protocol-v2.json').exists();save(OUT/'metadata-protocol-v2.json',protocol)
rows=[];bands=[]
for c in co:
 h,nir=parse(ROOT/c['released_path']);_,opt=parse(ROOT/c['raw_path']);z=float(c['zHEL']);t0=float(c['peak_header']);bp={}
 for band in 'grizJH':
  f=B/'release-filters'/f'RAISIN_DES_{band}.dat'
  if band=='J':
   with fits.open(ROOT/'sources/repos/djones1040__RAISIN_DataRelease@a383c4b/kcor/kcor_DES_NIR.fits') as hdu:a=np.column_stack([hdu['FilterTrans'].data.field(0),hdu['FilterTrans'].data['WFC3_IR_F125W-J']])
  else:a=np.loadtxt(f)
  l=a[:,0]/(1+z);w=a[:,0]*np.maximum(a[:,1],0);norm=np.trapezoid(w,a[:,0]);outside=(l<3000)|(l>18500);active=a[:,1]>.01*a[:,1].max()
  lo,hi=l[active][[0,-1]];strict=(l[a[:,1]>0].min()>=3000 and l[a[:,1]>0].max()<=18500)
  # integration with explicit endpoint interpolation at model-domain bounds
  g=np.unique(np.r_[a[:,0],3000*(1+z),18500*(1+z)]);v=np.interp(g,a[:,0],w,left=0,right=0);inside=(g>=3000*(1+z))&(g<=18500*(1+z));frac=1-np.trapezoid(v[inside],g[inside])/np.trapezoid(v,g)
  bp[band]=(lo>=3000 and hi<=18500);bands.append(dict(CID=c['CID'],band=band,zHEL=z,rest_1pct_min=float(lo),rest_1pct_max=float(hi),old_forward_gate=bool(bp[band]),positive_photon_weight_outside=max(0,float(frac)),all_positive_support_inside=bool(strict)))
 op=[r for r in opt if r['band'] in 'griz'];ni=[r for r in nir if r['band'] in 'JH'];r=dict(CID=c['CID'],zHEL=z,header_peak_feasibility_only=t0,N_optical=len(op),N_NIR=len(ni),NIR_header_phase=[(x['MJD']-t0)/(1+z) for x in ni],eligible_optical_bands=''.join(b for b in 'griz' if bp[b]))
 for low,high in [(-5,35),(0,30)]:
  sel=[x for x in op if bp[x['band']] and low<=(x['MJD']-t0)/(1+z)<=high];r[f'optical_{low}_{high}']=len(sel);r[f'optical_nights_{low}_{high}']=len(set(int(x['MJD']) for x in sel));r[f'NIR_{low}_{high}']=sum(low<=(x['MJD']-t0)/(1+z)<=high for x in ni)
 rows.append(r)
for name,data in [('cohort-metadata.csv',rows),('band-support.csv',bands)]:
 with (OUT/name).open('w') as f:w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
sortedco=sorted(rows,key=lambda r:(r['zHEL'],r['CID']));save(OUT/'metadata-result.json',dict(cohort=rows,pilot=[sortedco[0]['CID'],sortedco[-1]['CID']],N_optical=sum(r['N_optical'] for r in rows),N_NIR=sum(r['N_NIR'] for r in rows),no_flux_values_parsed=True))
print(json.dumps({'pilot':[sortedco[0]['CID'],sortedco[-1]['CID']],'cohort':rows},indent=2))
