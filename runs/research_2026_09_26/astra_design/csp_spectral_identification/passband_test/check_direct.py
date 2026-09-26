"""Separate direct trapezoid/reaggregation check; does not import runner."""
from pathlib import Path
import csv,json,hashlib
import numpy as np
from astropy.io import fits
P=Path(__file__).parent;ROOT=Path.cwd();S=P.parent
p=json.loads((P/'protocol.json').read_text());sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
for n,h in p['inputs_sha256'].items():assert sha(ROOT/n)==h
fd=ROOT/'sources/repos/PantheonPlusSH0ES__DataRelease@7fc6805/Pantheon+_Data/2_CALIBRATION/filters/CSP_TAMU_20180316'
k=ROOT/'sources/repos/djones1040__RAISIN_DataRelease@a383c4b/kcor/kcor_CSPDR3_BD17.fits'
with fits.open(k) as h:rw=np.asarray(h['PrimarySED'].data.field(0),float);rf=np.asarray(h['PrimarySED'].data['BD17'],float)
def integral(w,f,band):
 x=np.linspace(band[0,0],band[-1,0],int(np.ceil((band[-1,0]-band[0,0])/.25))+1)
 assert x[0]>=w[0] and x[-1]<=w[-1]
 v=x*np.interp(x,band[:,0],band[:,1])*np.interp(x,w,f)
 assert np.isfinite(v).all()
 return np.trapezoid(v,x)
bands={n:np.loadtxt(fd/n) for v in p['pairs'].values() for n in v};refs={n:integral(rw,rf,b) for n,b in bands.items()}
rows=list(csv.DictReader((P/'spectrum-results.csv').open()));out=[];cache={};maxgap=0
for r in rows:
 if r['target']!='observed':continue
 n=r['filename']
 if n not in cache:cache[n]=np.loadtxt(S/'spectra/observed_spectra/spec_txt'/f'{n}.txt')
 d=cache[n];a,b=p['pairs'][r['pair']];ia=integral(d[:,0]*1e4,d[:,1],bands[a]);ib=integral(d[:,0]*1e4,d[:,1],bands[b]);v=-2.5*np.log10(ib*refs[a]/(ia*refs[b]));gap=abs(v-float(r['delta_mag_equal_BD17']));maxgap=max(maxgap,gap)
 out.append({'name':r['name'],'pair':r['pair'],'phase':float(r['phase']),'v':v})
assert maxgap<1e-4
pairs=list(csv.DictReader((P/'paired-object-results.csv').open()));pairgaps=[]
for r in pairs:
 if r['target']!='observed' or r['clock_sigma_shift']!='0':continue
 q=[v for v in out if v['name']==r['name'] and v['pair']==r['pair']];e=[v['v'] for v in q if -7<=v['phase']<=7];l=[v['v'] for v in q if 10<=v['phase']<=20]
 assert len(e)==int(r['early_spectra']) and len(l)==int(r['late_spectra'])
 pairgaps.append(abs((np.mean(l)-np.mean(e))-float(r['delta_mag_equal_BD17'])))
# Exact aggregation of saved rows is independent from direct numerical quadrature.
sgaps=[]
for r in csv.DictReader((P/'summary.csv').open()):
 q=[v for v in pairs if v['pair']==r['pair'] and v['target']==r['target'] and v['clock_sigma_shift']==r['clock_sigma_shift'] and (r['exclude_2012fr']=='False' or v['name']!='SN2012fr')]
 assert len(q)==int(r['N_objects']);sgaps.append(abs(np.mean([float(v['delta_mag_equal_BD17']) for v in q])-float(r['equal_object_mean'])))
res={'pass':True,'source_sha256':sha(Path(__file__)),'verified_input_hashes':len(p['inputs_sha256']),'direct_observed_spectrum_pairs':len(out),'max_direct_quarterA_delta_mag_difference':maxgap,'paired_objects_compared':len(pairgaps),'max_direct_paired_phase_difference':max(pairgaps),'max_summary_reaggregation_difference':max(sgaps),'scope':'Independent direct interpolation/integration and aggregation, not independent photometric calibration or covariance.'}
(P/'independent-check.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
