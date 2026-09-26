"""Independent Simpson-polynomial primary spectral-integral/reaggregation check."""
from pathlib import Path
import csv,json,hashlib
from collections import defaultdict
import numpy as np
from astropy.io import fits

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'runs/research_2026_09_26/astra_design/csp_spectral_identification'
P=BASE/'passband_test'
OUT=ROOT/'runs/research_2026_09_26/csp_spectral_root_review'
FD=ROOT/'sources/repos/PantheonPlusSH0ES__DataRelease@7fc6805/Pantheon+_Data/2_CALIBRATION/filters/CSP_TAMU_20180316'
KC=ROOT/'sources/repos/djones1040__RAISIN_DataRelease@a383c4b/kcor/kcor_CSPDR3_BD17.fits'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()

def integral(wave,flux,band):
    lo,hi=band[0,0],band[-1,0]
    assert wave[0]<=lo and wave[-1]>=hi
    x=np.unique(np.concatenate([band[:,0],wave[(wave>lo)&(wave<hi)]]))
    mid=(x[:-1]+x[1:])/2
    # Simpson is exact for lambda * piecewise-linear T * piecewise-linear f.
    end=x*np.interp(x,band[:,0],band[:,1])*np.interp(x,wave,flux)
    val=mid*np.interp(mid,band[:,0],band[:,1])*np.interp(mid,wave,flux)
    assert np.all(np.isfinite(end)) and np.all(np.isfinite(val))
    answer=np.sum((x[1:]-x[:-1])*(end[:-1]+4*val+end[1:])/6)
    assert answer>0
    return answer

def main():
    OUT.mkdir(exist_ok=False)
    protocol=json.loads((P/'protocol.json').read_text())
    for path,digest in protocol['inputs_sha256'].items():assert sha(ROOT/path)==digest,path
    manifest=json.loads((P/'manifest.json').read_text())
    for path,digest in manifest['files_sha256'].items():assert sha(ROOT/path)==digest,path
    (OUT/'review-design.json').write_text(json.dumps({'status':'Independent verification after outcome report, not a new preregistered result. No executor imports.','scope':'Every3018 primary and shape-transport integral, independently integrate piecewise polynomials with endpoint/midpoint Simpson. Reaggregate503actual-z contrasts through frozen per-object phasewindows and primary summaries. Does not establish spectral physical calibration or full noise covariance.','protocol_sha256':sha(P/'protocol.json'),'source_sha256':sha(Path(__file__))},indent=2)+'\n')
    with fits.open(KC) as h:
        rw=np.array(h['PrimarySED'].data.field(0),float);rf=np.array(h['PrimarySED'].data['BD17'],float)
    bands={n:np.loadtxt(FD/n) for pp in protocol['pairs'].values() for n in pp}
    refs={n:integral(rw,rf,b) for n,b in bands.items()}
    cohort=list(csv.DictReader((P/'frozen-cohort.csv').open()))
    expected={(r['filename'],r['pair'],target):r for r in cohort for target in ['observed','0','0.01','0.03','0.05','0.08']}
    executed=list(csv.DictReader((P/'spectrum-results.csv').open()))
    assert len(expected)==len(executed)==3018
    cache={};gap=0.;independent=defaultdict(list);checked=[]
    for r in executed:
        key=(r['filename'],r['pair'],r['target']);meta=expected.pop(key)
        assert all(r[k]==meta[k] for k in ['name','phase','zhel','eTmax','EBV_MW'])
        if r['filename'] not in cache:cache[r['filename']]=np.loadtxt(BASE/'spectra/observed_spectra/spec_txt'/f"{r['filename']}.txt")
        d=cache[r['filename']];fac=(1+float(r['ztarget']))/(1+float(r['zhel']))
        w=d[:,0]*1e4*fac;f=d[:,1]/fac;a,b=protocol['pairs'][r['pair']]
        value=-2.5*np.log10(integral(w,f,bands[b])*refs[a]/(integral(w,f,bands[a])*refs[b]))
        error=abs(value-float(r['delta_mag_equal_BD17']));gap=max(gap,error);assert error<1e-10
        if r['target']=='observed':independent[(r['pair'],r['name'])].append((float(r['phase']),value))
    assert not expected
    saved={(r['pair'],r['name']):r for r in csv.DictReader((P/'paired-object-results.csv').open()) if r['target']=='observed' and r['clock_sigma_shift']=='0'}
    paired=defaultdict(list);pgap=0.
    for (pair,name),rr in independent.items():
        early=[v for p,v in rr if -7<=p<=7];late=[v for p,v in rr if 10<=p<=20]
        if not early or not late:continue
        v=float(np.mean(late)-np.mean(early));r=saved.pop((pair,name))
        assert len(early)==int(r['early_spectra']) and len(late)==int(r['late_spectra'])
        pgap=max(pgap,abs(v-float(r['delta_mag_equal_BD17'])));paired[pair].append(v)
    assert not saved and pgap<1e-10
    results=json.loads((P/'result.json').read_text());summary={}
    for r in results['primary_summaries']:
        values=np.array(paired[r['pair']]);assert len(values)==r['N_objects']
        ours={'equal_object_mean':float(values.mean()),'sample_SD_descriptive':float(values.std(ddof=1)),'median':float(np.median(values)),'min':float(values.min()),'max':float(values.max())}
        assert max(abs(v-r[k]) for k,v in ours.items())<1e-10
        summary[r['pair']]={'N':len(values),**ours}
    result={'pass':True,'verified_input_hashes':len(protocol['inputs_sha256']),'verified_manifest_hashes':len(manifest['files_sha256']),'independent_integrals':3018,'max_integral_mag_gap':gap,'max_paired_phase_gap':pgap,'primary_summary':summary,'scope':'Independent numerical verification only; unavailable full spectral covariance and independent phase/population calibration remain open.'}
    (OUT/'result.json').write_text(json.dumps(result,indent=2)+'\n');(OUT/'executed-source.py').write_bytes(Path(__file__).read_bytes())
    (OUT/'manifest.json').write_text(json.dumps({p.name:sha(p) for p in OUT.iterdir() if p.is_file() and p.name!='manifest.json'},indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
