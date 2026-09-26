"""Frozen descriptive passband experiment; no distances or population significance."""
from pathlib import Path
import csv,json,hashlib,sys
import numpy as np
from astropy.io import fits
from numpy.polynomial.legendre import leggauss
import extinction
P=Path(__file__).parent;S=P.parent;ROOT=Path.cwd()
FD=ROOT/'sources/repos/PantheonPlusSH0ES__DataRelease@7fc6805/Pantheon+_Data/2_CALIBRATION/filters/CSP_TAMU_20180316'
KC=ROOT/'sources/repos/djones1040__RAISIN_DataRelease@a383c4b/kcor/kcor_CSPDR3_BD17.fits'
PAIRS={'J_WIRC_minus_RC1':('Jrc1_SWO_TAM_scan_atm.dat','J_DUP_TAM_scan_atm.dat'),'J_WIRC_minus_RC2':('Jrc2_SWO_TAM_scan_atm.dat','J_DUP_TAM_scan_atm.dat'),'H_WIRC_minus_RetroCam':('H_SWO_TAM_scan_atm.dat','H_DUP_TAM_scan_atm.dat'),'Y_WIRC_minus_RetroCam':('Y_SWO_TAM_scan_atm.dat','Y_DUP_TAM_scan_atm.dat')}
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def dump(n,x): (P/n).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def csvout(n,rows):
 with (P/n).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def prepare():
 assert not (P/'protocol.json').exists()
 meta=list(csv.DictReader((S/'spectrum-schema-ledger.csv').open()));cohort=[]
 for r in meta:
  if r['host_contamination']!='No' or not .75<=float(r['sBV'])<=1.18 or float(r['e(Tmax)'])>.5 or not -7<=float(r['epoch'])<=30:continue
  for pair in PAIRS:
   col='SNRH_first' if pair[0]=='H' else 'SNR'+pair[0]
   if float(r[col])>=10:cohort.append({'name':r['name'],'filename':r['filename'],'pair':pair,'phase':float(r['epoch']),'zhel':float(r['zhel']),'eTmax':float(r['e(Tmax)']),'EBV_MW':float(r['EBV_MW'])})
 csvout('frozen-cohort.csv',cohort)
 paths=[Path(__file__),S/'metadata_protocol.json',S/'metadata-result.json',S/'spectra-acquisition.json',S/'spectrum-schema-ledger.csv',P/'frozen-cohort.csv',KC,FD/'README']+[FD/n for n in sorted({n for v in PAIRS.values() for n in v})]
 paths += [S/'spectra/observed_spectra/spec_txt'/f'{n}.txt' for n in sorted({r['filename'] for r in cohort})]
 dump('protocol.json',{'status':'Frozen before spectral passband integrals. Follow-up motivated by metadata filter-label merging and earlier Hsiao synthetic response.','primary_pairs':['J_WIRC_minus_RC1','J_WIRC_minus_RC2'],'secondary':'H_WIRC_minus_RetroCam','control':'Y_WIRC_minus_RetroCam','pairs':PAIRS,'eligible':'Host flagNo,sBV .75–1.18,eTmax≤.5d, per-band tableSNR≥10, restphase[-7,30]. All eligible individual spectra; paired phaseearly[-7,7], late[10,20]. Equal spectra per window then equal objects. No outcome-based exclusion/weights.','counts':{'spectrum_pair_rows':len(cohort),'spectra':len({r['filename'] for r in cohort}),'objects':len({r['name'] for r in cohort})},'units':'TXT observer wavelengthum→A; fluxf_lambda (FITSYUNITS erg/s/cm2/A). Main integrates delivered observer flux directly at actual z; no assumed dust correction. For targetz, wavelength multiplied (1+ztarget)/(1+zsource); overall fluxfactor irrelevant to normalized contrast and retained explicitly.','target_z_secondary':[0,.01,.03,.05,.08],'reference':'BD17 exactKCOR, equal assigned reference magnitude. Absolute natural-system calibration offsets not included. Same-object phase differences cancel static offsets.','gates':'No spectral extrapolation. No clipping negative flux/transmission. Require all nonzero weight spectral pixels finite and both countspositive. Missing errors do not block descriptive mean but diagonal error diagnostic is null if any nonzero derivative pixel lacks validerror. No covariance fills or discovery significance.','quadrature':'Exact Gauss2 on union of linear spectrum/filter knots, Gauss4 agreement≤1e-10mag, independent .5A trapezoid≤1e-4mag; if fails, predeclared .1A check before certification. Gray×7 test≤1e-12mag.','telluric':'Delivered archive lacks per-region visual quality mask. Full integral is descriptive. Fixed observed-frame intervals[13500,14500] and[18000,19500]A excluded in BOTH compared filters and BD17 reference for a separate common-support diagnostic. This changes estimand; not a filled full-band estimate. Record absolute lambda*T weight lost and contrast shift; >.005mag labels sensitivity, never removes primary row.','sensitivities':'At actualz only: conditional MWderedden factor10^(.4 F99(RV3.1,EBV_table)); no claim archive was corrected/un-corrected. Smooth multiplicative tilt corresponding±.2mag across8000–24000A. Same-object peak±quotederror reassignment for phase-window sensitivity, explicitly model-derived not independent timing. RemoveSN2012fr summary fixed in advance.','interpretation':'Measures passband-dependent shape of observed calibrated spectra. Gray normalization/distance cancels. Dust/color, telluric/slit response, phase labels and population/selection transport remain. No inference of gray luminosity, actual photometric calibration correction, corrected likelihood or cosmology.','inputs_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths}})
 print(json.dumps({'protocol_sha256':sha(P/'protocol.json'),'counts':len(cohort)}))

def weights(wave,band,order=2,masks=()):
 lo,hi=band[0,0],band[-1,0]
 if wave[0]>lo or wave[-1]<hi:return None
 if order==0:
  x=np.linspace(lo,hi,int(np.ceil((hi-lo)/.5))+1);q=np.full(len(x),(hi-lo)/(len(x)-1));q[[0,-1]]*=.5
 elif order==-1:
  x=np.linspace(lo,hi,int(np.ceil((hi-lo)/.1))+1);q=np.full(len(x),(hi-lo)/(len(x)-1));q[[0,-1]]*=.5
 else:
  bounds=[v for ab in masks for v in ab if lo<v<hi]
  k=np.unique(np.r_[band[:,0],wave[(wave>lo)&(wave<hi)],bounds]);mid=(k[1:]+k[:-1])/2;half=(k[1:]-k[:-1])/2;n,v=leggauss(order)
  x=(mid[:,None]+half[:,None]*n).ravel();q=(half[:,None]*np.broadcast_to(v,(len(half),len(v)))).ravel()
 t=np.interp(x,band[:,0],band[:,1]);u=q*x*t
 for a,b in masks:u[(x>=a)&(x<=b)]=0
 ix=np.minimum(np.searchsorted(wave,x,side='right')-1,len(wave)-2);a=(x-wave[ix])/(wave[ix+1]-wave[ix]);w=np.zeros(len(wave));np.add.at(w,ix,u*(1-a));np.add.at(w,ix+1,u*a)
 return w

def count(w,f):
 if w is None:return None
 use=w!=0
 if not np.isfinite(f[use]).all():return None
 v=float(w[use]@f[use]);return v if v>0 and np.isfinite(v) else None

def dm(ca,cb,ra,rb):
 if any(x is None for x in [ca,cb,ra,rb]):return None
 return float(-2.5*np.log10((cb/rb)/(ca/ra)))
def run():
 p=json.loads((P/'protocol.json').read_text())
 for n,h in p['inputs_sha256'].items():assert sha(ROOT/n)==h,n
 assert not (P/'spectrum-results.csv').exists()
 cohort=list(csv.DictReader((P/'frozen-cohort.csv').open()));bands={n:np.loadtxt(FD/n) for v in PAIRS.values() for n in v}
 with fits.open(KC) as h:rw=np.array(h['PrimarySED'].data.field(0),float);rf=np.array(h['PrimarySED'].data['BD17'],float)
 refs={o:{n:count(weights(rw,b,o),rf) for n,b in bands.items()} for o in [2,4,0,-1]}
 cache={};result=[]
 for ri,r in enumerate(cohort):
  if r['filename'] not in cache:cache[r['filename']]=np.loadtxt(S/'spectra/observed_spectra/spec_txt'/f"{r['filename']}.txt")
  d=cache[r['filename']];w0=d[:,0]*1e4;f0=d[:,1];e0=d[:,2];z=float(r['zhel']);a,b=PAIRS[r['pair']]
  for target in ['observed']+p['target_z_secondary']:
   zt=z if target=='observed' else target;fac=(1+zt)/(1+z);wave=w0*fac;flux=f0/fac;err=e0/fac
   wa=weights(wave,bands[a]);wb=weights(wave,bands[b]);ca=count(wa,flux);cb=count(wb,flux);v=dm(ca,cb,refs[2][a],refs[2][b]);rec={**r,'target':target,'ztarget':zt,'delta_mag_equal_BD17':v,'status':'PASS' if v is not None else 'unsupported_or_nonfinite_flux_or_nonpositive_counts'}
   rec.update({k:None for k in ['gauss2_4_gap','trap_halfA_gap','trap_tenthA_gap','gray_gap','diagonal_only_error_not_total','missing_error_derivative_pixels','telluric_cut_delta_mag','telluric_cut_change','telluric_absweight_a','telluric_absweight_b','foreground_removal_change','tilt_plus_change','tilt_minus_change']})
   if v is not None:
    alt={o:dm(count(weights(wave,bands[a],o),flux),count(weights(wave,bands[b],o),flux),refs[o][a],refs[o][b]) for o in [4,0]}
    rec['gauss2_4_gap']=abs(v-alt[4]);rec['trap_halfA_gap']=abs(v-alt[0]);rec['gray_gap']=abs(v-dm(count(wa,7*flux),count(wb,7*flux),refs[2][a],refs[2][b]))
    assert rec['gauss2_4_gap']<1e-10 and rec['gray_gap']<1e-12
    if rec['trap_halfA_gap']>1e-4:
     vv=dm(count(weights(wave,bands[a],-1),flux),count(weights(wave,bands[b],-1),flux),refs[-1][a],refs[-1][b]);rec['trap_tenthA_gap']=abs(v-vv);assert abs(v-vv)<1e-4
    grad=-2.5/np.log(10)*(wb/cb-wa/ca);use=grad!=0;bad=use&(~np.isfinite(err)|(err<=0));rec['missing_error_derivative_pixels']=int(bad.sum())
    if not bad.any():rec['diagonal_only_error_not_total']=float(np.sqrt(np.sum((grad[use]*err[use])**2)))
    masks=np.array([[13500,14500],[18000,19500]])*fac
    ma=weights(wave,bands[a],masks=masks);mb=weights(wave,bands[b],masks=masks)
    va=dm(count(ma,flux),count(mb,flux),count(weights(rw,bands[a],masks=masks),rf),count(weights(rw,bands[b],masks=masks),rf));rec['telluric_cut_delta_mag']=va;rec['telluric_cut_change']=None if va is None else va-v
    for label,name,full,cut in [('a',a,wa,ma),('b',b,wb,mb)]:
     ab=bands[name].copy();ab[:,1]=abs(ab[:,1]);q=weights(wave,ab);qc=weights(wave,ab,masks=masks);rec['telluric_absweight_'+label]=float(1-qc.sum()/q.sum())
    if target=='observed':
     ef=extinction.fitzpatrick99(np.ascontiguousarray(w0),3.1*float(r['EBV_MW']),3.1)
     corr=10**(.4*ef);rec['foreground_removal_change']=dm(count(wa,flux*corr),count(wb,flux*corr),refs[2][a],refs[2][b])-v
     for sign,lab in [(1,'plus'),(-1,'minus')]:
      mag=sign*.2*np.log(w0/16000)/np.log(3);ff=flux*10**(-.4*mag);rec['tilt_'+lab+'_change']=dm(count(wa,ff),count(wb,ff),refs[2][a],refs[2][b])-v
   result.append(rec)
  if ri%50==0:print('completed spectrum-pair',ri,flush=True)
 csvout('spectrum-results.csv',result)
 paired=[]
 for pair in PAIRS:
  for target in ['observed']+p['target_z_secondary']:
   for name in sorted({r['name'] for r in result}):
    rr=[r for r in result if r['name']==name and r['pair']==pair and r['target']==target and r['status']=='PASS']
    for clock in [0,-1,1]:
     early=[r for r in rr if -7<=float(r['phase'])-clock*float(r['eTmax'])/(1+float(r['zhel']))<=7];late=[r for r in rr if 10<=float(r['phase'])-clock*float(r['eTmax'])/(1+float(r['zhel']))<=20]
     if not early or not late:continue
     row={'pair':pair,'target':target,'name':name,'clock_sigma_shift':clock,'early_spectra':len(early),'late_spectra':len(late)}
     for k in ['delta_mag_equal_BD17','telluric_cut_delta_mag','foreground_removal_change','tilt_plus_change','tilt_minus_change']:
      vals=[r[k] for r in early+late];row[k]=None if any(v is None for v in vals) else float(np.mean([r[k] for r in late])-np.mean([r[k] for r in early]))
     paired.append(row)
 csvout('paired-object-results.csv',paired)
 summaries=[]
 for pair in PAIRS:
  for target in ['observed']+p['target_z_secondary']:
   for clock in [0,-1,1]:
    for exclude in [False,True]:
     rr=[r for r in paired if r['pair']==pair and r['target']==target and r['clock_sigma_shift']==clock and(not exclude or r['name']!='SN2012fr')]
     if not rr:continue
     x=np.array([r['delta_mag_equal_BD17'] for r in rr]);summaries.append({'pair':pair,'target':target,'clock_sigma_shift':clock,'exclude_2012fr':exclude,'N_objects':len(rr),'equal_object_mean':float(x.mean()),'median':float(np.median(x)),'sample_SD_descriptive':float(x.std(ddof=1)) if len(x)>1 else None,'min':float(x.min()),'max':float(x.max()),'max_telluric_change_in_phase_contrast':max([abs(r['telluric_cut_delta_mag']-r['delta_mag_equal_BD17']) for r in rr if r['telluric_cut_delta_mag'] is not None],default=None),'max_foreground_change_in_phase_contrast':max([abs(r['foreground_removal_change']) for r in rr if r['foreground_removal_change'] is not None],default=None),'max_tilt_change_in_phase_contrast':max([abs(r[k]) for r in rr for k in ['tilt_plus_change','tilt_minus_change'] if r[k] is not None],default=None)})
 csvout('summary.csv',summaries)
 good=[r for r in result if r['status']=='PASS'];obs=[r for r in good if r['target']=='observed']
 answer={'protocol_sha256':sha(P/'protocol.json'),'rows':len(result),'passed':len(good),'failed':len(result)-len(good),'quadrature_gauss_max':max(r['gauss2_4_gap'] for r in good),'quadrature_halfA_max':max(r['trap_halfA_gap'] for r in good),'tenthA_checks':sum(r['trap_tenthA_gap'] is not None for r in good),'gray_max':max(r['gray_gap'] for r in good),'observed_pairs':len(obs),'observed_missing_error':sum(r['missing_error_derivative_pixels']>0 for r in obs),'observed_max_telluric_shift':max(abs(r['telluric_cut_change']) for r in obs if r['telluric_cut_change'] is not None),'observed_max_foreground_shift':max(abs(r['foreground_removal_change']) for r in obs),'observed_max_tilt_shift':max(abs(r[k]) for r in obs for k in ['tilt_plus_change','tilt_minus_change']),'primary_summaries':[r for r in summaries if r['target']=='observed' and r['clock_sigma_shift']==0 and not r['exclude_2012fr']],'scope':p['interpretation']}
 dump('result.json',answer);dump('manifest.json',{'files_sha256':{str(f.relative_to(ROOT)):sha(f) for f in sorted(P.iterdir()) if f.is_file() and f.name not in ['manifest.json','run.log']}});print(json.dumps(answer,indent=2))
if __name__=='__main__':{'prepare':prepare,'run':run}[sys.argv[1]]()
