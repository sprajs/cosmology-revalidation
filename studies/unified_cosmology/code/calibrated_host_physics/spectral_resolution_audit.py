"""Read-only inventory of native DESI spectra, no age fits or new FSPS calls."""
import os
os.environ['OMP_NUM_THREADS']='1';os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import json,hashlib,sys,importlib.metadata
import numpy as np
from scipy import sparse
ROOT=Path(__file__).resolve().parents[4];INPUT=ROOT/'.work/unified-cosmology/calibrated-hosts/spectra';OUT=ROOT/'.work/unified-cosmology/full-spectrum-design';OUT.mkdir(parents=True,exist_ok=True)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
libpath=ROOT/'.work/unified-cosmology/calibrated-host-physics/ssp-library.npz'
origin=ROOT/'studies/unified_cosmology/results/calibrated_hosts/spectra.json'
origin_info=json.loads(origin.read_text())
for name,expected in origin_info['retrieved_files_sha256'].items():assert sha(ROOT/name)==expected,name
lib_origin=ROOT/'studies/unified_cosmology/results/calibrated_host_physics/stellar-library.json'
assert sha(libpath)==json.loads(lib_origin.read_text())['output_sha256'][str(libpath.relative_to(ROOT))]
lib=np.load(libpath);lw=lib['wavelength_A'];lr=lib['resolution_sigma_km_s']
regions={'common_3700_4900':(3700.,4900.),'old_indices_3850_4161':(3850.,4161.),'Hgamma_G4300':(4250.,4450.),'Hbeta':(4800.,4910.),'Mgb':(5100.,5250.),'Fe5270_5335':(5250.,5400.),'Halpha':(6500.,6630.)}
rows=[];sources={str(libpath.relative_to(ROOT)):sha(libpath)}
for mp in sorted(INPUT.glob('*.json')):
 m=json.loads(mp.read_text());p=ROOT/m['array_file'];assert sha(p)==m['array_sha256'];a=np.load(p);z=m['desi_z'];parts=[]
 for cam,lo,hi in [('B',-np.inf,5780.),('R',5780.,7570.),('Z',7570.,np.inf)]:
  wave=a[cam+'_WAVELENGTH'];rest=wave/(1+z);flux=a[cam+'_FLUX'];ivar=a[cam+'_IVAR'];mask=a[cam+'_MASK'];d=a[cam+'_RESOLUTION'].astype(float)
  valid=(wave>=lo)&(wave<hi)&(mask==0)&np.isfinite(flux)&np.isfinite(ivar)&(ivar>0)
  n=len(wave);k=d.shape[0]//2;R=sparse.dia_matrix((d,np.arange(k,-k-1,-1)),shape=(n,n)).tocsr()
  sums=np.asarray(R.sum(axis=0)).ravel();log=np.log(wave)-np.median(np.log(wave));mean=np.divide(R.T@log,sums,out=np.zeros(n),where=sums>0);var=np.divide(R.T@(log*log),sums,out=np.zeros(n),where=sums>0)-mean*mean
  sig=np.where((var>=0)&(sums>0),299792.458*np.sqrt(np.maximum(var,0)),np.nan);intr=np.interp(rest,lw,lr)
  parts.append(dict(cam=cam,wave=wave,rest=rest,flux=flux,ivar=ivar,valid=valid,sigma=sig,intrinsic=intr,sums=sums))
 row={'targetid':m['targetid'],'ozdes_id':m['oz_OzDES_ID'],'z':z,'coadd_numexp':m['fibermap']['COADD_NUMEXP'],'regions':{}}
 for name,(l,h) in regions.items():
  fields={k:[] for k in ['valid','flux','ivar','sigma','intrinsic','rest']};available=0
  for part in parts:
   sel=(part['rest']>=l)&(part['rest']<h)&((part['wave']<5780) if part['cam']=='B' else ((part['wave']>=5780)&(part['wave']<7570)) if part['cam']=='R' else(part['wave']>=7570))
   available+=int(sel.sum())
   for key in fields:fields[key].extend(part[key][sel].tolist())
  fields={k:np.array(v) for k,v in fields.items()};v=fields['valid'].astype(bool);f=fields['flux'][v];iv=fields['ivar'][v];ss=fields['sigma'][v];si=fields['intrinsic'][v]
  native_step=.8/(1+z);totalwidth=native_step*available;validwidth=native_step*v.sum()
  row['regions'][name]={'native_pixels_in_region':available,'valid_pixels':int(v.sum()),'approx_native_geometric_coverage':min(1.,float(totalwidth/(h-l))),'approx_valid_coverage':min(1.,float(validwidth/(h-l))),'signed_pixel_SNR_median':None if not len(f) else float(np.median(f*np.sqrt(iv))),'signed_mean_SNR':None if not len(f) else float(f.sum()/np.sqrt((1/iv).sum())),'max_abs_pixel_SNR':None if not len(f) else float(max(abs(f*np.sqrt(iv)))),'negative_flux_pixels':int((f<0).sum()),'instrument_sigma_kms_percentiles':None if not len(ss) else np.nanpercentile(ss,[0,5,50,95,100]).tolist(),'library_sigma_kms_percentiles':None if not len(si) else np.percentile(si,[0,50,100]).tolist(),'undefined_signed_response_second_moments':int((~np.isfinite(ss)).sum()),'instrument_narrower_than_library_valid_pixels':int((ss<si).sum()),'instrument_narrower_fraction':None if not len(ss) else float(np.sum(ss<si)/np.isfinite(ss).sum())}
 sources[str(mp.relative_to(ROOT))]=sha(mp);sources[str(p.relative_to(ROOT))]=sha(p);rows.append(row)
# Additional mask diagnostics use all available FSPS nebular centres, rather than
# the old Library.lines subset restricted to five diagnostic bands.
all_lines=lib['emission_line_wavelength_A']; extra_rows=[]
for row in rows:
 a=np.load(INPUT/(row['targetid']+'.npz'));z=row['z'];before=after=high_before=high_after=0;neg_elems=all_elems=0
 for cam,lo,hi in [('B',-np.inf,5780.),('R',5780.,7570.),('Z',7570.,np.inf)]:
  w=a[cam+'_WAVELENGTH'];r=w/(1+z);f=a[cam+'_FLUX'];iv=a[cam+'_IVAR'];d=a[cam+'_RESOLUTION'].astype(float)
  good=(w>=lo)&(w<hi)&(r>=3700)&(r<4900)&(a[cam+'_MASK']==0)&np.isfinite(f)&np.isfinite(iv)&(iv>0)
  keep=good.copy()
  for line in all_lines[(all_lines>3690)&(all_lines<4910)]:keep &= np.abs(np.log(r/line))*299792.458>400
  snr=np.abs(f*np.sqrt(np.maximum(iv,0)))
  before+=int(good.sum());after+=int(keep.sum());high_before+=int(np.sum(good&(snr>=20)));high_after+=int(np.sum(keep&(snr>=20)))
  neg_elems+=int(np.sum(d[:,good]<0));all_elems+=int(d[:,good].size)
 row['common_nebular_mask400']={'valid_pixels_before':before,'retained_pixels_after':after,'abs_SNR_ge20_before':high_before,'abs_SNR_ge20_after':high_after,'resolution_negative_coefficient_fraction':neg_elems/all_elems}
summary={}
for name in regions:
 r=[x['regions'][name] for x in rows];snr=[x['signed_pixel_SNR_median'] for x in r if x['signed_pixel_SNR_median'] is not None]
 summary[name]={'hosts_with_geometric_coverage95':sum(x['approx_native_geometric_coverage']>=.95 for x in r),'hosts_with_valid_coverage95':sum(x['approx_valid_coverage']>=.95 for x in r),'hosts_with_valid_coverage95_and_signed_mean_SNR5':sum(x['approx_valid_coverage']>=.95 and x['signed_mean_SNR']>=5 for x in r),'valid_pixel_count_min_median_max':np.percentile([x['valid_pixels'] for x in r],[0,50,100]).tolist(),'median_signed_pixel_SNR_min_median_max':np.percentile(snr,[0,50,100]).tolist(),'hosts_with_any_instrument_narrower_than_library':sum(x['instrument_narrower_than_library_valid_pixels']>0 for x in r),'narrower_fraction_all_valid_pixels':sum(x['instrument_narrower_than_library_valid_pixels'] for x in r)/sum(x['valid_pixels']-x['undefined_signed_response_second_moments'] for x in r),'undefined_signed_response_second_moments':sum(x['undefined_signed_response_second_moments'] for x in r),'max_abs_pixel_SNR':max(x['max_abs_pixel_SNR'] for x in r if x['max_abs_pixel_SNR'] is not None)}
result={'status':'input_availability_and_response_moment_inventory_only_no_fits','hosts':len(rows),'multi_exposure':sum(x['coadd_numexp']>1 for x in rows),'redshift_range':[min(x['z'] for x in rows),max(x['z'] for x in rows)],'library_shape':list(lib['spectra_Lsun_per_A_per_formed_Msun'].shape),'library_wavelength_min_max':[float(min(lw)),float(max(lw))],'library_sigma_optical':np.unique(lr[(lw>3700)&(lw<6800)]).tolist(),'regions':summary,'all_FSPS_line_centres_common_region':all_lines[(all_lines>3690)&(all_lines<4910)].tolist(),'nebular_mask400':{'retained_pixels_min_median_max':np.percentile([x['common_nebular_mask400']['retained_pixels_after'] for x in rows],[0,50,100]).tolist(),'hosts_with_abs_SNR_ge20_before':sum(x['common_nebular_mask400']['abs_SNR_ge20_before']>0 for x in rows),'hosts_with_abs_SNR_ge20_after':sum(x['common_nebular_mask400']['abs_SNR_ge20_after']>0 for x in rows),'pixels_abs_SNR_ge20_before':sum(x['common_nebular_mask400']['abs_SNR_ge20_before'] for x in rows),'pixels_abs_SNR_ge20_after':sum(x['common_nebular_mask400']['abs_SNR_ge20_after'] for x in rows)},'rows':rows,'input_sha256':sources,'source_sha256':{str(Path(__file__).relative_to(ROOT)):sha(__file__)},'limits':['Coverage uses native0.8A pixels and fixed existing camera boundaries; no imputation.','SNR is a signed formal-IVAR measurement; not a host selection or age detection.','Instrument sigma is a signed-column second moment of native R; negative response wings are preserved, not a fitted Gaussian LSF or full matching operator.','No emission masking beyond release mask in this availability inventory; feature presence does not imply stellar absorption is identified.']}
(OUT/'input-inventory.json').write_text(json.dumps(result,indent=2)+'\n');print('Native inventory complete:',len(rows),'hosts',flush=True)

# No age fit: test whether a covariance-preserving response remedy is numerically
# adequate for the fixed available family and a separate narrow-line stress.
from scipy.linalg import solve_triangular
import extinction
DESIGN=Path(__file__).with_name('full-spectrum-audit-design.json')
design=json.loads(DESIGN.read_text());windows=np.array(design['windows_vacuum_rest_A'],dtype=float)
reference=design['reference_window_index'];sigma=design['kernel_sigma_km_s']/299792.458
stellar=lib['spectra_Lsun_per_A_per_formed_Msun'].reshape(-1,len(lw))
dust=[(0.,-.7)]+[(tau,power) for tau in [.5,1.,2.] for power in [-.7,-1.]]

def blur(wave):
    # Native observed bins are uniform0.8A; retain the explicit quadrature
    # normalization rather than normalizing Gaussian rows after truncation.
    step=float(np.median(np.diff(wave)));n=len(wave)
    radius=int(np.ceil(design['kernel_tail_sigma']*sigma*wave[-1]/step))+1
    rr=[];cc=[];vv=[]
    for off in range(-radius,radius+1):
        ii=np.arange(max(0,-off),min(n,n-off));jj=ii+off
        x=np.log(wave[ii]/wave[jj])/sigma;keep=np.abs(x)<=design['kernel_tail_sigma']
        ii,jj,x=ii[keep],jj[keep],x[keep]
        rr.extend(ii);cc.extend(jj);vv.extend(step/wave[ii]/(np.sqrt(2*np.pi)*sigma)*np.exp(-.5*x*x))
    return sparse.csr_matrix((vv,(rr,cc)),shape=(n,n))

remedy=[]
for record in rows:
    target=record['targetid'];meta=json.loads((INPUT/(target+'.json')).read_text());arrays=np.load(ROOT/meta['array_file']);z=meta['desi_z'];parts=[]
    for cam,lo,hi in [('B',-np.inf,5780.),('R',5780.,7570.),('Z',7570.,np.inf)]:
        wave=arrays[cam+'_WAVELENGTH'].astype(float);rest=wave/(1+z)
        flux=arrays[cam+'_FLUX'].astype(float);ivar=arrays[cam+'_IVAR'].astype(float)
        valid=(arrays[cam+'_MASK']==0)&np.isfinite(flux)&np.isfinite(ivar)&(ivar>0)
        for line in all_lines[(all_lines>rest[0]*.997)&(all_lines<rest[-1]*1.003)]:
            valid &= np.abs(np.log(rest/line))*299792.458>400.
        G=blur(wave);support=G.copy();support.data[:]=1.
        usable=np.asarray(support@(~valid).astype(float)).ravel()==0
        usable &= (wave>wave[0]*np.exp(design['kernel_tail_sigma']*sigma))&(wave<wave[-1]*np.exp(-design['kernel_tail_sigma']*sigma))
        cuts=(wave>=lo)&(wave<hi)
        geom=np.array([np.maximum(0.,np.minimum(rest+.4/(1+z),h)-np.maximum(rest-.4/(1+z),l))*cuts for l,h in windows])
        widths=geom*usable
        d=arrays[cam+'_RESOLUTION'].astype(float);k=d.shape[0]//2
        R=sparse.dia_matrix((d,np.arange(k,-k-1,-1)),shape=(len(wave),len(wave))).tocsr()
        f=np.array([np.interp(rest,lw,s) for s in stellar]).T
        fore=10**(-.4*extinction.fitzpatrick99(wave,3.1*meta['fibermap']['EBV']*.86,3.1))/(1+z)
        parts.append(dict(wave=wave,rest=rest,flux=np.where(valid,flux,0.),variance=np.divide(1.,ivar,out=np.zeros_like(ivar),where=valid),G=G,R=R,f=f,fore=fore,widths=widths,geom=geom))
    retained=sum(p['widths'].sum(axis=1) for p in parts);coverage=sum(p['geom'].sum(axis=1) for p in parts)/np.diff(windows,axis=1).ravel()
    keep=(retained>0)&(coverage>=.95);indices=np.flatnonzero(keep);ref=int(np.flatnonzero(indices==reference)[0])
    y=np.zeros(len(indices));C=np.zeros((len(indices),len(indices)));left=np.zeros((len(indices),len(stellar)*len(dust)));right=left.copy();direct_cov=np.zeros_like(C)
    stress_ops=[]
    for p in parts:
        W=p['widths'][keep]/retained[keep,None];WG=np.asarray(p['G'].T@W.T).T
        WR=np.asarray(p['R'].T@W.T).T;WGR=np.asarray(p['R'].T@WG.T).T
        y+=WG@p['flux'];C+=(WG*p['variance'])@WG.T
        # Independent rank-one accumulation across native pixels verifies every
        # induced cross-band covariance entry; no diagonal substitute.
        direct_cov+=np.einsum('ik,k,jk->ij',WG,p['variance'],WG,optimize=False)
        for j,(tau,power) in enumerate(dust):
            attenuation=p['fore']*np.exp(-tau*(p['rest']/5500.)**power)
            A=WGR*attenuation
            B=np.asarray(p['G'].T@(WR*attenuation).T).T
            sl=slice(j*len(stellar),(j+1)*len(stellar))
            left[:,sl]+=A@p['f'];right[:,sl]+=B@p['f']
            # Scale line stress at observed reference continuum; factor below
            # removes the smooth screen at its origin, so the stated EW is
            # observed rather than intrinsic. It is only a numerical stress.
            stress_ops.append((A-B)/attenuation[None,:]/float(np.median(np.diff(p['wave']))))
    chol=np.linalg.cholesky(C);diag=np.sqrt(np.diag(C));assert y[ref]>0 and np.all(right[ref]>0)
    amplitude=y[ref]/right[ref];error=(left-right)*amplitude[None,:]
    white=solve_triangular(chol,error,lower=True);single=np.max(abs(error)/diag[:,None]);joint=np.max(np.linalg.norm(white,axis=0))
    stress_joint=0.;stress_single=0.
    for op in stress_ops:
        delta=op*y[ref]
        stress_joint=max(stress_joint,float(np.max(np.linalg.norm(solve_triangular(chol,delta,lower=True),axis=0))))
        stress_single=max(stress_single,float(np.max(abs(delta)/diag[:,None])))
    correlation=C/np.outer(diag,diag);np.fill_diagonal(correlation,0)
    eig=np.linalg.eigvalsh(C/np.outer(diag,diag));cov_error=np.max(abs(C-direct_cov))/np.max(abs(C))
    assert cov_error<1e-11 and eig[0]>0
    remedy.append({'targetid':target,'windows_retained':indices.tolist(),'geometric_coverage':coverage.tolist(),'retained_rest_width_A':retained.tolist(),'signed_broad_flux':y.tolist(),'full_broad_covariance':C.tolist(),'test_columns':left.shape[1],'minimum_correlation_eigenvalue':float(eig[0]),'maximum_abs_offdiagonal_correlation':float(np.max(abs(correlation))),'independent_covariance_relative_error':float(cov_error),'maximum_column_error_per_band_sigma':float(single),'maximum_column_joint_whitened_error':float(joint),'finite_family_gate_passed':bool(single<=.1 and joint<=.1),'unit_observed_EW_line_max_joint_error':stress_joint,'unit_observed_EW_line_max_single_error':stress_single,'age_fit_performed':False})
    print(target,'operator sigma',float(single),'joint',float(joint),'line1A',stress_joint,flush=True)
ledger=OUT/'response-remedy-ledger.json';ledger.write_text(json.dumps(remedy,indent=2)+'\n')
inventory=OUT/'input-inventory.json'
summary={k:v for k,v in result.items() if k not in ['rows','input_sha256','source_sha256']}
summary.update(status='input_and_cross_convolution_audit_no_age_fit',versions={name:importlib.metadata.version(name) for name in ['numpy','scipy','extinction']},python_version=sys.version,origin_sha256={str(p.relative_to(ROOT)):sha(p) for p in [origin,lib_origin,Path(__file__).with_name('requirements-lock.txt')]},design_sha256=sha(DESIGN),source_sha256={str(Path(__file__).relative_to(ROOT)):sha(__file__)},input_sha256=result['input_sha256'],output_sha256={str(p.relative_to(ROOT)):sha(p) for p in [inventory,ledger]},response_remedy={'hosts':len(remedy),'finite_family_columns_per_host':490,'passed_hosts':sum(r['finite_family_gate_passed'] for r in remedy),'maximum_error_per_band_sigma':max(r['maximum_column_error_per_band_sigma'] for r in remedy),'maximum_joint_whitened_error':max(r['maximum_column_joint_whitened_error'] for r in remedy),'hosts_1A_line_stress_above_same_01_joint_threshold':sum(r['unit_observed_EW_line_max_joint_error']>.1 for r in remedy),'worst_1A_observed_EW_line_joint_error':max(r['unit_observed_EW_line_max_joint_error'] for r in remedy),'maximum_cross_band_correlation':max(r['maximum_abs_offdiagonal_correlation'] for r in remedy),'independent_covariance_max_relative_error':max(r['independent_covariance_relative_error'] for r in remedy),'interpretation':'Finite-family response-order validation only, with all failure counts preserved. It does not validate the actual unknown library kernel, unmodelled sharper features, calibration uncertainty, coadd response truth or stellar-age identification.'})
public=ROOT/'studies/unified_cosmology/results/calibrated_host_physics/full-spectrum-audit.json'
public.write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary['response_remedy'],indent=2),flush=True)
