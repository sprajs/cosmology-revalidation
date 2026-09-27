"""Signed DESI native-pixel indices with joint formal statistical covariance."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from acquire import ROOT, WORK, RESULT, CODE, sha

BANDS=[('Dn_blue',3850.,3950.),('Dn_red',4000.,4100.),('Hd_blue',4041.60,4079.75),('Hd_feature',4083.50,4122.25),('Hd_red',4128.50,4161.)]


def native_vector(arrays):
    vectors=[]
    for b,lo,hi in [('B',-np.inf,5780.),('R',5780.,7570.),('Z',7570.,np.inf)]:
        w=arrays[b+'_WAVELENGTH'].astype(float);keep=(w>=lo)&(w<hi)
        vectors.append([w[keep],arrays[b+'_FLUX'][keep].astype(float),arrays[b+'_IVAR'][keep].astype(float),arrays[b+'_MASK'][keep],np.full(keep.sum(),b)])
    return [np.concatenate([x[j]for x in vectors]) for j in range(5)]


def measure(wave,flux,ivar,mask,arm,z):
    rest=wave/(1+z)
    # Native 0.8A sampling: each camera has contiguous pixel centers.
    assert np.allclose(np.diff(wave),.8,rtol=0,atol=1e-6)
    half=.4/(1+z)
    valid=np.isfinite(flux)&np.isfinite(ivar)&(ivar>0)&(mask==0)
    variance=np.zeros_like(ivar);variance[valid]=1/ivar[valid]
    safe=np.where(valid,flux,0.)
    weights=[];cover=[];widths=[]
    for k,lo,hi in BANDS:
        overlap=np.maximum(0.,np.minimum(rest+half,hi)-np.maximum(rest-half,lo))*valid
        width=overlap.sum();cover.append(float(width/(hi-lo)));widths.append(overlap)
        weight=overlap/width if width else overlap
        if k.startswith('Dn'):weight=weight*wave**2/2.99792458e18*1e12 # microJy per native flux unit
        weights.append(weight)
    weights=np.array(weights);means=weights@safe;cov=(weights*variance)@weights.T
    support=np.array(cover)>=.95
    se=np.sqrt(cov.diagonal());snr=np.divide(means,se,out=np.full(5,np.nan),where=se>0)
    result={'band_order':[x[0]for x in BANDS],'band_units':['microJy','microJy','1e-17 erg/s/cm2/A','1e-17 erg/s/cm2/A','1e-17 erg/s/cm2/A'],'band_flux':means.tolist(),'band_covariance_formal_diagonal_ivar':cov.tolist(),'band_coverage':cover,'band_snr_nominal':snr.tolist(),'Dn4000':None,'Dn4000_formal_sigma':None,'Hdelta_native_A':None,'Hdelta_formal_sigma_A':None,'index_covariance_formal':None}
    gradients=[];values=[]
    if support[:2].all() and snr[0]>5:
        dn=means[1]/means[0]
        gdn=weights[1]/means[0]-weights[0]*means[1]/means[0]**2
        result['Dn4000']=float(dn);result['Dn4000_formal_sigma']=float(np.sqrt(np.dot(gdn*gdn,variance)))
        gradients.append(gdn);values.append(dn)
    else:gradients.append(None)
    # Continuum anchors are fixed band midpoints, not fitted from absorption features.
    t=(rest-(4041.60+4079.75)/2)/((4128.50+4161.)/2-(4041.60+4079.75)/2)
    continuum=means[2]*(1-t)+means[4]*t
    feature=widths[3]>0
    if support[2:].all() and snr[2]>5 and snr[4]>5 and np.all(continuum[feature]>0):
        ew=float(np.sum(widths[3][feature]*(1-safe[feature]/continuum[feature])))
        gew=np.zeros_like(flux);gew[feature]=-widths[3][feature]/continuum[feature]
        dB=np.sum(widths[3][feature]*safe[feature]/continuum[feature]**2*(1-t[feature]))
        dR=np.sum(widths[3][feature]*safe[feature]/continuum[feature]**2*t[feature])
        gew+=dB*weights[2]+dR*weights[4]
        result['Hdelta_native_A']=ew;result['Hdelta_formal_sigma_A']=float(np.sqrt(np.dot(gew*gew,variance)))
        gradients.append(gew);values.append(ew)
    else:gradients.append(None)
    if all(g is not None for g in gradients):
        g=np.array(gradients);result['index_covariance_formal']=((g*variance)@g.T).tolist()
    relevant=(rest>=3850)&(rest<=4161)&valid
    pixel_snr=flux[relevant]*np.sqrt(ivar[relevant])
    result['relevant_pixels']=int(relevant.sum())
    result['signed_pixel_snr_percentiles']=np.percentile(pixel_snr,[0,5,50,95,100]).tolist() if len(pixel_snr) else None
    result['pixels_abs_snr_ge20']=int((abs(pixel_snr)>=20).sum())
    result['pixels_abs_snr_ge30']=int((abs(pixel_snr)>=30).sum())
    result['negative_flux_pixels']=int((flux[relevant]<0).sum())
    result['cameras_in_index_region']=sorted(set(arm[relevant]))
    result['gradient_arrays']=gradients
    return result


def main():
    records=[]
    for p in sorted((WORK/'spectra').glob('*.json')):
        meta=json.loads(p.read_text());arrays=np.load(ROOT/meta['array_file'])
        assert sha(ROOT/meta['array_file'])==meta['array_sha256']
        v=native_vector(arrays);r=measure(*v,meta['desi_z']);r.pop('gradient_arrays')
        r.update(targetid=meta['targetid'],oz_OzDES_ID=meta['oz_OzDES_ID'],z=meta['desi_z'],native_coordinate_gate=meta['native_coordinate_gate'],native_EBV=meta['fibermap'].get('EBV'),coadd_numexp=meta['fibermap']['COADD_NUMEXP'])
        # A declared response tilt is a sensitivity, not an empirical calibration error.
        tilted=[]
        for sign in [-1,1]:
            change=np.exp(sign*.02*(v[0]/(1+meta['desi_z'])-4000)/100)
            rr=measure(v[0],v[1]*change,v[2]/change**2,v[3],v[4],meta['desi_z'])
            tilted.append({k:rr[k]for k in ['Dn4000','Hdelta_native_A']})
        r['response_tilt_sensitivity_plus_minus_2percent_per_100_rest_A']=tilted
        records.append(r)
    rp=WORK/'band-records.json';rp.write_text(json.dumps(records,indent=2,allow_nan=False)+'\n')
    dn=[r for r in records if r['Dn4000'] is not None];hd=[r for r in records if r['Hdelta_native_A'] is not None]
    summary={'schema':'calibrated-host-band-measurements-v1','status':'complete','spectra':len(records),'Dn4000_supported':len(dn),'Hdelta_supported':len(hd),'both_indices_supported':sum(r['index_covariance_formal'] is not None for r in records),'hosts_with_any_abs_pixel_snr_ge20':sum(r['pixels_abs_snr_ge20']>0 for r in records),'hosts_with_any_abs_pixel_snr_ge30':sum(r['pixels_abs_snr_ge30']>0 for r in records),'spectra_with_negative_flux_pixels_in_index_region':sum(r['negative_flux_pixels']>0 for r in records),'multi_camera_index_regions':sum(len(r['cameras_in_index_region'])>1 for r in records),'multiple_exposure_coadds':sum(r['coadd_numexp']>1 for r in records),'Dn4000_range':[min(r['Dn4000']for r in dn),max(r['Dn4000']for r in dn)] if dn else None,'Hdelta_native_A_range':[min(r['Hdelta_native_A']for r in hd),max(r['Hdelta_native_A']for r in hd)] if hd else None,'interpretation':'Observed host spectral indices, not independently measured progenitor ages. Formal errors/covariance omit calibration, DR1 variance bug and residual interpixel terms. No population age or cosmological correction fit.','code_sha256':{str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),CODE/'design.json']},'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in sorted((WORK/'spectra').glob('*'))},'output_sha256':{str(rp.relative_to(ROOT)):sha(rp)}}
    (RESULT/'bands.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items()if not k.endswith('sha256')},indent=2))


if __name__=='__main__':main()
