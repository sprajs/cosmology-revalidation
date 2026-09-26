#!/usr/bin/env python3
"""Small, independent calibrated detector-image diagnostic; not a scene-model refit."""
from pathlib import Path
import json,numpy as np,pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import gaussian_filter,binary_dilation
R=Path(__file__).resolve().parents[3];O=R/'phase2/assumptions';I=O/'image_access';assert (O/'pixel_noise_preregistration.json').exists()
results=[];apertures=[]
def stat(x):
 x=np.asarray(x);m=np.median(x)
 return {'n':int(len(x)),'mean':float(np.mean(x)),'median':float(m),'rms':float(np.sqrt(np.mean(x*x))),'mad_sigma':float(1.4826*np.median(abs(x-m))),'abs_gt3_fraction':float(np.mean(abs(x)>3)),'abs_gt5_fraction':float(np.mean(abs(x)>5))} if len(x) else {'n':0}
for sci in sorted(I.glob('*.ext1.cutout.fits')):
 full=sci.with_name(sci.name.replace('.cutout.','.completeccd.'));data,hdr=fits.getdata(sci,header=True);fd,fh=fits.getdata(full,header=True);mask=fits.getdata(str(sci).replace('ext1','ext2'));weight=fits.getdata(str(sci).replace('ext1','ext3'))
 assert fd.shape==weight.shape==(4096,2048)
 ox=int(round(fh['CRPIX1']-hdr['CRPIX1']));oy=int(round(fh['CRPIX2']-hdr['CRPIX2']));ny,nx=data.shape
 exact=np.array_equal(fd[oy:oy+ny,ox:ox+nx],data);assert exact,'Science cutout does not match full CCD crop';w=weight[oy:oy+ny,ox:ox+nx]
 snx,sny=WCS(hdr).all_world2pix([[54.836567,-26.640186]],0)[0];yy,xx=np.indices(data.shape);outside=((xx-snx)**2+(yy-sny)**2)>40**2;valid=(mask==0)&(w>0)&np.isfinite(data)&outside
 baseline=float(np.median(data[valid]));sigma=float(hdr['SKYSIGMA']);res=data-baseline;pull=res*np.sqrt(np.maximum(w,0));smoothsigma=2.;kernel=np.zeros((41,41));kernel[20,20]=1;kernel=gaussian_filter(kernel,smoothsigma);smooth=gaussian_filter(res,smoothsigma);detected=abs(smooth)>4*sigma*np.sqrt(np.sum(kernel*kernel));screen=~binary_dilation(detected,iterations=4);core=valid&screen
 neighbor={}
 for axis in [0,1]:
  a=(slice(None,-1),slice(None)) if axis==0 else (slice(None),slice(None,-1));b=(slice(1,None),slice(None)) if axis==0 else (slice(None),slice(1,None));k=core[a]&core[b];l=pull[a][k];r=pull[b][k];neighbor[str(axis)]={'pairs':len(l),'normalized_product':float(np.sum(l*r)/np.sqrt(np.sum(l*l)*np.sum(r*r)))}
 aa=[]
 for cy in range(20,ny-12,30):
  for cx in range(20,nx-12,30):
   if (cx-snx)**2+(cy-sny)**2<40**2:continue
   rr=(xx-cx)**2+(yy-cy)**2;ap=rr<=9;an=(rr>=64)&(rr<=144);needed=ap|an;ok=(mask==0)&(w>0)&np.isfinite(data)
   if not np.all(ok[needed]):continue
   bg=float(np.median(data[an]));flux=float(np.sum(data[ap]-bg));var=float(np.sum(1/w[ap])+(ap.sum()**2)*np.pi/2*np.sum(1/w[an])/(an.sum()**2));row={'exposure':int(hdr['EXPNUM']),'band':hdr['FILTER'][0],'cx':cx,'cy':cy,'aperture_pixels':int(ap.sum()),'annulus_pixels':int(an.sum()),'background_electrons_per_pixel':bg,'aperture_flux_electrons':flux,'diagonal_sigma_electrons':np.sqrt(var),'pull':flux/np.sqrt(var),'source_screened':bool(np.all(screen[needed]))};aa.append(row);apertures.append(row)
 A=pd.DataFrame(aa);r={'exposure':int(hdr['EXPNUM']),'CCDNUM':int(hdr['CCDNUM']),'band':hdr['FILTER'][0],'MJD_OBS':float(hdr['MJD-OBS']),'SCI_BUNIT':hdr['BUNIT'],'shape':list(data.shape),'crop_origin_xy_zero_based':[ox,oy],'exact_science_crop_verified':exact,'SN_xy_zero_based':[float(snx),float(sny)],'median_weight_sigma_electrons':float(np.median(1/np.sqrt(w[valid]))),'header_SKYSIGMA_electrons':sigma,'sky_median_electrons':baseline,'pixel_pull_all_fixed_mask':stat(pull[valid]),'pixel_pull_source_screened':stat(pull[core]),'source_screened_neighbor_products':neighbor,'aperture_pull_all_fixed_controls':stat(A.pull),'aperture_pull_source_screened_controls':stat(A.loc[A.source_screened,'pull']),'caution':'Source-screened selection can truncate noise tails; all controls contain astrophysical contamination. One SN field and 3 exposures cannot establish survey-level calibration or reproduce SMP. Aperture variance includes annulus-median background estimate with Gaussian pi/2 approximation, but no measured pixel covariance.'};results.append(r)
 fits.PrimaryHDU(w,header=hdr).writeto(I/(sci.name.replace('ext1.cutout','weight.aligned_stamp')),overwrite=True)
pd.DataFrame(apertures).to_csv(O/'pixel_control_apertures.csv',index=False);(O/'pixel_noise_results.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(results,indent=2))
