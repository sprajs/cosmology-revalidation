#!/usr/bin/env python3
"""Independent uniform-grid integration review of the physical-age operator."""
from pathlib import Path
import sys,json,hashlib,datetime
import numpy as np
from scipy.integrate import trapezoid
import extinction
ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/'studies/host_ages/code/physical_ages'))
from model import Library,Grid,C_A_S
R=ROOT/'.work/host-transport';S=ROOT/'.work/physical-ages';O=ROOT/'studies/host_ages/results/host_transport'

def main():
    results=[]
    for family in ['sdss','des']:
     lib=Library(response_path=R/'filters/des-deep-responses.npz') if family=='des' else Library()
     grid=Grid(lib)
     for z,ebv in [(.1,0.),(.6,.035),(1.05,0.)]:
      ages,a,_,_=grid.at_redshift(z,ebv)
      # Independently construct interpolated SSPs and attenuation from stored spectra.
      keep=lib.ages<=ages[-1];sp=lib.spectra[:,keep,:]
      if len(ages)>keep.sum():
       right=np.searchsorted(lib.ages,ages[-1]);frac=np.log(ages[-1]/lib.ages[right-1])/np.log(lib.ages[right]/lib.ages[right-1]);sp=np.concatenate([sp,((1-frac)*lib.spectra[:,right-1]+frac*lib.spectra[:,right])[:,None,:]],axis=1)
      for j in [0,9,16,34,51,67]:
       par=grid.parameters[j];iz=list(lib.data['metallicity']).index(par['metallicity']);indices=sorted(set([0,3,7,len(ages)-1]));nu=sp[iz,indices]*lib.wave**2/C_A_S*np.exp(-par['tau']*(lib.wave/5500)**par['power'])
       flux=[]
       for band in lib.bands:
        fw,ft=lib.filters[band];wave=np.arange(fw.min(),fw.max()+.05,.1);throughput=np.interp(wave,fw,ft,left=0,right=0);foreground=10**(-.4*extinction.fitzpatrick99(wave,3.1*ebv,3.1));observed=np.array([np.interp(wave/(1+z),lib.wave,row) for row in nu]);flux.append(trapezoid(observed*throughput*foreground/wave,wave,axis=1)/trapezoid(throughput/wave,wave))
       flux=np.array(flux);operator=a[j][:,indices];color_error=-2.5*np.log10((operator/operator[2])/(flux/flux[2]));k=np.unravel_index(abs(color_error).argmax(),color_error.shape);results.append({'family':family,'z':z,'ebv':ebv,'cell':j,'max_abs_colour_error_mag':float(abs(color_error).max()),'band':lib.bands[k[0]],'age_Gyr':float(ages[indices[k[1]]])})
    maximum=max(x['max_abs_colour_error_mag'] for x in results)
    assert maximum<1e-6, maximum
    r={'passed':True,'tolerance_mag':1e-6,'checked_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'method':'Independent0.1Angstrom uniform observer grid trapezoid integral, native stored SSPs and independently constructed clock-interpolated ceiling. Compare photon-weighted AB colours against Grid operator over36 dust/metal/z/foreground cases and4 ages per case.','cases':results,'max_abs_colour_error_mag':max(x['max_abs_colour_error_mag'] for x in results),'review_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [S/'ssp-library.npz',R/'filters/des-deep-responses.npz']},'model_sha256':hashlib.sha256((ROOT/'studies/host_ages/code/physical_ages/model.py').read_bytes()).hexdigest()}
    (O/'physical-operator-review.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))

if __name__=="__main__":main()
