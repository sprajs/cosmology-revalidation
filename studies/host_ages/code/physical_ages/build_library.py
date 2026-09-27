#!/usr/bin/env python3
"""Generate physical SSP spectra per solar mass formed, with native filter checks."""
import os,json,time,hashlib,argparse
from pathlib import Path
from datetime import datetime,timezone
import numpy as np

ROOT=Path(__file__).resolve().parents[4]
HERE=Path(__file__).parent
WORK=ROOT/'.work/physical-ages'
OUT=ROOT/'studies/host_ages/results/physical_ages'
os.environ['SPS_HOME']=str(WORK/'fsps')
import fsps

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--refined',action='store_true');args=parser.parse_args()
    plan=json.loads((HERE/'design.json').read_text());start=time.perf_counter()
    ages=np.asarray(plan['ages_gyr']);metals=np.asarray(plan['metallicity_logZ_solar'])
    if args.refined:
        ages=np.sort(np.r_[ages,np.sqrt(ages[:-1]*ages[1:])])
        metals=np.sort(np.r_[metals,(metals[:-1]+metals[1:])/2])
    pop=fsps.StellarPopulation(zcontinuous=1,sfh=0,imf_type=1,add_dust_emission=False,add_neb_emission=True,gas_logu=-2.,dust1=0,dust2=0)
    bands=['sdss_'+b for b in 'ugriz']
    specs=[];mass=[];mags=[]
    for metal in metals:
        pop.params['logzsol']=metal;pop.params['gas_logz']=metal
        ss=[];mm=[];native=[]
        for age in ages:
            wave,flux=pop.get_spectrum(tage=age,peraa=True)
            assert abs(pop.formed_mass-1)<1e-12
            ss.append(flux);mm.append(pop.stellar_mass)
            native.append([pop.get_mags(tage=age,bands=bands,redshift=z).tolist() for z in [0.,.1,.3]])
        specs.append(ss);mass.append(mm);mags.append(native)
        print('Completed metallicity',metal,flush=True)
    arrays=dict(wavelength_A=wave,spectra_Lsun_per_A_per_formed_Msun=np.array(specs),surviving_mass_fraction=np.array(mass),ages_Gyr=ages,metallicity=metals,native_redshifts=[0.,.1,.3],native_sdss_mags=np.array(mags))
    for band in bands:
        f=fsps.get_filter(band);arrays[band+'_wavelength']=f.transmission[0];arrays[band+'_throughput']=f.transmission[1]
    suffix='-refined' if args.refined else ''
    path=WORK/('ssp-library'+suffix+'.npz');np.savez_compressed(path,**arrays)
    versions={'numpy':np.__version__,'fsps':fsps.__version__}
    result=dict(created_utc=datetime.now(timezone.utc).isoformat(),seconds=time.perf_counter()-start,refined=args.refined,ages_gyr=ages.tolist(),metallicities=metals.tolist(),libraries=[x.decode() for x in pop.libraries],versions=versions,age_semantics='Spectra and mixture coefficients use total formed mass. Surviving stellar+remnant mass fractions are separate, never silently exchanged.',nebular='Cloudy table gas metallicity equals stellar metallicity; gas_logU=-2 fixed; dust emission disabled; foreground and host attenuation applied in forward photometry.',output=str(path.relative_to(ROOT)),output_sha256=sha(path),code_sha256=sha(Path(__file__)),design_sha256=sha(HERE/'design.json'),acquisition_sha256=sha(OUT/'acquisition.json'))
    (OUT/('stellar-library'+suffix+'.json')).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
