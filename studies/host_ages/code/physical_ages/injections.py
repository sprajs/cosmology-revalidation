#!/usr/bin/env python3
"""Recovery from known SFHs, including populations outside the fitted family."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import json,hashlib,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd
from scipy.stats import norm
import fit
from model import Library, Grid, fieller
ROOT=Path(__file__).resolve().parents[4];WORK=ROOT/'.work/physical-ages';OUT=ROOT/'studies/host_ages/results/physical_ages'

def truths(sample):
    response=ROOT/'.work/host-transport/filters/des-deep-responses.npz' if sample=='des' else None
    z=.6 if sample=='des' else .07
    ebv=0. if sample=='des' else .025
    lib=Library(response_path=response);grid=Grid(lib);ages,a,_,_=grid.at_redshift(z,ebv);dn=grid.observed_dn(z,ebv)
    weights=np.zeros(len(ages));weights[np.argmin(abs(ages-.1))]=.1;weights[np.argmin(abs(ages-4))]=.9
    # One exact model cell: nontrivial SFH, fixed known nuisance.
    examples=[dict(kind='included_family',mean=a[38]@weights,dn=float(dn[38,1]@weights/(dn[38,0]@weights)),age=float(ages@weights))]
    refined=Library(WORK/'ssp-library-refined.npz',response)
    design=json.loads(Path(__file__).with_name('design.json').read_text());design['host_dust_tauV']=[.3]
    fine=Grid(refined,design);fages,fa,_,_=fine.at_redshift(z,ebv);fdn=fine.observed_dn(z,ebv)
    metal=next(i for i,p in enumerate(fine.parameters) if p['metallicity']==-.25 and p['power']==-.7)
    index=int(np.argmin(abs(fages-np.sqrt(2))))
    examples.append(dict(kind='native_off_grid_SSP',mean=fa[metal,:,index],dn=float(fdn[metal,1,index]/fdn[metal,0,index]),age=float(fages[index])))
    # Two ages experience different dust AND metallicity; no common-screen
    # fitting cell can exactly represent this general population by construction.
    temp=Grid(lib);ia=int(np.argmin(abs(lib.ages-.02)));ib=int(np.argmin(abs(lib.ages-6)))
    spectra=.15*lib.spectra[2,ia]*np.exp(-3.*(lib.wave/5500.)**-1.)+.85*lib.spectra[0,ib]*np.exp(-.2*(lib.wave/5500.)**-.7)
    mean=lib.photometry(spectra[None,:],z,ebv)[:,0]
    # Evaluate the same direct spectral-band integral using this new truth.
    temp.nu_flat=(spectra*lib.wave**2/2.99792458e18)[None,:]
    # Direct rest-band quadrature includes the foreground in observed wavelength.
    import extinction
    band=[]
    for lo,hi in [(3850.,3950.),(4000.,4100.)]:
        wave=np.linspace(lo,hi,1001);fnu=np.interp(wave,lib.wave,spectra*lib.wave**2/2.99792458e18)
        fnu*=10**(-.4*extinction.fitzpatrick99(wave*(1+z),3.1*ebv,3.1));band.append(np.trapz(fnu,wave)/(hi-lo))
    examples.append(dict(kind='differential_dust_and_metallicity',mean=mean,dn=band[1]/band[0],age=.15*lib.ages[ia]+.85*lib.ages[ib]))
    return z,ebv,examples

def worker(job):
    sample,item=job;fit.initialize(sample)
    return fit.fit_one(item)

def main():
    start=time.perf_counter();rng=np.random.default_rng(27912026);jobs=[];truth_record=[]
    for sample in ['sdss','des']:
        z,ebv,examples=truths(sample)
        for example in examples:
            mean=np.asarray(example['mean']);mean/=np.median(mean)
            sigma=np.maximum(.02*mean,.002*mean.max())
            truth_record.append(dict(sample=sample,kind=example['kind'],age_Gyr=example['age'],Dn4000=example['dn'],redshift=z))
            for k in range(40):
                flux=mean+rng.normal(size=len(mean))*sigma
                item=dict(id=f'{sample}-{example["kind"]}-{k}',sample=sample,generator=example['kind'],truth_age_Gyr=example['age'],z=z,ebv=ebv,flux=flux.tolist(),error=sigma.tolist())
                if sample=='sdss':
                    red,blue=np.array([example['dn'],1.])+rng.normal(size=2)*.015
                    for level in [.95,.975]:
                        item['dn_interval'+('95' if level==.95 else '975')]=fieller(red,blue,.015**2,.015**2,0.,norm.ppf((1+level)/2))
                jobs.append((sample,item))
    rows=[]
    with ProcessPoolExecutor(4) as pool:
        for i,result in enumerate(pool.map(worker,jobs,chunksize=1)):
            rows+=result
            if (i+1)%40==0:print(f'{i+1}/{len(jobs)} injections',flush=True)
    table=pd.DataFrame(rows);path=WORK/'physical-age-injections.csv';table.to_csv(path,index=False);summary=[]
    for keys,group in table.groupby(['sample','generator','model_discrepancy_mag']):
        truth=group.truth_age_Gyr;covered=(group.age_low_Gyr<=truth)&(group.age_high_Gyr>=truth)
        joint=(group.joint_age_low_Gyr<=truth)&(group.joint_age_high_Gyr>=truth)
        summary.append(dict(sample=keys[0],generator=keys[1],model_discrepancy_mag=keys[2],draws=len(group),photometry_feasible=int(group.age_low_Gyr.notna().sum()),
            truth_covered=int(covered.sum()),joint_feasible=int(group.joint_age_low_Gyr.notna().sum()),joint_truth_covered=int(joint.sum()),numerical_failures=int(group.numerical_failures.sum())))
    dependencies=[Path(__file__).with_name(name) for name in ['model.py','fit.py','design.json']]+[WORK/'ssp-library.npz',WORK/'ssp-library-refined.npz',ROOT/'.work/host-transport/filters/des-deep-responses.npz']
    result=dict(seed=27912026,draws=240,truths=truth_record,summaries=summary,seconds=time.perf_counter()-start,
        dependencies_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in dependencies},
        status='Targeted post-observation model-family checks, not a blind independent validation sample. Forty draws per generator are a Monte Carlo diagnostic, not a precision coverage certificate.',
        output=str(path.relative_to(ROOT)),output_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (OUT/'physical-injections.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
