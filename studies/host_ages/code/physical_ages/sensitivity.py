#!/usr/bin/env python3
"""Refinement, age-definition, foreground and cosmological-clock sensitivity.

The subset is selected by SHA256(ID), never by a fitted age or residual.
These comparisons diagnose assumptions; none is a selected preferred model.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import json,time,hashlib,copy
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd
from astropy.cosmology import FlatLambdaCDM
import model
import fit
ROOT=Path(__file__).resolve().parents[4];WORK=ROOT/'.work/physical-ages';OUT=ROOT/'studies/host_ages/results/physical_ages'

class SurvivingGrid(model.Grid):
    def at_redshift(self,z,ebv):
        age,a,dn,mass=super().at_redshift(z,ebv);self.current_mass=mass
        return age,a/mass[:,None,:],dn/mass[:,None,:],mass
    def observed_dn(self,z,ebv):
        return super().observed_dn(z,ebv)/self.current_mass[:,None,:]

def worker(job):
    sample,variant,item=job
    response=ROOT/'.work/host-transport/filters/des-deep-responses.npz' if sample=='des' else None
    path=WORK/('ssp-library-refined.npz' if variant=='refined_grid' else 'ssp-library.npz')
    lib=model.Library(path,response)
    design=json.loads(Path(__file__).with_name('design.json').read_text())
    if variant=='refined_grid':
        tau=np.array(design['host_dust_tauV']);design['host_dust_tauV']=np.sort(np.r_[tau,(tau[:-1]+tau[1:])/2]).tolist()
    model.COSMO=FlatLambdaCDM(H0=70,Om0=1. if variant=='nonaccelerating_EdS_clock' else .3)
    fit.GRID=(SurvivingGrid if variant=='surviving_mass_age' else model.Grid)(lib,design)
    item=copy.deepcopy(item)
    if variant=='SFD_scale_1':item['ebv']/=.86
    if variant=='native_SDSS_zero_points':
        ab=10**(-.4*np.array([-.04,0,0,0,.02]));item['flux']=(np.array(item['flux'])/ab).tolist();item['error']=(np.array(item['error'])/ab).tolist()
    if variant=='approximate_SDSS_to_AB_offsets':
        ab=10**(-.4*np.array([-.04,0,0,0,.02]));item['flux']=(np.array(item['flux'])*ab).tolist();item['error']=(np.array(item['error'])*ab).tolist()
    rows=fit.fit_one(item)
    for row in rows:row.update(sample=sample,variant=variant)
    return rows

def main():
    start=time.perf_counter();jobs=[];selection={}
    for sample in ['sdss','des','roman']:
        items=sorted(fit.inputs(sample)[0],key=lambda x:hashlib.sha256(x['id'].encode()).hexdigest())[:8]
        selection[sample]=[x['id'] for x in items]
        variants=['base','refined_grid','surviving_mass_age','nonaccelerating_EdS_clock']
        if sample=='sdss':variants+=['SFD_scale_1','native_SDSS_zero_points']
        if sample=='roman':variants+=['approximate_SDSS_to_AB_offsets']
        for variant in variants:
            for item in items:jobs.append((sample,variant,item))
    rows=[]
    with ProcessPoolExecutor(4) as pool:
        for i,result in enumerate(pool.map(worker,jobs,chunksize=1)):
            rows+=result
            if (i+1)%16==0:print(f'{i+1}/{len(jobs)} sensitivity jobs',flush=True)
    table=pd.DataFrame(rows);path=WORK/'age-bound-sensitivity.csv';table.to_csv(path,index=False)
    summary=[]
    for (sample,variant,floor),group in table.groupby(['sample','variant','model_discrepancy_mag']):
        base=table[(table['sample']==sample)&(table.variant=='base')&(table.model_discrepancy_mag==floor)]
        joined=group.merge(base,on='id',suffixes=('_new','_base'),validate='one_to_one')
        summary.append(dict(sample=sample,variant=variant,model_discrepancy_mag=floor,n=len(group),
            feasible=int(group.age_low_Gyr.notna().sum()),joint_feasible=int(group.joint_age_low_Gyr.notna().sum()),numerical_failures=int(group.numerical_failures.sum()),
            age_width_Gyr=fit.quantiles(group.age_high_Gyr-group.age_low_Gyr),
            lower_bound_shift_Gyr=fit.quantiles(joined.age_low_Gyr_new-joined.age_low_Gyr_base),
            upper_bound_shift_Gyr=fit.quantiles(joined.age_high_Gyr_new-joined.age_high_Gyr_base)))
    dependencies=[Path(__file__).with_name(name) for name in ['model.py','fit.py','design.json']]+[WORK/'ssp-library.npz',WORK/'ssp-library-refined.npz']
    dependencies += [p for sample in ['sdss','des','roman'] for p in fit.inputs(sample)[1]]
    result=dict(selection='First8IDs in SHA256 order within each eligible sample, no fitted-result selection',ids=selection,
        dependencies_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(set(dependencies))},
        variants='Refined includes additional native FSPS ages/metallicities and midpoint dust optical depths. Surviving mass changes the estimand, not merely the label. EdS changes only the stellar-age clock and free-mass stellar fit, not the supernova distance likelihood.',
        seconds=time.perf_counter()-start,summaries=summary,output=str(path.relative_to(ROOT)),output_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (OUT/'sensitivity.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
