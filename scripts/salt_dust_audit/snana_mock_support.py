#!/usr/bin/env python3
"""Audit actual written simulation RV/AV support directly from released FITS.

Evaluates model attenuation at fixed rest wavelengths; does not infer observed
broadband attenuation or redo SALT fitting/BBC. No real-SN dust is measured here.
"""
from __future__ import annotations
import csv,hashlib,json
from pathlib import Path
import numpy as np
from astropy.io import fits
from scipy.optimize import brentq
from snana_extinction import ROOT,SnanaExtinction

OUT=ROOT/'runs/salt_dust_audit/snana_mock_support'
BASE=ROOT/'sources/repos/des-science__DES-SN5YR@1.3/1_SIMULATIONS/SNIa_SIMULATIONS'


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    ext=SnanaExtinction(ROOT/'runs/salt_dust_audit/snana_extinction/libsnana_extinction.so')
    paths=sorted(BASE.glob('*/*_HEAD.FITS.gz'))
    arrays=[];sources=[]
    for path in paths:
        with fits.open(path,memmap=False)as hdul:
            d=hdul[1].data
            arrays.append({'path':str(path.relative_to(ROOT)), 'SNID':np.char.strip(np.asarray(d['SNID']).astype(str)), 'RV':np.asarray(d['SIM_RV'],dtype=float), 'AV':np.asarray(d['SIM_AV'],dtype=float), 'z':np.asarray(d['SIM_REDSHIFT_CMB'],dtype=float), 'logmass':np.asarray(d['HOSTGAL_LOGMASS'],dtype=float),'N':len(d)})
        sources.append({'path':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'N':arrays[-1]['N']})
    samples=[];negrows=[];allrv=[];allav=[];allz=[];allmass=[];alla=[]
    threshold_old=brentq(lambda r:float(ext(8000,r,option=-99)),.1,2)
    threshold_new=brentq(lambda r:float(ext(8000,r,option=99)),.1,2)
    for arr in arrays:
        rv=arr['RV'];av=arr['AV'];assert np.all(rv>0) and np.all(av>=0)
        ebv=av/rv
        # Historical approximation is the released mocks' chosen 99 implementation.
        old=np.array([float(ext(8000,r,e,option=-99))for r,e in zip(rv,ebv)])
        new=np.array([float(ext(8000,r,e,option=99))for r,e in zip(rv,ebv)])
        assert np.array_equal(old<0,(rv<threshold_old)&(av>0))
        assert np.array_equal(new<0,(rv<threshold_new)&(av>0))
        samples.append({'realization':Path(arr['path']).parent.name,'N_written_simulated':len(rv),'RV_min':float(rv.min()),'N_RV_lt1':int((rv<1).sum()),'N_RV_lt2':int((rv<2).sum()),'N_A8000_old_lt0':int((old<0).sum()),'N_A8000_exact_lt0':int((new<0).sum()),'min_A8000_old_mag':float(old.min()),'min_A8000_exact_mag':float(new.min())})
        for i in np.where((old<0)|(new<0))[0]:
            negrows.append({'realization':Path(arr['path']).parent.name,'SNID':arr['SNID'][i],'SIM_RV':rv[i],'SIM_AV':av[i],'SIM_EBV':ebv[i],'SIM_zCMB':arr['z'][i],'HOSTGAL_LOGMASS':arr['logmass'][i],'A8000_old_mag':old[i],'A8000_exact_mag':new[i]})
        allrv.append(rv);allav.append(av);allz.append(arr['z']);allmass.append(arr['logmass']);alla.append(np.column_stack([old,new]))
    for name,rows in [('per_realization.csv',samples),('negative_A8000_simulated_rows.csv',negrows)]:
        with (OUT/name).open('w')as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    rv=np.concatenate(allrv);av=np.concatenate(allav);z=np.concatenate(allz);mass=np.concatenate(allmass);atten=np.concatenate(alla)
    summaries=[]
    for label,sel in [('all',np.ones(len(rv),bool)),('logmass_lt10',mass<10),('logmass_ge10',mass>=10),('z_lt0p2',z<.2),('z_ge0p2_lt0p5',(z>=.2)&(z<.5)),('z_ge0p5',z>=.5)]:
        a=atten[sel]
        summaries.append({'subset':label,'N':int(sel.sum()),'N_RV_lt1':int((rv[sel]<1).sum()),'N_RV_lt2':int((rv[sel]<2).sum()),'N_A8000_old_lt0':int((a[:,0]<0).sum()),'N_A8000_exact_lt0':int((a[:,1]<0).sum()),'mean_A8000_exact_minus_old_mag':float(np.mean(a[:,1]-a[:,0]))})
    result={'classification':'Observed truth metadata of released simulated (not real) objects, after the generation/write selection. Further light-curve-fitting, classification and BBC cuts are not asserted.', 'N_realizations':len(paths),'N_written_simulated':len(rv),'RV_min':float(rv.min()),'zero_extinction_RV_at_8000A':{'old':threshold_old,'exact':threshold_new}, 'summaries':summaries,'caution':'A8000 is the host extinction component at a fixed rest wavelength evaluated on the stored SIM_AV/SIM_RV, independent of filters. Negative A implies that this component of the generator multiplier exceeds unity. Net broadband flux changes, fitted colors, selection changes and cosmological biases require matched reruns and are not established by this table. These counts are not error bars or frequencies of real interstellar dust.', 'input_files':sources}
    (OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    extras=[Path(__file__).resolve(),ROOT/'scripts/salt_dust_audit/snana_extinction.py',ROOT/'sources/repos/RickKessler__SNANA/src/genmag_SEDtools.c',ROOT/'sources/repos/RickKessler__SNANA/src/genmag_SALT2.c',ROOT/'phase2/official/build/SNANA-2fe0f56/src/genmag_SEDtools.c',BASE/'PIP_D5YR_SIM_V2_DATADESSIM_4D_P21.input']+sorted(OUT.glob('*.csv'))+[OUT/'results.json']
    (OUT/'manifest.json').write_text(json.dumps({'files':sources+[{'path':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}for p in extras]},indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items()if k!='input_files'},indent=2))

if __name__=='__main__':main()
