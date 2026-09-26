#!/usr/bin/env python3
"""Bounded actual-SNANA fit of a noise-preserving mock dust-floor injection.

This is a deliberately selected stress sample; no survey/host-redshift selection
or BBC regeneration. It cannot supply a population correction or cosmology bias.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,os,re,subprocess,time
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.io import fits
import sncosmo
from flux_response import ROOT,build_model
from snana_extinction import SnanaExtinction

OUT=ROOT/'runs/salt_dust_audit/snana_stress'
SOURCE=ROOT/'sources/repos/des-science__DES-SN5YR@1.3/1_SIMULATIONS/SNIa_SIMULATIONS/PIP_D5YR_SIM_V2_DATADESSIM_4D_P21-0001'
REF=ROOT/'phase2/official/results/snana_mock0001.FITRES.TEXT'
BINARY=ROOT/'phase2/official/build/SNANA-current/bin/snlc_fit.exe'


def read_fit(path):
    lines=path.read_text().splitlines();names=next(l.split()[1:]for l in lines if l.startswith('VARNAMES:'))
    rows=[l.split()[1:]for l in lines if l.startswith('SN:')]
    frame=pd.DataFrame(rows,columns=names);frame.CID=frame.CID.astype(str)
    for col in frame.columns:
        if col!='CID':
            try:frame[col]=pd.to_numeric(frame[col])
            except(ValueError,TypeError):pass
    return frame


class NativeDust(sncosmo.PropagationEffect):
    _param_names=['ebv','rv'];_minwave=1000.;_maxwave=25000.
    def __init__(self,ext,floor=False):
        self._parameters=np.array([0.,3.1]);self.ext=ext;self.floor=floor;self.cache={}
    def propagate(self,wave,flux,phase=None):
        ebv,rv=self._parameters
        key=(float(ebv),float(rv),wave.tobytes())
        if key not in self.cache:
            a=self.ext(wave,float(rv),float(ebv),option=-99)
            if self.floor:a=np.maximum(a,0)
            self.cache[key]=10**(-.4*a)
        return flux*self.cache[key]


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def prepare():
    OUT.mkdir(parents=True,exist_ok=True)
    hp=next(SOURCE.glob('*_HEAD.FITS.gz'));pp=next(SOURCE.glob('*_PHOT.FITS.gz'))
    ref=read_fit(REF)
    # Select on the already stored baseline for stable paired fitter diagnostics.
    good=ref[(ref.ERRFLAG_FIT==0)&(ref.x1.abs()<3)&(ref.c.abs()<.3)&(ref.x1ERR<1)&(ref.PKMJDERR<2)&(ref.FITPROB>.001)].copy()
    low=good[(good.SIM_RV<.65)&(good.SIM_AV>0)].sort_values(['zHEL','CID'])
    assert len(low)>=12
    # Twelve evenly spaced redshift order statistics, fixed before stress outcomes.
    chosen=low.iloc[np.rint(np.linspace(0,len(low)-1,12)).astype(int)].copy();chosen['role']='low_RV_stress'
    controls=[];pool=good[(good.SIM_RV>=2)&(good.SIM_RV<=4)].copy()
    for _,row in chosen.iterrows():
        i=(pool.zHEL-row.zHEL).abs().idxmin();v=pool.loc[i].copy();v['role']='comparison_control';v['matched_CID']=row.CID
        controls.append(v);pool=pool.drop(i)
    selected=pd.concat([chosen,pd.DataFrame(controls)],ignore_index=True)
    selected[['CID','role','matched_CID','zHEL','SIM_RV','SIM_AV','SIM_x1','SIM_c','SIM_PEAKMJD','mB','x1','c','PKMJD','FITPROB','NDOF']].to_csv(OUT/'selection.csv',index=False)
    (OUT/'selected.cid').write_text('\n'.join(selected.CID)+'\n')
    ext=SnanaExtinction(ROOT/'runs/salt_dust_audit/snana_extinction/libsnana_extinction.so')
    ordinary,bands,modelpaths,_=build_model()
    # Both models use the historical SNANA extinction law for generator truth.
    baseline=sncosmo.Model(source=ordinary.source,effects=[NativeDust(ext),NativeDust(ext)],effect_names=['mw','host'],effect_frames=['obs','rest'])
    floored=sncosmo.Model(source=ordinary.source,effects=[NativeDust(ext),NativeDust(ext,True)],effect_names=['mw','host'],effect_frames=['obs','rest'])
    ratios=[];changed=[]
    with fits.open(hp,memmap=False)as hh, fits.open(pp,memmap=False)as ph:
        heads=hh[1].data;phot=ph[1].data
        ids=np.char.strip(heads['SNID'].astype(str))
        modified=phot.copy()
        for _,sel in selected.iterrows():
            h=heads[np.where(ids==sel.CID)[0][0]]
            start=int(h['PTROBS_MIN'])-1;stop=int(h['PTROBS_MAX']);part=phot[start:stop]
            params=dict(z=float(h['SIM_REDSHIFT_HELIO']),t0=float(h['SIM_PEAKMJD']),x0=float(h['SIM_SALT2x0']),x1=float(h['SIM_SALT2x1']),c=float(h['SIM_SALT2c']),mwebv=float(h['SIM_MWEBV']),mwrv=3.1,hostebv=float(h['SIM_AV']/h['SIM_RV']),hostrv=float(h['SIM_RV']))
            baseline.set(**params);floored.set(**params)
            nchange=0;nskip=0;max_delta=0.;max_signal_shift=0.;max_noise_error=0.
            for j,pt in enumerate(part):
                t=float(pt['MJD']);band=str(pt['BAND']).strip();mag=float(pt['SIM_MAGOBS'])
                if band not in bands or not baseline.bandoverlap(bands[band]) or not baseline.mintime()<=t<=baseline.maxtime() or not -50<mag<50:
                    nskip+=1;continue
                before=float(baseline.bandflux(bands[band],t));after=float(floored.bandflux(bands[band],t))
                if before<=0 or after<=0:raise ValueError(f'Nonpositive model flux {sel.CID} {j}')
                ratio=after/before;trueflux=10**(.4*(27.5-mag));delta=trueflux*(ratio-1)
                # Underlying float32 serialization necessarily rounds the preserved residual.
                value=np.float32(float(pt['FLUXCAL'])+delta)
                modified['FLUXCAL'][start+j]=value
                modified['SIM_MAGOBS'][start+j]=np.float32(mag-2.5*np.log10(ratio))
                nchange+=abs(delta)>1e-10;max_delta=max(max_delta,abs(delta))
                if pt['FLUXCALERR']>0:max_signal_shift=max(max_signal_shift,abs(delta/pt['FLUXCALERR']))
                residual_error=float(value)-(float(pt['FLUXCAL'])+delta);max_noise_error=max(max_noise_error,abs(residual_error))
                ratios.append(dict(CID=sel.CID,role=sel.role,source_phot_row=start+j,band=band,MJD=t,phase=(t-params['t0'])/(1+params['z']),SIM_RV=params['hostrv'],SIM_EBV=params['hostebv'],flux_ratio_floor_over_old=ratio,delta_fluxcal=delta,fluxcalerr=float(pt['FLUXCALERR']),rounding_error=residual_error))
            changed.append(dict(CID=sel.CID,role=sel.role,N_observations=len(part),N_changed=nchange,N_not_modified_outside_model_support=nskip,max_abs_delta_fluxcal=max_delta,max_abs_delta_over_sigma=max_signal_shift,max_abs_noise_residual_rounding=max_noise_error))
        assert np.array_equal(modified['FLUXCALERR'],phot['FLUXCALERR'])
        assert np.array_equal(modified['MJD'],phot['MJD'])
        for variant in ['baseline','noop','floor']:
            version='SALT_DUST_STRESS_'+variant.upper();dest=OUT/version;dest.mkdir(exist_ok=True)
            # Each variant owns complete files; source row pointers remain unchanged.
            hcopy=fits.HDUList([h.copy()for h in hh]);pcopy=fits.HDUList([h.copy()for h in ph])
            hname=version+'_HEAD.FITS';pname=version+'_PHOT.FITS'
            hcopy[0].header['VERSION']=version;hcopy[0].header['PHOTFILE']=pname
            pcopy[0].header['VERSION']=version
            if variant=='floor':pcopy[1].data=modified
            hcopy.writeto(dest/hname,overwrite=True);pcopy.writeto(dest/pname,overwrite=True)
            (dest/(version+'.LIST')).write_text(hname+'\n')
            (dest/(version+'.README')).write_text('Local fixed-selected-sample dust stress fixture; see ../contract.json.\n')
            text=(ROOT/'phase2/official/inputs/snana_mock0001.nml').read_text()
            text=re.sub(r"PRIVATE_DATA_PATH\s*=\s*'[^']*'",f"PRIVATE_DATA_PATH = '{OUT}'",text)
            text=re.sub(r"VERSION_PHOTOMETRY\s*=\s*'[^']*'",f"VERSION_PHOTOMETRY = '{version}'",text)
            text=re.sub(r"TEXTFILE_PREFIX\s*=\s*'[^']*'",f"TEXTFILE_PREFIX = '{OUT/variant}'",text)
            text=text.replace("'FITRES(text:host)'","'FITRES(text:host) LCPLOT(text:col)'")
            text=text.replace('&SNLCINP',f"&SNLCINP\n MXLC_PLOT = 100\n SNCID_LIST_FILE = '{OUT/'selected.cid'}'",1)
            (OUT/(variant+'.nml')).write_text(text)
        # No-op photometry is exactly identical to baseline, all columns/records.
        with fits.open(OUT/'SALT_DUST_STRESS_BASELINE/SALT_DUST_STRESS_BASELINE_PHOT.FITS')as a,fits.open(OUT/'SALT_DUST_STRESS_NOOP/SALT_DUST_STRESS_NOOP_PHOT.FITS')as b:
            assert a[1].data.tobytes()==b[1].data.tobytes()
    pd.DataFrame(ratios).to_csv(OUT/'epoch_injections.csv',index=False);pd.DataFrame(changed).to_csv(OUT/'object_injections.csv',index=False)
    contract={'status':'Prepared bounded counterfactual stress, not a dust correction or end-to-end simulation.', 'N_selected':len(selected),'selection':'12 baseline-quality-passing RV<.65 objects, evenly spaced in sorted redshift in original mock realization0001; 12 unique baseline-quality comparison objects with RV in[2,4], matched by nearest zHEL. Purposeful sample, not representative.', 'intervention':'Only host A_lambda<0 is floored to0; all remaining historic dust curve values unchanged. MW uses historical curve in both ratio calculations. This is an ad hoc nonnegative stress, not a proposed physical law.', 'source_phot_row_convention':'zero-based offset in original PHOT FITS table; not one-based SNANA PTROBS', 'preserved':'Original epoch times, observed bandpasses, headers, redshifts, errors, flags, and noise residuals (up to float32 rounding); original pre-fitting selection fixed.', 'formula':'F_injected=F_original+10**(.4*(27.5-SIM_MAGOBS))*(F_model_floor/F_model_old-1). Band ratio uses SNcosmo SALT3 with compiled SNANA historical dust routine, at stored generator truth; actual snlc_fit.exe fits both.', 'limits':'No rerun of detection/host-redshift/classifier selection, no population retuning, no Poisson-error recomputation, no BBC correction. Epochs outside SNcosmo source support or with invalid SIM_MAGOBS are unchanged and counted. Ratio uses released SALT3 deterministic SED; minute C11 scatter or late-time extrapolation is not reimplemented. Actual fitter cuts may move.', 'noop_photometry_identical':True,'inputs':[{'path':str(p.relative_to(ROOT)),'sha256':sha(p)}for p in [hp,pp,REF,BINARY,Path(__file__).resolve(),ROOT/'scripts/salt_dust_audit/flux_response.py',ROOT/'scripts/salt_dust_audit/snana_extinction.py']+modelpaths]}
    (OUT/'contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    print(json.dumps({'prepared':len(selected),'objects':changed},indent=2,default=lambda x:x.item()))


def run():
    env=os.environ.copy();env['SNANA_DIR']=str(BINARY.parents[1]);env['SNDATA_ROOT']=str(ROOT/'phase2/official/inputs/SNDATA_ROOT');env['LD_LIBRARY_PATH']=str(ROOT/'phase2/official/build/sysroot/usr/lib')+(':'+env['LD_LIBRARY_PATH']if env.get('LD_LIBRARY_PATH')else '')
    for variant in ['baseline','noop','floor']:
        stamp={'variant':variant,'command':[str(BINARY),str(OUT/(variant+'.nml'))],'started_unix':time.time()}
        with (OUT/(variant+'.log')).open('w')as f:
            process=subprocess.Popen(stamp['command'],env=env,stdout=f,stderr=subprocess.STDOUT,cwd=OUT)
            stamp['pid']=process.pid;(OUT/'live_run.json').write_text(json.dumps(stamp,indent=2)+'\n');print(json.dumps(stamp),flush=True)
            rc=process.wait()
        stamp.update(exit_code=rc,ended_unix=time.time());(OUT/(variant+'-run.json')).write_text(json.dumps(stamp,indent=2)+'\n')
        if rc!=0 or not (OUT/(variant+'.FITRES.TEXT')).exists():raise RuntimeError(f'{variant} fit failed: see log')
    analyze()


def audit_epochs(selected,baseline_fit,floor_fit):
    hp=next(SOURCE.glob('*_HEAD.FITS.gz'));pp=next(SOURCE.glob('*_PHOT.FITS.gz'))
    heads=fits.getdata(hp,1);phot=fits.getdata(pp,1);ids=np.char.strip(heads['SNID'].astype(str))
    injection=pd.read_csv(OUT/'epoch_injections.csv',dtype={'CID':str}).set_index('source_phot_row')
    flags={};max_match=0.;ambiguous=0
    for variant in ['baseline','floor']:
        lclist=(OUT/(variant+'.LCLIST.TEXT')).read_text()
        cols=re.findall(r'^# col\s+\d+ -> (\w+)\s+:',lclist,re.M)
        plot=pd.read_csv(OUT/(variant+'.LCPLOT.TEXT'),sep=r'\s+',header=None,comment='#',dtype={0:str})
        plot.columns=cols[:plot.shape[1]];plot=plot[plot.DATAFLAG!=0]
        mapping={}
        variant_flux=fits.getdata(OUT/('SALT_DUST_STRESS_'+variant.upper())/('SALT_DUST_STRESS_'+variant.upper()+'_PHOT.FITS'),1)['FLUXCAL']
        for cid,group in plot.groupby('CID'):
            h=heads[np.where(ids==cid)[0][0]];start=int(h['PTROBS_MIN'])-1;stop=int(h['PTROBS_MAX']);part=phot[start:stop]
            used=set()
            for _,point in group.iterrows():
                candidates=np.flatnonzero(np.char.strip(part['BAND'].astype(str))==point.BAND)
                distances=np.abs(np.asarray(part['MJD'][candidates])-point.MJD)
                valid=candidates[distances<.0021]
                ambiguous+=int(len(valid)>1)
                valid=np.array([x for x in valid if x not in used],dtype=int)
                assert len(valid)>0,(cid,point.MJD)
                # Repeated same-time/filter exposures are distinguished by measured flux.
                nearest=int(valid[np.argmin(np.abs(variant_flux[start+valid]-point.FLUXCAL))])
                distance=float(abs(part['MJD'][nearest]-point.MJD));max_match=max(max_match,distance)
                assert abs(variant_flux[start+nearest]-point.FLUXCAL)<max(.0001,abs(point.FLUXCAL)*5e-5),(cid,nearest,point.FLUXCAL)
                used.add(nearest);mapping[start+nearest]=int(point.DATAFLAG)
        flags[variant]=mapping
    rows=[]
    means={'g':4827.7,'r':6434.8,'i':7828.1,'z':9181.2} # actual SNANA filtdump, 0.1A precision
    for _,s in selected.iterrows():
        cid=s.CID;h=heads[np.where(ids==cid)[0][0]];start=int(h['PTROBS_MIN'])-1;stop=int(h['PTROBS_MAX'])
        for idx in range(start,stop):
            pt=phot[idx];band=str(pt['BAND']).strip()
            inj=injection.loc[idx]if idx in injection.index else None
            delta=float(inj.delta_fluxcal)if inj is not None else 0.;ratio=float(inj.flux_ratio_floor_over_old)if inj is not None else 1.
            rows.append(dict(CID=cid,role=s.role,source_phot_row_zero_based=idx,band=band,MJD=float(pt['MJD']),baseline_phase=(float(pt['MJD'])-baseline_fit.loc[cid,'PKMJD'])/(1+baseline_fit.loc[cid,'zHEL']),band_mean_rest_wavelength_A=means[band]/(1+baseline_fit.loc[cid,'zHEL']),baseline_dataflag=flags['baseline'].get(idx),floor_dataflag=flags['floor'].get(idx),injection_model_supported=inj is not None,changed=abs(delta)>1e-10,delta_fluxcal=delta,delta_broadband_mag=-2.5*np.log10(ratio),delta_over_sigma=delta/float(pt['FLUXCALERR'])if pt['FLUXCALERR']>0 else np.nan))
    frame=pd.DataFrame(rows);frame.to_csv(OUT/'epoch_acceptance.csv',index=False)
    summaries=[]
    for cid,q in frame.groupby('CID'):
        changed=q.changed;accepted=q.baseline_dataflag==1
        summaries.append(dict(CID=cid,N_all_epochs=len(q),N_changed_all=int(changed.sum()),N_accepted_baseline=int(accepted.sum()),N_accepted_floor=int((q.floor_dataflag==1).sum()),N_changed_accepted_baseline=int((changed&accepted).sum()),N_changed_rejected_baseline=int((changed&(q.baseline_dataflag==-1)).sum()),N_changed_absent_from_baseline_plot=int((changed&q.baseline_dataflag.isna()).sum()),N_accepted_but_injection_unsupported=int((accepted&~q.injection_model_supported).sum()),N_acceptance_changes=int(((q.baseline_dataflag==1)!=(q.floor_dataflag==1)).sum()),max_abs_injected_mag_accepted=float(q.loc[accepted,'delta_broadband_mag'].abs().max()),max_abs_delta_over_sigma_accepted=float(q.loc[accepted,'delta_over_sigma'].abs().max())))
    summary=pd.DataFrame(summaries);summary.to_csv(OUT/'object_epoch_acceptance.csv',index=False)
    return {'max_abs_MJD_matching_difference_day':max_match,'same_time_candidates_resolved_by_flux_and_uniqueness':ambiguous,'dataflag_meaning':'From actual LCPLOT:1 accepted,-1 rejected; missing means not present in fitter plot data. Source SNLC SNLCPLOT packs NFITDATA+NFITDATA_REJECT and IEP_REJECT. All source row references in this file zero-based.','N_accepted_but_injection_unsupported':int(summary.N_accepted_but_injection_unsupported.sum()),'N_acceptance_changes':int(summary.N_acceptance_changes.sum()),'N_objects_changed_in_accepted_epochs':int((summary.N_changed_accepted_baseline>0).sum())}


def analyze():
    selected=pd.read_csv(OUT/'selection.csv',dtype={'CID':str});a=read_fit(OUT/'baseline.FITRES.TEXT');b=read_fit(OUT/'noop.FITRES.TEXT');f=read_fit(OUT/'floor.FITRES.TEXT');ref=read_fit(REF)
    a=a.set_index('CID');b=b.set_index('CID');f=f.set_index('CID');ref=ref.set_index('CID')
    columns=['mB','x1','c','PKMJD','x1ERR','cERR','mBERR','FITCHI2','NDOF','FITPROB','ERRFLAG_FIT','CUTFLAG_SNANA']
    assert set(a.index)==set(b.index)
    residual=float(np.max(np.abs(a[columns].to_numpy(float)-b.loc[a.index,columns].to_numpy(float))))
    assert residual==0,residual
    baseline_reference=float(np.max(np.abs(a[columns].to_numpy(float)-ref.loc[a.index,columns].to_numpy(float))))
    assert baseline_reference<1e-6,baseline_reference
    def passes(r):return bool(r.ERRFLAG_FIT==0 and abs(r.x1)<3 and abs(r.c)<.3 and r.x1ERR<1 and r.PKMJDERR<2 and r.FITPROB>.001)
    rows=[]
    for _,s in selected.iterrows():
        cid=s.CID;r=dict(CID=cid,role=s.role,zHEL=s.zHEL,SIM_RV=s.SIM_RV,SIM_AV=s.SIM_AV,baseline_written=cid in a.index,floor_written=cid in f.index)
        if cid in a.index:r['baseline_quality_proxy']=passes(a.loc[cid])
        if cid in f.index:r['floor_quality_proxy']=passes(f.loc[cid])
        if cid in a.index and cid in f.index:
            for c in columns:r['delta_'+c]=float(f.loc[cid,c]-a.loc[cid,c])
            # Fixed reported original-DES nuisance values; no fit of nuisance parameters or BBC.
            r['delta_mu_fixed_alpha0p16087_beta3p11780_noBBC']=r['delta_mB']+.16087*r['delta_x1']-3.11780*r['delta_c']
        rows.append(r)
    pd.DataFrame(rows).to_csv(OUT/'fit_response.csv',index=False)
    epoch_audit=audit_epochs(selected,a,f)
    result={'epoch_audit':epoch_audit,'standardization_coefficients':{'alpha':.16087,'beta':3.11780,'source':'original DES tag1.3 4_DISTANCES_COVMAT/README.md; held fixed, no BBC'},'classification':'Actual SNANA nonlinear paired refit of selected mock counterfactuals, using hybrid noise-preserving flux injection. No final BBC/cosmology bias estimated.', 'N_requested':len(rows),'N_baseline_written':len(a),'N_floor_written':len(f),'baseline_vs_noop_max_abs_numeric_difference':residual,'baseline_vs_existing_reference_max_abs_numeric_difference':baseline_reference,'quality_proxy':'ERRFLAG_FIT=0, |x1|<3, |c|<.3, x1ERR<1, PKMJDERR<2, FITPROB>.001; not complete BBC/classifier selection.','summaries':[]}
    df=pd.DataFrame(rows)
    for role,q in df.groupby('role'):
        v=q.get('delta_mu_fixed_alpha0p16087_beta3p11780_noBBC',pd.Series(dtype=float)).dropna()
        result['summaries'].append({'role':role,'N':len(q),'N_both_fits':len(v),'N_baseline_quality_proxy':int(q.baseline_quality_proxy.fillna(False).sum()),'N_floor_quality_proxy':int(q.floor_quality_proxy.fillna(False).sum()),'fixed_standardization_delta_mu_median':float(v.median())if len(v)else None,'fixed_standardization_delta_mu_min':float(v.min())if len(v)else None,'fixed_standardization_delta_mu_max':float(v.max())if len(v)else None,'max_abs_delta_c':float(q.delta_c.abs().max()),'max_abs_delta_mB':float(q.delta_mB.abs().max())})
    (OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    files=sorted(p for p in OUT.rglob('*')if p.is_file() and p.name not in ['manifest.json','live_run.json'])
    files+=[Path(__file__).resolve(),ROOT/'sources/repos/RickKessler__SNANA/src/snlc_fit.F90',ROOT/'sources/repos/des-science__DES-SN5YR@1.3/4_DISTANCES_COVMAT/README.md']
    (OUT/'manifest.json').write_text(json.dumps({'files':[{'path':str(p.relative_to(ROOT)),'sha256':sha(p),'bytes':p.stat().st_size}for p in files]},indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run','analyze']);args=p.parse_args();globals()[args.action]()
