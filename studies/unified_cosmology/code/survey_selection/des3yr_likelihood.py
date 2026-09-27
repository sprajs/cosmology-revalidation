#!/usr/bin/env python3
"""Exact released DES3YR 20-bin CosmoMC likelihood adapter, not 329-row covariance.

JLA source adds dmb^2 once to the supplied magnitude systematic covariance;
pecz=intrinsicdisp=0 in the released dataset. Model evaluated at author bin z.
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.linalg import cho_factor,cho_solve
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
from common import ROOT,WORK,RESULTS,sha
from acquire import get
COMMIT='83f00ac4742c7b48563204a22d7a7d3bb597a3e1'


def covariance(path):
    vals=np.loadtxt(path);n=int(vals[0]);assert vals.size==1+n*n
    c=vals[1:].reshape(n,n);assert np.max(abs(c-c.T))<1e-10
    return c


def main():
    b=WORK/'des3yr/05-COSMOLOGY/COSMOLOGY_INPUTS';out=WORK/'normalized';sources=[];records={}
    for name in ['supernovae_JLA.f90','supernovae.f90']:
        sources.append(get(f'https://raw.githubusercontent.com/cmbant/CosmoMC/{COMMIT}/source/{name}',WORK/'des3yr/cosmomc-source'/name))
    for label,events in [('DES+LOWz',329),('DESonly',207)]:
        name='combined' if label=='DES+LOWz' else 'des-only'
        data=b/f'lcparam_{label}.txt';raw=np.loadtxt(data);assert raw.shape[0]==20
        ids=raw[:,0].astype(int);z,zhel,mb,err=raw[:,1],raw[:,2],raw[:,4],raw[:,5]
        assert np.array_equal(ids,np.arange(20)) and np.all(err>0)
        nominal=np.diag(err**2);sys=covariance(b/f'sys_{label}_ALLSYS.txt');sys0=covariance(b/f'sys_{label}_STATONLY.txt');assert np.all(sys0==0)
        assembly={}
        for typ,extra in [('total',sys),('stat',sys0)]:
            cov=nominal+extra;fac=cho_factor(cov,lower=True);p=cho_solve(fac,np.eye(20));np.linalg.cholesky(cov)
            path=out/f'des3yr-{name}-{typ}.npz';np.savez_compressed(path,CID=np.array([f'DES3YR_{name}_BIN_{i:02d}' for i in ids]),zHD=z,zHEL=zhel,MU=mb+19.3,covariance=cov,precision=p)
            # Independent published flat reference geometry, evaluated at bin z,
            # one free global M. Arbitrary+19.3 convention cancels exactly.
            dl=(1+zhel)*299792.458/70*np.array([quad(lambda t:1/np.sqrt(.3*(1+t)**3+.7),0,float(q),epsabs=1e-12,epsrel=1e-12)[0] for q in z]);model=5*np.log10(dl)+25
            r=mb+19.3-model;one=np.ones(20);u=p@one;best=float(r@u/u.sum());chi=float(r@p@r-(r@u)**2/u.sum())
            brute=minimize_scalar(lambda m:(r-m)@p@(r-m),bracket=(-2,2),method='brent',options={'xtol':1e-12});assert abs(brute.fun-chi)<1e-7
            assembly[typ]={'NPZ':str(path.relative_to(ROOT)),'sha256':sha(path),'minimum_covariance_eigenvalue':float(np.linalg.eigvalsh(cov).min()),'inverse_identity_max':float(np.max(abs(cov@p-np.eye(20)))),'reference_Omega_m':.3,'reference_w':-1.,'reference_H0':70.,'reference_profiled_chi2':chi,'reference_profiled_global_offset':best,'direct_scalar_offset_minimization_delta_chi2':abs(float(brute.fun)-chi)}
        frame=pd.DataFrame({'likelihood_row':ids,'bin_ID':[f'DES3YR_{name}_BIN_{i:02d}' for i in ids],'z_author':z,'zHEL_author':zhel,'published_standardized_magnitude':mb,'MU_arbitrary_offset_19_3':mb+19.3,'statistical_error':err,'empty_bin_sentinel':err>=999})
        path=out/f'des3yr-{name}-bins.csv';frame.to_csv(path,index=False,float_format='%.17g')
        records[name]={'likelihood_rows':20,'physical_SNe_represented':events,'empty_bins_with_999mag_error':int((err>=999).sum()),'nonempty_bins':int((err<999).sum()),'assembly':assembly,'ledger':{'path':str(path.relative_to(ROOT)),'sha256':sha(path)},'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [data,b/f'sys_{label}_ALLSYS.txt',b/f'sys_{label}_STATONLY.txt',b/f'{label}_ALLSYS.dataset']}}
    result={'status':'passed_exact_released_binned_interface','code_sha256':sha(__file__),'author_CosmoMC_source_commit':COMMIT,'author_CosmoMC_source':sources,'samples':records,'equations':'C_total=diag(dmb^2)+C_mag_sys; othercovariancesdisabled,pecz=0,intrinsicdisp=0. dL=(1+zHEL)(1+zCMB)DA(zCMB) evaluated at eachreleasedbin redshift. Onegloballyfree magnitude marginalized/profiled.','author_source_lines':{'dmb_squared':476,'add_prevars_once':[813,828,913],'redshift_geometry':[1193,1198]},'limits':['20Gaussian likelihoodrows represent329physical SNe forcombinedsample, not20new events; two rows have999magemptybinsentinels andnegligibleweight.','No329x329eventsystematiccovariance wasrecovered; donotexpand20bincovariancetoindividualSNe.','The+19.3mag conversion is onlyaninterfaceconvention; itisnotanexternalabsolute-calibration measurement.','Modelis evaluatedat releasedbinz, not newlyaveragedoverunavailable individualwithin-bin responseweights.','Alternative toDovekie/Pantheon+, never independentjoint input because events/calibration overlap.','HistoricalSALT2,selection,scattermodelandBBCcompression assumptionsremainconditional.']}
    (RESULTS/'des3yr-interface.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(records,indent=2))
if __name__=='__main__':main()
