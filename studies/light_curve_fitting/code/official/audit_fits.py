#!/usr/bin/env python3
"""Audit freshly fitted observables; export exact Minuit Hessian covariances."""
from pathlib import Path
import argparse,hashlib,json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[3]; OUT=ROOT/'phase2/official'
RELEASE=ROOT/'sources/repos/des-science__DES-SN5YR@1.3'

def read_fit(path):
    return pd.read_csv(path,sep=r'\s+',comment='#',dtype={'CID':str}).drop(columns='VARNAMES:')

def hessian(path):
    rows=[]
    for line in path.read_text().splitlines():
        if line.startswith('PHASE2_HESSIAN:'):
            q=line.split(); vals=np.array([float(x) for x in q[2:]])
            assert len(vals)==20
            rows.append((q[1],vals[:4],vals[4:].reshape(4,4)))
    assert len(set(x[0] for x in rows))==len(rows)
    return {x[0]:(x[1],x[2]) for x in rows}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--sample',choices=['pilot','all','recovered','conditioned','conditioned_multistart_best'],default='pilot');args=ap.parse_args()
    labels=['pilot12_hessian'] if args.sample=='pilot' else [f'des_hessian{i}' for i in range(4)]+['lowz_hessian','foundation_hessian']
    if args.sample=='recovered':labels=[f'des_mask32{i}' for i in range(4)]+['lowz_consistent','foundation_hessian']
    if args.sample=='conditioned':labels=[f'des_conditioned{i}' for i in range(4)]+['lowz_consistent','foundation_hessian']
    if args.sample=='conditioned_multistart_best':
        labels=[f'des_conditioned{i}' for i in range(4)]+['lowz_consistent','foundation_hessian']
        manifest=json.loads((OUT/'inputs/conditioned-multistart-contract.json').read_text())
        labels += ['conditioned_seed_'+t['label'] for t in manifest['tasks']]
    frames=[];pars=[];covs=[];logs=[]
    for label in labels:
        log=OUT/'diagnostics'/('snana-'+label.replace('_','-')+'.log')
        if label.startswith('conditioned_seed_'):
            log=OUT/'diagnostics'/('snana-conditioned-seed-'+label.removeprefix('conditioned_seed_')+'.log')
        # DES labels contain numeric suffix without hyphen, consistently with runner.
        assert 'ENDING PROGRAM GRACEFULLY.' in log.read_text(),f'Incomplete run: {log}'
        frame=read_fit(OUT/f'results/snana_{label}.FITRES.TEXT');raw=hessian(log)
        assert set(frame.CID)==set(raw), (label,len(frame),len(raw))
        for cid in frame.CID:
            p,c=raw[cid];pars.append(p);covs.append(c)
        frame['source_fit_label']=label
        frames.append(frame);logs.append(log)
    fits=pd.concat(frames,ignore_index=True);params=np.array(pars);cov=np.array(covs)
    if args.sample=='conditioned_multistart_best':
        selected=fits.groupby(['CID','IDSURVEY'],sort=False).FITCHI2.idxmin().to_numpy()
        fits=fits.loc[selected].reset_index(drop=True);params=params[selected];cov=cov[selected]
    assert not fits.duplicated(['CID','IDSURVEY']).any()
    # Exact covariance changes basis from x0 to mx=-2.5log10(x0); additive mB
    # convention is 10.635 and therefore has identical covariance.
    jac=np.ones((len(fits),4));jac[:,0]=-2.5/np.log(10)/params[:,0]
    mcov=cov*jac[:,:,None]*jac[:,None,:]
    mpar=params.copy();mpar[:,0]=-2.5*np.log10(params[:,0])
    corr=mcov/np.sqrt(np.diagonal(mcov,axis1=1,axis2=2))[:,:,None]/np.sqrt(np.diagonal(mcov,axis1=1,axis2=2))[:,None,:]
    eig=np.linalg.eigvalsh(corr)
    fits['x0_double']=params[:,0];fits['x1_double']=params[:,1];fits['c_double']=params[:,2];fits['t0_double']=params[:,3]
    fits['mx_double']=mpar[:,0];fits['mB_double']=mpar[:,0]+10.635
    fits['hessian_correlation_min_eigenvalue']=eig[:,0]
    for i,k in enumerate(['mx','x1','c','t0']):fits[f'{k}_hessian_sigma']=np.sqrt(mcov[:,i,i])
    fits.to_csv(OUT/f'results/{args.sample}_refit_observables.csv',index=False,float_format='%.17g')
    np.savez_compressed(OUT/f'results/{args.sample}_refit_covariance.npz',CID=fits.CID.astype(str).to_numpy(dtype='U32'),IDSURVEY=fits.IDSURVEY.to_numpy(),parameter_order=np.array(['x0','x1','c','t0']),parameters=params,covariance=cov,mag_parameter_order=np.array(['mx','x1','c','t0']),mag_parameters=mpar,mag_covariance=mcov)
    ref=pd.read_csv(RELEASE/'4_DISTANCES_COVMAT/DES-SN5YR_HD+MetaData.csv',dtype={'CID':str})
    comp=fits.merge(ref,on=['CID','IDSURVEY'],suffixes=('_refit','_published'),validate='one_to_one')
    diag={'sample':args.sample,'n_fits':len(fits),'n_matched_reference':len(comp),'n_published':len(ref),'hessian_nonpositive_correlation_eigenvalues':int(np.sum(eig[:,0]<=0)), 'hessian_min_correlation_eigenvalue':float(np.min(eig)), 'hessian_max_asymmetry':float(np.max(np.abs(cov-cov.swapaxes(1,2)))),'fit_error_flags':{str(k):int(v) for k,v in fits.ERRFLAG_FIT.value_counts().items()},'differences':{},'surveys':{},'missing_reference_keys':[list(x) for x in set(zip(ref.CID,ref.IDSURVEY))-set(zip(fits.CID,fits.IDSURVEY))], 'covariance_definition':'Minuit full Hessian FITERRMAT; 4x4 order x0,x1,c,t0; transformed mx,x1,c,t0. MINOS scalar errors retained separately in CSV; not used as Hessian diagonals. No PSD repair.'}
    for key in ['mB','x1','c','PKMJD','NDOF','FITPROB']:
        delta=comp[key+'_refit']-comp[key+'_published'];comp['delta_'+key]=delta
        diag['differences'][key]={'median':float(np.median(delta)),'rms':float(np.sqrt(np.mean(delta**2))),'q90_abs':float(np.quantile(np.abs(delta),.9)),'max_abs':float(np.max(np.abs(delta)))}
    comp['delta_Tripp_fixed_coefficients']=comp.delta_mB+.16087*comp.delta_x1-3.1178*comp.delta_c
    for survey,g in comp.groupby('IDSURVEY'):
        diag['surveys'][str(survey)]={'n':len(g),'mB_median_shift':float(g.delta_mB.median()),'Tripp_median_shift':float(g.delta_Tripp_fixed_coefficients.median()),'n_different_NDOF':int((g.delta_NDOF!=0).sum())}
    comp.to_csv(OUT/f'results/{args.sample}_comparison.csv',index=False,float_format='%.17g')
    diag['source_hashes']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in logs+[OUT/'build/SNANA-current/bin/snlc_fit.exe']}
    if args.sample in ['conditioned','conditioned_multistart_best']:
        diag['arm_definition']='DES author-accepted-mask conditioning with no new clipping/phase cuts; lowz and Foundation retain recovered baseline fits. Not an independent recovery of author masking.'
    if args.sample=='conditioned_multistart_best':
        diag['branch_selection']='Lowest reported SNANA FITCHI2 among default plus3x3starting grid and reference diagnostic for10prespecified DES anomalies. No choice based on closeness to published values; no omitted objects. Finite-start best, not proof of global minimum; iterative covariance differs across branch solutions.1307748 branch ranking also confirmed under both common frozen covariances.'
    (OUT/f'diagnostics/{args.sample}_refit_summary.json').write_text(json.dumps(diag,indent=2)+'\n');print(json.dumps(diag,indent=2))

if __name__=='__main__':main()
