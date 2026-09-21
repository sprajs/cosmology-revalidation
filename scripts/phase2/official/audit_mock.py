#!/usr/bin/env python3
"""Trace first nominal Ia mock through actual light-curve fitting and quality cuts."""
from pathlib import Path
import argparse,json,hashlib,re
import numpy as np
import pandas as pd
from astropy.io import fits

ROOT=Path(__file__).resolve().parents[3];OUT=ROOT/'phase2/official'
VERSION='PIP_D5YR_SIM_V2_DATADESSIM_4D_P21-0001'
SOURCE=ROOT/'sources/repos/des-science__DES-SN5YR@1.3/1_SIMULATIONS/SNIa_SIMULATIONS'/VERSION
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--label',default='mock0001');ap.add_argument('--source',type=Path,default=SOURCE);args=ap.parse_args()
    hp=next(args.source.glob('*HEAD.FITS*'));h=fits.getdata(hp,1)
    log=OUT/f'diagnostics/snana-{args.label.replace("_","-")}.log';logtext=log.read_text();assert 'ENDING PROGRAM GRACEFULLY.' in logtext
    a=pd.read_csv(OUT/f'results/snana_{args.label}.FITRES.TEXT',sep=r'\s+',comment='#',dtype={'CID':str}).drop(columns='VARNAMES:')
    a['BBC_quality_only']=(a.x1.abs()<3)&(a.c.abs()<.3)&(a.x1ERR<1)&(a.PKMJDERR<2)&(a.cERR<1.5)&(a.FITPROB>.001)&(a.zHD>.025)&(a.zHD<1.2)
    from audit_fits import hessian
    if 'PHASE2_HESSIAN:' in logtext:
        raw=hessian(log);assert set(raw)==set(a.CID)
        params=np.array([raw[c][0] for c in a.CID]);cov=np.array([raw[c][1] for c in a.CID])
        jac=np.ones((len(a),4));jac[:,0]=-2.5/np.log(10)/params[:,0];mcov=cov*jac[:,:,None]*jac[:,None,:]
        mpar=params.copy();mpar[:,0]=-2.5*np.log10(params[:,0])
        np.savez_compressed(OUT/f'results/{args.label}_refit_covariance.npz',CID=a.CID.to_numpy(dtype='U32'),IDSURVEY=a.IDSURVEY.to_numpy(),parameters=params,covariance=cov,parameter_order=np.array(['x0','x1','c','t0']),mag_parameters=mpar,mag_covariance=mcov,mag_parameter_order=np.array(['mx','x1','c','t0']))
        for i,k in enumerate(['x0','x1','c','t0']):a[k+'_double']=params[:,i]
        a['mx_double']=mpar[:,0];a['mB_double']=mpar[:,0]+10.635
    a.to_csv(OUT/f'results/{args.label}_fit_and_quality.csv.gz',index=False,float_format='%.17g')
    status=pd.DataFrame({'CID':[str(x).strip() for x in h['SNID']]})
    status['lightcurve_fit_pass']=status.CID.isin(a.CID)
    status['BBC_quality_only']=status.CID.isin(a.loc[a.BBC_quality_only,'CID'])
    status.to_csv(OUT/f'results/{args.label}_all_detected_status.csv',index=False)
    npre=int(re.search(r'Finished processing\s+(\d+) SN after snana\s+cuts',logtext)[1])
    result={'version':args.source.name,'detected_head_n':len(h),'SNANA_preselection_pass_n':npre,'lightcurve_fit_pass_n':len(a),'BBC_quality_only_n':int(a.BBC_quality_only.sum()),'not_complete_BBC_selection':['bias-correction map coverage','BBC chi2max=16 rejection','consistent sample across systematic variants','classifier weighting and contaminant-prior likelihood not reapplied'], 'heads_sha256':hashlib.file_digest(hp.open('rb'),'sha256').hexdigest(),'fitter':'SNANA 886408a4 with output-only patch; official original DES namelist; no scientific source modification','full_hessian_exported':'PHASE2_HESSIAN:' in logtext}
    (OUT/f'diagnostics/{args.label}_summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
