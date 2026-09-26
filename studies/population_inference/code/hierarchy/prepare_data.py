"""Join audited reference membership to coherent pre-BBC refit observables.

Preparation only: this neither validates the survey selection nor fits a model.
All outputs identify their particular refit arm and preserve reference values.
"""
import argparse,datetime,hashlib,json
from pathlib import Path
import numpy as np,pandas as pd
from scipy.special import ndtr

ROOT=Path(__file__).resolve().parents[3]
def main():
    p=argparse.ArgumentParser();p.add_argument('--refit-prefix',default='all_refit');p.add_argument('--name',default='literal-baseline');a=p.parse_args()
    indir=ROOT/'phase2/official/results';refitpath=indir/(a.refit_prefix+'_observables.csv');covpath=indir/(a.refit_prefix+'_covariance.npz')
    metapath=ROOT/'phase2/data_audit/original_metadata.csv.gz';foldpath=ROOT/'phase2/hierarchy/folds.csv'
    refit=pd.read_csv(refitpath,dtype={'CID':str});meta=pd.read_csv(metapath,dtype={'CID':str});folds=pd.read_csv(foldpath,dtype={'CID':str});cv=np.load(covpath,allow_pickle=False)
    keys=['IDSURVEY','CID'];assert not meta.duplicated(keys).any() and not refit.duplicated(keys).any()
    cov_index=pd.DataFrame({'IDSURVEY':cv['IDSURVEY'].astype(int),'CID':cv['CID'].astype(str),'cov_row':np.arange(len(cv['CID']))})
    assert not cov_index.duplicated(keys).any();assert list(cv['mag_parameter_order'])==['mx','x1','c','t0']
    merged=meta.merge(refit,on=keys,how='left',suffixes=('_reference',''),validate='one_to_one').merge(cov_index,on=keys,how='left',validate='one_to_one').merge(folds[keys+['fold','field','depth']],on=keys,how='left',validate='one_to_one')
    excluded=merged.loc[merged.cov_row.isna(),keys].copy();excluded['reason']='No successful coherent-covariance fit in this particular arm; never replaced by a rounded metadata covariance.'
    use=merged.loc[merged.cov_row.notna()].copy();idx=use.cov_row.astype(int).to_numpy();pars=cv['mag_parameters'][idx,:3].copy();pars[:,0]+=10.635;cov=cv['mag_covariance'][idx,:3,:3]
    expected=use[['mB_double','x1_double','c_double']].to_numpy();np.testing.assert_allclose(pars,expected,atol=2e-12,rtol=1e-12)
    assert np.isfinite(cov).all() and (np.linalg.eigvalsh(cov).min(axis=1)>0).all()
    mass=use.HOST_LOGMASS.to_numpy();masserr=use.HOST_LOGMASS_ERR.to_numpy();assert (masserr>0).all() and (mass>0).all()
    hostprob=ndtr((mass-10.)/masserr)
    pia=use.PROB_SNNV19.to_numpy().copy();lowz=use.IDSURVEY.to_numpy()!=10
    assert (use.loc[lowz,'TYPE_reference']==1).all();pia[lowz]=1.
    assert ((pia>=0)&(pia<=1)).all()
    out=ROOT/'phase2/hierarchy/data'/a.name;out.mkdir(parents=True,exist_ok=True)
    data={'y':pars,'cov':cov,'z':use.zHD_reference.to_numpy(),'zhel':use.zHEL_reference.to_numpy(),'zerr':use.zHDERR_reference.to_numpy(),'host_prob':hostprob,'host_mass':mass,'host_mass_error':masserr,'pIa':pia,'survey':use.IDSURVEY.to_numpy().astype(int),'cid_index':use.CIDint.to_numpy().astype(int),'fold':use.fold.to_numpy().astype(int),'reference_y':use[['mB_reference','x1_reference','c_reference']].to_numpy()}
    np.savez_compressed(out/'data.npz',**data)
    table=use[keys+['CIDint','TYPE_reference','field','depth','fold','cov_row']].copy();table['pIa_input']=pia;table['host_prob_above_10']=hostprob;table['delta_mB']=pars[:,0]-data['reference_y'][:,0];table['delta_x1']=pars[:,1]-data['reference_y'][:,1];table['delta_c']=pars[:,2]-data['reference_y'][:,2];table.to_csv(out/'rows.csv',index=False);excluded.to_csv(out/'excluded.csv',index=False)
    qualification={'status':'prepared_not_inferred','selection':'Not specified or validated here. The synthetic probit selector is not a DES selector. No real-population model ranking has been run by this script.','observables':'Actual independently refitted pre-BBC mB,x1,c and coherent Hessian covariance. No BBC bias correction or released MU error inserted.','reference':'Original published membership fixed before alternatives. Missing fits are explicit; matched-comparison rules still apply.','host_membership':'Probability above logmass10 from Gaussian measured mass error with locally flat mass prior. This is a declared approximation, not a measured physical binary host property.','classification':'DES: published SNNV19 probabilities, which require calibration and a contaminant model. Spectroscopically confirmed TYPE1 external low-z sample: probability set to1 explicitly. No BEAMS posterior probability reused as a prior.','uncertainties_not_yet_added':'Shared calibration, source-model uncertainty beyond the fitter, redshift/peculiar-velocity/lensing contribution to distance, host prior/matching uncertainty, selection and contamination. A positive Hessian is not a validation of the flux likelihood.','refit_gate':'Read official reproduction report for approximation and clipping status of this named arm.'}
    (out/'contract.json').write_text(json.dumps({'arguments':vars(a),'n_reference':len(meta),'n_prepared':len(use),'n_excluded':len(excluded),'qualification':qualification},indent=2)+'\n')
    files=[refitpath,covpath,metapath,foldpath,Path(__file__)];sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
    (out/'manifest.json').write_text(json.dumps({'time_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'inputs_sha256':{str(f.relative_to(ROOT)):sha(f) for f in files},'outputs_sha256':{str(f.relative_to(ROOT)):sha(f) for f in out.iterdir() if f.name!='manifest.json'}},indent=2)+'\n')
    print(json.dumps({'n_prepared':len(use),'n_excluded':len(excluded),'path':str(out)}))
if __name__=='__main__':main()
