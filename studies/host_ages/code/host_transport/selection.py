#!/usr/bin/env python3
"""Field-held-out audit of actual eight-band host-measurement availability."""
from pathlib import Path
import json,hashlib,datetime,sys
import numpy as np,pandas as pd
from scipy.optimize import minimize
from scipy.special import expit
from scipy.stats import rankdata
ROOT=Path(__file__).resolve().parents[4];sys.path.insert(0,str(ROOT));from lib.records import fitres
W=ROOT/'.work/host-transport';O=ROOT/'studies/host_ages/results/host_transport';H=Path(__file__).parent

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    sp=W/'des-smp-hosts.csv';qp=W/'des-deep-crosswalk.csv';mp=W/'des/DES-Dovekie_Metadata.csv'
    s=pd.read_csv(sp,dtype={'SNID':str});q=pd.read_csv(qp,dtype={'SNID':str});m=fitres(mp);m=m[m.IDSURVEY==10].reset_index()
    d=m[['CID','FIELD','zHEL']].merge(s,left_on='CID',right_on='SNID',validate='one_to_one');d['field']=d.FIELD.map(lambda value: next((part for part in value.replace('SN-','').split('+') if part in ['C3','E2','X3']), 'outside'))
    d=d[d.field.isin(['C3','E2','X3'])].copy();full={k:len(g) for k,g in d.groupby('field')}
    avail=set(q.loc[q.primary_match&q.selected_Dovekie&q.deep_quality&q.host_dlr_lt4,'SNID']);d['available']=d.CID.isin(avail)
    d['host_i']=d.HOSTGAL_MAG_i;d['host_g_minus_i']=d.HOSTGAL_MAG_g-d.HOSTGAL_MAG_i;d['redshift']=d.zHEL
    ok=np.isfinite(d[['host_i','host_g_minus_i','redshift']]).all(axis=1)&d.HOSTGAL_MAG_i.between(0,50)&d.HOSTGAL_MAG_g.between(0,50)
    excluded=d[~ok][['CID','field','available']].to_dict('records');d=d[ok].reset_index(drop=True)
    cols=['redshift','host_i','host_g_minus_i'];x=d[cols].to_numpy();y=d.available.to_numpy(float);pred=np.full(len(d),np.nan);validation=[]
    for field in ['C3','E2','X3']:
        train=d.field.to_numpy()!=field;test=~train;mean=x[train].mean(0);std=x[train].std(0);a=np.column_stack([np.ones(train.sum()),(x[train]-mean)/std]);b=np.column_stack([np.ones(test.sum()),(x[test]-mean)/std]);t=y[train]
        def loss(coef):
            eta=a@coef;return float(np.sum(np.logaddexp(0,eta)-t*eta)+.5*np.sum(coef[1:]**2)),a.T@(expit(eta)-t)+np.r_[0,coef[1:]]
        opt=minimize(loss,np.zeros(4),jac=True,method='BFGS',options={'gtol':1e-8});assert np.max(abs(loss(opt.x)[1]))<1e-5
        p=expit(b@opt.x);pred[test]=p;truth=y[test];n1=truth.sum();n0=len(truth)-n1;auc=(rankdata(p)[truth==1].sum()-n1*(n1+1)/2)/(n1*n0)
        validation.append({'heldout_field':field,'n':int(test.sum()),'available':int(n1),'observed_fraction':float(truth.mean()),'mean_predicted_fraction':float(p.mean()),'Brier':float(np.mean((p-truth)**2)),'constant_training_fraction_Brier':float(np.mean((t.mean()-truth)**2)),'AUC':float(auc),'coefficients':opt.x.tolist(),'training_mean':mean.tolist(),'training_std':std.tolist(),'p_range': [float(p.min()),float(p.max())]})
    assert np.isfinite(pred).all();d['heldout_availability_probability']=pred;d.to_csv(W/'des-eightband-availability.csv',index=False)
    w=1/pred[y==1];w/=w.sum();means={c:{'full_selected_field_mean':float(d[c].mean()),'available_unweighted_mean':float(d.loc[d.available,c].mean()),'available_IPW_mean':float(w@d.loc[d.available,c])} for c in cols}
    result={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'method':json.loads((H/'selection-design.json').read_text()),'full_selected_field_counts':full,'excluded_missing_corrected_g_or_i':excluded,'analysed_selected_events':len(d),'available_events':int(y.sum()),'heldout_fields':validation,'descriptive_IPW_ESS':float(1/(w@w)),'descriptive_IPW_largest_weight':float(w.max()),'means':means,'inference_limit':'IPW is only a measured-covariate diagnostic: no missing-at-random assumption established and no true-age distribution transported. Holdout failure/miscalibration must not be hidden by refitting target fields.','code_sha256':{str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),H/'selection-design.json']},'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in [sp,qp,mp]},'output_sha256':sha(W/'des-eightband-availability.csv')}
    (O/'des-selection.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
