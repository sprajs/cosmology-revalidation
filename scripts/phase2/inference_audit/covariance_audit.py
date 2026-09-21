#!/usr/bin/env python3
from pathlib import Path
import json,numpy as np,pandas as pd
from scipy.linalg import eigvalsh
R=Path(__file__).resolve().parents[3];O=R/'phase2/inference_audit'
S=json.load(open(R/'phase2/official/portable_pilot/double_parameters_and_hessian.json'));F=json.load(open(R/'phase2/independent_flux/independent_refits.json'));rows=[]
for r in F:
 if r.get('interpolation')!='sncosmo':continue
 cid=str(r['CID']);p=np.array(S[cid]['parameters_x0_x1_c_t0']);sv=np.array(S[cid]['hessian_covariance']);iv=np.array(r['hessian_covariance']);js=np.array([-2.5/np.log(10)/p[0],1,1,1]);ji=np.array([-2.5/np.log(10),1,1,1]);A=(sv*js[:,None]*js[None,:])[:3,:3];B=(iv*ji[:,None]*ji[None,:])[:3,:3];e=eigvalsh(B,A)
 row={'CID':cid,'family':r['family'],'generalized_eigen_min':float(e.min()),'generalized_eigen_max':float(e.max()),'marginal_sd_ratios':np.sqrt(np.diag(B)/np.diag(A)).tolist(),'projections':[]}
 for alpha in [.05,.15,.3]:
  for beta in [1.,3.,5.]:
   # Mean m=-alpha*x+beta*c => residual vector[1,+alpha,-beta].
   v=np.array([1,alpha,-beta]);row['projections'].append({'alpha':alpha,'beta':beta,'variance_ratio_independent_over_SNANA':float((v@B@v)/(v@A@v))})
 rows.append(row)
cohort=[]
for name in ['conditioned-multistart-best','recovered-mask']:
 p=R/'phase2/hierarchy/data'/name;a=np.load(p/'data.npz');d=pd.read_csv(p/'rows.csv',dtype={'CID':str});keep=(a['survey']==10)&(a['pIa']>.999);ids=dict(zip(d.loc[keep,'CID'],a['fold'][keep].astype(int).tolist()));canonical=json.load(open(R/'phase2/hierarchy/conditional-cohort.json'))['cid_to_fold'];assert ids==canonical
 v=a['cov'][keep];assert np.all(np.linalg.eigvalsh(v)>0);cohort.append({'arm':name,'n':int(keep.sum()),'n_test':int((keep&(a['fold']==0)).sum()),'frozen_membership_identical':True,'measurement_sd_quantiles':{k:np.quantile(np.sqrt(v[:,j,j]),[.05,.5,.95]).tolist() for j,k in enumerate(['mB','x1','c'])},'near_colour_hinge_within_1sigma_count':int(np.sum(abs(a['y'][keep,2])<np.sqrt(v[:,2,2]))),'training_multimodal_retained':bool(((d.CID=='1307748')&keep&(a['fold']!=0)).any())})
result={'pilot_geometry':rows,'cohorts':cohort,'global_generalized_eigen_range':[min(r['generalized_eigen_min'] for r in rows),max(r['generalized_eigen_max'] for r in rows)],'interpretation':'Independent matrix transform/recalculation. Pilot covariance agreement is directional and local; scalar1.20 stress is not a calibrated surveywide covariance prior. No real heldout model scores read.'}
(O/'covariance_audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
