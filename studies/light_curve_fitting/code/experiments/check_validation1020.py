"""Independent profile-chi-square and QR predictive checks for frozen transfer."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from scipy.linalg import solve_triangular
P=Path(__file__).resolve().parent/'validation1020';O=P/'analysis';m=pd.read_csv(P/'cohort.csv',dtype={'CID':str});r=json.loads((O/'result.json').read_text());dc=np.load(P/'frozen-discovery-coefficients.npz');dr=np.load(P/'frozen-rest-discovery-coefficients.npz');c,V=dc['basis_mean'],dc['basis_covariance'];cr,Vr=dr['basis_mean'],dr['basis_covariance'];K=.4*np.log(10)
rows=[];stores={arm:{'r':[],'T':[]} for arm in r['arms']}
for cid in m.CID:
    x=np.load(O/'objects'/f'{cid}.npz');f=x['official_flux_model'];b=x['band'];off=np.array([c[0],-sum(c),c[1],c[2]])
    mode=-K*f*np.array([off['griz'.index(bb)] for bb in b])
    for arm,ids in [('published_mask',np.arange(len(f))),('first_exact_duplicate',x['first_duplicate_keep'])]:
        C=x['exact_covariance'][np.ix_(ids,ids)];J=x['jacobian_flux'][ids];L=np.linalg.cholesky(C)
        A=solve_triangular(L,J,lower=True);y=solve_triangular(L,(x['observed_flux']-f)[ids],lower=True);d=solve_triangular(L,mode[ids],lower=True)
        baseline=y-A@np.linalg.lstsq(A,y,rcond=None)[0];changed=y-d-A@np.linalg.lstsq(A,y-d,rcond=None)[0]
        gain=.5*(baseline@baseline-changed@changed);saved=x[arm+'_T'][:,:3]@c;yy=x[arm+'_r'];expected=yy@saved-.5*saved@saved
        rows.append({'CID':cid,'arm':arm,'gain_error':float(gain-expected),'chi_error':float(baseline@baseline-yy@yy)})
        stores[arm]['r'].append(yy);stores[arm]['T'].append(x[arm+'_T'])
checks={}
for arm,d in stores.items():
    yy=np.concatenate(d['r']);tt=np.vstack(d['T'])
    for fam,ii,cc,vv,saved in [('observer',slice(0,3),c,V,r['arms'][arm]['joint_discovery_posterior_predictive_gain']),('rest',slice(3,6),cr,Vr,r['arms'][arm]['frozen_rest_secondary']['joint_discovery_posterior_predictive_gain'])]:
        Q,R=np.linalg.qr(tt[:,ii],mode='reduced');yc=Q.T@yy;cov=np.eye(3)+R@vv@R.T;L=np.linalg.cholesky(cov);w=solve_triangular(L,yc-R@cc,lower=True)
        gain=-.5*w@w-np.log(np.diag(L)).sum()+.5*yc@yc;checks[arm+'_'+fam]=float(gain-saved)
assert max(abs(x['gain_error']) for x in rows)<1e-8
assert max(abs(x['chi_error']) for x in rows)<1e-8
assert max(abs(v) for v in checks.values())<1e-8
out={'objects':len(m),'paired_profile_checks':len(rows),'max_profile_gain_error':max(abs(x['gain_error']) for x in rows),'max_profile_chi2_error':max(abs(x['chi_error']) for x in rows),'joint_predictive_QR_errors':checks,'method':'Independent raw-C weighted least-squares profile difference for every object and duplicate arm; orthogonal3D Gaussian covariance Cholesky for each shared predictive marginal.','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(O/'independent-check.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
