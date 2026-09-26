"""Discovery-conditioned shared latent predictive audit; no validation fitting."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from scipy.linalg import solve_triangular
from scipy.stats import norm
P=Path(__file__).resolve().parent;O=P/'shared1020';VPATH=P/'validation1020';m=pd.read_csv(VPATH/'cohort.csv',dtype={'CID':str});disc=np.load(P/'shared43'/'projected-modes.npz');dc=np.load(VPATH/'frozen-discovery-coefficients.npz');c=dc['basis_mean']
assert not (O/'result.json').exists()
assert len(m)==1020 and m.CID.is_unique and set(m.CID).isdisjoint(set(disc['CID']))
discovery64=pd.read_csv(P.parents[2]/'runs/salt_dust_audit/flux_response/selected_objects.csv',dtype={'CID':str})
assert set(m.CID).isdisjoint(set(discovery64.CID))
checks=[];res=[];cal=[];obs=[];Fs=[];us=[];offsets=[0]
for count,row in enumerate(m.itertuples()):
 cid=row.CID;nom=np.load(VPATH/'objectives'/f'objective_{cid}.npz');cached=np.load(VPATH/'analysis'/'objects'/f'{cid}.npz')
 order=pd.DataFrame({'MJD':nom['MJD'],'band':nom['band']}).sort_values(['MJD','band'],kind='stable').index.to_numpy()
 assert len(order)==row.expected_epochs
 assert np.array_equal(nom['model_flux'][order],cached['official_flux_model'])
 assert np.array_equal(nom['data_flux'][order],cached['observed_flux'])
 assert np.array_equal(nom['MJD'][order],cached['MJD']) and np.array_equal(nom['band'][order],cached['band'])
 assert np.array_equal(nom['data_fluxerr'][order],cached['quoted_error'])
 assert np.array_equal(nom['frozen_flux_covariance'][np.ix_(order,order)],cached['exact_covariance'])
 C=cached['exact_covariance'];J=cached['jacobian_flux'];L=np.linalg.cholesky(C);jw=solve_triangular(L,J,lower=True);U,s,_=np.linalg.svd(jw,full_matrices=True);assert np.sum(s>s[0]*1e-10)==4;Q=U[:,4:]
 r=Q.T@solve_triangular(L,(nom['data_flux']-nom['model_flux'])[order],lower=True);assert np.allclose(r,cached['published_mask_r'],rtol=0,atol=1e-10)
 T=cached['published_mask_T'][:,:3];columns=[]
 for k in range(1,10):
  q=np.load(O/f'model{k:03d}'/'objectives'/f'objective_{cid}.npz')
  for key in ['MJD','band','data_flux','data_fluxerr','zHEL','MWEBV','parameters_x0_x1_c_t0']:assert np.array_equal(q[key],nom[key]),(cid,k,key)
  assert np.isfinite(q['model_flux']).all() and np.isfinite(q['frozen_flux_covariance']).all()
  delta=q['model_flux'][order]-nom['model_flux'][order];columns.append(delta)
  checks.append(dict(CID=cid,model=k,epochs=len(delta),parameter_data_epoch_exact_identity=True,max_absolute_fractional_prediction_shift=float(max(abs(delta/nom['model_flux'][order])))))
 W=.3*Q.T@solve_triangular(L,np.column_stack(columns),lower=True);X=np.column_stack([W,T]);assert np.isfinite(X).all() and np.isfinite(r).all();Fs.append(X.T@X);us.append(X.T@r)
 res.append(r);cal.append(W);obs.append(T);offsets.append(offsets[-1]+len(r))
 if count%100==0:print('mode projection',count+1,'of',len(m),flush=True)
y=np.concatenate(res);A=np.vstack(cal);T=np.vstack(obs);FrawD=disc['joint_F'].sum(0);grawD=disc['joint_u'].sum(0);FrawV=np.array(Fs).sum(0);grawV=np.array(us).sum(0)
R=lambda F,g:float(.5*g@np.linalg.solve(np.eye(len(g))+F,g)-.5*np.linalg.slogdet(np.eye(len(g))+F)[1])
model_results={};save={};posterior={}
bandmap=np.array([[1,0,0],[-1,-1,-1],[0,1,0],[0,0,1.]])
physical_cov=2*.02**2*np.linalg.inv(bandmap.T@bandmap);physical_L=np.linalg.cholesky(physical_cov)
M0=np.eye(12)[:,:9];M1=np.diag(np.r_[np.ones(9),np.ones(3)*.02]);M2=np.eye(12);M2[9:,9:]=physical_L
for name,M in [('calibration_only',M0),('calibration_plus_observer',M1),('calibration_plus_isotropic_griz_observer',M2)]:
 inds=np.arange(M.shape[1]);FD=M.T@FrawD@M;gD=M.T@grawD;FV=M.T@FrawV@M;gV=M.T@grawV
 VD=np.linalg.inv(np.eye(len(inds))+FD);mu=VD@gD;X=np.column_stack([A,T])@M
 gain=R(FD+FV,gD+gV)-R(FD,gD)
 # Independent direct posterior predictive density in QR observation subspace.
 Q,B=np.linalg.qr(X,mode='reduced');yp=Q.T@y;C=np.eye(len(inds))+B@VD@B.T;L=np.linalg.cholesky(C);w=solve_triangular(L,yp-B@mu,lower=True)
 direct=float(-.5*w@w-np.log(np.diag(L)).sum()+.5*yp@yp)
 assert abs(gain-direct)<1e-8,(name,gain,direct)
 prediction=X@mu;postdf=pd.DataFrame({'CID':disc['CID']})
 model_results[name]={'joint_validation_predictive_gain_over_zero':gain,'independent_QR_predictive_gain':direct,'independent_error':direct-gain,
   'discovery_latent_posterior_mean':mu.tolist(),'discovery_latent_posterior_sd':np.sqrt(np.diag(VD)).tolist(),
   'fixed_discovery_posterior_mean_validation_gain':float(y@prediction-.5*prediction@prediction),
   'discovery_log_marginal_gain_over_zero':R(FD,gD),'validation_noise_design_information_trace':float(np.trace(FV)),
   'conditioning':'Shared joint latent posterior learned only on43 discovery objects; shared posterior uncertainty integrated once on1020 validation objects.'}
 save[name+'_discovery_posterior_covariance']=VD;save[name+'_discovery_posterior_mean']=mu;save[name+'_validation_prediction']=prediction
 posterior[name]=(mu,VD)
mu,VD=posterior['calibration_only'];v=T@c;e=y-A@mu;B=np.linalg.inv(np.linalg.inv(VD)+A.T@A)
apply=lambda w:w-A@B@(A.T@w)
I=float(v@apply(v));av=float(v@apply(e));gain=av-.5*I;z=av/np.sqrt(I)
# Physical-space coefficient posterior under combined model: last3 remain original contrast basis.
combined=np.array(model_results['calibration_plus_observer']['discovery_latent_posterior_mean']);bandmap=np.array([[1,0,0],[-1,-1,-1],[0,1,0],[0,0,1.]])
result={'scope':'Post-validation exploratory attribution through nine native coupled calibration/SALT finite variants; discovery-conditioned shared Gaussian second-moment approximation, not full release-systematic posterior.',
 'discovery_objects':43,'validation_objects':len(m),'validation_epochs':int(m.expected_epochs.sum()),'validation_projected_dimension':len(y),'nine_mode_parameter_data_epoch_identity':True,
 'models':model_results,'conditional_predictive_log_ratio_calibration_plus_observer_vs_calibration':model_results['calibration_plus_observer']['joint_validation_predictive_gain_over_zero']-model_results['calibration_only']['joint_validation_predictive_gain_over_zero'],
 'original_frozen_observer_direction_after_calibration_null_conditioning':{'score_gain_for_conditional_mean_shift':gain,'matched_filter_centered':av,'information':I,'conditional_Gaussian_Z':z,'conditional_Gaussian_one_sided_tail':float(norm.sf(z)),'qualification':'Direction fixed by discovery, residual centred on calibration-only discovery posterior and shared predictive covariance used. Partial nine-mode model; not original validation significance or full-budget rejection.'},
 'joint_model_discovery_observer_griz_mean_mag':(bandmap@(.02*combined[-3:])).tolist(),
 'isotropic_prior_sensitivity':{'coefficient_covariance':physical_cov.tolist(),'physical_band_covariance':(bandmap@physical_cov@bandmap.T).tolist(),'inherited_physical_prior_total_variance':float(6*.02**2),'isotropic_physical_prior_total_variance':float(np.trace(bandmap@physical_cov@bandmap.T)),'conditional_predictive_log_ratio_vs_calibration':model_results['calibration_plus_isotropic_griz_observer']['joint_validation_predictive_gain_over_zero']-model_results['calibration_only']['joint_validation_predictive_gain_over_zero']},
 'limitations':['Nine finite variations interpreted as Gaussian modes with source amplitude scale .3, no renormalization. This is a new epoch-level second-moment model.','The same latent variables connect discovery and validation; no independent validation covariance inflation is used.','CALSPEC and other released systematics omitted. Fragilistic sample coverage is anisotropic; wavelength-shift covariance/seed mapping not established.','Nominal fixed masks, C and local tangent; no selection/classification/BBC/SALT retraining likelihood regeneration.','The observer prior remains .02^2I in inherited [g-r,i-r,z-r] coordinates, not isotropic in an orthonormal zero-sum4-band gauge.','All model definitions here were frozen after the earlier validation result, so this is exploratory attribution rather than pristine hypothesis confirmation.']}
np.savez_compressed(O/'projected-modes.npz',calibration_design=A,observer_design=T,residual=y,joint_F=np.array(Fs),joint_u=np.array(us),CID=m.CID.to_numpy(dtype='U'),field=m.field.to_numpy(dtype='U'),row_offsets=np.array(offsets),**save)
pd.DataFrame(checks).to_csv(O/'variant-gates.csv',index=False)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();result['source_sha256']=sha(__file__);result['protocol_sha256']=sha(P/'shared1020-protocol.md')
result['input_sha256']={str(p.relative_to(P)):sha(p) for p in [O/'input-manifest.json',O/'prior-amendment-manifest.json',P/'shared43'/'projected-modes.npz',VPATH/'frozen-discovery-coefficients.npz']+[O/f'model{k:03d}'/'objectives'/'manifest.json' for k in range(1,10)]}
(O/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
