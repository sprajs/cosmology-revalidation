"""Native12-mode shared predictive sensitivity with explicit CALSPEC unit map."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from scipy.linalg import solve_triangular
from scipy.stats import norm
P=Path(__file__).resolve().parent;O=P/'expanded12';assert not (O/'result.json').exists()
units=json.loads((O/'calspec-unit-map.json').read_text())['bands'];d9=np.load(P/'shared43'/'projected-modes.npz');v9=np.load(P/'shared1020'/'projected-modes.npz');dmat=np.load(P/'exact43'/'comparison'/'matched-matrices.npz');c=np.load(P/'validation1020'/'frozen-discovery-coefficients.npz')['basis_mean'];checks=[];samples={}
for cohort,old in [('discovery',d9),('validation',v9)]:
 ids=[str(s) for s in old['CID']];rs=[];As=[];Ts=[];Fs=[];gs=[];offsets=[0]
 for i,cid in enumerate(ids):
  if cohort=='discovery':
   q=np.load(P/'exact43'/'objectives'/f'objective_{cid}.npz');pp=cid+'__';C=dmat[pp+'exact_covariance'];J=dmat[pp+'jacobian_flux'].copy();A0=old[pp+'calibration_design'];T=old[pp+'observer_design'];r=old[pp+'r']
  else:
   q=np.load(P/'validation1020'/'objectives'/f'objective_{cid}.npz');ca=np.load(P/'validation1020'/'analysis'/'objects'/f'{cid}.npz');C=ca['exact_covariance'];J=ca['jacobian_flux'].copy();lo,hi=old['row_offsets'][i:i+2];A0=old['calibration_design'][lo:hi];T=old['observer_design'][lo:hi];r=old['residual'][lo:hi]
  order=pd.DataFrame({'MJD':q['MJD'],'band':q['band']}).sort_values(['MJD','band'],kind='stable').index.to_numpy();f=q['model_flux'][order];J[:,0]=-.4*np.log(10)*f
  assert np.array_equal(C,q['frozen_flux_covariance'][np.ix_(order,order)])
  L=np.linalg.cholesky(C);U,s,_=np.linalg.svd(solve_triangular(L,J,lower=True),full_matrices=True);assert np.sum(s>s[0]*1e-10)==4;Q=U[:,4:]
  assert np.allclose(Q.T@solve_triangular(L,(q['data_flux']-q['model_flux'])[order],lower=True),r,rtol=0,atol=1e-10)
  columns=[]
  for name,amp in [('CALSPEC',1.),('MWEBV',1.),('COLORLAW',.3)]:
   a=np.load(O/cohort/name/'objectives'/f'objective_{cid}.npz')
   for key in ['MJD','band','zHEL','parameters_x0_x1_c_t0']:assert np.array_equal(a[key],q[key]),(cohort,cid,name,key)
   unit_scale=np.ones(len(f));d={}
   if name=='CALSPEC':
    unit_scale=np.array([units[b]['flux_scale_float32'] for b in q['band']]);data_expected=(q['data_flux'].astype(np.float32)*unit_scale.astype(np.float32)).astype(float);err_expected=(q['data_fluxerr'].astype(np.float32)*unit_scale.astype(np.float32)).astype(float)
    data_diff=a['data_flux']-data_expected;err_diff=a['data_fluxerr']-err_expected
    assert np.allclose(a['data_flux'],data_expected,rtol=2e-7,atol=1e-8),(cid,'CALSPEC flux unit map')
    assert np.allclose(a['data_fluxerr'],err_expected,rtol=2e-7,atol=1e-8),(cid,'CALSPEC error unit map')
    assert np.array_equal(a['MWEBV'],q['MWEBV'])
    d={'data_source_map_max_absolute_error':float(max(abs(data_diff))),'quoted_error_source_map_max_absolute_error':float(max(abs(err_diff))),'inverse_map_data_rounding_max_measurement_sigma':float(max(abs(a['data_flux']/unit_scale-q['data_flux'])/q['data_fluxerr'])),'model_flux_changed_before_unit_mapping':bool(not np.array_equal(a['model_flux'],q['model_flux']))}
   else:
    assert np.array_equal(a['data_flux'],q['data_flux']) and np.array_equal(a['data_fluxerr'],q['data_fluxerr'])
    if name=='MWEBV':
     expected=(q['MWEBV'].astype(np.float32)*np.float32(.95)).astype(float);assert np.allclose(a['MWEBV'],expected,rtol=1e-7,atol=1e-10),(cid,'MWEBV intended input')
     d={'max_MWEBV_source_map_error':float(max(abs(a['MWEBV']-expected)))}
    else:assert np.array_equal(a['MWEBV'],q['MWEBV'])
   mapped=a['model_flux']/unit_scale;delta=mapped[order]-f;assert np.isfinite(delta).all();columns.append(amp*Q.T@solve_triangular(L,delta,lower=True))
   checks.append(dict(cohort=cohort,CID=cid,mode=name,amplitude_scale=amp,epochs=len(f),parameter_epoch_band_z_identity=True,observed_data_units='CALSPEC source scale mapped back' if name=='CALSPEC' else 'exact nominal identity',**d))
  A=np.column_stack([A0,np.column_stack(columns)]);X=np.column_stack([A,T]);assert np.isfinite(X).all();rs.append(r);As.append(A);Ts.append(T);Fs.append(X.T@X);gs.append(X.T@r);offsets.append(offsets[-1]+len(r))
  if i%200==0:print(cohort,'processed',i+1,'of',len(ids),flush=True)
 samples[cohort]=dict(CID=np.array(ids,dtype='U'),field=old['field'],residual=np.concatenate(rs),calibration_design=np.vstack(As),observer_design=np.vstack(Ts),joint_F=np.array(Fs),joint_u=np.array(gs),row_offsets=np.array(offsets))
 np.savez_compressed(O/f'{cohort}-projected-modes.npz',**samples[cohort])
assert set(samples['discovery']['CID']).isdisjoint(set(samples['validation']['CID']))
D=samples['discovery'];V=samples['validation'];FD0=D['joint_F'].sum(0);gD0=D['joint_u'].sum(0);FV0=V['joint_F'].sum(0);gV0=V['joint_u'].sum(0);A=V['calibration_design'];T=V['observer_design'];y=V['residual'];B=np.array([[1,0,0],[-1,-1,-1],[0,1,0],[0,0,1.]])
prior=2*.02**2*np.linalg.inv(B.T@B);M0=np.eye(15)[:,:12];M1=np.diag(np.r_[np.ones(12),np.ones(3)*.02]);M2=np.eye(15);M2[12:,12:]=np.linalg.cholesky(prior)
R=lambda F,g:float(.5*g@np.linalg.solve(np.eye(len(g))+F,g)-.5*np.linalg.slogdet(np.eye(len(g))+F)[1]);models={};post={};save={}
for name,M in [('systematics_only',M0),('systematics_plus_observer',M1),('systematics_plus_isotropic_griz_observer',M2)]:
 FD=M.T@FD0@M;gD=M.T@gD0;FV=M.T@FV0@M;gV=M.T@gV0;Cpost=np.linalg.inv(np.eye(len(gD))+FD);mean=Cpost@gD;X=np.column_stack([A,T])@M
 gain=R(FD+FV,gD+gV)-R(FD,gD);Q,H=np.linalg.qr(X,mode='reduced');yp=Q.T@y;Cp=np.eye(len(gD))+H@Cpost@H.T;L=np.linalg.cholesky(Cp);z=solve_triangular(L,yp-H@mean,lower=True);direct=float(-.5*z@z-np.log(np.diag(L)).sum()+.5*yp@yp);assert abs(gain-direct)<1e-8
 models[name]={'joint_validation_predictive_gain_over_zero':gain,'independent_QR_gain':direct,'independent_error':direct-gain,'discovery_latent_posterior_mean':mean.tolist(),'discovery_latent_posterior_sd':np.sqrt(np.diag(Cpost)).tolist(),'discovery_log_marginal_gain_over_zero':R(FD,gD)};post[name]=(mean,Cpost);save[name+'_discovery_posterior_mean']=mean;save[name+'_discovery_posterior_covariance']=Cpost
mean,Cpost=post['systematics_only'];inv=np.linalg.inv(np.linalg.inv(Cpost)+A.T@A);apply=lambda v:v-A@inv@(A.T@v);v=T@c;e=y-A@mean;I=float(v@apply(v));score=float(v@apply(e));Z=score/np.sqrt(I)
checks=pd.DataFrame(checks);checks.to_csv(O/'variant-gates.csv',index=False);np.savez_compressed(O/'posterior-models.npz',**save)
result={'scope':'Expanded12-mode exploratory sensitivity to existing native calibration/training/CALSPEC/MWEBV/COLORLAW variations, correct shared discovery conditioning, not a complete released-systematic model.',
 'discovery_objects':43,'validation_objects':1020,'validation_epochs':39606,'mode_order':[f'cal_{j}' for j in range(1,10)]+['CALSPEC','MWEBV','COLORLAW'],'mode_amplitude_scales':[.3]*9+[1.,1.,.3],
 'models':models,'conditional_predictive_log_ratio_extra_observer':models['systematics_plus_observer']['joint_validation_predictive_gain_over_zero']-models['systematics_only']['joint_validation_predictive_gain_over_zero'],
 'isotropic_prior_conditional_predictive_log_ratio_extra_observer':models['systematics_plus_isotropic_griz_observer']['joint_validation_predictive_gain_over_zero']-models['systematics_only']['joint_validation_predictive_gain_over_zero'],
 'original_fixed_observer_after_systematic_null_conditioning':{'centered_matched_product':score,'information':I,'conditional_Gaussian_Z':Z,'conditional_Gaussian_one_sided_tail':float(norm.sf(Z)),'full_fixed_direction_shift_gain':score-.5*I,'qualification':'Discovery-selected fixed direction, conditional partial12-mode null mean/covariance; exploratory not universal significance.'},
 'CALSPEC_unit_map_gates':{'max_data_source_map_absolute_error':float(checks.data_source_map_max_absolute_error.max()),'max_error_source_map_absolute_error':float(checks.quoted_error_source_map_max_absolute_error.max()),'max_inverse_data_rounding_in_measurement_sigma':float(checks.inverse_map_data_rounding_max_measurement_sigma.max()),'model_changed_before_mapping':bool(checks.model_flux_changed_before_unit_mapping.dropna().any()),'formula':'delta_m=.00714*transmission_weighted_mean_A/10000; a=10^(-.4delta_m); native data/error multiplied by a; mean divided by a to nominal units.'},
 'MWEBV_source_map_max_error':float(checks.max_MWEBV_source_map_error.max()),'limitations':['Finite native mean changes projected with nominal local tangent and C; variant nonlinear coordinate transport and covariance changes not regenerated.','Gaussian second-moment interpretation of saved scales is an epoch-level approximation, not complete original latent posterior.','Other published systematic modes, simulation/BBC, classifier/selection and SALT-training likelihood remain outside scope.','No physical zero-point/extinction correction or changed cosmology is identified by these conditional predictions.','Inherited and equal-total-variance isotropic observer priors both retained without selecting the favourable result.']}
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();result['source_sha256']=sha(__file__);result['protocol_sha256']=sha(P/'expanded12-protocol.md');result['inputs_sha256']={str(p.relative_to(P)):sha(p) for p in [O/'input-manifest.json',O/'calspec-unit-map.json',P/'shared43'/'projected-modes.npz',P/'shared1020'/'projected-modes.npz']+[O/cohort/name/'objectives'/'manifest.json' for cohort in ['discovery','validation'] for name in ['CALSPEC','MWEBV','COLORLAW']]}
(O/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
