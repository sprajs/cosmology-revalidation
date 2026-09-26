"""Actual native coupled mode geometry and a labelled9-mode covariance stress."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from scipy.linalg import solve_triangular,eigvalsh
P=Path(__file__).resolve().parent;ROOT=P.parents[2];O=P/'shared43';D=P/'exact43'/'comparison';x=np.load(D/'matched-matrices.npz');m=pd.read_csv(D/'object-checks.csv',dtype={'CID':str});dc=np.load(P/'validation1020'/'frozen-discovery-coefficients.npz');c=dc['basis_mean'];save={};checks=[];rs=[];Ws=[];Ts=[];coeff=[];gram=[];uu=[]
for row in m.itertuples():
 cid=row.CID;p=cid+'__';q=np.load(O/'model000'/'objectives'/f'objective_{cid}.npz');order=pd.DataFrame({'MJD':q['MJD'],'band':q['band']}).sort_values(['MJD','band'],kind='stable').index.to_numpy()
 f=q['model_flux'][order];C=x[p+'exact_covariance'];J=x[p+'jacobian_flux'].copy();J[:,0]=-.4*np.log(10)*f;L=np.linalg.cholesky(C);jw=solve_triangular(L,J,lower=True);U,s,_=np.linalg.svd(jw,full_matrices=True);Q=U[:,4:]
 r=Q.T@solve_triangular(L,q['data_flux'][order]-f,lower=True);assert np.allclose(r,x[p+'official_mean_exactC_r'],atol=1e-10,rtol=0)
 T=x[p+'official_mean_exactC_T'][:,:3];columns=[]
 for k in range(1,10):
  a=np.load(O/f'model{k:03d}'/'objectives'/f'objective_{cid}.npz')
  check=dict(CID=cid,model=k,epochs=len(f),parameters_max_absolute=float(max(abs(a['parameters_x0_x1_c_t0']-q['parameters_x0_x1_c_t0']))),data_flux_max_absolute=float(max(abs(a['data_flux']-q['data_flux']))),data_error_max_absolute=float(max(abs(a['data_fluxerr']-q['data_fluxerr']))),mean_fractional_max=float(max(abs(a['model_flux']/q['model_flux']-1))),covariance_relative_frobenius=float(np.linalg.norm(a['frozen_flux_covariance']-q['frozen_flux_covariance'])/np.linalg.norm(q['frozen_flux_covariance'])))
  for key in ['MJD','band','data_flux','data_fluxerr','zHEL','MWEBV','parameters_x0_x1_c_t0']:assert np.array_equal(a[key],q[key]),(cid,k,key)
  columns.append(a['model_flux'][order]-f);checks.append(check)
 raw=np.column_stack(columns);W=.3*Q.T@solve_triangular(L,raw,lower=True);joint=np.column_stack([W,T]);gram.append(joint.T@joint);uu.append(joint.T@r)
 save[p+'r']=r;save[p+'calibration_design']=W;save[p+'observer_design']=T;save[p+'raw_coupled_flux_differences']=raw
 rs.append(r);Ws.append(W);Ts.append(T)
r=np.concatenate(rs);W=np.vstack(Ws);T=np.vstack(Ts);v=T@c;F=W.T@W;post=np.linalg.inv(np.eye(9)+F);a=post@W.T@r
kinv=lambda v:v-W@post@(W.T@v)
kv=kinv(v);kr=kinv(r);I=float(v@v);Is=float(v@kv);ap=float(v@r);aps=float(v@kr)
A,s,B=np.linalg.svd(W,full_matrices=False);rank=int(np.sum(s>s[0]*1e-10));vp=A[:,:rank]@(A[:,:rank].T@v)
# Original .3 scales are squared in the outer-product covariance. No centring/normalization.
modelbase=ROOT/'sources/repos/des-science__DES-SN5YR@1.3/2_LCFIT_MODEL/SALT3.DES5YR-SYS';shifts=[];magoffset=[]
for k in range(10):
 info=(modelbase/f'SALT3.MODEL{k:03d}'/'SALT3.INFO').read_text().splitlines();ms={};ws={}
 for ln in info:
  z=ln.split()
  if not z:continue
  if z[0]=='MAG_OFFSET:':magoffset.append({'model':k,'MAG_OFFSET':float(z[1])})
  if len(z)>=4 and z[1]=='DES' and z[0] in ['MAGSHIFT:','WAVESHIFT:']:(ms if z[0]=='MAGSHIFT:' else ws)[z[2][-1]]=float(z[3])
 if k:
  assert set(ms)==set(ws)==set('griz');shifts.append([ms[b] for b in 'griz'])
  log=(O/f'model{k:03d}'/'fit.log').read_text()
  for b in 'griz':assert f'Update DES-{b}' in log
shift=np.array(shifts);sample_cov=(.3*shift).T@(.3*shift)
fpath=ROOT/'sources/repos/PantheonPlusSH0ES__DataRelease@7fc6805/Pantheon+_Data/2_CALIBRATION/FRAGILISTIC_COVARIANCE.npz';fc=np.load(fpath,allow_pickle=True);cov=fc['cov'];labels=[str(a).strip() for a in fc['labels']];idx=[labels.index('DES5YR '+b) for b in 'griz'];block=cov[np.ix_(idx,idx)];eig=eigvalsh(sample_cov,block)
result={'scope':'Retrospective43-object native finite coupled calibration+surface attribution and approximate shared9-mode covariance stress, not full systematic-budget likelihood.',
 'objects':len(m),'epochs':sum(v['epochs'] for v in checks if v['model']==1),'nominal_parameter_data_epoch_identity':True,'native_MAG_OFFSET':magoffset,'signed_variants':9,'amplitude_scales':[.3]*9,'weight_normalization':'None; covariance=.09*sum(delta_j delta_j^T), uncentred.',
 'calibration_design_singular_values':s.tolist(),'rank':rank,'fixed_observer_prediction_norm_squared':I,'fixed_observer_span_explained_fraction':float(vp@vp/I),'unpenalized_observer_span_residual_norm_squared':float((v-vp)@(v-vp)),
 'nominal_fixed_observer_gain':ap-.5*I,'nine_mode_covariance_fixed_observer_gain':aps-.5*Is,'nominal_signed_matched_filter':ap/np.sqrt(I),'nine_mode_covariance_signed_matched_filter':aps/np.sqrt(Is),'nine_mode_information':Is,'information_retained_fraction':Is/I,
 'nine_mode_posterior_mean':a.tolist(),'nine_mode_posterior_sd':np.sqrt(np.diag(post)).tolist(),'nine_mode_log_marginal_gain_over_fixed_null':float(.5*r@W@post@W.T@r-.5*np.linalg.slogdet(np.eye(9)+F)[1]),
 'zeropoint_covariance_coverage':{'band_order':'griz','nine_scaled_uncentred_covariance':sample_cov.tolist(),'Fragilistic_DES_block':block.tolist(),'nine_to_Fragilistic_generalized_eigenvalues':eig.tolist(),'nine_sd':np.sqrt(np.diag(sample_cov)).tolist(),'Fragilistic_sd':np.sqrt(np.diag(block)).tolist(),'source_sha256':hashlib.sha256(fpath.read_bytes()).hexdigest(),'qualification':'Coverage diagnostic only: finite realized draws vs upstream marginal; no prior addition or renormalization.'},
 'limitations':['Epoch calibration covariance stress is a new Gaussian second-moment approximation using released finite variations and distance-level scale convention.','No CALSPEC or other systematics, nonlinear nuisance refits, selection/BBC or full upstream prior.','Discovery and validation share latent coordinates; validation requires discovery-conditioned predictive covariance, not independent inflation.','Finite native predictions use fixed nominal coordinates/C/J; variant covariance changes are recorded but not mixed into mean attribution.']}
save['joint_F']=np.array(gram);save['joint_u']=np.array(uu);save['CID']=m.CID.to_numpy(dtype='U');save['field']=m.field.to_numpy(dtype='U');save['coefficient_calibration_posterior_mean']=a;save['coefficient_calibration_posterior_covariance']=post
np.savez_compressed(O/'projected-modes.npz',**save);pd.DataFrame(checks).to_csv(O/'variant-gates.csv',index=False)
result['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest();result['protocol_sha256']=hashlib.sha256((P/'shared43-protocol.md').read_bytes()).hexdigest()
(O/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
