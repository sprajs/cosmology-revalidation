"""First-party GLS of the full released SH0ES linear system; no H0 prior export."""
import json,importlib.metadata
import numpy as np
from astropy.io import fits
from scipy import linalg,special,stats
from common import ROOT,HERE,WORK,OUT,sha,relative,write

NAMES=[f'host_distance_modulus_{i+1}'for i in range(37)]+[
 'delta_mu_NGC4258','Cepheid_MH_10day_solar','delta_mu_LMC','mu_M31',
 'Cepheid_slope_delta_from_minus3p285','SN_MB','Cepheid_metallicity_slope',
 'released_zero_constraint_auxiliary','ground_HST_zeropoint_offset','five_log10_H0']

def load():
 acquisition=json.loads((OUT/'acquisition.json').read_text())
 for row in acquisition['files']:assert sha(ROOT/row['path'])==row['sha256']
 data={c:fits.getdata(next(WORK.glob('all'+c+'_*.fits'))).astype(np.float64)for c in'ylc'}
 y,A,C=data['y'],data['l'].T,data['c']
 assert y.ndim==1 and A.shape==(len(y),len(NAMES))and C.shape==(len(y),len(y))
 assert all(np.isfinite(x).all()for x in[y,A,C])and np.array_equal(C,C.T)
 return y,A,C

def gls(y,A,C):
 chol=linalg.cholesky(C,lower=True)
 wy=linalg.solve_triangular(chol,y,lower=True);WA=linalg.solve_triangular(chol,A,lower=True)
 fit,_,rank,singular=linalg.lstsq(WA,wy,lapack_driver='gelsd')
 assert rank==A.shape[1],'No silent rank reduction or regularization.'
 Q,R=linalg.qr(WA,mode='economic');qrfit=linalg.solve_triangular(R,Q.T@wy)
 Rinv=linalg.solve_triangular(R,np.eye(len(R)));cov=Rinv@Rinv.T
 normal=WA.T@WA;direct=linalg.solve(normal,WA.T@wy,assume_a='pos')
 directcov=linalg.solve(normal,np.eye(len(normal)),assume_a='pos')
 residual=y-A@fit;white=linalg.solve_triangular(chol,residual,lower=True)
 relative_fit=max(np.max(abs(fit-qrfit)),np.max(abs(fit-direct)))
 cov_error=np.max(abs(cov-directcov))/np.max(abs(cov))
 normal_residual=np.max(abs(WA.T@white))/max(1,np.max(abs(WA.T@wy)))
 assert relative_fit<1e-8 and cov_error<1e-9 and normal_residual<1e-10
 return fit,cov,{'rank':int(rank),'design_condition_whitened':float(singular[0]/singular[-1]),
  'fit_cross_solver_max_absolute':float(relative_fit),'covariance_cross_solver_relative_max':float(cov_error),
  'residual_normal_equation_relative_max':float(normal_residual),'chi2':float(white@white),
  'degrees_of_freedom':len(y)-rank,'logdet_covariance':float(2*np.log(np.diag(chol)).sum())}

def h0_summary(value,var):
 t=np.log(10)/5;mu=t*value;variance=t*t*var
 return {'median_km_s_Mpc':float(np.exp(mu)),'mean_km_s_Mpc':float(np.exp(mu+variance/2)),
  'sd_km_s_Mpc':float(np.sqrt(np.expm1(variance)*np.exp(2*mu+variance))),
  'linear_propagation_sd_km_s_Mpc':float(np.exp(mu)*t*np.sqrt(var)),
  'quantile_probabilities':[.025,.16,.5,.84,.975],
  'equal_tail_quantiles_km_s_Mpc':np.exp(mu+np.sqrt(variance)*stats.norm.ppf([.025,.16,.5,.84,.975])).tolist(),
  'measure':'Gaussian linear-coordinate likelihood with flat q measure; exactly lognormal H0. This is not a new Gaussian H0 likelihood for use with overlapping SN.'}

def main():
 y,A,C=load();fit,cov,validation=gls(y,A,C);sigma=np.sqrt(np.diag(cov))
 reference=np.loadtxt(WORK/'lstsq_results.txt');assert reference.shape==(47,2)
 h0=h0_summary(fit[-1],cov[-1,-1]);ndof=validation['degrees_of_freedom']
 singleton=[]
 for row in np.flatnonzero(np.count_nonzero(A,axis=1)==1):
  col=int(np.flatnonzero(A[row])[0]);singleton.append({'row_zero_based':int(row),'parameter':NAMES[col],
   'coefficient':float(A[row,col]),'released_y':float(y[row]),'released_sigma':float(np.sqrt(C[row,row])),
   'offdiagonal_covariance_nonzero_count':int(np.count_nonzero(C[row])-1)})
 sn=A[:,42]!=0;hf=A[:,46]!=0;cal=sn&~hf
 blocks={'Cepheid_and_external':np.flatnonzero(~sn),'SN_calibrators':np.flatnonzero(cal),'SN_Hubble_flow':np.flatnonzero(hf)}
 blockcov=[]
 for i,(left,il)in enumerate(blocks.items()):
  for right,ir in list(blocks.items())[i+1:]:
   sub=C[np.ix_(il,ir)];den=np.sqrt(np.diag(C)[il,None]*np.diag(C)[None,ir])
   blockcov.append({'left':left,'right':right,'entries':int(sub.size),'nonzero':int(np.count_nonzero(sub)),
    'maximum_absolute_correlation':float(np.max(abs(sub/den)))})
 output=WORK/'gls-result.npz';np.savez_compressed(output,q=fit,covariance=cov,parameter_names=np.array(NAMES),
  row_SN_calibrator=cal,row_SN_HF=hf)
 result={'status':'passed_released_linear_system_reconstruction','rows':len(y),'parameters':len(fit),'H0':h0,
  'validation':validation,'conditional_chi2_per_dof':validation['chi2']/ndof,
  'conditional_chi2_upper_tail':float(stats.chi2.sf(validation['chi2'],ndof)),
  'parameters_by_index':[{'index_zero_based':i,'name':name,'fit':float(fit[i]),'sigma':float(sigma[i]),
   'author_initialization_value':float(reference[i,0]),'author_initialization_sigma':float(reference[i,1]),
   'minus_author_initialization':float(fit[i]-reference[i,0])}for i,name in enumerate(NAMES)],
  'Cepheid_period_slope':{'fit':float(fit[41]-3.285),'sigma':float(sigma[41]),'release_reference_offset':-3.285},
  'published_comparison':{'Riess2022_baseline_H0':73.04,'Riess2022_fit_sigma':1.01,'Riess2022_including_analysis_variant_sigma':1.04,
   'reconstructed_median_minus_baseline':h0['median_km_s_Mpc']-73.04,
   'reconstructed_propagated_sigma_minus_fit_sigma':h0['linear_propagation_sd_km_s_Mpc']-1.01,
   'additional_analysis_variant_systematic_not_added':True},
  'row_blocks':{name:len(ix)for name,ix in blocks.items()},'released_cross_covariance_blocks':blockcov,'singleton_constraints':singleton,
  'source_sha256':{relative(p):sha(p)for p in[HERE/'solve.py',HERE/'common.py',HERE/'design.json']},
  'acquisition_sha256':sha(OUT/'acquisition.json'),'output_sha256':{relative(output):sha(output)},
  'versions':{p:importlib.metadata.version(p)for p in['numpy','scipy','astropy']},
  'limits':['Conditional on released selected/corrected Cepheid and SN products, fixed covariance, standardization and low-redshift expansion; not an upstream photometric reconstruction.',
   'Author initialization values are broad-box reference values, not additional measurements; zero rounded uncertainty at coordinate44 is retained as the actual finite released Gaussian constraint.',
   'No covariance rescaling, data clipping, new anchor prior or cosmological H0 prior is applied.',
   'No independence from Dovekie or Pantheon+ is asserted; their shared events and calibration require a joint covariance or latent response model.']}
 write(OUT/'reconstruction.json',result)
 print(json.dumps({'status':result['status'],'H0':h0,'validation':validation,'published_comparison':result['published_comparison']},indent=2))
if __name__=='__main__':main()
