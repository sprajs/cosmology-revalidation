"""Conditional recovery/PPC and independent-ensemble Monte Carlo diagnostics."""
import json
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.linalg import cholesky
from scipy.interpolate import PchipInterpolator
from scipy.stats import chi2
import pandas as pd
from core import ROOT, Pantheon, mu, qvalue, manifest

out=ROOT/'runs/cosmology/diagnostics';out.mkdir(parents=True,exist_ok=True)
sn=Pantheon();rng=np.random.default_rng(61020);inputs=sn.inputs.copy();results={}
n=200;truth=.33
noise=rng.normal(size=(n,len(sn.z)))@cholesky(sn.cov,lower=True).T
fake=mu(sn.z,[truth],'lcdm',sn.zhel)[0]+noise
y=fake-sn.fid
B=sn.P;BAB=B.T@sn.A@B;BAy=y@sn.A@B;yAy=np.einsum('bi,ij,bj->b',y,sn.A,y)
est=[];err=[]
for i in range(n):
    def chi(om):
        g=mu(sn.nodes,[om],'lcdm')[0]-sn.fidnodes
        return float(yAy[i]-2*g@BAy[i]+g@BAB@g)
    fit=minimize_scalar(chi,bounds=(.05,.7),method='bounded',options={'xatol':1e-10})
    h=.0001;curv=(chi(fit.x+h)+chi(fit.x-h)-2*chi(fit.x))/h**2
    est.append(float(fit.x));err.append(float(np.sqrt(2/curv)))
est=np.array(est);err=np.array(err)
results['lcdm_conditional_recovery']={'n':n,'true_Om':truth,'estimate_mean':float(est.mean()),
    'bias':float(est.mean()-truth),'mean_mcse':float(est.std(ddof=1)/np.sqrt(n)),
    'empirical_sd':float(est.std(ddof=1)),'mean_curvature_sd':float(err.mean()),
    'one_sigma_coverage':float(np.mean(abs(est-truth)<err)),
    'coverage_binomial_se_at_nominal':float(np.sqrt(.6827*(1-.6827)/n)),
    'qualification':'Simulated Gaussian released-distance likelihood, not survey/photometry/selection recovery.'}
pd.DataFrame({'Om_hat':est,'curvature_sd':err}).to_csv(out/'lcdm_injections.csv',index=False)

results['posterior_predictive']={}
for name in ['pantheon-cpl','pantheon-cpl-c14fixed','pantheon-bao-cpl-c14slope']:
    base=ROOT/'runs/cosmology'/name
    config=json.loads((base/'configuration.json').read_text())
    path=base/'chains.npz';inputs.append(path)
    a=np.load(path,allow_pickle=False);chain=a['chain'].reshape(-1,a['chain'].shape[-1])
    chosen=chain[rng.choice(len(chain),4000,replace=False)]
    amp=0.
    if config['correction']:
        p=ROOT/config['correction'];inputs.append(p);d=pd.read_csv(p)
        sn.reset_template(PchipInterpolator(d.z,d.delta_mu)(sn.z))
        amp=chosen[:,-1] if config['amplitude']=='normal' else 1.
    else:sn.reset_template(np.zeros(len(sn.z)))
    observed=sn.chisq(chosen[:,:3],'cpl',amp)
    replicated=rng.chisquare(len(sn.z)-1,len(chosen))
    results['posterior_predictive'][name]={'replicate_gaussian_df':len(sn.z)-1,
        'observed_chi2_median':float(np.median(observed)),'replicate_chi2_median':float(np.median(replicated)),
        'P_replicate_chi2_at_least_observed':float(np.mean(replicated>=observed)),
        'P_replicate_chi2_at_most_observed':float(np.mean(replicated<=observed)),
        'analytic_P_replicate_chi2_at_most_observed':float(np.mean(chi2.cdf(observed,len(sn.z)-1))),
        'qualification':'Posterior predictive Gaussian distance-product check; lower observed scatter can indicate conservative covariance or non-Gaussian/selected residuals, not validation of all physics.'}

results['replicate_ensembles']={}
for first,second in [('pantheon-bao-cpl','pantheon-bao-cpl-repeat'),
                     ('pantheon-bao-cpl-c14slope','pantheon-bao-cpl-c14slope-repeat')]:
    if not (ROOT/'runs/cosmology'/second/'chains.npz').exists():continue
    pair=[]
    for name in [first,second]:
        p=ROOT/'runs/cosmology'/name/'chains.npz';inputs.append(p)
        a=np.load(p,allow_pickle=False)['chain'];q=.5+1.5*a[:,:,1]*(1-a[:,:,0])
        blocks=np.array_split(q,25,axis=0)
        mean_se=np.std([b.mean() for b in blocks],ddof=1)/np.sqrt(len(blocks))
        prob_se=np.std([(b<0).mean() for b in blocks],ddof=1)/np.sqrt(len(blocks))
        pair.append({'name':name,'q0_mean':float(q.mean()),'q0_sd':float(q.std()),
                     'mean_block_mcse':float(mean_se),'P_q0_negative_sample_fraction':float((q<0).mean()),
                     'probability_block_mcse':float(prob_se),'opposite_sign_draws':int(np.count_nonzero(q>=0)),
                     'zero_tail_notice':'zero block MCSE when all draws share a sign is not a resolved posterior tail'})
    results['replicate_ensembles'][first]={'runs':pair,
        'mean_difference_in_combined_block_mcse':float(abs(pair[0]['q0_mean']-pair[1]['q0_mean'])/
                 np.hypot(pair[0]['mean_block_mcse'],pair[1]['mean_block_mcse']))}
(out/'results.json').write_text(json.dumps(results,indent=2)+'\n')
manifest(out,'Conditional Gaussian recovery, posterior predictions, and independent ensemble comparison',inputs,
         {'seed':61020,'recovery_draws':200,'PPC_draws':4000,'ensemble_blocks':25},
         [out/'results.json',out/'lcdm_injections.csv'])
print(json.dumps(results,indent=2))
