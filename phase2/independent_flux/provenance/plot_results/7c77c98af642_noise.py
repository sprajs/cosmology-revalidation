"""Separate fixed/recomputed noise and normalized/unnormalized objectives."""
from fit import *


class NoiseCase(Case):
    def __init__(self,engine,cid):
        super().__init__(engine,cid)
        f=pd.read_csv(BUNDLE/'exact_covariance_components.csv',dtype={'CID':str})
        f=f[f.CID==cid].sort_values('fit_row_one_based')
        assert len(f)==len(self.t)
        self.mw=f.frozen_MWXT_FLUXERR.values
        self.oldflux=f.frozen_MODELFLUX.values
        self.dc=np.diag(f.DATAFLUX_ERR.values**2+f.FUDGEFLUX_ERR.values**2)
        _,sc=engine.sncosmo_flux(self.p0,self.b,self.t,self.z,self.ebv,True)
        self.native=self.dc+sc+np.outer(self.mw,self.mw)

    def recomputed(self,x,convention):
        p=x.copy();p[3]+=self.tref
        f=self.flux(x)
        if convention=='native':
            _,cov=self.engine.sncosmo_flux(p,self.b,self.t,self.z,self.ebv,True)
        else:
            _,cov,_,_=self.engine.model_covariance(p,self.b,self.t,self.z,self.grid)
            sigma_relative=np.sqrt(np.diag(cov))/abs(f)
            nonlinear=-np.expm1(-sigma_relative)/sigma_relative
            cov=cov*np.outer(nonlinear,nonlinear)
        mw=self.mw*f/self.oldflux
        return cov+self.dc+np.outer(mw,mw)

    def dynamic_fit(self,convention,logdet):
        reference_logdet=np.linalg.slogdet(self.cov)[1]
        def objective(x):
            cv=self.recomputed(x,convention)
            r=self.y-self.flux(x);l=np.linalg.cholesky(cv)
            w=solve_triangular(l,r,lower=True)
            value=w@w+np.sum(self.prior_residuals(x)**2)
            if logdet:value+=2*np.log(np.diag(l)).sum()-reference_logdet
            return .5*value
        bounds=[(self.x0[0]-3,self.x0[0]+3),(-5,5),(-.5,.5),(-5,5)]
        bounds[-1]=(max(-5,float(np.max(self.t-self.tref-(1+self.z)*self.engine.tmax)+1e-6)),min(5,float(np.min(self.t-self.tref-(1+self.z)*self.engine.tmin)-1e-6)))
        # Use explicit central derivatives in dimensionless/time-offset coordinates.
        grad=lambda x:jacobian(lambda a:np.array([objective(a)]),x).ravel()
        res=minimize(objective,self.x0,jac=grad,method='L-BFGS-B',bounds=bounds,options={'ftol':1e-12,'gtol':1e-6,'maxiter':250})
        h=hessian(objective,res.x);eig=np.linalg.eigvalsh(h)
        return {'CID':self.cid,'noise_convention':convention,'includes_logdet':logdet,
                'success':bool(res.success),'message':str(res.message),'objective':float(res.fun),
                'x':res.x.tolist(),'delta_mB_x1_c_t0':((res.x-self.x0)*np.array([-2.5/np.log(10),1,1,1])).tolist(),
                'hessian_min_eigenvalue':float(eig[0]),'hessian_covariance':np.linalg.inv(h).tolist() if eig[0]>0 else None}


def main():
    engine=Engine();_,_,pars=read_inputs();out=[];held=[];components=[];dynamic=[]
    for cid in sorted(pars):
        c=NoiseCase(engine,cid)
        fit,mu,j,hc=c.fit(cov=c.native,two_starts=True)
        fit['delta_mB_x1_c_t0']=(np.array(fit['x'])-c.x0)*np.array([-2.5/np.log(10),1,1,1])
        fit['delta_mB_x1_c_t0']=fit['delta_mB_x1_c_t0'].tolist()
        out.append(fit)
        model=c.cov-c.dc-np.outer(c.mw,c.mw)
        pold=np.array(c.priors['previous_iteration_FITVAL_x0_x1_c_t0']);pold[0]=np.log(pold[0])
        pf,pc,snake,cd=c.engine.model_covariance(pold,c.b,c.t,c.z,c.grid)
        rel=np.sqrt(np.diag(pc))/abs(pf);fac=-np.expm1(-rel)/rel
        nonlinear=pc*np.outer(fac,fac)
        comps={'CID':cid,'model_variance_fraction_total_median':float(np.median(np.diag(model)/np.diag(c.cov))),
               'mw_variance_fraction_total_median':float(np.median(c.mw**2/np.diag(c.cov))),
               'colour_dispersion_fraction_model_variance_median':float(np.median(cd**2/(snake+cd**2))),
               'reconstructed_model_cov_relative_frobenius_error':float(np.linalg.norm(nonlinear-model)/np.linalg.norm(model)),
               'native_to_official_total_sigma_median':float(np.median(np.sqrt(np.diag(c.native)/np.diag(c.cov))))}
        components.append(comps)
        testmask=np.array([int(hashlib.sha256(f'independent-flux-v1|{cid}|{int(np.floor(t+.5))}'.encode()).hexdigest(),16)%5==0 for t in c.t])
        train,test=np.where(~testmask)[0],np.where(testmask)[0]
        for covname,cov in [('official',c.cov),('native',c.native)]:
            for family in ['salt','f99','phase_colour']:
                f,mu,j,hc=c.fit(family,train,cov=cov,two_starts=True)
                f['noise_convention']=covname
                f['prediction']=c.predictive(f,mu,j,hc,train,test,cov=cov)
                held.append(f)
        for convention in ['native','snana_formula']:
            for logdet in [False,True]:dynamic.append(c.dynamic_fit(convention,logdet))
        print(cid,'noise arms complete',flush=True)
    summary={'components':components,'frozen_native_refits':out,'dynamic_refits':dynamic,'heldout':held}
    (OUT/'noise_results.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    quick={'dynamic_successes':sum(q['success'] for q in dynamic),'dynamic_total':len(dynamic),
           'max_native_frozen_abs_delta_mB_x1_c_t0':np.max(abs(np.array([q['delta_mB_x1_c_t0'] for q in out])),axis=0).tolist(),
           'epoch_comparisons':{}}
    indexed={(q['CID'],q['noise_convention'],q['family']):q for q in held}
    for noise,family in [('native','salt'),('native','f99'),('native','phase_colour'),('official','f99'),('official','phase_colour')]:
        diff=[];chis=[];n=0
        for cid in sorted(pars):
            a=indexed[(cid,noise,family)]['prediction'];b=indexed[(cid,'official','salt')]['prediction']
            if a['valid'] and b['valid']:
                diff.append(a['logpdf_laplace']-b['logpdf_laplace']);chis.append(a['test_conditional_chi2']);n+=a['n_test']
        quick['epoch_comparisons'][noise+'_'+family]={'sum_delta_logpdf_vs_official_salt':float(np.sum(diff)),
            'object_cluster_standard_error_sum':float(np.std(diff,ddof=1)*np.sqrt(len(diff))),
            'conditional_chi2_sum':float(np.sum(chis)),'n_test':n,'n_valid_objects':len(diff)}
    salt=[q['prediction'] for q in held if q['noise_convention']=='official' and q['family']=='salt']
    quick['official_salt_chi2']=float(sum(p['test_conditional_chi2'] for p in salt));quick['n_test']=sum(p['n_test'] for p in salt)
    for convention in ['native','snana_formula']:
        for ld in [False,True]:
            r=[q for q in dynamic if q['noise_convention']==convention and q['includes_logdet']==ld]
            quick[f'{convention}_logdet{ld}_max_abs_parameter_shift']=np.max(abs(np.array([q['delta_mB_x1_c_t0'] for q in r])),axis=0).tolist()
    (OUT/'noise_summary.json').write_text(json.dumps(quick,indent=2,allow_nan=False)+'\n')
    provenance('noise',[OUT/'noise_results.json',OUT/'noise_summary.json'],[OUT/'PLAN.json',OUT/'AMENDMENTS.md'])
    print(json.dumps(quick,indent=2),flush=True)


if __name__=='__main__':main()
