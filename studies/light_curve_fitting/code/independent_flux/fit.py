"""Independent optimizer, frozen flux likelihood and conditional predictions."""
from engine import *
from scipy.linalg import solve_triangular
from scipy.optimize import least_squares, minimize


def jacobian(fun,x):
    h=np.array([2e-5,2e-4,2e-5,2e-3]+([2e-5] if len(x)==5 else []))
    return np.column_stack([(fun(x+np.eye(len(x))[j]*h[j])-fun(x-np.eye(len(x))[j]*h[j]))/(2*h[j]) for j in range(len(x))])


def hessian(fun,x,scale=1.):
    h=np.array([2e-4,2e-3,2e-4,1e-2]+([2e-4] if len(x)==5 else []))*scale
    f=fun(x); n=len(x); out=np.empty((n,n)); eye=np.eye(n)
    for i in range(n):
        out[i,i]=(fun(x+h[i]*eye[i])-2*f+fun(x-h[i]*eye[i]))/h[i]**2
        for j in range(i):
            out[i,j]=out[j,i]=(fun(x+h[i]*eye[i]+h[j]*eye[j])-fun(x+h[i]*eye[i]-h[j]*eye[j])-fun(x-h[i]*eye[i]+h[j]*eye[j])+fun(x-h[i]*eye[i]-h[j]*eye[j]))/(4*h[i]*h[j])
    return out


def normal_logpdf(residual,cov):
    l=np.linalg.cholesky(cov)
    white=solve_triangular(l,residual,lower=True)
    return float(-.5*(white@white+2*np.log(np.diag(l)).sum()+len(residual)*np.log(2*np.pi)))


class Case:
    def __init__(self,engine,cid):
        self.engine,self.cid=engine,cid
        self.q={k:v for k,v in np.load(BUNDLE/f'objective_{cid}.npz').items()}
        self.z=float(self.q['zHEL'][0]);self.ebv=float(self.q['MWEBV'][0])
        self.p0=self.q['parameters_x0_x1_c_t0'].copy();self.p0[0]=np.log(self.p0[0])
        self.tref=self.p0[3];self.x0=self.p0.copy();self.x0[3]=0.
        self.t=self.q['MJD'];self.b=self.q['band'];self.y=self.q['data_flux']
        self.grid=engine.prepare(self.z,self.ebv,5.)
        self.priors=json.loads((BUNDLE/'exact_prior_setup.json').read_text())['objects'][cid]
        self.priorcentre=self.priors['actual_t0_prior_centre']-self.tref
        self.cov=self.q['frozen_flux_covariance']

    def flux(self,x,family='salt',interpolation='sncosmo'):
        p=x.copy();p[3]+=self.tref
        return self.engine.flux(p,self.b,self.t,self.z,self.grid,family,interpolation)

    def prior_residuals(self,x):
        out=[(x[3]-self.priorcentre)/10.]
        if len(x)==5:out.append(x[4]/.1)
        return np.array(out)

    def fit(self,family='salt',train=None,cov=None,y=None,start=None,interpolation='sncosmo',two_starts=False):
        n=len(self.y);train=np.arange(n) if train is None else np.asarray(train)
        cov=self.cov if cov is None else cov; y=self.y if y is None else y
        ctrain=cov[np.ix_(train,train)];l=np.linalg.cholesky(ctrain)
        pinit=np.r_[self.x0,0.] if family=='phase_colour' else self.x0.copy()
        if start is not None:pinit=np.array(start)
        low=np.array([self.x0[0]-3.,-5.,-.5,-10.]+([-.3] if family=='phase_colour' else []))
        high=np.array([self.x0[0]+3.,5.,.5,10.]+([.3] if family=='phase_colour' else []))
        # Preserve model domain during all optimizer evaluations.
        low[3]=max(low[3],float(np.max(self.t-self.tref-(1+self.z)*self.engine.tmax)+1e-6))
        high[3]=min(high[3],float(np.min(self.t-self.tref-(1+self.z)*self.engine.tmin)-1e-6))
        def residual(x):
            r=solve_triangular(l,y[train]-self.flux(x,family,interpolation)[train],lower=True)
            return np.r_[r,self.prior_residuals(x)]
        starts=[np.clip(pinit,low+1e-8,high-1e-8)]
        if two_starts:starts.append(np.clip(pinit+np.array([.05,.15,-.01,.15]+([.02] if family=='phase_colour' else [])),low+1e-8,high-1e-8))
        sols=[least_squares(residual,s,bounds=(low,high),diff_step=1e-4,xtol=1e-10,ftol=1e-10,gtol=1e-8,max_nfev=250) for s in starts]
        sol=min(sols,key=lambda s:s.fun@s.fun);x=sol.x
        objective=lambda a:.5*np.sum(residual(a)**2)
        h=hessian(objective,x);h2=hessian(objective,x,.5)
        eig=np.linalg.eigvalsh(h);hc=None
        if eig[0]>0:hc=np.linalg.inv(h)
        boundary=bool(np.any(np.minimum(x-low,high-x)<1e-4))
        j=jacobian(lambda a:self.flux(a,family,interpolation),x)
        out={'CID':self.cid,'family':family,'ntrain':len(train),'success':bool(sol.success),'message':sol.message,
             'nfev':int(sol.nfev),'x':x.tolist(),'chi2_total':float(sol.fun@sol.fun),
             'chi2_data':float(sol.fun[:-len(self.prior_residuals(x))]@sol.fun[:-len(self.prior_residuals(x))]),
             'hessian_min_eigenvalue':float(eig[0]),'boundary':boundary,
             'hessian_stepsize_relative_difference':float(np.linalg.norm(h-h2)/np.linalg.norm(h)),
             'hessian_covariance':None if hc is None else hc.tolist(),
             'gauss_newton_covariance':np.linalg.inv(sol.jac.T@sol.jac).tolist(),
             'two_start_chi2_gap':float(max(s.fun@s.fun for s in sols)-min(s.fun@s.fun for s in sols))}
        return out,self.flux(x,family,interpolation),j,hc

    def predictive(self,fit,mu,j,hc,train,test,cov=None,y=None):
        cov=self.cov if cov is None else cov;y=self.y if y is None else y
        train,test=np.asarray(train),np.asarray(test)
        if hc is None or fit['boundary'] or not fit['success']:
            return {'valid':False}
        a=np.linalg.solve(cov[np.ix_(train,train)],cov[np.ix_(train,test)]).T
        condmean=mu[test]+a@(y[train]-mu[train])
        noise=cov[np.ix_(test,test)]-a@cov[np.ix_(train,test)]
        jp=j[test]-a@j[train]
        predictive_cov=noise+jp@hc@jp.T
        r=y[test]-condmean
        return {'valid':True,'n_test':len(test),'logpdf_plugin':normal_logpdf(r,noise),
                'logpdf_laplace':normal_logpdf(r,predictive_cov),
                'test_conditional_chi2':float(r@np.linalg.solve(predictive_cov,r)),
                'test_marginal_standardized_residuals':(r/np.sqrt(np.diag(predictive_cov))).tolist(),
                'test_model':condmean.tolist(),'test_variance':np.diag(predictive_cov).tolist()}


def run():
    engine=Engine();_,_,original=read_inputs(); results=[];full=[];arrays={}
    for cid in sorted(original):
        c=Case(engine,cid)
        for interp in ['sncosmo','scipy']:
            f,mu,j,hc=c.fit(interpolation=interp,two_starts=True)
            f['interpolation']=interp
            transform=np.diag([-2.5/np.log(10),1,1,1])
            delta=(np.array(f['x'])-c.x0)@transform
            ref=np.array(original[cid]['hessian_covariance'])
            tr=np.diag([-2.5/np.log(10)/np.exp(c.x0[0]),1,1,1]);ref=tr@ref@tr
            f['delta_mB_x1_c_t0_from_SNANA']=delta.tolist()
            f['delta_in_SNANA_hessian_sigma']=(delta/np.sqrt(np.diag(ref))).tolist()
            if hc is not None:
                own=transform@hc@transform
                f['hessian_sigma_ratio_to_SNANA']=(np.sqrt(np.diag(own)/np.diag(ref))).tolist()
                f['hessian_covariance_relative_difference']=float(np.linalg.norm(own-ref)/np.linalg.norm(ref))
            full.append(f)
            arrays[cid+'_'+interp+'_prediction']=mu
        groups=[('epoch',np.array([int(hashlib.sha256(f'independent-flux-v1|{cid}|{int(np.floor(t+.5))}'.encode()).hexdigest(),16)%5==0 for t in c.t]))]
        for band in np.unique(c.b):
            if len(np.unique(c.b[c.b!=band]))>=2:groups.append(('band_'+band,c.b==band))
        for split,testmask in groups:
            test,train=np.where(testmask)[0],np.where(~testmask)[0]
            for family in ['salt','f99','phase_colour']:
                f,mu,j,hc=c.fit(family,train,two_starts=True)
                f['split']=split;f['train_rows']=train.tolist();f['test_rows']=test.tolist()
                f['prediction']=c.predictive(f,mu,j,hc,train,test)
                results.append(f)
        print(cid,'completed',len(results),'heldout fits',flush=True)
    (OUT/'independent_refits.json').write_text(json.dumps(full,indent=2,allow_nan=False)+'\n')
    (OUT/'heldout_predictions.json').write_text(json.dumps(results,indent=2,allow_nan=False)+'\n')
    np.savez_compressed(OUT/'independent_refit_predictions.npz',**arrays)
    summarize(full,results)
    provenance('fit',[OUT/x for x in ['independent_refits.json','heldout_predictions.json','independent_refit_predictions.npz','fit_summary.json']],[OUT/'PLAN.json',OUT/'AMENDMENTS.md'])


def summarize(full,results):
    summary={'refits':{},'heldout':{}}
    for interp in ['sncosmo','scipy']:
        q=[r for r in full if r['interpolation']==interp]
        summary['refits'][interp]={'max_abs_delta_mB_x1_c_t0':np.max(abs(np.array([r['delta_mB_x1_c_t0_from_SNANA'] for r in q])),axis=0).tolist(),
                                   'max_abs_delta_in_SNANA_hessian_sigma':np.max(abs(np.array([r['delta_in_SNANA_hessian_sigma'] for r in q])),axis=0).tolist(),
                                   'hessian_sigma_ratio_range':[np.min(np.array([r['hessian_sigma_ratio_to_SNANA'] for r in q]),axis=0).tolist(),np.max(np.array([r['hessian_sigma_ratio_to_SNANA'] for r in q]),axis=0).tolist()],
                                   'failures':sum(not r['success'] or r['boundary'] or r['hessian_covariance'] is None for r in q),
                                   'max_two_start_chi2_gap':max(r['two_start_chi2_gap'] for r in q),
                                   'max_hessian_stepsize_difference':max(r['hessian_stepsize_relative_difference'] for r in q)}
    for splitkind in ['epoch','band']:
        q=[r for r in results if r['split'].startswith(splitkind)]
        groups={}
        for r in q:groups.setdefault((r['CID'],r['split']),{})[r['family']]=r
        per=[]
        for key,d in groups.items():
            if not all(r['prediction']['valid'] for r in d.values()):continue
            row={'CID':key[0],'split':key[1],'n_test':d['salt']['prediction']['n_test']}
            for fam in ['f99','phase_colour']:
                for score in ['logpdf_laplace','logpdf_plugin']:
                    row[fam+'_'+score+'_difference']=d[fam]['prediction'][score]-d['salt']['prediction'][score]
            per.append(row)
        frame=pd.DataFrame(per)
        item={'paired_valid_folds':len(per),'total_folds':len(groups),'folds':per}
        for fam in ['f99','phase_colour']:
            col=fam+'_logpdf_laplace_difference'
            totals=frame.groupby('CID')[col].sum()
            item[fam]={'sum_delta_logpdf':float(totals.sum()),'n_objects':len(totals),
                       'object_cluster_standard_error_sum':float(totals.std(ddof=1)*np.sqrt(len(totals))),
                       'positive_objects':int(sum(totals>0)),
                       'worst_object':str(totals.idxmin()),'worst_delta':float(totals.min()),
                       'best_object':str(totals.idxmax()),'best_delta':float(totals.max())}
        summary['heldout'][splitkind]=item
    (OUT/'fit_summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in summary['refits'].items()},indent=2),flush=True)
    print(json.dumps({k:{x:y for x,y in v.items() if x!='folds'} for k,v in summary['heldout'].items()},indent=2),flush=True)


if __name__=='__main__':run()
