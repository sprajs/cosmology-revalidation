"""Covariance geometry, independent propagation and synthetic recovery tests."""
from noise import *


def run():
    rng=np.random.default_rng(9212026);engine=Engine();_,_,pars=read_inputs()
    full=json.loads((OUT/'independent_refits.json').read_text())
    geom=[];profiles=[];linear=[];nonlinear=[];injection=[];noiseless=[]
    for cid in sorted(pars):
        c=NoiseCase(engine,cid)
        fit=next(r for r in full if r['CID']==cid and r['interpolation']=='sncosmo')
        conv=np.diag([-2.5/np.log(10),1,1,1]);cov=conv@np.array(fit['hessian_covariance'])@conv
        ref=np.array(pars[cid]['hessian_covariance']);tr=np.diag([-2.5/np.log(10)/np.exp(c.x0[0]),1,1,1]);ref=tr@ref@tr
        a=np.linalg.cholesky(ref[:3,:3]);white=np.linalg.solve(a,cov[:3,:3])@np.linalg.inv(a.T)
        corr=lambda v:v/np.sqrt(np.outer(np.diag(v),np.diag(v)))
        geom.append({'CID':cid,'generalized_eigenvalues_mB_x1_c':np.linalg.eigvalsh(white).tolist(),
                     'max_absolute_correlation_difference_3x3':float(np.max(abs(corr(cov[:3,:3])-corr(ref[:3,:3])))),
                     'max_absolute_correlation_difference_4x4':float(np.max(abs(corr(cov)-corr(ref))))})
        if cid in ['1896213','1442085','1318689','1652139']:
            x=np.array(fit['x']);sig=np.sqrt(np.array(fit['hessian_covariance'])[3,3]);l=np.linalg.cholesky(c.cov)
            for q in [-2.,-1.,-.5,0.,.5,1.,2.]:
                dt=x[3]+q*sig
                def residual(a):
                    xx=np.r_[a,dt]
                    return np.r_[solve_triangular(l,c.y-c.flux(xx),lower=True),c.prior_residuals(xx)]
                r=least_squares(residual,x[:3],bounds=([x[0]-3,-5,-.5],[x[0]+3,5,.5]),xtol=1e-11,ftol=1e-11,gtol=1e-8)
                profiles.append({'CID':cid,'delta_t0':q*sig,'independent_sigma':sig,
                                 'snana_hessian_sigma':float(np.sqrt(ref[3,3])),
                                 'delta_chi2_profile':float(r.fun@r.fun-fit['chi2_total']),
                                 'quadratic_prediction':q*q})
        # Local synthetic generator: a Cholesky draw is independent of fitting.
        truth=c.x0.copy();c.priorcentre=truth[3];mu=c.flux(truth)
        nf,_,_,_=c.fit(y=mu,start=truth+np.array([.03,.2,.01,.2]),two_starts=True)
        noiseless.append({'CID':cid,'max_abs_scaled_parameter_error':float(np.max(abs(np.array(nf['x'])-truth)/np.array([.1,1,.1,1]))),
                          'chi2':nf['chi2_total'],'success':nf['success']})
        j=jacobian(c.flux,truth);ci=np.linalg.inv(c.cov);pr=np.diag([0,0,0,.01]);h=j.T@ci@j+pr
        A=np.linalg.solve(h,j.T@ci);actual=A@c.cov@A.T
        noise=rng.standard_normal((10000,len(mu)))@np.linalg.cholesky(c.cov).T
        deviations=noise@A.T;samplecov=np.cov(deviations,rowvar=False)
        linear.append({'CID':cid,'diagonal_variance_ratio_MC_to_analytic':(np.diag(samplecov)/np.diag(actual)).tolist(),
                       'analytic_sampling_vs_inverse_hessian_sigma_ratio':np.sqrt(np.diag(actual)/np.diag(np.linalg.inv(h))).tolist()})
        for rep in range(10):
            y=mu+np.linalg.cholesky(c.cov)@rng.standard_normal(len(mu))
            for name,cv in [('correct_frozen',c.cov),('native_frozen',c.native)]:
                f,_,_,hc=c.fit(y=y,cov=cv)
                entry={'CID':cid,'replicate':rep,'noise_used':name,'success':f['success'],'boundary':f['boundary'],
                       'hessian_positive':hc is not None,'delta':(np.array(f['x'])-truth).tolist()}
                if hc is not None:
                    pull=(np.array(f['x'])-truth)/np.sqrt(np.diag(hc));entry['pulls']=pull.tolist()
                    entry['covers_95']=(abs(pull)<1.95996398454).tolist()
                nonlinear.append(entry)
        testmask=np.array([int(hashlib.sha256(f'independent-flux-v1|{cid}|{int(np.floor(t+.5))}'.encode()).hexdigest(),16)%5==0 for t in c.t])
        train,test=np.where(~testmask)[0],np.where(testmask)[0]
        y=c.flux(np.r_[truth,.06],family='phase_colour')+np.linalg.cholesky(c.cov)@rng.standard_normal(len(mu))
        for family in ['salt','phase_colour']:
            f,ff,jj,hc=c.fit(family,train,y=y)
            f['prediction']=c.predictive(f,ff,jj,hc,train,test,y=y)
            injection.append(f)
        print(cid,'synthetic validation complete',flush=True)
    # Exact scalar-spectral-error transport counterexample, independent of real maps.
    basis=np.array([2.,.5]);bcov=np.array([[.04,.001],[.001,.0025]]);x1=.4;v=np.array([1.,x1]);truthvar=float(v@bcov@v);mean=float(v@basis)
    perturb=rng.multivariate_normal(np.zeros(2),bcov,size=100000)@v
    transport=[]
    for z in [0.,.5,1.]:
        for color in [0.,.1]:
            colorfactor=10**(-.4*color*1.3);scale=colorfactor/(1+z)
            truefracvar=truthvar/mean**2
            mcfracvar=np.var(scale*perturb,ddof=1)/(scale*mean)**2
            transport.append({'z':z,'color':color,'mc_relative_variance':float(mcfracvar),
                              'propagated_relative_variance':truefracvar,
                              'native_effective_wavelength_limit':truefracvar,
                              'snana_denominator_limit':truthvar/(scale*mean)**2})
    out={'geometry':geom,'t0_profiles':profiles,'linearized_covariance':linear,'noiseless':noiseless,
         'nonlinear_synthetic':nonlinear,'phase_colour_injection':injection,'variance_transport_toy':transport}
    (OUT/'validation_results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    summary={'generalized_covariance_eigenvalue_min':min(min(g['generalized_eigenvalues_mB_x1_c']) for g in geom),
             'generalized_covariance_eigenvalue_max':max(max(g['generalized_eigenvalues_mB_x1_c']) for g in geom),
             'max_correlation_difference_3x3':max(g['max_absolute_correlation_difference_3x3'] for g in geom),
             'max_correlation_difference_4x4':max(g['max_absolute_correlation_difference_4x4'] for g in geom),
             'noiseless_max_parameter_error_scaled':max(g['max_abs_scaled_parameter_error'] for g in noiseless),
             'linear_MC_sigma_ratio_range':[float(np.sqrt(min(min(g['diagonal_variance_ratio_MC_to_analytic']) for g in linear))),float(np.sqrt(max(max(g['diagonal_variance_ratio_MC_to_analytic']) for g in linear)))],
             'synthetic':{}}
    for name in ['correct_frozen','native_frozen']:
        allq=[q for q in nonlinear if q['noise_used']==name]
        q=[q for q in allq if q['success'] and not q['boundary'] and q['hessian_positive']]
        pulls=np.array([r['pulls'] for r in q])
        summary['synthetic'][name]={'attempted':len(allq),'interior_success_positive_hessian':len(q),'mean_pulls':pulls.mean(axis=0).tolist(),
                                   'std_pulls':pulls.std(axis=0,ddof=1).tolist(),'coverage95':np.mean(abs(pulls)<1.95996398454,axis=0).tolist()}
    scores=[]
    for cid in pars:
        d={r['family']:r for r in injection if r['CID']==cid}
        if all(r['prediction']['valid'] for r in d.values()):scores.append(d['phase_colour']['prediction']['logpdf_laplace']-d['salt']['prediction']['logpdf_laplace'])
    summary['phase_colour_injection']={'n':len(scores),'sum_delta_logpdf':float(sum(scores)),
                                       'object_cluster_standard_error_sum':float(np.std(scores,ddof=1)*np.sqrt(len(scores)))}
    (OUT/'validation_summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    provenance('validate',[OUT/'validation_results.json',OUT/'validation_summary.json'],[OUT/'PLAN.json',OUT/'AMENDMENTS.md',OUT/'independent_refits.json'])
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':run()
