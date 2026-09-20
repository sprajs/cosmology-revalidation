"""Run preregistered corrected-distance cosmology experiments."""
import argparse
import json
import os
from pathlib import Path
import sys
import time

os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import numpy as np
import emcee
from scipy.optimize import minimize
from scipy.interpolate import PchipInterpolator
from core import ROOT, Pantheon, BAO, mu, qvalue, integral, efunc, manifest, sha


def summarize(a):
    a=np.asarray(a)
    return dict(mean=float(a.mean()),sd=float(a.std()),
                q025=float(np.quantile(a,.025)),q16=float(np.quantile(a,.16)),
                median=float(np.median(a)),q84=float(np.quantile(a,.84)),
                q975=float(np.quantile(a,.975)))


def fit(args):
    out=ROOT/'runs/cosmology'/args.name; out.mkdir(parents=True,exist_ok=True)
    sn=Pantheon(magnitude_table=ROOT/args.magnitude_table if args.magnitude_table else None) if args.data in ('sn','joint') else None
    bao=BAO() if args.data in ('bao','joint') else None
    config=vars(args).copy()
    inputs=[]
    if sn: inputs+=sn.inputs
    if bao: inputs+=bao.inputs
    if args.correction:
        import pandas as pd
        p=ROOT/args.correction
        df=pd.read_csv(p)
        f=PchipInterpolator(df[args.zcolumn],df[args.column],extrapolate=False)
        correction=f(sn.z)*args.scale
        assert np.isfinite(correction).all(), 'Correction does not cover SN redshifts'
        sn.reset_template(correction)
        inputs.append(p)
    assert args.data=='sn' or args.model=='cpl'
    base={'lcdm':([.33],['Om'],[(.01,.99)]),
          'wcdm':([.3,-.9],['Om','w'],[(.01,.99),(-3.,1.)]),
          'cpl':([.3,-.8,-.5],['Om','w0','wa'],[(.01,.99),(-3.,1.),(-3.,2.)]),
          'kinematic':([-.4,.8],['q0','q1'],[(-3.,2.),(-8.,8.)]),
          'qbins':([-.4,-.4,-.2,.1,.3],[f'q{i}' for i in range(5)],[(-3.,2.)]*5)}
    start,names,bounds=base[args.model]; ncos=len(start)
    if bao:
        start=start+[10000.]; names=names+['H0_rd']; bounds=bounds+[(5000.,15000.)]
    if args.amplitude in ('normal','uniform'):
        start=start+[1.]; names=names+['evolution_amplitude']; bounds=bounds+[(-2.,4.)]
    bounds=np.asarray(bounds)
    config['actual_priors']={n:list(b) for n,b in zip(names,bounds)}
    config['early_DE_prior']='w0+wa<0' if args.model=='cpl' else None
    config['normal_amplitude_prior']=[1.,4/30] if args.amplitude=='normal' else None
    config['qbin_smoothing_sd']=args.tau if args.model=='qbins' else None
    (out/'configuration.json').write_text(json.dumps(config,indent=2)+'\n')

    def chisq(t):
        t=np.atleast_2d(t); ans=np.zeros(len(t))
        amp=t[:,-1] if args.amplitude in ('normal','uniform') else (1. if args.amplitude=='fixed' else 0.)
        if sn: ans+=sn.chisq(t[:,:ncos],args.model,amp)
        if bao: ans+=bao.chisq(t[:,:ncos],t[:,ncos])
        return ans

    def logpost(t):
        t=np.atleast_2d(t); lp=np.full(len(t),-np.inf)
        ok=np.all((t>bounds[:,0])&(t<bounds[:,1]),axis=1)
        if args.model=='cpl': ok &= t[:,1]+t[:,2]<0
        if np.any(ok):
            v=t[ok]; val=-.5*chisq(v)
            if args.amplitude=='normal': val-=.5*((v[:,-1]-1)/(4/30))**2
            if args.model=='qbins': val-=.5*np.sum(np.diff(v[:,:ncos],axis=1)**2,axis=1)/args.tau**2
            lp[ok]=val
        return lp

    rng=np.random.default_rng(args.seed)
    np.random.seed(args.seed)
    opts=[]
    for k in range(4):
        x=np.array(start)+(rng.normal(size=len(start))*.02*(bounds[:,1]-bounds[:,0]) if k else 0)
        x=np.clip(x,bounds[:,0]+1e-5,bounds[:,1]-1e-5)
        r=minimize(lambda x:-logpost(x)[0],x,method='Nelder-Mead',
                   options={'maxiter':15000,'xatol':1e-7,'fatol':1e-7})
        opts.append(dict(x=r.x.tolist(),objective=float(r.fun),success=bool(r.success),message=r.message))
    best=min(opts,key=lambda r:r['objective']); x=np.array(best['x'])
    walkers=[]
    while len(walkers)<args.walkers:
        v=x+rng.normal(size=len(x))*.005*(bounds[:,1]-bounds[:,0])
        if np.isfinite(logpost(v)[0]): walkers.append(v)
    sampler=emcee.EnsembleSampler(args.walkers,len(x),logpost,vectorize=True)
    t0=time.time()
    sampler.run_mcmc(np.array(walkers),args.steps,progress=False)
    chain=sampler.get_chain(discard=args.burn)
    flat=chain.reshape(-1,len(x))
    tau=sampler.get_autocorr_time(discard=args.burn,tol=0)
    q=qvalue([0.],flat[:,:ncos],args.model)[:,0]
    # Monte Carlo error: block means across retained time, preserving walker interaction.
    qb=q.reshape(chain.shape[:2])
    chunks=np.array_split((qb<0).mean(axis=1),20)
    mcse=float(np.std([c.mean() for c in chunks],ddof=1)/np.sqrt(len(chunks)))
    summary=dict(parameters={n:summarize(flat[:,i]) for i,n in enumerate(names)},
        q0=summarize(q),P_q0_negative=float(np.mean(q<0)),P_q0_negative_block_mcse=mcse,
        label_q0=('q averaged over 0<z<0.1; not point q(0)' if args.model=='qbins' else 'model-extrapolated q(0)'),
        optimizer_starts=opts,posterior_mode_chisq=float(chisq(x)[0]),
        acceptance_fraction=float(sampler.acceptance_fraction.mean()),autocorrelation_time=tau.tolist(),
        ESS_per_parameter=(len(flat)/tau).tolist(),retained_steps_per_tau=((args.steps-args.burn)/tau).tolist(),
        runtime_seconds=time.time()-t0,
        prior_boundary_fraction={n:float(np.mean((flat[:,i]-bounds[i,0]<.01*(bounds[i,1]-bounds[i,0]))|
                                               (bounds[i,1]-flat[:,i]<.01*(bounds[i,1]-bounds[i,0])))) for i,n in enumerate(names)})
    if sn:
        amp=x[-1] if args.amplitude in ('normal','uniform') else (1. if args.amplitude=='fixed' else 0.)
        summary['mode_sn_chisq']=float(sn.chisq(x[:ncos],args.model,amp)[0])
        summary['mode_full_minus_compressed_chisq']=float(sn.chisq_full(x[:ncos],args.model,amp)[0]-summary['mode_sn_chisq'])
        summary['data_diagnostics']=sn.diagnostics()
    if bao: summary['mode_bao_chisq']=float(bao.chisq(x[:ncos],x[ncos])[0])
    z=np.linspace(0,2.3,231)
    select=rng.choice(len(flat),min(5000,len(flat)),replace=False)
    qgrid=qvalue(z,flat[select,:ncos],args.model)
    import pandas as pd
    pd.DataFrame({'z':z,'q025':np.quantile(qgrid,.025,axis=0),'q16':np.quantile(qgrid,.16,axis=0),
                  'median':np.median(qgrid,axis=0),'q84':np.quantile(qgrid,.84,axis=0),
                  'q975':np.quantile(qgrid,.975,axis=0),'P_q_negative':(qgrid<0).mean(axis=0)}).to_csv(out/'q_history.csv',index=False)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    np.savez_compressed(out/'chains.npz',chain=chain,log_prob=sampler.get_log_prob(discard=args.burn),names=names)
    manifest(out,'Conditional released-distance likelihood fit; not raw-photometry reconstruction',inputs,config,
             [out/f for f in ['summary.json','configuration.json','q_history.csv','chains.npz']])
    print(json.dumps({'name':args.name,'q0':summary['q0'],'P_q0_negative':summary['P_q0_negative'],
                      'min_retained_tau':min(summary['retained_steps_per_tau']),'seconds':summary['runtime_seconds']}),flush=True)


def validate():
    from astropy.cosmology import Flatw0waCDM
    from scipy.integrate import quad
    out=ROOT/'runs/cosmology/validation'; out.mkdir(parents=True,exist_ok=True)
    sn=Pantheon(); rng=np.random.default_rng(94314)
    z=np.array([.01,.1,.8,2.3])
    tests={}
    for label,t in [('lcdm',[.3,-1,0]),('cpl',[.353,-.42,-1.75]),('extreme',[.02,-2.9,1.9])]:
        cos=Flatw0waCDM(H0=70,Om0=t[0],w0=t[1],wa=t[2],Tcmb0=0)
        ref=cos.distmod(z).value-5*np.log10(299792.458/70)-25
        tests['astropy_'+label+'_max_mag_error']=float(np.max(np.abs(mu(z,t)[0]-ref)))
        qnum=(1+z)*(efunc(z+1e-5,t)[0]-efunc(z-1e-5,t)[0])/(2e-5)/efunc(z,t)[0]-1
        tests['q_derivative_'+label+'_error']=float(np.max(abs(qnum-qvalue(z,t)[0])))
    for q in [-1.,0.,.5]:
        direct=integral(z,[q]*5,'qbins')[0]
        ref=np.log1p(z) if q==0 else ((1+z)**(-q)-1)/(-q)
        tests['analytic_constant_q_'+str(q)]=float(np.max(abs(direct-ref)))
    for model,theta in [('cpl',np.column_stack([rng.uniform(.01,.99,100),rng.uniform(-3,1,100),rng.uniform(-3,2,100)])),
                        ('qbins',rng.uniform(-3,2,(100,5))),('kinematic',rng.uniform([-3,-8],[2,8],(100,2)))]:
        approx=sn.chisq(theta,model); full=sn.chisq_full(theta,model)
        resid=mu(sn.z,theta,model,sn.zhel)-sn.fid-(mu(sn.nodes,theta,model)-sn.fidnodes)@sn.P.T
        tests[model+'_compression_max_mag_error']=float(abs(resid).max())
        tests[model+'_compression_max_chisq_error']=float(abs(approx-full).max())
    # The source's tiny nonsymmetry: direct general inverse versus symmetric likelihood.
    r=np.loadtxt(sn.inputs[1]); n=int(r[0]); c=r[1:].reshape(n,n)[np.ix_(sn.original_indices,sn.original_indices)]
    inv=np.linalg.inv(c); u=inv@np.ones(len(c)); v=np.ones(len(c))@inv
    A=inv-np.outer(u,v)/v.sum()
    residual=sn.mag-mu(sn.z,[.3316],'lcdm',sn.zhel)[0]
    tests['source_general_inverse_minus_symmetric_chisq']=float(residual@A@residual-sn.chisq_full([.3316],'lcdm')[0])
    tests['offset_invariance_chisq_error']=float((residual+7.13)@sn.A@(residual+7.13)-residual@sn.A@residual)
    # Free dimming function exactly exchanges coasting and LCDM distances.
    m1=mu(z,[.3],'lcdm')[0]; m2=mu(z,[0.]*5,'qbins')[0]; evolution=m1-m2
    tests['exact_evolution_degeneracy_mag_error']=float(abs(m1-(m2+evolution)).max())
    tests['LCDM_minus_flat_coasting_magnitude_at_z']=dict(zip(map(str,z),evolution.tolist()))
    tests['data']=sn.diagnostics()
    copies=[]
    for p in (ROOT/'sources/repos').glob('PantheonPlusSH0ES__DataRelease*/Pantheon+_Data/4_DISTANCES_AND_COVAR/*'):
        if p.name in {x.name for x in sn.inputs}:
            copies.append({'path':str(p.relative_to(ROOT)),'sha256':sha(p)})
    tests['full_release_copy_hashes']=copies
    tests['compact_input_hashes']={str(p.relative_to(ROOT)):sha(p) for p in sn.inputs}
    (out/'results.json').write_text(json.dumps(tests,indent=2)+'\n')
    manifest(out,'Independent analytic and Astropy checks plus full covariance/compression audit',sn.inputs,
             {'seed':94314,'draws_per_model':100,'Gauss_Legendre_order':48,'nodes_per_q_interval':40},[out/'results.json'])
    print(json.dumps(tests,indent=2))
    assert max(tests[k] for k in tests if k.startswith('astropy_'))<1e-7
    assert max(tests[k] for k in tests if 'compression_max_mag' in k)<1e-5
    assert max(tests[k] for k in tests if 'compression_max_chisq' in k)<.01


if __name__=='__main__':
    p=argparse.ArgumentParser(); sub=p.add_subparsers(dest='action',required=True)
    sub.add_parser('validate')
    f=sub.add_parser('fit')
    f.add_argument('--name',required=True); f.add_argument('--model',choices=['lcdm','wcdm','cpl','kinematic','qbins'],default='cpl')
    f.add_argument('--data',choices=['sn','bao','joint'],default='sn')
    f.add_argument('--amplitude',choices=['none','fixed','normal','uniform'],default='none')
    f.add_argument('--correction'); f.add_argument('--column',default='delta_mu'); f.add_argument('--zcolumn',default='z')
    f.add_argument('--magnitude-table',help='Explicitly revised m_b_corr table; identifiers/redshifts must match original row order')
    f.add_argument('--scale',type=float,default=1.); f.add_argument('--tau',type=float,default=.5)
    f.add_argument('--steps',type=int,default=7000); f.add_argument('--burn',type=int,default=1500)
    f.add_argument('--walkers',type=int,default=40); f.add_argument('--seed',type=int,default=2092026)
    a=p.parse_args()
    if a.action=='validate': validate()
    else: fit(a)
