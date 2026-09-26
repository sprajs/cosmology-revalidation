#!/usr/bin/env python3
"""Independent numerical checks of distances, ages, BAO, populations and dust.

Run from any directory with the clean edition's locked Python environment.
The tested workflow kernels are not imported. Checks start from frozen inputs;
the fresh workflow outputs are comparison targets, never the calculation source.
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
import json
import argparse
import hashlib
import zipfile
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.linalg import cholesky, solve_triangular, cho_solve
from scipy.optimize import minimize_scalar, minimize, LinearConstraint, brentq
from scipy.integrate import quad, solve_ivp, cumulative_trapezoid
from scipy.interpolate import PchipInterpolator
from scipy.stats import chi2
from astropy.cosmology import FlatLambdaCDM

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data'
OUT = ROOT / 'validation/reports'
PREFIX = 'baseline/'
SECOND_SEED = 'robustness/lcdm'


def summary(name):
    return json.loads((ROOT / 'results' / (PREFIX + name) / 'summary.json').read_text())


def load_distances():
    table = pd.read_csv(DATA / 'distances/Pantheon+SH0ES.dat', sep=r'\s+', dtype={'CID': str})
    raw = np.loadtxt(DATA / 'distances/Pantheon+SH0ES_STAT+SYS.cov')
    n = int(raw[0])
    assert n == len(table) and len(raw) == 1+n*n
    cov = raw[1:].reshape(n, n)
    assert np.max(abs(cov-cov.T)) < 1e-7
    return table, (cov+cov.T)/2


def cosmology(table, covariance):
    mask = table.zHD.to_numpy() > .01
    d = table.loc[mask]
    C = covariance[np.ix_(mask, mask)]
    L = cholesky(C, lower=True)
    one = solve_triangular(L, np.ones(len(d)), lower=True)
    y = d.m_b_corr.to_numpy()
    z, zh = d.zHD.to_numpy(), d.zHEL.to_numpy()

    def model(om):
        # Independent Astropy distance implementation and direct full covariance.
        return 5*np.log10((1+zh)*FlatLambdaCDM(H0=70, Om0=om, Tcmb0=0).comoving_distance(z).value)

    def objective(om):
        w = solve_triangular(L, y-model(om), lower=True)
        w -= one * (one @ w)/(one @ one)
        return float(w @ w)

    fit = minimize_scalar(objective, bounds=(.1, .6), method='bounded', options={'xatol':1e-12})
    assert fit.success
    grid = np.linspace(.01, .99, 1401)
    prediction = np.array([model(x) for x in grid])
    white = solve_triangular(L, (y[None, :]-prediction).T, lower=True)
    white -= one[:, None] * (one @ white)[None, :]/(one @ one)
    chis = np.sum(white*white, axis=0)
    posterior = np.exp(-.5*(chis-chis.min()))
    posterior /= np.trapezoid(posterior, grid)
    mean = np.trapezoid(grid*posterior, grid)
    sd = np.sqrt(np.trapezoid((grid-mean)**2*posterior, grid))
    cdf = cumulative_trapezoid(posterior, grid, initial=0)
    ci = np.interp([.025,.5,.975], cdf, grid)
    # Grid halving is an independent integration-resolution check.
    coarse = posterior[::2]/np.trapezoid(posterior[::2],grid[::2])
    coarse_mean = np.trapezoid(grid[::2]*coarse,grid[::2])
    assert abs(coarse_mean-mean)<1e-7
    old = summary('cosmology')
    assert abs(fit.x-old['mode'][0])<2e-6
    assert abs(fit.fun-old['mode_chisq'])<1e-5
    assert abs(mean-old['posterior']['Om']['mean'])<.001
    pd.DataFrame({'Omega_m':grid,'posterior_density':posterior,'chi2':chis}).to_csv(OUT/'lcdm-grid.csv',index=False)
    # Do not infer absolute H0 from an intercept-marginalized SN likelihood.
    result = {'rows':len(d),'unique_objects':int(d.CID.nunique()),'direct_full_covariance_mode':float(fit.x),
        'direct_chi2':float(fit.fun),'posterior_mean':float(mean),'posterior_sd':float(sd),
        'posterior_95_interval':ci[[0,2]].tolist(),'posterior_q0_mean':float(1.5*mean-1),
        'posterior_q0_sd':float(1.5*sd),'grid_halving_mean_difference':float(coarse_mean-mean),
        'direct_vs_workflow_mode_difference':float(fit.x-old['mode'][0]),
        'direct_vs_workflow_chi2_difference':float(fit.fun-old['mode_chisq']),
        'grid_vs_mcmc_mean_difference':float(mean-old['posterior']['Om']['mean']),
        'chi2_per_approximate_dof':float(fit.fun/(len(d)-2)),
        'conditional_gaussian_lower_tail':float(chi2.cdf(fit.fun,len(d)-2)),
        'lower_tail_scope':'Reference for fixed Gaussian covariance and locally linear two-parameter mean; upstream fitted covariance and selection are not regenerated. Not a permission to rescale errors.',
        'published_reference':{'url':'https://arxiv.org/abs/2202.04077','Omega_m':.334,'reported_sd':.018,
            'difference_in_published_sd_units':float((mean-.334)/.018),
            'interpretation':'Descriptive difference, not an independent tension significance; overlapping data and differing release/pipeline details.'}}
    second_path=ROOT/'results'/SECOND_SEED/'summary.json'
    if second_path.exists():
        second=json.loads(second_path.read_text())
        result['second_seed']={'mean':second['posterior']['Om']['mean'],'sd':second['posterior']['Om']['sd'],
            'within_ensemble_diagnostics_pass':second['valid_for_posterior_summary'],
            'mean_minus_grid':second['posterior']['Om']['mean']-float(mean)}
    # Fixed-truth parametric recovery, not prior-based simulation calibration.
    # Refit every realization through the independently constructed full likelihood.
    truth=.33; simulations=200; rng=np.random.default_rng(260926731)
    h=solve_triangular(L,(prediction-model(truth)).T,lower=True)
    h-=one[:,None]*(one@h)[None,:]/(one@one)
    noise=rng.normal(size=(simulations,len(d)))
    noise-=np.outer(noise@one,one)/(one@one)
    simchi=np.sum(noise*noise,axis=1)[:,None]-2*noise@h+np.sum(h*h,axis=0)[None,:]
    density=np.exp(-.5*(simchi-simchi.min(axis=1)[:,None]))
    density/=np.trapezoid(density,grid,axis=1)[:,None]
    simmean=np.trapezoid(density*grid,grid,axis=1)
    simcdf=cumulative_trapezoid(density,grid,axis=1,initial=0)
    intervals=np.array([np.interp([.158655,.841345],cdfrow,grid) for cdfrow in simcdf])
    covered=(intervals[:,0]<=truth)&(truth<=intervals[:,1])
    result['fixed_truth_recovery']={'draws':simulations,'seed':260926731,'true_Omega_m':truth,
        'mean_posterior_mean':float(simmean.mean()),'bias':float(simmean.mean()-truth),
        'bias_monte_carlo_se':float(simmean.std(ddof=1)/np.sqrt(simulations)),
        'coverage_68_percent_interval':float(covered.mean()),'covered_draws':int(covered.sum()),
        'coverage_binomial_se':float(np.sqrt(.68269*(1-.68269)/simulations)),
        'scope':'Conditional Gaussian fixed-covariance recovery; no upstream selection or error-model validation. Fixed truth rather than a draw from the prior.'}
    pd.DataFrame({'posterior_mean':simmean,'lower_68':intervals[:,0],'upper_68':intervals[:,1],'covers_truth':covered}).to_csv(OUT/'lcdm-recovery.csv',index=False)
    return result


def ages(table,covariance):
    rows=[]
    with zipfile.ZipFile(DATA/'ages/supplement.zip') as z:
        for basename,label in [('table1.dat','G11'),('table2.dat','R19')]:
            names=[x for x in z.namelist() if Path(x).name==basename]
            assert len(names)==1
            text=z.read(names[0]).decode()
            for line in text.splitlines():
                fields=[x.strip().rstrip('\\').strip() for x in line.split('&')]
                if len(fields)!=5 or not fields[0].isdigit():continue
                residual=fields[3].split('(')
                rows.append({'CID':str(int(fields[0])),'age_source':label,'age':float(fields[1]),'age_err':float(fields[2]),
                    'author_hr':float(residual[0]),'author_hr_original':float(residual[-1].strip(') ')),'author_hr_err':float(fields[4])})
    original=pd.DataFrame(rows)
    table=table.copy();table['pp_row']=np.arange(len(table))
    merged=original.merge(table[table.IDSURVEY==1],on='CID',how='left',validate='many_to_one',indicator=True)
    matched=merged[merged._merge=='both'].drop_duplicates('CID',keep='first')
    selected=matched[(matched.zHD>.06)&(matched.zHD<.42)].copy()
    expected=pd.read_csv(ROOT/'results'/f'{PREFIX}ages'/'selected.csv',dtype={'CID':str})
    assert list(selected.CID)==list(expected.CID)
    reference=FlatLambdaCDM(H0=70,Om0=.3,Tcmb0=0).distmod(selected.zHD.to_numpy()).value
    reference+=5*np.log10((1+selected.zHEL.to_numpy())/(1+selected.zHD.to_numpy()))
    outcome={'hr_corrected':selected.MU_SH0ES.to_numpy()-reference,
        'hr_no_bias':selected.MU_SH0ES.to_numpy()-reference+selected.biasCor_m_b.to_numpy(),
        'hr_tripp':selected.mB.to_numpy()+.148*selected.x1.to_numpy()-3.112*selected.c.to_numpy()+19.253-reference}
    x=selected.age.to_numpy();e=selected.MU_SH0ES_ERR_DIAG.to_numpy();w=1/e**2
    xc=x-np.average(x,weights=w);den=np.sum(w*xc**2)
    indices=selected.pp_row.to_numpy(dtype=int)
    C=covariance[np.ix_(indices,indices)];L=cholesky(C,lower=True)
    X=np.column_stack([np.ones(len(x)),x]);wx=solve_triangular(L,X,lower=True)
    regressions=[]
    ref=summary('ages')
    for name,y in outcome.items():
        slope=float(np.sum(w*xc*y)/den);se=float(1/np.sqrt(den))
        target=next(j for j in ref['wls'] if j['outcome']==name)
        assert abs(slope-target['slope'])<1e-12 and abs(se-target['slope_se'])<1e-12
        wy=solve_triangular(L,y,lower=True);beta=np.linalg.lstsq(wx,wy,rcond=None)[0]
        v=np.linalg.inv(wx.T@wx)
        # Influence diagnostic, not removal of influential objects from primary fit.
        leave=[]
        for i in range(len(x)):
            keep=np.arange(len(x))!=i;xx=x[keep];ww=w[keep];yy=y[keep];center=xx-np.average(xx,weights=ww)
            leave.append(float(np.sum(ww*center*yy)/np.sum(ww*center**2)))
        regressions.append({'outcome':name,'wls_slope':slope,'wls_se':se,'full_covariance_gls_slope':float(beta[1]),
            'full_covariance_gls_se':float(np.sqrt(v[1,1])),'leave_one_object_slope_range':[min(leave),max(leave)]})
    # Independent joint-normal block construction checks conditional Gaussian algebra.
    latent_errors=[]
    for fit in ref['gaussian_latent_age_mle']:
        mean,logtau,intercept,slope,logscatter=fit['theta'];t2=np.exp(2*logtau);n=len(x)
        joint=np.block([[np.diag(t2+selected.age_err.to_numpy()**2),np.eye(n)*slope*t2],
            [np.eye(n)*slope*t2,C+np.eye(n)*(slope*slope*t2+np.exp(2*logscatter))]])
        chol=cholesky(joint,lower=True);r=np.r_[x-mean,outcome[fit['outcome']]-intercept-slope*mean]
        v=solve_triangular(chol,r,lower=True);nll=float(np.log(np.diag(chol)).sum()+.5*(v@v))
        latent_errors.append({'outcome':fit['outcome'],'joint_minus_conditional_nll':nll-fit['nll'],'boundary':fit['boundary']})
        assert abs(nll-fit['nll'])<1e-8
    author_slopes=[]
    for label,frame in [('full_G11',original[original.age_source=='G11']),('full_R19',original[original.age_source=='R19']),('matched_union',selected)]:
        xx=frame.age.to_numpy();ww=1/frame.author_hr_err.to_numpy()**2;center=xx-np.average(xx,weights=ww)
        denom=np.sum(ww*center**2)
        for column in ['author_hr','author_hr_original']:
            author_slopes.append({'sample':label,'outcome':column,'n':len(frame),'fixed_age_wls_slope':float(np.sum(ww*center*frame[column])/denom),'se':float(denom**-.5)})
    sigmas=np.sqrt(np.diag(C))
    return {'literal_source_rows':len(original),'source_counts':original.age_source.value_counts().to_dict(),
        'unmatched_rows':int((merged._merge!='both').sum()),'selected_unique_objects':len(selected),
        'independent_regressions':regressions,'latent_block_likelihood_checks':latent_errors,
        'bias_term_slope_contribution':float(np.sum(w*xc*selected.biasCor_m_b.to_numpy())/den),
        'author_residual_diagnostics':author_slopes,
        'uncertainty_column_check':{'median_quoted_mu_error':float(np.median(e)),
            'median_full_covariance_diagonal_error':float(np.median(sigmas)),
            'median_ratio':float(np.median(e/sigmas)),
            'scope':'Tabulated plotting errors and supplied covariance diagonal differ in the frozen author release; this is not introduced by the clean extraction. Cosmology uses the full covariance.'},
        'scope':'Diagnostics on selected age summaries; GLS still treats ages as exact. Leave-one-object is influence sensitivity, not selection correction. No original age PDFs or exact published LINMIX reconstruction.'}


def populations():
    h0,om,w0,wa=63.6,.353,-.42,-1.75
    def e_la(la):
        a=np.exp(la)
        return np.sqrt(om*a**-3+(1-om)*a**(-3*(1+w0+wa))*np.exp(-3*wa*(1-a)))
    start=977.7922216807892/h0*2*np.exp(-24)/(3*np.sqrt(om))
    sol=solve_ivp(lambda t,y:[977.7922216807892/h0/e_la(t)],(-16,0),[start],rtol=2e-12,atol=2e-13,dense_output=True)
    assert sol.success
    la=np.linspace(-16,0,60001);times=sol.sol(la)[0];redshift=PchipInterpolator(times,np.expm1(-la),extrapolate=False)
    table=pd.read_csv(ROOT/'results'/f'{PREFIX}populations'/'delay_curves.csv');checks=[]
    for z in [0.,1.,2.5]:
        age=float(sol.sol(-np.log1p(z))[0])
        for name,tp,power,alpha in [('C14_smooth',.3,-1.,20.),('cut300',.3,-1.,None),('W26_cut40',.04,-1.13,None)]:
            low=1e-7 if alpha is not None else tp
            def density(delay):
                birth=float(redshift(age-delay))
                if not np.isfinite(birth):return 0.
                sf=.180/(10**(-.997*(birth-1.243))+10**(.241*(birth-1.243))) if birth<100 else 0.
                phi=(delay/tp)**power/(1+(delay/tp)**(power-alpha)) if alpha is not None else delay**power
                return sf*phi
            high=age*(1-1e-9)
            mass=quad(density,low,high,epsabs=1e-10,epsrel=2e-8,points=[tp] if alpha else None)[0]
            mean=quad(lambda d:d*density(d),low,high,epsabs=1e-10,epsrel=2e-8,points=[tp] if alpha else None)[0]/mass
            median=brentq(lambda v:quad(density,low,v,epsabs=1e-10,epsrel=2e-8,points=[tp] if alpha and v>tp else None)[0]/mass-.5,low,high,xtol=1e-9)
            target=table[(np.isclose(table.z,z))&(table.dtd==name)].iloc[0]
            gap=max(abs(mean-target.mean_delay_gyr),abs(median-target.median_delay_gyr))
            assert gap<3e-5,(z,name,gap)
            checks.append({'z':z,'dtd':name,'independent_mean_gyr':mean,'independent_median_gyr':median,'max_difference_gyr':gap})
    return {'method':'Independent adaptive ODE clock, inverse interpolation and adaptive scalar integration; all three DTDs at z=0,1,2.5. Other redshifts are not independently sampled.','checks':checks,'max_delay_difference_gyr':max(x['max_difference_gyr'] for x in checks)}


def bao():
    table=pd.read_csv(DATA/'bao/desi_gaussian_bao_ALL_GCcomb_mean.txt',sep=r'\s+',comment='#',names=['z','value','kind'])
    raw=np.loadtxt(DATA/'bao/desi_gaussian_bao_ALL_GCcomb_cov.txt');zs=np.sort(table[table.kind=='DH_over_rs'].z.unique())
    ids=[int(table.index[(table.z==z)&(table.kind==kind)][0]) for kind in ['DM_over_rs','DH_over_rs'] for z in zs]
    scale=np.r_[np.ones(len(zs)),1+zs];y=table.value.to_numpy()[ids]*scale;C=raw[np.ix_(ids,ids)]*np.outer(scale,scale)
    L=cholesky(C,lower=True);wy=solve_triangular(L,y,lower=True);t=np.log1p(zs);n=len(zs);constraints=[]
    for i in range(n):
        unit=np.zeros(2*n);unit[n+i]=1;constraints.append(unit)
        row=np.zeros(2*n);row[i]=1;row[n+i]=-(t[i]-(t[i-1] if i else 0))
        if i:row[i-1]=-1
        constraints.append(row)
        if i:
            row=np.zeros(2*n);row[n+i-1]=1;row[n+i]=-1;constraints.append(row)
            row=np.zeros(2*n);row[i-1]=1;row[i]=-1;row[n+i-1]=t[i]-t[i-1];constraints.append(row)
    B=np.array(constraints)@L;B/=np.linalg.norm(B,axis=1)[:,None]
    f=minimize(lambda v:.5*np.sum((v-wy)**2),np.zeros(2*n),jac=lambda v:v-wy,method='SLSQP',constraints=[LinearConstraint(B,0,np.inf)],options={'ftol':1e-12,'maxiter':2000})
    assert f.success and np.min(B@f.x)>-1e-7
    target=summary('bao-shape');assert abs(2*f.fun-target['statistic'])<1e-6
    return {'independent_inequality_projection':float(2*f.fun),'workflow_difference':float(2*f.fun-target['statistic']),
        'anisotropic_rows':2*n,'excluded_isotropic_rows':int((table.kind=='DV_over_rs').sum()),
        'interpretation':'Conservative cone-tail result under flat geometry. No CMB likelihood; no interpretation as a present-day q0 probability.'}


def dust():
    rows=pd.read_csv(ROOT/'results'/f'{PREFIX}dust'/'geometry.csv')
    # Independent integral over uniformly distributed source depth.
    errors=[]
    for tau in [0.,1e-10,.01,.1,1.,10.]:
        numeric=quad(lambda d:np.exp(-tau*d),0,1,epsabs=1e-13)[0]
        analytic=1 if tau==0 else -np.expm1(-tau)/tau
        errors.append(abs(numeric-analytic))
    assert max(errors)<1e-12
    return {'quadrature_max_transmission_error':max(errors),'stored_geometry_rows':len(rows),
        'interpretation':'Constructed absorption and latent-mean identity, not an empirical dust-population estimate.'}


def main():
    global PREFIX, SECOND_SEED
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', default=PREFIX.rstrip('/'))
    parser.add_argument('--second-seed', default=SECOND_SEED)
    args=parser.parse_args()
    PREFIX=args.results.rstrip('/')+'/'
    SECOND_SEED=args.second_seed
    OUT.mkdir(parents=True,exist_ok=True)
    table,cov=load_distances()
    result={'cosmology':cosmology(table,cov),'ages':ages(table,cov),'populations':populations(),'bao':bao(),'dust':dust()}
    result['audit_code_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result['input_sha256']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [DATA/'distances/Pantheon+SH0ES.dat',DATA/'distances/Pantheon+SH0ES_STAT+SYS.cov',DATA/'ages/supplement.zip',DATA/'bao/desi_gaussian_bao_ALL_GCcomb_mean.txt',DATA/'bao/desi_gaussian_bao_ALL_GCcomb_cov.txt']}
    (OUT/'core.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'passed','groups':list(result)[:5]},indent=2))

if __name__=='__main__':main()
