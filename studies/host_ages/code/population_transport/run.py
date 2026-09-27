#!/usr/bin/env python3
"""Observed transport support and conditional host-to-progenitor sensitivity.

Run from anywhere with the repository's locked Python environment. Numerical
arrays are regenerated in .work; only compact scientific records are retained.
No likelihood for true host age or survey selection is inferred by this script.
"""
from pathlib import Path
import datetime
import hashlib
import json
import sys

import numpy as np
import pandas as pd
from astropy.cosmology import Flatw0waCDM
from scipy.integrate import cumulative_trapezoid, quad
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import brentq, minimize_scalar
from scipy.stats import norm, ks_2samp

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
from lib.ages import read_age
from lib.populations import COSMOS, DTDS, csfh, clock
from lib.cosmology import Pantheon, mu

HERE = Path(__file__).parent
OUT = ROOT / 'studies/host_ages/results/population_transport'
WORK = ROOT / '.work/population-transport'
CFG = json.loads((HERE / 'scenarios.json').read_text())


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def jwrite(p, obj):
    p.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')


def weighted_mean_range(x, ratio):
    """Exact extrema over independent density weights with max/min <= ratio."""
    x = np.sort(np.asarray(x))
    if ratio == 1:
        return [float(x.mean()), float(x.mean())]
    candidates = []
    for k in range(len(x) + 1):
        w = np.where(np.arange(len(x)) < k, ratio, 1)
        candidates.append(float(w @ x / w.sum()))
        candidates.append(float((ratio + 1 - w) @ x / (ratio + 1 - w).sum()))
    return [min(candidates), max(candidates)]


def entropy_weights(x, target):
    if not min(x) < target < max(x):
        return {'feasible': False, 'reason': 'target outside convex support'}
    def weights(lam):
        a = lam * (x - x.mean())
        w = np.exp(a - a.max())
        return w / w.sum()
    lam = brentq(lambda a: weights(a) @ x - target, -100, 100)
    w = weights(lam)
    return dict(feasible=True, target_mean_gyr=float(target),
                achieved_mean_gyr=float(w @ x), effective_size=float(1 / (w @ w)),
                max_weight=float(w.max()), exponential_tilt_per_gyr=float(lam))


def empirical(pp, full_cov):
    frames = [read_age(ROOT / f'data/ages/table{i}.dat', name)
              for i, name in [(1, 'G11'), (2, 'R19')]]
    ages = pd.concat(frames, ignore_index=True)
    joined = ages.merge(pp[pp.IDSURVEY == 1], on='CID', validate='many_to_one')
    rng = np.random.default_rng(CFG['seed'])
    all_results = {}
    retained = {}
    for policy in CFG['age_catalogues']:
        x = joined.drop_duplicates('CID', keep='first' if policy == 'G11_first' else 'last')
        x = x[(x.zHD > .06) & (x.zHD < .42)].copy()
        assert len(x) == 196 and x.CID.nunique() == 196
        x.to_csv(WORK / f'{policy}-observed-hosts.csv', index=False)
        a = x.age.to_numpy()
        low = x.zHD.to_numpy() < .2
        hi = ~low
        contrast = float(a[low].mean() - a[hi].mean())
        boot = []
        for _ in range(CFG['bootstrap_replicates']):
            al = rng.choice(a[low], low.sum(), replace=True)
            ah = rng.choice(a[hi], hi.sum(), replace=True)
            boot.append([al.mean()-ah.mean(), np.median(al)-np.median(ah)])
        boot = np.asarray(boot)
        bins = CFG['age_histogram_edges_gyr']
        hl, hh = [np.histogram(a[mask], bins)[0]/mask.sum() for mask in [low, hi]]
        central = [max(np.quantile(a[low], .05), np.quantile(a[hi], .05)),
                   min(np.quantile(a[low], .95), np.quantile(a[hi], .95))]
        counts = []
        for mask in [low, hi]:
            counts.append(dict(n=int(mask.sum()), z_min=float(x.zHD.to_numpy()[mask].min()),
                z_max=float(x.zHD.to_numpy()[mask].max()), age_mean_gyr=float(a[mask].mean()),
                age_median_gyr=float(np.median(a[mask])), age_sd_gyr=float(a[mask].std(ddof=1)),
                age_error_median_gyr=float(np.median(x.age_err.to_numpy()[mask])),
                age_range_gyr=[float(a[mask].min()),float(a[mask].max())],
                outside_shared_central_support=int(((a[mask]<central[0]) | (a[mask]>central[1])).sum())))
        selection = []
        for ratio in CFG['selection_density_ratio_limits']:
            lrange, hrange = weighted_mean_range(a[low], ratio), weighted_mean_range(a[hi], ratio)
            delta = [lrange[0]-hrange[1], lrange[1]-hrange[0]]
            selection.append(dict(max_to_min_density_weight_ratio=ratio,
                low_mean_range_gyr=lrange, high_mean_range_gyr=hrange,
                low_minus_high_contrast_range_gyr=delta,
                conditional_high_minus_low_brightness_mag=[.03*v for v in delta]))
        ids = x.pp_row.to_numpy(dtype=int)
        c = full_cov[np.ix_(ids, ids)]
        cf = cho_factor(c)
        one = np.ones(len(x))
        nuisance = {'intercept': np.column_stack([one]),
            'intercept_z': np.column_stack([one,x.zHD]),
            'intercept_z_mass_colour_width': np.column_stack([one,x.zHD,x.HOST_LOGMASS,x.c,x.x1])}
        info = []
        for name, design in nuisance.items():
            icx = cho_solve(cf, design)
            fitted = design @ np.linalg.solve(design.T @ icx, icx.T @ a)
            residual = a - fitted
            information = float(residual @ cho_solve(cf, residual))
            info.append(dict(nuisance=name, known_age_information_gyr2_per_mag2=information,
                optimistic_slope_se_mag_per_gyr=float(1/np.sqrt(information)),
                two_sided_5percent_power_for_minus_0p03=float(norm.cdf(-norm.ppf(.975)-.03*np.sqrt(information))+norm.sf(norm.ppf(.975)-.03*np.sqrt(information))),
                fraction_information_relative_intercept=None))
        for d in info:
            d['fraction_information_relative_intercept']=d['known_age_information_gyr2_per_mag2']/info[0]['known_age_information_gyr2_per_mag2']
        all_results[policy] = dict(sample_bins=counts,
            age_mean_low_minus_high_gyr=contrast,
            age_mean_contrast_bootstrap_95=np.quantile(boot[:,0],[.025,.975]).tolist(),
            age_median_low_minus_high_gyr=float(np.median(a[low])-np.median(a[hi])),
            age_median_contrast_bootstrap_95=np.quantile(boot[:,1],[.025,.975]).tolist(),
            bootstrap_scope='Physical SN resampling of published summary values only; does not draw latent ages or propagate host-inference priors.',
            age_histogram_overlap=float(np.minimum(hl,hh).sum()),
            shared_central_age_support_gyr=[float(v) for v in central],
            age_ks_statistic=float(ks_2samp(a[low],a[hi]).statistic),
            age_ks_pvalue=float(ks_2samp(a[low],a[hi]).pvalue),
            match_high_observed_mean_by_weighting_low=entropy_weights(a[low],a[hi].mean()),
            match_assumed_3p1Gyr_mean_by_weighting_low=entropy_weights(a[low],3.1),
            conditional_selection_sensitivity=selection,
            known_age_design_information=info,
            naive_classical_error_variance_gyr2=float(np.var(a,ddof=1)-np.mean(x.age_err.to_numpy()**2)))
        retained[policy] = x
    coverage=[]
    physical = pp.drop_duplicates('CID')
    known = set(retained['G11_first'].CID)
    for lo,hi in zip([.01,.06,.2,.42,1], [.06,.2,.42,1,2.5]):
        sel=physical[(physical.zHD>lo)&(physical.zHD<=hi)]
        coverage.append(dict(z_interval=[lo,hi], pantheon_unique_CID=int(len(sel)),
            age_study_matched_CID=int(sel.CID.isin(known).sum())))
    return dict(catalogues=all_results,coverage_of_existing_age_study_in_pantheon=coverage),retained


def mapping_tests(frame):
    age=frame.age.to_numpy();z=frame.zHD.to_numpy();low=z<.2;high=~low
    b=CFG['reference_host_slope_mag_per_gyr']
    mean=.2+.6*age
    agef,_=clock(COSMOS['lcdm_H70'])
    upper=agef(-np.log1p(z))
    assert np.all((mean>.04)&(mean<upper))
    width=np.minimum.reduce([.45*mean,.9*(mean-.04),.9*(upper-mean)])
    tau=mean[:,None]+width[:,None]*np.array([-1.,1.])
    gw=np.array([.5,.5])
    assert np.all((tau>=.04)&(tau<=upper[:,None]))
    cases={
        'affine': (mean, np.ones(len(age))),
        'quadratic': (.1+.15*age+.06*age**2,np.ones(len(age))),
        'redshift_dependent': (.2+.6*age-1.2*(z-.12),np.ones(len(age))),
        'mean_preserving_scatter': (tau@gw,np.ones(len(age))),
    }
    weights=np.exp(-tau/4)
    selection=weights@gw
    cases['scatter_with_delay_selection']=((tau*weights)@gw/selection,selection)
    answer=[]
    olddelta=age[low].mean()-age[high].mean()
    for name,(mapped,w) in cases.items():
        al=age[low]; tl=mapped[low]; wl=w[low]/w[low].sum()
        cov=wl@((al-wl@al)*(tl-wl@tl)); var=wl@((al-wl@al)**2)
        k=float(cov/var)
        slope=b/k
        tlo=float(w[low]@mapped[low]/w[low].sum())
        thi=float(w[high]@mapped[high]/w[high].sum())
        delta=tlo-thi
        high_minus_low=slope*(thi-tlo)
        answer.append(dict(mapping=name, low_sample_regression_mapping_gyr_per_gyr=k,
            jointly_calibrated_progenitor_slope=slope, low_minus_high_mean_delay_gyr=delta,
            predicted_high_minus_low_brightness_mag=high_minus_low,
            original_host_linear_prediction_mag=-b*olddelta,
            departure_from_simple_compensation_mag=high_minus_low+b*olddelta,
            low_observed_host_slope_reconstructed_mag_per_gyr=slope*k,
            scope='Constructed map on measured hosts, not inferred progenitor ages; last case additionally reweights hosts under hypothetical delay selection.'))
    assert abs(answer[0]['departure_from_simple_compensation_mag'])<1e-14
    assert abs(answer[3]['departure_from_simple_compensation_mag'])<1e-14
    assert max(abs(d['low_observed_host_slope_reconstructed_mag_per_gyr']-b) for d in answer)<1e-14
    return answer


def distribution(c,sfh,dtd,z,n,cn,tilt=0):
    agef,zf=clock(c,cn)
    cosmic=float(agef(-np.log1p(z)))
    tp,s,alpha=DTDS[dtd]
    lo=1e-7 if alpha is not None else tp
    delay=np.geomspace(lo,cosmic*(1-1e-9),n)
    form=np.nan_to_num(csfh(zf(cosmic-delay),sfh))
    x=np.log(delay/tp)
    phi=delay**s if alpha is None else np.exp(alpha*x-np.logaddexp((alpha-s)*x,0))
    density=form*phi*np.exp(tilt*delay)
    cumulative=cumulative_trapezoid(density,delay,initial=0)
    area=cumulative[-1]
    return dict(mean=float(np.trapezoid(delay*density,delay)/area),
                median=float(np.interp(area/2,cumulative,delay)),
                second_moment=float(np.trapezoid(delay**2*density,delay)/area),
                cosmic_age=cosmic)


def cosmic_scenarios():
    records=[]
    for cname in CFG['cosmologies']:
        c=COSMOS[cname]
        for sfh in CFG['sfhs']:
            for dtd in CFG['dtds']:
                for tilt in CFG['selection_log_weight_per_gyr']:
                    series=[distribution(c,sfh,dtd,z,CFG['delay_points'],CFG['clock_points'],tilt) for z in CFG['population_z']]
                    zero=series[0]
                    for z,d in zip(CFG['population_z'],series):
                        records.append(dict(cosmology=cname,sfh=sfh,dtd=dtd,selection_tilt_per_gyr=tilt,z=z,
                            mean_delay_gyr=d['mean'],median_delay_gyr=d['median'],cosmic_age_gyr=d['cosmic_age'],
                            conditional_linear_mean_brightness_drift_mag=.03*(zero['mean']-d['mean']),
                            conditional_linear_median_template_mag=.03*(zero['median']-d['median']),
                            nonlinear_mean_minus_function_of_median_mag=(-.03*d['mean']+.002*d['second_moment'])-(-.03*d['median']+.002*d['median']**2)))
    pd.DataFrame(records).to_csv(WORK/'cosmic-scenarios.csv',index=False)
    # Every DTD/clock receives a finer-grid check, with a second quadrature scheme.
    errors=[]
    for cname in CFG['cosmologies']:
        c=COSMOS[cname]
        _,zf=clock(c,2*CFG['clock_points']-1)
        agef,_=clock(c,CFG['clock_points'])
        for sfh in CFG['sfhs']:
            for dtd in CFG['dtds']:
                for z in [0.,1.,2.]:
                    coarse=distribution(c,sfh,dtd,z,CFG['delay_points'],CFG['clock_points'])
                    fine=distribution(c,sfh,dtd,z,2*CFG['delay_points']-1,2*CFG['clock_points']-1)
                    t=fine['cosmic_age'];tp,s,alpha=DTDS[dtd]
                    lo=1e-7 if alpha is not None else tp
                    def f(v,power):
                        delay=np.exp(v); xx=np.log(delay/tp)
                        phi=delay**s if alpha is None else np.exp(alpha*xx-np.logaddexp((alpha-s)*xx,0))
                        return float(np.nan_to_num(csfh(zf(t-delay),sfh))*phi*delay**(power+1))
                    area=quad(lambda v:f(v,0),np.log(lo),np.log(t*(1-1e-9)),epsabs=1e-9,epsrel=1e-8,limit=200)[0]
                    mean=quad(lambda v:f(v,1),np.log(lo),np.log(t*(1-1e-9)),epsabs=1e-9,epsrel=1e-8,limit=200)[0]/area
                    errors.append(dict(cosmology=cname,sfh=sfh,dtd=dtd,z=z,
                        grid_mean_difference_gyr=abs(fine['mean']-coarse['mean']),
                        grid_median_difference_gyr=abs(fine['median']-coarse['median']),
                        independent_log_quad_mean_difference_gyr=abs(fine['mean']-mean),
                        independent_astropy_clock_difference_gyr=abs(float(agef(-np.log1p(z)))-Flatw0waCDM(**c,Tcmb0=0).age(z).value)))
    assert max(x['grid_median_difference_gyr'] for x in errors)<1e-4
    assert max(x['independent_log_quad_mean_difference_gyr'] for x in errors)<1e-4
    assert max(x['independent_astropy_clock_difference_gyr'] for x in errors)<2e-6
    return dict(at_z1=[x for x in records if x['z']==1], integration_checks=errors,
        scope='Parent cosmic SFH convolution and arbitrary exp(tilt*delay) selection stress, not survey selection likelihood or observed host reconstruction. Fixed .030 coefficient compares imposed delay models; it is not a refitted progenitor slope.')


def precision():
    p=Pantheon()
    om=.33225848
    h=1e-4
    derivative=(mu(p.z,[om+h],'lcdm',p.zhel)[0]-mu(p.z,[om-h],'lcdm',p.zhel)[0])/(2*h)
    shape=np.minimum(p.z,1)
    response=float(derivative@p.A@shape/(derivative@p.A@derivative))
    sigma_om=float(1/np.sqrt(derivative@p.A@derivative))
    baseline=json.loads((ROOT/'results/robustness/lcdm/summary.json').read_text())
    reference_sigma=float(baseline['q0']['sd'])
    target_q=CFG['tolerance_fraction_of_baseline_q0_sigma']*reference_sigma
    a=target_q/abs(1.5*response)
    true=mu(p.z,[om],'lcdm',p.zhel)[0]
    outcomes=[]
    for amp in [-a,a,.01]:
        y=true+amp*shape
        def chi(v):
            r=y-mu(p.z,[v],'lcdm',p.zhel)[0]
            return float(r@p.A@r)
        fit=minimize_scalar(chi,bounds=(.15,.55),method='bounded',options={'xatol':1e-12})
        outcomes.append(dict(bias_amplitude_mag=amp, recovered_omega_m=float(fit.x),
            delta_q0=float(1.5*(fit.x-om)),linear_delta_q0=1.5*response*amp))
    return dict(shape=CFG['bias_shape'], reference_q0_sigma=reference_sigma,
        target_max_abs_delta_q0=target_q, covariance_fisher_sigma_q0=1.5*sigma_om,
        linear_q0_response_per_mag=1.5*response,
        shape_specific_amplitude_tolerance_mag=a, independent_nonlinear_distance_fits=outcomes,
        scope='Known-covariance flat-LCDM design response around Omega_m=.33225848, with free magnitude intercept. No cosmological measurement is changed, and this does not give a shape-independent magnitude tolerance.')


def main():
    OUT.mkdir(parents=True,exist_ok=True);WORK.mkdir(parents=True,exist_ok=True)
    inputs=[ROOT/'data/ages/table1.dat',ROOT/'data/ages/table2.dat',ROOT/'data/ages/supplement.zip',
        ROOT/'data/distances/Pantheon+SH0ES.dat',ROOT/'data/distances/Pantheon+SH0ES_STAT+SYS.cov',
        ROOT/'results/robustness/lcdm/summary.json']
    frozen={d['path']:d['sha256'] for d in json.loads((ROOT/'provenance/inputs.json').read_text())['files']}
    for p in inputs[:5]:
        assert sha(p)==frozen[str(p.relative_to(ROOT))],f'Frozen input mismatch: {p}'
    pp=pd.read_csv(inputs[3],sep=r'\s+',dtype={'CID':str});pp['pp_row']=np.arange(len(pp))
    raw=np.loadtxt(inputs[4]);n=int(raw[0]);assert n==len(pp)
    covariance=raw[1:].reshape(n,n);covariance=(covariance+covariance.T)/2
    empirical_result,frames=empirical(pp,covariance)
    results=dict(experiment='Host-to-progenitor mapping and selected-population transport',
        empirical=empirical_result, mapping_scenarios={k:mapping_tests(v) for k,v in frames.items()},
        cosmic_population=cosmic_scenarios(),precision=precision(),
        conclusion='Observed summaries permit limited within-sample transport diagnostics. The full selected high-redshift age distribution and calibrated host-to-progenitor map remain unidentified; no empirical physical bound on an additional correction is established.')
    jwrite(OUT/'summary.json',results)
    code=[Path(__file__),HERE/'scenarios.json',ROOT/'lib/ages.py',ROOT/'lib/populations.py',ROOT/'lib/cosmology.py']
    record=dict(completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        command='OPENBLAS_NUM_THREADS=1 .venv/bin/python studies/host_ages/code/population_transport/run.py',
        input_sha256={str(p.relative_to(ROOT)):sha(p) for p in inputs},
        code_sha256={str(p.relative_to(ROOT)):sha(p) for p in code},
        output_sha256={str(p.relative_to(ROOT)):sha(p) for p in [OUT/'summary.json',*WORK.glob('*.csv')]},
        numerical_scope='New conditional experiment; not exact author reproduction. Historical results unchanged.')
    jwrite(OUT/'run.json',record)
    print(json.dumps(dict(observed=empirical_result,precision=results['precision'],mapping=results['mapping_scenarios']['G11_first']),indent=2))

if __name__=='__main__':
    main()
