"""Synthetic transport counterexamples, not author-data reproductions."""
from common import ROOT, OUT, manifest
import json
import numpy as np
import pandas as pd

SEED, N, S = 20260920, 200_000, -.03

def slope(x, y):
    return np.mean((x-x.mean())*(y-y.mean()))/np.var(x)

def evaluate(label, h0, h1, t0, t1, y0, y1, extra=None):
    bh, bt = slope(h0,y0), slope(t0,y0)
    truth = y1.mean()-y0.mean()
    pred_h, pred_t = bh*(h1.mean()-h0.mean()), bt*(t1.mean()-t0.mean())
    # Monte Carlo precision for the whole pipeline from 20 independent blocks.
    err_blocks=[]
    for ih, il in zip(np.array_split(np.arange(len(h0)),20), np.array_split(np.arange(len(h1)),20)):
        p=slope(h0[ih], y0[ih])*(h1[il].mean()-h0[ih].mean())
        err_blocks.append(p-(y1[il].mean()-y0[ih].mean()))
    return dict(case=label, n0=len(h0),n1=len(h1), host_slope=bh, progenitor_slope=bt,
                observed_mean_shift=truth, slope_times_supplied_delay_mean_shift=S*(t1.mean()-t0.mean()),
                predicted_host_product=pred_h, predicted_progenitor_product=pred_t,
                host_prediction_error=pred_h-truth, progenitor_prediction_error=pred_t-truth,
                host_error_mc_se=np.std(err_blocks,ddof=1)/np.sqrt(20), **(extra or {}))

def main():
    rng=np.random.default_rng(SEED)
    h0=rng.uniform(2,10,N); h1=rng.uniform(1,5,N)
    t0=.5+.3*h0; t1=.5+.3*h1
    rows=[evaluate('affine_deterministic',h0,h1,t0,t1,S*t0,S*t1)]
    assert abs(rows[-1]['host_prediction_error']) < 1e-14
    assert abs(rows[-1]['predicted_host_product']-rows[-1]['predicted_progenitor_product']) < 1e-14

    # Classical host measurement error: fixed affine re-expression preserves the wrong prediction.
    hn0=h0+rng.normal(0,1.5,N); hn1=h1+rng.normal(0,1.5,N)
    y0,y1=S*t0,S*t1
    row=evaluate('noisy_host_classical',hn0,hn1,t0,t1,y0,y1)
    row['analytic_host_slope']=S*.3*((8**2/12)/(8**2/12+1.5**2))
    tn0=.5+.3*hn0; tn1=.5+.3*hn1
    row['same_noisy_proxy_remap_product']=slope(tn0,y0)*(tn1.mean()-tn0.mean())
    assert abs(row['same_noisy_proxy_remap_product']-row['predicted_host_product']) < 1e-14
    # Correct errors-in-variables denominator when known variance is available.
    row['known_error_variance_corrected_product']=slope(hn0,y0)*np.var(hn0)/(np.var(hn0)-1.5**2)*(hn1.mean()-hn0.mean())
    rows.append(row)

    # Broad progenitor distribution is not by itself failure of a stable conditional-mean map.
    u0=rng.uniform(-.4,.4,N);u1=rng.uniform(-.4,.4,N)
    ts0=t0+u0;ts1=t1+u1
    rows.append(evaluate('latent_scatter_stable_conditional_mean',h0,h1,ts0,ts1,S*ts0,S*ts1))
    # Replacing unknown delay by an independent draw adds variance but no covariance with observed Y.
    tr0=t0+rng.uniform(-.4,.4,N); tr1=t1+rng.uniform(-.4,.4,N)
    row=evaluate('independent_delay_imputation',h0,h1,tr0,tr1,S*ts0,S*ts1)
    row['analytic_imputed_slope']=S*(.3**2*(8**2/12))/(.3**2*(8**2/12)+.8**2/12)
    rows.append(row)

    # Nonlinear known map: fit on actual transformed quantity works, one host slope need not transport.
    tq0=.1*h0**2;tq1=.1*h1**2
    rows.append(evaluate('nonlinear_true_delay_response',h0,h1,tq0,tq1,S*tq0,S*tq1))
    # Exact reparameterization of a host-linear response is nonlinear in transformed age.
    row=evaluate('nonlinear_reparameterization_host_response',h0,h1,tq0,tq1,S*h0,S*h1)
    row['exact_function_pushforward_shift']=S*(np.sqrt(10*tq1).mean()-np.sqrt(10*tq0).mean())
    assert abs(row['exact_function_pushforward_shift']-row['observed_mean_shift']) < 1e-14
    rows.append(row)

    # A conditional offset can change even at the same host-age value.
    tz1=t1+.8
    row=evaluate('population_dependent_mapping',h0,h1,t0,tz1,S*t0,S*tz1)
    row['analytic_host_prediction_error']=-S*.8
    rows.append(row)

    # Survey inclusion changes the conditional residual distribution as well as the age mix.
    e0=rng.normal(0,.1,N);e1=rng.normal(0,.1,N)
    sy0=S*t0+e0;sy1=S*t1+e1; keep=sy1<-.05
    rows.append(evaluate('brightness_selected_high_z',h0,h1[keep],t0,t1[keep],sy0,sy1[keep]))

    df=pd.DataFrame(rows)
    p=OUT/'controlled-mapping.csv';df.to_csv(p,index=False)
    mixture=dict(label='synthetic_two_point_mixture', age_support_gyr=[1,9],
        low_z_weights=[.6,.4],high_z_weights=[.9,.1],
        low_z_mean_gyr=4.2,high_z_mean_gyr=1.8,low_z_median_gyr=1,high_z_median_gyr=1,
        mean_luminosity_shift_mag=S*(1.8-4.2),median_product_mag=0)
    q=OUT/'mean-median-counterexample.json';q.write_text(json.dumps(mixture,indent=2)+'\n')
    rconf=np.random.default_rng(20260921)
    u=rconf.normal(size=N)
    ha=5+1.5*u+rconf.normal(0,.2,N)
    mass=10+.5*u+rconf.normal(0,.2,N)
    lum=-.05*u+rconf.normal(0,.05,N)
    hi=mass>=10
    def step(y): return y[hi].mean()-y[~hi].mean()
    age_corrected=lum-slope(ha,lum)*(ha-ha.mean())
    mass_corrected=lum-step(lum)*(hi.astype(float)-.5)
    asym=dict(data_type='synthetic common-driver model, no age->luminosity arrow',seed=20260921,n=N,
        age_slope_before=slope(ha,lum),age_slope_after_mass_step=slope(ha,mass_corrected),
        mass_step_before=step(lum),mass_step_after_age_correction=step(age_corrected))
    a=OUT/'causal-asymmetry-counterexample.json';a.write_text(json.dumps(asym,indent=2)+'\n')
    manifest('controlled-mapping',dict(seed=SEED,n_per_population=N,coefficient_mag_gyr=S,
        data_type='synthetic mathematical controlled experiments; no observational claims'),
        [__file__,ROOT/'scripts/mapping/common.py',ROOT/'docs/experiments/mapping-plan.md',ROOT/'uv.lock'],[p,q,a])
    print(df[['case','host_prediction_error','progenitor_prediction_error','host_error_mc_se']].to_string(index=False))
    print(json.dumps(asym,indent=2))

if __name__=='__main__':main()
