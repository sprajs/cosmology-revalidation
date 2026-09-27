#!/usr/bin/env python3
"""Held-out mapping diagnostics of an older public TITAN author snapshot.

All ages here are posterior summaries from one author's SED/SFH/DTD model.
Held-out prediction of those summaries does not independently validate physics.
"""
from pathlib import Path
import datetime
import hashlib
import json
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT=Path(__file__).resolve().parents[4]
HERE=Path(__file__).parent
OUT=ROOT/'studies/host_ages/results/population_transport'
WORK=ROOT/'.work/population-transport'
SOURCE=WORK/'titan-author'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def write(p,data):
    p.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')


def main():
    registry=json.loads((OUT/'author-inputs.json').read_text())
    for row in registry['files']:
        assert sha(ROOT/row['path'])==row['sha256'],row['path']
    plan=json.loads((HERE/'titan-scenarios.json').read_text())
    d=pd.read_csv(SOURCE/'host_props_with_SN_age_good_Mar18.csv')
    raw=pd.read_csv(SOURCE/'host_props.csv')
    reliable=pd.read_csv(SOURCE/'all_parameters_reliable.csv')
    ddlr=pd.read_csv(SOURCE/'titan_ddlr.csv')
    merged=raw.merge(ddlr,left_on='transient',right_on='objname',validate='one_to_one')
    recovered=merged[merged.transient.isin(reliable.objname)&merged.d_dlr.le(10)]
    assert set(recovered.transient)==set(d.transient)
    assert len(d)==8610 and d.transient.nunique()==len(d)
    common=[x for x in raw if x!='transient']
    x=d.set_index('transient')[common].sort_index()
    y=recovered.set_index('transient')[common].sort_index()
    source_column_max_difference=float(np.nanmax(np.abs(x.to_numpy()-y.to_numpy())))
    assert source_column_max_difference<1e-13
    age=d.mass_weighted_age_50.to_numpy()
    delay=d.prog_age_50.to_numpy()
    covariates=d[['stellar_mass_50','ssfr_50','dust:Av_50','iyer2019:metallicity_50']].to_numpy()
    assert np.isfinite(covariates).all() and np.isfinite(age).all() and np.isfinite(delay).all()
    test=np.array([int(hashlib.sha256(s.encode()).hexdigest()[:8],16)%5==0 for s in d.transient])
    train=~test
    mean=age[train].mean();scale=age[train].std();a=(age-mean)/scale
    cm=covariates[train].mean(axis=0);cs=covariates[train].std(axis=0)
    c=(covariates-cm)/cs
    designs={'affine_host_age':np.column_stack([np.ones(len(d)),a]),
        'cubic_host_age':np.column_stack([np.ones(len(d)),a,a*a,a*a*a]),
        'cubic_host_age_plus_mass_ssfr_dust_metallicity':np.column_stack([np.ones(len(d)),a,a*a,a*a*a,c])}
    fits=[];predictions={}
    for name,design in designs.items():
        coeff=np.linalg.lstsq(design[train],delay[train],rcond=None)[0]
        pred=design@coeff
        predictions[name]=pred
        residual=delay[test]-pred[test]
        fits.append(dict(model=name,train_objects=int(train.sum()),heldout_objects=int(test.sum()),
            coefficients_on_standardized_predictors=coeff.tolist(), heldout_rmse_gyr=float(np.sqrt(np.mean(residual**2))),
            heldout_mae_gyr=float(np.mean(np.abs(residual))),
            heldout_r_squared=float(1-np.sum(residual**2)/np.sum((delay[test]-delay[train].mean())**2))))
    rng=np.random.default_rng(plan['seed'])
    indices=np.flatnonzero(test)
    contrasts=[]
    for name in list(designs)[1:]:
        original=(delay[test]-predictions['affine_host_age'][test])**2
        alternative=(delay[test]-predictions[name][test])**2
        boot=[]
        for _ in range(2000):
            sample=rng.integers(len(indices),size=len(indices))
            boot.append(np.sqrt(original[sample].mean())-np.sqrt(alternative[sample].mean()))
        contrasts.append(dict(comparison='affine RMSE minus '+name+' RMSE',
            reduction_gyr=float(np.sqrt(original.mean())-np.sqrt(alternative.mean())),
            paired_host_bootstrap_95_gyr=np.quantile(boot,[.025,.975]).tolist(),
            uncertainty_scope='Resampling held-out physical-host posterior centres; shared SED model, calibration and joint posterior uncertainty are not included.'))
    uv=d.uv_colour_50.to_numpy();vj=d.vj_colour_50.to_numpy()
    quiescent=(vj<1.6)&(uv>1.3)&(uv>.88*vj+.69)
    starforming=~(uv>.88*vj+.69)
    low=d.stellar_mass_50.to_numpy()<=10
    groups={ 'low_mass_quiescent':low&quiescent,'low_mass_starforming':low&starforming,
             'high_mass_quiescent':(~low)&quiescent,'high_mass_starforming':(~low)&starforming,
             'outside_both_notebook_UVJ_groups':~(quiescent|starforming)}
    groupstats=[]
    for name,mask in groups.items():
        coeff=np.polyfit(age[mask],delay[mask],1)
        groupstats.append(dict(group=name,n=int(mask.sum()),mean_host_age_gyr=float(age[mask].mean()),
            mean_posterior_centre_expected_delay_gyr=float(delay[mask].mean()),
            host_to_delay_linear_slope=float(coeff[0]),
            affine_prediction_mean_residual_gyr=float(np.mean(delay[mask]-predictions['affine_host_age'][mask]))))
    pp=pd.read_csv(ROOT/'data/distances/Pantheon+SH0ES.dat',sep=r'\s+',dtype={'CID':str})
    pp['pp_row']=np.arange(len(pp))
    matched=d.merge(pp,left_on='transient',right_on='CID',validate='one_to_many')
    matched.to_csv(WORK/'titan-pantheon-literal-crosswalk.csv',index=False)
    unique=matched.drop_duplicates('transient')
    above=unique[unique.zHD>.01]
    finite_order=[]
    for prefix in ['mass_weighted_age','prog_age','dust:Av','stellar_mass']:
        finite_order.append(dict(parameter=prefix,disordered_quantiles=int(((d[prefix+'_16']>d[prefix+'_50'])|(d[prefix+'_50']>d[prefix+'_84'])).sum())))
    pd.DataFrame(dict(transient=d.transient,heldout=test,host_age=age,progenitor_expected_delay_summary=delay,
        **{f'prediction_{k}':v for k,v in predictions.items()})).to_csv(WORK/'titan-mapping-predictions.csv',index=False)
    result=dict(source_commit=registry['commit'],source_status=registry['commit_snapshot_status'],
        rows=int(len(d)),unique_transients=int(d.transient.nunique()),
        public_notebook_sample_reconstruction='Exact equality of names; all 57 original host-property columns agree to CSV floating-point round-trip precision after reliable-ID and d_dlr<=10 filter.',
        source_column_max_abs_difference=source_column_max_difference,
        published_final_sample_rows=6983,final_sample_reproduced=False,quantile_order_checks=finite_order,
        age_semantics='prog_age_50 is the posterior median of the per-SFH mean progenitor delay. It is not the median of the progenitor-delay PDF.',
        host_age_mean_gyr=float(age.mean()),expected_delay_posterior_centres_mean_gyr=float(delay.mean()),
        expected_delay_posterior_centres_median_gyr=float(np.median(delay)),
        mean_delay_16_to_84_posterior_width_gyr=float(np.mean(d.prog_age_84-d.prog_age_16)),
        spearman_host_age_vs_delay=float(spearmanr(age,delay).statistic),
        standardization=dict(host_age_train_mean=mean,host_age_train_sd=scale,covariate_train_mean=cm.tolist(),covariate_train_sd=cs.tolist()),
        heldout_fits=fits,paired_comparisons=contrasts,subpopulations=groupstats,
        literal_distance_linkage=dict(pantheon_rows=int(len(matched)),physical_names=int(unique.transient.nunique()),
            above_z0p01_objects=int(len(above)),z_range=[float(unique.zHD.min()),float(unique.zHD.max())],
            names=sorted(unique.transient.tolist()), high_z_observed_ages_added=0,
            limitation='No sky coordinates or redshifts in this author host-summary table; no fuzzy/sky matches attempted. Pantheon repeats are not independent ages. These few very low-z matches do not identify high-z transport.'),
        unresolved_inputs=registry['missing_for_full_reproduction'],
        inference_limits=['Posterior centres share SFH/dust/metallicity/clock assumptions, with no joint likelihood or priors supplied here.',
            'Train/test split tests predictive shape of published model-derived centres, not independent stellar/progenitor ages. Host identifiers are absent, so distinct SNe in the same galaxy cannot be grouped; the split is by transient name, not guaranteed by unique galaxy.',
            'Prepublication 8610-host snapshot cannot replace final 6983-host catalogue.',
            'Author notebook high-z curve is entered from W26; it is not high-z observed host ages.',
            'Age support and fitted map must be combined with selection and a refitted luminosity slope before cosmological use.'])
    write(OUT/'titan-summary.json',result)
    record=dict(completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        command='OPENBLAS_NUM_THREADS=1 .venv/bin/python studies/host_ages/code/population_transport/titan.py',
        input_sha256={r['path']:r['sha256'] for r in registry['files']},
        code_sha256={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'titan-scenarios.json',HERE/'acquire_titan.py']},
        additional_input_sha256={'data/distances/Pantheon+SH0ES.dat':sha(ROOT/'data/distances/Pantheon+SH0ES.dat')},
        output_sha256={str(p.relative_to(ROOT)):sha(p) for p in [OUT/'titan-summary.json',WORK/'titan-pantheon-literal-crosswalk.csv',WORK/'titan-mapping-predictions.csv']})
    write(OUT/'titan-run.json',record)
    print(json.dumps(result,indent=2))

if __name__=='__main__': main()
