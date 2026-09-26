"""Independent array arithmetic check of completed nonlinear validation."""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
from scipy.linalg import solve_triangular

ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent/'resolved'
BASE=ROOT/'runs/research_2026_09_26/astra_design/validation1020'
cohort=pd.read_csv(HERE/'cohort.csv',dtype={'CID':str})
membership=pd.read_csv(HERE/'contrast-membership.csv',dtype={'CID':str})
old_membership=pd.read_csv(ROOT/'runs/research_2026_09_26/shared_distance_response_12/contrast-membership.csv',dtype={'CID':str})
assert membership.set_index('CID')[['low_quartile','high_quartile']].sort_index().equals(old_membership.set_index('CID')[['low_quartile','high_quartile']].sort_index())
root=json.loads((HERE/'result.json').read_text())
assert root['all_gates_pass'] and root['objects']==1020
objective_error=0.;shift_error=0.;comparisons={};scores={};rows=[]
for cid in cohort.CID:
    record=json.loads((HERE/'objects'/f'{cid}.json').read_text())
    assert record['all_gates_pass'] and record['status']=='complete'
    with np.load(BASE/'objectives'/f'objective_{cid}.npz') as d:
        data=d['data_flux'];official=d['model_flux'];cov=d['frozen_flux_covariance'];time=d['MJD'];band=d['band']
    order=pd.DataFrame({'t':time,'b':band}).sort_values(['t','b'],kind='stable').index.to_numpy()
    with np.load(BASE/'analysis/objects'/f'{cid}.npz') as d:
        native=np.empty(len(time));native[order]=d['native_flux_model']
    targets={'native_noiseless':native,'observed_flux':data,'official_mean_sensitivity':official}
    chol=np.linalg.cholesky(cov)
    with np.load(HERE/'objects'/f'{cid}-predictions.npz') as prediction:
        changes={}
        for response in record['responses']:
            target=response['target'];mode=response['mode']
            base=prediction[target+'__nominal'];alt=prediction[target+'__'+mode]
            rb=solve_triangular(chol,base-targets[target],lower=True)
            ra=solve_triangular(chol,alt-targets[target],lower=True)
            objective_error=max(objective_error,abs(float(ra@ra)-response['alternative_chi2']),abs(float(rb@rb)-response['baseline_chi2']))
            delta=np.array(response['alternative_theta'])-np.array(response['baseline_theta'])
            mu=float(np.array([1.,.16087,-3.1178,0.])@delta)
            shift_error=max(shift_error,abs(mu-response['nonlinear_standardized_mag']))
            correction=solve_triangular(chol,alt-base,lower=True)
            changes.setdefault(target,{})[mode]=correction
            if target=='observed_flux':scores[mode]=scores.get(mode,0)+float(.5*(rb@rb-ra@ra))
            rows.append(dict(CID=cid,target=target,mode=mode,mu=mu))
        for target,values in changes.items():
            a=values['observer'];b=values['sed']
            result=comparisons.setdefault(target,dict(observer=0.,sed=0.,cross=0.,difference=0.))
            result['observer']+=float(a@a);result['sed']+=float(b@b)
            result['cross']+=float(a@b);result['difference']+=float((a-b)@(a-b))
assert objective_error<1e-8 and shift_error<1e-10
for mode,value in scores.items():assert abs(value-root['descriptive_observed_local_fit_scores'][mode]['descriptive_profile_log_score_change'])<1e-8
low=set(membership.loc[membership.low_quartile,'CID']);high=set(membership.loc[membership.high_quartile,'CID'])
frame=pd.DataFrame(rows);contrast_error=0.
for target in ['native_noiseless','observed_flux']:
    for mode in ['observer','sed']:
        d=frame[(frame.target==target)&(frame['mode']==mode)]
        contrast=float(d.loc[d.CID.isin(high),'mu'].mean()-d.loc[d.CID.isin(low),'mu'].mean())
        saved=root['exact_high255_low255_contrasts'][target+'__'+mode]['nonlinear_standardized_mag']
        contrast_error=max(contrast_error,abs(contrast-saved))
        assert contrast_error<1e-10
    c=comparisons[target];calc=c['difference']/c['observer']
    assert abs(calc-root['nonlinear_mean_change_comparisons'][target]['relative_squared_mismatch'])<1e-10
six=json.loads((ROOT/'runs/research_2026_09_26/sed_nonlinear_refit/result.json').read_text())
six_error=0.
for old in six['responses']:
    if old['mask']!='accepted_full':continue
    d=frame[(frame.CID==old['CID'])&(frame.target==old['target'])&(frame['mode']==old['mode'])]
    assert len(d)==1
    six_error=max(six_error,abs(float(d.mu.iloc[0])-old['nonlinear_standardized_mag']))
assert six_error<1e-8
report=dict(objects=1020,all_gates_pass=True,exact_membership_preserved=True,max_objective_error=objective_error,max_standardized_shift_error=shift_error,max_contrast_error=contrast_error,max_prior_six_object_response_error=six_error,descriptive_scores_recomputed=scores,inputs_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),HERE/'result.json',HERE/'freeze.json']})
(HERE/'independent-verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
