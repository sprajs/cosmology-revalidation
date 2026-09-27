#!/usr/bin/env python3
"""Selected-sample redshift responses and paired Monte Carlo covariance.

Postprocessing of the fully retained flux injections. Bins are descriptive,
chosen after the pilot; they are not new independent hypothesis tests.
"""
from pathlib import Path
import hashlib,json
from datetime import datetime,timezone
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[4]
OUT=ROOT/'studies/host_ages/results/age_recovery'
WORK=ROOT/'.work/age-recovery'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    run=json.loads((OUT/'flux-recovery.json').read_text());p=WORK/'flux-event-ledger.csv'
    assert sha(p)==run['provenance']['event_ledger_sha256']
    data=pd.read_csv(p);data=data[(data.flux_scale==1)&data.selected&data.success&~data.boundary].copy()
    edges=[0,.3,.6,.9,1.3];data['zbin']=pd.cut(data.z,edges,labels=False)
    curves={}
    for arm,d in data.groupby('arm'):
        a=d.pivot_table(index='draw',columns='zbin',values='delta_Tripp',aggfunc='mean').reindex(index=range(12),columns=range(4))
        assert a.notna().all().all()
        values=a.to_numpy();values-=values[:,[0]]
        curves[arm]=values
    rows={}
    for arm,v in curves.items():
        diff=v-curves['zero'];cov=np.cov(diff,rowvar=False,ddof=1)
        rows[arm]={'mean_response_relative_low_z_mag':diff.mean(axis=0).tolist(),'paired_noise_covariance_mag2':cov.tolist(),'covariance_of_estimated_mean_mag2':(cov/len(diff)).tolist(),'per_bin_selected_event_draws':data[data.arm==arm].groupby('zbin').size().tolist(),'distinct_cadences':data[data.arm==arm].groupby('zbin').CID.nunique().tolist()}
    report={'created_utc':datetime.now(timezone.utc).isoformat(),'scope':'Each arm recomputes synthetic detection on the same 12 eligible observed cadences and 12 paired noise draws. This is a response to specified perturbations relative to zero injection, with equal selected-event weights and reference z<.3. It is not the full-population cosmological B_std; training, actual selection, covariance regeneration and host-age likelihood are missing.','bin_edges':edges,'reference_bin':0,'bin_choice':'descriptive post-pilot binning, not a prespecified significance test','all_draws_have_all_bins':True,'draws':12,'arms':rows,'code_sha256':sha(Path(__file__)),'input_record_sha256':sha(OUT/'flux-recovery.json'),'event_ledger_sha256':sha(p)}
    (OUT/'flux-redshift-response.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    for k,v in rows.items():print(k,v['mean_response_relative_low_z_mag'])
if __name__=='__main__':main()
