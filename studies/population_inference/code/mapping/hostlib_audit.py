"""Summaries of actual author-supplied W22 simulated host rows, NOT real SN observations.

Equal row weights: generating code has already sampled rate×mass-function.
No inferred W26 selection function is imposed and no low-z extrapolation is made.
"""
from common import ROOT, OUT, manifest
import numpy as np
import pandas as pd
import json

def main():
    repo=ROOT/'sources/repos/wisemanp__des_sn_hosts'
    path=repo/'simulations/data/Hostlib_Allz_Quenched_Bursts_16062022_combined.csv.gz'
    cols=['GALID','ZTRUE','LOGMASS','LOG_sSFR','mean_age','pred_rate_total','SN_age','Av']
    df=pd.read_csv(path,usecols=cols)
    assert not df[cols].isna().any().any()
    df['host_age_gyr']=df.mean_age/1000
    rows=[]
    for z,d in df.groupby('ZTRUE',sort=True):
        lo=d.LOGMASS<10;hi=~lo
        rows.append(dict(z=z,n=len(d),mean_host_age_gyr=d.host_age_gyr.mean(),
            mean_delay_gyr=d.SN_age.mean(),median_delay_gyr=d.SN_age.median(),
            sd_delay_gyr=d.SN_age.std(),high_mass_fraction=hi.mean(),
            mean_host_lowmass=d.loc[lo,'host_age_gyr'].mean(),mean_host_highmass=d.loc[hi,'host_age_gyr'].mean(),
            mean_delay_lowmass=d.loc[lo,'SN_age'].mean(),mean_delay_highmass=d.loc[hi,'SN_age'].mean(),
            naive_delay_on_host_slope=(np.cov(d.host_age_gyr,d.SN_age,ddof=0)[0,1]/np.var(d.host_age_gyr)
                                      if np.var(d.host_age_gyr)>1e-15 else np.nan)))
    res=pd.DataFrame(rows)
    # This is a source-product descriptive quantity only; no fit or observational uncertainty.
    res['illustrative_host_age_mass_step_mag']=.03*(res.mean_host_highmass-res.mean_host_lowmass)
    res['illustrative_delay_mass_step_mag']=.03*(res.mean_delay_highmass-res.mean_delay_lowmass)
    p=OUT/'author-hostlib-summary.csv';res.to_csv(p,index=False)
    q=OUT/'author-hostlib-audit.json'
    summary=dict(data_type='author-supplied simulated host library, not observed galaxies',
        rows=len(df), redshift_count=len(res), min_z=float(df.ZTRUE.min()),max_z=float(df.ZTRUE.max()),
        min_host_gyr=float(df.host_age_gyr.min()),max_host_gyr=float(df.host_age_gyr.max()),
        min_delay_gyr=float(df.SN_age.min()),max_delay_gyr=float(df.SN_age.max()),
        unique_galid=int(df.GALID.nunique()),
        raw_unique_redshifts=int(df.ZTRUE.nunique()),rounded_12dp_unique_redshifts=int(df.ZTRUE.round(12).nunique()),
        groups_not_10000_rows=int((res.n!=10000).sum()),
        weights='equal supplied rows; rate already sampled upstream; no extra rate weighting',
        status='genuine source-product summary; W26 figures not reproduced; selection and exact run bundle unavailable')
    q.write_text(json.dumps(summary,indent=2)+'\n')
    manifest('author-hostlib',dict(seed=None,host_age_unit_input='Myr',delay_unit_input='Gyr',selection_applied=False),
        [__file__,ROOT/'scripts/mapping/common.py',ROOT/'docs/experiments/mapping-plan.md',ROOT/'uv.lock',path,
         repo/'simulations/data/README.txt',repo/'simulations/scripts/make_hostlib_snana.py',
         repo/'simulations/aura.py',repo/'simulations/utils/gal_functions.py'],[p,q])
    print(json.dumps(summary,indent=2));print(res.loc[(res.z<.002)|np.isclose(res.z,.5)|np.isclose(res.z,1)|np.isclose(res.z,1.2)].to_string(index=False))

if __name__=='__main__':main()
