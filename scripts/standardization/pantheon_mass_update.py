#!/usr/bin/env python3
"""STD-03 independent Vandermonde reconstruction of a public mass revision."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,subprocess
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[2]
SRC=ROOT/'sources/updates/2026-09-20-standardization/pantheonplus-hostmass-correction-reconstruction-c5583379eb06f9a4b045353bb39a0991f8c15e4c'
OLD=ROOT/'sources/repos/PantheonPlusSH0ES__DataRelease@7fc6805/Pantheon+_Data/4_DISTANCES_AND_COVAR/Pantheon+SH0ES.dat'
CURRENT=ROOT/'sources/repos/PantheonPlusSH0ES__DataRelease/Pantheon+_Data/4_DISTANCES_AND_COVAR/Pantheon+SH0ES.dat'
OUT=ROOT/'runs/standardization/pantheon_mass_update';OUT.mkdir(parents=True,exist_ok=True)
DER=ROOT/'data/derived/standardization';DER.mkdir(parents=True,exist_ok=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return pd.read_csv(p,sep=r'\s+',dtype={'CID':str})
def main():
    original=SRC/'data/Pantheon+SH0ES.dat'; published=SRC/'data/Pantheon+SH0ES_HostMassCorrected_Reconstructed.dat'
    d=read(original); target=read(published)
    z=d.zCMB.ge(.01)&d.zCMB.lt(.15);lo=z&d.HOST_LOGMASS.lt(10);hi=z&d.HOST_LOGMASS.ge(10)
    affected=lo&d.HOST_LOGMASS.ge(9.4)
    x=np.vander(d.c,4)
    coeff_lo=np.linalg.lstsq(x[lo],d.loc[lo,'biasCor_m_b'],rcond=None)[0]
    coeff_hi=np.linalg.lstsq(x[hi],d.loc[hi,'biasCor_m_b'],rcond=None)[0]
    delta=np.zeros(len(d));delta[affected]=x[affected]@(coeff_lo-coeff_hi)
    updated=d.copy();updated['m_b_corr']+=delta
    path=DER/'Pantheon_W26_lowz_mass_revision.dat';updated.to_csv(path,sep=' ',index=False,float_format='%.17g')
    vector=d[['CID','IDSURVEY']].copy()
    vector.insert(0,'row_index_zero_based',np.arange(len(d)))
    vector['delta_m_b_corr']=delta
    vector.to_csv(DER/'Pantheon_W26_lowz_mass_revision_delta.csv',index=False,float_format='%.17g')
    cols=[c for c in d.columns if c not in ['CID','m_b_corr']]
    assert np.max(np.abs(updated.m_b_corr-target.m_b_corr))<1e-12
    assert np.array_equal(d.CID,target.CID)
    assert np.allclose(d[cols].astype(float),target[cols].astype(float),atol=0,rtol=0)
    # Original author inputs must agree before using inherited covariance/order.
    comparisons={}
    for name,p in [('W26_snapshot',OLD),('current_acquisition',CURRENT)]:
        a=read(p)
        comparisons[name]={'sha256_equal':sha(p)==sha(original),'all_values_equal':bool(a.equals(d)),
                           'CID_survey_order_equal':bool(a[['CID','IDSURVEY']].equals(d[['CID','IDSURVEY']]))}
        assert comparisons[name]['all_values_equal']
    audit=d.loc[affected,['CID','IDSURVEY','zCMB','zHD','HOST_LOGMASS','c','m_b_corr']].copy()
    audit['delta_m_b_corr']=delta[affected];audit['m_b_corr_new']=updated.loc[affected,'m_b_corr']
    audit.to_csv(OUT/'affected_rows.csv',index=False)
    res={'release_comparisons':comparisons,'rows':len(d),'low_z_rows':int(z.sum()),'low_mass_fit_rows':int(lo.sum()),'high_mass_fit_rows':int(hi.sum()),
         'changed_rows':int(affected.sum()),'changed_unique_CIDs':int(d.loc[affected,'CID'].nunique()),
         'changed_calibrator_rows':int(d.loc[affected,'IS_CALIBRATOR'].sum()),'changed_SDSS_rows':int((d.loc[affected,'IDSURVEY']==1).sum()),
         'mean_delta_m_b_corr_changed':float(delta[affected].mean()),'mean_delta_m_b_corr_all_lowz':float(delta[z].mean()),
         'range_delta_m_b_corr_changed':[float(delta[affected].min()),float(delta[affected].max())],
         'max_abs_difference_to_public_reconstructed_table':float(np.max(np.abs(updated.m_b_corr-target.m_b_corr))),
         'descending_cubic_coefficients_low_minus_high':(coeff_lo-coeff_hi).tolist(),
         'unchanged_columns':'Every column except m_b_corr (including MU_SH0ES, masses, errors, biasCor_m_b). Use m_b_corr for this branch.',
         'uncertainty':'No covariance update. Does not propagate host-mass shift or curve-fit uncertainty; deterministic sensitivity only.',
         'sign':'NEW minus OLD standardized apparent magnitude is positive on average (+0.07220 mag); public file and -biasCor sign support this. Hoyt Appendix F quotes Delta mu=-0.073; retain sign discrepancy/convention explicitly.',
         'reproduction':'Exact numerical reproduction of public Roy Choudhury branch-transfer table; comparison to private Hoyt file only asserted by public author, not independently verified here.'}
    (OUT/'results.json').write_text(json.dumps(res,indent=2)+'\n')
    inputs=[original,published,OLD,CURRENT,SRC/'scripts/build_corrected_pantheonplus.py',ROOT/'sources/updates/2026-09-20-standardization/2601.19424.pdf',ROOT/'sources/updates/2026-09-20-standardization/2607.24443.pdf']
    manifest={'experiment':'STD-03','utc':datetime.now(timezone.utc).isoformat(),'code_revision':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
              'script_sha256':sha(Path(__file__)),'plan_sha256':sha(ROOT/'docs/experiments/standardization-plan.md'),
              'inputs':{str(p.relative_to(ROOT)):sha(p) for p in inputs},'seed':None,'environment':{'numpy':np.__version__,'pandas':pd.__version__},
              'outputs':[str(path.relative_to(ROOT)),'data/derived/standardization/Pantheon_W26_lowz_mass_revision_delta.csv',str((OUT/'results.json').relative_to(ROOT)),str((OUT/'affected_rows.csv').relative_to(ROOT))]}
    (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(res,indent=2))
if __name__=='__main__':main()
