#!/usr/bin/env python3
"""Observed-light-curve closure of published SNN weights, not classifier retraining.

Before predictions select256 secure-host-z objects by SHA256(CID), independent
of class probabilities. Use public PEAKMJD, published flag/window convention;
preserve discrepancies with released rounded probabilities rather than tune.
"""
import hashlib,json,time
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from astropy.io import fits
from supernnova.validation.validate_onthefly import classify_lcs
from common import ROOT,WORK,RESULTS,sha

def main():
    torch.set_num_threads(1);torch.set_num_interop_threads(1);torch.manual_seed(92731);np.random.seed(92731)
    base=WORK/'release/0_DATA/DES-SN5YR_DES';hp=base/'DES-SN5YR_DES_HEAD.FITS.gz';pp=base/'DES-SN5YR_DES_PHOT.FITS.gz'
    cp=next((WORK/'release/3_CLASSIFICATION').glob('*.csv'));ref=pd.read_csv(cp,dtype={'SNID':str}).set_index('SNID');reference=ref.iloc[:,0].to_dict()
    with fits.open(hp) as h:heads=h[1].data.copy()
    eligible=[str(r['SNID']).strip() for r in heads if str(r['SNID']).strip() in reference and r['REDSHIFT_FINAL']>0 and r['PEAKMJD']>50000]
    ids=set(sorted(eligible,key=lambda c:hashlib.sha256(c.encode()).hexdigest())[:256]);assert len(ids)==256
    design={'code_sha256':sha(__file__),'sample_size':256,'eligible':len(eligible),'CIDs':sorted(ids),'selection':'lowestSHA256CID among released-classified positive measuredredshift andPEAKMJD>50000','peak':'publicHEADPEAKMJD, not regeneratedoriginalPKMJDINI','window_days':[-30,100],'badflag_bitmask':1016,'seed':92731}
    (WORK/'classifier-closure-design.json').write_text(json.dumps(design,indent=2)+'\n')
    rows=[];meta=[]
    with fits.open(pp) as f:
        phot=f[1].data
        for h in heads:
            name=str(h['SNID']).strip()
            if name not in ids:continue
            q=phot[int(h['PTROBS_MIN'])-1:int(h['PTROBS_MAX'])];window=(q['MJD']-h['PEAKMJD']>-30)&(q['MJD']-h['PEAKMJD']<100);q=q[window];good=(q['PHOTFLAG']&1016)==0;q=q[good]
            meta.append({'CID':name,'epochs_retained':len(q),'epochs_window_rejected':int((~window).sum()),'epochs_flag_rejected':int((~good).sum())})
            assert len(q)>0
            for o in q:
                band=str(o['BAND']).strip().replace('DES-','');assert band in 'griz'
                rows.append({'SNID':name,'MJD':float(o['MJD']),'FLT':band,'FLUXCAL':float(o['FLUXCAL']),'FLUXCALERR':float(o['FLUXCALERR']),'HOSTGAL_SPECZ':float(h['REDSHIFT_FINAL']),'HOSTGAL_SPECZ_ERR':float(h['REDSHIFT_FINAL_ERR'])})
    frame=pd.DataFrame(rows)
    model=WORK/'SNDATA_ROOT_2026-04-10/models/classifiers/DES-SN5YR/SNNTRAINV19_z_TRAINDES_V19/model.pt'
    start=time.monotonic();out=[]
    for chunk in np.array_split(sorted(ids),4):
        names,pred=classify_lcs(frame[frame.SNID.isin(chunk)].copy(),str(model),'cpu');assert pred.shape==(len(chunk),1,2)
        out.extend({'CID':str(n),'pIa_reconstructed':float(p[0,0]),'pIa_released':reference[str(n)]} for n,p in zip(names,pred))
    d=pd.DataFrame(out).merge(pd.DataFrame(meta),on='CID',validate='one_to_one');d['delta']=d.pIa_reconstructed-d.pIa_released;path=WORK/'classifier-closure.csv';d.to_csv(path,index=False)
    maxdiff=float(d.delta.abs().max());samegate=int((d.pIa_reconstructed.gt(.5)==d.pIa_released.gt(.5)).sum())
    result={'status':'executed_exact_closure_not_established' if maxdiff>5.1e-5 else 'rounding_closure_passed','code_sha256':sha(__file__),'N':len(d),'seconds':time.monotonic()-start,'max_abs_probability_difference':maxdiff,'mean_abs_probability_difference':float(d.delta.abs().mean()),'median_abs_probability_difference':float(d.delta.abs().median()),'pIa_0_5_gate_agree':samegate,'pIa_0_5_gate_disagree':len(d)-samegate,'within_released_4decimal_rounding':int(d.delta.abs().le(5.1e-5).sum()),'design_sha256':sha(WORK/'classifier-closure-design.json'),'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [hp,pp,cp,model]},'output':{str(path.relative_to(ROOT)):sha(path)},'limitations':['Observed HEADPEAKMJD and REDSHIFT_FINAL replace unavailableoriginalcompletepreprocessing state; this is a declared closuretest, not tuned reconstruction.','Classifier uses released pretrainedweights, not newtraining; probability discrepancy blocks a claimofexactoriginalclassificationclosure.','Published probabilities can be safely consumedaspartoftheconditionalreleaseddistance likelihood without this classifier substituting newprobabilities.']}
    (RESULTS/'classifier-closure.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
