#!/usr/bin/env python3
"""Full released-candidate SNN closure with declared timing-source comparison.

Preserve HEAD timing baseline; compare independently regenerated documented
OPT_SETPKMJD16 timing. No probability-based tuning, threshold changes or
retraining. Every released classifier CID receives a prediction or explicit
unsupported reason in every arm.
"""
import json,sys,time
import numpy as np
import pandas as pd
import torch
from astropy.io import fits
from supernnova.validation.validate_onthefly import classify_lcs
from common import ROOT,WORK,RESULTS,sha
sys.path.insert(0,str(ROOT));from lib.records import fitres

def summary(d):
    good=d.status.eq('predicted');q=d[good];delta=q.pIa_reconstructed-q.pIa_released
    return {'released_CIDs':len(d),'predicted':int(good.sum()),'unsupported':int((~good).sum()),'within_4decimal_rounding':int(delta.abs().le(5.1e-5).sum()),'mean_absolute_delta':float(delta.abs().mean()),'median_absolute_delta':float(delta.abs().median()),'max_absolute_delta':float(delta.abs().max()),'gate_0_5_disagreement':int((q.pIa_reconstructed.gt(.5)!=q.pIa_released.gt(.5)).sum()),'absolute_delta_quantiles':delta.abs().quantile([.5,.9,.95,.99,.999,1]).to_dict()}

def main():
    torch.set_num_threads(1);torch.set_num_interop_threads(1);torch.manual_seed(92731);np.random.seed(92731)
    b=WORK/'release/0_DATA/DES-SN5YR_DES';hp=b/'DES-SN5YR_DES_HEAD.FITS.gz';pp=b/'DES-SN5YR_DES_PHOT.FITS.gz';cp=next((WORK/'release/3_CLASSIFICATION').glob('*.csv'))
    ref=pd.read_csv(cp,dtype={'SNID':str}).set_index('SNID').iloc[:,0];assert len(ref)==17733
    npth=WORK/'native-dataprep/dataprep.SNANA.TEXT';native=fitres(npth);assert native.index.is_unique
    ledger=pd.read_csv(WORK/'normalized/dovekie-ledger.csv',dtype={'CID':str});accepted=set(ledger[ledger.survey_family.eq('DES')].CID)
    model=WORK/'SNDATA_ROOT_2026-04-10/models/classifiers/DES-SN5YR/SNNTRAINV19_z_TRAINDES_V19/model.pt'
    with fits.open(hp) as h:heads={str(r['SNID']).strip():r for r in h[1].data.copy()}
    plan={'code_sha256':sha(__file__),'CIDs':sorted(ref.index),'modes':['public_HEAD_PEAKMJD','native_OPT_SETPKMJD16_PKMJDINI'],'window':[-30,100],'photflag_mask':1016,'redshift_feature':'REDSHIFT_FINAL plusREDSHIFT_FINAL_ERR, documentedPippin override','native_source_sha256':sha(npth),'comparison_after_initial256_discrepancies':True,'no_tuning_to_probabilities':True}
    (WORK/'classifier-full-design.json').write_text(json.dumps(plan,indent=2)+'\n')
    results={};start=time.monotonic()
    with fits.open(pp) as ph:
        phot=ph[1].data
        for mode in plan['modes']:
            out=[]
            for chunk in np.array_split(sorted(ref.index),int(np.ceil(len(ref)/128))):
                frames=[];metadata=[]
                for name in chunk:
                    h=heads[name];peak=float(h['PEAKMJD']) if mode.startswith('public') else float(native.loc[name,'PKMJDINI']) if name in native.index else np.nan
                    meta={'CID':name,'pIa_released':float(ref[name]),'in_Dovekie':name in accepted,'spectroscopic_type':int(h['SNTYPE']),'public_peak':float(h['PEAKMJD']),'used_peak':peak,'z_feature':float(h['REDSHIFT_FINAL'])}
                    if not np.isfinite(peak) or peak<50000:meta['status']='unsupported_peak';out.append(meta);continue
                    q=phot[int(h['PTROBS_MIN'])-1:int(h['PTROBS_MAX'])];win=(q['MJD']-peak>-30)&(q['MJD']-peak<100);q=q[win];valid=(q['PHOTFLAG']&1016)==0;q=q[valid]
                    meta.update(epochs_retained=len(q),epochs_window_rejected=int((~win).sum()),epochs_flag_rejected=int((~valid).sum()))
                    if not len(q):meta['status']='unsupported_no_epochs';out.append(meta);continue
                    bands=np.char.strip(np.asarray(q['BAND']).astype(str));assert set(bands)<=set('griz')
                    frame=pd.DataFrame({'SNID':name,'MJD':q['MJD'].astype(float),'FLT':bands,'FLUXCAL':q['FLUXCAL'].astype(float),'FLUXCALERR':q['FLUXCALERR'].astype(float),'HOSTGAL_SPECZ':float(h['REDSHIFT_FINAL']),'HOSTGAL_SPECZ_ERR':float(h['REDSHIFT_FINAL_ERR'])})
                    frames.append(frame);metadata.append(meta)
                if not frames:continue
                names,pred=classify_lcs(pd.concat(frames,ignore_index=True),str(model),'cpu');assert pred.shape==(len(metadata),1,2)
                lookup={str(n):float(v[0,0]) for n,v in zip(names,pred)}
                for m in metadata:m.update(status='predicted',pIa_reconstructed=lookup[m['CID']]);out.append(m)
            d=pd.DataFrame(out);assert len(d)==len(ref) and d.CID.is_unique
            path=WORK/f'classifier-full-{mode}.csv';d.to_csv(path,index=False)
            results[mode]={'all':summary(d),'Dovekie':summary(d[d.in_Dovekie]),'spectroscopic_Ia':summary(d[d.spectroscopic_type.eq(1)]),'unsupported_reasons':d.status.value_counts().to_dict(),'rows_sha256':sha(path),'path':str(path.relative_to(ROOT))}
            print(mode,json.dumps(results[mode]),flush=True)
    record={'status':'executed_exact_closure_unless_all_rows_rounding_pass','code_sha256':sha(__file__),'seconds':time.monotonic()-start,'plan_sha256':sha(WORK/'classifier-full-design.json'),'results':results,'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [hp,pp,cp,npth,model,model.parent/'cli_args.json',model.parent/'data_norm.json'] if p.exists()},'interpretation':'Released-model conditionalinference withdocumentedmeasuredinputs; fullproductionclaim requires allrequestedrows andpreprocessing/equivalence gates, notmerely highclassificationagreement.'}
    (RESULTS/'classifier-full.json').write_text(json.dumps(record,indent=2)+'\n')
if __name__=='__main__':main()
