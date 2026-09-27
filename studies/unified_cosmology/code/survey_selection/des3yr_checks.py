#!/usr/bin/env python3
"""Executed support audit of a separate spectroscopically typed DES3YR route."""
import gzip,json,sys
from collections import Counter
import numpy as np
import pandas as pd
from common import ROOT,WORK,RESULTS,sha
sys.path.insert(0,str(ROOT));from lib.records import fitres

def main():
    base=WORK/'des3yr';record={'code_sha256':sha(__file__)};frames={};sources=[]
    ledger=pd.read_csv(WORK/'normalized/dovekie-ledger.csv',dtype={'CID':str}).set_index('CID')
    for label in ['G10','C11']:
        path=base/f'04-BBCFITS/SALT2mu_DES+LOWZ_{label}.FITRES';d=fitres(path);sources.append(path);assert d.index.is_unique
        des=d[d.IDSURVEY.eq(10)];assert len(des)==207
        phot=[]
        for cid in d.index:
            group='DES' if d.loc[cid,'IDSURVEY']==10 else 'LOWZ'
            files=[p for p in (base/'02-DATA_PHOTOMETRY'/f'DES-SN3YR_{group}').iterdir() if p.suffix.lower()=='.dat']
            matched=[]
            for p in files:
                for line in p.open():
                    if line.startswith('SNID:'):
                        if line.split()[1]==cid:matched.append(p)
                        break
            assert len(matched)==1,(cid,len(matched));phot.append(matched[0])
        closure=float(np.max(abs(d.MU-d.MUMODEL-d.MURES-d.M0DIF)));assert closure<.0002
        shared=[c for c in des.index if c in ledger.index]
        pos=np.hypot((des.loc[shared,'RA'].to_numpy()-ledger.loc[shared,'RA'].to_numpy())*np.cos(np.deg2rad(des.loc[shared,'DECL'])),des.loc[shared,'DECL'].to_numpy()-ledger.loc[shared,'DEC'].to_numpy())*3600
        assert max(pos)<1
        record[label]={'rows':len(d),'DES':len(des),'low_z':len(d)-len(des),'all_rows_have_exact_SNID_photometry':True,'MU_MUMODEL_MURES_M0DIF_identity_max':closure,'shared_DES_Dovekie_events':len(shared),'shared_max_sky_arcsec':float(max(pos)),'DES_redshift_range':[float(des.zHD.min()),float(des.zHD.max())]}
        frames[label]=d
        path=base/f'03-SIM_BIASCOR/BIASCOR_{label}.FITRES.gz';sources.append(path)
        with gzip.open(path,'rt') as f:
            for line in f:
                if line.startswith('VARNAMES:'):names=line.split()[1:];break
            ia=names.index('SIM_alpha')+1;ib=names.index('SIM_beta')+1;it=names.index('SIM_TYPE_INDEX')+1;isur=names.index('IDSURVEY')+1
            grid=Counter();types=Counter();surveys=Counter();n=0
            for line in f:
                if not line.startswith('SN:'):continue
                row=line.split();n+=1;grid[(float(row[ia]),float(row[ib]))]+=1;types[row[it]]+=1;surveys[row[isur]]+=1
        assert len(grid)==4 and set(types)=={'1'}
        record[label]['bias_training']={'rows':n,'intrinsic_scatter':label,'alpha_beta_grid_counts':{f'{a},{b}':v for (a,b),v in grid.items()},'type_counts':dict(types),'survey_counts':dict(surveys),'contains_physical_age_or_extinction_likelihood':False}
    assert np.array_equal(frames['G10'].index,frames['C11'].index)
    for col in ['x0','x1','c','mB','PKMJD']:assert np.array_equal(frames['G10'][col],frames['C11'][col])
    specpath=base/'02-DATA_SPECTRA/DES-SN3YR_CLASSIFICATIONS.LIST';sources.append(specpath)
    spec=pd.read_csv(specpath,sep=r'\s+',comment='#',header=None,names=['name','CID','RA','DEC','telescope','date','type','z'],dtype={'CID':str})
    selected=set(frames['G10'][frames['G10'].IDSURVEY.eq(10)].index);q=spec[spec.CID.isin(selected)]
    labels=q.groupby('CID')['type'].agg(set)
    definite=int(sum('SNIa' in s for s in labels));provisional=int(sum(s=={'SNIa?'} for s in labels))
    assert definite+provisional==207
    record['spectroscopic_evidence']={'classification_rows':len(spec),'classified_transients':int(spec.CID.nunique()),'DES207_with_spectrum_classification':int(q.CID.nunique()),'accepted_spectral_types':q.type.value_counts().to_dict(),'DES_with_at_least_one_SNIa_label':definite,'DES_with_only_provisional_SNIa_question_mark_labels':provisional,'labels_are_released_classifications_not_independently_certified_truth':True,'repeated_exposures_not_new_SNe':True}
    record['status']='passed_historical_support_checks'
    record['inputs']={str(p.relative_to(ROOT)):sha(p) for p in sources}
    record['limits']=['207typedDESobjects have a matchinghistoricalspecselection model; do not applyit to353typed5yrSNe.','These publicbias-training outputs contain aboutonemillion events on a2x2alpha/beta grid, not newlygeneratedphysicalage/dustresponses.','G10/C11 are alternativeintrinsicscatter models, not independentobservations.','Original releasedraw/fit/cosmology products arehistoricalcalibration; fullnewphysicalinference stillneeds population/selection normalization, detector/systematicresponse and hostlikelihood.','Dovekie sharesphysicalevents; two compilations cannotbe multipliedas independentlikelihoods.']
    (RESULTS/'des3yr-support.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))
if __name__=='__main__':main()
