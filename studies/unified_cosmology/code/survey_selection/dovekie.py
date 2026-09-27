#!/usr/bin/env python3
"""Strict adapter for the released corrected-distance likelihood, not raw flux.

Published NPZ arrays contain upper-triangular TOTAL PRECISION. For a subset,
invert the full precision, subset covariance, then invert that covariance.
"""
import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.linalg import cho_factor,cho_solve
from astropy.io import fits
from common import ROOT,HERE,WORK,RESULTS,sha
sys.path.insert(0,str(ROOT))
from lib.records import fitres
RELEASE=WORK/'release'


def unpack(path):
    with np.load(path,allow_pickle=False) as d:
        n=int(d['nsn'][0]);v=np.asarray(d['cov'],dtype=float)
    assert v.shape==(n*(n+1)//2,)
    p=np.zeros((n,n));p[np.triu_indices(n)]=v;p=p+p.T-np.diag(np.diag(p))
    assert np.isfinite(p).all()
    factor=cho_factor(p,lower=True)
    cov=cho_solve(factor,np.eye(n))
    cov=(cov+cov.T)/2
    return p,cov


def load(mask=None,stat_only=False):
    folder=RELEASE/'4_DISTANCES_COVMAT'
    hd=fitres(folder/'DES-Dovekie_HD.csv')
    precision,covariance=unpack(folder/('STATONLY.npz' if stat_only else 'STAT+SYS.npz'))
    assert len(hd)==len(precision) and hd.zHD.gt(0).all()
    if mask is not None:
        mask=np.asarray(mask)
        assert mask.dtype==bool and mask.shape==(len(hd),)
        hd=hd.iloc[np.flatnonzero(mask)]
        covariance=covariance[np.ix_(mask,mask)]
        precision=cho_solve(cho_factor(covariance,lower=True),np.eye(len(hd)))
    return {'cid':hd.index.to_numpy(),'survey':hd.IDSURVEY.to_numpy(),
            'zHD':hd.zHD.to_numpy(),'zHEL':hd.zHEL.to_numpy(),'MU':hd.MU.to_numpy(),
            'precision':precision,'covariance':covariance,
            'scope':'Released corrected distances, conditional on author light-curve, bias/CC treatment and calibration covariance; no independent raw-survey likelihood.'}


def marginalized_loglike(mu_model, data, include_constant=False):
    residual=np.asarray(mu_model)-data['MU'];p=data['precision']
    one=np.ones(len(residual));u=p@one
    chi2=residual@p@residual-(residual@u)**2/(one@u)
    if include_constant:
        chi2+=np.log(one@u/(2*np.pi))+np.linalg.slogdet(data['covariance'])[1]
    return float(-.5*chi2)


def main():
    folder=RELEASE/'4_DISTANCES_COVMAT'
    hd=fitres(folder/'DES-Dovekie_HD.csv')
    meta=fitres(folder/'DES-Dovekie_Metadata.csv').loc[hd.index]
    assert hd.index.is_unique and np.array_equal(hd.IDSURVEY,meta.IDSURVEY)
    assert np.max(abs(hd.zHD-meta.zHD))<1e-8 and np.max(abs(hd.MU-meta.MU))<1e-8
    heads={}
    for group in ['DES','Foundation','LOWZ']:
        path=RELEASE/'0_DATA'/f'DES-SN5YR_{group}'/f'DES-SN5YR_{group}_HEAD.FITS.gz'
        with fits.open(path) as h:
            rows=h[1].data
            heads[group]={str(r['SNID']).strip():{k:r[k].item() if isinstance(r[k],np.generic) else r[k]
                         for k in ['RA','DEC','SNTYPE','HOSTGAL_OBJID','NAME_IAUC','NAME_TRANSIENT','IAUC'] if k in rows.names} for r in rows}
        assert len(heads[group])==len(rows)
    ledger=[]
    for row,(cid,d) in enumerate(hd.iterrows()):
        group='DES' if d.IDSURVEY==10 else 'Foundation' if d.IDSURVEY==150 else 'LOWZ'
        head=heads[group].get(cid)
        assert head is not None,(cid,group)
        iau=str(head.get('NAME_IAUC',head.get('IAUC',''))).strip()
        goodiau=iau not in ['', 'UNKNOWN','NULL','-9']
        event='IAU:'+iau if goodiau else group+':'+cid
        ledger.append({'row':row,'physical_ID':event,'CID':cid,'IDSURVEY':int(d.IDSURVEY),
                       'survey_family':group,'zHD':d.zHD,'zHEL':d.zHEL,'MU':d.MU,
                       'MUERR':d.MUERR,'RA':head['RA'],'DEC':head['DEC'],
                       'IAU_alias':iau if goodiau else '',
                       'transient_alias':str(head.get('NAME_TRANSIENT','')).strip(),
                       'host_identifier':str(head['HOSTGAL_OBJID']),
                       'spectroscopic_type':int(head['SNTYPE']),
                       'host_mass_summary':meta.loc[cid,'HOST_LOGMASS'],
                       'calibration_release':'DES-Dovekie; calibrated SMP flux with calibration offsets in KCOR/model, total covariance already includes released calibration systematics'})
    frame=pd.DataFrame(ledger);assert frame.physical_ID.is_unique
    out=WORK/'normalized';out.mkdir(exist_ok=True)
    frame.to_csv(out/'dovekie-ledger.csv',index=False,float_format='%.17g')
    checks={}
    for label,stat in [('total',False),('stat',True)]:
        d=load(stat_only=stat);p,c=d['precision'],d['covariance'];n=len(p)
        checks[label]={'n':n,'precision_min_eigenvalue':float(np.linalg.eigvalsh(p).min()),
                       'covariance_min_eigenvalue':float(np.linalg.eigvalsh(c).min()),
                       'inverse_identity_max':float(np.max(abs(p@c-np.eye(n)))),
                       'covariance_diagonal_error_range':np.sqrt(np.diag(c))[[np.argmin(np.diag(c)),np.argmax(np.diag(c))]].tolist()}
        assert checks[label]['inverse_identity_max']<1e-10
        np.savez_compressed(out/f'dovekie-{label}.npz',precision=p,covariance=c,
                            CID=d['cid'].astype(str),zHD=d['zHD'],zHEL=d['zHEL'],MU=d['MU'])
    data=load();trial=data['MU']+np.linspace(-.1,.1,len(frame))
    checks['global_offset_loglike_invariance']=abs(marginalized_loglike(trial,data)-marginalized_loglike(trial+3,data))
    assert checks['global_offset_loglike_invariance']<1e-7
    mask=frame.survey_family.eq('DES').to_numpy();sub=load(mask=mask)
    checks['DES_subset_precision_submatrix_difference']=float(np.max(abs(sub['precision']-data['precision'][np.ix_(mask,mask)])))
    checks['zHD_to_zHEL_convention']='Luminosity distance from background zHD is multiplied by(1+zHEL)/(1+zHD), matching released DA(zHD)*(1+zHD)*(1+zHEL).'
    record={'code_sha256':sha(__file__),'release_commit':json.loads((WORK/'release-commit.json').read_text())['sha'],
            'rows':len(frame),'survey_counts':frame.survey_family.value_counts().to_dict(),
            'IDSURVEY_counts':frame.IDSURVEY.value_counts().to_dict(),'checks':checks,
            'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [folder/'DES-Dovekie_HD.csv',folder/'DES-Dovekie_Metadata.csv',folder/'STAT+SYS.npz',folder/'STATONLY.npz',folder/'DES-Dovekie-SN_Likelihood.py']},
            'outputs':{str(p.relative_to(ROOT)):sha(p) for p in out.iterdir()},
            'interface':{'ledger':str(out/'dovekie-ledger.csv'),'total':str(out/'dovekie-total.npz'),'stat':str(out/'dovekie-stat.npz')},
            'limits':['Covariance already contains released calibration/selection/CC uncertainty; do not add quoted MUERR variance again.',
                      'Do not combine with Pantheon+ as independent: shared events/calibrations require explicit crosscovariance.',
                      'Physical_ID is exact within released author names; unknownIAU events retain survey-qualified identifier, not a fabricated global alias.',
                      'README1635DES claim differs from actual1623DES rows; actual data ordering/counts govern this adapter.']}
    (RESULTS/'dovekie-interface.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({'rows':len(frame),'counts':record['survey_counts'],'checks':checks},indent=2))
if __name__=='__main__':main()
