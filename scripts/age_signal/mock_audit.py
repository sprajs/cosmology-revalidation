#!/usr/bin/env python3
"""Paired controlled surrogate tests of the C26 bin-centering operation."""
import sys,json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).resolve().parent))
from analyse import OUT, ROOT, digest, save_manifest, mu
SEED=20260920; NREP=1000; NPER=200; B=-.034; NOISE=.12
ZCENTERS=np.arange(.025,.4,.05); CUTS=[.2,.25,.3,.35,.42]
def slope(x,y):
    xc=x-x.mean(axis=-1,keepdims=True); yc=y-y.mean(axis=-1,keepdims=True)
    return np.sum(xc*yc,axis=-1)/np.sum(xc*xc,axis=-1)
def run():
    rng=np.random.default_rng(SEED); records=[]; identities=[]
    z=np.repeat(ZCENTERS,NPER); dm=mu(z)
    for gradient in [0.,5.,10.,15.]:
        age=7-gradient*z+rng.normal(0,1.5,size=(NREP,len(z)))
        age=np.clip(age,.05,13.5)
        epsilon=rng.normal(0,NOISE,age.shape)
        y=B*(age-5)+epsilon
        # Apparent magnitudes generated and baseline subtracted explicitly.
        apparent=dm[None,:]+(-19.3)+y
        hr=apparent-dm[None,:]+19.3
        assert np.allclose(hr,y,atol=1e-12)
        xbin=age.reshape(NREP,len(ZCENTERS),NPER)
        ybin=hr.reshape(NREP,len(ZCENTERS),NPER)
        xc=(xbin-xbin.mean(axis=-1,keepdims=True)).reshape(NREP,-1)
        yc=(ybin-ybin.mean(axis=-1,keepdims=True)).reshape(NREP,-1)
        for cut in CUTS:
            m=z<cut; x=age[:,m]; yy=hr[:,m]; xx=xc[:,m]; cc=yc[:,m]
            # exact finite-sample prediction for noiseless bin centering
            attenuation=np.sum(xx*xx,axis=1)/np.sum((x-x.mean(axis=1,keepdims=True))**2,axis=1)
            clean=(B*(xbin-xbin.mean(axis=-1,keepdims=True))).reshape(NREP,-1)[:,m]
            actual=slope(x,clean)
            identities.append({'gradient':gradient,'cut':cut,'max_identity_error':float(np.max(np.abs(actual-B*attenuation))),'mean_attenuation':float(attenuation.mean())})
            for case,xv,yv in [('physical_LCDM',x,yy),('center_y_only',x,cc),('within_bins',xx,cc),('null',x,epsilon[:,m])]:
                beta=slope(xv,yv)
                expected=0 if case=='null' else B
                err=yv-(yv.mean(axis=1,keepdims=True)+beta[:,None]*(xv-xv.mean(axis=1,keepdims=True)))
                se=np.sqrt(np.sum(err*err,axis=1)/(len(xv[0])-2)/np.sum((xv-xv.mean(axis=1,keepdims=True))**2,axis=1))
                records.append({'gradient_Gyr_per_z':gradient,'zcut':cut,'case':case,'n':int(m.sum()),'mean_slope':beta.mean(),'sd_slope':beta.std(ddof=1),'mean_mc_se':beta.std(ddof=1)/np.sqrt(NREP),'q025':np.quantile(beta,.025),'q975':np.quantile(beta,.975),'coverage_true_95':np.mean(np.abs(beta-expected)<1.96*se),'reject_zero_95':np.mean(np.abs(beta)>1.96*se)})
    df=pd.DataFrame(records);df.to_csv(OUT/'mock-results.csv',index=False)
    (OUT/'mock-identities.json').write_text(json.dumps(identities,indent=2))
    fig,ax=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    for case,label in [('physical_LCDM','Fixed true cosmology'),('center_y_only','Residuals bin-centered'),('within_bins','Ages and residuals bin-centered')]:
        q=df[(df.gradient_Gyr_per_z==10)&(df.case==case)];ax[0].errorbar(q.zcut,q.mean_slope,yerr=q.sd_slope,marker='o',capsize=2,label=label)
    ax[0].axhline(B,color='black',ls='--',lw=1);ax[0].set(xlabel='Maximum redshift',ylabel='Recovered slope (mag/Gyr)',title='Synthetic populations; error bars = trial scatter');ax[0].legend(fontsize=8)
    for case,label in [('physical_LCDM','Fixed true cosmology'),('center_y_only','Residuals bin-centered'),('within_bins','Within-bin slope')]:
        q=df[(df.zcut==.42)&(df.case==case)];ax[1].plot(q.gradient_Gyr_per_z,q.mean_slope,'o-',label=label)
    ax[1].axhline(B,color='black',ls='--',lw=1);ax[1].set(xlabel='Imposed age evolution (Gyr per redshift)',ylabel='Recovered slope (mag/Gyr)',title='Full sample: sensitivity to age evolution')
    fig.savefig(OUT/'mock-audit.png',dpi=170);plt.close(fig)
    # save_manifest is imported; include this actual script explicitly as input as well.
    save_manifest('A4-mock-audit',[Path(__file__)],{'seed':SEED,'nrep':NREP,'n_per_bin':NPER,'zcenters':ZCENTERS.tolist(),'b':B,'noise_mag':NOISE,'age_sigma_Gyr':1.5,'age_mean':'7-gradient*z','age_clip':[.05,13.5],'gradients':[0,5,10,15],'status':'controlled surrogate, not exact author mock'},[OUT/'mock-results.csv',OUT/'mock-identities.json',OUT/'mock-audit.png'])
    print(df[(df.gradient_Gyr_per_z==10)&(df.zcut.isin([.2,.42]))].to_string(index=False))
if __name__=='__main__':run()
