#!/usr/bin/env python3
"""Publication-style evidence figure; derived plot rows stay outside Git."""
import json,hashlib
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[4];WORK=ROOT/'.work/physical-ages'

def main():
    source=ROOT/'.work/galaxy-validation/host-spectrum-age-ledger.csv'
    d=pd.read_csv(source,dtype={'specObjID':str,'physical_host':str})
    d=d[(d['sample']=='ZTF')&d.index_quality&np.isfinite(d.sed_age)].drop_duplicates('physical_host')
    assert len(d)==609
    spatial=ROOT/'.work/galaxy-validation/manga-central-SN-indices.csv'
    s=pd.read_csv(spatial);s=s[s['sample']=='ZTF'].drop_duplicates('MANGAID');s=s[s.central_valid&s.SN_site_valid]
    assert len(s)==34
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.labelcolor':'#292d33','text.color':'#292d33','xtick.color':'#292d33','ytick.color':'#292d33','figure.facecolor':'white','axes.facecolor':'white'})
    fig,ax=plt.subplots(1,3,figsize=(13.6,4.3),layout='constrained',gridspec_kw={'width_ratios':[1,1,1.2]})
    blue='#256797';orange='#b85b2d'
    from scipy.stats import spearmanr
    assert abs(spearmanr(d.sed_age,d.d4000_n).statistic-.659)<.0005
    ax[0].scatter(d.sed_age,d.d4000_n,s=12,c=blue,alpha=.43,edgecolors='none')
    ax[0].set(xlabel='Published global age estimate (Gyr)',ylabel=r'Observed narrow 4000 Å break, $D_n$',xlim=(2,10.5),ylim=(.85,2.5))
    ax[0].set_title('A   Age ordering is observed',loc='left',fontweight='bold',pad=13)
    ax[0].text(.04,.97,'609 distinct hosts\nSpearman ρ = 0.659 [0.605, 0.706]',transform=ax[0].transAxes,va='top',fontsize=9)
    lim=(.85,2.3);ax[1].plot(lim,lim,':',color='#444444',lw=1.2)
    ax[1].errorbar(s.central_Dn4000,s.SN_site_Dn4000,xerr=s.central_Dn4000_error,yerr=s.SN_site_Dn4000_error,fmt='o',ms=4,color=blue,ecolor='#bbbbbb',elinewidth=.8)
    ax[1].set(xlabel=r'Galaxy-centre $D_n$',ylabel=r'Supernova-site $D_n$',xlim=lim,ylim=lim)
    ax[1].set_title('B   Location matters',loc='left',fontweight='bold',pad=13)
    ax[1].text(.05,.97,'34 distinct galaxies\nMedian site − centre = −0.123',transform=ax[1].transAxes,va='top',fontsize=9)
    rows=[]
    for i,(sample,label) in enumerate([('sdss','SDSS fibres'),('roman','Roman local apertures'),('des','DES deep hosts')]):
        p=ROOT/f'studies/host_ages/results/physical_ages/{sample}-age-bounds.json';j=json.loads(p.read_text());r=j['summaries'][0];width=r['age_interval_width_Gyr']
        q=[width['p05'],width['median'],width['p95']];rows.append(dict(sample=sample,low=q[0],median=q[1],high=q[2],admitted=r['photometry_feasible'],total=r['objects']))
        ax[2].errorbar(q[1],2-i,xerr=[[q[1]-q[0]],[q[2]-q[1]]],fmt='o',color=orange,capsize=4,lw=2,ms=5)
        ax[2].text(.15,2-i+.17,f'{label}: {r["photometry_feasible"]}/{r["objects"]} compatible',fontsize=9)
    ax[2].set(xlabel='Width of conditional age region (Gyr)',yticks=[],xlim=(0,14),ylim=(-.5,2.55))
    ax[2].set_title('C   Absolute ages remain broad',loc='left',fontweight='bold',pad=13)
    ax[2].text(.01,-.24,'Dots: median width; bars: 5th–95th object percentiles.\nDifferent redshifts impose different cosmic-age ceilings.',transform=ax[2].transAxes,fontsize=8.8)
    for a in ax:a.grid(alpha=.15);a.set_axisbelow(True)
    fig.suptitle('Stellar-population evidence and the limits of host-age inference',fontsize=14,fontweight='bold')
    output=ROOT/'docs/figures/host-age-physics.png';fig.savefig(output,dpi=180,bbox_inches='tight');plt.close(fig)
    pd.DataFrame(rows).to_csv(WORK/'figure-age-widths.csv',index=False)
    record={'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [source,spatial]+[ROOT/f'studies/host_ages/results/physical_ages/{s}-age-bounds.json' for s in ['sdss','roman','des']]},'figure':str(output.relative_to(ROOT)),'figure_sha256':hashlib.sha256(output.read_bytes()).hexdigest(),'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'A: catalogue indices and photometric age summaries. B: paired observed spatial indices with released measurement errors. C: statistical-only, model-conditional95% age-region widths; displayed whiskers are object percentiles, not uncertainties on medians. No aggregate cosmological correction is inferred.'}
    (ROOT/'studies/host_ages/results/physical_ages/figure.json').write_text(json.dumps(record,indent=2)+'\n')

if __name__=='__main__':main()
