"""Export diagnostic figure, never use cosmology as a scoring criterion."""
from engine import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    f=pd.read_csv(OUT/'fixed_predictions_exact.csv.gz')
    v=json.loads((OUT/'validation_results.json').read_text())
    score=json.loads((OUT/'fit_summary.json').read_text())['heldout']['epoch']['folds']
    score=pd.DataFrame(score).sort_values('CID')
    fig,ax=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for band,col in zip('griz',['#2a9d8f','#e76f51','#457b9d','#805eac']):
        d=f[f.BAND==band]
        ax[0,0].scatter(d.phase,1000*d.sncosmo_flux_fracdiff,s=12,alpha=.6,label=band,color=col)
    ax[0,0].axhline(0,color='0.4',lw=.8)
    ax[0,0].set(xlabel='Rest-frame phase (day)',ylabel='Mean flux difference (parts per thousand)',title='Independent sncosmo minus SNANA: fixed parameters')
    ax[0,0].legend(ncol=4,frameon=False)
    for cid,d in f.groupby('CID'):
        ax[0,1].scatter(d.phase,d.sncosmo_model_magerr/d.snana_model_magerr,s=11,alpha=.55)
    ax[0,1].axhline(1,color='0.4',lw=.8)
    ax[0,1].set(xlabel='Rest-frame phase (day)',ylabel='Native sncosmo / SNANA model sigma',title='Same runtime maps; different uncertainty transport')
    profile=pd.DataFrame([p for p in v['t0_profiles'] if p['CID']=='1896213'])
    x=np.linspace(-1.2,1.2,201)
    ax[1,0].plot(x,(x/profile.independent_sigma.iloc[0])**2,label='Independent Hessian',color='#2a9d8f')
    ax[1,0].plot(x,(x/profile.snana_hessian_sigma.iloc[0])**2,label='SNANA local Hessian',color='#e76f51',ls='--')
    ax[1,0].scatter(profile.delta_t0,profile.delta_chi2_profile,label='Independent profile',color='k',s=20,zorder=4)
    ax[1,0].set(xlabel='Peak-time displacement (observer day)',ylabel='Profile Δχ²',ylim=(0,4.8),title='CID 1896213: local Hessian understates time uncertainty')
    ax[1,0].legend(frameon=False,fontsize=8)
    y=np.arange(len(score))
    ax[1,1].barh(y-.17,score.f99_logpdf_laplace_difference,height=.32,label='F99 law − SALT',color='#457b9d')
    ax[1,1].barh(y+.17,score.phase_colour_logpdf_laplace_difference,height=.32,label='Phase colour − SALT',color='#e9c46a')
    ax[1,1].set_yticks(y,score.CID,fontsize=8)
    ax[1,1].axvline(0,color='0.4',lw=.8)
    ax[1,1].set(xlabel='Held-out epoch Δ log predictive density',title='12-SN pilot: weak, heterogeneous mean-model evidence')
    ax[1,1].legend(frameon=False,fontsize=8)
    fig.savefig(OUT/'independent_flux_diagnostics.png',dpi=170)
    fig.savefig(OUT/'independent_flux_diagnostics.pdf')
    provenance('plot_results',[OUT/'independent_flux_diagnostics.png',OUT/'independent_flux_diagnostics.pdf'],[OUT/'fixed_predictions_exact.csv.gz',OUT/'validation_results.json',OUT/'fit_summary.json'])


if __name__=='__main__':main()
