"""Build source-linked numerical summaries and standalone research figures."""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from core import ROOT,manifest,mu

OUT=ROOT/'runs/cosmology/summary';OUT.mkdir(parents=True,exist_ok=True)
rows=[];inputs=[]
for p in sorted((ROOT/'runs/cosmology').glob('*/summary.json')):
    d=json.loads(p.read_text())
    if 'q0' not in d:continue
    c=json.loads((p.parent/'configuration.json').read_text())
    row={'experiment':p.parent.name,'model':c['model'],'data':c['data'],'amplitude_prior':c['amplitude'],
         'q0_mean':d['q0']['mean'],'q0_sd':d['q0']['sd'],'q0_q16':d['q0']['q16'],'q0_q84':d['q0']['q84'],
         'P_q0_negative':d['P_q0_negative'],'P_block_mcse':d['P_q0_negative_block_mcse'],
         'minimum_steps_per_tau':min(d['retained_steps_per_tau']),'minimum_ESS':min(d['ESS_per_parameter']),
         'chi2_at_posterior_mode':d['posterior_mode_chisq'],'sn_chi2_at_mode':d.get('mode_sn_chisq'),
         'bao_chi2_at_mode':d.get('mode_bao_chisq'),'max_boundary_fraction':max(d['prior_boundary_fraction'].values())}
    for k,v in d['parameters'].items():row[k+'_mean']=v['mean'];row[k+'_sd']=v['sd']
    rows.append(row);inputs.extend([p,p.parent/'configuration.json'])
df=pd.DataFrame(rows);df.to_csv(OUT/'all_fits.csv',index=False)

plt.rcParams.update({'font.size':10,'font.family':'DejaVu Sans','axes.spines.top':False,'axes.spines.right':False})
palette=['#225ea8','#d17a13','#7b5293','#497c59']
chosen=[('pantheon-bao-cpl','Original corrected distances'),
        ('pantheon-bao-cpl-hostmass','114-row host-mass revision'),
        ('pantheon-bao-cpl-shortdelay','Short-delay DTD, mean'),
        ('pantheon-bao-cpl-c14mean','C14 DTD, mean'),
        ('pantheon-bao-cpl-c14fixed','C14 DTD, median, fixed slope'),
        ('pantheon-bao-cpl-c14slope','C14 median, shared slope uncertainty'),
        ('pantheon-bao-cpl-c14free','C14 median, broad amplitude prior')]
fig,ax=plt.subplots(figsize=(11,6.3))
fig.subplots_adjust(left=.43,right=.98,bottom=.18,top=.88)
for i,(key,label) in enumerate(chosen):
    r=df[df.experiment==key]
    if r.empty:continue
    r=r.iloc[0];x=r.q0_mean
    ax.errorbar(x,i,xerr=[[x-r.q0_q16],[r.q0_q84-x]],fmt='o',capsize=4,
                color=palette[0] if i<2 else palette[1],ms=5)
    probability_label='>0.999*' if r.P_q0_negative>=.9995 else f'{r.P_q0_negative:.3f}'
    ax.text(.37,i,probability_label,va='center',fontsize=9)
ax.axvline(0,color='#444444',lw=1,ls='--');ax.set_yticks(range(len(chosen)),[x[1] for x in chosen]);ax.invert_yaxis()
ax.set_xlabel('Present deceleration parameter q₀ (negative = acceleration)')
ax.set_title('Same Pantheon+ and DESI BAO data; different correction assumptions')
ax.text(.37,-.65,'P(q₀ < 0)',fontsize=9);ax.set_xlim(-.75,.55);ax.grid(axis='x',alpha=.2)
fig.text(.02,.015,'Dots: posterior means; bars: central 68% intervals. Flat CPL; no CMB. Age templates use fixed CPL63.6/B13 clocks.\nThe broad amplitude prior is a sensitivity test. * No opposite-sign draws in these runs; not a precisely resolved tail probability.',fontsize=8)
fig.savefig(OUT/'correction-sensitivity.png',dpi=180,bbox_inches='tight');plt.close(fig)

fig,axs=plt.subplots(1,2,figsize=(11,5.4))
fig.subplots_adjust(left=.07,right=.98,bottom=.21,top=.88,wspace=.23)
for i,(name,label) in enumerate([('pantheon-cpl','CPL, original distances'),('pantheon-cpl-c14fixed','CPL, fixed median correction'),
                                ('pantheon-qbins-tau15','q bins, original distances'),('pantheon-qbins-c14fixed','q bins, fixed median correction')]):
    p=ROOT/'runs/cosmology'/name/'q_history.csv';inputs.append(p);d=pd.read_csv(p)
    ax=axs[i//2];color=palette[i%2]
    ax.plot(d.z,d['median'],color=color,label=label,lw=1.5)
    ax.fill_between(d.z,d.q16,d.q84,color=color,alpha=.18)
for ax in axs:
    ax.axhline(0,color='#444444',ls='--',lw=1);ax.set_xlim(0,1.1)
    ax.set_xlabel('Redshift z');ax.set_ylabel('q(z)');ax.legend(fontsize=8,loc='upper right');ax.grid(alpha=.15)
axs[0].set_title('SN-only, flat CPL');axs[1].set_title('SN-only, five q bins; smoothness scale 1.5')
axs[0].set_ylim(-.7,.6);axs[1].set_ylim(-1.7,2.1)
fig.text(.02,.025,'Lines: posterior medians; bands: central 68% intervals. Panels have different vertical scales. Fixed correction is a conditional sensitivity.\nThe first q bin covers 0 < z < 0.1; it is not a direct point measurement of q(0).',fontsize=8)
fig.savefig(OUT/'q-history-model-sensitivity.png',dpi=180,bbox_inches='tight');plt.close(fig)

manifest(OUT,'Comparisons of independently fitted conditional distance likelihoods',inputs,
         {'interval':'central 68%','plot_probabilities':'posterior sample fractions, finite MC precision'},
         [OUT/'all_fits.csv',OUT/'correction-sensitivity.png',OUT/'q-history-model-sensitivity.png'])
print(df[['experiment','q0_mean','q0_sd','P_q0_negative','minimum_steps_per_tau']].to_string(index=False))
