from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).parent
rows=list(csv.DictReader((P/'paired-object-results.csv').open()));summ=list(csv.DictReader((P/'summary.csv').open()))
fig,ax=plt.subplots(1,2,figsize=(11.4,5.6),gridspec_kw={'width_ratios':[1.25,1]})
cols=['#ac2c3d','#227d70'];names=sorted({r['name'] for r in rows if r['pair']=='J_WIRC_minus_RC1' and r['target']=='observed' and r['clock_sigma_shift']=='0'})
for i,name in enumerate(names):
 vals=[1000*float(next(r for r in rows if r['name']==name and r['pair']==pair and r['target']=='observed' and r['clock_sigma_shift']=='0')['delta_mag_equal_BD17']) for pair in ['J_WIRC_minus_RC1','J_WIRC_minus_RC2']]
 ax[0].plot(vals,[i,i],color='.7',lw=1,zorder=1)
 for j,v in enumerate(vals):ax[0].scatter(v,i,c=cols[j],s=28,label=['WIRC − RC1','WIRC − RC2'][j] if i==0 else None,zorder=2)
ax[0].axvline(0,c='.4',lw=.8);ax[0].set_yticks(range(len(names)),names,fontsize=8);ax[0].invert_yaxis();ax[0].set_xlabel('Late − early change in passband contrast [mmag]');ax[0].set_title('Measured spectra at each object’s redshift',fontsize=11);fig.legend(*ax[0].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.30,.94),ncol=2,frameon=False,fontsize=9)
for pair,c,label in [('J_WIRC_minus_RC1',cols[0],'J: WIRC − RC1 (13 objects)'),('J_WIRC_minus_RC2',cols[1],'J: WIRC − RC2 (13)'),('H_WIRC_minus_RetroCam','#5c5796','H (12)'),('Y_WIRC_minus_RetroCam','#9f7725','Y (17)')]:
 rr=[r for r in summ if r['pair']==pair and r['target']!='observed' and r['clock_sigma_shift']=='0' and r['exclude_2012fr']=='False'];rr.sort(key=lambda r:float(r['target']));ax[1].plot([float(r['target']) for r in rr],[1000*float(r['equal_object_mean']) for r in rr],'-o',c=c,label=label,ms=4)
ax[1].axhline(0,c='.4',lw=.8);ax[1].set_xlabel('Synthetic target redshift');ax[1].set_ylabel('Equal-object mean phase contrast [mmag]');ax[1].set_title('Same delivered spectral shapes, redshifted',fontsize=11);ax[1].legend(fontsize=8)
for a in ax:a.grid(alpha=.15)
fig.suptitle('CSP-II FIRE spectra: early [−7, 7] to late [10, 20] rest-frame days',fontsize=13)
fig.text(.5,.028,'Descriptive contrasts; no uncertainty bars or cosmological inference. Static zero points and gray scaling cancel.\nTelluric, calibration, timing and population-transfer limits remain; target-z transport is not an observed high-z sample.',ha='center',fontsize=9)
fig.tight_layout(rect=[0,.1,1,.94]);fig.savefig(P/'spectral-phase-contrast.png',dpi=160);fig.savefig(P/'spectral-phase-contrast.pdf')
