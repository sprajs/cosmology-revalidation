#!/usr/bin/env python3
from pathlib import Path
import sys,json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).resolve().parent))
from analyse import ROOT,OUT,DER,save_manifest

def run():
    j=json.loads((OUT/'linmix-summary.json').read_text());results={q['name']:q for q in j};cuts=[.2,.25,.3,.35,.42];paper=np.array([-.034,-.027,-.026,-.025,-.024]);pe=np.array([.012,.013,.012,.012,.012]);comp=[]
    fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained');d=pd.read_csv(DER/'matched_G11_first.csv')
    for y,label,color in [('hr_corrected','Corrected Pantheon+','tab:blue'),('hr_no_bias','Bias correction reversed','tab:orange')]:
        axes[0].scatter(d.age,d[y]-d[y].mean(),s=10,alpha=.4,color=color)
        chain=np.load(OUT/'linmix'/f'G11_first_{y}_error_diag_0.42.npz')['chain'];xs=np.linspace(0,10,100);ys=chain['alpha'][:,None]+chain['beta'][:,None]*xs-d[y].mean();mean=ys.mean(axis=0);lo,hi=np.quantile(ys,[.16,.84],axis=0)
        axes[0].plot(xs,mean,color=color,label=label);axes[0].fill_between(xs,lo,hi,color=color,alpha=.15)
    axes[0].set(xlabel='Inferred host age (Gyr)',ylabel='Mean-centered Hubble residual (mag)',ylim=(-.6,.6),title='196 unique SNe, same ages and redshift range');axes[0].legend(fontsize=8)
    axes[1].errorbar(cuts,paper,pe,fmt='o--',color='black',label='C26 published',capsize=2)
    for pol,color,offset in [('G11_first','tab:blue',-.002),('R19_first','tab:orange',.002)]:
        mean=[];sd=[]
        for cut,p,e in zip(cuts,paper,pe):
            r=results[f'{pol}_hr_tripp_error_no_covadd_{cut}'];mean.append(r['mean']);sd.append(r['sd']);comp.append({'age_policy':pol,'zcut':cut,'n':r['n'],'published_slope':p,'published_sd':e,'our_slope':r['mean'],'our_sd':r['sd'],'status':'approximate: author age overlap/uncertainty settings unavailable'})
        axes[1].errorbar(np.array(cuts)+offset,mean,sd,fmt='s-',color=color,capsize=2,label=pol.replace('_',' '))
    axes[1].axhline(0,color='grey',lw=.8);axes[1].set(xlabel='Maximum redshift',ylabel='Age slope (mag/Gyr)',title='Uncorrected variants; cuts share observations');axes[1].legend(fontsize=8)
    fig.savefig(OUT/'matched-age-summary.png',dpi=170);plt.close(fig);csv=OUT/'published-cut-comparison.csv';pd.DataFrame(comp).to_csv(csv,index=False)
    save_manifest('summary',[OUT/'linmix-summary.json',DER/'matched_G11_first.csv',Path(__file__)],{'published_source':'doi:10.1093/mnras/stag1513 published Figure 2','plot_note':'Left bands posterior16–84%; no uncertainty bars on individual ages for legibility. Curves do not measure progenitor ages.'},[csv,OUT/'matched-age-summary.png'])
if __name__=='__main__':run()
