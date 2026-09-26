"""Scientific figure for independently verified deterministic filter responses."""
from pathlib import Path
import json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'runs/research_2026_09_26/csp_native_filter_response/root-state-review'

def main():
    path=BASE/'response-decomposition/result.json';data=json.loads(path.read_text())
    verified=json.loads((BASE/'stable-response/result.json').read_text())
    assert data['pass_arithmetic'] and verified['pass_root_numerical_and_physical_row_checks']
    out=BASE/'figure';out.mkdir(exist_ok=False)
    rows=data['objects'];n=len(rows);assert n==42
    y=np.arange(n);native=np.array([r['total_delta_D']*1000 for r in rows])
    oldC=np.array([r['fixed_old_C_mean_response']*1000 for r in rows])
    affected=np.array([r['affected'] for r in rows]);mean=np.mean(native)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(figsize=(10,11.5));fig.subplots_adjust(left=.16,right=.97,top=.865,bottom=.16)
    fig.text(.16,.965,'CSP J-filter relabel: fitted distance response',ha='left',fontsize=18,weight='bold')
    fig.text(.16,.933,'42 supernovae · 188 labels changed · WIRC J mapped to the RC2 approximation',ha='left',fontsize=10,color='#444444')
    fig.text(.16,.906,f'Mean across all 42: {mean:+.3f} mmag     |     10 unchanged controls: exactly zero',fontsize=11)
    ax.axvline(0,color='#858585',lw=1)
    ax.axvline(mean,color='#ab3f38',lw=1.2,ls=(0,(4,3)))
    for i in range(n):
        if i%2==0:ax.axhspan(i-.5,i+.5,color='#f4f6f7',zorder=0)
        if affected[i]:
            ax.plot([oldC[i],native[i]],[i,i],color='#91aabd',lw=1.3,zorder=2)
    ax.scatter(oldC[affected],y[affected],s=32,facecolor='white',edgecolor='#5b768b',lw=1,zorder=3)
    ax.scatter(native[affected],y[affected],s=29,color='#166c91',zorder=4)
    ax.scatter(native[~affected],y[~affected],s=24,marker='D',color='#555555',zorder=4)
    ax.set_yticks(y,[r['CID'] for r in rows]);ax.set_ylim(n-.4,-.8)
    lim=max(np.max(abs(native)),np.max(abs(oldC)))+.6;ax.set_xlim(-max(4,lim*.65),lim)
    ax.set_xlabel('Change in distance modulus (mmag; 1 mmag = 0.001 mag)',labelpad=9)
    ax.tick_params(axis='y',length=0,pad=8,labelsize=9)
    ax.grid(axis='x',color='#dbe2e6',lw=.6);ax.set_axisbelow(True)
    handles=[Line2D([],[],marker='o',ls='',color='#166c91',label='Full native response'),
             Line2D([],[],marker='o',ls='',markerfacecolor='white',markeredgecolor='#5b768b',label='Changed mean, original covariance'),
             Line2D([],[],marker='D',ls='',color='#555555',label='Unchanged control')]
    fig.legend(handles=handles,loc='lower left',bbox_to_anchor=(.155,.065),ncol=3,frameon=False,fontsize=9,handletextpad=.5,columnspacing=1.4)
    fig.text(.16,.044,'Templates, peak times, stretch and extinction held fixed. Points are processing shifts, not uncertainty intervals.',fontsize=9,color='#444444')
    fig.text(.16,.025,'Bias correction, selection and training have not been regenerated. No corrected cosmology is inferred.',fontsize=9,color='#444444')
    for ext in ['png','pdf']:fig.savefig(out/f'csp-filter-response.{ext}',dpi=180,facecolor='white')
    plt.close(fig)
    (out/'executed-source.py').write_bytes(Path(__file__).read_bytes())
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    (out/'manifest.json').write_text(json.dumps(dict(inputs={str(path.relative_to(ROOT)):sha(path)},outputs={p.name:sha(p) for p in out.iterdir() if p.is_file() and p.name!='manifest.json'}),indent=2)+'\n')
    print(out/'csp-filter-response.png')

if __name__=='__main__':main()
