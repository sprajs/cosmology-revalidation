"""Draw only qualified late-time posterior intervals, retaining two CPL seeds."""
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[4]
RESULTS=ROOT/'studies/unified_cosmology/results/inference'


def main():
    filenames=['late-lcdm-none-dovekie.json','late-wcdm-none-dovekie.json',
               'late-nested-refined-927088.json','late-nested-refined-927089.json']
    data=[json.loads((RESULTS/name).read_text()) for name in filenames]
    assert all(d['status']=='passed' for d in data[:2])
    comparison=RESULTS/'late-nested-refined-comparison.json'
    assert json.loads(comparison.read_text())['status']=='passed_independent_seed_diagnostics'
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,
                         'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(10,4.7),sharey=True)
    colors=['#225f9b','#b66719','#278175','#278175']
    positions=[2.,1.,.035,-.035]
    for axis,key,title in zip(axes,['q0','j0'],['Present acceleration','Change in scale-factor acceleration']):
        axis.axvline(0,color='#606060',lw=1,ls='--',zorder=0)
        for row,y,color in zip(data,positions,colors):
            lo,_,middle,_,hi=row['posterior'][key]['quantiles_025_16_50_84_975']
            axis.plot([lo,hi],[y,y],color=color,lw=2)
            axis.scatter(middle,y,s=37,color=color,zorder=3)
        axis.set_title(title,pad=14,fontsize=12)
        axis.set_xlabel('$q_0$' if key=='q0' else '$j_0$')
        axis.set_yticks([0,1,2],['CPL','Constant $w$','Flat ΛCDM'])
        axis.set_ylim(-.5,2.5)
        axis.grid(axis='x',alpha=.15)
        axis.spines['left'].set_visible(False)
        axis.tick_params(axis='y',length=0)
    axes[0].set_xlim(-.65,.04)
    axes[1].set_xlim(-1.3,1.8)
    fig.suptitle('The same distances under three expansion models',fontsize=16,y=.99)
    fig.text(.5,.90,'Dovekie supernovae + DESI DR2 BAO · no CMB or additional age correction',ha='center',fontsize=11,color='#444444')
    fig.text(.27,.16,'Negative q: accelerating expansion',ha='center',fontsize=10)
    fig.text(.75,.16,'Negative j: acceleration decreases with time',ha='center',fontsize=10)
    fig.text(.5,.065,'Markers: posterior medians. Lines: 95% equal-tail intervals. ΛCDM fixes j = 1.\nCPL shows two independent calculations; its remote tails remain less precise.',ha='center',fontsize=9,color='#444444')
    fig.subplots_adjust(top=.78,bottom=.29,left=.13,right=.97,wspace=.27)
    target=ROOT/'studies/unified_cosmology/figures/late-expansion.png'
    target.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(target,dpi=180,facecolor='white');plt.close(fig)
    inputs=[RESULTS/name for name in filenames]+[comparison,Path(__file__)]
    record={'status':'generated_from_qualified_intervals','figure':str(target.relative_to(ROOT)),
            'sha256':hashlib.sha256(target.read_bytes()).hexdigest(),
            'inputs_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
            'scope':'Late-time released distances; no CMB. Separate CPL runs retained, no tail significance inferred.'}
    (RESULTS/'late-figure.json').write_text(json.dumps(record,indent=2)+'\n')


if __name__=='__main__':main()
