"""Describe saved restricted profiles; thresholds are not confidence levels."""
from pathlib import Path
import csv,json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path(__file__).resolve().parent;G=O/'fixed-c-profile/global-profile'
sha=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()

def main():
    rows=list(csv.DictReader((G/'shape-profiles.csv').open()));res=json.loads((G/'result.json').read_text());out={};curves={}
    for label in res['metrics']:
      rr=[{k:(v if k in ['metric','stage'] else float(v)) for k,v in x.items()} for x in rows if x['metric']==label and x['stage']=='fine']
      q0=res['metrics'][label]['finest']['Q'];mi=[]
      for j in range(1,len(rr)-1):
        if rr[j]['Q']<=rr[j-1]['Q'] and rr[j]['Q']<=rr[j+1]['Q']:
          mi.append({k:rr[j][k] for k in ['shape','AV','Q','DLMAG']}|{'delta_Q':rr[j]['Q']-q0})
      ranges={}
      for threshold in [.1,1.,4.,9.]:
        keep=[x for x in rr if x['Q']<=q0+threshold]
        ranges[str(threshold)]={k:[min(x[k] for x in keep),max(x[k] for x in keep)] for k in ['shape','AV','DLMAG']}
      out[label]=dict(fine_grid_local_minima=mi,profiled_shape_ridge_ranges_by_deltaQ=ranges,
        warning='Ranges follow AV/amplitude-profiled shape ridge; not full distance-profile confidence or posterior intervals.')
      curves[label]=rr
    (G/'branch-geometry.json').write_text(json.dumps(out,indent=2)+'\n')
    fig,ax=plt.subplots(2,1,figsize=(7.3,6.2),sharex=True,layout='constrained')
    for arm,color,name in [('A','#33658a','A: 65 positive epochs'),('B','#ba5a31','B: 74 signed epochs')]:
      label='Banchor_'+arm;rr=curves[label];s=np.array([r['shape'] for r in rr]);q=np.array([r['Q'] for r in rr]);d=np.array([r['DLMAG'] for r in rr]);best=res['metrics'][label]['finest']
      ax[0].plot(s,q-best['Q'],c=color,label=name,lw=1.8)
      ax[1].plot(s,d,c=color,lw=1.8)
      ax[0].scatter([best['shape']],[0],c=color,s=32,zorder=4)
      ax[1].scatter([best['shape']],[best['DLMAG']],c=color,s=32,zorder=4)
    ax[0].set(ylim=(-.07,4.1),ylabel='Quadratic above each arm minimum')
    ax[1].set(xlim=(.9,1.3),xlabel='SNooPy stretch, fixed header peak',ylabel='DLMAG along profiled shape ridge')
    ax[0].legend(frameon=False,loc='upper left')
    ax[0].annotate('Signed-data competing branch\nΔQ = 0.116; DLMAG = 42.768',xy=(1.18923,.11574),xytext=(1.08,2.8),fontsize=9,arrowprops={'arrowstyle':'->','color':'0.35'})
    for a in ax:a.grid(alpha=.18);a.spines[['top','right']].set_visible(False)
    fig.suptitle('DES16C1cim: stable minima, weakly separated distance branches',fontsize=12)
    fig.savefig(G/'profile-branches.png',dpi=180);fig.savefig(G/'profile-branches.pdf');plt.close(fig)
    manifest=dict(source_sha256=sha(__file__),input_hashes={f:sha(G/f) for f in ['shape-profiles.csv','result.json']},
      output_hashes={f:sha(G/f) for f in ['branch-geometry.json','profile-branches.png','profile-branches.pdf']})
    (G/'branch-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')

if __name__=='__main__':main()
