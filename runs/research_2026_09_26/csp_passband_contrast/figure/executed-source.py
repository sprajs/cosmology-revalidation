"""Descriptive figure of frozen template grid; derived WIRC−RC2 is post hoc."""
from pathlib import Path
import csv,json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'runs/research_2026_09_26/csp_passband_contrast'
OUT=BASE/'figure'
def main():
    OUT.mkdir(exist_ok=False)
    rows=list(csv.DictReader((BASE/'grid.csv').open()))
    values={(r['pair'],float(r['phase']),float(r['z'])):float(r['delta_mag_equal_BD17']) for r in rows}
    phases=sorted(set(k[1] for k in values));zs=sorted(set(k[2] for k in values))
    derived=[]
    for p in phases:
        for z in zs:
            d=values['J_WIRC_minus_RC1',p,z]-values['J_RC2_minus_RC1',p,z]
            d0=values['J_WIRC_minus_RC1',0,z]-values['J_RC2_minus_RC1',0,z]
            derived.append(dict(phase=p,z=z,delta_mag_equal_BD17=d,phase_difference=d-d0))
    fig,ax=plt.subplots(1,3,figsize=(11.4,4),sharey=True)
    for a,pair,title in zip(ax,['J_WIRC_minus_RC1','J_WIRC_minus_RC2','H_WIRC_minus_RetroCam'],['J: WIRC − RC1','J: WIRC − RC2','H: WIRC − RetroCam']):
        for z,color in zip([0,.03,.08],['#006d77','#b35c00','#6441a5']):
            if pair=='J_WIRC_minus_RC2':y=[next(r['phase_difference'] for r in derived if r['phase']==p and r['z']==z) for p in phases]
            else:y=[values[pair,p,z]-values[pair,0,z] for p in phases]
            a.plot(phases,np.array(y)*1000,marker='o',ms=3,color=color,label=f'z = {z:.2f}')
        a.axhline(0,color='#999999',lw=.7);a.set_title(title);a.set_xlabel('Rest-frame phase (days)');a.grid(alpha=.18)
    ax[0].set_ylabel('Contrast change from phase 0 (mmag)');ax[-1].legend(frameon=False,fontsize=9)
    fig.suptitle('Conditional response of the archived Hsiao spectrum to CSP filters',fontsize=13,y=.99)
    fig.text(.5,.025,'Template calculation, not measured bias. Phase subtraction removes constant zero-point offsets. WIRC−RC2 derived after the primary grid.',ha='center',fontsize=8.5)
    fig.tight_layout(rect=(0,.07,1,.95));fig.savefig(OUT/'phase-contrast.png',dpi=180);fig.savefig(OUT/'phase-contrast.pdf');plt.close(fig)
    with (OUT/'derived-wirc-rc2.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=derived[0]);w.writeheader();w.writerows(derived)
    s={k:[min(r[k] for r in derived),max(r[k] for r in derived)] for k in ['delta_mag_equal_BD17','phase_difference']}
    (OUT/'derived-result.json').write_text(json.dumps({'status':'Post-primary-grid algebraic contrast; not a newly frozen primary outcome.','WIRC_minus_RC2':s},indent=2)+'\n')
    (OUT/'executed-source.py').write_bytes(Path(__file__).read_bytes())
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    (OUT/'manifest.json').write_text(json.dumps({'parent_grid_sha256':sha(BASE/'grid.csv'),'files_sha256':{p.name:sha(p) for p in OUT.iterdir() if p.is_file() and p.name!='manifest.json'}},indent=2)+'\n')
    print(s)
if __name__=='__main__':main()
