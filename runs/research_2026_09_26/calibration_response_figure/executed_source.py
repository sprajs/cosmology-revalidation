"""Scientific response figure; hypothetical offsets are explicitly labelled."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]


def main():
    out=ROOT/'runs/research_2026_09_26/calibration_response_figure'
    out.mkdir(exist_ok=False)
    src=Path(__file__).read_bytes();(out/'executed_source.py').write_bytes(src)
    paths=[ROOT/'runs/research_2026_09_26/astra_design/exact43/comparison/fitted-distance-response.csv',
        ROOT/'runs/research_2026_09_26/calibration_pattern_refit/paired_refits.csv',
        ROOT/'runs/research_2026_09_26/calibration_pattern_refit/summary.json']
    lin=pd.read_csv(paths[0]);lin=lin[lin.arm=='official_mean_exactC']
    full=pd.read_csv(paths[1]);mag=json.loads(paths[2].read_text())['delta_observed_magnitude_griz']
    plt.rcParams.update({'font.size':12})
    fig,axes=plt.subplots(1,2,figsize=(11.8,4.8),gridspec_kw={'width_ratios':[1,2]})
    ax=axes[0];ax.bar(list('griz'),np.array(mag)*1000,color=['#64748b','#64748b','#64748b','#64748b'])
    ax.axhline(0,color='black',lw=.8);ax.set_ylabel('Added observed magnitude (mmag)')
    ax.set_xlabel('DES filter');ax.set_title('Hypothetical filter correction')
    ax=axes[1];ax.scatter(lin.zHEL,lin.fixed_reference_preBBC,s=25,color='#94a3b8',label='Local response: 43 pilot objects')
    ax.scatter(full.zHEL,full.nonlinear_delta_standardized,s=50,color='#b45309',marker='D',label='Nonlinear refits: 6 checks',zorder=3)
    ax.axhline(0,color='black',lw=.8);ax.set_xlabel('Redshift');ax.set_ylabel('Pre-BBC standardized-distance change (mag)')
    ax.set_title('Distance response with fixed training and α, β')
    high=full.sort_values('zHEL').iloc[-1]
    ax.annotate('Δχ² = −0.008 for this refit',xy=(high.zHEL,high.nonlinear_delta_standardized),
                xytext=(.35,.145),arrowprops={'arrowstyle':'->','color':'#475569'},fontsize=10)
    ax.legend(loc='upper left',fontsize=10,frameon=False)
    fig.suptitle('Similar light-curve fits can imply different distances',fontsize=16)
    fig.text(.5,.015,'A response experiment, not an estimated physical bias. Training, selection and BBC are not rerun.',ha='center',fontsize=10)
    fig.tight_layout(rect=[0,.04,1,.94])
    for ext in ['png','pdf']:fig.savefig(out/f'calibration-response.{ext}',dpi=180)
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    (out/'manifest.json').write_text(json.dumps({'source_sha256':hashlib.sha256(src).hexdigest(),
        'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},
        'outputs_sha256':{p.name:sha(p) for p in out.glob('*') if p.is_file()}},indent=2)+'\n')


if __name__=='__main__':main()
