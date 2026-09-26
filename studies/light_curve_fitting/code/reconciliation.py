"""E05 screening: SN luminosity correction required by BAO-only histories.

Histories were fitted without these SN magnitudes. Template amplitudes here are
conditional descriptive regressions, not independent measurements of age physics.
"""
from pathlib import Path
import sys
import json
import hashlib
from datetime import datetime,timezone

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize_scalar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from bao_recent_counterexample import geometry
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts/cosmology'))
from core import Pantheon,mu

SOURCE=Path(__file__).read_bytes()
CORE=(ROOT/'scripts/cosmology/core.py').read_bytes()
GEOMETRY=(Path(__file__).with_name('bao_recent_counterexample.py')).read_bytes()


def main():
    out=ROOT/'runs/research_2026_09_26/reconciliation'
    if out.exists():raise FileExistsError(out)
    out.mkdir()
    for name,b in [('executed_source.py',SOURCE),('imported_core.py',CORE),('imported_geometry.py',GEOMETRY)]:
        (out/name).write_bytes(b)
    sn=Pantheon()
    hp=ROOT/'runs/research_2026_09_26/bao_shape/recent_counterexample/results.json'
    histories=json.loads(hp.read_text())['fits']
    tp=[ROOT/'runs/cosmology/templates'/f for f in ['c14-cpl63-median.csv','c14-cpl63-mean.csv','cut40-cpl63-mean.csv']]
    templates={p.stem:PchipInterpolator(pd.read_csv(p).z,pd.read_csv(p).delta_mu)(sn.z) for p in tp}
    inv=np.linalg.inv(sn.cov);u=inv.sum(axis=1)
    center=lambda v:v-float(u@v/u.sum())
    lcdm=minimize_scalar(lambda om:sn.chisq_full([om],'lcdm')[0],bounds=(.05,.8),method='bounded',options={'xatol':1e-11})
    ref=mu(sn.z,[lcdm.x],'lcdm',sn.zhel)[0]
    fitrows=[];vectors=[]
    fig,axs=plt.subplots(1,2,figsize=(11,4.8),layout='constrained')
    colors={'main__free':'tab:purple','main__recent_q_plus_half':'tab:orange','main__always_nonnegative':'tab:gray'}
    for name in ['main__free','main__recent_nonnegative','main__recent_q_plus_half','main__always_nonnegative',
                 'coarser__recent_q_plus_half']:
        h=histories[name];pars=np.array(h['best']['x']);dm,_=geometry(sn.z,pars[0],pars[1:],np.array(h['edges']))
        mh=5*np.log10((1+sn.zhel)*dm)
        # m=M+mu+deltaM, so deltaM required to make history h match reference is mu_ref-mu_h.
        required=center(ref-mh)
        residual=sn.mag-mh
        chi0=float(residual@sn.A@residual)
        base={'history':name,'BAO_only_chi2':h['best']['chi2'],'SN_conditional_chi2_no_evolution':chi0,
              'required_shape_min_mag':float(required.min()),'required_shape_max_mag':float(required.max()),
              'required_shape_peak_to_peak_mag':float(np.ptp(required))}
        for tn,t in templates.items():
            info=float(t@sn.A@t); amp=float(t@sn.A@residual/info)
            delta=residual-amp*t;chi=float(delta@sn.A@delta)
            shapeamp=float(t@sn.A@required/info)
            shaperes=required-shapeamp*center(t)
            fitrows.append({**base,'template':tn,'fitted_amplitude_fixed_BAO_history':amp,
                'conditional_amplitude_sd':float(1/np.sqrt(info)),
                'equivalent_fixed_template_slope_mag_per_Gyr':.03*amp,
                'SN_conditional_chi2_after_template':chi,
                'required_shape_template_projection':shapeamp,
                'remaining_shape_chi2_metric':float(shaperes@sn.A@shaperes)})
        for i in range(len(sn.z)):
            vectors.append({'history':name,'CID':str(sn.df.CID.iloc[i]),'IDSURVEY':int(sn.df.IDSURVEY.iloc[i]),
                            'original_row':int(sn.original_indices[i]),'zHD':sn.z[i],
                            'required_relative_luminosity_mag':required[i],
                            'observed_residual_centered_mag':center(residual)[i]})
        if name in colors:
            idx=np.argsort(sn.z);label=name.replace('main__','').replace('_',' ')
            axs[0].plot(sn.z[idx],required[idx],label=label,color=colors[name])
            q=np.array(pars[1:]);e=h['edges']
            axs[1].stairs(q,e,label=label,color=colors[name])
    pd.DataFrame(fitrows).to_csv(out/'template_screen.csv',index=False)
    pd.DataFrame(vectors).to_csv(out/'required_vectors.csv.gz',index=False)
    axs[0].axhline(0,color='black',lw=.5);axs[0].set(xlabel='Redshift z',ylabel='Required relative luminosity term (mag)',
        title='Correction to match fitted SN ΛCDM shape',xlim=(.01,1.5))
    axs[1].axhline(0,color='black',lw=.5);axs[1].set(xlabel='Redshift z',ylabel='Piecewise q(z)',
        title='Histories fitted to BAO alone',xlim=(0,1.5));axs[1].legend(fontsize=8)
    fig.suptitle('Conditional required corrections; their physical amplitude is not measured here')
    fig.savefig(out/'required-corrections.png',dpi=170);fig.savefig(out/'required-corrections.pdf');plt.close(fig)
    report={'scope':'E05 screening only; not selection-normalized inference or an identified physical correction',
            'SN_rows':len(sn.z),'unique_CID':int(sn.df.CID.nunique()),'z_cut':'zHD>0.01',
            'SN_reference_LCDM_Omega_m':float(lcdm.x),'SN_reference_chi2':float(lcdm.fun),
            'magnitude_convention':'m=M+mu+deltaM; deltaM_required=mu_SNreference-mu_BAO_history minus GLS constant',
            'template_convention':'Saved positive fixed-clock .030 mag/Gyr templates; multiplying by amplitude is a stipulated law, not remeasured host/progenitor ages.',
            'uncertainty_boundary':'Histories are fixed at BAO optima and their uncertainty is not propagated. Amplitude SD conditions on released total covariance and each fixed template. No independent age prior or cross-probe covariance is added.',
            'selection_boundary':'Uses already standardized/released distances. Cannot establish the proper raw-flux correction or absence of double-counted physical channels.',
            'calculation_validation':'Direct covariance/profile solve; original exact row order; no likelihood interpolation for BAO histories',
            'histories_do_not_constrain_q0_model_independently':True}
    (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    inputs=sn.inputs+[hp]+tp
    outputs=[out/f for f in ['summary.json','template_screen.csv','required_vectors.csv.gz','required-corrections.png','required-corrections.pdf']]
    (out/'manifest.json').write_text(json.dumps({'utc':datetime.now(timezone.utc).isoformat(),
       'source_sha256':hashlib.sha256(SOURCE).hexdigest(),'core_sha256':hashlib.sha256(CORE).hexdigest(),
       'geometry_sha256':hashlib.sha256(GEOMETRY).hexdigest(),
       'inputs_sha256':{str(p.relative_to(ROOT)):sha(p) for p in inputs},
       'outputs_sha256':{str(p.relative_to(ROOT)):sha(p) for p in outputs}},indent=2)+'\n')
    print(pd.DataFrame(fitrows).to_string(index=False))


if __name__=='__main__':main()
