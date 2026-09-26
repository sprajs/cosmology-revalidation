"""Conditional scale/curvature sensitivity, not an inferred systematic correction."""
from pathlib import Path
import hashlib
import json
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.stats import norm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from bao_shape import ROOT, OUT, load_data

SOURCE=Path(__file__).read_bytes()


def main():
    out=OUT/'assumption_budget'
    if out.exists():raise FileExistsError(out)
    out.mkdir();(out/'executed_source.py').write_bytes(SOURCE)
    z,y,C,hashes,_,raw,rawcov=load_data();n=len(z)
    idx=int(np.flatnonzero(raw.kind=='DV_over_rs')[0]);b=float(raw.z.iloc[idx]);tb=np.log1p(b)
    fac=(b/(1+b)*tb**2)**(1/3)
    s=float(raw.value.iloc[idx]/fac);ss=float(np.sqrt(rawcov[idx,idx])/fac)
    gj=y[n:];gs=np.sqrt(np.diag(C)[n:]);ratio=s/gj
    allrows=[]
    for j in range(n):
        # All possible correlations preserve Var(s-g)<= (sd_s+sd_g)^2.
        deficit=(gj[j]-s)/(ss+gs[j])
        curvature=0.
        if ratio[j]<1:
            x=brentq(lambda x:np.sinc(x/np.pi)-ratio[j]**1.5,1e-10,np.pi-1e-8)
            curvature=(x/tb)**2
        allrows.append({'z':z[j],'s_over_g':ratio[j],
            'BGS_fractional_increase_to_equality':gj[j]/s-1,
            'DH_fractional_decrease_to_equality':1-ratio[j],
            'necessary_abs_closed_Omega_k_central_values':curvature,
            'worst_correlation_standardized_deficit':deficit,
            'six_test_Bonferroni_upper_p_any_correlation':min(1.,6*norm.sf(max(0,deficit)))})
    pd.DataFrame(allrows).to_csv(out/'per_redshift_budget.csv',index=False)
    best=allrows[int(np.argmax(gj-s))]
    kappas=np.array([0,.01,.1,1,3,5,8.291275330017672,10])
    table=[]
    for kappa in kappas:
        u=np.sqrt(kappa)*tb
        f=np.sinc(u/np.pi)**(2/3)
        zz=(f*gj-s)/np.sqrt(ss**2+f*f*gs**2)
        # Bonferroni is conservative with arbitrary correlation among the six contrasts.
        bound=min(1.,6*norm.sf(max(0,float(zz.max()))))
        table.append({'Omega_k_lower_bound':-kappa,'BGS_bound_factor':f,
                      'maximum_deficit_standardized':float(zz.max()),
                      'six_test_Bonferroni_upper_p_under_released_crosscov':bound})
    pd.DataFrame(table).to_csv(out/'curvature_sensitivity.csv',index=False)
    report={'status':'conditional sensitivity under noacceleration, not a fitted physical bias or curvature posterior',
            'worst_central_contrast':best,
            'curvature_derivation':'q>=0 implies E>=1+z and I_b<=ln(1+b). For Omega_k=-kappa, DM/DC=sinc(sqrt(kappa)*I_b)>=sinc(sqrt(kappa)*ln(1+b)) while sqrt(kappa)*ln(1+b)<pi. Therefore s_b/g_j>=sinc(...)^(2/3).',
            'geometry_domain':'standard FLRW Omega_k; positive expansion and same constant ruler; fixed point-effective-redshift compressed likelihood',
            'ruler_interpretation':'At coasting equality, rd(z_j)/rd(z_b)=s_b/g_j. This is a required relative ruler/extraction scale, not measured ruler evolution.',
            'correlation_bound':'Uses only the quoted two marginal variances; six-test Bonferroni remains valid for any cross-correlation. It does not cover unmodeled extra variance or mean bias.',
            'curvature_scan':table}
    (out/'results.json').write_text(json.dumps(report,indent=2)+'\n')
    fig,ax=plt.subplots(figsize=(8,4.8),layout='constrained')
    ax.axhspan(s-ss,s+ss,color='tab:orange',alpha=.25,label='BGS upper bound: 1σ measurement band')
    ax.axhline(s,color='tab:orange',lw=1.5)
    ax.errorbar(z,gj,yerr=gs,fmt='o',capsize=3,label=r'Radial BAO: $(1+z)D_H/r_d$')
    ax.set(xlabel='Redshift z',ylabel='Dimensionless BAO quantity',
           title='No acceleration + flat geometry requires every blue mean below orange')
    ax.legend(fontsize=9,loc='lower left')
    fig.savefig(out/'bgs-radial-bound.png',dpi=180);fig.savefig(out/'bgs-radial-bound.pdf');plt.close(fig)
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    outputs=[p for p in out.iterdir() if p.name!='executed_source.py']
    (out/'manifest.json').write_text(json.dumps({'utc':datetime.now(timezone.utc).isoformat(),
        'inputs_sha256':hashes,'source_sha256':hashlib.sha256(SOURCE).hexdigest(),
        'outputs_sha256':{str(p.relative_to(ROOT)):sha(p) for p in outputs}},indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
