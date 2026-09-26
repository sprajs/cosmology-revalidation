"""Fixed stellar-branch vectors projected into SN data; amplitudes never fitted."""
from pathlib import Path
import json
import hashlib
import numpy as np
import pandas as pd
from scipy.linalg import solve_triangular

ROOT=Path(__file__).resolve().parents[2]


def main():
    out=ROOT/'runs/research_2026_09_26/calibration_stars/sn_transfer'
    out.mkdir(exist_ok=False)
    src=Path(__file__).read_bytes();(out/'executed_source.py').write_bytes(src)
    starfile=out.parent/'extinction_branch_changes.csv'
    meta=ROOT/'runs/research_2026_09_26/foreground_response/map_samples.csv'
    matrixfile=ROOT/'runs/salt_dust_audit/flux_response/matrices.npz'
    d=pd.read_csv(meta,dtype={'CID':str});s=pd.read_csv(starfile);s=s[s.arm=='fixed_intersection']
    a=np.load(matrixfile);D=np.array([1,.16087,-3.1178,0]);rows=[]
    for cat,p in s.groupby('library'):
        delta=p.set_index('band').loc[list('griz'),'delta_mean_correction_stored_minus_zero'].to_numpy()
        for row in d.itertuples():
            pre=row.CID+'__';J=a[pre+'jacobian_flux'];G=a[pre+'nuisance_flux'][:,5:]
            residual=a[pre+'flux_observed']-a[pre+'flux_model']
            for weighting in ['measurement_only','measurement_plus_model']:
                L=np.linalg.cholesky(a[pre+weighting+'_covariance'])
                jw=solve_triangular(L,J,lower=True);rw=solve_triangular(L,residual,lower=True)
                # Positive correction adds magnitude to observed photometry.
                v=solve_triangular(L,G@delta,lower=True)
                Q=np.linalg.qr(jw,mode='reduced')[0];vp=v-Q@(Q.T@v);rp=rw-Q@(Q.T@rw)
                shift=float(D@a[pre+weighting+'_response'][:,5:]@delta)
                rows.append(dict(library=cat,CID=row.CID,zHEL=row.zHEL,field=row.field,highIa=row.PROB_SNNV19>.999,
                    weighting=weighting,delta_standardized_preBBC=shift,
                    residual_loglike_gain_if_applied=float(-rp@vp-.5*(vp@vp)),template_information=float(vp@vp)))
    r=pd.DataFrame(rows);r.to_csv(out/'object_responses.csv',index=False)
    summaries=[]
    for (cat,w),p in r.groupby(['library','weighting']):
        for subset in ['all64','highIa43']:
            t=p if subset=='all64' else p[p.highIa]
            v=t.delta_standardized_preBBC.to_numpy();z=t.zHEL.to_numpy();q=np.quantile(z,[.25,.75])
            summaries.append(dict(library=cat,weighting=w,subset=subset,objects=len(t),
                fixed_score_gain=float(t.residual_loglike_gain_if_applied.sum()),template_SNR=float(np.sqrt(t.template_information.sum())),
                shift_mean_mag=float(v.mean()),shift_centered_RMS_mag=float(v.std()),
                high_minus_low_z_quartile_shift_mag=float(v[z>=q[1]].mean()-v[z<=q[0]].mean()),
                shift_min_mag=float(v.min()),shift_max_mag=float(v.max())))
    report={'scope':'Fixed descriptive stellar-extinction branch offsets applied to observed magnitude convention using original trained SALT and fixed native covariance; no fitted amplitude or calibration/selection/BBC closure.',
        'sign':'Positive delta adds to observed magnitude; fitted parameters shift +R_zp delta. Residual log likelihood changes -r_perp dot v_perp minus norm(v_perp)^2/2.',
        'summaries':summaries}
    (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    (out/'manifest.json').write_text(json.dumps({'source_sha256':hashlib.sha256(src).hexdigest(),
        'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in [starfile,meta,matrixfile]}},indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
