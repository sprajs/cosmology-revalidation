"""Independent second-start check of the saved six-object pattern refit."""
from datetime import datetime,timezone
from pathlib import Path
import hashlib,json,sys
import numpy as np,pandas as pd
from scipy.linalg import solve_triangular
from scipy.optimize import least_squares
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from scripts.salt_dust_audit.flux_response import build_model,K
OUT=ROOT/'runs/research_2026_09_26/calibration_pattern_verification'
ORIG=ROOT/'runs/research_2026_09_26/calibration_pattern_refit'
COEF=ROOT/'runs/research_2026_09_26/astra_design/validation1020/frozen-discovery-coefficients.npz'
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
manifest=json.loads((ORIG/'manifest.json').read_text())
assert sha(ROOT/'scripts/research_2026_09_26/calibration_pattern_refit.py')==manifest['source_sha256']
assert sha(ORIG/'executed_source.py')==manifest['source_sha256']
for rel,digest in {**manifest['inputs_sha256'],**manifest['outputs_sha256']}.items():assert sha(ROOT/rel)==digest,rel
coeff=np.load(COEF);residual_dimming=coeff['gauge_griz']@coeff['basis_mean'];dmag=-residual_dimming
summary=json.loads((ORIG/'summary.json').read_text())
np.testing.assert_allclose(dmag,summary['delta_observed_magnitude_griz'],atol=1e-15,rtol=0)
assert abs(dmag.sum())<1e-14
original=pd.read_csv(ORIG/'paired_refits.csv',dtype={'CID':str})
cohort=pd.read_csv(ORIG/'cohort.csv',dtype={'CID':str})
assert len(original)==len(cohort)==6 and cohort.CID.is_unique
model,bands,paths,zp=build_model();offsets={str(r['Filter Name'])[-1]:float(r['Primary Mag']) for r in zp}
D=np.array([1,.16087,-3.1178,0.]);second_start=np.array([.05,.25,.03,2.0]);records=[]
for row in original.itertuples():
    p=ROOT/f'runs/research_2026_09_26/astra_design/exact43/objectives/objective_{row.CID}.npz';o=np.load(p)
    t=o['MJD'];band=o['band'];bs=np.array([bands[b] for b in band],dtype=object)
    fc=np.array([10**(-.4*(.27+offsets[b])) for b in band]);correction=np.array([dmag['griz'.index(b)] for b in band]);scale=np.exp(-K*correction)
    x0,x1,c,t0=o['parameters_x0_x1_c_t0'];z=float(o['zHEL'][0]);ebv=float(o['MWEBV'][0]);y=o['data_flux'];C=o['frozen_flux_covariance'];L=np.linalg.cholesky(C)
    transformed_C=scale[:,None]*C*scale[None,:];transformed_L=np.linalg.cholesky(transformed_C)
    chol_error=float(np.max(abs(transformed_L-scale[:,None]*L)))
    assert np.allclose(transformed_L,scale[:,None]*L,atol=1e-8,rtol=1e-10)
    def flux(theta):
        model.set(z=z,t0=t0+theta[3],x0=x0*np.exp(-K*theta[0]),x1=x1+theta[1],c=c+theta[2],mwebv=ebv,mwrv=3.1,hostebv=0,hostrv=3.1)
        return fc*model.bandflux(bs,t,zp=27.5,zpsys='ab')
    probe=flux(second_start)
    lhs=solve_triangular(transformed_L,probe-scale*y,lower=True)
    rhs=solve_triangular(L,probe/scale-y,lower=True)
    whiten_error=float(np.max(abs(lhs-rhs)));assert whiten_error<1e-8
    fun=lambda theta,s:solve_triangular(L,flux(theta)/s-y,lower=True)
    fit=lambda s:least_squares(lambda theta:fun(theta,s),second_start,xtol=1e-12,ftol=1e-12,gtol=1e-10,max_nfev=400,bounds=([-2,-8,-1,-30],[2,8,1,30]))
    old=fit(np.ones(len(y)));new=fit(scale)
    shift=float(D@(new.x-old.x));old_chi2=float(old.fun@old.fun);new_chi2=float(new.fun@new.fun)
    records.append({'CID':row.CID,'zHEL':z,'second_start':second_start.tolist(),'baseline_success':bool(old.success),'alternative_success':bool(new.success),
       'baseline_chi2_second':old_chi2,'alternative_chi2_second':new_chi2,'baseline_chi2_original':row.baseline_chi2,'alternative_chi2_original':row.alternative_chi2,
       'delta_standardized_second':shift,'delta_standardized_original':row.nonlinear_delta_standardized,
       'baseline_chi2_difference':old_chi2-row.baseline_chi2,'alternative_chi2_difference':new_chi2-row.alternative_chi2,'delta_standardized_difference':shift-row.nonlinear_delta_standardized,
       'baseline_optimality_second':float(old.optimality),'alternative_optimality_second':float(new.optimality),'cholesky_max_abs_error':chol_error,'whitened_residual_max_abs_error':whiten_error})
    print(row.CID,'shift difference',records[-1]['delta_standardized_difference'],flush=True)
frame=pd.DataFrame(records)
assert frame.baseline_success.all() and frame.alternative_success.all()
max_shift=float(frame.delta_standardized_difference.abs().max());max_chi=float(max(frame.baseline_chi2_difference.abs().max(),frame.alternative_chi2_difference.abs().max()))
assert max_shift<1e-4 and max_chi<1e-4,(max_shift,max_chi)
high=frame.loc[frame.zHEL.idxmax()]
assert abs(high.delta_standardized_second-.17887031456855215)<1e-4
OUT.mkdir(exist_ok=False)
frame.to_csv(OUT/'second_start_refits.csv',index=False)
report={'scope':'Post-run numerical verification of a hypothetical fixed-vector response; not a measured calibration bias.',
 'data_magnitude_correction_griz':dmag.tolist(),'frozen_residual_dimming_griz':residual_dimming.tolist(),'correction_is_negative_frozen_pattern':bool(np.allclose(dmag,-residual_dimming,atol=1e-15,rtol=0)),
 'covariance_rule':'C_new=S C S, chol(C_new)=S chol(C) for positive diagonal S; whitening L^-1(f/S-y) is exact',
 'local_sign':'For S=exp(-K dmag), f/S=f+K f dmag+O(dmag^2); local refit shift solves J delta=-K f dmag',
 'second_start':second_start.tolist(),'fit_pairs':len(frame),'all_success':True,'max_abs_shift_difference_mag':max_shift,'max_abs_chi2_difference':max_chi,
 'max_abs_cholesky_error':float(frame.cholesky_max_abs_error.max()),'max_abs_whitened_residual_error':float(frame.whitened_residual_max_abs_error.max()),
 'high_z_object':{'CID':high.CID,'zHEL':float(high.zHEL),'delta_standardized_mag':float(high.delta_standardized_second),'chi2_change':float(high.alternative_chi2_second-high.baseline_chi2_second)},
 'interpretation':'Large high-z fitted-distance change with negligible residual chi2 change is weak identification under this conditional frozen pattern; not evidence of an observed bias.'}
(OUT/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
(OUT/'manifest.json').write_text(json.dumps({'verified_utc':datetime.now(timezone.utc).isoformat(),'verifier_source_sha256':sha(Path(__file__)),'original_manifest_sha256':sha(ORIG/'manifest.json'),'frozen_coefficients_sha256':sha(COEF),'outputs_sha256':{str(p.relative_to(ROOT)):sha(p) for p in OUT.iterdir() if p.is_file() and p.name!='manifest.json'}},indent=2)+'\n')
print(json.dumps({k:report[k] for k in ('fit_pairs','max_abs_shift_difference_mag','max_abs_chi2_difference','max_abs_cholesky_error','max_abs_whitened_residual_error','high_z_object')},indent=2))
