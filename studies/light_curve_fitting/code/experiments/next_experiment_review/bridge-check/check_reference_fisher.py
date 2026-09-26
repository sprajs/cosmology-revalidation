"""Independent raw-reference sign and saved local-Fisher arithmetic; no flux outcomes."""
from pathlib import Path
import csv, hashlib, json
import numpy as np
from astropy.io import fits
from scipy.linalg import null_space

ROOT=Path(__file__).resolve().parents[5]
OUT=Path(__file__).resolve().parent
SRC=ROOT/'runs/research_2026_09_26/bayesn_distance_identification'
KCOR=ROOT/'sources/repos/djones1040__RAISIN_DataRelease@a383c4b/kcor/kcor_CSPDR3_BD17.fits'
NATIVE=SRC/'official-code/bayesn/bayesn-filters'
FILTER=NATIVE/'filters/LCO/Swope/B_tel_ccd_atm_ext_1.2.dat'
STAR=NATIVE/'standards/bd_17d4708_stisnic_007.dat'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def integral(curve,star):
    """Exact integral of lambda times piecewise-linear T,F (two Gauss nodes/interval)."""
    lo=max(curve[0,0],star[0,0]);hi=min(curve[-1,0],star[-1,0])
    grid=np.unique(np.r_[lo,hi,curve[:,0],star[:,0]])
    grid=grid[(grid>=lo)&(grid<=hi)]
    mid=(grid[1:]+grid[:-1])/2;half=np.diff(grid)/2
    x=np.column_stack([mid-half/np.sqrt(3),mid+half/np.sqrt(3)])
    t=np.interp(x,curve[:,0],curve[:,1]);f=np.interp(x,star[:,0],star[:,1])
    return float(np.sum(half[:,None]*x*t*f)),float(np.sum(half[:,None]*x*t))

def zero(curve,star,mag):
    a,b=integral(curve,star);return a/b*10**(.4*mag)

with fits.open(KCOR) as h:
    c=np.column_stack([h['FilterTrans'].data.field(0),h['FilterTrans'].data['CSP-B']]).astype(float)
    s=np.column_stack([h['PrimarySED'].data.field(0),h['PrimarySED'].data['BD17']]).astype(float)
    z=h['ZPoff'].data
    row=z[np.array([str(x).strip()=='CSP-B' for x in z['Filter Name']])][0]
    mag=float(row['Primary Mag'])
cn=np.loadtxt(FILTER);sn=np.loadtxt(STAR);mn=9.896
old=zero(c,s,mag);native_same=zero(c,sn,mn);native_own=zero(cn,sn,mn)
# The released primary itself has magnitude m_old. Through the same curve,
# m_native=-2.5log10(count_star/count_native_ref)+m_native_ref.
count_release,_=integral(c,s);count_native,_=integral(c,sn)
star_m_native=-2.5*np.log10(count_release/count_native)+mn
delta=star_m_native-mag
assert abs(delta-2.5*np.log10(native_same/old))<1e-12
assert delta>0
bridge=next(r for r in csv.DictReader((SRC/'calibration-bridge.csv').open()) if r['survey']=='CSP' and r['raw_band']=='B')
reference=dict(release_primary_mag=mag,native_primary_mag=mn,release_primary_native_magnitude_same_curve=float(star_m_native),
    same_curve_delta_mag=float(delta),own_curve_zero_ratio_delta_mag=float(2.5*np.log10(native_own/old)),
    same_curve_scale_native_over_release=float(old/native_same),
    scale_from_magnitude_delta=float(10**(-.4*delta)),
    agent_simpson_same_curve_delta_mag=float(bridge['same_curve_delta_mag']),
    independent_linear_integral_minus_agent_cubic_simpson_mag=float(delta-float(bridge['same_curve_delta_mag'])),
    native_own_vs_same_curve_delta_mag=float(2.5*np.log10(native_own/native_same)),
    interpretation='Positive native-minus-release magnitude means smaller native FLUXCAL. Same-curve sign directly verified using the released standard as the test star. Independent piecewise-linear merged-grid integration; difference from source cubic/Simpson is numerical/interpolation sensitivity, not measured calibration bias.')

a=np.load(SRC/'forward-geometry-arrays.npz',allow_pickle=False)
decl=json.loads((SRC/'forward-geometry-result.json').read_text())
results=[]
for i,r in enumerate(decl['objects']):
    keep=a['mask'][:,i]>0
    A=a['Jhalf'][keep,i,:]/a['errors'][keep,i,None]
    dust=A[:,1:3];g=A[:,0:1];E=A[:,3:]
    # Independent data-space marginal covariance, not parameter-space Schur subtraction.
    K=np.eye(len(A))+E@E.T
    Wd=np.linalg.solve(K,dust);Wg=np.linalg.solve(K,g)
    no=dust.T@Wd-(dust.T@Wg)@np.linalg.solve(g.T@Wg,g.T@Wd)
    Kext=K+(r['sigma_ext_mag']**2+.088**2)*(g@g.T)
    ext=dust.T@np.linalg.solve(Kext,dust)
    N=np.column_stack([g,E]);Q=null_space(N.T)
    free=(Q.T@dust).T@(Q.T@dust)
    sing=np.linalg.svd(N,compute_uv=False)
    prior_info=np.asarray(r['no_distance']['information'])
    invsqrt=np.linalg.inv(np.linalg.cholesky((no+no.T)/2))
    relative_free=invsqrt@free@invsqrt.T
    item=dict(CID=r['CID'],epochs=len(A),nuisance_rank=len(A)-Q.shape[1],singular_values=sing.tolist(),
       no_distance_information=no.tolist(),external_distance_information=ext.tolist(),unconstrained_information=free.tolist(),
       no_distance_max_error=float(np.max(abs(no-prior_info))),
       external_distance_max_error=float(np.max(abs(ext-np.asarray(r['external_distance']['information'])))),
       unconstrained_max_error=float(np.max(abs(free-np.asarray(r['unregularized_dust_information'])))),
       generalized_free_vs_trained_information_eigenvalues=np.linalg.eigvalsh(relative_free).tolist())
    results.append(item)
assert max(r['no_distance_max_error'] for r in results)<1e-7
assert max(r['external_distance_max_error'] for r in results)<1e-7
files=[KCOR,FILTER,STAR,SRC/'calibration-bridge.csv',SRC/'forward-geometry-arrays.npz',SRC/'forward-geometry-result.json',SRC/'forward_geometry.py',Path(__file__)]
out=dict(status='Arithmetic and reference sign PASS; physical forward resolution gate currently fails, so information numbers remain provisional.',
    reference=reference,local_fisher=results,parent_forward_gates_pass=decl['gates_pass'],
    interpretation='The Fisher calculation is local joint likelihood plus prior curvature, with the trained epsilon and theta priors held fixed. Unconstrained epsilon is an identifiability stress, not a proposed physical population. No measured SN flux enters this check.',
    inputs={str(p.relative_to(ROOT)):sha(p) for p in files})
(OUT/'reference-fisher-check.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(dict(reference=reference,max_no_distance_error=max(r['no_distance_max_error'] for r in results),max_external_error=max(r['external_distance_max_error'] for r in results),free_fraction_eigenvalues={r['CID']:r['generalized_free_vs_trained_information_eigenvalues'] for r in results},forward_gates_pass=decl['gates_pass']),indent=2))
