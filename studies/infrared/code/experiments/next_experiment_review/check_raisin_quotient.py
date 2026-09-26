"""Independent NIR28 raw-file and free-intercept covariance quotient audit."""
from pathlib import Path
from functools import lru_cache
from decimal import Decimal
import json,hashlib
import numpy as np
import pandas as pd
from scipy.integrate import quad

P=Path(__file__).resolve().parent;R=P.parents[3];O=P/'quotient-check';O.mkdir(exist_ok=True)
Q=R/'runs/research_2026_09_26/raisin_covariance_quotient'
protocol=json.loads((Q/'protocol.json').read_text());sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
for name,h in protocol['input_sha256'].items():assert sha(R/name)==h,name
result=json.loads((Q/'result.json').read_text());saved=np.load(Q/'arrays.npz')
assert sha(Q/'protocol.json')==result['protocol_sha256'] and sha(Q/'arrays.npz')==result['arrays_sha256']
release=R/'sources/repos/djones1040__RAISIN_DataRelease@a383c4b/distances/w'
members=pd.read_csv(R/'runs/research_2026_09_26/raisin_differential/frozen-membership.csv',dtype={'CID':str})
ids=list(members.CID);n=len(ids);assert n==79
def fit(k):
    p=release/f'nir_dist/RAISIN_combined_FITOPT{k:03d}.FITRES'
    lines=p.read_text().splitlines();names=next(x.split()[1:] for x in lines if x.startswith('VARNAMES:'))
    rows=[dict(zip(names,x.split()[1:])) for x in lines if x.startswith('SN:')]
    assert [v['CID'] for v in rows]==ids
    return {c:np.array([float(v[c]) for v in rows]) for c in ['DLMAG','DLMAGERR','zHD']}
base=fit(0)
# QR of differences with the last object is independent of the saved Helmert basis.
diff=np.c_[np.eye(n-1),-np.ones(n-1)];B=np.linalg.qr(diff.T,mode='reduced')[0].T
assert np.max(abs(B@np.ones(n)))<1e-14 and np.max(abs(B@B.T-np.eye(n-1)))<1e-14
H=B.T@B
@lru_cache(None)
def mu(z):
    integral=quad(lambda x:1/np.sqrt(.3*(1+x)**3+.7),0,float(z),epsabs=1e-13,epsrel=1e-13)[0]
    return 5*np.log10(299792.458/70*(1+z)*integral)+25
mu0=np.array([mu(float(z)) for z in base['zHD']]);vectors=[];fullvectors=[];zchanged=[]
for k in range(1,29):
    t=fit(k);v=t['DLMAG']-base['DLMAG']-(np.array([mu(float(z)) for z in t['zHD']])-mu0)
    vectors.append(B@v)
    off=np.average(t['DLMAG']-base['DLMAG'],weights=1/t['DLMAGERR']**2)
    fullvectors.append(v-off)
    zchanged.append(int(np.count_nonzero(t['zHD']!=base['zHD'])))
V=np.column_stack(vectors);source=V@V.T
FV=np.column_stack(fullvectors);rawsource=FV@FV.T
assert np.max(abs(source-B@rawsource@B.T))<1e-12
sys=release/'nir_syst';tokens=(sys/'RAISIN_all.covmat').read_text().split();assert int(tokens[0])==n
off=np.array(list(map(float,tokens[1:]))).reshape(n,n)
assert np.array_equal(off,off.T) and np.count_nonzero(np.diag(off))==0
allerr=np.loadtxt(sys/'RAISIN_all_lcparams_cosmosis.txt')[:,5];staterr=np.loadtxt(sys/'RAISIN_stat_lcparams_cosmosis.txt')[:,5]
total=off+np.diag(allerr**2);published=off+np.diag(allerr**2-staterr**2)
delta=source-B@published@B.T
# Literal precision of scientific-notation exports; zeros are exactly written zeros.
T=np.array([float(Decimal(10)**Decimal(Decimal(s).as_tuple().exponent)/2) if float(s)!=0 else 0. for s in tokens[1:]]).reshape(n,n)
np.fill_diagonal(T,(allerr+staterr)*1e-6)
T+=1e-12;bound=abs(B)@T@abs(B.T)+1e-12
rootB=saved['B'];rootdelta=saved['nir_28_contrast_difference']
recon_error=float(np.max(abs(B.T@delta@B-rootB.T@rootdelta@rootB)))
assert recon_error<1e-11
invc=np.linalg.inv(total);one=np.ones(n);t=invc@one
profile=invc-np.outer(t,t)/(one@t)
quotient=B.T@np.linalg.inv(B@total@B.T)@B
precision_error=float(np.max(abs(profile-quotient)));assert precision_error<1e-9
ld=np.linalg.slogdet(total)[1];lq=np.linalg.slogdet(B@total@B.T)[1]
det_error=float(abs(lq-ld-np.log(one@t)+np.log(n)));assert det_error<1e-10
u=np.linspace(-.03,.04,n);gray=np.outer(one,u)+np.outer(u,one)
gray_error=float(np.max(abs(B@gray@B.T)));assert gray_error<1e-14
relative=float(np.linalg.norm(delta)/np.linalg.norm(B@published@B.T))
assert abs(relative-result['branches']['nir']['options']['28']['relative_centered_systematic_frobenius_error'])<1e-10
# Validate the parent's rounding expression in its own basis, separately from ours.
rootTol=abs(rootB)@T@abs(rootB.T)+1e-12
roundratio=float(np.max(abs(rootdelta)/rootTol))
rootratio=result['branches']['nir']['options']['28']['max_error_over_propagated_rounding_bound']
assert abs(roundratio/rootratio-1)<1e-7
paths=[Path(__file__),Q/'protocol.json',Q/'result.json',Q/'arrays.npz',R/'scripts/research_2026_09_26/raisin_covariance_quotient.py']
out=dict(status='PASS independent NIR FITOPT1..28 reconstruction and quotient algebra; published/source nonclosure persists',verified_frozen_input_hashes=len(protocol['input_sha256']),basis='QR of differences to last object; comparison with saved Helmert via full centered representation',cosmology='Independent adaptive quadrature of frozen flat Omega_m=.3,H0=70 convention; exact historical imported cosmo module not separately recovered in this check',changed_redshift_counts_by_option=zchanged,source_vs_saved_centered_matrix_max_error=recon_error,profile_precision_identity_max_error=precision_error,marginal_intercept_determinant_identity_error=det_error,common_offset_covariance_projection_max=gray_error,relative_centered_systematic_error_NIR28=relative,parent_rounding_ratio_reproduced=roundratio,our_basis_any_component_exceeds_rounding=bool(np.any(abs(delta)>bound)),rounding_scope='Literal printed covariance/dmb export rounding of the frozen recipe. FITRES DLMAG/DLMAGERR printed at near full precision; this does not bound missing historical upstream products, which are a provenance gap rather than print rounding.',likelihood_scope='B C B^T gives the Gaussian contrast likelihood and, up to fixed sqrt(n), flat-intercept-marginal likelihood. Equality of quadratic precision alone does not justify dropping determinant or prior factors in comparisons with changing C.',scientific_scope='Does not prove the published covariance is wrong; no fitted weights, repaired all matrix or new cosmology.',sha256={str(p.relative_to(R)):sha(p) for p in paths})
(O/'result.json').write_text(json.dumps(out,indent=2)+'\n');np.savez_compressed(O/'arrays.npz',B=B,source_contrast=source,published_contrast=B@published@B.T,difference=delta,rounding_bound=bound)
print(json.dumps({k:v for k,v in out.items() if k not in ['sha256','changed_redshift_counts_by_option']},indent=2))
