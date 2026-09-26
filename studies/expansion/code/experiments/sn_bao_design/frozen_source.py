#!/usr/bin/env python3
"""Outcome-free design for a relative SN-versus-transverse-BAO comparison.

The SN magnitude column is deliberately not read. This is a compressed-data
feasibility and conditional power calculation, not a cosmic-opacity result.
"""
import argparse
import hashlib
import json
from pathlib import Path
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy.linalg import cho_factor, cho_solve
from scipy.stats import norm
from numpy.polynomial.legendre import leggauss

ROOT = Path(__file__).resolve().parents[2]
SN = ROOT/'sources/repos/CobayaSampler__sn_data/PantheonPlus'
BAO = ROOT/'sources/repos/CobayaSampler__bao_data/desi_bao_dr2'
OUT = ROOT/'runs/research_2026_09_26/sn_bao_design'
SCRIPT = Path(__file__).resolve()
INPUTS = [SN/'Pantheon+SH0ES.dat', SN/'Pantheon+SH0ES_STAT+SYS.cov',
          BAO/'desi_gaussian_bao_ALL_GCcomb_mean.txt',
          BAO/'desi_gaussian_bao_ALL_GCcomb_cov.txt']


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def write(p, x):
    p.write_text(json.dumps(x, indent=2, allow_nan=False)+'\n')


def freeze():
    OUT.mkdir(parents=True,exist_ok=True)
    p = OUT/'protocol.json'
    if p.exists():
        raise FileExistsError(p)
    write(p,dict(created_utc=datetime.now(timezone.utc).isoformat(),
        stage='Before design results, no SN apparent-magnitude column read by this workstream. These releases were previously analysed elsewhere in the project.',
        estimand='Relative total standardized-luminosity/opacity/calibration/selection drift: r_SN(z)-5log10[(DM/rd)/z], with one free common magnitude/ruler intercept. r_SN=m_b_corr-5log10[(1+zHEL)zHD].',
        assumptions='Published standardized magnitudes and total covariance retained. Constant comoving ruler, standard distance duality, released Gaussian effective-redshift BAO approximation, unknown cross-probe covariance taken zero for forecast only. No LCDM fit, H0 or rd calibration.',
        selection='SN zHD>.01 as existing released-product branch; preserve duplicate observation rows and full covariance. BAO DM nodes inside SN redshift support; exclude DV and DH, use marginal DM covariance submatrix.',
        interpolation='Local GLS polynomial in [ln(1+z)-ln(1+z_node)]/halfwidth, primary degree2 halfwidth.10. Frozen sensitivities degree2 widths.075,.125 and degree3 width.10. Marginal local C in each GLS; cross-node errors from W Cfull W^T. No kernel tuning with outcomes.',
        admission='For every frozen interpolation design: >=20 distinct CID, >=5 distinct CID strictly each side, full polynomial rank, condition number<1e8, max named-test-curve interpolation bias<=.005mag. Interpolation tests are a finite sensitivity set, not a bound over arbitrary smooth histories.',
        recovery_curves='Flat constant-w Om in[.15,.3,.5], w in[-1.5,-1,-.5]; positive kinematic E=(1+z)^(1+q0+q1)exp[-q1*z/(1+z)], q0 in[-1,0,.5], q1 in[0,1,2]. Test r=5log10[(integral dz/E)/z]. Curves test approximation only, not priors or fits.',
        power='If >=2 nodes admitted, estimate Gaussian forecast SE for drift linear in ln(1+z) normalized to1 between min/max admitted nodes. Profile intercept. Power at two-sided5% for endpoint differences .02,.05,.10mag. Retain inter-node SN covariance. BAO log delta-method approximation verified with50000 Gaussian draws, seed260926117.',
        stop='Do not score observed drift in this script. Broad BAO redshift kernels are unresolved; interpolation recovery alone does not validate compressed likelihood under arbitrary histories. No opacity attribution or cosmology correction.',
        input_sha256={str(p.relative_to(ROOT)):sha(p) for p in INPUTS+[SCRIPT]}))
    print(json.dumps({'protocol_sha256':sha(p)}))


def curve(z, family, p1, p2, nodes=96):
    gx,gw=leggauss(nodes)
    zz=z[:,None]*(gx+1)/2
    if family=='w':
        e=np.sqrt(p1*(1+zz)**3+(1-p1)*(1+zz)**(3*(1+p2)))
    else:
        e=(1+zz)**(1+p1+p2)*np.exp(-p2*zz/(1+zz))
    d=z/2*np.sum(gw/e,axis=1)
    return 5*np.log10(d/z)


def design():
    protocol=json.loads((OUT/'protocol.json').read_text())
    for name,h in protocol['input_sha256'].items():
        assert sha(ROOT/name)==h,name
    if (OUT/'result.json').exists():
        raise FileExistsError('Preserve original design')
    df=pd.read_csv(INPUTS[0],sep=r'\s+',usecols=['CID','zHD','zHEL'])
    raw=np.loadtxt(INPUTS[1]); n=int(raw[0]); assert n==len(df)
    full=raw[1:].reshape(n,n)
    asym=float(abs(full-full.T).max());assert asym<1e-7
    full=(full+full.T)/2
    take=np.flatnonzero(df.zHD.to_numpy()>.01)
    df=df.iloc[take].reset_index(drop=True); C=full[np.ix_(take,take)]
    z=df.zHD.to_numpy(); t=np.log1p(z)
    bd=pd.read_csv(INPUTS[2],sep=r'\s+',comment='#',names=['z','value','kind'])
    bc=np.loadtxt(INPUTS[3])
    bix=np.flatnonzero((bd.kind=='DM_over_rs')&(bd.z>=z.min())&(bd.z<=z.max()))
    zb=bd.z.to_numpy()[bix]; d=bd.value.to_numpy()[bix]
    cb=bc[np.ix_(bix,bix)]
    cbmag=cb*np.outer(5/np.log(10)/d,5/np.log(10)/d)
    families=[('w',om,w)for om in [.15,.3,.5]for w in[-1.5,-1.,-.5]]
    families += [('q',q0,q1)for q0 in[-1.,0.,.5]for q1 in[0.,1.,2.]]
    zz=np.r_[z,zb]
    curves=np.array([curve(zz,*f)for f in families])
    quaderr=float(max(np.max(abs(curve(zz,*f,192)-curves[i]))for i,f in enumerate(families)))
    assert quaderr<1e-10
    settings=[('primary',2,.10),('narrow',2,.075),('wide',2,.125),('cubic',3,.10)]
    rows=[]; matrices={}; admissible=[]
    for label,degree,hw in settings:
        W=np.zeros((len(zb),len(z))); gate=[]
        for j,znode in enumerate(zb):
            ix=np.flatnonzero(abs(t-np.log1p(znode))<=hw)
            x=(t[ix]-np.log1p(znode))/hw
            X=np.stack([x**k for k in range(degree+1)],axis=1)
            Cp=C[np.ix_(ix,ix)]
            cx=cho_solve(cho_factor(Cp),X)
            gram=X.T@cx; condition=float(np.linalg.cond(gram))
            assert np.linalg.matrix_rank(gram)==degree+1
            w=np.linalg.solve(gram,cx.T)[0]
            W[j,ix]=w
            bias=curves[:,:len(z)]@W[j]-curves[:,len(z)+j]
            counts=[len(set(df.CID.iloc[ix])),len(set(df.CID.iloc[ix[x<0]])),len(set(df.CID.iloc[ix[x>0]]))]
            passed=counts[0]>=20 and min(counts[1:])>=5 and condition<1e8 and np.max(abs(bias))<=.005
            gate.append(bool(passed))
            assert np.max(abs(w@X-np.r_[1.,np.zeros(degree)]))<1e-10
            rows.append(dict(design=label,degree=degree,halfwidth_log1pz=hw,z_node=float(znode),
                observation_rows=len(ix),distinct_CID=counts[0],distinct_below=counts[1],distinct_above=counts[2],
                min_z=float(z[ix].min()),max_z=float(z[ix].max()),condition=condition,
                sum_abs_weights=float(abs(w).sum()),largest_abs_weight=float(abs(w).max()),
                sn_sigma_mag=float(np.sqrt(w@Cp@w)),bao_sigma_mag=float(np.sqrt(cbmag[j,j])),
                max_finite_family_bias_mag=float(abs(bias).max()),admitted=bool(passed)))
        matrices[label]=W;admissible.append(gate)
    admitted=np.all(admissible,axis=0);ix=np.flatnonzero(admitted)
    forecast={}
    for label,W in matrices.items():
        csn=W@C@W.T; V=csn+cbmag
        if len(ix)>=2:
            g=(np.log1p(zb[ix])-np.log1p(zb[ix[0]]))/(np.log1p(zb[ix[-1]])-np.log1p(zb[ix[0]]))
            X=np.column_stack([np.ones(len(ix)),g]);vi=V[np.ix_(ix,ix)]
            pcov=np.linalg.inv(X.T@cho_solve(cho_factor(vi),X)); se=np.sqrt(pcov[1,1])
            power={str(a):float(norm.cdf(-norm.ppf(.975)-a/se)+norm.sf(norm.ppf(.975)-a/se))for a in[.02,.05,.10]}
            forecast[label]=dict(endpoint_drift_se_mag=float(se),two_sided_5pct_power=power,
                point_sigma_mag=np.sqrt(np.diag(V)).tolist())
        matrices[label+'_sn_cov']=csn
        matrices[label+'_total_cov']=V
    rng=np.random.default_rng(260926117)
    draws=rng.multivariate_normal(d,cb,size=50000)
    assert np.min(draws)>0
    logdelta=5*np.log10(draws/d)
    mc=np.cov(logdelta,rowvar=False)
    # MC errors are recorded rather than interpreted as physical model accuracy.
    delta_check=dict(draws=50000,seed=260926117,
        mean_bias_mag=logdelta.mean(axis=0).tolist(),analytic_second_order_bias_mag=(-2.5/np.log(10)*np.diag(cb)/d**2).tolist(),
        covariance_relative_frobenius_error=float(np.linalg.norm(mc-cbmag)/np.linalg.norm(cbmag)))
    pd.DataFrame(rows).to_csv(OUT/'design-table.csv',index=False)
    np.savez_compressed(OUT/'weights.npz',**matrices,admitted=admitted,bao_z=zb,bao_magnitude_cov=cbmag,
        original_sn_rows=take,sn_z=z)
    write(OUT/'result.json',dict(protocol_sha256=sha(OUT/'protocol.json'),sn_rows=len(z),sn_distinct_CID=df.CID.nunique(),
        bao_nodes=zb.tolist(),admitted_nodes=zb[ix].tolist(),input_covariance_asymmetry_mag2=asym,
        quadrature_check_max_mag=quaderr,forecast=forecast,bao_log_delta_method_check=delta_check,
        no_observed_SN_magnitudes_read=True,no_outcome_scoring=True,
        remaining_gate='Exact BAO compression/redshift-kernel response or an explicitly limited smooth-history validity domain. Finite toy-family interpolation recovery is not proof for arbitrary histories.',
        artifacts_sha256={p.name:sha(p)for p in[OUT/'design-table.csv',OUT/'weights.npz']}))
    (OUT/'executed_source.py').write_bytes(SCRIPT.read_bytes())
    print((OUT/'result.json').read_text())


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['freeze','design'])
    a=p.parse_args();globals()[a.mode]()
