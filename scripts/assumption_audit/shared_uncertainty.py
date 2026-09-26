"""Audit shared correction uncertainty without altering the original fits.

The diagonal variant is an intentionally wrong counterfactual for a single
global slope. Importance weighting uses the *marginal* Gaussian likelihood,
integrating that slope exactly, not the joint likelihood at a sampled slope.
"""
from pathlib import Path
import hashlib
import json
import sys

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts/cosmology'))
from core import Pantheon, BAO, mu


def digest(p):
    return hashlib.file_digest(open(p, 'rb'), 'sha256').hexdigest()


def projected(c):
    inv = cho_solve(cho_factor(c, lower=True), np.eye(len(c)))
    u = inv.sum(axis=1)
    return inv - np.outer(u, u) / u.sum()


def main():
    out = ROOT / 'runs/assumption_audit'
    out.mkdir(parents=True, exist_ok=True)
    sn, bao = Pantheon(), BAO()
    template_path = ROOT / 'runs/cosmology/templates/c14-cpl63-median.csv'
    d = pd.read_csv(template_path)
    t = PchipInterpolator(d.z, d.delta_mu)(sn.z)
    s = 4 / 30  # exactly the original shared amplitude prior
    covs = {'fixed': sn.cov, 'shared': sn.cov + s*s*np.outer(t,t),
            'incorrect_independent_rows': sn.cov + np.diag((s*t)**2)}
    projectors = {k: projected(c) for k,c in covs.items()}
    y = sn.mag - sn.fid - t
    comp = {k: (y@a@y, sn.P.T@a@y, sn.P.T@a@sn.P) for k,a in projectors.items()}

    def chi(theta, treatment):
        theta = np.atleast_2d(theta)
        g = mu(sn.nodes, theta[:,:3]) - sn.fidnodes
        yy, by, bb = comp[treatment]
        return yy - 2*g@by + np.sum((g@bb)*g,axis=1) + bao.chisq(theta[:,:3],theta[:,3])

    base = ROOT/'runs/cosmology/pantheon-bao-cpl-c14slope'
    chain_path = base/'chains.npz'
    a = np.load(chain_path)
    draws = a['chain'].reshape(-1,5)
    rng = np.random.default_rng(920612)
    ix = rng.choice(len(draws),min(80000,len(draws)),replace=False)
    draws = draws[ix]
    q = .5 + 1.5*draws[:,1]*(1-draws[:,0])
    reference = np.concatenate([chi(c[:,:4],'shared') for c in np.array_split(draws,160)])
    result = {'interpretation': 'Conditional fixed-template Gaussian sensitivity; not new measured age uncertainty.',
              'normal_slope_mag_per_Gyr': [0.030, 0.004], 'source_draws': len(draws), 'treatments': {}}
    for label in covs:
        c = np.concatenate([chi(v[:,:4],label) for v in np.array_split(draws,160)])
        lw = -.5*(c-reference)
        w = np.exp(lw-lw.max()); w /= w.sum()
        m = np.sum(w*q)
        order = np.argsort(q); cdf=np.cumsum(w[order])-.5*w[order]
        x0 = np.array([.35,-.4,-1.7,10000.])
        scales = np.array([1,1,1,10000.])
        def fun(x):
            v=x*scales
            if not (.01<v[0]<.99 and -3<v[1]<1 and -3<v[2]<2 and v[1]+v[2]<0 and 5000<v[3]<15000):
                return 1e20
            return float(chi(v,label)[0])
        opt = minimize(fun,x0/scales,method='Nelder-Mead',options={'maxiter':8000,'xatol':1e-9,'fatol':1e-8})
        x=opt.x*scales
        assert opt.success
        result['treatments'][label] = {
            'q0_mean':float(m), 'q0_sd':float(np.sqrt(np.sum(w*(q-m)**2))),
            'q0_q025_q975':np.interp([.025,.975],cdf,q[order]).tolist(),
            'P_q0_negative':float(w[q<0].sum()), 'importance_ESS':float(1/np.sum(w*w)),
            'max_normalized_weight':float(w.max()), 'marginal_mode':x.tolist(),
            'mode_q0':float(.5+1.5*x[1]*(1-x[0])), 'mode_chisq':float(opt.fun),
        }
    # Algebra check: integrate the global amplitude after profiling the flat
    # intercept. Sherman-Morrison must give the same quadratic form as C+s²tt'.
    A=sn.A; At=A@t
    analytic=A-np.outer(At,At)/(1/s**2+t@At)
    tests={'shared_projector_max_difference':float(abs(analytic-projectors['shared']).max())}
    check=draws[:30,:4]
    residual=sn.mag-mu(sn.z,check[:,:3],zhel=sn.zhel)-t
    direct=np.einsum('bi,ij,bj->b',residual,projectors['shared'],residual)+bao.chisq(check[:,:3],check[:,3])
    tests['full_minus_compressed_max_chisq']=float(abs(direct-chi(check,'shared')).max())
    # Explicit MAP amplitude elimination independently checks the sign.
    explicit=[]
    for r in residual:
        delta=(t@A@r)/(1/s**2+t@A@t)
        explicit.append((r-delta*t)@A@(r-delta*t)+(delta/s)**2)
    tests['explicit_amplitude_minimum_vs_marginal_quadratic_max']=float(abs(np.array(explicit)-(direct-bao.chisq(check[:,:3],check[:,3]))).max())
    assert tests['shared_projector_max_difference']<1e-8
    assert tests['full_minus_compressed_max_chisq']<1e-4
    assert tests['explicit_amplitude_minimum_vs_marginal_quadratic_max']<1e-7
    result['checks']=tests
    result['sample_correlations']={k:float(np.corrcoef(draws[:,4],v)[0,1]) for k,v in
                                    {'q0':q,'Om':draws[:,0],'w0':draws[:,1],'wa':draws[:,2]}.items()}
    result['limitations']=[
        'Importance ESS is weight concentration, not independent-chain ESS; inherited autocorrelation remains.',
        'The shared case was already treated correctly in the main analysis.',
        'This integrates only the shared slope; clock, DTD, SFH, host-age/dust, overlap and selection uncertainty remain conditional.',
        'External slope prior independence from the reused supernova data is assumed, not demonstrated.',
        'The Gaussian covariance equivalence uses a fixed template and its normal prior; a cosmology-dependent template would require its changing normalization.',
    ]
    target=out/'shared-uncertainty.json';target.write_text(json.dumps(result,indent=2)+'\n')
    inputs=sn.inputs+bao.inputs+[template_path,chain_path,ROOT/'scripts/cosmology/core.py',Path(__file__)]
    (out/'shared-uncertainty-manifest.json').write_text(json.dumps({
        'inputs_sha256':{str(p.relative_to(ROOT)):digest(p) for p in inputs},
        'outputs_sha256':{str(target.relative_to(ROOT)):digest(target)},
        'seed':920612,'method':'Exact Gaussian marginalization plus importance reweighting of existing shared-slope posterior',
    },indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
