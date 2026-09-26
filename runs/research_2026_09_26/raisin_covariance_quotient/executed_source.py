#!/usr/bin/env python3
"""Check RAISIN covariance closure after eliminating a free common magnitude.

This is a source-reproduction diagnostic, not an outcome-fitted covariance.
Run freeze before score; no covariance is repaired by this script.
"""
import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.cosmology import FlatLambdaCDM
from scipy.linalg import helmert, eigvalsh, cho_factor, cho_solve

ROOT = Path(__file__).resolve().parents[2]
REL = ROOT / 'sources/repos/djones1040__RAISIN_DataRelease@a383c4b'
OLD = ROOT / 'runs/research_2026_09_26/raisin_differential'
OUT = ROOT / 'runs/research_2026_09_26/raisin_covariance_quotient'
SOURCE = OLD / 'mass-threshold/code/raisin_cosmo/cosmo_sys.py'
SCRIPT = Path(__file__).resolve()


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def write(p, x):
    p.write_text(json.dumps(x, indent=2, allow_nan=False) + '\n')


def readfit(p):
    lines = p.read_text().splitlines()
    cols = next(x.split()[1:] for x in lines if x.startswith('VARNAMES:'))
    rows = [x.split()[1:] for x in lines if x.startswith('SN:')]
    d = pd.DataFrame(rows, columns=cols).set_index('CID')
    for k in ['zHD', 'DLMAG', 'DLMAGERR']:
        d[k] = d[k].astype(float)
    assert d.index.is_unique
    return d


def inputs():
    paths = [SOURCE, OLD/'frozen-membership.csv', SCRIPT]
    for b in ['nir', 'optical', 'opticalnir']:
        paths += list((REL/f'distances/w/{b}_dist').glob('*.FITRES'))
        paths += [REL/f'distances/w/{b}_syst/RAISIN_all.covmat']
        paths += [REL/f'distances/w/{b}_syst/RAISIN_{g}_lcparams_cosmosis.txt'
                  for g in ['all', 'stat']]
        paths += [OLD/f'mass-threshold/code/output/cosmo_fitres_{b}/RAISIN_combined_FITOPT029_new.FITRES']
    return sorted(set(paths))


def freeze():
    OUT.mkdir(parents=True, exist_ok=True)
    p = OUT/'protocol.json'
    if p.exists():
        raise FileExistsError(p)
    write(p, dict(created_utc=datetime.now(timezone.utc).isoformat(),
        stage='Before intercept-quotient outcomes; earlier raw all-matrix failure known.',
        source_formula='Unit-weight sum of FITOPT1:28 or1:29 signed outer products; variant inverse-variance centering and fiducial mu70,.3 redshift adjustment, exactly as pinned source.',
        identification='Common SN magnitude is free. Full Gaussian shape likelihood depends on B y and B C B^T for any rank(n-1) B with B1=0. Absolute calibration and marginal-likelihood normalization are not inferred.',
        primary='For each branch and each documented option range, compare B(C_source-C_published)B^T with propagated literal export rounding bounds, B=orthonormal Helmert contrasts. No fit of weights, selection of variants, centering prescription or tolerance.',
        rounding='Off-diagonal half last printed %.5e digit; reconstructed diagonal bound (dmb_all+dmb_stat)*1e-6+1e-12; propagate absolute element bounds as abs(B) T abs(B)^T plus1e-12 numerical tolerance.',
        checks='Helmert orthonormality/null; agreement with n-by-n centered H representation; equality of profiled-intercept inverse quadratic form with contrast likelihood; common-offset injection leaves contrast covariance unchanged.',
        descriptive='Relative centered systematic Frobenius mismatch, change in frozen high37-minus-low42 variance, and generalized eigenvalues of mismatch relative to total published contrast covariance. These diagnose nonclosure, not real covariance error.',
        prohibition='No repaired all matrix, no cosmology fit, no inference that source reconstruction is more correct than the published export.',
        input_sha256={str(x.relative_to(ROOT)): sha(x) for x in inputs()}))
    print(json.dumps({'protocol_sha256': sha(p)}))


def score():
    p = OUT/'protocol.json'
    protocol = json.loads(p.read_text())
    for name, h in protocol['input_sha256'].items():
        assert sha(ROOT/name) == h, name
    if (OUT/'result.json').exists():
        raise FileExistsError('Keep original score')
    members = pd.read_csv(OLD/'frozen-membership.csv').set_index('CID')
    n = len(members)
    B = helmert(n)
    H = np.eye(n)-np.ones((n,n))/n
    high = (members.stratum == 'high').to_numpy()
    a = high/high.sum() - (~high)/(~high).sum()
    assert n == 79 and high.sum() == 37
    assert np.max(abs(B@B.T-np.eye(n-1))) < 1e-14
    assert np.max(abs(B.T@B-H)) < 1e-14
    cosm = FlatLambdaCDM(H0=70, Om0=.3)
    results, arrays = {}, {'B': B, 'a': a}
    for branch in ['nir', 'optical', 'opticalnir']:
        def tab(k):
            q = REL/f'distances/w/{branch}_dist/RAISIN_combined_FITOPT{k:03d}.FITRES'
            if k == 29:
                q = OLD/f'mass-threshold/code/output/cosmo_fitres_{branch}/RAISIN_combined_FITOPT029_new.FITRES'
            return readfit(q).loc[members.index]
        base = tab(0)
        mu0 = cosm.distmod(base.zHD.to_numpy()).value
        sysp = REL/f'distances/w/{branch}_syst'
        raw = np.loadtxt(sysp/'RAISIN_all.covmat')
        assert int(raw[0]) == n
        off = raw[1:].reshape(n,n)
        allerr = np.loadtxt(sysp/'RAISIN_all_lcparams_cosmosis.txt')[:,5]
        staterr = np.loadtxt(sysp/'RAISIN_stat_lcparams_cosmosis.txt')[:,5]
        pub = off + np.diag(allerr**2-staterr**2)
        total = off + np.diag(allerr**2)
        tol = np.full((n,n), 1e-12)
        nz = off != 0
        tol[nz] += .5*10.**(np.floor(np.log10(abs(off[nz])))-5)
        np.fill_diagonal(tol, (allerr+staterr)*1e-6+1e-12)
        btol = abs(B)@tol@abs(B.T) + 1e-12
        Cb = B@total@B.T
        inv = cho_solve(cho_factor(total), np.eye(n))
        u = inv@np.ones(n)
        P = inv-np.outer(u,u)/u.sum()
        P2 = B.T@cho_solve(cho_factor(Cb),B)
        profile_err = float(np.max(abs(P-P2)))
        assert profile_err < 1e-9
        v = np.linspace(-.01,.01,n)
        injected = total + np.outer(np.ones(n), v)+np.outer(v,np.ones(n))
        offset_err = float(np.max(abs(B@injected@B.T-Cb)))
        assert offset_err < 1e-14
        recon = np.zeros((n,n))
        br = dict(profile_precision_identity_max=profile_err,offset_injection_error=offset_err,options={})
        for k in range(1,30):
            d = tab(k)
            dm = d.DLMAG.to_numpy()-base.DLMAG.to_numpy()
            r = dm-np.average(dm,weights=1/d.DLMAGERR.to_numpy()**2)
            r -= cosm.distmod(d.zHD.to_numpy()).value-mu0
            recon += np.outer(r,r)
            if k not in [28,29]:
                continue
            delta = recon-pub
            db = B@delta@B.T
            hdelta = H@delta@H
            identity = float(np.max(abs(B.T@db@B-hdelta)))
            assert identity < 1e-12
            ev = eigvalsh(db,Cb)
            br['options'][str(k)] = dict(pass_rounding_gate=bool(np.all(abs(db)<=btol)),
                max_absolute_contrast_covariance_error_mag2=float(abs(db).max()),
                max_error_over_propagated_rounding_bound=float(np.max(abs(db)/btol)),
                relative_centered_systematic_frobenius_error=float(np.linalg.norm(db)/np.linalg.norm(B@pub@B.T)),
                published_highlow_sys_variance_mag2=float(a@pub@a),
                reconstructed_highlow_sys_variance_mag2=float(a@recon@a),
                highlow_variance_difference_mag2=float(a@delta@a),
                generalized_total_covariance_error_eigenvalue_range=[float(ev[0]),float(ev[-1])],
                centered_representation_identity_max_error=identity)
            arrays[f'{branch}_{k}_contrast_difference']=db
        results[branch] = br
    np.savez_compressed(OUT/'arrays.npz', **arrays)
    write(OUT/'result.json', dict(protocol_sha256=sha(p), branches=results,
        scientific_scope='Failed source reproduction even in observable distance-shape subspace is not evidence that the supplied covariance itself is wrong.',
        arrays_sha256=sha(OUT/'arrays.npz')))
    (OUT/'executed_source.py').write_bytes(SCRIPT.read_bytes())
    print(json.dumps(results,indent=2))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('mode',choices=['freeze','score'])
    args = ap.parse_args(); globals()[args.mode]()
