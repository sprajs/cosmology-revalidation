"""Read-only BAO release audit and posterior covariance sensitivity.

Writes only runs/assumption_audit/bao/. Existing fits and source files are untouched.
Correlation stress tests preserve the supplied diagonal and within-tracer blocks;
they add ONLY cross-tracer theoretical systematic covariance at published scales.
"""
from pathlib import Path
import hashlib
import json
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.linalg import cho_factor, cho_solve
import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts/cosmology'))
from core import BAO

OUT = ROOT / 'runs/assumption_audit/bao'
OUT.mkdir(parents=True, exist_ok=True)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def independent_prediction(zs, kinds, pars):
    om, w0, wa, hrd = pars
    def ez(z):
        return np.sqrt(om * (1 + z)**3 + (1 - om) *
                       (1 + z)**(3 * (1 + w0 + wa)) * np.exp(-3 * wa * z / (1 + z)))
    ans = []
    for z, kind in zip(zs, kinds):
        dm = 299792.458 / hrd * quad(lambda zz: 1 / ez(zz), 0, z,
                                    epsabs=1e-12, epsrel=1e-12)[0]
        dh = 299792.458 / hrd / ez(z)
        ans.append({'DM_over_rs': dm, 'DH_over_rs': dh,
                    'DV_over_rs': (z * dm**2 * dh)**(1/3)}[kind])
    return np.array(ans)


def q_summary(x, weight):
    q = .5 + 1.5 * x[:, 1] * (1 - x[:, 0])
    w = weight / weight.sum()
    mean = float(w @ q)
    return dict(q0_mean=mean, q0_sd=float(np.sqrt(w @ ((q - mean)**2))),
                p_q0_negative=float(np.clip(w @ (q < 0), 0., 1.)),
                weight_ess_fraction=float(1 / (len(w) * (w @ w))))


def main():
    b = BAO()
    inputs = [*b.inputs, ROOT / 'scripts/cosmology/core.py']
    cov = np.loadtxt(b.inputs[1])
    kinds = b.data.kind.to_numpy()
    std = np.sqrt(np.diag(cov))
    corr = cov / np.outer(std, std)
    zs = np.unique(b.z)
    blocks = []
    for z in zs:
        idx = np.flatnonzero(b.z == z)
        blocks.append(dict(z=float(z), indices=idx.tolist(), kinds=kinds[idx].tolist(),
                           values=b.value[idx].tolist(), std=std[idx].tolist(),
                           correlation=corr[np.ix_(idx, idx)].tolist()))
    assert set(kinds) == {'DV_over_rs', 'DM_over_rs', 'DH_over_rs'}
    assert cov.shape == (13, 13)
    assert np.array_equal(cov, cov.T)
    cross = b.z[:, None] != b.z[None, :]
    assert not np.any(cov[cross])
    # Compare individual release blocks after identifying their observable rows.
    subchecks = []
    for path in sorted(b.inputs[0].parent.glob('*_mean.txt')):
        if '_ALL_' in path.name:
            continue
        d = pd.read_csv(path, sep=r'\s+', comment='#', names=['z', 'value', 'kind'])
        cc = path.with_name(path.name.replace('_mean', '_cov'))
        inputs += [path, cc]
        indices = [int(np.flatnonzero((b.z == row.z) & (kinds == row.kind))[0])
                   for row in d.itertuples()]
        value_error = float(np.max(np.abs(b.value[indices] - d.value.to_numpy())))
        cov_error = float(np.max(np.abs(cov[np.ix_(indices, indices)] - np.atleast_2d(np.loadtxt(cc)))))
        # ALL prints eight decimal digits; the Ly-alpha-only file prints more.
        assert value_error == 0 and cov_error < 3e-10
        subchecks.append(dict(file=str(path.relative_to(ROOT)), indices=indices,
                              max_value_difference=value_error, max_covariance_difference=cov_error))
    checks = []
    rng = np.random.default_rng(221091)
    pars = [[.3, -1, 0, 10000], [.353, -.42, -1.75, 9500]]
    for _ in range(12):
        pars.append([rng.uniform(.1, .6), rng.uniform(-1.5, -.1),
                     rng.uniform(-2, 0), rng.uniform(8000, 11000)])
    for p in pars:
        pred = independent_prediction(b.z, kinds, p)
        residual = b.value - pred
        direct_chisq = float(residual @ cho_solve(cho_factor(cov, lower=True), residual))
        checks.append(dict(theta=p, max_prediction_difference=float(np.max(np.abs(pred - b.prediction(p[:3], p[3])[0]))),
                           chisq_difference=float(b.chisq(p[:3], p[3])[0] - direct_chisq)))
    assert max(x['max_prediction_difference'] for x in checks) < 1e-8
    assert max(abs(x['chisq_difference']) for x in checks) < 1e-7
    # A redshift-independent BAO scale shift is absorbed by the free H0*rd.
    theta, hrd, scale = [.3, -.8, -.5], 10000., 1.015
    residual = b.value - b.prediction(theta, hrd)[0]
    shifted = scale * b.value - b.prediction(theta, hrd / scale)[0]
    rescale_error = float(shifted @ np.linalg.solve(scale**2 * cov, shifted) -
                          residual @ np.linalg.solve(cov, residual))
    assert abs(rescale_error) < 1e-9
    # In log-distance variables: ln DM = ln alpha_iso - ln alpha_AP / 3;
    # ln DH = ln alpha_iso + 2 ln alpha_AP / 3. BGS has no AP observation.
    galaxy = b.z < 2
    uiso = .001 * b.value * galaxy
    uap = .002 * b.value * np.select([kinds == 'DM_over_rs', kinds == 'DH_over_rs'], [-1/3, 2/3], default=0) * galaxy
    cross_theory = (np.outer(uiso, uiso) + np.outer(uap, uap)) * cross
    covariants = {'baseline': cov}
    for rho in [.5, 1.]:
        cc = cov + rho * cross_theory
        assert np.all(np.linalg.eigvalsh(cc) > 0)
        assert np.array_equal(np.diag(cc), np.diag(cov))
        covariants[f'galaxy_theory_rho_{rho:g}'] = cc
    # Deliberately incorrect diagonal-only control, clearly separate from valid alternatives.
    covariants['diagonal_only_control_NOT_recommended'] = np.diag(np.diag(cov))
    sensitivity = {}
    for name in ['bao-cpl', 'pantheon-bao-cpl', 'pantheon-bao-cpl-c14fixed', 'pantheon-bao-cpl-c14slope']:
        path = ROOT / 'runs/cosmology' / name / 'chains.npz'
        inputs.append(path)
        chain = np.load(path)['chain'][::10]
        flat = chain.reshape(-1, chain.shape[-1])
        residual = np.concatenate([b.value - b.prediction(t[:, :3], t[:, 3])
                                   for t in np.array_split(flat, 50)])
        old_chi = np.einsum('bi,ij,bj->b', residual, np.linalg.inv(cov), residual)
        sensitivity[name] = {'retained_shape': list(chain.shape), 'variants': {}}
        base = q_summary(flat, np.ones(len(flat)))
        for variant, cc in covariants.items():
            new_chi = np.einsum('bi,ij,bj->b', residual, np.linalg.inv(cc), residual)
            logw = -.5 * (new_chi - old_chi)
            weights = np.exp(logw - logw.max())
            summary = q_summary(flat, weights)
            summary['q0_mean_shift'] = summary['q0_mean'] - base['q0_mean']
            # Paired 20-time-block difference MCSE (walker dependence retained).
            shifts = []
            probabilities = []
            for ids in np.array_split(np.arange(len(chain)), 20):
                sel = (ids[:, None] * chain.shape[1] + np.arange(chain.shape[1])).ravel()
                aa = q_summary(flat[sel], weights[sel]); bb = q_summary(flat[sel], np.ones(len(sel)))
                shifts.append(aa['q0_mean'] - bb['q0_mean'])
                probabilities.append(aa['p_q0_negative'] - bb['p_q0_negative'])
            summary['paired_q0_shift_block_mcse'] = float(np.std(shifts, ddof=1) / np.sqrt(20))
            summary['p_shift'] = summary['p_q0_negative'] - base['p_q0_negative']
            summary['paired_p_shift_block_mcse'] = float(np.std(probabilities, ddof=1) / np.sqrt(20))
            sensitivity[name]['variants'][variant] = summary
        q = .5 + 1.5 * flat[:, 1] * (1 - flat[:, 0])
        # Exact prior-coordinate sensitivity inside the same finite hrd support.
        for exponent, label in [(-1, 'uniform_log_H0rd'), (-2, 'uniform_inverse_H0rd')]:
            sensitivity[name][label] = q_summary(flat, flat[:, 3]**exponent)
        sensitivity[name]['wa_below_minus_2_8'] = float(np.mean(flat[:, 2] < -2.8))
        sensitivity[name]['q0_H0rd_correlation'] = float(np.corrcoef(q, flat[:, 3])[0, 1])
    # The reference likelihood configurations encode foreground/nuisance assumptions.
    cmb = []
    for d in sorted((ROOT / 'data/bao/desi-dr2-reference/cobaya/base_w_wa').iterdir()):
        config = d / 'chain.updated.yaml'
        inputs.append(config)
        y = yaml.safe_load(config.read_text())
        pieces = []
        for path in sorted(d.glob('chain.[1-4].txt')):
            inputs.append(path)
            columns = path.open().readline().strip().lstrip('#').split()
            a = np.loadtxt(path); a = a[int(.3 * len(a)):]
            pieces.append(pd.DataFrame(a, columns=columns))
        df = pd.concat(pieces)
        # Today massive neutrinos are approximated as matter as in original q derivation.
        orad = 1 - df.omegam.to_numpy() - df.omegal.to_numpy()
        dq = (.5 - 1.5 * df.w.to_numpy()) * orad
        qapprox = .5 + 1.5 * df.w.to_numpy() * (1 - df.omegam.to_numpy())
        weights = df.weight.to_numpy()
        cmb.append(dict(directory=str(d.relative_to(ROOT)),
                        likelihoods=list(y['likelihood']),
                        fixed_parameters={k: y['params'][k] for k in ['mnu', 'nnu', 'omk']},
                        theory=y['theory'],
                        q0_radiation_correction_mean=float(np.average(dq, weights=weights)),
                        q0_radiation_correction_max=float(dq.max()),
                        weighted_sign_probability_change=float(np.average((qapprox + dq) < 0, weights=weights) - np.average(qapprox < 0, weights=weights))))
    result = dict(created_utc=datetime.now(timezone.utc).isoformat(),
                  scope='Released compressed BAO likelihood and existing posterior sensitivity; no raw-catalog refit, no CMB likelihood rerun.',
                  data_integrity=dict(n=13, blocks=blocks, min_covariance_eigenvalue=float(np.linalg.eigvalsh(cov).min()),
                                      cross_redshift_nonzero_entries=int(np.count_nonzero(cov[cross])), individual_block_matches=subchecks),
                  independent_quadrature_and_cholesky=checks,
                  common_1p5_percent_scale_chisq_invariance_error=rescale_error,
                  covariance_stress_definition='Only off-block theoretical systematics correlated: 0.1% alpha_iso and 0.2% alpha_AP in the six galaxy/quasar blocks; Ly-alpha has a separate theory budget and is left untouched. Same diagonal errors. Linearized Jacobian at release means.',
                  posterior_sensitivity=sensitivity, cmb_configurations=cmb)
    dest = OUT / 'release_check.json'
    dest.write_text(json.dumps(result, indent=2) + '\n')
    archive = OUT / 'provenance'; archive.mkdir(exist_ok=True)
    for source in [Path(__file__), ROOT / 'scripts/cosmology/core.py']:
        (archive / (sha(source) + '.py')).write_bytes(source.read_bytes())
    (OUT / 'manifest.json').write_text(json.dumps(dict(
        created_utc=result['created_utc'], script_sha256=sha(__file__),
        inputs_sha256={str(p.relative_to(ROOT)): sha(p) for p in dict.fromkeys(inputs)},
        outputs_sha256={str(dest.relative_to(ROOT)): sha(dest)}), indent=2) + '\n')
    print(json.dumps(dict(n=13, max_quadrature_error=max(x['max_prediction_difference'] for x in checks),
                          covariance_sensitivity=sensitivity), indent=2))


if __name__ == '__main__':
    main()
