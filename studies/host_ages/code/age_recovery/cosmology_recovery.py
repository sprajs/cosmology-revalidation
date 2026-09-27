#!/usr/bin/env python3
"""Observed-design injection into actual nonlinear cosmology fits, not survey simulation."""
from pathlib import Path
import hashlib
import json
import sys
from datetime import datetime, timezone
import numpy as np
import pandas as pd
from scipy.linalg import solve_triangular
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
from lib.cosmology import mu, Pantheon, BAO, qvalue

HERE = Path(__file__).resolve().parent
OUT = ROOT / 'studies/host_ages/results/age_recovery'
WORK = ROOT / '.work/age-recovery'


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def save(p, x):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(x, indent=2, allow_nan=False) + '\n')


def design():
    rows = []
    for i in (1, 2):
        for line in (ROOT / f'data/ages/table{i}.dat').read_text().splitlines()[1:]:
            c = [v.strip().replace('\\', '').strip() for v in line.split('&')]
            if len(c) == 5:
                rows.append(dict(CID=str(int(c[0])), age=float(c[1]), age_err=float(c[2]), age_source=i))
    ages = pd.DataFrame(rows).drop_duplicates('CID', keep='first')
    pp = pd.read_csv(ROOT / 'data/distances/Pantheon+SH0ES.dat', sep=r'\s+', dtype={'CID': str})
    pp['source_row'] = np.arange(len(pp))
    d = ages.merge(pp[pp.IDSURVEY == 1], on='CID', validate='one_to_one')
    d = d[(d.zHD > .06) & (d.zHD < .42)].copy()
    assert len(d) == d.CID.nunique() == 196
    raw = np.loadtxt(ROOT / 'data/distances/Pantheon+SH0ES_STAT+SYS.cov')
    n = int(raw[0]); cov = raw[1:].reshape(n, n)
    cov = (cov + cov.T) / 2
    ix = d.source_row.to_numpy()
    return d, cov[np.ix_(ix, ix)]


def summary(v):
    v = np.asarray(v)
    return dict(mean=float(v.mean()), sd=float(v.std(ddof=1)), mcse=float(v.std(ddof=1) / np.sqrt(len(v))), q025_q975=np.quantile(v, [.025, .975]).tolist())


class Fit:
    """Profile the linear predictors for every actual nonlinear FLRW distance."""
    def __init__(self, z, zhel, L, X):
        self.L = L
        self.X = solve_triangular(L, X, lower=True)
        self.Q = np.linalg.qr(self.X, mode='reduced')[0]
        self.Xplus = np.linalg.solve(self.X.T @ self.X, self.X.T)
        self.grid = np.linspace(.001, .999, 1001)
        actual = mu(z, self.grid[:, None], 'lcdm', zhel).T
        self.spline = CubicSpline(self.grid, solve_triangular(L, actual, lower=True), axis=1)
        # Cubic interpolation is solely numerical acceleration of the FLRW integral.
        probe = np.array([.0017, .21357, .67311, .9982])
        direct = solve_triangular(L, mu(z, probe[:, None], 'lcdm', zhel).T, lower=True)
        self.interpolation_error = float(np.max(abs(self.spline(probe) - direct)))
        assert self.interpolation_error < 1e-8

    def residual(self, om, y):
        r = y - self.spline(om)
        b = self.Xplus @ r
        return r - self.X @ b, b

    def fit(self, y):
        # Grid ensures global branch selection; refine the best interval.
        r = y[:, None] - self.spline(self.grid)
        r -= self.Q @ (self.Q.T @ r)
        k = int(np.argmin(np.sum(r*r, axis=0)))
        lo, hi = self.grid[max(0, k-1)], self.grid[min(len(self.grid)-1, k+1)]
        sol = minimize_scalar(lambda o: np.sum(self.residual(o, y)[0]**2), bounds=(lo, hi), method='bounded', options={'xatol': 1e-10})
        assert sol.success
        resid, b = self.residual(sol.x, y)
        J = np.column_stack([self.spline(sol.x, 1), self.X])
        C = np.linalg.inv(J.T @ J)
        return sol.x, b, C, resid, bool(sol.x < .00101 or sol.x > .99899)


def precision():
    sn, bao = Pantheon(), BAO()
    # Fixed fiducial for a prospective response calculation, not fitted correction.
    t = np.array([.3, -1., 0., 10000.])
    def derivative(fun, p, h):
        return np.column_stack([(fun(p + np.eye(len(p))[j]*h[j]) - fun(p - np.eye(len(p))[j]*h[j]))/(2*h[j]) for j in range(len(p))])
    h = np.array([1e-4, 1e-4, 1e-4, 1.])
    D = derivative(lambda p: mu(sn.z, p[:3], 'cpl', sn.zhel)[0], t, h)
    B = derivative(lambda p: bao.prediction(p[:3], p[3])[0], t, h)
    C = np.linalg.inv(D.T @ sn.A @ D + B.T @ bao.inv @ B)
    # g has unit contrast from .05 to 1; intercept is removed by sn.A.
    g = (sn.z/(1+sn.z) - .05/1.05) / (.5-.05/1.05)
    response = C @ D.T @ sn.A @ g
    Jq = np.array([1.5, 1.5*(1-t[0]), 0, 0])
    qsig = float(np.sqrt(Jq@C@Jq)); qr = float(Jq@response)
    report = dict(model='local Fisher SN+BAO flat CPL; no CMB and no posterior priors', fiducial=t.tolist(), parameter_order=['Omega_m','w0','wa','H0_rd'], covariance=C.tolist(), bias_response_per_mag=response.tolist(), q0_sigma=qsig, q0_response_per_mag=qr, allowed_template_contrast_mag_for_0_1sigma={'q0': .1*qsig/abs(qr), 'w0': .1*np.sqrt(C[1,1])/abs(response[1]), 'wa': .1*np.sqrt(C[2,2])/abs(response[2])})
    d = (mu(sn.z, [.3001], 'lcdm', sn.zhel)[0] - mu(sn.z, [.2999], 'lcdm', sn.zhel)[0])/.0002
    var = 1/(d@sn.A@d); r = var*d@sn.A@g
    report['SN_flatLCDM'] = dict(Om_sigma=float(np.sqrt(var)), Om_response_per_mag=float(r), q0_sigma=float(1.5*np.sqrt(var)), q0_response_per_mag=float(1.5*r), allowed_template_contrast_mag_for_0_1sigma=float(.1*np.sqrt(var)/abs(r)))
    report['warning'] = 'Shape-specific local precision budget, not a data-supported bound on age evolution or a full posterior result.'
    return report


def main():
    p = json.loads((HERE/'cosmology-protocol.json').read_text())
    d, C = design(); L = np.linalg.cholesky(C)
    d.to_csv(WORK/'observed-design.csv', index=False)
    z, zh = d.zHD.to_numpy(), d.zHEL.to_numpy()
    A = d.age.to_numpy(); A -= A.mean()
    aw = solve_triangular(L, A, lower=True)
    onew = solve_triangular(L, np.ones(len(d)), lower=True)
    ar = aw - onew*(onew@aw)/(onew@onew)
    slope_row = ar/(ar@ar)
    rng = np.random.default_rng(p['seed'])
    eps = rng.normal(size=(p['realizations'], len(d)))
    baseline = solve_triangular(L, mu(z, [.3], 'lcdm', zh)[0], lower=True)
    out = dict(scope=p['selection_scope'], known_injection_age_not_host_likelihood=True, n=len(d), seed=p['seed'], realizations=len(eps), experiments=[], precision_budget=precision())
    saved = []
    for kind in ['intercept_only', 'refit_width_colour_mass']:
        X = np.ones((len(d), 1))
        if kind != 'intercept_only':
            X = np.column_stack([X, d.x1, d.c, (d.HOST_LOGMASS > 10).astype(float)])
        std = Fit(z, zh, L, X)
        joint = Fit(z, zh, L, np.column_stack([X, A]))
        # Reference age/low-z zero point cancels with the freely fitted intercept.
        for btrue in p['injected_slopes_mag_per_Gyr']:
            draws = []
            for i, e in enumerate(eps):
                y = baseline + btrue*aw + e
                om, coeff, cc, r, edge = std.fit(y)
                oj, cj, cv, rj, edgej = joint.fit(y)
                bs = float(slope_row@r)
                # Appending full historical correction after standardization, without refitting nuisances.
                adjusted = y - std.X@coeff - (-.03)*aw
                # Refit intercept and cosmology only, because already-standardized nuisances are frozen.
                # Faster equivalent to a prebuilt intercept-only profile.
                if kind == 'intercept_only':
                    template_fit = std
                else:
                    if i == 0 and btrue == p['injected_slopes_mag_per_Gyr'][0]:
                        template_fit = Fit(z, zh, L, np.ones((len(d),1)))
                ot, _, _, _, edget = template_fit.fit(adjusted)
                se_age = float(np.sqrt(cv[-1,-1])); se_om = float(np.sqrt(cv[0,0]))
                draws.append([om, bs, oj, cj[-1], se_age, se_om, float(abs(cj[-1]-btrue)<=1.96*se_age), float(abs(oj-.3)<=1.96*se_om), float(edge), float(edgej), ot, float(edget)])
            a = np.asarray(draws)
            row = dict(model=kind, injected_slope=btrue, standard_Om=summary(a[:,0]), postfit_residual_slope=summary(a[:,1]), joint_Om=summary(a[:,2]), joint_age=summary(a[:,3]), nominal_age_interval_mean_se=float(a[:,4].mean()), age_coverage=float(a[:,6].mean()), Om_coverage=float(a[:,7].mean()), standard_boundary_rate=float(a[:,8].mean()), joint_boundary_rate=float(a[:,9].mean()), appended_full_template_Om=summary(a[:,10]), template_boundary_rate=float(a[:,11].mean()), false_positive_or_power=float(np.mean(abs(a[:,3])>1.96*a[:,4])))
            row['coverage_mcse'] = float(np.sqrt(row['age_coverage']*(1-row['age_coverage'])/len(a)))
            # Noiseless response separates absorption from random realization changes.
            om0, _, _, r0, _ = std.fit(baseline + btrue*aw)
            row['noiseless'] = dict(standard_Om=om0, postfit_residual_slope=float(slope_row@r0))
            out['experiments'].append(row)
            for i, values in enumerate(draws): saved.append([kind,btrue,i]+values)
            print(kind, btrue, row['postfit_residual_slope']['mean'], row['joint_age']['mean'], row['age_coverage'], flush=True)
        # Linear geometry is a separate independent check of the small-injection limit.
        der = std.spline(.3, 1); Xw=std.X
        J = np.column_stack([Xw,der]); Q=np.linalg.qr(J,mode='reduced')[0]
        ageperp=aw-Q@(Q.T@aw)
        gain=float(slope_row@ageperp)
        interval=np.sqrt(1/(ageperp@ageperp))
        out.setdefault('geometry',{})[kind]=dict(postfit_slope_gain=gain, joint_known_age_se=float(interval), power_minus_0_03=float(norm.cdf(-1.96+.03/interval)+norm.sf(1.96+.03/interval)), interpolation_error=std.interpolation_error)
    bins = np.digitize(z, [.10,.15,.20,.25,.30,.35])
    V = np.column_stack([bins==i for i in np.unique(bins)]).astype(float)
    Q=np.linalg.qr(solve_triangular(L,V,lower=True),mode='reduced')[0]
    out['bin_center_surrogate_gain']=float(slope_row@(aw-Q@(Q.T@aw)))
    out['validation']={'correlated_noise_generation':'whitened draws equivalent to L@epsilon; C is full released covariance', 'inverse_covariance_identity_max':float(np.max(abs(np.linalg.solve(C,C)-np.eye(len(d)))))}
    cols=['model','injected_slope','draw','standard_Om','postfit_slope','joint_Om','joint_age','joint_age_se','joint_Om_se','age_coverage','Om_coverage','standard_boundary','joint_boundary','appended_Om','appended_boundary']
    pd.DataFrame(saved,columns=cols).to_csv(WORK/'cosmology-draws.csv',index=False)
    inputs=[ROOT/'data/ages/table1.dat', ROOT/'data/ages/table2.dat',ROOT/'data/distances/Pantheon+SH0ES.dat',ROOT/'data/distances/Pantheon+SH0ES_STAT+SYS.cov']+list((ROOT/'data/bao').glob('*.txt'))
    out['provenance']=dict(time_utc=datetime.now(timezone.utc).isoformat(), input_sha256={str(x.relative_to(ROOT)):digest(x) for x in inputs}, code_sha256={str(x.relative_to(ROOT)):digest(x) for x in [Path(__file__),ROOT/'lib/cosmology.py']}, protocol_sha256=digest(HERE/'cosmology-protocol.json'), draw_sha256=digest(WORK/'cosmology-draws.csv'))
    save(OUT/'cosmology-recovery.json',out)

if __name__=='__main__':
    WORK.mkdir(parents=True,exist_ok=True)
    main()
