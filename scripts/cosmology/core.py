"""Independent, low-redshift metric-distance and released-product likelihoods.

No CMB, raw photometry, luminosity evolution, or early-universe model is implicit.
Distances here are in units c/H0; a common SN magnitude offset is marginalized.
"""
from pathlib import Path
import hashlib
import json
import platform
import subprocess
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
from scipy.linalg import cho_factor, cho_solve
from numpy.polynomial.legendre import leggauss

ROOT = Path(__file__).resolve().parents[2]
EDGES = np.array([0., .1, .3, .6, 1., 2.5])
GX, GW = leggauss(48)
# Capture source at import, not at a long fit's completion while other work proceeds.
CODE_AT_START = {str(p.relative_to(ROOT)): p.read_bytes() for p in Path(__file__).parent.glob('*.py')}
REVISION_AT_START = subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
DIRTY_AT_START = bool(subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip())


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def manifest(directory, purpose, inputs, configuration, outputs):
    files = {str(Path(p).relative_to(ROOT)): sha(p) for p in inputs}
    code = {p: hashlib.sha256(b).hexdigest() for p,b in CODE_AT_START.items()}
    archive=ROOT/'runs/cosmology/provenance'; archive.mkdir(parents=True,exist_ok=True)
    for p,b in CODE_AT_START.items():
        dest=archive/(code[p]+'.py')
        if not dest.exists(): dest.write_bytes(b)
        assert sha(dest)==code[p]
    record = dict(created_utc=datetime.now(timezone.utc).isoformat(), purpose=purpose,
                  inputs_sha256=files, code_sha256=code, configuration=configuration,
                  git_revision=REVISION_AT_START, git_dirty=DIRTY_AT_START,
                  code_capture='at process import; exact bytes under runs/cosmology/provenance/<sha256>.py',
                  python=platform.python_version(), lock_sha256=sha(ROOT/'uv.lock'),
                  outputs_sha256={str(Path(p).relative_to(ROOT)): sha(p) for p in outputs})
    (directory/'manifest.json').write_text(json.dumps(record, indent=2)+'\n')


def efunc(z, theta, model='cpl'):
    """Vectorized E(z)=H(z)/H0. Last axis of theta contains parameters."""
    t = np.atleast_2d(theta)
    z = np.asarray(z)
    if model == 'kinematic':
        q0, q1 = [t[:, k].reshape((-1,)+(1,)*z.ndim) for k in range(2)]
        return (1+z)**(1+q0+q1)*np.exp(-q1*z/(1+z))
    om = t[:, 0].reshape((-1,)+(1,)*z.ndim)
    w0 = -1. if model == 'lcdm' else t[:, 1].reshape((-1,)+(1,)*z.ndim)
    wa = t[:, 2].reshape((-1,)+(1,)*z.ndim) if model == 'cpl' else 0.
    return np.sqrt(om*(1+z)**3 + (1-om)*(1+z)**(3*(1+w0+wa))*np.exp(-3*wa*z/(1+z)))


def integral(z, theta, model='cpl'):
    """Integral dz/E, with exact power-law integration for piecewise q."""
    z = np.atleast_1d(z)
    t = np.atleast_2d(theta)
    if model == 'qbins':
        ans = np.zeros((len(t),len(z)))
        estart = np.ones(len(t))
        for k, (lo, hi) in enumerate(zip(EDGES[:-1], EDGES[1:])):
            v = np.maximum(np.minimum(z,hi),lo)
            lr = np.log((1+v)/(1+lo))
            q = t[:, k, None]
            # expm1(-q*x)/(-q), analytic q=0 continuation.
            f = np.empty((len(t),len(z)))
            np.divide(np.expm1(-q*lr), -q, out=f, where=np.abs(q)>1e-10)
            f = np.where(np.abs(q)>1e-10, f, lr)
            ans += (1+lo)/estart[:,None]*f
            estart *= ((1+hi)/(1+lo))**(1+t[:,k])
        return ans
    zz = z[:,None]*(GX+1)/2
    return np.sum(GW/efunc(zz,t,model),axis=-1)*z/2


def mu(z, theta, model='cpl', zhel=None):
    if zhel is None:
        zhel = z
    return 5*np.log10((1+np.asarray(zhel))*integral(z,theta,model))


def qvalue(z, theta, model='cpl'):
    z = np.atleast_1d(z)
    t = np.atleast_2d(theta)
    if model == 'kinematic':
        return t[:,0,None]+t[:,1,None]*z/(1+z)
    if model == 'qbins':
        return t[:,np.minimum(np.searchsorted(EDGES[1:],z,side='right'),4)]
    om = t[:,0,None]
    w0 = -1. if model == 'lcdm' else t[:,1,None]
    wa = t[:,2,None] if model == 'cpl' else 0.
    de = (1-om)*(1+z)**(3*(1+w0+wa))*np.exp(-3*wa*z/(1+z))
    return .5*(om*(1+z)**3+(1+3*(w0+wa*z/(1+z)))*de)/(om*(1+z)**3+de)


def interpolation_basis(z, per_segment=40):
    """Separate cubic interpolants at q-bin edges avoid crossing a q jump."""
    nodes, mats = [], []
    for lo, hi in zip(EDGES[:-1], EDGES[1:]):
        active = (z >= lo) & (z <= hi)
        if not np.any(active):
            continue
        a = max(float(z.min()),lo,1e-5)
        b = min(float(z.max()),hi)
        g = np.linspace(a,b,per_segment)
        mat = np.zeros((len(z),per_segment))
        # Assign exact bin boundaries once, to the interval on their right.
        active &= (z < hi) | (hi==EDGES[-1])
        mat[active] = CubicSpline(g,np.eye(per_segment))(z[active])
        nodes.extend(g)
        mats.append(mat)
    return np.asarray(nodes), np.concatenate(mats,axis=1)


class Pantheon:
    def __init__(self, zmax=None, per_segment=40, magnitude_table=None):
        base = ROOT/'sources/repos/CobayaSampler__sn_data/PantheonPlus'
        self.inputs = [base/'Pantheon+SH0ES.dat',base/'Pantheon+SH0ES_STAT+SYS.cov']
        df = pd.read_csv(self.inputs[0],sep=r'\s+')
        self.magnitude_revision_rows = 0
        if magnitude_table is not None:
            new = pd.read_csv(magnitude_table,sep=r'\s+')
            assert len(new)==len(df)
            for col in ['CID','IDSURVEY']:
                assert np.array_equal(new[col].to_numpy(),df[col].to_numpy()), col
            for col in ['zHD','zHEL']:
                assert np.allclose(new[col].to_numpy(),df[col].to_numpy(),rtol=0,atol=1e-14), col
            self.magnitude_revision_rows = int(np.count_nonzero(np.abs(new.m_b_corr.to_numpy()-df.m_b_corr.to_numpy())>1e-10))
            df['m_b_corr'] = new.m_b_corr.to_numpy()
            self.inputs.append(Path(magnitude_table))
        raw = np.loadtxt(self.inputs[1])
        n = int(raw[0]); assert n==len(df) and raw.size==1+n*n
        cov = raw[1:].reshape(n,n)
        self.input_asymmetry = float(np.max(np.abs(cov-cov.T)))
        assert self.input_asymmetry < 1e-7, 'Inspect material covariance asymmetry'
        # Printed source entries have <=3e-8 mag^2 antisymmetric round-off.
        cov = (cov+cov.T)/2
        mask = df.zHD.to_numpy()>.01
        if zmax is not None:
            mask &= df.zHEL.to_numpy()<=zmax
        self.original_indices = np.flatnonzero(mask)
        self.df = df.loc[mask].copy()
        self.z = self.df.zHD.to_numpy()
        self.zhel = self.df.zHEL.to_numpy()
        self.mag = self.df.m_b_corr.to_numpy()
        self.cov = cov[np.ix_(mask,mask)]
        self.chol = cho_factor(self.cov,lower=True)
        inv = cho_solve(self.chol,np.eye(len(self.z)))
        u = inv.sum(axis=1)
        self.offset_precision = u.sum()
        self.A = inv-np.outer(u,u)/u.sum()
        self.nodes, self.P = interpolation_basis(self.z,per_segment)
        self.fid = mu(self.z,[.3],'lcdm',self.zhel)[0]
        self.fidnodes = mu(self.nodes,[.3],'lcdm')[0]
        self.template = np.zeros(len(self.z))
        self.reset_template(self.template)

    def reset_template(self,template):
        self.template = np.asarray(template)
        B = np.column_stack([self.P,self.template])
        y = self.mag-self.fid
        self.yAy = y@self.A@y
        self.BAy = B.T@self.A@y
        self.BAB = B.T@self.A@B

    def chisq(self,theta,model='cpl',amplitude=0.):
        g = mu(self.nodes,theta,model)-self.fidnodes
        amp = np.broadcast_to(np.asarray(amplitude),len(g))
        v = np.column_stack([g,amp])
        return self.yAy-2*v@self.BAy+np.sum((v@self.BAB)*v,axis=1)

    def chisq_full(self,theta,model='cpl',amplitude=0.):
        g = mu(self.z,theta,model,self.zhel)
        amp = np.broadcast_to(np.asarray(amplitude),len(g))
        r = self.mag-g-amp[:,None]*self.template
        return np.einsum('bi,ij,bj->b',r,self.A,r)

    def diagnostics(self):
        return dict(rows=len(self.z),unique_CID=int(self.df.CID.nunique()),
                    input_rows=int(np.loadtxt(self.inputs[1],max_rows=1)),
                    zHD_min=float(self.z.min()),zHD_max=float(self.z.max()),
                    covariance_min_eigenvalue=float(np.linalg.eigvalsh(self.cov)[0]),
                    input_covariance_max_asymmetry_mag2=self.input_asymmetry,
                    magnitude_revision_rows=self.magnitude_revision_rows,
                    offset_projection_error=float(np.max(np.abs(self.A@np.ones(len(self.z))))),
                    interpolation_nodes=len(self.nodes))


class BAO:
    def __init__(self):
        p=ROOT/'sources/repos/CobayaSampler__bao_data/desi_bao_dr2'
        self.inputs=[p/'desi_gaussian_bao_ALL_GCcomb_mean.txt',p/'desi_gaussian_bao_ALL_GCcomb_cov.txt']
        self.data=pd.read_csv(self.inputs[0],sep=r'\s+',comment='#',names=['z','value','kind'])
        self.z=self.data.z.to_numpy()
        self.inv=np.linalg.inv(np.loadtxt(self.inputs[1]))
        self.value=self.data.value.to_numpy()

    def prediction(self,theta,hrd):
        dm=299792.458/np.atleast_1d(hrd)[:,None]*integral(self.z,theta,'cpl')
        dh=299792.458/np.atleast_1d(hrd)[:,None]/efunc(self.z,theta,'cpl')
        dv=(self.z*dm**2*dh)**(1/3)
        return np.where(self.data.kind.to_numpy()=='DM_over_rs',dm,
                        np.where(self.data.kind.to_numpy()=='DH_over_rs',dh,dv))

    def chisq(self,theta,hrd):
        r=self.value-self.prediction(theta,hrd)
        return np.einsum('bi,ij,bj->b',r,self.inv,r)
