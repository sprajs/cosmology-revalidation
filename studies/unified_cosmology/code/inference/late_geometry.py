"""Independent late-time SN+BAO geometry with a free H0*rdrag scale.

No CMB information or early-time w0+wa cut. Matter plus CPL dark energy is
used over z<=2.33; radiation is omitted here, unlike the full CAMB analysis.
Chebyshev compression approximates only the smooth *theory* distance curve;
all 1820 objects and the full covariance enter its exact quadratic form.
"""
from pathlib import Path
import argparse
import json
import hashlib
import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.interpolate import CubicSpline
from scipy.integrate import quad
from scipy.optimize import differential_evolution, minimize
from numpy.polynomial import chebyshev as cheb
from numpy.polynomial.legendre import leggauss

ROOT = Path(__file__).resolve().parents[4]
WORK = ROOT/'.work/unified-cosmology/inference/late-geometry'
BAO = ROOT/'.work/unified-cosmology/external-probes/packages/data/bao_data/desi_bao_dr2'
SN = ROOT/'.work/unified-cosmology/survey-selection/normalized/dovekie-total.npz'


def sample_path(sample):
    if sample=='dovekie':return SN
    if sample=='pantheon':return ROOT/'.work/unified-cosmology/inference/pantheon/pantheon-total.npz'
    if sample=='des3yr':return SN.parent/'des3yr-combined-total.npz'
    raise ValueError(sample)


class Geometry:
    def __init__(self, model='cpl', evolution='none', sn=True, order=40, quadrature=128, sample='dovekie'):
        self.model,self.evolution,self.sn = model,evolution,sn
        self.names = ['Omega_m','H0_rdrag']+({'lcdm':[], 'wcdm':['w0'], 'cpl':['w0','wa']}[model])+(['epsilon'] if evolution=='linear' else [])
        self.bounds = [[.01,.99],[5000.,15000.]]+({'lcdm':[], 'wcdm':[[-3.,1.]],'cpl':[[-3.,1.],[-5.,3.]]}[model])+([[-.5,.5]] if evolution=='linear' else [])
        raw = np.loadtxt(BAO/'desi_gaussian_bao_ALL_GCcomb_mean.txt',dtype=str)
        self.bz,self.bmean,self.btype = raw[:,0].astype(float),raw[:,1].astype(float),raw[:,2]
        self.bprecision = np.linalg.inv(np.loadtxt(BAO/'desi_gaussian_bao_ALL_GCcomb_cov.txt'))
        sn_path = sample_path(sample)
        with np.load(sn_path) as f:
            self.sz,self.zhel,self.observed,self.cov = (f[k] for k in ['zHD','zHEL','MU','covariance'])
        # Bounded by observed max; z=0 is regular in log(DM/z).
        self.scale = np.log1p(max(self.sz))
        self.xnodes = cheb.chebpts2(order)
        self.nodes = np.expm1((self.xnodes+1)*self.scale/2)
        v = cheb.chebvander(2*np.log1p(self.sz)/self.scale-1,order-1)
        self.W = np.linalg.solve(cheb.chebvander(self.xnodes,order-1).T,v.T).T
        if evolution in ['smooth01','smooth03']:
            knots = [0.,.1,.4,.8,1.3]
            basis = CubicSpline(knots,np.eye(5)[:,1:],bc_type='natural')(self.sz)
            width = .1 if evolution=='smooth01' else .3
            self.cov = self.cov+width*width*basis@basis.T
        p = cho_solve(cho_factor(self.cov,lower=True),np.eye(len(self.sz)))
        u = p@np.ones(len(p)); self.p = p-np.outer(u,u)/u.sum()
        self.data = self.observed-5*np.log10(self.sz*(1+self.zhel))
        self.data -= self.data.mean()
        self.Q = self.W.T@self.p@self.W
        self.b = self.W.T@self.p@self.data
        self.dd = self.data@self.p@self.data
        gx,gw = leggauss(quadrature)
        self.gx,self.gw = (gx+1)/2,gw/2

    def physical(self, theta):
        theta = np.atleast_2d(theta)
        om,scale = theta[:,0],theta[:,1]
        w = np.full(len(theta),-1.) if self.model=='lcdm' else theta[:,2]
        wa = theta[:,3] if self.model=='cpl' else np.zeros(len(theta))
        eps = theta[:,-1] if self.evolution=='linear' else np.zeros(len(theta))
        return om,scale,w,wa,eps

    @staticmethod
    def expansion(z, om, w, wa):
        l = np.log1p(z)
        return np.sqrt(om[:,None,None]*np.exp(3*l)+(1-om[:,None,None])*np.exp(3*(1+w[:,None,None]+wa[:,None,None])*l-3*wa[:,None,None]*z/(1+z)))

    def integrals(self,z,om,w,wa):
        z = np.asarray(z)
        return np.sum(self.gw/self.expansion(z[None,:,None]*self.gx,om,w,wa),axis=-1)

    def components(self,theta):
        om,scale,w,wa,eps = self.physical(theta)
        ratio = self.integrals(self.nodes,om,w,wa)
        curve = 5*np.log10(ratio)+eps[:,None]*np.log1p(self.nodes)/np.log(2)
        schi = np.einsum('bi,ij,bj->b',curve,self.Q,curve)-2*curve@self.b+self.dd
        ii = self.bz*self.integrals(self.bz,om,w,wa)
        ee = self.expansion(self.bz[None,:,None],om,w,wa)[:,:,0]
        dh = 299792.458/scale[:,None]/ee
        dm = 299792.458/scale[:,None]*ii
        dv = np.cbrt(self.bz*dh*dm*dm)
        pred = np.where(self.btype=='DH_over_rs',dh,np.where(self.btype=='DM_over_rs',dm,dv))
        r = pred-self.bmean
        bchi = np.einsum('bi,ij,bj->b',r,self.bprecision,r)
        return schi,bchi

    def logp(self,theta):
        theta = np.atleast_2d(theta)
        low,high = np.asarray(self.bounds).T
        keep = ((theta>low)&(theta<high)).all(axis=1)
        values = np.full(len(theta),-np.inf)
        if keep.any():
            ss,bb = self.components(theta[keep])
            values[keep] = -.5*(bb+(ss if self.sn else 0))
        return values

    def direct_chi2(self,theta):
        om,_,w,wa,eps = [p[0] for p in self.physical(theta)]
        def inv_e(z):
            return 1/np.sqrt(om*(1+z)**3+(1-om)*np.exp(3*(1+w+wa)*np.log1p(z)-3*wa*z/(1+z)))
        ratios = np.array([quad(inv_e,0,z,epsabs=1e-11,epsrel=1e-11)[0]/z for z in self.sz])
        residual = 5*np.log10(ratios)+eps*np.log1p(self.sz)/np.log(2)-self.data
        return float(residual@self.p@residual)


def validate():
    rng = np.random.default_rng(92627)
    checks = []
    for evolution in ['none','linear','smooth01','smooth03']:
        g = Geometry(evolution=evolution)
        probes = [np.array([.31,9900,-.8,-.7]+([.08] if evolution=='linear' else []))]
        for _ in range(12):
            probes.append(np.array([rng.uniform(*bounds) for bounds in g.bounds]))
        errors = []
        for theta in probes:
            compressed = g.components(theta)[0][0]
            errors.append(abs(compressed-g.direct_chi2(theta)))
        assert max(errors)<1e-4,(evolution,max(errors))
        checks.append({'evolution':evolution,'prior_and_reference_probes':len(probes),'maximum_SN_chi2_absolute_error':max(errors)})
    out = ROOT/'studies/unified_cosmology/results/inference'
    out.mkdir(parents=True,exist_ok=True)
    (out/'late-geometry-validation.json').write_text(json.dumps({'checks':checks,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n')
    print(json.dumps(checks,indent=2))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--validate',action='store_true')
    p.add_argument('--model',choices=['lcdm','wcdm','cpl'],default='cpl')
    p.add_argument('--evolution',choices=['none','linear','smooth01','smooth03'],default='none')
    p.add_argument('--no-sn',action='store_true')
    p.add_argument('--sample',choices=['dovekie','pantheon','des3yr'],default='dovekie')
    p.add_argument('--steps',type=int,default=10000)
    p.add_argument('--ensembles',type=int,default=4)
    p.add_argument('--resume',action='store_true')
    p.add_argument('--move',choices=['stretch','mixed'],default='mixed')
    a = p.parse_args()
    if a.validate:
        validate();return
    import emcee
    g = Geometry(a.model,a.evolution,not a.no_sn,sample=a.sample)
    name = f'{a.model}-{a.evolution}-'+('nonsn' if a.no_sn else a.sample)
    out = WORK/name;out.mkdir(parents=True,exist_ok=True)
    if a.resume:
        record = json.loads((out/'fit.json').read_text())
        assert record['parameter_names']==g.names and record['bounds']==g.bounds
    else:
        if (out/'fit.json').exists() and any(out.glob('ensemble-*.h5')):
            raise RuntimeError('Existing run: use --resume to extend explicitly.')
        fit = differential_evolution(lambda t:-2*g.logp(t)[0],g.bounds,seed=272609,tol=1e-8,polish=True)
        record = {'parameter_names':g.names,'bounds':g.bounds,'optimum':fit.x.tolist(),'chi2':float(fit.fun),'optimizer_success':bool(fit.success),'components_SN_BAO':[float(v[0]) for v in g.components(fit.x)],'arguments':vars(a),'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
        (out/'fit.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record),flush=True)
    width = np.array([.005,20]+({'lcdm':[],'wcdm':[.02],'cpl':[.02,.04]}[a.model])+([.01] if a.evolution=='linear' else []))
    for ensemble in range(a.ensembles):
        backend = emcee.backends.HDFBackend(str(out/f'ensemble-{ensemble}.h5'))
        previous = backend.iteration if backend.initialized else 0
        if previous and not a.resume:
            raise RuntimeError('Existing run: refuse silent overwrite or repeated chains.')
        np.random.seed(272609+ensemble)
        start = None if previous else np.array(record['optimum'])+np.random.normal(size=(48,len(g.names)))*width
        moves = ([(emcee.moves.DEMove(),.7),(emcee.moves.DESnookerMove(),.1),
                  (emcee.moves.StretchMove(),.2)] if a.move=='mixed' else emcee.moves.StretchMove())
        sampler = emcee.EnsembleSampler(48,len(g.names),g.logp,vectorize=True,backend=backend,moves=moves)
        if a.steps>previous:
            sampler.run_mcmc(start,a.steps-previous,progress=False)
        print(json.dumps({'ensemble':ensemble,'steps':a.steps,'acceptance':float(sampler.acceptance_fraction.mean()),'autocorrelation_time':sampler.get_autocorr_time(tol=0).tolist()}),flush=True)
    with (out/'run-history.jsonl').open('a') as file:
        file.write(json.dumps({'arguments':vars(a),'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})+'\n')


if __name__=='__main__':
    main()
