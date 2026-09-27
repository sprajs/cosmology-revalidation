"""Source-bound M20 forward operator and two proper distance densities."""
from __future__ import annotations
import hashlib
import inspect
import json
import os
from pathlib import Path
import sys
import textwrap

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
ARCHIVE = Path.home()/'.local/share/cosmology-revalidation/archive/17487bf659fcbdeeea072221492bac14b04a0a85'
RUN = Path('runs/research_2026_09_26')
GRAY = .088


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def configure_cpu(cpu=None):
    available = sorted(os.sched_getaffinity(0))
    os.sched_setaffinity(0, [available[0] if cpu is None else cpu])
    for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
        os.environ[key] = '1'
    os.environ['JAX_PLATFORM_NAME'] = 'cpu'
    os.environ['JAX_ENABLE_X64'] = 'true'
    os.environ['XLA_FLAGS'] = '--xla_cpu_multi_thread_eigen=false --xla_force_host_platform_device_count=1'


def dump(path, value):
    path = Path(path)
    assert not path.exists(), f'Never overwrite an execution record: {path}'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def source_hashes():
    return {p.name: sha(p) for p in sorted(HERE.iterdir()) if p.suffix in ('.py','.json')}


def check_parent(archive):
    design = json.loads((HERE/'design.json').read_text())
    result = ROOT/'studies/infrared/results/bayesn-nir-forward-gate.json'
    assert sha(result) == design['parent_forward_sha256']
    gate = json.loads(result.read_text())
    assert gate['status'] == 'passed'
    assert sha(ROOT/'studies/infrared/code/bayesn_nir_forward_gate.py') == gate['source_sha256']
    for name, digest in gate['execution_design']['input_sha256'].items():
        assert sha(Path(archive)/name) == digest, name
    return design, gate


def filter_config(archive, directory):
    from ruamel.yaml import YAML
    archive = Path(archive)
    cfg = YAML(typ='safe').load((archive/RUN/'bayesn_signed_optical_pilot/release-filter-config.yaml').read_text())
    for group in ('standards','filters'):
        for item in cfg[group].values():
            item['path'] = str(archive/item['path'].split('/supernova/',1)[1])
    path = Path(directory)/'release-filter-config.yaml'
    if not path.exists():
        with path.open('w') as stream:
            YAML().dump(cfg,stream)
    else:
        assert YAML(typ='safe').load(path.read_text()) == cfg
    return path


class Kernel:
    """All rows share one object's 47 latent coordinates; no data flux required."""
    def __init__(self, archive, config, metadata, rows, bins=2400):
        import jax
        import jax.numpy as jnp
        sys.path.insert(0,str(Path(archive)/RUN/'bayesn_distance_identification/official-code'))
        import bayesn.bayesn_model as module
        from bayesn import SEDmodel
        source = textwrap.dedent(inspect.getsource(SEDmodel._setup_band_weights))
        assert source.count('self.spectrum_bins = 300') == 1
        source = source.replace('self.spectrum_bins = 300',f'self.spectrum_bins = {int(bins)}')
        namespace = {}
        exec(compile(source,'<declared-resolution>', 'exec'),module.__dict__,namespace)
        class Refined(SEDmodel):
            _setup_band_weights = namespace['_setup_band_weights']
        self.model = m = Refined(load_model='M20_model',num_devices=1,filter_yaml=str(config))
        times = jnp.array([r['trigger_rest_time'] for r in rows])[:,None]
        assert bool(jnp.all((times>10)&(times<30)))
        nobs = len(rows)
        band = jnp.array([m.band_dict['RAISIN_DES_'+r['band']] for r in rows])[:,None]
        weight = m._calculate_band_weights(jnp.array([metadata['zHEL']]),jnp.array([metadata['MWEBV']]))
        def forward(p):
            p = p[None,:]
            phase = times-p[:,46][None,:]
            jt = m.J_t_map(phase.flatten(order='F'),m.tau_knots,m.KD_t).reshape((nobs,1,6),order='F').transpose(1,2,0)
            hs = jnp.array([19+jnp.floor(phase),19+jnp.ceil(phase),jnp.remainder(phase,1)])
            middle = (m.L_Sigma@p[:,4:46].T).T.reshape((1,7,6),order='F')
            eps = jnp.zeros((1,9,6)).at[:,1:-1,:].set(middle)
            return m.get_flux_batch(m.M0,p[:,3],p[:,1],m.W0,m.W1,eps,p[:,0],p[:,2],band,jnp.ones((nobs,1)),jt,hs,weight)[:,0]
        self.forward = jax.jit(forward)
        self.batch = jax.jit(jax.vmap(forward))


def broad_logpdf(d, lo=20., hi=50.):
    import jax.numpy as jnp
    from jax.scipy.special import log_ndtr
    a = jnp.where(d<(lo+hi)/2,(d-lo)/GRAY,(hi-d)/GRAY)
    b = jnp.where(d<(lo+hi)/2,(d-hi)/GRAY,(lo-d)/GRAY)
    la,lb = log_ndtr(a),log_ndtr(b)
    return la+jnp.log(-jnp.expm1(lb-la))-jnp.log(hi-lo)


def numpyro_model(kernel, data, arm):
    import jax
    import jax.numpy as jnp
    import numpyro
    import numpyro.distributions as dist
    from numpyro.distributions import constraints
    class BroadD(dist.Distribution):
        support = constraints.real
        arg_constraints = {}
        def sample(self,key,sample_shape=()):
            a,b = jax.random.split(key)
            return jax.random.uniform(a,sample_shape,minval=20.,maxval=50.)+GRAY*jax.random.normal(b,sample_shape)
        def log_prob(self,value):
            return broad_logpdf(value)
    def likelihood():
        meta = data['metadata']
        D = numpyro.sample('D',dist.Normal(meta['mu_LCDM'],(meta['sigma_external']**2+GRAY**2)**.5) if arm=='LCDM' else BroadD())
        AV = numpyro.sample('AV',dist.Exponential(1/.329))
        RV = numpyro.sample('RV',dist.Uniform(1.2,6.))
        theta = numpyro.sample('theta',dist.Normal(0.,1.))
        eps = numpyro.sample('epsilon_white',dist.Normal(jnp.zeros(42),jnp.ones(42)).to_event(1))
        tau = numpyro.sample('tau',dist.Uniform(-10.,20.))
        p = jnp.concatenate([jnp.array([D,AV,RV,theta]),eps,jnp.array([tau])])
        numpyro.sample('optical',dist.Normal(kernel.forward(p),jnp.array(data['errors'])).to_event(1),obs=jnp.array(data['flux']))
    return likelihood
