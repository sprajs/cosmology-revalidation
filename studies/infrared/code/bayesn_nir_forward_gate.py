"""Metadata-only continuation of the frozen signed-optical BayeSN gate.

No measured NIR flux is converted to a number, no sealed payload is opened,
and no likelihood, posterior or predictive outcome score is evaluated.
The archived author source and science inputs are read-only. This caller
remaps historical paths and changes only the integration resolution.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import inspect
import json
import os
from pathlib import Path
import sys
import textwrap
import time

ROOT = Path(__file__).resolve().parents[3]
ARCHIVE = Path('/home/szymon/.local/share/cosmology-revalidation/archive/17487bf659fcbdeeea072221492bac14b04a0a85')
RUN = Path('runs/research_2026_09_26')
CIDS = ('DES16E2clk', 'DES16X3cry')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, default=ARCHIVE)
    parser.add_argument('--work', type=Path, default=ROOT/'.work/infrared/nir-forward-gate')
    parser.add_argument('--result', type=Path, default=ROOT/'studies/infrared/results/bayesn-nir-forward-gate.json')
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    assert not args.result.exists(), 'Preserve every previous result; choose a new --result.'
    os.sched_setaffinity(0, [min(os.sched_getaffinity(0))])
    for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        os.environ[key] = '1'
    os.environ['JAX_PLATFORM_NAME'] = 'cpu'
    os.environ['JAX_ENABLE_X64'] = 'true'
    os.environ['XLA_FLAGS'] = '--xla_cpu_multi_thread_eigen=false --xla_force_host_platform_device_count=1'
    start = time.monotonic()
    import numpy as np
    from ruamel.yaml import YAML
    import jax
    import jax.numpy as jnp

    pilot = args.archive/RUN/'bayesn_signed_optical_pilot'
    old = args.archive/RUN/'bayesn_distance_identification'
    design_dir = args.archive/RUN/'astra_design/bayesn_signed_pilot'
    adapter = json.loads((pilot/'adapter-manifest.json').read_text())
    inputs = {}
    def pin(path, expected=None):
        digest = sha(path)
        if expected is not None:
            assert digest == expected, path
        inputs[str(path.relative_to(args.archive))] = digest
        return digest
    pin(pilot/'adapter-manifest.json')
    for name in ('frozen-cohort.csv', 'frozen-mask-metadata.csv', 'pilot-design-protocol.json'):
        pin(design_dir/name, adapter['inputs'][str((RUN/'astra_design/bayesn_signed_pilot')/name)])
    pin(pilot/'release-filter-config.yaml', adapter['outputs']['release-filter-config.yaml'])
    pin(pilot/'RAISIN_DES_J.dat', adapter['outputs']['RAISIN_DES_J.dat'])
    for path in sorted((old/'official-code/bayesn').rglob('*')):
        if path.is_file() and '__pycache__' not in path.parts and path.suffix != '.pyc':
            pin(path)
    assert inputs[str(RUN/'bayesn_distance_identification/official-code/bayesn/bayesn_model.py')] == adapter['inputs'][str(RUN/'bayesn_distance_identification/official-code/bayesn/bayesn_model.py')]
    cohort = {r['CID']: r for r in csv.DictReader((design_dir/'frozen-cohort.csv').open())}
    metadata = [r for r in csv.DictReader((design_dir/'frozen-mask-metadata.csv').open())
                if r['CID'] in CIDS and r['role'] == 'heldout_NIR' and r['primary_keep'] == 'True']
    rows = []
    # Only the quoted-error coordinate is added to frozen epoch/filter metadata.
    # The FLUXCAL token at index 4 is never converted, stored or inspected.
    for r in metadata:
        path = args.archive/r['source_path']
        pin(path, adapter['inputs'][r['source_path']])
        with path.open() as stream:
            for number, line in enumerate(stream, 1):
                if number == int(r['source_line']):
                    tokens = line.split()
                    assert tokens[0] == 'OBS:' and float(tokens[1]) == float(r['MJD']) and tokens[2] == r['band']
                    error = float(tokens[5])
                    assert error > 0 and np.isfinite(error)
                    break
            else:
                raise ValueError('Missing frozen metadata row')
        rows.append({k: r[k] for k in ('CID', 'source_path', 'source_line', 'MJD', 'band', 'trigger_rest_time')} | {'error': error})
    assert [sum(r['CID'] == cid for r in rows) for cid in CIDS] == [6, 5]
    # An executable input-read guard excludes the sealed outcome file.
    def no_sealed_payload(event, arguments):
        if event == 'open' and 'nir-sealed' in str(arguments[0]):
            raise RuntimeError('Sealed NIR outcome access forbidden')
    sys.addaudithook(no_sealed_payload)
    yaml = YAML(typ='safe')
    config = yaml.load((pilot/'release-filter-config.yaml').read_text())
    remaps = {}
    for group in ('standards', 'filters'):
        for entry in config[group].values():
            original = entry['path']
            suffix = original.split('/supernova/', 1)[1]
            path = args.archive/suffix
            pin(path)
            entry['path'] = str(path)
            remaps[original] = str(path)
    config_path = args.work/'release-filter-config.yaml'
    with config_path.open('w') as stream:
        YAML().dump(config, stream)
    plan = dict(scope='Metadata-only J/H forward validation; no posterior or observed-flux score.',
                source_sha256=sha(__file__), archive=str(args.archive), input_sha256=inputs,
                objects=list(CIDS), rows=rows, integration_bins=[1200, 2400, 4800],
                test_parameters='D=mu_LCDM; fiducial AV=.3,RV=3.1,theta=0,epsilon=0,tau=5; tau endpoints -10,20; AV=.6,RV=2,theta=.5,epsilon seed202; zero dust.',
                gates=dict(final_resolution_quoted_sigma=.1, static_dynamic_quoted_sigma=1e-8,
                           eager_jit_quoted_sigma=1e-8, gray_relative=1e-12, AV0_RV_quoted_sigma=1e-12,
                           knot_symmetric_step_quoted_sigma=1e-5), wall_cap_seconds=300,
                interpretation='One worker, no sampling. Finite knot tests do not assert differentiability at Hsiao integer-phase kinks.',
                historical_path_remapping=remaps)
    plan_path = args.work/'execution-design.json'
    assert not plan_path.exists(), 'Use a fresh work directory for every execution.'
    plan_path.write_text(json.dumps(plan, indent=2)+'\n')
    sys.path.insert(0, str(old/'official-code'))
    import bayesn.bayesn_model as module
    from bayesn import SEDmodel
    source = textwrap.dedent(inspect.getsource(SEDmodel._setup_band_weights))
    assert source.count('self.spectrum_bins = 300') == 1
    source = source.replace('self.spectrum_bins = 300', 'self.spectrum_bins = self.audit_spectrum_bins')
    namespace = {}
    exec(compile(source, '<resolution-only-override>', 'exec'), module.__dict__, namespace)
    class Refined(SEDmodel):
        _setup_band_weights = namespace['_setup_band_weights']
    n, nobs = 2, 6
    times = np.full((nobs, n), 15.)
    errors = np.ones((nobs, n))
    mask = np.zeros((nobs, n))
    bands = np.full((nobs, n), 'NULL_BAND', object)
    for j, cid in enumerate(CIDS):
        for i, row in enumerate(r for r in rows if r['CID'] == cid):
            times[i,j], errors[i,j], mask[i,j] = float(row['trigger_rest_time']), row['error'], 1
            bands[i,j] = 'RAISIN_DES_'+row['band']
    assert np.all((times[mask > 0] > 10) & (times[mask > 0] < 30))
    fid = np.zeros((n, 47))
    fid[:,0] = [float(cohort[c]['mu_LCDM']) for c in CIDS]
    fid[:,1], fid[:,2], fid[:,46] = .3, 3.1, 5.
    cases = {'fiducial': fid}
    for label, tau in [('tau_lower', -10.), ('tau_upper', 20.)]:
        q = fid.copy(); q[:,46] = tau; cases[label] = q
    q = fid.copy(); q[:,1:4] = [.6, 2., .5]; q[:,4:46] = np.random.default_rng(202).normal(size=(n,42)); cases['nonzero_SED'] = q
    q = fid.copy(); q[:,1] = 0.; cases['zero_dust'] = q
    z = np.array([float(cohort[c]['zHEL']) for c in CIDS])
    mw = np.array([float(cohort[c]['MWEBV']) for c in CIDS])
    records, previous = [], None
    for bins in plan['integration_bins']:
        Refined.audit_spectrum_bins = bins
        model = Refined(load_model='M20_model', num_devices=1, filter_yaml=str(config_path))
        weights = model._calculate_band_weights(jnp.array(z), jnp.array(mw))
        indices = jnp.array([[model.band_dict[b] for b in row] for row in bands])
        def epsilon(p):
            middle = (model.L_Sigma@p[:,4:46].T).T.reshape((n,7,6), order='F')
            return jnp.zeros((n,9,6)).at[:,1:-1,:].set(middle)
        def forward(p):
            phase = jnp.array(times)-p[:,46][None,:]
            jt = model.J_t_map(phase.flatten(order='F'), model.tau_knots, model.KD_t).reshape((nobs,n,6), order='F').transpose(1,2,0)
            hs = jnp.array([19+jnp.floor(phase),19+jnp.ceil(phase),jnp.remainder(phase,1)])
            return model.get_flux_batch(model.M0,p[:,3],p[:,1],model.W0,model.W1,epsilon(p),p[:,0],p[:,2],indices,jnp.array(mask),jt,hs,weights)
        compiled = jax.jit(forward)
        predictions, checks = {}, []
        for name, p in cases.items():
            prediction = np.asarray(compiled(jnp.array(p)))
            eager = np.asarray(forward(jnp.array(p)))
            assert np.isfinite(prediction).all()
            static = np.zeros_like(prediction)
            eps = np.asarray(epsilon(jnp.array(p)))
            for j in range(n):
                m = mask[:,j] > 0
                # Independent public API builds phases and row layout itself.
                got, _, _ = model.simulate_light_curve(times[m,j], 1, bands[m,j].tolist(),
                    yerr=0.,err_type='flux',z=z[j],mu=p[j,0],ebv_mw=mw[j],RV=p[j,2],
                    tmax=p[j,46],del_M=0.,AV=p[j,1],theta=p[j,3],eps=eps[j:j+1],mag=False)
                static[m,j] = np.asarray(got)[:,0]
            checks.append(dict(case=name, eager_jit_sigma=float(np.max(abs(prediction-eager)/errors*mask)),
                               static_dynamic_sigma=float(np.max(abs(prediction-static)/errors*mask)),
                               resolution_delta_sigma=None if previous is None else float(np.max(abs(prediction-previous[name])/errors*mask))))
            predictions[name] = prediction
        records.append(dict(bins=bins, checks=checks, elapsed_seconds=time.monotonic()-start))
        previous = predictions
        print(json.dumps(records[-1]), flush=True)
        if time.monotonic()-start > 300:
            raise TimeoutError('Prespecified 300-second gate cap')
    zero = cases['zero_dust'].copy(); zero[:,2] = 1.2
    f1 = np.asarray(compiled(jnp.array(zero))); zero[:,2] = 6.
    av0 = float(np.max(abs(f1-np.asarray(compiled(jnp.array(zero))))/errors*mask))
    q = fid.copy(); q[:,0] += .15
    expected = previous['fiducial']*10**(-.4*.15)
    gray = float(np.max(abs(np.asarray(compiled(jnp.array(q)))-expected)/np.maximum(abs(expected),1e-200)*mask))
    crossings = []
    for j,cid in enumerate(CIDS):
        for knot in np.asarray(model.tau_knots):
            tau = times[0,j]-knot
            if not -10 < tau < 20:
                continue
            q = fid.copy(); q[j,46] = tau
            center = np.asarray(compiled(jnp.array(q)))
            q[j,46] = tau-1e-6; left = np.asarray(compiled(jnp.array(q)))
            q[j,46] = tau+1e-6; right = np.asarray(compiled(jnp.array(q)))
            crossings.append(dict(cid=cid,tau=float(tau),finite=bool(np.isfinite(center+left+right).all()),
                                  symmetric_step_sigma=float(np.max(abs(left+right-2*center)/errors*mask))))
    gates = plan['gates']
    passed = all(c['eager_jit_sigma'] < gates['eager_jit_quoted_sigma'] and c['static_dynamic_sigma'] < gates['static_dynamic_quoted_sigma'] for r in records for c in r['checks'])
    passed &= all(c['resolution_delta_sigma'] < gates['final_resolution_quoted_sigma'] for c in records[-1]['checks'])
    passed &= gray < gates['gray_relative'] and av0 < gates['AV0_RV_quoted_sigma']
    passed &= all(c['finite'] and c['symmetric_step_sigma'] < gates['knot_symmetric_step_quoted_sigma'] for c in crossings)
    for path,digest in inputs.items():
        assert sha(args.archive/path) == digest, 'Source drift during gate'
    assert sha(__file__) == plan['source_sha256']
    result = dict(status='passed' if passed else 'failed', scope=plan['scope'],
                  source_sha256=sha(__file__), execution_design=plan, execution_design_sha256=sha(plan_path),
                  source_archive_commit='17487bf659fcbdeeea072221492bac14b04a0a85',
                  numerical_checks=records, AV0_RV_sigma=av0, gray_relative=gray, knot_crossings=crossings,
                  versions=dict(numpy=np.__version__,jax=jax.__version__,python=sys.version),
                  elapsed_seconds=time.monotonic()-start, CPU_affinity=sorted(os.sched_getaffinity(0)),
                  observed_NIR_flux_parsed=False, sealed_NIR_payload_opened=False, posterior_executed=False,
                  static_comparison='Official simulate_light_curve; explicit identical latents, zero random noise; same author spectral kernel.',
                  remaining_gates=['Full synthetic recovery and mixing','Observed optical posterior with frozen diagnostics','Held-out joint NIR predictive score and MC error'])
    args.result.parent.mkdir(parents=True, exist_ok=True)
    args.result.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k: result[k] for k in ('status','elapsed_seconds','AV0_RV_sigma','gray_relative')}),flush=True)
    assert passed


if __name__ == '__main__':
    main()
