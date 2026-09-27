"""Synthetic anchored-replacement mathematics and fail-closed consumer tests."""
import argparse
import ast
import copy
from contextlib import ExitStack
import importlib.metadata
import json
from pathlib import Path
from types import SimpleNamespace
import tempfile
from unittest.mock import patch

import numpy as np
from scipy.special import logsumexp
from scipy.stats import norm

import anchored_bridge as bridge


def rejected(call, contains=None):
    try:
        call()
    except (AssertionError, ValueError) as error:
        if contains:
            assert contains in str(error), str(error)
        return
    raise AssertionError('Invalid fixture accepted.')


def _validate():
    design = json.loads(bridge.DESIGN.read_text())
    gates = json.loads(bridge.GATES.read_text())['overlap_gates']
    model_path = bridge.ROOT/'.work/unified-cosmology/inference/quartic-spectral/surrogate-quartic-v1.npz'
    configurations = []
    for model in ['lcdm', 'cpl']:
        settings = dict(model=model, evolution='none', sample='dovekie', calibration='official_planck',
                        fast_lensing=True, gpu=False, surrogate=str(model_path))
        options = {k: settings[k] for k in ['model', 'evolution', 'sample', 'calibration']}
        frozen = {'configuration': bridge.canonical(bridge.parent_configuration(**options, surrogate=model_path))}
        old, new = bridge.configurations(settings, frozen, design)
        left, right = map(copy.deepcopy, [bridge.canonical(old), bridge.canonical(new)])
        left['likelihood'].pop('released_sn'); right['likelihood'].pop('released_sn')
        assert left == right
        configurations.append({'model': model, 'non_SN_configuration_identical': True})
        for key, value in [('model','other'), ('evolution','linear'), ('sample','pantheon'),
                           ('calibration','independent'), ('fast_lensing',False), ('gpu',True)]:
            rejected(lambda key=key, value=value: bridge.configurations(dict(settings, **{key:value}), frozen, design))
        bad = copy.deepcopy(frozen); bad['configuration']['params']['H0']['prior']['max'] += 1.
        rejected(lambda: bridge.configurations(settings, bad, design), 'configuration changed')
        changed = copy.deepcopy(new); changed['params']['H0']['prior']['max'] += 1.
        with patch.object(bridge.anchored, 'configuration', return_value=changed):
            rejected(lambda: bridge.configurations(settings, frozen, design), 'non-SN')

    rng = np.random.default_rng(design['synthetic_seed'])
    n = design['synthetic_points']; groups = np.repeat(np.arange(4), n//4)
    fixed = (-.1, 1.7); original = (.2, 1.3); replacement = (.45, 1.1)
    precision_old = 1+1/fixed[1]**2+1/original[1]**2
    mean_old = (fixed[0]/fixed[1]**2+original[0]/original[1]**2)/precision_old
    proposal_mean, proposal_sd = mean_old+.02, precision_old**-.5*1.03
    x = rng.normal(proposal_mean, proposal_sd, n); aux = rng.normal(size=n)
    prior = norm.logpdf(x)+norm.logpdf(aux)
    cmb = norm.logpdf(x, *fixed); source = norm.logpdf(x, *original)
    target = norm.logpdf(x, *replacement)
    proposal = norm.logpdf(x, proposal_mean, proposal_sd)+norm.logpdf(aux)
    exact = prior+cmb+source
    records, backgrounds, rows = [], [], []
    for i in range(n):
        record = {'status':'finite', 'point':{'w':float(x[i]), 'A_fg':float(aux[i])},
            'derived':{'q0':float(x[i]-.25), 'q05':float(x[i]*.8), 'q1':float(x[i]*.5),
                       'j0':1., 'sn_chi2':float(((x[i]-original[0])/original[1])**2)},
            'exact_loglikes':{'released_sn':float(source[i]), 'synthetic_CMB_BAO':float(cmb[i])},
            'proposal_loglikes':{'released_sn':float(source[i]+proposal[i]-exact[i]), 'synthetic_CMB_BAO':float(cmb[i])},
            'exact_logpost':float(exact[i]), 'proposal_logpost':float(proposal[i]),
            'log_weight':float(exact[i]-proposal[i])}
        background = {'old_SN_reconstructed_loglike':float(source[i]),
                      'anchored_SN':{'loglike':float(target[i]), 'chi2':float(((x[i]-replacement[0])/replacement[1])**2)},
                      'background_calls':1}
        records.append(record); backgrounds.append(background)
        rows.append(bridge.replace_SN(record, background, record['exact_loglikes'], design))
    output = bridge.summarize(records, rows, groups, gates)
    assert output['qualified_under_declared_numerical_gates'], output['failed_gates']
    precision = 1+1/fixed[1]**2+1/replacement[1]**2
    true_mean = (fixed[0]/fixed[1]**2+replacement[0]/replacement[1]**2)/precision
    true_sd = precision**-.5
    posterior = output['posterior']['w']
    assert abs(posterior['mean']-true_mean) < design['synthetic_mean_absolute_tolerance']
    assert abs(posterior['sd']-true_sd) < design['synthetic_sd_absolute_tolerance']
    independent_lw = prior+cmb+target-proposal
    error = float(np.max(abs(independent_lw-[r['target_logweight'] for r in rows])))
    assert error < 1e-12
    assert 'A_fg' not in output['posterior'] and 'sn_chi2' not in output['posterior']
    assert 'A_fg' in output['weight_diagnostics']['weighted_chain_stability']
    weights = np.exp(independent_lw-logsumexp(independent_lw))
    target_chi2 = ((x-replacement[0])/replacement[1])**2
    assert abs(output['weight_diagnostics']['weighted_summaries_for_diagnostics']['sn_chi2']['mean']-weights@target_chi2) < 1e-12

    null_error = 0.
    shifted = []
    for record, background in zip(records[:256], backgrounds[:256]):
        null = copy.deepcopy(background)
        null['anchored_SN']['loglike'] = null['old_SN_reconstructed_loglike']
        row = bridge.replace_SN(record, null, record['exact_loglikes'], design)
        null_error = max(null_error, abs(row['target_logweight']-record['log_weight']))
        altered = copy.deepcopy(background); altered['anchored_SN']['loglike'] += 700.
        shifted.append(bridge.replace_SN(record, altered, record['exact_loglikes'], design)['target_logweight'])
    assert null_error == 0.
    first = independent_lw[:256]; shifted = np.asarray(shifted)
    constant_error = float(np.max(abs(np.exp(first-logsumexp(first))-np.exp(shifted-logsumexp(shifted)))))
    assert constant_error < 1e-13
    altered = copy.deepcopy(records[0]); altered['proposal_logpost'] += .01; altered['log_weight'] -= .01
    rejected(lambda: bridge.replace_SN(altered, backgrounds[0], altered['exact_loglikes'], design), 'priors differ')
    mismatch = copy.deepcopy(backgrounds[0]); mismatch['old_SN_reconstructed_loglike'] += .001
    bad_row = bridge.replace_SN(records[0], mismatch, records[0]['exact_loglikes'], design)
    assert bad_row['status'] == 'source_SN_density_mismatch'
    assert bad_row['target_logweight'] == rows[0]['target_logweight'], 'Reconstructed old density must not enter the ratio.'
    failed = bridge.summarize(records, [bad_row]+rows[1:], groups, gates)
    assert failed['posterior'] is None and failed['status']=='failed_anchored_bridge_evaluation'
    collapsed = copy.deepcopy(rows); collapsed[0]['target_logweight'] += 1000.
    failed = bridge.summarize(records, collapsed, groups, gates)
    assert 'raw_weight_ESS' in failed['failed_gates']
    assert all(failed[k] is None for k in ['posterior','conditional_sign_fractions','weighted_covariance','paired_mean_changes_from_parent'])
    json.dumps(failed, allow_nan=False)
    selected_indices = np.concatenate([np.arange(200)+i*(n//4) for i in range(4)])
    short = bridge.summarize([records[i] for i in selected_indices], [rows[i] for i in selected_indices], np.repeat(np.arange(4),200), gates)
    assert 'minimum_exact_points' in short['failed_gates'] and short['posterior'] is None

    with tempfile.TemporaryDirectory(dir=bridge.ROOT/'.work') as directory:
        directory = Path(directory); cache = directory/'refused'
        with patch.object(bridge, 'summarize_run', side_effect=ValueError('unqualified-parent')):
            with patch.object(bridge, 'configurations', side_effect=AssertionError('too-early')):
                rejected(lambda: bridge.actual(directory, directory/'absent', cache), 'unqualified-parent')
        assert not cache.exists()

        # Explicitly synthetic consumer fixture. Qualification and target/data/
        # background interfaces are mocked; file seals and cache paths are real.
        folder = directory/'chain'; native_dir = folder/'exact'; native_dir.mkdir(parents=True)
        data = directory/'synthetic-data'; data.write_text('not observational data')
        selection_path = native_dir/'selection.json'; summary_path = directory/'summary.json'
        picked = np.concatenate([np.arange(100)+i*(n//4) for i in range(4)])
        selected = [records[i] for i in picked]
        chosen_backgrounds = {r['point']['w']: backgrounds[i] for r,i in zip(selected,picked)}
        settings = {'model':'lcdm','evolution':'none','sample':'dovekie','calibration':'official_planck'}
        selection_path.write_text(json.dumps({'settings':settings,'points':[r['point'] for r in selected], 'groups':np.repeat(np.arange(4),100).tolist()}))
        summary_path.write_text(json.dumps({'selection_path':bridge.relative(selection_path)}))
        frozen = {'identity':'synthetic-only', 'source_sha256':{}, 'versions':{'numpy':importlib.metadata.version('numpy')},
                  'assets':{'synthetic':'no-data'}, 'sample_sha256':bridge.digest(data)}
        manifest = folder/'run-0.json'; manifest.write_text(json.dumps({'target_identity':frozen}))
        paths = [manifest,selection_path,summary_path]
        for index,record in enumerate(selected):
            path=native_dir/f'{index:05d}.json'; path.write_text(json.dumps(record)); paths.append(path)
        parent = {'input_sha256':{bridge.relative(p):bridge.digest(p) for p in paths},'settings':settings}
        config = {'likelihood':{'released_sn':{'data_file':str(data)},'synthetic_CMB_BAO':{}},'theory':{'camb':{'extra_args':{}}}}
        target_identity = {'identity':'synthetic-target', 'assets':frozen['assets'], 'versions':frozen['versions'],
                           'source_sha256':{}, 'anchored_calibration':{'input_and_audit_sha256':{}}}
        with patch.object(bridge,'summarize_run',return_value=parent), patch.object(bridge,'configurations',return_value=(config,config)), \
             patch.object(bridge.anchored,'identify',return_value=target_identity), patch.object(bridge,'load_SN_interfaces',return_value=(None,None)), \
             patch.object(bridge,'background_densities',side_effect=lambda p,*a:chosen_backgrounds[p['w']]):
            cache=directory/'cache'
            first=bridge.actual(folder,summary_path,cache)
            assert first['posterior'] is None and 'minimum_exact_points' in first['failed_gates']
            with patch.object(bridge,'background_densities',side_effect=AssertionError('cache-replay-called-background')):
                second=bridge.actual(folder,summary_path,cache)
            assert first==second
            row_path=cache/'00000.json'; raw=row_path.read_bytes(); changed=json.loads(raw); changed['target_logweight']+=.1
            row_path.write_text(json.dumps(changed))
            rejected(lambda:bridge.actual(folder,summary_path,cache),'identity changed');row_path.write_bytes(raw)
            ledger=cache/'record-hashes.json'; ledger_raw=ledger.read_bytes();ledger.unlink()
            row_path.write_text(json.dumps(changed))
            rejected(lambda:bridge.actual(folder,summary_path,cache),'payload changed');row_path.write_bytes(raw);ledger.write_bytes(ledger_raw)
            old=manifest.read_bytes();manifest.write_bytes(old+b' ')
            rejected(lambda:bridge.actual(folder,summary_path,directory/'bad-parent'),'identity changed');manifest.write_bytes(old)
            old=data.read_bytes();data.write_bytes(b'altered')
            rejected(lambda:bridge.actual(folder,summary_path,directory/'bad-data'),'Parent SN data changed');data.write_bytes(old)
            assert not (directory/'bad-parent').exists() and not (directory/'bad-data').exists()

    # Exercise the production distance operator with a fake background. No CAMB
    # solver executes. Its distinct zHD/zHEL, repeats and reversed order make an
    # accidental heliocentric/cosmological redshift or ordering swap detectable.
    import camb
    zold=np.array([.32,.11,.32]); zhel=np.array([.321,.109,.322])
    znew=np.array([.27,.07,.27,.41])
    observed=np.array([40.,38.,40.01]); cov=np.array([[.04,.006,.002],[.006,.06,.004],[.002,.004,.05]])
    old=bridge.IntegratedLuminosity(zold,zhel,observed,cov)
    requested=[]; received=[]; kwargs_seen=[]; native_guard=[]
    def distance(z):return 4317*np.asarray(z)/(1+np.asarray(z))
    def get_distance(z):requested.append(np.array(z));return distance(z)
    def evaluate_new(da):
        received.append(np.array(da));return {'loglike':-float(np.sum(da/1000)**2), 'chi2':float(np.sum((da/1000)**2))}
    def background_only(parameters):
        assert parameters.WantTransfer is True
        for forbidden in [camb.get_results,camb.get_transfer_functions]:
            try:forbidden(parameters)
            except RuntimeError:native_guard.append(True)
            else:raise AssertionError('Native solver was not guarded.')
        return SimpleNamespace(angular_diameter_distance=get_distance)
    def parameters(**kwargs):kwargs_seen.append(kwargs);return SimpleNamespace(WantTransfer=False)
    fake_point=dict(H0=69.,ombh2=.0224,omch2=.12,ns=.97,tau=.055,logA=3.04,w=-.92,wa=.17)
    with patch.object(camb,'set_params',side_effect=parameters), patch.object(camb,'get_background',side_effect=background_only), \
         patch.object(bridge,'native_thermal_parameters',return_value=SimpleNamespace(WantTransfer=True)):
        synthetic_background=bridge.background_densities(fake_point,{'AccuracyBoost':1.},old,
            SimpleNamespace(z_hd_noncalibrator=znew,evaluate=evaluate_new))
    assert len(requested)==len(received)==len(kwargs_seen)==1 and len(native_guard)==2
    assert np.array_equal(requested[0],np.unique(np.r_[zold,znew]))
    assert np.array_equal(received[0],distance(znew))
    prediction=5*np.log10(distance(zold)*(1+zold)*(1+zhel))+25
    residual=prediction-observed; precision=np.linalg.inv(cov);one=np.ones(3)
    q=residual@precision@residual-(residual@precision@one)**2/(one@precision@one)
    direct_ll=-.5*(q+np.linalg.slogdet(cov)[1]+np.log(one@precision@one)+2*np.log(2*np.pi))
    background_error=abs(direct_ll-synthetic_background['old_SN_reconstructed_loglike'])
    assert background_error<1e-11
    assert abs(kwargs_seen[0]['As']-1e-10*np.exp(fake_point['logA']))<1e-25
    assert kwargs_seen[0]['w']==fake_point['w'] and kwargs_seen[0]['wa']==fake_point['wa']
    assert synthetic_background['distance_unit']=='Mpc' and synthetic_background['H0_unit']=='km/s/Mpc'

    sources = [Path(__file__), Path(bridge.__file__), bridge.DESIGN, bridge.GATES,
               bridge.HERE/'anchored_adapter.py', bridge.HERE/'expansion_history.py',
               bridge.HERE/'measurement_summary.py', bridge.HERE/'luminosity_sensitivity.py',
               bridge.HERE/'probe_omission.py', bridge.HERE/'modern_fast.py', bridge.HERE/'modern_run.py',
               bridge.HERE/'target_identity.py', bridge.HERE.parent/'distance_ladder/calibration_interface.py']
    for path in sources:
        if path.suffix=='.py':ast.parse(path.read_text())
    return {'status':'passed_synthetic_anchored_bridge_checks_no_observational_evaluation',
        'configuration_cases':configurations, 'synthetic_points':n,
        'known_Gaussian_recovery':{'expected_mean':true_mean,'mean':posterior['mean'],'expected_sd':true_sd,'sd':posterior['sd'],
                                   'maximum_density_error':error,'raw_weight_ESS':output['weight_diagnostics']['raw_weight_ESS']},
        'null_replacement_logweight_error':null_error,'additive_constant_normalized_weight_error':constant_error,
        'target_SN_chi2_replaces_old_diagnostic':True,'source_reconstruction_not_substituted_in_weight':True,
        'all_failed_gates_withhold_posterior':True,'unsupported_settings_and_altered_nonSN_priors_rejected':True,
        'consumer_fixture_points':400,'consumer_interfaces_and_qualification_explicitly_mocked':True,
        'cache_replay_identical_no_new_backgrounds':True,'cache_payload_and_manifest_and_parent_and_data_tampering_rejected':True,
        'synthetic_background_ordering_units_and_heliocentric_convention_passed':True,
        'synthetic_old_SN_independent_precision_loglike_error':background_error,
        'native_spectrum_and_transfer_calls_guarded':True,
        'unqualified_parent_refused_before_work':True,
        'observational_bridge_points':0,'CAMB_background_calls':0,'CMB_spectrum_calls':0,
        'source_sha256':{bridge.relative(p):bridge.digest(p) for p in sources},
        'limitations':'Synthetic algebra and consumer integrity only. Real every-point native SN/background closure and importance overlap are still required; no anchored cosmological measurement follows.'}


def validate():
    import camb
    import cobaya.model
    attempted=[]
    def forbidden(*args,**kwargs):
        attempted.append(True)
        raise AssertionError('Validation must not construct a model or execute native CAMB.')
    with ExitStack() as stack:
        for name in ['set_params','get_background','get_results','get_transfer_functions']:
            stack.enter_context(patch.object(camb,name,forbidden))
        stack.enter_context(patch.object(cobaya.model,'get_model',forbidden))
        stack.enter_context(patch.object(cobaya.model.Model,'__init__',forbidden))
        result=_validate()
    assert not attempted
    result['native_and_model_entrypoints_globally_forbidden']=True
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=validate()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':result['status'],'native_calls':0}))


if __name__=='__main__':main()
