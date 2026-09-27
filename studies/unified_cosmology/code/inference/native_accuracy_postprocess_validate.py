"""Known Gaussian targets and corruption controls for native2 aggregation."""
import argparse
from contextlib import ExitStack
import copy
from dataclasses import replace
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
from scipy.special import logsumexp
from scipy.stats import norm

import native_accuracy_postprocess as consumer
from luminosity_sensitivity import weighted_summary


def rejected(fn):
    try:
        fn()
    except (AssertionError, ValueError):
        return True
    raise AssertionError('Corrupted numerical input was accepted.')


def run():
    rng = np.random.default_rng(274019)
    count = 2000
    x = rng.normal(.1, 1.4, count)
    groups = np.repeat(np.arange(4), 500)
    # Deliberately retain two repeated selected slots, in different chains.
    x[167] = x[166]; x[1229] = x[1228]
    points = [{'H0':float(68+v)} for v in x]
    derived = [{'omegam':float(.3+.01*v), 'q0':float(-.5+.1*v),
                'q05':float(-.1+.1*v), 'q1':float(.2+.1*v),
                'j0':float(-.2+.5*v), 'sn_chi2':float((v-.5)**2/4)} for v in x]
    logp = norm.logpdf(x, .2, 1.)
    logq = norm.logpdf(x, .1, 1.4)
    lw = logp-logq
    locations = [{'chain':f'chain.{i//500+1}.txt','expanded_index':int(3*(i%500)),
                  'row':int(3*(i%500))} for i in range(count)]
    for i in [167, 1229]:
        locations[i] = copy.deepcopy(locations[i-1])
    chronology = consumer.chronology(points, groups, locations, logq)
    assert sum(row['repeated_slots'] for row in chronology.values()) == 2
    values, weights = consumer.arrays(points, derived, lw, groups)
    posterior = {k:weighted_summary(v,weights) for k,v in values.items()}
    summary = consumer.summaries(values, weights, groups, posterior)
    matrix = np.column_stack([values[k] for k in summary['weighted_covariance']['parameter_order']])
    direct_cov = sum(float(w)*np.outer(row-weights@matrix,row-weights@matrix)
                     for w,row in zip(weights,matrix))
    error_cov = float(np.max(abs(direct_cov-np.array(summary['weighted_covariance']['matrix']))))
    assert error_cov < 1e-12
    assert abs(weights@x-.2) < .08 and abs(np.sqrt(weights@(x-weights@x)**2)-1) < .08
    for key in ('q0','q05','q1','j0'):
        assert summary['conditional_sign_fractions'][key+'_below_zero']['fraction'] == float(weights@(values[key]<0))
    covariance_names = summary['weighted_covariance']['parameter_order']
    selected = {k:values[k] for k in covariance_names}
    boot = consumer.quantiles.block_quantile_precision(selected, lw, groups)
    for k in selected:
        assert boot['baseline_quantiles'][k] == posterior[k]['quantiles_025_16_50_84_975']
    # Probe removal has a separate analytic Gaussian target: subtracting
    # -(x-.5)^2/(2*2^2) from log N(.2,1) gives mean .1, variance 4/3.
    design = json.loads(consumer.omission.DESIGN.read_text())
    gates = json.loads(consumer.omission.GATES.read_text())['overlap_gates']
    expected = design['required_parent_components']
    removed = design['alternatives']['sn']['removed_components']
    views, rows = [], []
    for i,v in enumerate(x):
        components = {k:0. for k in expected}
        components['released_sn'] = float(-.5*((v-.5)/2)**2)
        components['planck_2018_lowl.TT'] = float(logp[i]-components['released_sn'])
        a,b = consumer.remove_components(points[i], derived[i], components, logp[i],
                    logq[i], lw[i], removed, expected, 1e-9)
        views.append(a); rows.append(b)
    result = consumer.omission.summarize_omission(views, rows, groups, gates)
    assert result['status'] == 'qualified_conditional_probe_omission', result['failed_gates']
    w = np.exp(np.array([r['target_logweight'] for r in rows])-logsumexp([r['target_logweight'] for r in rows]))
    mean, sd = float(w@x), float(np.sqrt(w@(x-w@x)**2))
    assert abs(mean-.1) < .08 and abs(sd-np.sqrt(4/3)) < .08
    independent_w = norm.pdf(x,.1,np.sqrt(4/3))/norm.pdf(x,.1,1.4)
    independent_w /= independent_w.sum()
    weight_error = float(np.max(abs(w-independent_w)))
    assert weight_error < 1e-14
    # Strong support failure remains a failure with no reported posterior.
    failed_rows = copy.deepcopy(rows)
    for i,row in enumerate(failed_rows):
        row['target_logweight'] += 20 if i < 500 else 0
    failed = consumer.omission.summarize_omission(views, failed_rows, groups, gates)
    assert failed['posterior'] is None and failed['conditional_sign_fractions'] is None
    assert failed['weighted_covariance'] is None and failed['failed_gates']
    controls = {}
    bad = copy.deepcopy(locations); bad[90]['expanded_index'] = 0
    controls['decreasing_chronology'] = rejected(lambda:consumer.chronology(points,groups,bad,logq))
    bad_points = copy.deepcopy(points); bad_points[167]['H0'] += 1
    controls['repeat_different_point'] = rejected(lambda:consumer.chronology(bad_points,groups,locations,logq))
    bad_q = logq.copy(); bad_q[167] += 1
    controls['repeat_different_proposal_density'] = rejected(lambda:consumer.chronology(points,groups,locations,bad_q))
    controls['interleaved_chains'] = rejected(lambda:consumer.chronology(points,np.tile(np.arange(4),500),locations,logq))
    bad_derived = copy.deepcopy(derived); bad_derived[5]['q0'] = float('nan')
    controls['nonfinite_native2_derived'] = rejected(lambda:consumer.arrays(points,bad_derived,lw,groups))
    bad_posterior = copy.deepcopy(posterior); bad_posterior['H0']['mean'] += .1
    controls['stale_parent_summary'] = rejected(lambda:consumer.summaries(values,weights,groups,bad_posterior))
    controls['wrong_density_weight'] = rejected(lambda:consumer.remove_components(points[-1],derived[-1],components,
        logp[-1],logq[-1],lw[-1]+.1,removed,expected,1e-9))
    controls['duplicate_removal'] = rejected(lambda:consumer.remove_components(points[-1],derived[-1],components,
        logp[-1],logq[-1],lw[-1],removed+removed,expected,1e-9))
    controls['missing_component'] = rejected(lambda:consumer.remove_components(points[-1],derived[-1],
        {k:v for k,v in components.items() if k!='released_sn'},logp[-1],logq[-1],lw[-1],removed,expected,1e-9))
    return {'status':'passed_native_accuracy2_postprocess_kernel_validation',
            'scope':'Pure kernels only. Typed observational boundary validation remains required before execution.',
            'points':count,'explicit_repeated_slots':2,'seed':274019,
            'parent_Gaussian_mean_error':float(weights@x-.2),
            'removed_Gaussian_mean':mean,'removed_Gaussian_SD':sd,
            'removed_Gaussian_truth':{'mean':.1,'SD':float(np.sqrt(4/3))},
            'independent_normalized_density_weight_error':weight_error,
            'independent_covariance_error':error_cov,
            'negative_controls':controls,'failed_overlap_posterior_withheld':True,
            'quantile_partitions':[10,20],'replicates_per_partition':500,
            'physical_calls':0,
            'source_sha256':{consumer.omission.relative(p):consumer.omission.digest(p)
                for p in [Path(__file__),Path(consumer.__file__),consumer.DESIGN,
                          Path(consumer.omission.__file__),consumer.omission.DESIGN,
                          consumer.omission.GATES,Path(consumer.quantiles.__file__)]}}


def boundary_tests():
    """Exercise the real consumer with typed synthetic parents, no physics.

    Only qualify_run is replaced. The real qualifier has a separate execution
    validator; this test checks that its consumer calls and respects that
    boundary, real file hashes, numerical accuracy and original pure gates.
    """
    import native_accuracy_qualify as qualifier
    rng = np.random.default_rng(274020)
    n = 2000
    x = rng.normal(.1, 1.4, n)
    epsilon = rng.uniform(-.5, .5, n)
    foreground = rng.uniform(0., 2., n)
    for i in [167, 1229]:
        x[i], epsilon[i], foreground[i] = x[i-1], epsilon[i-1], foreground[i-1]
    points = tuple({'H0':float(68+v), 'epsilon':float(e), 'A_fg':float(f)}
                   for v,e,f in zip(x,epsilon,foreground))
    derived = tuple({'omegam':float(.3+.01*v), 'q0':float(-.5+.1*v),
                     'q05':float(-.1+.1*v), 'q1':float(.2+.1*v),
                     'j0':float(-.2+.5*v), 'sn_chi2':float((v-.5)**2/4)} for v in x)
    logp, logq = norm.logpdf(x,.2,1.), norm.logpdf(x,.1,1.4)
    lw = logp-logq
    groups = np.repeat(np.arange(4),500)
    locations = tuple({'chain':f'chain.{i//500+1}.txt', 'expanded_index':int(3*(i%500)),
                       'row':int(3*(i%500))} for i in range(n))
    for i in [167,1229]: locations[i].update(locations[i-1])
    values,weights = consumer.arrays(points,derived,lw,groups)
    posterior = {k:weighted_summary(v,weights) for k,v in values.items()}
    omission_design = json.loads(consumer.omission.DESIGN.read_text())
    expected = omission_design['required_parent_components']
    records = []
    for i,v in enumerate(x):
        components = {k:0. for k in expected}
        components['released_sn'] = float(-.5*((v-.5)/2)**2)
        components['bao.desi_dr2'] = float(-.1*(v+.1)**2)
        components['act_dr6_lenslike.ACTDR6LensLike'] = float(-.01*v**2)
        components['SPT2023_lensing'] = float(-.015*(v-.1)**2)
        components['planck_2018_lowl.TT'] = float(logp[i]-sum(components.values()))
        records.append({'native2_derived':derived[i], 'native2_loglikes':components,
                        'native2_logpost':float(logp[i])})
    configuration = {'theory':{'camb':{'extra_args':{name:2 for name in
                         ('AccuracyBoost','lAccuracyBoost','lSampleBoost')}}},
                     'likelihood':{name:{} for name in expected},
                     'params':{'H0':{'prior':{'min':50,'max':90}},
                               'epsilon':{'prior':{'min':-.5,'max':.5}},
                               'A_fg':{'prior':{'min':0,'max':2}}}}
    cases = []
    with tempfile.TemporaryDirectory(prefix='native2-postprocess-',dir=consumer.ROOT/'.work') as tmp:
        root = Path(tmp)
        source = consumer.sources()
        input_file = root/'sealed-synthetic-parent.json'
        input_file.write_text('{"scope":"SYNTHETIC parent only, not an observed posterior"}\n')
        original_input = input_file.read_bytes()
        validation = root/'studies/unified_cosmology/results/inference/native-accuracy-postprocess-validation.json'
        validation.parent.mkdir(parents=True)
        proof = {'status':'passed_native_accuracy2_postprocess_validation',
                 'source_sha256':source,'physical_calls':0,'scope':'temporary fixture only'}
        def restore():
            input_file.write_bytes(original_input)
            validation.write_text(json.dumps(proof)+'\n')
        restore()
        original = qualifier.QualifiedNativeAccuracy2(
            numerical_accuracy=2,numerical_target_identity='a'*64,
            native_configuration=configuration,parent_settings={'model':'cpl','evolution':'linear',
              'sample':'synthetic','calibration':'official_planck'},
            parent_proposal_target_identity='b'*64,records=tuple(records),selected_points=points,
            groups=groups,locations=locations,proposal_logposts=logq,logweights=lw,
            normalized_weights=weights,posterior_summary={'posterior':posterior},
            input_sha256={consumer.omission.relative(input_file):consumer.omission.digest(input_file)},
            plan={'scope':'synthetic'},plan_path=root/'plan.json',summary_path=root/'summary.json')
        def invoke(target=original,action='summary',qualifier_error=None,mutate_summary=None):
            with ExitStack() as stack:
                stack.enter_context(qualifier.rt.guards_without_physics())
                stack.enter_context(patch.object(consumer,'ROOT',root))
                guard = stack.enter_context(patch.object(qualifier,'qualify_run',
                    side_effect=qualifier_error,return_value=target))
                if mutate_summary is not None:
                    real = consumer.summaries
                    def mutate(*args):
                        answer = real(*args); mutate_summary(); return answer
                    stack.enter_context(patch.object(consumer,'summaries',mutate))
                answer = consumer.actual(root/'work',root/'summary.json',action)
                assert guard.call_count == 1
                return answer
        for action in ['summary','quantiles','omit-sn','omit-bao','omit-lensing']:
            restore(); before = copy.deepcopy(configuration)
            answer = invoke(action=action)
            assert configuration == before, 'Consumer mutated parent configuration.'
            assert answer['numerical_accuracy']==2 and answer['physical_calls']==0
            assert answer['parent_numerical_target_identity']=='a'*64
            assert answer['parent_proposal_target_identity']=='b'*64
            assert sum(r['repeated_slots'] for r in answer['chronology'].values())==2
            assert answer['source_sha256']==source
            assert answer['input_sha256'][consumer.omission.relative(validation)]==consumer.omission.digest(validation)
            config = answer['target_definition']['configuration']
            assert config['params']==configuration['params'] and config['theory']==configuration['theory']
            if action.startswith('omit-'):
                removed = omission_design['alternatives'][action[5:]]['removed_components']
                assert set(config['likelihood'])==set(expected)-set(removed)
                assert answer['status']=='qualified_conditional_probe_omission',answer['failed_gates']
                direct = lw-np.array([sum(r['native2_loglikes'][k] for k in removed) for r in records])
                direct = np.exp(direct-logsumexp(direct))
                assert abs(answer['posterior']['H0']['mean']-direct@values['H0'])<1e-12
                if action=='omit-sn':
                    assert 'epsilon' not in answer['posterior']
                    assert 'epsilon' in answer['weight_diagnostics']['weighted_summaries_for_diagnostics']
                if action=='omit-lensing': assert answer['unused_auxiliary_coordinates']==['A_fg']
            else:
                assert config==configuration
            cases.append({'case':'typed_'+action,'passed':True,'qualified_parent_calls':1})
        def negative(name,target=original,setup=None,qualifier_error=None,mutate_summary=None):
            restore()
            if setup: setup()
            assert rejected(lambda:invoke(target=target,qualifier_error=qualifier_error,mutate_summary=mutate_summary))
            cases.append({'case':name,'passed':True,'rejected':True})
        negative('unqualified_parent_rejected',qualifier_error=AssertionError('synthetic parent fails gates'))
        restore()
        with patch.object(consumer,'arrays',side_effect=RuntimeError('Aggregation before qualification is forbidden.')):
            assert rejected(lambda:invoke(qualifier_error=AssertionError('synthetic parent fails gates')))
        cases.append({'case':'qualification_precedes_aggregation','passed':True,'rejected':True})
        negative('descriptor_is_not_typed_qualification',target=SimpleNamespace(**original.__dict__))
        negative('wrong_numerical_accuracy',target=replace(original,numerical_accuracy=1))
        negative('changed_parent_input',setup=lambda:input_file.write_text('changed'))
        negative('changed_parent_input_during_summary',mutate_summary=lambda:input_file.write_text('changed'))
        negative('changed_validation_during_summary',mutate_summary=lambda:validation.write_text('{}'))
        def alter_proof(key,value):
            bad=copy.deepcopy(proof);bad[key]=value;validation.write_text(json.dumps(bad))
        negative('failed_validation_proof',setup=lambda:alter_proof('status','failed'))
        negative('stale_validation_source',setup=lambda:alter_proof('source_sha256',{}))
        negative('proof_admits_physical_calls',setup=lambda:alter_proof('physical_calls',1))
        bad_config=copy.deepcopy(configuration);bad_config['theory']['camb']['extra_args']['AccuracyBoost']=1
        negative('wrong_native_configuration',target=replace(original,native_configuration=bad_config))
        negative('stale_normalized_weights',target=replace(original,normalized_weights=np.full(n,1/n)))
        bad=copy.deepcopy(posterior);bad['H0']['mean']+=.1
        negative('stale_qualified_posterior',target=replace(original,posterior_summary={'posterior':bad}))
        bad=copy.deepcopy(locations);bad[99]['expanded_index']=0
        negative('decreasing_typed_locations',target=replace(original,locations=bad))
        # A genuinely failed removal must preserve its diagnostics and withhold
        # the posterior even when its parent was qualified by the typed boundary.
        bad=copy.deepcopy(records)
        for i,row in enumerate(bad):
            if i<500:
                row['native2_loglikes']['released_sn']-=20
                row['native2_loglikes']['planck_2018_lowl.TT']+=20
        restore();failed=invoke(target=replace(original,records=bad),action='omit-sn')
        assert failed['posterior'] is None and failed['conditional_sign_fractions'] is None
        assert failed['weighted_covariance'] is None and failed['failed_gates']
        cases.append({'case':'failed_child_overlap_withholds_posterior','passed':True})
    return {'cases':cases,'count':len(cases),'physical_calls':0,
            'mock_scope':'Temporary synthetic typed parents only; real consumer actual() and all pure aggregation/hash checks run. This does not replace independent validation of the real qualifier.'}


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--boundary',action='store_true',help='Also test the actual consumer against the current typed qualifier API, using synthetic parents only.')
    a=p.parse_args(); assert not a.output.exists()
    result=run()
    if a.boundary:
        result['boundary_validation']=boundary_tests()
        result['status']='passed_native_accuracy2_postprocess_validation'
        result['scope']='Pure kernels and typed consumer boundary, synthetic only. No observational postprocessing or physical calls.'
        result['source_sha256']=consumer.sources()
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ['status','independent_normalized_density_weight_error',
                                          'independent_covariance_error','physical_calls']}))
