"""Verify and compact saved synthetic gradient checks, without model evaluations."""
import argparse
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[3]

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def summarize(folder,archive):
    summary_path=folder/'summary.json';raw=json.loads(summary_path.read_text())
    runner=Path(__file__).with_name('bayesn_synthetic_gradients.py')
    design=Path(__file__).with_name('bayesn-synthetic-gradient-design.json')
    assert raw['source_sha256']==sha(runner)and raw['design_sha256']==sha(design)
    assert sha(folder/'execution-design.json')==raw['execution_design_sha256']
    for p,h in raw['input_sha256'].items():assert sha(ROOT/p)==h
    for name,h in raw['reviewed_package_sha256'].items():assert sha(runner.parent/'bayesn_heldout'/name)==h
    registry_path=ROOT/'studies/infrared/results/bayesn-nir-forward-gate.json'
    registry=json.loads(registry_path.read_text())['execution_design']['input_sha256']
    expected=raw['archive_opened_input_sha256'];matched=[];unlisted=[]
    for p,h in expected.items():
        assert '/photometry/'not in p.lower()and not p.lower().endswith('truth.json')
        assert sha(archive/p)==h
        if p in registry:
            assert registry[p]==h
            matched.append(p)
        else:unlisted.append(p)
    records=[];smooth=0;intentional=0;failures=[];all_fd=0;native=[]
    for item in raw['states']:
        p=ROOT/item['raw_path'];assert sha(p)==item['raw_sha256']
        state=json.loads(p.read_text());groups={}
        for v in state['finite_differences']:
            groups.setdefault(v['coordinate'],[]).append(v)
            all_fd+=1
        assert len(groups)==47 and all(len(v)==11 for v in groups.values())
        for coordinate in state['coordinate_summary']:
            steps=groups[coordinate['coordinate']]
            ok=[not s['crosses_integer_phase']and abs(s['finite_difference']-s['autodiff'])/(1+abs(s['autodiff']))<1e-4 for s in steps]
            agreement=any(a and b for a,b in zip(ok[:-1],ok[1:]))
            assert agreement==coordinate['adjacent_steps_agree']
            if coordinate['status']=='intentional_phase_knot_nondifferentiability':
                assert state['state']['kind']=='phase_control'and state['state']['control_delta_tau_days']==0 and coordinate['coordinate']=='logit_tau'
                intentional+=1
            elif agreement:smooth+=1
            else:failures.append({'state':item['state_index'],'coordinate':coordinate['coordinate']})
        native.append(state['native_potential_gradient_max_error'])
        supported=[c['best_noncrossing_normalized_error']for c in state['coordinate_summary']if c['best_noncrossing_normalized_error']is not None]
        records.append({'state_index':item['state_index'],'selection':item['selection'],'raw_path':item['raw_path'],'raw_sha256':item['raw_sha256'],
            'nearest_integer_phase_distance_days':state['nearest_integer_phase_distance'],
            'independent_density_error':state['independent_density_error'],
            'native_potential_value_error':state['native_potential_value_error'],
            'native_potential_gradient_max_error':state['native_potential_gradient_max_error'],
            'maximum_best_noncrossing_normalized_error':max(supported),
            'maximum_smallest_step_normalized_error_excluding_exact_knot':max(c['smallest_step_error']for c in state['coordinate_summary']if c['status']!='intentional_phase_knot_nondifferentiability'),
            'physical_tau_sweep':state['physical_tau_sweep']})
    assert len(records)==9 and all_fd==4653
    # This exact-equality check is a reporting safeguard added after inspecting the
    # completed diagnostic, not a tuned scientific tolerance or a new model run.
    native_equal=all(x==0 for x in native)
    return {'status':'passed_local_smooth_gradient_checks_with_expected_phase_kink'if not failures and native_equal else 'gradient_checks_require_investigation',
        'scope':'Synthetic optical only, failed parent chains; not posterior qualification or divergence-cause identification.',
        'states':records,'smooth_coordinate_checks_passed':smooth,'intentional_exact_knot_checks':intentional,
        'finite_difference_step_records_verified':all_fd,'smooth_failures':failures,
        'native_potential_gradients_bitwise_equal':native_equal,'native_gradient_gate_timing':'Post-execution defensive equality check; main predeclared finite-difference thresholds unchanged.',
        'maximum_native_potential_gradient_difference':max(native),
        'source_sha256':{str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),runner,design]},
        'raw_summary':str(summary_path.relative_to(ROOT)),'raw_summary_sha256':sha(summary_path),
        'execution_design_sha256':raw['execution_design_sha256'],'input_sha256':raw['input_sha256'],
        'reviewed_package_sha256':raw['reviewed_package_sha256'],
        'parent_gate_registry_sha256':sha(registry_path),'opened_archive_files_verified':len(expected),
        'opened_archive_files_matching_parent_gate':len(matched),'opened_files_not_listed_in_parent_gate':unlisted,
        'archive_opened_input_sha256':expected,
        'interpretation':['Every predeclared smooth coordinate comparison has adjacent supported finite-difference steps agreeing with autodiff; this does not prove globally correct gradients.',
            'At the exact Hsiao integer-phase node, the two one-sided derivatives differ. floor=ceil makes the pointwise autodiff interpolation contribution special; a unique smooth derivative is not defined there.',
            'Both nearby phase controls pass. The two retained divergence-associated states also pass, but these are not the actual failed leapfrog locations.',
            'The original sampling failures are preserved. No new sampler, smoothing choice, target change, held-out score or cosmological inference is endorsed.']}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--audit-folder',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--archive',type=Path,default=Path.home()/'.local/share/cosmology-revalidation/archive/17487bf659fcbdeeea072221492bac14b04a0a85')
    args=p.parse_args();result=summarize(args.audit_folder.resolve(),args.archive.resolve());assert not args.output.exists()
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k]for k in ['status','smooth_coordinate_checks_passed','intentional_exact_knot_checks','native_potential_gradients_bitwise_equal','opened_archive_files_verified','opened_archive_files_matching_parent_gate','opened_files_not_listed_in_parent_gate']}))
