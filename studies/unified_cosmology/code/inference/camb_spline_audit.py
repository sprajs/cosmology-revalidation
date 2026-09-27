"""Prepare a typed, bounded CAMB1.6.6 source-control comparison. No physics."""
import argparse
import json
from pathlib import Path
import subprocess

import native_accuracy_runtime as runtime
from target_identity import canonical

ROOT=runtime.ROOT
HERE=Path(__file__).resolve().parent
DESIGN=HERE/'camb-spline-audit-design.json'
BUILD=HERE.parent/'external_probes/camb_spline_build.py'
FINALIZE=HERE.parent/'external_probes/camb_spline_finalize.py'
PROOF=HERE.parent/'external_probes/camb_spline_indexing.f90'


def sources():
    return {runtime.relative(p):runtime.digest(p) for p in [Path(__file__),DESIGN,BUILD,FINALIZE,PROOF]}


def verify_build(path):
    value=json.loads(Path(path).read_text())
    assert value['status']=='built_matched_controls_no_physical_calls' and value['physical_calls']==0
    package_path=Path(path).parent/'package-record.json'
    package=json.loads(package_path.read_text())
    assert package['status']=='matched_python_data_and_different_declared_binaries'
    assert package['original_build_sha256']==runtime.digest(path)
    runtime.verify_hashes(package['source_sha256'])
    assert package['source_sha256'][runtime.relative(FINALIZE)]==runtime.digest(FINALIZE)
    runtime.verify_hashes(value['source_sha256'])
    assert value['source_sha256'][runtime.relative(BUILD)]==runtime.digest(BUILD)
    for row in value['assets']:assert runtime.digest(ROOT/row['path'])==row['sha256']
    original=value['installed_provenance'];runtime.verify_hashes(original['installed_file_sha256'])
    assert original['matching_wheel_payload_files']==41
    record={runtime.relative(path):runtime.digest(path),runtime.relative(package_path):runtime.digest(package_path)}
    record.update(original['installed_file_sha256'])
    for arm,row in value['builds'].items():
        runtime.verify_hashes(package['packages'][arm]['file_sha256'])
        record.update(package['packages'][arm]['file_sha256'])
        assert package['packages'][arm]['binary_sha256']==row['binary_sha256']
        tree=(ROOT/row['binary_path']).parent.parent
        for name,expected in row['source_sha256'].items():assert runtime.digest(tree/name)==expected
        assert runtime.digest(ROOT/row['binary_path'])==row['binary_sha256']
        runtime.verify_hashes(row['environment_pth_sha256'])
        assert runtime.digest(ROOT/row['log_path'])==row['log_sha256']
        record.update({runtime.relative(tree/name):expected for name,expected in row['source_sha256'].items()})
        record[row['binary_path']]=row['binary_sha256'];record.update(row['environment_pth_sha256'])
    assert set(value['changed_source_files'])=={'fortran/cmbmain.f90'}
    old=value['patch']['old'];new=value['patch']['new']
    control=(ROOT/value['builds']['control']['binary_path']).parent.parent/'fortran/cmbmain.f90'
    fixed=(ROOT/value['builds']['fixed']['binary_path']).parent.parent/'fortran/cmbmain.f90'
    assert control.read_text().count(old)==1 and control.read_text().replace(old,new)==fixed.read_text()
    return value,record


def prepare(screen,review,build_path,work):
    assert not work.exists(),'Fresh preparation directory required.'
    design=json.loads(DESIGN.read_text());source=sources()
    with runtime.guards_without_physics():
        evidence,two=runtime.evidence(screen,review)
        one,two=runtime.contracts.configuration_pair(evidence['settings'])
    build,bound=verify_build(build_path)
    compiler=ROOT/build['compiler']['path'];assert runtime.digest(compiler)==build['compiler']['sha256']
    work.mkdir(parents=True)
    proof_binary=work/'indexing-proof'
    command=[str(compiler),'-O0','-fcheck=all',str(PROOF),'-o',str(proof_binary)]
    subprocess.run(command,check=True)
    proof_output=subprocess.check_output([str(proof_binary)],text=True)
    assert 'PASS: old first knot is sentinel; old remaining knots shift; fixed knots match.' in proof_output
    references=[evidence['references'][i] for i in design['reference_indices']]
    original_python=ROOT/'.work/unified-cosmology/external-probes/.modern-venv/bin/python'
    arms={'original_wheel':{'python':str(original_python.relative_to(ROOT)),'implementation':build['installed_provenance']},
          **{name:{'python':str(Path(row['environment_path'])/'bin/python'),'implementation':row}
             for name,row in [('source_control',build['builds']['control']),('source_fixed',build['builds']['fixed'])]}}
    # A symlinked interpreter resolves outside the venv. Preserve launch paths, not realpaths.
    arms['original_wheel']['python']=str(original_python.relative_to(ROOT))
    for name,key in [('source_control','control'),('source_fixed','fixed')]:
        arms[name]['python']=str(Path(build['builds'][key]['environment_path'])/'bin/python')
    requests=[]
    for arm in design['arms']:
        for accuracy in design['accuracies']:
            for ordinal,ref in enumerate(references):
                requests.append({'index':len(requests),'arm':arm,'accuracy':accuracy,'ordinal':ordinal,
                                 'audit_index':ref['audit_index'],'parent_index':ref['parent_index'],'point':ref['point']})
    assert len(requests)==design['requests']
    plan={'schema':'camb166-spline-comparison-prepared-v1','status':'prepared_no_physical_calls',
          'physical_calls':0,'design':design,'source_sha256':source,'evidence':evidence,
          'input_sha256':runtime.merge(evidence['input_sha256'],bound),
          'configurations':{'1':canonical(one),'2':canonical(two)},'arms':arms,'requests':requests,
          'indexing_proof':{'source_sha256':runtime.digest(PROOF),'compiler_sha256':runtime.digest(compiler),
                            'binary_sha256':runtime.digest(proof_binary),'output':proof_output},
          'scope':'Implementation audit only; no corrected posterior, no physical execution, no original target relabelled.'}
    assert sources()==source;verify_build(build_path)
    plan['identity']=runtime.identity(plan)
    runtime.write_new(work/'plan.json',plan)
    return plan


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--screen',type=Path,required=True)
    parser.add_argument('--review',type=Path,required=True)
    parser.add_argument('--build-record',type=Path,required=True)
    parser.add_argument('--work',type=Path,required=True)
    args=parser.parse_args()
    value=prepare(args.screen.resolve(),args.review.resolve(),args.build_record.resolve(),args.work.resolve())
    print(json.dumps({'status':value['status'],'requests':len(value['requests']),'identity':value['identity']}))
