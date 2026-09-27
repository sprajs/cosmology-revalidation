"""Eight predeclared native2 requests, two persistent workers; fresh work only."""
import argparse
from pathlib import Path
import native_accuracy_runtime as rt


def layout(plan):
    requests=[]
    for worker,sequence in enumerate(plan['design']['pilot_sequences']):
        for ordinal,audit in enumerate(sequence):
            requests.append({'index':len(requests),'worker':worker,'ordinal':ordinal,'audit_index':audit,
                'origin':'persistent_native','point':plan['evidence']['references'][audit]['point']})
    return requests


def prepare(screen,review,work):
    plan=rt.base_plan(screen,review,'persistent_accuracy2_pilot')
    plan.update(workers=plan['design']['pilot_workers'],requests=layout(plan))
    return rt.finalize_plan(plan,work)


def verify_plan(path):
    plan,info=rt.verify_base(path)
    assert plan['kind']=='persistent_accuracy2_pilot' and plan['workers']==2
    assert plan['requests']==layout(plan) and len(plan['requests'])==8
    return plan,info


def build_record(plan,request,value,raw):
    reference=plan['evidence']['references'][request['audit_index']]
    check=rt.compare(value,raw,*rt.old_reference(reference),plan['design']['pilot_tolerances'])
    return {'reference':reference,'reference_replay':check,'failed_checks':check['failed_checks'].copy()}


def replay(work):
    work=Path(work);plan,_=verify_plan(work/'plan.json');rows,files=rt.collect(work,plan)
    failures=rt.process_integrity(work,plan,rows)
    seen={}
    for request,row in zip(plan['requests'],rows):
        if row['status']!='finite_native_accuracy2':failures.append('request:'+str(request['index']));continue
        raw=rt.load_spectrum(rt.ROOT/row['spectrum_path'])
        expected=build_record(plan,request,row['value'],raw)
        assert all(row[k]==v for k,v in expected.items())
        assert not expected['failed_checks']
        key=(request['worker'],rt.identity(request['point']))
        repeat=rt.compare(row['value'],raw,*seen[key],plan['design']['pilot_tolerances']) if key in seen else None
        assert row['repeat_replay']==repeat
        if repeat:assert not repeat['failed_checks']
        seen[key]=(row['value'],raw)
    execution=work/'execution.json'
    if execution.exists():files[rt.relative(execution)]=rt.digest(execution)
    report={'schema':'persistent-native-accuracy2-pilot-v2',
        'status':'failed_or_incomplete_persistent_accuracy2_pilot' if failures else 'passed_persistent_accuracy2_replay_pilot',
        'plan_identity':plan['identity'],'plan_path':rt.relative(work/'plan.json'),'plan_sha256':rt.digest(work/'plan.json'),
        'parent_proposal_target_identity':plan['evidence']['parent_proposal_identity'],
        'native_configuration':plan['evidence']['native_configuration'],
        'failures':failures,'requests':rows,'source_sha256':plan['source_sha256'],
        'known_native_logposterior_invocations':sum(r['native_logposterior_invocations'] for r in rows),
        'unknown_count_attempts':sum(bool(r.get('unknown_native_invocation_count')) for r in rows),
        'posterior_qualified':False,'timing_interpretation':'Per-process construction and per-request wall time under actual contention; no solver-call equivalence.'}
    return report,files


def verify_completed(work):
    with rt.guards_without_physics():
        report,files=replay(work);saved=rt.verify_saved_summary(work,report,files)
        assert report['status']=='passed_persistent_accuracy2_replay_pilot'
        assert report['known_native_logposterior_invocations']==8 and report['unknown_count_attempts']==0
        plan,_=verify_plan(Path(work)/'plan.json')
        bindings=rt.merge(files,plan['evidence']['input_sha256'],plan['source_sha256'],plan['validation_sha256'],
                  {rt.relative(Path(work)/n):rt.digest(Path(work)/n) for n in ['plan.json','summary.json','record-hashes.json']})
        return {'plan':plan,'summary':rt.plain(saved),'input_sha256':bindings}


def main():
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','execute','report','worker'])
    p.add_argument('--work',required=True,type=Path);p.add_argument('--screen',type=Path);p.add_argument('--thermal-review',type=Path)
    p.add_argument('--worker',type=int);p.add_argument('--output',type=Path);a=p.parse_args()
    if a.action=='prepare':
        assert a.screen and a.thermal_review;result=prepare(a.screen,a.thermal_review,a.work)
        print({'status':'prepared_without_physical_calls','plan_identity':result['identity']});return
    if a.action=='worker':
        assert a.worker is not None
        if not rt.persistent_worker(a.work,a.worker,verify_plan,build_record):raise SystemExit(1)
        return
    plan,_=verify_plan(a.work/'plan.json')
    if a.action=='execute':
        rt.execute_workers(a.work,Path(__file__).name)
        report,files=replay(a.work);rt.seal_summary(a.work,plan,report,files,a.output)
    else:
        report,files=replay(a.work);rt.verify_saved_summary(a.work,report,files,a.output)
    print({'status':report['status'],'known_invocations':report['known_native_logposterior_invocations'],
           'unknown_attempts':report['unknown_count_attempts']})


if __name__=='__main__':main()
