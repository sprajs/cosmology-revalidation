"""Separate complete accuracy-two correction, retaining every original slot.

This executor never updates an accuracy-one record. Physical execution is an
explicit CLI action after fresh scientific gates and a successful pilot.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import native_accuracy_runtime as rt
import native_accuracy_pilot as pilot


def prerequisites(evidence,refinement_work,pilot_work):
    screen=rt.ROOT/evidence['screen_path'];review=rt.ROOT/evidence['review_path']
    refinement=rt.contracts.prerequisite_snapshot(screen,review,refinement_work)
    checked=pilot.verify_completed(pilot_work)
    assert checked['plan']['evidence']==evidence
    assert refinement['parent_target_identity']==evidence['parent_proposal_identity']
    bindings=rt.merge(checked['input_sha256'],
        {rt.relative(p):rt.digest(p) for p in sorted(Path(refinement_work).rglob('*')) if p.is_file()})
    return {'refinement_work':rt.relative(refinement_work),'pilot_work':rt.relative(pilot_work),
            'refinement':refinement,'pilot_plan_identity':checked['plan']['identity'],
            'pilot_summary_sha256':rt.digest(Path(pilot_work)/'summary.json'),'input_sha256':bindings}


def parent_slots(evidence,workers):
    """The evidence consumer has already freshly qualified this exact parent."""
    from exact_correction import verify_record
    old=rt.original.sealed_read(rt.ROOT/evidence['original_plan_path'])
    summary_path=rt.ROOT/old['correction_summary_path'];summary=json.loads(summary_path.read_text())
    selection_path=rt.ROOT/summary['selection_path'];selection=json.loads(selection_path.read_text())
    assert len(selection['points'])==2000
    assert np.array_equal(selection['groups'],np.repeat(np.arange(4),500))
    refs={r['parent_index']:r for r in evidence['references']};assert len(refs)==32
    requests=[];ordinals=[0]*workers;new_index=0
    for index,point in enumerate(selection['points']):
        path=selection_path.parent/f'{index:05d}.json';checksum=rt.digest(path)
        assert checksum==summary['native_record_sha256'][rt.relative(path)]==evidence['input_sha256'][rt.relative(path)]
        native=json.loads(path.read_text());verify_record(native)
        assert native['index']==index and native['point']==point and native['status']=='finite'
        worker=None;ordinal=None
        if index not in refs:
            worker=new_index%workers;ordinal=ordinals[worker];ordinals[worker]+=1;new_index+=1
        else:
            assert refs[index]['point']==point and refs[index]['native_record_sha256']==checksum
            assert refs[index]['original_location']==selection['locations'][index]
        requests.append({'index':index,'point':point,'group':selection['groups'][index],
             'location':selection['locations'][index],'parent_native_path':rt.relative(path),'parent_native_sha256':checksum,
             'origin':'reused_fixed32' if index in refs else 'persistent_native',
             'audit_index':refs[index]['audit_index'] if index in refs else None,'worker':worker,'ordinal':ordinal})
    assert new_index==1968
    # Nondecreasing chronology retains repeated selected RLE slots verbatim.
    for group in range(4):
        rows=[r for r in requests if r['group']==group]
        assert len({r['location']['chain'] for r in rows})==1
        for left,right in zip(rows,rows[1:]):
            assert left['location']['expanded_index']<=right['location']['expanded_index']
            assert left['location']['row']<=right['location']['row']
            if left['location']['expanded_index']==right['location']['expanded_index']:
                assert left['point']==right['point'] and left['location']['row']==right['location']['row']
    return {'selection_path':rt.relative(selection_path),'selection_sha256':rt.digest(selection_path),
            'correction_summary_path':rt.relative(summary_path),'correction_summary_sha256':rt.digest(summary_path),
            'requests':requests}


def target_identity(plan):
    return rt.identity({'kind':'native_accuracy2_full_importance_target',
        'native_configuration':plan['evidence']['native_configuration'],
        'original_parent_proposal_identity':plan['evidence']['parent_proposal_identity'],
        'original_selection':{k:v for k,v in plan['selection'].items() if k!='requests'},
        'selected_slots':plan['requests'],'evidence':plan['evidence'],
        'prerequisites':plan['prerequisites'],'source_sha256':plan['source_sha256'],
        'design':plan['design'],'converter_contract':plan['converter_contract']})


def prepare(screen,review,refinement_work,pilot_work,work,workers=4):
    assert isinstance(workers,int) and 1<=workers<=4
    plan=rt.base_plan(screen,review,'full_accuracy2_correction')
    plan['workers']=workers
    plan['prerequisites']=prerequisites(plan['evidence'],refinement_work,pilot_work)
    plan['selection']=parent_slots(plan['evidence'],workers)
    plan['requests']=plan['selection']['requests']
    metadata=[rt.old_reference(r)[0]['metadata'] for r in plan['evidence']['references']]
    assert all(x==metadata[0] for x in metadata)
    plan['expected_native_metadata']=metadata[0]
    plan['numerical_target_identity']=target_identity(plan)
    return rt.finalize_plan(plan,work)


def verify_plan(path):
    plan,info=rt.verify_base(path)
    assert plan['kind']=='full_accuracy2_correction' and 1<=plan['workers']<=4
    original=prerequisites(plan['evidence'],rt.ROOT/plan['prerequisites']['refinement_work'],
                           rt.ROOT/plan['prerequisites']['pilot_work'])
    assert original==plan['prerequisites']
    selection=parent_slots(plan['evidence'],plan['workers'])
    assert selection==plan['selection'] and plan['requests']==selection['requests']
    metadata=[rt.old_reference(r)[0]['metadata'] for r in plan['evidence']['references']]
    assert all(x==plan['expected_native_metadata'] for x in metadata)
    assert plan['numerical_target_identity']==target_identity(plan)
    return plan,info


def parent_record(request):
    from exact_correction import verify_record
    path=rt.ROOT/request['parent_native_path'];assert rt.digest(path)==request['parent_native_sha256']
    parent=json.loads(path.read_text());verify_record(parent)
    assert parent['point']==request['point'] and parent['index']==request['index']
    return parent


def build_record(plan,request,value,raw):
    assert rt.canonical(value['metadata'])==rt.canonical(plan['expected_native_metadata'])
    parent=parent_record(request)
    density=rt.contracts.accuracy2_density(parent,value['loglikes'],value['logpriors'],value['logpost'],value['derived'])
    return {'index':request['index'],'point':request['point'],'group':request['group'],'location':request['location'],
            'parent_native_path':request['parent_native_path'],'parent_native_sha256':request['parent_native_sha256'],
            'numerical_target_identity':plan['numerical_target_identity'],
            'native1_logpost':parent['exact_logpost'],'native1_loglikes':parent['exact_loglikes'],
            'native1_logweight':parent['log_weight'],'native1_derived':parent['derived'],
            'proposal_logpost':parent['proposal_logpost'],'proposal_loglikes':parent['proposal_loglikes'],
            'native2_logpost':density['accuracy2_logpost'],'native2_loglikes':density['accuracy2_loglikes'],
            'native2_logpriors':density['accuracy2_logpriors'],'native2_derived':density['accuracy2_derived'],
            'log_weight':density['accuracy2_logweight'],'log_weight_by_composition':density['accuracy2_logweight_by_composition'],
            'density_closure':density['closure'],'failed_checks':[]}


def reuse_fixed32(work,plan):
    for request in plan['requests']:
        if request['origin']!='reused_fixed32':continue
        paths=rt.request_paths(work,request['index'])
        assert not any(p.exists() for p in paths.values())
        ref=plan['evidence']['references'][request['audit_index']]
        assert ref['parent_index']==request['index'] and ref['point']==request['point']
        value,raw=rt.old_reference(ref)
        row=build_record(plan,request,value,raw)
        rt.save_spectrum(paths['spectrum'],raw)
        row.update(binding=rt.request_identity(plan,request),origin='reused_fixed32',value=value,
              status='finite_native_accuracy2',native_logposterior_invocations=0,reference=ref,
              original_reference_status=ref['original_status'],original_reference_failed_checks=ref['original_failed_checks'],
              spectrum_path=rt.relative(paths['spectrum']),spectrum_sha256=rt.digest(paths['spectrum']))
        rt.write_new(paths['record'],row)


def summarize_rows(rows,groups):
    from exact_correction import summarize
    if any(r['status']!='finite_native_accuracy2' for r in rows):
        return {'status':'failed_or_incomplete_native_accuracy2_execution','posterior':None,
                'failures':[{'index':i,'status':r['status']} for i,r in enumerate(rows) if r['status']!='finite_native_accuracy2']}
    pure=[{'status':'finite','point':r['point'],'derived':r['native2_derived'],'log_weight':r['log_weight']} for r in rows]
    diagnostics=summarize(pure,np.asarray(groups))
    passed=diagnostics['status']=='passed_importance_weight_gates' and not diagnostics['failed_gates']
    return {'status':'qualified_native_accuracy2_importance_posterior' if passed else 'failed_native_accuracy2_importance_gates',
            'posterior':diagnostics['posterior'] if passed else None,'importance_diagnostics':diagnostics}


def replay(work):
    work=Path(work);plan,_=verify_plan(work/'plan.json');rows,files=rt.collect(work,plan)
    process_failures=rt.process_integrity(work,plan,rows);seen={}
    for request,row in zip(plan['requests'],rows):
        if row['status']!='finite_native_accuracy2':continue
        assert row['origin']==request['origin']
        raw=rt.load_spectrum(rt.ROOT/row['spectrum_path'])
        expected=build_record(plan,request,row['value'],raw)
        assert all(row[k]==v for k,v in expected.items())
        if request['origin']=='reused_fixed32':
            ref=plan['evidence']['references'][request['audit_index']]
            assert row['reference']==ref and row['original_reference_status']==ref['original_status']
            assert row['original_reference_failed_checks']==ref['original_failed_checks']
            value,original_raw=rt.old_reference(ref)
            assert row['value']==value and all(np.array_equal(raw[k],v) for k,v in original_raw.items())
        else:
            key=(request['worker'],rt.identity(request['point']))
            repeat=rt.compare(row['value'],raw,*seen[key],plan['design']['pilot_tolerances']) if key in seen else None
            assert row['repeat_replay']==repeat
            if repeat:assert not repeat['failed_checks']
            seen[key]=(row['value'],raw)
    result=summarize_rows(rows,[r['group'] for r in plan['requests']])
    if process_failures:
        result.update(status='failed_or_incomplete_native_accuracy2_execution',posterior=None)
    result.update(schema='native-accuracy2-full-correction-v2',numerical_accuracy=2,
       numerical_target_identity=plan['numerical_target_identity'],
       parent_proposal_target_identity=plan['evidence']['parent_proposal_identity'],
       native_configuration=plan['evidence']['native_configuration'],
       plan_path=rt.relative(work/'plan.json'),plan_sha256=rt.digest(work/'plan.json'),plan_identity=plan['identity'],
       source_sha256=plan['source_sha256'],process_failures=process_failures,
       original_selected_slots=len(plan['requests']),
       reused_fixed32=sum(r.get('origin')=='reused_fixed32' for r in rows),
       known_native_logposterior_invocations=sum(r['native_logposterior_invocations'] for r in rows),
       unknown_count_attempts=sum(bool(r.get('unknown_native_invocation_count')) for r in rows),
       records=[{'index':i,'status':r['status'],'path':r.get('record_path'),'sha256':r.get('record_sha256')}
                for i,r in enumerate(rows)],
       limits='Qualified only at the declared accuracy2 configuration. Fixed32 accuracy2-to3 checks do not prove global numerical convergence, mode coverage, or empirical model adequacy.')
    execution=work/'execution.json'
    if execution.exists():files[rt.relative(execution)]=rt.digest(execution)
    return result,files,rows


def main():
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','execute','report','worker'])
    p.add_argument('--work',type=Path,required=True);p.add_argument('--screen',type=Path);p.add_argument('--thermal-review',type=Path)
    p.add_argument('--refinement-work',type=Path);p.add_argument('--pilot-work',type=Path)
    p.add_argument('--workers',type=int,default=4);p.add_argument('--worker',type=int);p.add_argument('--output',type=Path)
    a=p.parse_args()
    if a.action=='prepare':
        assert a.screen and a.thermal_review and a.refinement_work and a.pilot_work
        plan=prepare(a.screen,a.thermal_review,a.refinement_work,a.pilot_work,a.work,a.workers)
        print({'status':'prepared_without_physical_calls','numerical_target_identity':plan['numerical_target_identity']});return
    if a.action=='worker':
        if not rt.persistent_worker(a.work,a.worker,verify_plan,build_record):raise SystemExit(1)
        return
    plan,_=verify_plan(a.work/'plan.json')
    if a.action=='execute':
        rt.execute_workers(a.work,Path(__file__).name,before_workers=lambda:reuse_fixed32(a.work,plan))
        result,files,_=replay(a.work);rt.seal_summary(a.work,plan,result,files,a.output)
    else:
        result,files,_=replay(a.work);rt.verify_saved_summary(a.work,result,files,a.output)
    print({'status':result['status'],'known_invocations':result['known_native_logposterior_invocations'],
           'unknown_attempts':result['unknown_count_attempts']})


if __name__=='__main__':main()
