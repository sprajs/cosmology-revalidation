"""Fresh, theory-forbidden typed boundary for the distinct native2 posterior."""
from dataclasses import dataclass
from pathlib import Path
import numpy as np
from scipy.special import logsumexp
import native_accuracy_runtime as rt
import native_accuracy_execute as executor


@dataclass(frozen=True)
class QualifiedNativeAccuracy2:
    numerical_accuracy:int
    numerical_target_identity:str
    native_configuration:dict
    parent_settings:dict
    parent_proposal_target_identity:str
    records:tuple
    selected_points:tuple
    groups:np.ndarray
    locations:tuple
    proposal_logposts:np.ndarray
    logweights:np.ndarray
    normalized_weights:np.ndarray
    posterior_summary:dict
    input_sha256:dict
    plan:dict
    plan_path:Path
    summary_path:Path


def qualify_run(work,summary_path):
    """Reconstruct qualifications and verify all recorded bytes before returning."""
    work=Path(work).resolve();summary_path=Path(summary_path).resolve()
    with rt.guards_without_physics():
        plan,configuration=executor.verify_plan(work/'plan.json')
        report,files,rows=executor.replay(work)
        saved=rt.verify_saved_summary(work,report,files,summary_path)
        assert report['status']=='qualified_native_accuracy2_importance_posterior'
        assert report['original_selected_slots']==2000 and report['reused_fixed32']==32
        assert report['known_native_logposterior_invocations']==1968 and report['unknown_count_attempts']==0
        assert len(rows)==2000 and not report['process_failures']
        assert report['importance_diagnostics']['status']=='passed_importance_weight_gates'
        assert not report['importance_diagnostics']['failed_gates'] and report['posterior'] is not None
        inputs=rt.merge(plan['evidence']['input_sha256'],plan['prerequisites']['input_sha256'],
                        plan['source_sha256'],plan['validation_sha256'],files,
                        {rt.relative(p):rt.digest(p) for p in [work/'plan.json',work/'summary.json',work/'record-hashes.json',summary_path]})
        rt.verify_hashes(inputs)
        weights=np.array([r['log_weight'] for r in rows]);groups=np.array([r['group'] for r in plan['requests']])
        assert np.array_equal(groups,np.repeat(np.arange(4),500))
        return QualifiedNativeAccuracy2(numerical_accuracy=2,numerical_target_identity=plan['numerical_target_identity'],
            native_configuration=configuration,parent_settings=plan['evidence']['settings'],
            parent_proposal_target_identity=plan['evidence']['parent_proposal_identity'],records=tuple(rows),
            selected_points=tuple(r['point'] for r in plan['requests']),groups=groups,
            locations=tuple(r['location'] for r in plan['requests']),
            proposal_logposts=np.array([r['proposal_logpost'] for r in rows]),logweights=weights,
            normalized_weights=np.exp(weights-logsumexp(weights)),posterior_summary=rt.plain(saved),input_sha256=inputs,
            plan=plan,plan_path=work/'plan.json',summary_path=summary_path)
