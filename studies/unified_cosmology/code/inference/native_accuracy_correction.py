"""Preparation-only contracts for a separate native accuracy-2 correction.

No worker/model constructor or physical execution CLI is provided here. The
frozen design requires a separately authorized persistent-worker replay pilot
and a completed observational qualifier before a new posterior can be admitted.
"""
import copy
import hashlib
import inspect
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
DESIGN = HERE/'native-accuracy-correction-design.json'
SPECTRA = ('tt','ee','bb','te','pp')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def identity(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def sources():
    paths = [Path(__file__),DESIGN,HERE/'native_accuracy_measurement.py',
             HERE/'native_accuracy_validate.py',HERE/'native_posterior_precision.py',
             HERE/'native-posterior-precision-design.json',
             HERE/'native_precision_review_consumer.py',HERE/'native_precision_refinement.py',
             HERE/'native-precision-review-consumer-design.json',
             HERE/'native_precision_refinement_design.json',
             HERE/'native_precision_thermal_review.py',HERE/'native-precision-thermal-review-design.json',
             HERE/'exact_correction.py',HERE/'luminosity_sensitivity.py',
             HERE/'luminosity-sensitivity-design.json',HERE/'target_identity.py']
    return {str(p.relative_to(ROOT)):digest(p) for p in paths}


def configuration_pair(settings):
    """Only build native configuration dictionaries; never construct a model."""
    import native_posterior_precision as original
    design = json.loads(original.DESIGN.read_text())
    one,two = original.native_configurations(settings,design)
    assert_only_boost_change(one,two)
    return one,two


def assert_only_boost_change(one,two):
    from target_identity import canonical
    design = json.loads(DESIGN.read_text())
    assert set(one['theory']) == set(two['theory']) == {'camb'}
    changed = copy.deepcopy(two)
    for name in design['numerical_controls']:
        assert one['theory']['camb']['extra_args'][name] == design['parent_accuracy']
        assert changed['theory']['camb']['extra_args'][name] == design['target_accuracy']
        changed['theory']['camb']['extra_args'][name] = design['parent_accuracy']
    assert canonical(changed) == canonical(one), 'Numerical target changes more than its three boosts.'


def converter_contract():
    """Bind the installed provider conversion source without instantiating it."""
    import importlib.metadata
    from cobaya.theories.camb.camb import CAMB
    source = inspect.getsource(CAMB._get_Cl)
    return {'archive_ell_factor':True,'output_ell_factor':False,'units':'FIRASmuK2',
            'pp':'dimensionless C_phi_phi',
            'cobaya_version':importlib.metadata.version('cobaya'),
            'provider_conversion_source_sha256':hashlib.sha256(source.encode()).hexdigest(),
            'provider_unit_source_sha256':hashlib.sha256(inspect.getsource(CAMB._cmb_unit_factor).encode()).hexdigest(),
            'provider_entrypoint_source_sha256':hashlib.sha256(inspect.getsource(CAMB.get_Cl).encode()).hexdigest(),
            'converter_source_sha256':digest(__file__)}


def raw_from_dl(spectra):
    """Invert the installed CAMB provider convention at valid archived multipoles.

    The provider treats temperature/polarization ell=0 units differently across
    its two output modes. Ordinary physical spectra have zero there; reject a
    nonzero value instead of inventing unavailable underlying-unit information.
    """
    assert set(spectra) == {'ell',*SPECTRA}
    ell = np.asarray(spectra['ell'])
    assert ell.ndim == 1 and len(ell)>=2
    assert np.array_equal(ell,np.arange(len(ell)))
    out = {'ell':ell.copy()}
    ll = ell[1:].astype(float)*(ell[1:].astype(float)+1)
    for key in SPECTRA:
        value = np.asarray(spectra[key],dtype=float)
        assert value.shape == ell.shape and np.isfinite(value).all()
        if key != 'pp':
            assert value[0] == 0., 'Nonzero temperature/polarization ell0 cannot be inverted from this archive alone.'
        out[key] = value.copy()
        out[key][1:] *= 2*np.pi/(ll**(2 if key=='pp' else 1))
    return out


def dl_from_raw(spectra):
    """Synthetic round-trip helper; not an alternate physical spectrum source."""
    raw_from_dl(spectra) # Shared shape/finite/low-ell checks only.
    ell = np.asarray(spectra['ell'])
    ll = ell[1:].astype(float)*(ell[1:].astype(float)+1)
    result = {'ell':ell.copy()}
    for key in SPECTRA:
        result[key] = np.asarray(spectra[key],dtype=float).copy()
        result[key][1:] *= ll**(2 if key=='pp' else 1)/(2*np.pi)
    return result


def accuracy2_density(parent,loglikes,logpriors,logpost,derived):
    """Validate both raw weight equations and preserve distinct native targets."""
    limit = json.loads(DESIGN.read_text())['density_absolute_tolerance']
    assert parent['status']=='finite'
    assert set(loglikes)==set(parent['exact_loglikes'])
    assert set(derived)==set(parent['derived'])
    assert np.isfinite([parent['log_weight'],parent['exact_logpost'],parent['proposal_logpost'],
                       *parent['exact_loglikes'].values(),*loglikes.values(),*logpriors,
                       logpost,*derived.values()]).all()
    prior1=parent['exact_logpost']-sum(parent['exact_loglikes'].values())
    errors = {
        'parent_weight':parent['log_weight']-(parent['exact_logpost']-parent['proposal_logpost']),
        'prior_difference':sum(logpriors)-prior1,
        'native2_density':logpost-sum(loglikes.values())-sum(logpriors),
    }
    direct=logpost-parent['proposal_logpost']
    chained=parent['log_weight']+logpost-parent['exact_logpost']
    errors['weight_composition']=direct-chained
    assert all(abs(x)<=limit for x in errors.values()), 'Native accuracy2 density/prior accounting failed.'
    return {'status':'finite_native_accuracy2_contract',
            'accuracy1_logpost':float(parent['exact_logpost']),
            'proposal_logpost':float(parent['proposal_logpost']),
            'accuracy1_logweight':float(parent['log_weight']),
            'accuracy2_logpost':float(logpost),'accuracy2_logweight':float(direct),
            'accuracy2_logweight_by_composition':float(chained),
            'accuracy2_loglikes':dict(loglikes),'accuracy2_logpriors':list(logpriors),
            'accuracy2_derived':dict(derived),'closure':{k:float(v) for k,v in errors.items()},
            'observational_qualification':False}


def prerequisite_snapshot(screen_path,thermal_review_path,refinement_work):
    """Read-only evidence gate; unavailable/mismatched evidence fails closed.

    This does not qualify a full accuracy2 posterior and does not create a plan,
    change an original record, instantiate a model or execute a native point.
    """
    import native_precision_review_consumer as consumer
    import native_precision_refinement as refinement
    receipt=consumer.verify(Path(screen_path),Path(thermal_review_path))
    assert receipt['reviewed_screen_status']=='numerical_sensitivity_requires_followup'
    flags=receipt['reviewed_diagnostic_flags']
    assert flags and all(x in ('total_centered_RMS','total_centered_max_absolute') or
                         x.startswith('component_variation:') for x in flags)
    work=Path(refinement_work).resolve()
    plan=refinement.verify_plan(work/'plan.json')
    assert plan['verified_original_receipt_identity']==refinement.original.identity(receipt)
    assert plan['parent_target_identity']==receipt['parent_target_identity']
    replay_sha=refinement.require_replay(work,plan)
    report=refinement.run_stage(work,'refine',False,1)
    saved=refinement.original.sealed_read(work/'refine-summary.json')
    assert {k:v for k,v in saved.items() if k!='payload_sha256'}==report
    assert report['status']=='no_large_variation_detected_on_fixed32'
    assert not report.get('diagnostic_flags') and not report.get('failures')
    assert report['known_native_logposterior_invocations']==32
    assert report['unknown_count_attempts']==0 and report['cumulative_known_native_logposterior_invocations']==33
    assert report['replay_gate_sha256']==replay_sha
    assert len(report['points'])==32
    assert [r['point'] for r in report['points']]==[p['point'] for p in plan['points']]
    return {'status':'verified_numerical_screen_prerequisites_only',
            'parent_target_identity':receipt['parent_target_identity'],
            'original_screen_status':receipt['original_screen_status'],
            'reviewed_accuracy1_to2_status':receipt['reviewed_screen_status'],
            'accuracy2_to3_status':report['status'],'refinement_plan_identity':plan['identity'],
            'replay_sha256':replay_sha,'refinement_summary_sha256':digest(work/'refine-summary.json'),
            'verified_original_receipt_identity':refinement.original.identity(receipt),
            'physical_execution_authorized':False,'posterior_qualified':False}


if __name__=='__main__':
    raise SystemExit('Preparation-only contracts: no physical execution CLI. Run native_accuracy_validate.py for synthetic validation.')
