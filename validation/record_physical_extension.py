#!/usr/bin/env python3
"""Bind the second physical extension without rewriting the first execution.

Run only after physical_extension_checks.py reports a completed, passing audit.
An honestly recorded failed scientific-identification gate is not a missing run.
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'studies/host_ages'
PREVIOUS_COMMIT = 'ca8cfaabba5db7d438aa8ced24b1e947cf7b9eec'


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def read(path):
    return json.loads(path.read_text())


def inventory(paths):
    return {str(p.relative_to(ROOT)): sha(p) for p in sorted(set(paths))
            if p.is_file() and '__pycache__' not in p.parts}


def main():
    report_path = ROOT/'validation/reports/physical-extension.json'
    audit = read(report_path)
    assert audit['status'] == 'passed' and audit['passed']
    assert audit['checker_sha256'] == sha(ROOT/'validation/physical_extension_checks.py')
    now = datetime.now(timezone.utc).isoformat()
    code = [p for name in ['nebular_dust', 'infrared_resolution', 'infrared_photometry']
            for p in (BASE/'code'/name).rglob('*')]
    code += list((BASE/'code/survey_physics').glob('scaled*'))
    results = [p for name in ['nebular_dust', 'infrared_resolution', 'infrared_photometry']
               for p in (BASE/'results'/name).glob('*.json')]
    results += list((BASE/'results/survey_physics').glob('scaled*.json'))
    results += [BASE/'results/galaxy_validation/scaled-bbc-method-review.json']
    notes = [BASE/'notes'/name for name in ['nebular-dust-results.md',
             'infrared-resolution-results.md', 'survey-bbc-support-results.md',
             'survey-bbc-method-review.md', 'infrared-photometry-results.md']]
    assert all(path.exists() for path in notes)
    documents = [ROOT/'README.md', ROOT/'docs/physical-program-results.md', ROOT/'docs/workflows.md', BASE/'notes/host-transport-results.md', ROOT/'studies/infrared/notes/raisin-signed-cohort.md',
                 BASE/'README.md', ROOT/'validation/README.md', ROOT/'provenance/README.md']+notes
    previous = ROOT/'provenance/physical-program-execution.json'
    campaign = read(BASE/'results/survey_physics/scaled-campaign.json')
    record = {
        'schema': 'cosmology-physical-program-extension-v1',
        'recorded_utc': now, 'base_repository_commit': PREVIOUS_COMMIT,
        'previous_execution': {'path': str(previous.relative_to(ROOT)), 'sha256': sha(previous)},
        'scope': 'Observed nebular attenuation comparisons, empirical infrared source information, shorter-wavelength infrared catalogue recovery and enlarged native survey-correction experiments. No additional cosmological correction is inferred.',
        'branches': {
            'nebular_dust': 'Exact observed line joins, common-cohort brightness comparisons, signed weak-line sensitivities and held-out prediction.',
            'infrared_resolution': 'Recovered empirical beams, actual optical source positions and image sampling; conditional linear information, independent numerical checks and explicit source/beam/pixel sensitivities.',
            'infrared_photometry': 'Public shorter-wavelength infrared catalogue observations and explicit source-association, coverage and measurement limits.',
            'survey_physics_scaled': 'Independent native training shards, fixed scientific BBC gates, literal versus retained-age training targets, map-geometry/nuisance identification and supported paired contrasts.'},
        'sample_boundaries': {
            'nebular_parent': 165, 'valid_signed_Balmer': 163,
            'Balmer_formal_SNR3': 138, 'infrared_hosts': 265,
            'infrared_host_band_configurations': 795,
            'new_native_training_attempts': campaign['attempts_total'],
            'counting_note': 'Galaxy cohorts are inherited subsets, not additional distinct-galaxy discoveries; enlarged training is new and evaluation uses existing independent paired simulations.'},
        'validation': {
            'path': str(report_path.relative_to(ROOT)), 'sha256': sha(report_path),
            'meaning': 'Completed numerical and provenance audit. Native or physical identification failures remain explicitly failed in the component scientific records.'},
        'baseline_cosmology': 'No cosmology kernel, frozen data, likelihood, posterior or baseline summary changed. Prior scientific and validation identities are checked against the ca8cfaab execution record.',
        'historical_label_clarification': 'The older DES deep source record labels all 1906261 catalogue rows galaxies. The current text distinguishes 1127570 KNN_CLASS=1 objects before quality cuts. Existing host and neighbour calculations already select that class; source/result bytes and selected measurements are unchanged.',
        'unresolved_physical_requirements': [
            'Unblinded brightness verification and joint age-dust-metallicity-SFH measurement distributions.',
            'Local progenitor/environment mapping rather than assigning central-fibre properties to the explosion.',
            'Calibrated infrared image covariance, exact coadd beams, source completeness and deblended dust-heating inference.',
            'Supported survey correction across the intended population, with training uncertainty, nuisance identifiability and a correctly defined luminosity target.',
            'Validated follow-up selection, light-curve-surface/classifier retraining and contamination or BEAMS coverage.',
            'An observationally supported residual distance correction before a new joint cosmological measurement.'],
        'code_and_design_sha256': inventory(code),
        'compact_result_sha256': inventory(results),
        'documentation_sha256': inventory(documents),
        'figure_sha256': inventory([ROOT/'docs/figures/infrared-source-separation.png']),
        'validation_sha256': inventory([report_path, ROOT/'validation/physical_extension_checks.py', Path(__file__)]),
        'environment_sha256': inventory([ROOT/'pyproject.toml', ROOT/'uv.lock']),
        'publication_policy': 'Only authored code/designs, scientific documents and compact numerical/provenance records. Downloads, native builds and bulk generated tables remain ignored.',
        'new_physical_cosmology_claim': False, 'author_contact_made': False}
    target = ROOT/'provenance/physical-program-extension.json'
    target.write_text(json.dumps(record, indent=2, allow_nan=False)+'\n')
    edition_path = ROOT/'provenance/edition.json'
    history = ROOT/'provenance/history/edition-before-physical-extension.json'
    if not history.exists():
        history.write_bytes(edition_path.read_bytes())
    edition = read(edition_path)
    edition['updated_utc'] = now
    edition['scope'] = 'Current manuscript, original studies, age-correction comparisons, galaxy/survey physics, observed nebular dust, infrared source information and enlarged native correction experiments. A complete physical cosmology measurement remains unidentified.'
    edition['physical_program_extension_manifest'] = {'path': str(target.relative_to(ROOT)), 'sha256': sha(target), 'scope': record['scope']}
    edition['documentation_sha256'] = inventory([ROOT/'README.md']+list((ROOT/'docs').rglob('*.md'))+list((ROOT/'studies').rglob('*.md'))+[ROOT/'provenance/README.md', ROOT/'validation/README.md'])
    edition['validation_tools_sha256'] = inventory(list((ROOT/'validation').glob('*.py')))
    edition_path.write_text(json.dumps(edition, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'manifest': str(target.relative_to(ROOT)), 'code_and_design_files': len(record['code_and_design_sha256']), 'result_records': len(results), 'new_native_attempts': campaign['attempts_total']}))


if __name__ == '__main__':
    main()
