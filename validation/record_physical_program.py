#!/usr/bin/env python3
"""Record the completed physical experiments without rewriting past executions.

Run after the component analyses and validation/physical_program_checks.py.
The earlier distance likelihood and its frozen numerical outputs are unchanged.
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'studies/host_ages'
TOPICS = ['galaxy_validation', 'host_transport', 'physical_ages', 'survey_physics']


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def read(path):
    return json.loads(path.read_text())


def inventory(paths):
    return {str(p.relative_to(ROOT)): sha(p) for p in sorted(set(paths))}


def main():
    now = datetime.now(timezone.utc).isoformat()
    report = ROOT / 'validation/reports/physical-program.json'
    assert read(report)['status'] == 'passed', 'Run the physical programme audit first'
    previous = read(ROOT / 'validation/manifest.json')
    assert previous['input_files_verified'] == 391 and len(previous['runs']) == 29
    assert all(x['summary_matches_pre_reorganization'] for x in previous['baseline_comparisons'].values())
    sources = [p for t in TOPICS for p in (BASE/'code'/t).rglob('*')
               if p.is_file() and '__pycache__' not in p.parts]
    results = [p for t in TOPICS for p in (BASE/'results'/t).glob('*.json')]
    docs = [ROOT/'README.md', ROOT/'docs/physical-program-results.md',
            ROOT/'docs/age-correction-plan.md', ROOT/'docs/age-correction-results.md',
            ROOT/'validation/README.md', ROOT/'provenance/README.md', BASE/'README.md']
    docs += [BASE/'notes'/n for n in ['galaxy-validation-results.md',
            'host-transport-results.md', 'physical-ages-results.md', 'survey-physics-results.md']]
    native = read(BASE/'results/survey_physics/validation.json')
    bounds = {n: read(BASE/f'results/physical_ages/{n}-age-bounds.json') for n in ['sdss', 'roman', 'des']}
    record = {
        'schema': 'cosmology-physical-program-execution-v1', 'recorded_utc': now,
        'base_repository_commit': '93c2026016120e0a2e4ce9305a3b315760f86746',
        'manuscript': 'docs/physical-program-results.md',
        'question': 'What observed galaxy light and physically propagated survey interventions establish about residual age effects after standardization.',
        'scope': 'Executed observational and conditional physical tests. No complete selected-population luminosity model or new corrected cosmological measurement is identified.',
        'branches': {
            'galaxy_validation': 'Recovered and independently integrated spectra, spatial environments, same-object age comparisons and incremental brightness predictions.',
            'host_transport': 'Recovered local and high-redshift photometry and SPIRE images; tested calibration, association, signed measurements, confusion and cross-survey transport.',
            'physical_ages': 'Physical nonnegative stellar mixtures, formed-mass age bounds, held-out and joint spectral-index checks, numerical recovery and CUDA comparison.',
            'survey_physics': 'Native photon generation, selection and fitting, reconstructed classifier signal acceptance, separate frozen/refitted corrections and native BBC support tests.'},
        'physical_flux_samples': {n: d['objects'] for n, d in bounds.items()},
        'counting_boundary': 'The three flux samples can overlap and measure different apertures; their sum is not a distinct-galaxy count.',
        'native_generated_attempts': native['generated_attempts'],
        'native_scientific_jobs': native['native_run_count'],
        'corrections_during_validation': [
            'Retained all finite signed Roman local measurements, including five formerly excluded nondetections.',
            'Replaced insufficient DES observer-grid quadrature by Gauss integration over both spectral and passband knots; reran all physical fits.',
            'Replaced defective SciPy 1.15.3 NNLS on the new stellar design with independently checked column-scaled BVLS.',
            'Required original-unit feasibility and duality checks for every conic retry candidate.',
            'Kept historical solver and survey engineering failures as explicitly scoped regression records.'],
        'unresolved_gates': [
            'Absolute/progenitor ages and realistic joint age-dust-metallicity-SFH likelihoods are not uniquely measured.',
            'Fibre and local environments differ; foreground, photometric and spectral covariance remain partly unavailable.',
            'Far-infrared host luminosities require source deblending and dust-heating assumptions.',
            'Host follow-up/measurement availability is not a validated cosmic-population selection function.',
            'Production-dimensional BBC lacks interpolation support; original LFS mocks remain inaccessible.',
            'Light-curve surface/classifier retraining, contamination/BEAMS and full survey coverage are not established.',
            'No observationally validated residual distance correction is available for a new cosmological likelihood.'],
        'protocol_history': 'Initial choices and dated amendments remain in component designs. Exploratory extensions and engineering repairs are identified; none is retroactively called blinded.',
        'code_and_design_sha256': inventory(sources),
        'compact_result_sha256': inventory(results),
        'documentation_sha256': inventory(docs),
        'figure_sha256': inventory([ROOT/'docs/figures/host-age-physics.png']),
        'environment_sha256': inventory([ROOT/'pyproject.toml', ROOT/'uv.lock']),
        'validation_sha256': inventory([report, ROOT/'validation/manifest.json',
            ROOT/'validation/reports/age-correction-execution.json',
            ROOT/'validation/physical_program_checks.py', Path(__file__)]),
        'baseline_cosmology': 'All 391 frozen input identities and 29 reference run records are reverified, including 19 unchanged default summaries. These numerical distance runs were not rerun in this physical campaign.',
        'publication_policy': 'First-party code, design definitions, scientific notes, compact numerical/provenance JSON and manuscript figures. Downloaded libraries/data, generated per-object tables, chains and native builds remain local.',
        'author_contact_made': False, 'new_physical_cosmology_claim': False}
    target = ROOT/'provenance/physical-program-execution.json'
    target.write_text(json.dumps(record, indent=2)+'\n')
    edition_path = ROOT/'provenance/edition.json'
    history = ROOT/'provenance/history/edition-before-physical-program.json'
    if not history.exists():
        history.write_bytes(edition_path.read_bytes())
    edition = read(edition_path)
    edition['updated_utc'] = now
    edition['scope'] = 'Current manuscript, original studies, age-correction comparisons and galaxy/survey physical extension. Successful numerical checks do not establish all physical hypotheses or a complete unified cosmology measurement.'
    edition['documentation_sha256'] = inventory([ROOT/'README.md'] + list((ROOT/'docs').rglob('*.md')) + list((ROOT/'studies').rglob('*.md')) + [ROOT/'provenance/README.md', ROOT/'validation/README.md'])
    edition['physical_program_execution_manifest'] = {'path': str(target.relative_to(ROOT)), 'sha256': sha(target), 'scope': record['scope']}
    edition['validation_tools_sha256'] = inventory(list((ROOT/'validation').glob('*.py')))
    edition_path.write_text(json.dumps(edition, indent=2)+'\n')
    print(json.dumps({'manifest': str(target.relative_to(ROOT)), 'code_and_design_files': len(sources), 'result_records': len(results), 'physical_samples': record['physical_flux_samples'], 'native_attempts': record['native_generated_attempts']}))


if __name__ == '__main__':
    main()
