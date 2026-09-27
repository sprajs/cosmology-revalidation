"""Read-only aggregate identity and exact-spectrum acquisition audit.

Run in the modern environment after the individual numerical validators.
The only output written is aggregate-validation.json. Qualified historical
outcomes are preserved; this is not a posterior or emulator validation.
"""
import ast
from collections import Counter
import importlib
import json
from pathlib import Path
import subprocess

import numpy as np

from acquire import HERE, ROOT, WORK, PACKAGES, RESULTS, sha, tree
from modern_acquire import REPOS


PRODUCERS = {
    'acquisition': 'acquire', 'author-acquisition': 'author_acquire',
    'author-configuration': 'author_config',
    'author-density-check': 'author_density_check',
    'author-density-reconstruction': 'author_evaluate',
    'modern-acquisition': 'modern_acquire', 'modern-proposal': 'modern_proposal',
    'modern-release-checks': 'modern_release_checks',
    'modern-restoration': 'restore_modern', 'proposal': 'proposal',
    'sn-proposal': 'sn_proposal', 'spectral-training': 'spectral_training',
    'table-cpl-sensitivity': 'table_cpl_sensitivity',
    'theory-domain': 'theory_domain', 'theory-trim-check': 'theory_trim_check',
}


def main():
    records = {p.stem: json.loads(p.read_text()) for p in RESULTS.glob('*.json')
               if p.stem != 'aggregate-validation'}
    initial = {str(p.relative_to(ROOT)): sha(p) for p in sorted(HERE.glob('*'))
               if p.is_file()}
    for p in HERE.glob('*.py'):
        ast.parse(p.read_text(), filename=str(p))
    source_checks = 0
    for name, record in records.items():
        identity = record.get('code_sha256')
        if isinstance(identity, str):
            assert sha(HERE/(PRODUCERS[name]+'.py')) == identity, name
            source_checks += 1
        elif isinstance(identity, dict):
            for filename, expected in identity.items():
                assert sha(HERE/filename) == expected, (name, filename)
                source_checks += 1
        if isinstance(record.get('source_sha256'), dict):
            for filename, expected in record['source_sha256'].items():
                assert sha(ROOT/filename) == expected, (name, filename)
                source_checks += 1

    asset_checks = {}
    acquisition = records['acquisition']
    assert sha(HERE/'requirements-lock.txt') == acquisition['requirements_lock_sha256']
    assert sha(ROOT/acquisition['file_inventory']) == acquisition['file_inventory_sha256']
    for name, expected in acquisition['assets'].items():
        path = (PACKAGES/'code/planck/clipy/clipy' if name == 'clipy_source'
                else PACKAGES/'data'/name)
        actual, _ = tree(path)
        assert actual == {k: expected[k] for k in actual}, name
        asset_checks['baseline/'+name] = actual
    modern = records['modern-acquisition']
    assert sha(HERE/'modern-requirements-lock.txt') == modern['requirements_lock_sha256']
    assert sha(ROOT/modern['inventory_path']) == modern['inventory_sha256']
    paths = {'ACTDR6_CMB': PACKAGES/'data/ACTDR6CMBonly',
             'ACT_Planck_lensing': PACKAGES/'data/ACT_dr6_likelihood/v1.2'}
    paths.update({name: Path(importlib.import_module(name).__file__).parent
                  for name in REPOS})
    for name, path in paths.items():
        actual, _ = tree(path)
        assert actual == {k: modern['assets'][name][k] for k in actual}, name
        asset_checks['modern/'+name] = actual
    author = records['author-acquisition']
    for expected in author['assets'].values():
        path = ROOT/expected['path']
        assert sha(path) == expected['sha256'] and path.stat().st_size == expected['bytes']
    actual, _ = tree(WORK/'author-chains/act-lensing/v1.1')
    assert actual == author['ACT_lensing_extracted_tree']
    assert sha(ROOT/author['lensing_inventory_path']) == author['lensing_inventory_sha256']
    asset_checks['author/ACT_Planck_lensing_v1.1'] = actual

    dependency_checks = {}
    for env, locks in [('.venv', ['requirements-lock.txt']),
                       ('.modern-venv', ['modern-requirements-lock.txt', 'sampling-runtime-lock.txt'])]:
        # Ask each environment itself: no mutation and no scientific imports.
        probe = ('import importlib.metadata as m,json;'
                 'print(json.dumps({d.metadata["Name"].lower().replace("_","-"): '
                 'd.version for d in m.distributions()}))')
        installed = json.loads(subprocess.check_output(
            [str(WORK/env/'bin/python'), '-c', probe], text=True))
        checked = {}
        for lock in locks:
            for line in (HERE/lock).read_text().splitlines():
                if '==' not in line:
                    continue  # VCS source identity is checked by the trees above.
                name, version = line.split('==')
                assert installed[name.lower().replace('_','-')] == version, (env, name)
                checked[name] = version
        dependency_checks[env] = checked

    def link(record, field, filename):
        assert records[record][field] == sha(RESULTS/filename), (record, field)

    link('validation', 'acquisition_sha256', 'acquisition.json')
    link('modern-restoration', 'frozen_acquisition_sha256', 'modern-acquisition.json')
    link('fast-lensing-validation', 'acquisition_sha256', 'modern-acquisition.json')
    link('author-density-check', 'source_result_sha256', 'author-density-reconstruction.json')
    for filename, expected in records['author-density-check']['dependency_records'].items():
        assert sha(RESULTS/filename) == expected

    design = json.loads((HERE/'spectral-training-design.json').read_text())
    training = records['spectral-training']
    assert training['status'] == 'complete_exact_acquisition'
    assert training['design_sha256'] == sha(HERE/'spectral-training-design.json')
    assert design['modern_acquisition_sha256'] == sha(RESULTS/'modern-acquisition.json')
    assert design['adapter_sha256'] == sha(HERE/'modern_adapter.py')
    assert design['proposal_sha256'] == sha(WORK/'proposal-cpl-dovekie.covmat')
    assert design['proposal_transform_sha256'] == sha(RESULTS/'proposal.json')
    centre = np.array(design['centre']); chol = np.array(design['coordinate_cholesky'])
    np.testing.assert_allclose(chol@chol.T, design['coordinate_covariance'], rtol=1e-12, atol=1e-18)
    expected_points = {}
    for split, n in design['counts'].items():
        rng = np.random.default_rng(design['seeds'][split])
        draws = rng.standard_normal((n, 8))
        widths = np.where(rng.random(n)<0.1, 2.5, 1.5)
        expected_points[split] = centre + (draws*widths[:,None])@chol.T
    counts = {key: Counter() for key in design['counts']}
    exclusions = []; seen = set(); max_angle_difference = 0.
    for row in training['rows']:
        split, index = row['split'], row['index']
        assert (split,index) not in seen
        seen.add((split,index))
        path = WORK/'spectral-training'/split/f'{index:04d}.json'
        assert sha(path) == row['manifest_sha256']
        rec = json.loads(path.read_text()); assert rec['status'] == row['status']
        assert rec['adapter_sha256'] == design['adapter_sha256']
        assert rec['training_code_sha256'] == training['code_sha256']
        np.testing.assert_array_equal(rec['requested_coordinates'], expected_points[split][index])
        counts[split][rec['status']] += 1
        if rec['status'] != 'finite_exact':
            assert rec['status'] == 'nonfinite_prior_theory_or_likelihood'
            assert rec['physical_point']['tau'] < .01
            assert not path.with_suffix('.npz').exists()
            exclusions.append({'split':split,'index':index,'tau':rec['physical_point']['tau'],
                               'reason':'below_declared_tau_prior'})
            continue
        spectrum = ROOT/rec['file']; assert sha(spectrum) == rec['sha256']
        with np.load(spectrum) as data:
            assert data['spectra'].shape == (5,10152)
            assert np.isfinite(data['spectra']).all() and np.isfinite(data['loglikes']).all()
            np.testing.assert_array_equal(data['ell'], np.arange(10152))
            np.testing.assert_array_equal(data['coordinates'], rec['actual_coordinates'])
            np.testing.assert_array_equal(data['requested_coordinates'], rec['requested_coordinates'])
            np.testing.assert_array_equal(data['loglikes'], rec['loglikes'])
            np.testing.assert_array_equal(data['physical_point'], [rec['physical_point'][k] for k in design['physical_names']])
            assert data['derived_rdrag'] > 0
            max_angle_difference = max(max_angle_difference, abs(float(data['coordinates'][0]-data['requested_coordinates'][0])))
        for key, expected in training['theory_metadata'].items():
            assert rec[key] == expected, (split,index,key)
    assert seen == {(split,i) for split,n in design['counts'].items() for i in range(n)}
    assert counts['train'] == {'finite_exact':509,'nonfinite_prior_theory_or_likelihood':3}
    assert counts['holdout'] == {'finite_exact':96}
    assert training['theory_metadata']['CAMB_Params_max_l'] == 10251
    assert training['theory_metadata']['finalized_theory_extra_args']['lmax'] == 9001

    for name in ['validation','modern-validation','modern-release-checks',
                 'modern-restoration','fast-lensing-validation','fast-lensing-independent',
                 'sn-proposal']:
        assert records[name]['status'] == 'passed', name
    assert records['theory-trim-check']['status'] == 'completed_not_adopted'
    assert records['theory-trim-check']['target_unchanged']
    assert not records['author-density-check']['exact_author_runtime_recovered']
    assert not records['author-density-check']['can_claim_independent_cosmological_measurement']
    final = {str(p.relative_to(ROOT)): sha(p) for p in sorted(HERE.glob('*')) if p.is_file()}
    assert final == initial, 'Source changed during the read-only audit.'
    output = {
        'status':'passed', 'scope':'Identity, exact-acquisition integrity and previously executed numerical checks only; not posterior convergence, emulator validation or missing-covariance certification.',
        'source_sha256':final, 'result_sha256':{str((RESULTS/(k+'.json')).relative_to(ROOT)):sha(RESULTS/(k+'.json')) for k in sorted(records)},
        'producer_source_references_verified':source_checks,
        'verified_asset_trees':asset_checks, 'author_single_files_verified':len(author['assets']),
        'dependency_versions_verified':dependency_checks,
        'exact_spectral_acquisition':{'attempts':len(seen),'counts':{k:dict(v) for k,v in counts.items()},
            'prior_exclusions':exclusions,'maximum_requested_actual_acoustic_angle_difference':max_angle_difference,
            'theory_metadata':training['theory_metadata'],'all_row_json_and_npz_hashes_verified':True,
            'all_seeded_design_rows_reconstructed':True},
        'qualified_outcomes_retained':{name:records[name].get('status',records[name].get('effective_primary_CPL_domain')) for name in ['author-configuration','author-density-reconstruction','author-density-check','table-cpl-sensitivity','theory-domain','theory-trim-check']},
        'remaining_scientific_limits':['Analytic CAMB CPL domain w0+wa<=0.',
            'Missing primary/lensing and inter-experiment sampling cross-covariance beyond released ACT/Planck lensing covariance.',
            'Author runtime not byte-identical: five-point conditional density reconstruction only.',
            'Numerical accuracy sensitivity is measured at reference, not fully across posterior.',
            'Spectral emulator held-out validation and exact posterior correction belong to the inference lane.'],
    }
    (RESULTS/'aggregate-validation.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({'status':'passed','attempts':len(seen),'counts':output['exact_spectral_acquisition']['counts'],
                      'asset_trees':len(asset_checks),'source_references':source_checks}))


if __name__ == '__main__':
    main()
