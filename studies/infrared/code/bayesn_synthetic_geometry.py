"""Inspect preserved synthetic posterior attempts without evaluating a model.

This diagnostic reads optical metadata and finished chains only. Correlations
from a failed attempt are sampler diagnostics, not qualified posterior results.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path

os.environ.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
PACKAGE = Path(__file__).resolve().parent/'bayesn_heldout'
NAMES = ['D', 'AV', 'RV', 'theta']+[f'epsilon_white[{i}]' for i in range(42)]+['tau']


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def quantiles(values):
    values = np.asarray(values).ravel()
    finite = values[np.isfinite(values)]
    names = ['min', 'q025', 'q16', 'median', 'q84', 'q975', 'max']
    result = dict(zip(names, map(float, np.quantile(finite, [0,.025,.16,.5,.84,.975,1])))) if len(finite) else dict.fromkeys(names)
    return result|dict(total_count=len(values), finite_count=len(finite),
                       status='finite' if len(finite)==len(values)>0 else 'empty_or_nonfinite_diagnostic_subset')


def finite_or_null(value):
    return float(value) if np.isfinite(value) else None


def energy_bfmi(energy):
    variance = np.var(energy,ddof=1)
    return finite_or_null(np.mean(np.diff(energy)**2)/variance) if variance>0 else None


def geometry(matrix, names):
    if not np.isfinite(matrix).all():
        return dict(status='nonfinite_coordinates', condition_number=None,
                    explanation='No rows removed or interpreted as a qualified posterior.')
    standard_deviation = matrix.std(axis=0, ddof=1)
    active = (standard_deviation>0)&(np.ptp(matrix,axis=0)>0)
    constants = [name for name,keep in zip(names,active) if not keep]
    names = [name for name,keep in zip(names,active) if keep]
    if len(names)<2:
        return dict(status='fewer_than_two_varying_coordinates', condition_number=None,
                    constant_coordinates=constants, active_parameter_order=names)
    matrix = matrix[:,active]; standard_deviation = standard_deviation[active]
    standardized = (matrix-matrix.mean(axis=0))/standard_deviation
    correlation = standardized.T@standardized/(len(matrix)-1)
    eigenvalues, eigenvectors = np.linalg.eigh(correlation)
    pairs = [(abs(correlation[i,j]), i, j) for i in range(len(names)) for j in range(i)]
    pairs.sort(reverse=True)
    cutoff = np.finfo(float).eps*len(names)*eigenvalues[-1]
    full_rank = eigenvalues[0]>cutoff
    return dict(status='full_rank' if full_rank and not constants else 'constant_or_rank_deficient',
                active_parameter_order=names, constant_coordinates=constants,
                correlation=correlation.tolist(), eigenvalues=eigenvalues.tolist(),
                numerical_rank=int(np.sum(eigenvalues>cutoff)), rank_cutoff=float(cutoff),
                condition_number=float(eigenvalues[-1]/eigenvalues[0]) if full_rank else None,
                top_absolute_correlations=[dict(first=names[i], second=names[j], correlation=float(correlation[i,j]))
                                           for _,i,j in pairs[:15]],
                narrowest_standardized_direction=dict(zip(names,map(float,eigenvectors[:,0]))))


def analyze(work, case_index, external_outcome):
    prepared_path = work/'prepared.json'
    prepared = json.loads(prepared_path.read_text())
    current = {path.name:sha(path) for path in PACKAGE.iterdir() if path.suffix in ('.py','.json')}
    assert prepared['source_sha256']==current
    case = work/f'case-{case_index:02d}'
    optical_path = case/'optical.json'
    assert sha(optical_path)==prepared['cases'][case_index]['optical_sha256']
    optical = json.loads(optical_path.read_text())
    assert optical['scope']=='synthetic' and all(row['role']=='optical' for row in optical['rows'])
    times = np.asarray([row['trigger_rest_time'] for row in optical['rows']])
    files = {str(prepared_path):sha(prepared_path), str(optical_path):sha(optical_path)}
    arms = {}
    for arm in ('LCDM', 'broad'):
        folder = case/arm
        chain_records, arrays = [], []
        for record_path in sorted(folder.glob('chain-*.json')):
            record = json.loads(record_path.read_text())
            path = folder/record['archive']
            assert sha(path)==record['sha256']
            files[str(record_path)] = sha(record_path); files[str(path)] = sha(path)
            with np.load(path) as raw:
                parameters = np.column_stack([raw[name] for name in ('D','AV','RV','theta','epsilon_white','tau')])
                assert parameters.shape==(1000,47) and record['warmup']==record['draws']==1000
                energy = raw['stat_energy']; potential = raw['stat_potential_energy']
                divergent = raw['stat_diverging'].astype(bool)
                assert int(divergent.sum())==record['divergences']
                phase = times[None,:]-parameters[:,-1,None]
                distance_to_integer = np.min(abs(phase-np.round(phase)),axis=1)
                divergence_states = []
                for index in np.flatnonzero(divergent):
                    divergence_states.append(dict(draw=int(index), parameters=[finite_or_null(v) for v in parameters[index]],
                        energy=finite_or_null(energy[index]), potential_energy=finite_or_null(potential[index]),
                        num_steps=int(raw['stat_num_steps'][index]),
                        nearest_Hsiao_integer_phase_distance=finite_or_null(distance_to_integer[index])))
                chain_records.append(record|dict(
                    all_retained_parameters_finite=bool(np.isfinite(parameters).all()),
                    EBFMI=energy_bfmi(energy),
                    energy=quantiles(energy), potential_energy=quantiles(potential),
                    absolute_adjacent_energy_difference=quantiles(abs(np.diff(energy))),
                    num_steps=quantiles(raw['stat_num_steps']),
                    mean_accept_probability=finite_or_null(raw['stat_accept_prob'].mean()),
                    maximum_tree_fraction=float(np.mean(raw['stat_num_steps']>=1023)),
                    divergence_associated_retained_states=divergence_states,
                    nondivergent_nearest_integer_phase=quantiles(distance_to_integer[~divergent]),
                    nondivergent_within_0_001_rest_days_of_integer_phase=int(np.sum((distance_to_integer<.001)&~divergent)),
                    physical_coordinate_quantiles={NAMES[index]:quantiles(parameters[:,index]) for index in (0,1,2,3,46)}))
                arrays.append(parameters)
        entry = dict(completed_chains=len(arrays), chains=chain_records,
                     complete_four_chain_attempt=len(arrays)==4,
                     divergence_gate_passed=bool(arrays) and sum(row['divergences'] for row in chain_records)==0)
        manifest_path = folder/'posterior-manifest.json'
        if manifest_path.exists():
            manifest = json.loads(manifest_path.read_text())
            assert manifest['worker_spec']['source_sha256']==current
            assert manifest['worker_spec']['optical_sha256']==sha(optical_path)
            assert len(manifest['chains'])==len(chain_records)==4
            for recorded, analyzed in zip(manifest['chains'],chain_records):
                assert all(recorded[key]==analyzed[key] for key in recorded)
            files[str(manifest_path)] = sha(manifest_path)
            entry['worker_wall_seconds'] = manifest['seconds']
        diagnostic_path = folder/'diagnostics.json'
        if diagnostic_path.exists():
            files[str(diagnostic_path)] = sha(diagnostic_path)
            entry['official_optical_diagnostics'] = json.loads(diagnostic_path.read_text())
            assert entry['official_optical_diagnostics']['posterior_manifest_sha256']==sha(manifest_path)
            assert entry['official_optical_diagnostics']['divergences']==sum(row['divergences'] for row in chain_records)
        if arrays:
            pooled = np.concatenate(arrays)
            transformed = pooled.copy()
            with np.errstate(divide='ignore',invalid='ignore'):
                transformed[:,1] = np.log(pooled[:,1])
                transformed[:,2] = np.log((pooled[:,2]-1.2)/(6.-pooled[:,2]))
                transformed[:,46] = np.log((pooled[:,46]+10.)/(20.-pooled[:,46]))
            transformed_names = NAMES.copy()
            transformed_names[1]='log_AV'; transformed_names[2]='logit_RV'; transformed_names[46]='logit_tau'
            entry['physical_geometry'] = geometry(pooled, NAMES)
            entry['unconstrained_geometry'] = geometry(transformed, transformed_names)
            within = np.concatenate([array-array.mean(axis=0) for array in arrays])
            entry['within_chain_physical_geometry'] = geometry(within, NAMES)
        arms[arm] = entry
    return dict(schema='bayesn-synthetic-geometry-v1', case_index=case_index,
                external_execution_outcome=external_outcome,
                parameter_order=NAMES, arms=arms,
                qualified_scientific_posterior=False,
                source_sha256={str(Path(__file__).resolve().relative_to(ROOT)):sha(__file__)},
                reviewed_package_sha256=current, input_sha256=files,
                caveats=[
                    'This is a diagnostic of preserved chains, including failed or incomplete attempts; no posterior inference or NIR score.',
                    'A divergence flag labels a retained transition state, not the actual failing leapfrog location. Proximity to a phase knot cannot establish a cause.',
                    'Hsiao interpolation is piecewise linear at integer rest phases. The recorded distances identify a diagnostic lead only.',
                    'The final adapted metric, step size and full divergent trajectories were not archived; their exact Hamiltonian error cannot be reconstructed.',
                    'Unconstrained coordinates use the same positive/interval transforms as the declared sampler; geometry does not change its target or qualify a replacement.'])


def compact(result, full_path):
    arms = {}
    incomplete = False
    failed = False
    for arm, entry in result['arms'].items():
        diagnostic = entry.get('official_optical_diagnostics')
        incomplete |= diagnostic is None
        failed |= diagnostic is not None and not diagnostic['sampler_passed']
        states = [state for chain in entry['chains'] for state in chain['divergence_associated_retained_states']]
        parameters = np.asarray([state['parameters'] for state in states],dtype=float)
        arms[arm] = dict(completed_chains=entry['completed_chains'],
                        worker_wall_seconds=entry.get('worker_wall_seconds'),
                        chain_wall_seconds=[chain['seconds'] for chain in entry['chains']],
                        divergences_by_chain=[chain['divergences'] for chain in entry['chains']],
                        mean_accept_probability_by_chain=[chain['mean_accept_probability'] for chain in entry['chains']],
                        median_leapfrog_steps_by_chain=[chain['num_steps']['median'] for chain in entry['chains']],
                        diagnostics={key:value for key,value in (diagnostic or {}).items() if key!='scalar_parameters'},
                        geometry={label:{key:entry[label].get(key) for key in ('status','condition_number','top_absolute_correlations')}
                                  for label in ('physical_geometry','unconstrained_geometry','within_chain_physical_geometry') if label in entry})
        if states:
            valid_phases = [state['nearest_Hsiao_integer_phase_distance'] for state in states
                            if state['nearest_Hsiao_integer_phase_distance'] is not None]
            arms[arm]['divergence_retained_states'] = dict(
                count=len(states), unique_coordinate_vectors=len({tuple(row) for row in parameters}) if np.isfinite(parameters).all() else None,
                physical_coordinate_ranges={NAMES[index]:[finite_or_null(np.min(parameters[:,index])),finite_or_null(np.max(parameters[:,index]))]
                                            for index in (0,1,2,3,46)},
                within_0_001_rest_days_of_any_integer_phase=sum(value<.001 for value in valid_phases),
                nondivergent_within_same_phase_threshold=sum(chain['nondivergent_within_0_001_rest_days_of_integer_phase'] for chain in entry['chains']),
                nondivergent_states=sum(chain['draws']-chain['divergences'] for chain in entry['chains']),
                phase_threshold_scope='Exploratory description chosen after the failed sampling attempt; repeated states are not independent trials.',
                interpretation='Retained states flagged by a divergent transition, not actual failed leapfrog states; no causal knot diagnosis.')
    return dict(schema='bayesn-synthetic-case-diagnostic-v1', case_index=result['case_index'],
                status='incomplete_attempt' if incomplete else 'failed_sampling_gate' if failed else 'sampling_gates_passed_prediction_not_assessed',
                external_execution_outcome=result['external_execution_outcome'], arms=arms,
                scope='Synthetic-only numerical sampling diagnostic; no qualified posterior or infrared predictive comparison.',
                full_geometry_path=str(full_path.relative_to(ROOT)), full_geometry_sha256=sha(full_path),
                source_sha256=result['source_sha256'], reviewed_package_sha256=result['reviewed_package_sha256'],
                input_sha256={str(Path(name).relative_to(ROOT)):value for name,value in result['input_sha256'].items()},
                no_model_or_observed_photometry_calls=True, caveats=result['caveats'])


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--case', type=int, default=0)
    parser.add_argument('--external-outcome', choices=('completed','timeout','failed'), required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--summary-output', type=Path)
    args = parser.parse_args()
    assert not args.output.exists(), 'Preserve previous diagnostic records.'
    result = analyze(args.work.resolve(),args.case,args.external_outcome)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if args.summary_output:
        assert not args.summary_output.exists(), 'Preserve previous compact records.'
        args.summary_output.parent.mkdir(parents=True,exist_ok=True)
        args.summary_output.write_text(json.dumps(compact(result,args.output.resolve()),indent=2,allow_nan=False)+'\n')
    print(json.dumps({arm:{'completed_chains':entry['completed_chains'],
                           'divergences':sum(row['divergences'] for row in entry['chains'])}
                      for arm,entry in result['arms'].items()},indent=2))
