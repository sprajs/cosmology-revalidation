"""Diagnose proposal fallback costs without calculating native CMB spectra.

Use --chain-folder to snapshot completed chain rows and the current covariance,
or --frozen-design to replay the numerical design embedded in a previous compact
result. Only the requested ignored output directory and compact summary change.
This is a fixed-centre sensitivity study, not a Markov chain or ESS estimate.
"""
import argparse
from collections import Counter
from contextlib import contextmanager
import importlib.metadata
import json
from pathlib import Path
import time

import camb
import numpy as np
from cobaya.model import get_model
from cobaya.samplers.mcmc.proposal import BlockedProposer
import cobaya.samplers.mcmc.proposal as native_proposal

from modern_fast import configuration
from spectral_surrogate import coordinates, polynomial, SPECTRA, COORDINATES
from target_identity import ROOT, digest


LIMITS = [
    'Fixed accepted centres cycle without accept/reject updates; rates are not actual chain telemetry.',
    'The native radial distribution and random rotations are used, with paired displacements across scales.',
    'The small number of tail proposals does not establish a precise fallback probability.',
    'Hypothetical envelope counts do not certify interpolation accuracy outside the adopted envelope.',
    'The 30-second native-fallback cost is an explicit scenario, not measured by this diagnostic.',
    'Acceptance probabilities, accepted jump distances and ESS per time are not measured.',
    'No active chain, target, interpolator or proposal setting is changed.',
]


def relative(path):
    path = Path(path).resolve()
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2)+'\n')


def capture_design(args):
    if args.frozen_design:
        record = json.loads(args.frozen_design.read_text())
        design = record.get('diagnostic_design', record)
        assert design['model_sha256'] == digest(args.model_file)
        return design, {'kind':'frozen_numerical_design',
                        'path':relative(args.frozen_design), 'sha256':digest(args.frozen_design)}
    snapshots = args.output/'snapshots'
    snapshots.mkdir()
    covariance_file = args.chain_folder/'chain.covmat'
    covariance_bytes = covariance_file.read_bytes()
    (snapshots/'chain.covmat').write_bytes(covariance_bytes)
    names = covariance_bytes.decode().splitlines()[0].lstrip('# ').split()
    covariance = np.loadtxt(snapshots/'chain.covmat')
    inputs = {relative(covariance_file):digest(snapshots/'chain.covmat')}
    row_snapshots = {}; centres = []
    for path in sorted(args.chain_folder.glob('chain.*.txt')):
        raw = path.read_bytes(); lines = raw.decode().splitlines()
        columns = lines[0].lstrip('# ').split()
        # A running process may have an incomplete final line. Never retain it.
        complete = [line for line in lines[1:] if len(line.split()) == len(columns)]
        if not complete:
            continue
        snapshot = snapshots/path.name
        snapshot.write_text(lines[0]+'\n'+'\n'.join(complete)+'\n')
        import hashlib
        inputs[relative(path)] = hashlib.sha256(raw).hexdigest()
        row_snapshots[path.name] = {'sha256':digest(snapshot), 'rows':len(complete),
                                   'omitted_incomplete_rows':len(lines)-1-len(complete)}
        data = np.atleast_2d(np.loadtxt(snapshot))
        for index in np.unique(np.floor(np.array([.25,.75])*(len(data)-1)).astype(int)):
            centres.append({'chain':path.name,'row':int(index),
                'point':{key:float(data[index,columns.index(key)]) for key in names}})
    if len(centres) < 4:
        raise ValueError('At least four distinct accepted-row centres are required.')
    design = {
        'seed':args.seed, 'proposals_per_scale':args.draws, 'scales':args.scales,
        'centres':centres, 'parameter_names':names, 'proposal_covariance':covariance.tolist(),
        'covariance_sha256':digest(snapshots/'chain.covmat'),
        'model_sha256':digest(args.model_file), 'source_input_sha256':inputs,
        'chain_snapshots':row_snapshots,
        'target':{'model':args.model,'evolution':args.evolution,'sample':args.sample,
                  'calibration':args.calibration},
        'proposal':'Native Cobaya single correlated block, unit-scale random rotation/radial mixture; common displacements across scales; cycle through frozen centres.',
    }
    return design, {'kind':'chain_snapshots','chain_folder':relative(args.chain_folder)}


@contextmanager
def prohibit_native_spectra():
    """Fail closed if the diagnostic accidentally asks CAMB for full spectra."""
    original = camb.get_results
    def forbidden(*args, **kwargs):
        raise RuntimeError('Native spectral calls are forbidden in this diagnostic.')
    camb.get_results = forbidden
    try:
        yield
    finally:
        camb.get_results = original


def classify(model, design):
    theory = model.theory['spectral_surrogate']
    extra = dict(theory.extra_args)
    extra['lmax'] = max(extra.get('lmax',0), max(theory.required_cls.values(),default=0))
    names = design['parameter_names']
    assert list(model.parameterization.sampled_params()) == names
    proposer = BlockedProposer([list(range(len(names)))],
        np.random.default_rng(design['seed']),oversampling_factors=[1],proposal_scale=1.)
    proposer.set_covariance(np.array(design['proposal_covariance']))
    displacements = []
    for _ in range(design['proposals_per_scale']):
        delta = np.zeros(len(names)); proposer.get_proposal(delta)
        displacements.append(delta)
    records = []
    for scale in design['scales']:
        for index, delta in enumerate(displacements):
            centre_index = index % len(design['centres'])
            centre = design['centres'][centre_index]['point']
            point = {key:float(centre[key]+scale*delta[j]) for j,key in enumerate(names)}
            record = {'scale':scale,'draw':index,'centre':centre_index,'point':point}
            if not np.isfinite(model.logpriors(point)).all():
                record['category'] = 'prior_reject'; records.append(record); continue
            inputs = model.parameterization.to_input(point)
            pars_input = {key:inputs[key] for key in theory.input_params}
            start = time.perf_counter()
            try:
                parameters = camb.set_params(**pars_input, **extra)
                background = camb.get_background(parameters)
            except camb.CAMBError as error:
                record.update(category='background_failure',error=str(error))
                records.append(record); continue
            record['background_seconds'] = time.perf_counter()-start
            coord = coordinates(parameters,background)
            white = np.linalg.solve(theory.chol,coord-theory.centre)
            record.update(max_whitened=float(max(abs(white))),
                          max_coordinate=COORDINATES[int(np.argmax(abs(white)))],
                          whitened=white.tolist())
            start = time.perf_counter()
            spectra = ((polynomial(white,theory.exponents)@theory.coefficients)
                       *theory.output_scale).reshape(5,theory.length)
            record['polynomial_seconds'] = time.perf_counter()-start
            record['nonfinite_spectra'] = bool(not np.isfinite(spectra).all())
            record['negative_spectra'] = [key for j,key in enumerate(SPECTRA)
                if key in ['tt','ee','pp'] and np.any(spectra[j,2:]<=0)]
            record['category'] = ('envelope_fallback' if record['max_whitened']>theory.envelope
                else 'spectrum_fallback' if record['nonfinite_spectra'] or record['negative_spectra']
                else 'fast')
            records.append(record)
        print(json.dumps({'scale':scale,'counts':dict(Counter(
            r['category'] for r in records if r['scale']==scale))}),flush=True)
    return records


def benchmark(model, records, scale):
    """Local in-memory instrumentation; restore methods even on failure."""
    originals = []; current = {}; measurements = []
    try:
        for name, component in list(model.theory.items())+list(model.likelihood.items()):
            original = component.calculate; originals.append((component,original))
            def wrapped(*args, _original=original, _name=name, **kwargs):
                start = time.perf_counter()
                try:
                    return _original(*args,**kwargs)
                finally:
                    current.setdefault(_name,[]).append(time.perf_counter()-start)
            component.calculate = wrapped
        fast = [r for r in records if r['scale']==scale and r['category']=='fast']
        if len(fast)<9:
            raise ValueError('At least nine fast proposals are required for timing.')
        for index in np.linspace(0,len(fast)-1,9).astype(int):
            current = {}; start = time.perf_counter()
            result = model.logposterior(fast[index]['point'],cached=False)
            elapsed = time.perf_counter()-start
            assert np.isfinite(result.logpost)
            measurements.append({'draw':fast[index]['draw'],'seconds':elapsed,
                'logpost':float(result.logpost),
                'component_seconds':{key:sum(value) for key,value in current.items()}})
    finally:
        for component, original in originals:
            component.calculate = original
    return measurements


def summarize(records, design, measurements, envelope):
    median = float(np.median([r['seconds'] for r in measurements[1:]]))
    result = {}
    for scale in design['scales']:
        rows = [r for r in records if r['scale']==scale]
        counts = Counter(r['category'] for r in rows)
        finite = [r for r in rows if 'max_whitened' in r]
        fallbacks = counts['envelope_fallback']+counts['spectrum_fallback']
        item = {'counts':dict(counts),'attempts':len(rows),
            'fallback_fraction_all_proposals':fallbacks/len(rows),
            'fallback_fraction_prior_valid':fallbacks/(len(rows)-counts['prior_reject']),
            'envelope_excursion_coordinate_counts':dict(Counter(
                r['max_coordinate'] for r in finite if r['max_whitened']>envelope)),
            'negative_polynomial_spectra_anywhere_count':sum(
                bool(r['negative_spectra']) or r['nonfinite_spectra'] for r in finite),
            'max_whitened_quantiles':np.quantile(
                [r['max_whitened'] for r in finite],[0,.5,.9,.95,1]).tolist(),
            'hypothetical_envelope_counts':{str(e):sum(
                r['max_whitened']>e for r in finite) for e in [3,4,5,6]},
            # Retain the original conservative simple scenario: even prior
            # rejections are charged a complete fast evaluation here.
            'conditional_mean_seconds_if_fallback_30s':
                (1-fallbacks/len(rows))*median+fallbacks/len(rows)*30,
            'scenario_prior_rejections_charged_as_fast':True,
        }
        result[str(scale)] = item
    return result, median


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument('--chain-folder',type=Path)
    inputs.add_argument('--frozen-design',type=Path)
    parser.add_argument('--model-file',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True,help='New ignored working directory')
    parser.add_argument('--summary',type=Path,default=ROOT/'studies/unified_cosmology/results/inference/proposal-throughput.json')
    parser.add_argument('--historical-result',type=Path)
    parser.add_argument('--seed',type=int,default=2727766)
    parser.add_argument('--draws',type=int,default=192)
    parser.add_argument('--scales',type=float,nargs='+',default=[2.4,1.6,1.2])
    parser.add_argument('--model',choices=['lcdm','cpl'],default='cpl')
    parser.add_argument('--sample',choices=['dovekie','pantheon','des3yr'],default='dovekie')
    parser.add_argument('--evolution',choices=['none','linear','smooth01','smooth03'],default='none')
    parser.add_argument('--calibration',choices=['official_planck','paper_literal'],default='official_planck')
    args = parser.parse_args()
    args.output = args.output.resolve(); args.model_file = args.model_file.resolve()
    if not args.output.is_relative_to(ROOT/'.work') or args.output.exists():
        raise ValueError('--output must be a new directory beneath repository .work.')
    if args.draws<9 or min(args.scales)<=0:
        raise ValueError('Positive scales and at least nine draws are required.')
    args.output.mkdir(parents=True)
    design, acquisition = capture_design(args)
    source_paths = [Path(__file__),Path(__file__).with_name('spectral_surrogate.py'),
        Path(__file__).with_name('modern_fast.py'),Path(__file__).with_name('modern_run.py'),
        Path(__file__).parent.parent/'external_probes/fast_lensing.py',
        Path(__file__).parent.parent/'external_probes/modern_adapter.py',
        Path(native_proposal.__file__)]
    source_hashes = {relative(p):digest(p) for p in source_paths}
    write_json(args.output/'design.json',design)  # Freeze before evaluating outcomes.
    start = time.perf_counter()
    with prohibit_native_spectra(), get_model(configuration(
            surrogate=args.model_file,**design['target'])) as model:
        construction = time.perf_counter()-start
        records = classify(model,design)
        measurements = benchmark(model,records,design['scales'][0])
        theory = model.theory['spectral_surrogate']
        assert theory.exact_calls == 0
        summary, median = summarize(records,design,measurements,theory.envelope)
        envelope = float(theory.envelope)
    full = {'records':records,'benchmark':measurements}
    write_json(args.output/'full-records.json',full)
    history = None
    if args.historical_result:
        previous = json.loads(args.historical_result.read_text())
        for scale, record in summary.items():
            assert record['counts'] == previous['summary'][scale]['counts']
            np.testing.assert_allclose(record['max_whitened_quantiles'],
                previous['summary'][scale]['max_whitened_quantiles'],rtol=0,atol=1e-10)
        history = {'path':relative(args.historical_result),'sha256':digest(args.historical_result),
            'original_source_sha256':previous['script_sha256'],
            'original_design_sha256':previous['design_sha256'],
            'original_warm_median_seconds':previous['warm_median_seconds'],
            'original_summary':previous['summary'],
            'classification_reproduced':True,'whitened_quantiles_agree_atol':1e-10}
    for path in source_paths:
        assert digest(path) == source_hashes[relative(path)], 'Source changed during diagnostic.'
    result = {'status':'completed','diagnostic_design':design,'acquisition':acquisition,
        'source_sha256':source_hashes,'design_sha256':digest(args.output/'design.json'),
        'full_records_path':relative(args.output/'full-records.json'),
        'full_records_sha256':digest(args.output/'full-records.json'),
        'versions':{key:importlib.metadata.version(key) for key in ['camb','cobaya','numpy','scipy']},
        'construction_seconds':construction,'summary':summary,'warm_median_seconds':median,
        'warm_component_median_seconds':{key:float(np.median(
            [r['component_seconds'].get(key,0) for r in measurements[1:]]))
            for key in measurements[-1]['component_seconds']},
        'historical_outcome':history,'exact_spectrum_calls':0,'adopted_envelope':envelope,
        'limitations':LIMITS}
    args.summary.parent.mkdir(parents=True,exist_ok=True)
    write_json(args.summary,result)
    print(json.dumps({'status':'completed','summary':relative(args.summary),
                      'warm_median_seconds':median,'counts':{k:v['counts'] for k,v in summary.items()}}))


if __name__ == '__main__':
    main()
