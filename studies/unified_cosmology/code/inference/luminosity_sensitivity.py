"""Integrate declared SN luminosity modes and audit importance-weight overlap.

The CMB/BAO factors are unchanged. Only exact background distances and the
released SN likelihood are recomputed. Provisional chains never qualify for
posterior claims, even when their conditional weight diagnostics look benign.
"""
import argparse
from collections import Counter
from contextlib import contextmanager
import json
from pathlib import Path
import time

import numpy as np
from scipy.linalg import cholesky, solve_triangular, cho_factor, cho_solve
from scipy.interpolate import CubicSpline
from scipy.special import log_ndtr, ndtr, logsumexp

from late_geometry import sample_path
from likelihood import expansion_diagnostics
from target_identity import ROOT, digest

HERE = Path(__file__).resolve().parent
DESIGN_PATH = HERE/'luminosity-sensitivity-design.json'
RESULTS = HERE.parents[1]/'results/inference'


def log_normal_interval(lower, upper):
    """Stable log[Phi(upper)-Phi(lower)], including both remote tails."""
    if not lower < upper:
        raise ValueError('Normal interval must have positive width.')
    if upper <= 0:
        large, small = log_ndtr(upper), log_ndtr(lower)
        return float(large+np.log(-np.expm1(small-large)))
    if lower >= 0:
        large, small = log_ndtr(-lower), log_ndtr(-upper)
        return float(large+np.log(-np.expm1(small-large)))
    return float(np.log1p(-ndtr(lower)-ndtr(-upper)))


class IntegratedLuminosity:
    """Gaussian distance likelihood after the same improper flat-M integral."""
    def __init__(self, z, zhel, observed, covariance):
        self.z = np.asarray(z); self.zhel = np.asarray(zhel)
        self.observed = np.asarray(observed); self.n = len(self.z)
        self.factor = cholesky(covariance, lower=True)
        self.constant = solve_triangular(self.factor,np.ones(self.n),lower=True)
        self.A = float(self.constant@self.constant)
        self.unit_constant = self.constant/np.sqrt(self.A)
        self.lognorm = (2*np.log(np.diag(self.factor)).sum()+np.log(self.A)
                        +(self.n-1)*np.log(2*np.pi))
        self.linear = self.project(np.log1p(self.z)/np.log(2))
        self.linear_information = float(self.linear@self.linear)
        self.knots = np.array([0.,.1,.4,.8,1.3])
        self.basis = CubicSpline(self.knots,np.eye(5)[:,1:],bc_type='natural')(self.z)
        self.spline = self.project(self.basis)
        gram = self.spline.T@self.spline
        self.gaussian = {}
        for name, sigma in [('smooth01',.1),('smooth03',.3)]:
            matrix = np.eye(4)/sigma**2+gram
            factor = cho_factor(matrix,lower=True)
            logdet = np.linalg.slogdet(np.eye(4)+sigma**2*gram)[1]
            self.gaussian[name] = (factor,float(logdet))

    def project(self, vector):
        # Remove an arbitrary constant first to protect large-magnitude offsets.
        vector = np.asarray(vector)
        centered = vector-vector.mean(axis=0)
        white = solve_triangular(self.factor,centered,lower=True)
        if white.ndim == 1:
            return white-self.unit_constant*(self.unit_constant@white)
        return white-np.outer(self.unit_constant,self.unit_constant@white)

    def from_prediction(self, prediction):
        residual = self.project(np.asarray(prediction)-self.observed)
        baseline = -.5*(float(residual@residual)+self.lognorm)
        a = self.linear_information; b = float(self.linear@residual)
        if a <= np.finfo(float).eps:
            if abs(b)>1e-10:
                raise ArithmeticError('Degenerate linear mode has nonzero score.')
            linear = 0.
        else:
            mean = -b/a; sigma = 1/np.sqrt(a)
            probability = log_normal_interval((-.5-mean)/sigma,(.5-mean)/sigma)
            # Prior width is exactly 1 mag; retain its normalization explicitly.
            linear = .5*b*b/a+.5*np.log(2*np.pi/a)-np.log(1.)+probability
        ratios = {'linear':float(linear)}
        v = self.spline.T@residual
        for name, (factor,logdet) in self.gaussian.items():
            ratios[name] = float(.5*v@cho_solve(factor,v)-.5*logdet)
        return {'baseline_SN_loglike':float(baseline),'log_likelihood_ratios':ratios,
                'linear_conditional_untruncated_mean':float(-b/a) if a>0 else None,
                'linear_conditional_sigma':float(1/np.sqrt(a)) if a>0 else None}


@contextmanager
def background_only():
    import camb
    original = camb.get_results
    def forbidden(*args,**kwargs):
        raise RuntimeError('Native CMB spectrum calls are forbidden in this lane.')
    camb.get_results = forbidden
    try:
        yield
    finally:
        camb.get_results = original


def exact_background_record(point, settings, sn):
    import camb
    import sys
    sys.path.insert(0,str(HERE.parent/'external_probes'))
    from modern_adapter import modern_info
    extra = modern_info(settings['model'])['theory']['camb']['extra_args']
    cosmology = {key:point[key] for key in ['H0','ombh2','omch2','ns','tau']}
    cosmology.update(As=1e-10*np.exp(point['logA']),
        w=point.get('w',-1.),wa=point.get('wa',0.))
    with background_only():
        parameters = camb.set_params(**cosmology,**extra)
        background = camb.get_background(parameters)
    unique, inverse = np.unique(sn.z,return_inverse=True)
    da = background.angular_diameter_distance(unique)[inverse]
    prediction = 5*np.log10(da*(1+sn.z)*(1+sn.zhel))+25
    record = sn.from_prediction(prediction)
    record['spectral_coordinates'] = [100*background.cosmomc_theta(),parameters.ombh2,
        parameters.omch2,np.log(1e10*parameters.InitPower.As),parameters.InitPower.ns,
        parameters.Reion.optical_depth,parameters.DarkEnergy.w,parameters.DarkEnergy.wa]
    grid = np.concatenate([np.arange(5)*.001+z for z in [0.,.5,1.]])
    record['derived'] = dict(expansion_diagnostics(grid,background.hubble_parameter(grid)),
        omegam=float(parameters.omegam),rdrag=float(background.get_derived_params()['rdrag']))
    return record


def weighted_summary(values, weights):
    order = np.argsort(values); x = np.asarray(values)[order]; w = weights[order]
    mean = float(w@x); sd = float(np.sqrt(w@((x-mean)**2)))
    return {'mean':mean,'sd':sd,'quantiles_025_16_50_84_975':
            np.interp([.025,.16,.5,.84,.975],np.cumsum(w)-w/2,x).tolist()}


def weight_diagnostics(logweights, values, groups, gates):
    import arviz as az
    lw = np.asarray(logweights); groups = np.asarray(groups)
    if not np.isfinite(lw).all():
        return {'status':'nonfinite_weights','qualified_overlap':False}
    weights = np.exp(lw-logsumexp(lw))
    constant_weights = np.ptp(lw)<1e-14
    _, shape = az.psislw(lw) if not constant_weights else (None,np.nan)
    k = float(shape) if np.isfinite(shape) else None
    ess = float(1/(weights@weights))
    failures = []
    if ess < gates['minimum_raw_weight_ESS']: failures.append('raw_weight_ESS')
    if not constant_weights and (k is None or k >= gates['maximum_Pareto_k_exclusive']): failures.append('Pareto_k')
    chains = np.unique(groups)
    if len(chains)<gates['minimum_independent_chains']: failures.append('independent_chains')
    per_chain = {}; posterior = {}; stability = {}
    for group in chains:
        keep = groups==group; fraction = float(weights[keep].sum())
        cw = weights[keep]/fraction if fraction>0 else np.zeros(sum(keep))
        ce = float(1/(cw@cw)) if fraction>0 else 0.
        per_chain[str(group)] = {'points':int(sum(keep)),'weight_fraction':fraction,'raw_weight_ESS':ce}
        if ce<gates['minimum_per_chain_raw_weight_ESS']: failures.append(f'chain_{group}_weight_ESS')
        if not gates['minimum_chain_weight_fraction']<=fraction<=gates['maximum_chain_weight_fraction']:
            failures.append(f'chain_{group}_weight_fraction')
    for name, x in values.items():
        x = np.asarray(x); summary = weighted_summary(x,weights); posterior[name] = summary
        means = {str(group):(float(weights[groups==group]@x[groups==group]/weights[groups==group].sum())
                 if weights[groups==group].sum()>0 else None) for group in chains}
        batches = {}
        for count in gates['batch_partitions_per_chain']:
            numerators = []; denominators = []; sizes = []
            for group in chains:
                for indices in np.array_split(np.flatnonzero(groups==group),count):
                    sizes.append(len(indices)); denominators.append(float(weights[indices].sum()))
                    numerators.append(float(weights[indices]@(x[indices]-summary['mean'])))
            mcse = float(np.sqrt(np.var(numerators,ddof=1)/len(numerators))/np.mean(denominators))
            batches[str(count)] = {'mean_mcse':mcse,'minimum_points_per_batch':min(sizes),
                                  'minimum_weight_per_batch':min(denominators),'maximum_weight_per_batch':max(denominators)}
        sd = summary['sd']; constant = np.ptp(x)<1e-12
        support_degenerate = not constant and sd<=np.finfo(float).eps
        deviation = (None if any(v is None for v in means.values()) or support_degenerate else
                     0. if constant else max(abs(v-summary['mean']) for v in means.values())/sd)
        max_mcse = (None if support_degenerate else 0. if constant else
                    max(v['mean_mcse'] for v in batches.values())/sd)
        stability[name] = {'independent_chain_weighted_means':means,
            'maximum_chain_mean_deviation_in_pooled_SD':deviation,
            'maximum_batch_MCSE_in_pooled_SD':max_mcse,'batch_sensitivity':batches}
        if deviation is None or deviation>gates['maximum_chain_mean_deviation_in_pooled_SD']: failures.append(name+'_chain_means')
        if max_mcse is None or max_mcse>gates['maximum_weighted_batch_MCSE_in_pooled_SD']: failures.append(name+'_batch_MCSE')
    if stability and min(v['minimum_points_per_batch'] for s in stability.values() for v in s['batch_sensitivity'].values())<gates['minimum_points_per_batch']:
        failures.append('batch_size')
    return {'status':'passed_conditional_overlap_gates' if not failures else 'insufficient_overlap_or_stability',
        'qualified_overlap':not failures,'raw_weight_ESS':ess,'Pareto_k':k,
        'Pareto_k_status':'constant_weights_no_tail_to_fit' if constant_weights else 'finite' if k is not None else 'nonfinite_tail_fit',
        'largest_normalized_weight':float(weights.max()),'log_weight_quantiles':np.quantile(lw,[0,.025,.5,.975,1]).tolist(),
        'per_chain':per_chain,'weighted_summaries_for_diagnostics':posterior,
        'weighted_chain_stability':stability,'failed_gates':failures,
        'weights_used':'raw, untrimmed; Pareto smoothing is diagnostic only'}


def select_provisional(folder, output, count):
    paths = sorted(folder.glob('chain.[0-9]*.txt'))
    if len(paths)!=4: raise ValueError('Provisional comparison requires four chain snapshots.')
    manifest = json.loads((folder/'run-0.json').read_text())
    settings = manifest['arguments']; assert settings['evolution']=='none'
    for rank in range(4):
        assert json.loads((folder/f'run-{rank}.json').read_text())['target_identity']==manifest['target_identity']
    names = [key for key,value in manifest['target_identity']['configuration']['params'].items()
             if isinstance(value,dict) and 'prior' in value]
    snapshots = output/'chain-snapshots'; snapshots.mkdir()
    rows = []; origins = []; columns = None; hashes = {}
    for path in paths:
        raw = path.read_bytes(); lines = raw.decode().splitlines()
        columns_here = lines[0].lstrip('# ').split()
        if columns is None: columns = columns_here
        assert columns==columns_here
        complete = [line for line in lines[1:] if len(line.split())==len(columns)]
        if not complete: raise ValueError('A provisional chain has no complete rows.')
        dest = snapshots/path.name; dest.write_text(lines[0]+'\n'+'\n'.join(complete)+'\n')
        hashes[str(dest.relative_to(ROOT))] = digest(dest)
        data = np.loadtxt(dest,ndmin=2); frequency = data[:,0]
        assert np.allclose(frequency,np.round(frequency)) and np.all(frequency>0)
        origin = np.repeat(np.arange(len(data)),frequency.astype(int))
        origin = origin[int(.3*len(origin)):]
        rows.append(data); origins.append(origin)
    length = min(map(len,origins)); per_chain = min(count//4,length)
    if per_chain<40: raise ValueError('Too few provisional iterations for diagnostic batching.')
    rng = np.random.default_rng(2727869); points = []; groups = []; baseline = []; locations = []
    for group,(data,origin) in enumerate(zip(rows,origins)):
        origin = origin[-length:]
        indices = np.floor((np.arange(per_chain)+rng.random(per_chain))*length/per_chain).astype(int)
        for index in indices:
            row = data[origin[index]]
            points.append({key:float(row[columns.index(key)]) for key in names})
            groups.append(group); baseline.append(float(-.5*row[columns.index('chi2__released_sn')]))
            locations.append({'chain':paths[group].name,'row':int(origin[index]),'common_expanded_index':int(index)})
    return {'mode':'provisional','settings':settings,'points':points,'groups':groups,
        'base_logweights':[0.]*len(points),'stored_baseline_SN_loglikes':baseline,'locations':locations,
        'chain_sha256':hashes,'parent_target_identity':manifest['target_identity'],
        'parent_gates_passed':False,'qualification':'No exact CMB/proposal correction and no parent convergence qualification; overlap diagnosis only.'}


def select_exact(folder, exact_summary, diagnostics):
    selection = json.loads((folder/'selection.json').read_text())
    check = json.loads(diagnostics.read_text()); summary = json.loads(exact_summary.read_text())
    assert selection['settings']['evolution']=='none'
    assert selection['diagnostics_sha256']==digest(diagnostics)
    assert summary['selection_sha256']==digest(folder/'selection.json')
    assert summary['code_sha256']==digest(HERE/'exact_correction.py')
    records = []; identities = {}
    for index, point in enumerate(selection['points']):
        path = folder/f'{index:05d}.json'; record = json.loads(path.read_text())
        assert record['index']==index and record['point']==point
        assert record['target_identity']==selection['correction_identity']
        if record['status']!='finite': raise ValueError(f'Exact correction failed at {index}.')
        assert abs(record['log_weight']-(record['exact_logpost']-record['proposal_logpost']))<1e-10
        records.append(record); identities[str(path.relative_to(ROOT))]=digest(path)
    parent = json.loads((folder.parent/'run-0.json').read_text())['target_identity']
    return {'mode':'exact_corrected','settings':selection['settings'],'points':selection['points'],
        'groups':selection['groups'],'base_logweights':[r['log_weight'] for r in records],
        'stored_baseline_SN_loglikes':[r['exact_loglikes']['released_sn'] for r in records],
        'parent_target_identity':parent,'selection_sha256':digest(folder/'selection.json'),
        'exact_summary_sha256':digest(exact_summary),'parent_diagnostics_sha256':digest(diagnostics),
        'exact_summary_path':str(exact_summary.resolve().relative_to(ROOT)),
        'parent_diagnostics_path':str(diagnostics.resolve().relative_to(ROOT)),
        'exact_record_sha256':identities,
        'parent_gates_passed':check['status']=='passed' and summary['status']=='passed_importance_weight_gates',
        'qualification':'Requires additional alternative-specific overlap and weighted-stability checks.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--provisional-folder',type=Path); mode.add_argument('--exact-folder',type=Path)
    mode.add_argument('--frozen-selection',type=Path,help='Replay a previously recorded point selection')
    parser.add_argument('--exact-summary',type=Path); parser.add_argument('--parent-diagnostics',type=Path)
    parser.add_argument('--points',type=int,default=2000)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--summary',type=Path,default=RESULTS/'luminosity-sensitivity.json')
    args = parser.parse_args(); args.output = args.output.resolve()
    if not args.output.is_relative_to(ROOT/'.work') or args.output.exists():
        raise ValueError('Output must be a new ignored directory under .work.')
    args.output.mkdir(parents=True)
    design = json.loads(DESIGN_PATH.read_text()); gates = design['overlap_gates']
    if args.frozen_selection:
        selection = json.loads(args.frozen_selection.read_text())
        assert selection['settings']['evolution']=='none'
        for field in ['chain_sha256','exact_record_sha256']:
            for path, expected in selection.get(field,{}).items():
                assert digest(ROOT/path)==expected
        if selection['mode']=='exact_corrected':
            for stem in ['exact_summary','parent_diagnostics']:
                assert digest(ROOT/selection[stem+'_path'])==selection[stem+'_sha256']
    elif args.exact_folder:
        if not args.exact_summary or not args.parent_diagnostics: parser.error('Exact mode requires --exact-summary and --parent-diagnostics.')
        selection = select_exact(args.exact_folder.resolve(),args.exact_summary,args.parent_diagnostics)
    else:
        selection = select_provisional(args.provisional_folder.resolve(),args.output,args.points)
    source_paths = [Path(__file__),DESIGN_PATH,HERE/'likelihood.py',HERE/'late_geometry.py',
        HERE.parent/'external_probes/modern_adapter.py',HERE.parent/'external_probes/adapter.py']
    source = {str(p.relative_to(ROOT)):digest(p) for p in source_paths}
    dataset = sample_path(selection['settings']['sample'])
    assert digest(dataset)==selection['parent_target_identity']['sample_sha256']
    for filename, expected in selection['parent_target_identity']['source_sha256'].items():
        assert digest(ROOT/filename)==expected, filename
    (args.output/'selection.json').write_text(json.dumps(selection,indent=2)+'\n')
    with np.load(dataset) as data:
        sn = IntegratedLuminosity(*(data[key] for key in ['zHD','zHEL','MU','covariance']))
    records = []; start = time.monotonic()
    for index, point in enumerate(selection['points']):
        try:
            record = exact_background_record(point,selection['settings'],sn)
            difference = record['baseline_SN_loglike']-selection['stored_baseline_SN_loglikes'][index]
            tolerance = .005 if selection['mode']=='provisional' else 1e-5
            record.update(index=index,status='finite',baseline_SN_difference=float(difference))
            if abs(difference)>tolerance: record['status']='baseline_density_mismatch'
        except Exception as error:
            record = {'index':index,'status':'exception','exception':repr(error)}
        records.append(record)
        if (index+1)%200==0: print(json.dumps({'evaluated':index+1,'total':len(selection['points'])}),flush=True)
    rows_path = args.output/'records.json'; rows_path.write_text(json.dumps(records,indent=2)+'\n')
    failures = [r for r in records if r['status']!='finite']
    out = {'scope':design['scope'],'mode':selection['mode'],'points':len(records),
        'design_sha256':digest(DESIGN_PATH),'source_sha256':source,'SN_input_sha256':digest(dataset),
        'selection_path':str((args.output/'selection.json').relative_to(ROOT)),
        'selection_sha256':digest(args.output/'selection.json'),'records_path':str(rows_path.relative_to(ROOT)),
        'records_sha256':digest(rows_path),'parent_gates_passed':selection['parent_gates_passed'],
        'seconds':time.monotonic()-start,'native_spectrum_calls':0,'failures':failures,
        'status':'failed_evaluation' if failures else 'provisional_overlap_diagnosis'}
    if args.frozen_selection:
        out['replayed_selection_sha256'] = digest(args.frozen_selection)
    if not failures:
        merged = [dict(point,**r['derived']) for point,r in zip(selection['points'],records)]
        names = sorted(set.intersection(*(set(r) for r in merged)))
        values = {name:np.array([r[name] for r in merged]) for name in names}
        base = np.asarray(selection['base_logweights']); groups = np.asarray(selection['groups'])
        model_file = Path(selection['settings']['surrogate'])
        assert digest(model_file)==selection['parent_target_identity']['surrogate_sha256']
        with np.load(model_file) as model:
            spectral_x = np.linalg.solve(model['coordinate_cholesky'],
                (np.array([r['spectral_coordinates'] for r in records])-model['centre']).T).T
        maximum_coordinate = np.max(abs(spectral_x),axis=1)
        out['spectral_model_sha256'] = digest(model_file)
        out['observed_maximum_absolute_whitened_coordinate'] = float(max(maximum_coordinate))
        out['observed_points_outside_envelope4'] = int(sum(maximum_coordinate>4))
        out['maximum_baseline_SN_loglike_difference'] = max(abs(r['baseline_SN_difference']) for r in records)
        out['baseline_weight_diagnostics'] = weight_diagnostics(base,values,groups,gates)
        out['alternatives'] = {}
        for alternative in ['linear','smooth01','smooth03']:
            ratios = np.array([r['log_likelihood_ratios'][alternative] for r in records])
            result = weight_diagnostics(base+ratios,values,groups,gates)
            result['SN_log_likelihood_ratio_quantiles'] = np.quantile(ratios,[0,.025,.5,.975,1]).tolist()
            normalized = np.exp(base+ratios-logsumexp(base+ratios))
            result['proposal_cost_diagnostic'] = {
                'weighted_fraction_of_observed_points_outside_envelope4':float(normalized@(maximum_coordinate>4)),
                'weighted_fraction_of_observed_points_above3p5':float(normalized@(maximum_coordinate>3.5)),
                'maximum_coordinate_weighted_summary':weighted_summary(maximum_coordinate,normalized),
                'scope':'Observed provisional/exact-weighted point cloud only; no coverage claim about unsampled alternative tails or interpolation accuracy.'}
            result['posterior_qualified'] = bool(selection['mode']=='exact_corrected'
                and selection['parent_gates_passed'] and len(records)>=gates['minimum_exact_points']
                and result['qualified_overlap'] and out['baseline_weight_diagnostics']['qualified_overlap'])
            out['alternatives'][alternative] = result
        if selection['mode']=='exact_corrected':
            out['status'] = ('qualified_sensitivity' if all(r['posterior_qualified'] for r in out['alternatives'].values())
                             else 'insufficient_overlap_or_parent_qualification')
        out['interpretation'] = ('Diagnostic-only weighted summaries; do not report alternative posterior intervals unless that alternative posterior_qualified is true. '
            'These luminosity priors are declared sensitivities, never measured age corrections. Missing remote support cannot be repaired by smoothing weights.')
    for path in source_paths: assert digest(path)==source[str(path.relative_to(ROOT))]
    args.summary.parent.mkdir(parents=True,exist_ok=True)
    args.summary.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':out['status'],'points':len(records),'failures':len(failures),
        'alternatives':{k:{q:v[q] for q in ['raw_weight_ESS','Pareto_k','posterior_qualified']} for k,v in out.get('alternatives',{}).items()}}))


if __name__=='__main__':
    main()
