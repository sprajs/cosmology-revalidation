"""Replay saved mean predictions, compare nonlinear correction shapes."""
from pathlib import Path
import importlib.util,json,hashlib
import numpy as np
from scipy.linalg import solve_triangular

ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('gate',ROOT/'scripts/research_2026_09_26/sed_nonlinear_refit.py')
gate=importlib.util.module_from_spec(spec);spec.loader.exec_module(gate)
r=json.loads((HERE/'result.json').read_text())
pred=np.load(HERE/'fitted-model-means.npz')
with np.load(ROOT/r['sed_coefficient_source']) as data: coeff=data['broad_drift_discovery_rms0.05'][:,3]
observer=np.array(r['frozen_observer_griz_mag'])
model,bands,paths,zp=gate.sed.fr.build_model();offsets={str(x['Filter Name'])[-1]:float(x['Primary Mag']) for x in zp}
groups={};replay_error=0.;objective_error=0.;record=[]
for row in r['cohort']:
    cid=row['CID'];path=gate.BASE/'objectives'/f'objective_{cid}.npz'
    with np.load(path) as data:
        nom={k:data[k] for k in ['MJD','band','model_flux','data_flux','data_fluxerr','zHEL','MWEBV','parameters_x0_x1_c_t0','frozen_flux_covariance']}
    engine=gate.Engine(nom,model,bands,offsets,coeff,observer)
    chol=np.linalg.cholesky(nom['frozen_flux_covariance'])
    target_values={'native_noiseless':engine.flux(np.zeros(4)),'official_mean_sensitivity':nom['model_flux'],'observed_flux':nom['data_flux']}
    for response in [x for x in r['responses'] if x['CID']==cid and x['mask']=='accepted_full']:
        target=response['target'];mode=response['mode']
        baseline=engine.flux(np.array(response['baseline_theta']))
        alternative=engine.flux(np.array(response['alternative_theta']),mode)
        saved=pred[f'{cid}__accepted_full__{target}__{mode}']
        replay_error=max(replay_error,float(np.max(abs(alternative-saved))))
        resid=solve_triangular(chol,alternative-target_values[target],lower=True)
        objective_error=max(objective_error,abs(float(resid@resid)-response['alternative_chi2']))
        change=solve_triangular(chol,alternative-baseline,lower=True)
        groups.setdefault(target,{}).setdefault(mode,[]).append(change)
        record.append(dict(CID=cid,target=target,mode=mode,correction_norm=float(np.linalg.norm(change))))
comparisons={}
for target,modes in groups.items():
    a=np.concatenate(modes['observer']);b=np.concatenate(modes['sed'])
    comparisons[target]={'observer_correction_norm':float(np.linalg.norm(a)),'sed_correction_norm':float(np.linalg.norm(b)),
                         'difference_norm':float(np.linalg.norm(a-b)),
                         'relative_squared_mismatch_to_observer':float(np.sum((a-b)**2)/np.sum(a*a)),
                         'cosine':float(a@b/(np.linalg.norm(a)*np.linalg.norm(b)))}
assert replay_error<1e-10 and objective_error<1e-8
out={'replayed_mean_max_abs_error':replay_error,'recomputed_objective_max_abs_error':objective_error,'nonlinear_correction_shape_comparison':comparisons,'object_norms':record,'inputs_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),HERE/'result.json',HERE/'fitted-model-means.npz']}}
(HERE/'saved-fit-verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ['object_norms','inputs_sha256']},indent=2))
