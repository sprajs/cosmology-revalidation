"""Read existing response matrices only; quantify conditional design information.

This does not fit observed residuals, select a correction, or regenerate audits.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
MATRIX = ROOT / 'runs/salt_dust_audit/flux_response/matrices.npz'
META = ROOT / 'runs/salt_dust_audit/flux_response/selected_objects.csv'
SUMMARY = ROOT / 'runs/salt_dust_audit/flux_response/summary.json'

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    arrays = np.load(MATRIX)
    meta = pd.read_csv(META, dtype={'CID': str}).set_index('CID')
    rows = []
    statistics = {}
    for weighting in ('measurement_only', 'measurement_plus_model'):
        blocks, targets, nuisances = [], [], []
        for cid, row in meta.iterrows():
            k = arrays[f'{cid}__{weighting}_residual_gram']
            blocks.append(k)
            cov = arrays[f'{cid}__{weighting}_covariance']
            chol = np.linalg.cholesky(cov)
            jw = np.linalg.solve(chol, arrays[f'{cid}__jacobian_flux'])
            gw = np.linalg.solve(chol, arrays[f'{cid}__nuisance_flux'])
            r = gw - jw @ np.linalg.lstsq(jw, gw, rcond=None)[0]
            assert np.allclose(k, r.T @ r, atol=1e-7, rtol=1e-7)
            # One common added host E; calibration columns are shared across SNe.
            targets.append(r[:, 2])
            nuisances.append(r[:, 4:])
            if weighting == 'measurement_plus_model':
                bands = arrays[f'{cid}__band']
                phase = (arrays[f'{cid}__mjd'] - row['PKMJD']) / (1 + row['zHEL'])
                rows.append(dict(CID=cid, zHEL=float(row['zHEL']),
                                 nepoch=int(len(bands)),
                                 bands=''.join(sorted(set(bands))),
                                 phase_min=float(min(phase)), phase_max=float(max(phase)),
                                 host_E_001_residual_snr=float(.01*np.sqrt(k[2,2]))))
        gram = sum(blocks)
        target = np.concatenate(targets)
        nuisance = np.vstack(nuisances)
        after = target - nuisance @ np.linalg.lstsq(nuisance, target, rcond=None)[0]
        singular = np.linalg.svd(nuisance, compute_uv=False)
        info = float(target @ target)
        residual_info = float(after @ after)
        # A calibration-like alternative with the same perturbation unit is shown
        # only to compare design sensitivity, never as an empirical prior.
        statistics[weighting] = dict(
            common_host_E_001_residual_snr=float(.01*np.sqrt(info)),
            common_host_E_001_residual_snr_after_shared_MW_and_griz=float(.01*np.sqrt(residual_info)),
            information_retained_after_shared_MW_and_griz=residual_info/info,
            common_host_E_conditional_sigma_mag=float(1/np.sqrt(info)),
            common_host_E_sigma_after_shared_MW_and_griz_mag=float(1/np.sqrt(residual_info)),
            shared_nuisance_singular_values=singular.tolist(),
            common_zp_i_001_residual_snr=float(.01*np.sqrt(gram[7,7])),
            common_zp_z_001_residual_snr=float(.01*np.sqrt(gram[8,8])),
            residual_gram=gram.tolist())
    pd.DataFrame(rows).to_csv(OUT/'object-design.csv', index=False)
    result = dict(
        purpose='Design-only information from local derivatives; no observed residual fit.',
        conditioning='64 deterministic redshift-ranked published accepted light curves; frozen per-object covariance; no selection, clipping, training, or cross-object model covariance.',
        statistic='Known same-signed host E(B-V) addition to all 64 objects after independent SALT coordinates. The sigma is inverse square root of local information, not a real reddening uncertainty.',
        limitation='Shared training covariance is absent; unknown E_i and population dust parameters are harder. Pooling common E is an optimistic design calculation under these fixed covariance assumptions.',
        count=len(meta), epochs=sum(r['nepoch'] for r in rows),
        results=statistics,
        bands_counts=pd.Series([r['bands'] for r in rows]).value_counts().to_dict(),
        inputs_sha256={str(p.relative_to(ROOT)):digest(p) for p in (MATRIX,META,SUMMARY)},
        source_sha256=digest(Path(__file__)),
        numpy_version=np.__version__, pandas_version=pd.__version__)
    (OUT/'design-geometry.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:{kk:vv for kk,vv in v.items() if kk not in ('residual_gram','shared_nuisance_singular_values')} for k,v in statistics.items()}, indent=2))
    print('Bands',result['bands_counts'])

if __name__ == '__main__':
    main()
