"""Frozen-vector marginal wavelength/phase diagnostic; no model selection."""
from pathlib import Path
import argparse
import hashlib
import json
import sys

import numpy as np
import pandas as pd
from scipy.linalg import solve_triangular

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.salt_dust_audit.flux_response import build_model, K
from sncosmo.utils import integration_grid

BASE = ROOT / 'runs/research_2026_09_26/astra_design/validation1020'
OUT = ROOT / 'runs/research_2026_09_26/residual_localization'
PROTOCOL = ROOT / 'docs/research-2026-09-26/residual-localization-protocol.md'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--stage', choices=['design', 'score'], required=True)
    stage = parser.parse_args().stage
    dest = OUT / stage
    dest.mkdir(parents=True, exist_ok=False)
    inputs = {}

    def record(path):
        path = Path(path)
        key = str(path.relative_to(ROOT))
        if key not in inputs:
            inputs[key] = hashlib.sha256(path.read_bytes()).hexdigest()
        return path

    record(PROTOCOL); record(__file__)
    (dest / 'executed_source.py').write_bytes(Path(__file__).read_bytes())
    coeff = np.load(record(BASE / 'frozen-discovery-coefficients.npz'), allow_pickle=False)
    vector = coeff['gauge_griz'] @ coeff['basis_mean']
    cohort = pd.read_csv(record(BASE / 'cohort.csv'), dtype={'CID': str})
    _, bands, paths, _ = build_model()
    for p in paths:
        record(p)
    grids = {}
    for b in 'griz':
        bp = bands[b]
        wave, dw = integration_grid(bp.minwave(), bp.maxwave(), 5.)
        weight = wave * bp(wave) * dw
        assert np.all(weight >= 0) and weight.sum() > 0
        grids[b] = (wave, weight / weight.sum())
    rows, masks = [], {}
    if stage == 'score':
        frozen = pd.read_csv(record(OUT / 'design/object-information.csv'), dtype={'CID': str})
        frozen = frozen.set_index(['CID', 'partition'])
        original = pd.read_csv(record(BASE / 'analysis/object-scores.csv'), dtype={'CID': str})
        original = original[original.arm == 'published_mask'].set_index('CID')
        frozen_masks = np.load(record(OUT / 'design/epoch-indices.npz'), allow_pickle=False)
    for item in cohort.itertuples():
        cid = str(item.CID)
        with np.load(record(BASE / f'objectives/objective_{cid}.npz'), allow_pickle=False) as nominal:
            z = float(nominal['zHEL'][0])
            t0 = float(nominal['parameters_x0_x1_c_t0'][3])
        with np.load(record(BASE / f'analysis/objects/{cid}.npz'), allow_pickle=False) as cache:
            time, band = cache['MJD'], cache['band']
            flux, cov = cache['official_flux_model'], cache['exact_covariance']
            jac = cache['jacobian_flux'].copy()
            jac[:, 0] = -K * flux
            observed = cache['observed_flux'] if stage == 'score' else None
        epoch = np.arange(len(time))
        partitions = {'all': np.ones(len(time), dtype=bool),
                      'pre_peak': time <= t0, 'post_peak': time > t0}
        for cutoff in [3000, 3500, 4000]:
            fractions = {b: float(w[lam / (1+z) < cutoff].sum()) for b, (lam, w) in grids.items()}
            keep = np.array([fractions[b] <= .01 for b in band])
            partitions[f'optical_{cutoff}'] = keep
            partitions[f'uv_exposed_{cutoff}'] = ~keep
        for label, keep in partitions.items():
            idx = epoch[keep]
            key = cid + '__' + label
            masks[key] = idx
            if stage == 'score':
                assert np.array_equal(idx, frozen_masks[key])
            row = dict(CID=cid, partition=label, zHEL=z, epochs=len(idx), rank=0,
                       dimension=0, information=0., matched_product=0., gain=0., chi2=0.)
            if len(idx):
                chol = np.linalg.cholesky(cov[np.ix_(idx, idx)])
                jw = solve_triangular(chol, jac[idx], lower=True)
                u, singular, _ = np.linalg.svd(jw, full_matrices=False)
                rank = int(np.sum(singular > singular[0] * 1e-10)) if singular[0] else 0
                q = u[:, :rank]
                pred = solve_triangular(chol, -K * flux[idx] * np.array([vector['griz'.index(b)] for b in band[idx]]), lower=True)
                pred -= q @ (q.T @ pred)
                row.update(rank=rank, dimension=len(idx)-rank, information=float(pred @ pred))
                assert np.linalg.norm(q.T @ pred) < 1e-8
                if stage == 'score':
                    rw = solve_triangular(chol, observed[idx] - flux[idx], lower=True)
                    rw -= q @ (q.T @ rw)
                    row.update(matched_product=float(pred @ rw), gain=float(pred @ rw - pred @ pred / 2), chi2=float(rw @ rw))
                    expected = frozen.loc[(cid, label)]
                    assert abs(row['information'] - expected.information) < 1e-10
                    assert row['rank'] == expected['rank'] and row['epochs'] == expected.epochs
                    if label == 'all':
                        ref = original.loc[cid]
                        assert abs(row['information'] - ref.information) < 1e-10
                        assert abs(row['matched_product'] - ref.matched_filter) < 1e-10
                        assert abs(row['gain'] - ref.fixed_prediction_gain) < 1e-10
            if stage == 'design':
                for k in ['matched_product', 'gain', 'chi2']:
                    row.pop(k)
            rows.append(row)
    table = pd.DataFrame(rows)
    table.to_csv(dest / ('object-information.csv' if stage == 'design' else 'object-scores.csv'), index=False)
    if stage == 'design':
        np.savez_compressed(dest / 'epoch-indices.npz', **masks)
    results = {}
    for label, group in table.groupby('partition'):
        info = float(group.information.sum())
        value = dict(objects=len(group), objects_with_rows=int((group.epochs>0).sum()),
                     rank4_objects=int((group['rank']==4).sum()), residual_dimension=int(group.dimension.sum()),
                     epochs=int(group.epochs.sum()), information=info,
                     rank_deficient_information=float(group.loc[group['rank']<4, 'information'].sum()))
        if stage == 'score':
            total = float(group.matched_product.sum())
            value.update(matched_product=total, gain=float(group.gain.sum()),
                         descriptive_amplitude=total/info if info > 0 else None,
                         chi2=float(group.chi2.sum()),
                         rank_deficient_matched_product=float(group.loc[group['rank']<4,'matched_product'].sum()))
        results[label] = value
    result = dict(stage=stage, scope='Post-discovery mechanism diagnostic; fixed C and local nuisance projection; no p-values or physical correction.',
                  observed_flux_read=stage=='score', partitions=results)
    (dest / 'result.json').write_text(json.dumps(result, indent=2)+'\n')
    outputs = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in dest.iterdir() if p.is_file()}
    (dest / 'manifest.json').write_text(json.dumps(dict(inputs_sha256=inputs, outputs_sha256=outputs), indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
