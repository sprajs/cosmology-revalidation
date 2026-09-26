"""Independent least-squares projection and nested-mask covariance checks."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.linalg import solve_triangular

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / 'runs/research_2026_09_26/astra_design/validation1020'
OUT = ROOT / 'runs/research_2026_09_26/residual_localization'


def main():
    scores = pd.read_csv(OUT / 'score/object-scores.csv', dtype={'CID': str})
    indices = np.load(OUT / 'design/epoch-indices.npz', allow_pickle=False)
    coefficients = np.load(BASE / 'frozen-discovery-coefficients.npz', allow_pickle=False)
    vector = coefficients['gauge_griz'] @ coefficients['basis_mean']
    hashes = 0
    for stage in ['design', 'score']:
        manifest = json.loads((OUT / stage / 'manifest.json').read_text())
        for group in ['inputs_sha256', 'outputs_sha256']:
            for name, expected in manifest[group].items():
                # Design source/protocol were explicitly amended after geometry,
                # before scores. Preserved executed source retains its exact hash.
                target = ROOT / name
                if stage == 'design' and name == 'scripts/research_2026_09_26/residual_localization.py':
                    target = OUT / 'design/executed_source.py'
                if stage == 'design' and name == 'docs/research-2026-09-26/residual-localization-protocol.md':
                    target = OUT / 'design/protocol-before-amendment.md'
                assert hashlib.sha256(target.read_bytes()).hexdigest() == expected, name
                hashes += 1
    errors = dict(score=0., cross_covariance=0., orthogonality=0., complement_variance=0.)
    complements = []
    for cid, group in scores.groupby('CID'):
        data = np.load(BASE / f'analysis/objects/{cid}.npz', allow_pickle=False)
        cov = data['exact_covariance']; mean = data['official_flux_model']; y = data['observed_flux']
        jac = data['jacobian_flux'].copy(); jac[:,0] = -.4*np.log(10)*mean
        delta = -.4*np.log(10)*mean*np.array([vector['griz'.index(b)] for b in data['band']])
        kernels = {}
        for row in group.itertuples():
            idx = indices[cid+'__'+row.partition]
            kernel = np.zeros(len(mean))
            if len(idx):
                chol = np.linalg.cholesky(cov[np.ix_(idx,idx)])
                jw = solve_triangular(chol, jac[idx], lower=True)
                pw = solve_triangular(chol, delta[idx], lower=True)
                nuisance = np.linalg.lstsq(jw, pw, rcond=1e-10)[0]
                projected = pw-jw@nuisance
                kernel[idx] = solve_triangular(chol.T, projected, lower=False)
            m = float(kernel@(y-mean)); info = float(kernel@cov@kernel)
            errors['score'] = max(errors['score'], abs(m-row.matched_product), abs(info-row.information))
            kernels[row.partition] = kernel
        for cutoff in [3000,3500,4000]:
            label = f'optical_{cutoff}'; full = kernels['all']; retained = kernels[label]; rest = full-retained
            i_full = float(full@cov@full); i_retained = float(retained@cov@retained)
            cross = float(full@cov@retained); i_rest = float(rest@cov@rest)
            errors['cross_covariance'] = max(errors['cross_covariance'],abs(cross-i_retained))
            errors['orthogonality'] = max(errors['orthogonality'],abs(rest@cov@retained))
            errors['complement_variance'] = max(errors['complement_variance'],abs(i_rest-(i_full-i_retained)))
            complements.append(dict(CID=cid,cutoff=cutoff,information=i_rest,matched_product=float(rest@(y-mean))))
    assert max(errors.values()) < 1e-9, errors
    frame = pd.DataFrame(complements)
    frame.to_csv(OUT / 'nested-complement-object-scores.csv',index=False)
    aggregate = {}
    for cutoff, group in frame.groupby('cutoff'):
        m = float(group.matched_product.sum()); info = float(group.information.sum())
        aggregate[str(cutoff)] = dict(information=info,matched_product=m,
                                    descriptive_amplitude=m/info, gain=m-info/2,
                                    fixed_nominal_noise_sd_amplitude=1/np.sqrt(info))
    result = dict(status='PASS',checked_hashes=hashes,score_rows=len(scores),errors=errors,
                  interpretation='Nested removed information includes cross-band constraints; not a UV-only score. Conditional Gaussian covariance algebra is not physical-null validation.',
                  nested_complements=aggregate,
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (OUT / 'independent-verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
