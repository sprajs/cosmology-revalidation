"""Independent principal-minor and reference-row checks of the grouped SN audit."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / 'sources/repos/PantheonPlusSH0ES__DataRelease/Pantheon+_Data/4_DISTANCES_AND_COVAR'


def load(path):
    a = np.loadtxt(path); n = int(a[0]); c = a[1:].reshape(n, n)
    return (c + c.T) / 2


def main():
    paths = [BASE/'Pantheon+SH0ES_STATONLY.cov', BASE/'Pantheon+SH0ES.dat']
    stat = load(paths[0]); data = pd.read_csv(paths[1], sep=r'\s+')
    valid = (data.zHD.to_numpy() > .1) & ~data.CID.duplicated(keep=False).to_numpy()
    ans = {}
    for label in ['MWEBV', 'MWCOLORLAW']:
        path = BASE/f'sytematic_groupings/Pantheon+SH0ES_122221_{label}.cov'; paths.append(path)
        w = load(path); d = w - stat
        ks = np.argsort(np.where(valid, np.diag(d), -np.inf))[-4:]
        bg = []
        for k in ks:
            v = d[k] / np.sqrt(d[k, k]); bg.append(w - np.outer(v, v))
        pair = d[np.ix_([47, 48], [47, 48])]
        ans[label] = dict(alternative_reference_indices=ks.tolist(),
                          max_baseline_difference_across_four_references=max(float(abs(b-bg[0]).max()) for b in bg),
                          two_row_difference_47_48=pair.tolist(),
                          two_row_min_eigenvalue=float(np.linalg.eigvalsh(pair)[0]))
        assert ans[label]['two_row_min_eigenvalue'] < -.05
        assert ans[label]['max_baseline_difference_across_four_references'] < 1e-7
    ans['inputs_sha256'] = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                            for p in [*paths, Path(__file__)]}
    out = ROOT / 'runs/assumption_audit/bao/group-peer-review.json'
    out.write_text(json.dumps(ans, indent=2) + '\n')
    print(json.dumps(ans, indent=2))


if __name__ == '__main__':
    main()
