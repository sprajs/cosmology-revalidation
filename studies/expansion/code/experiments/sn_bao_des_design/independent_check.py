"""Independent outcome-free DES forecast check; never loads MU."""
import os
for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[key] = "1"

import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import solve_triangular
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
DES = ROOT / "sources/repos/des-science__DES-SN5YR@1.3/4_DISTANCES_COVMAT"
BAO = ROOT / "sources/repos/CobayaSampler__bao_data/desi_bao_dr2"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def flattened_gzip(path):
    with gzip.open(path, "rt") as f:
        values = np.fromstring(f.read(), sep=" ")
    n = int(values[0])
    assert len(values) == 1 + n * n
    return values[1:].reshape((n, n))


def main():
    frozen = json.loads((OUT / "protocol.json").read_text())
    result = json.loads((OUT / "result.json").read_text())
    assert result["protocol_sha256"] == sha(OUT / "protocol.json")
    assert sha(OUT / "frozen_source.py") == sha(OUT / "executed_source.py") == sha(ROOT / "scripts/research_2026_09_26/sn_bao_des_design.py")
    for rel, digest in frozen["input_sha256"].items():
        assert sha(ROOT / rel) == digest
    for filename, digest in result["artifacts_sha256"].items():
        assert sha(OUT / filename) == digest
    source = (OUT / "executed_source.py").read_text()
    design_body = source[source.index("def design():"):source.index("if __name__")]
    assert "usecols=['CID','IDSURVEY','zHD','zHEL','MUERR_FINAL']" in design_body
    assert "['MU']" not in design_body and "usecols=None" not in design_body

    # Independent release-row and covariance construction, without the MU column.
    hd = pd.read_csv(DES / "DES-SN5YR_HD.csv", usecols=["CID", "IDSURVEY", "zHD", "zHEL", "MUERR_FINAL"], dtype={"CID": str})
    meta = pd.read_csv(DES / "DES-SN5YR_HD+MetaData.csv", usecols=["CID", "IDSURVEY", "zHD", "zHEL", "MUERR_FINAL"], dtype={"CID": str})
    assert list(hd.columns) != list(meta.columns)  # files order their headers differently
    assert hd.equals(meta[list(hd.columns)])
    sys = flattened_gzip(DES / "STAT+SYS.txt.gz")
    statonly = flattened_gzip(DES / "STATONLY.txt.gz")
    assert len(hd) == len(sys) == 1829 and np.count_nonzero(statonly) == 0
    assert np.max(abs(sys - sys.T)) == 0
    des_rows = np.flatnonzero(hd.IDSURVEY.to_numpy() == 10)
    assert len(des_rows) == 1635 and np.array_equal(des_rows, result["des_file_row_indices"])
    selected = hd.iloc[des_rows].reset_index(drop=True)
    C = sys[np.ix_(des_rows, des_rows)].copy()
    C[np.diag_indices_from(C)] += selected.MUERR_FINAL.to_numpy() ** 2
    np.linalg.cholesky(C)

    bao = pd.read_csv(BAO / "desi_gaussian_bao_ALL_GCcomb_mean.txt", sep=r"\s+", comment="#", names=["z", "value", "kind"])
    bao_cov = np.loadtxt(BAO / "desi_gaussian_bao_ALL_GCcomb_cov.txt")
    bz = selected.zHD.to_numpy()
    bi = np.flatnonzero((bao.kind == "DM_over_rs") & (bao.z >= bz.min()) & (bao.z <= bz.max()))
    nodes = bao.z.to_numpy()[bi]
    values = bao.value.to_numpy()[bi]
    J = np.diag(5 / np.log(10) / values)
    bao_mag_cov = J @ bao_cov[np.ix_(bi, bi)] @ J
    table = pd.read_csv(OUT / "design-table.csv")
    with np.load(OUT / "weights.npz", allow_pickle=False) as data:
        weights = {k: data[k] for k in data.files}
    assert np.array_equal(nodes, weights["bao_z"])
    assert np.array_equal(bz, weights["sn_z"])
    assert np.array_equal(selected.zHEL.to_numpy(), weights["sn_zHEL"])
    assert np.max(abs(bao_mag_cov - weights["bao_magnitude_cov"])) < 1e-15
    assert len(table) == 12
    all_admitted = np.ones(len(nodes), dtype=bool)
    for node_index, node in enumerate(nodes):
        subset = table[abs(table.z_node - node) < 1e-12]
        assert len(subset) == 4
        for _, row in subset.iterrows():
            prescribed = row.distinct_CID >= 20 and row.distinct_below >= 5 and row.distinct_above >= 5 and row.condition < 1e8 and row.max_finite_family_bias_mag <= .005
            assert bool(row.admitted) == bool(prescribed)
            all_admitted[node_index] &= prescribed
    assert np.array_equal(all_admitted, weights["admitted"])
    assert np.array_equal(nodes[all_admitted], result["admitted_nodes"])

    # Primary local GLS recomputed from a whitened QR, independent of the
    # primary code's normal-equation solve.
    W = np.zeros((len(nodes), len(bz)))
    local_details = []
    for j, node in enumerate(nodes):
        halfwidth = .10
        local = np.flatnonzero(abs(np.log1p(bz) - np.log1p(node)) <= halfwidth)
        x = (np.log1p(bz[local]) - np.log1p(node)) / halfwidth
        X = np.column_stack([np.ones(len(x)), x, x*x])
        subC = C[np.ix_(local, local)]
        L = np.linalg.cholesky(subC)
        Q, R = np.linalg.qr(solve_triangular(L, X, lower=True), mode="reduced")
        direction = solve_triangular(R.T, [1., 0., 0.], lower=True)
        W[j, local] = solve_triangular(L.T, Q @ direction, lower=False)
        row = table[(table.design == "primary") & (abs(table.z_node - node) < 1e-12)].iloc[0]
        assert len(local) == int(row.observation_rows)
        assert selected.CID.iloc[local].nunique() == int(row.distinct_CID)
        assert selected.CID.iloc[local[x < 0]].nunique() == int(row.distinct_below)
        assert selected.CID.iloc[local[x > 0]].nunique() == int(row.distinct_above)
        assert np.max(abs(W[j, local] @ X - [1, 0, 0])) < 1e-12
        local_details.append({"z": float(node), "rows": int(len(local)), "cid": int(selected.CID.iloc[local].nunique())})
    max_weight_difference = float(np.max(abs(W - weights["primary"])))
    sn_cov = W @ C @ W.T
    V = sn_cov + bao_mag_cov
    assert max_weight_difference < 1e-10
    assert np.max(abs(sn_cov - weights["primary_sn_cov"])) < 1e-12
    assert np.max(abs(V - weights["primary_total_cov"])) < 1e-12

    # Helmert contrasts explicitly remove the unknown common intercept in
    # the three-node forecast; no standardised magnitude is required.
    H = np.array([[1., -1., 0.], [1., 1., -2.]]) / np.array([[np.sqrt(2.)], [np.sqrt(6.)]])
    assert np.max(abs(H @ np.ones(3))) < 1e-15
    g = (np.log1p(nodes) - np.log1p(nodes[0])) / (np.log1p(nodes[-1]) - np.log1p(nodes[0]))
    hg = H @ g
    hv = H @ V @ H.T
    se = float(1 / np.sqrt(hg @ np.linalg.solve(hv, hg)))
    assert abs(se - result["forecast"]["primary"]["endpoint_drift_se_mag"]) < 1e-12
    power = {}
    for a in (.02, .05, .10):
        critical = norm.ppf(.975)
        power[str(a)] = float(norm.cdf(-critical-a/se) + norm.sf(critical-a/se))
        assert abs(power[str(a)] - result["forecast"]["primary"]["two_sided_5pct_power"][str(a)]) < 1e-12
    check = {
        "gate": "pass",
        "scope": "Outcome-free DES release covariance/GLS/Helmert forecast verification; no MU loaded.",
        "protocol_sha256": sha(OUT / "protocol.json"),
        "source_sha256": sha(Path(__file__)),
        "source_and_input_hashes_pass": True,
        "release_rows": len(hd),
        "des_rows": len(selected),
        "des_rows_contiguous_first": bool(np.array_equal(des_rows, np.arange(1635))),
        "statonly_nonzero_entries": int(np.count_nonzero(statonly)),
        "systematic_covariance_asymmetry": float(np.max(abs(sys - sys.T))),
        "covariance_rule": "STAT+SYS systematic submatrix plus MUERR_FINAL squared diagonal once",
        "bao_marginal_dm_nodes": nodes.tolist(),
        "admitted_all_settings": nodes[all_admitted].tolist(),
        "primary_local_support": local_details,
        "primary_qr_weight_max_abs_difference": max_weight_difference,
        "primary_sn_cov_max_abs_difference": float(np.max(abs(sn_cov - weights["primary_sn_cov"]))),
        "primary_total_cov_max_abs_difference": float(np.max(abs(V - weights["primary_total_cov"]))),
        "primary_helmert_free_intercept_se_mag": se,
        "primary_two_sided_power": power,
    }
    (OUT / "independent-check.json").write_text(json.dumps(check, indent=2) + "\n")
    print(json.dumps(check, indent=2))


if __name__ == "__main__":
    main()
