"""Frozen synthetic CALWF3 3.7.3 linfit variance algebra, no image input."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    protocol = json.loads((HERE / "protocol.json").read_text())
    for record in protocol["source_files"]:
        path = Path(record["path"])
        assert sha(path) == record["sha256"] and path.stat().st_size == record["bytes"]
    t = np.asarray(protocol["times_seconds"], dtype=np.float64)
    assert len(t) == 7 and np.allclose(np.diff(t), 50)
    sigma = float(protocol["sigma_CDS_electrons"])
    T = t[-1] - t[0]
    cov_poisson_per_rate = np.minimum.outer(t, t)
    eye = np.eye(len(t))
    one = np.ones(len(t))
    e0 = eye[0]
    # Cov(N_i - N_0, N_j - N_0) for independent single-read N_i.
    cov_read = sigma**2 / 2 * (eye + np.outer(one, one)
                               - np.outer(e0, one) - np.outer(one, e0))
    result = []
    for power in protocol["power_branches"]:
        # C pow(0,0) and NumPy both yield one in the equal-weight branch.
        u = np.abs((np.arange(len(t)) - 3) / 3) ** power
        weighted_mean_t = np.dot(u, t) / u.sum()
        Q = np.dot(u, (t - weighted_mean_t)**2)
        h = u * (t - weighted_mean_t) / Q
        assert abs(h.sum()) < 1e-14 and abs(np.dot(h, t)-1) < 1e-14
        # Exact source normalizer with W_i = u_i / sigma_CDS^2, and no
        # denominator clamp for any of the six seven-read cases.
        W = u / sigma**2
        S, Sx, Sxx = W.sum(), np.dot(W,t), np.dot(W,t*t)
        denom = S*Sxx-Sx*Sx
        assert denom > 1e-6
        source_fit_uncert_sq = S/denom
        reported_read = sigma**2 / Q
        assert np.isclose(source_fit_uncert_sq,reported_read,rtol=1e-13)
        actual_read = float(h @ cov_read @ h)
        independent_single_read = sigma**2 / 2 * float(h @ h)
        assert np.isclose(actual_read,independent_single_read,rtol=1e-13,atol=1e-15)
        actual_poisson_per_rate = float(h @ cov_poisson_per_rate @ h)
        increment_route = float(sum((t[k]-t[k-1]) * np.sum(h[k:])**2
                                    for k in range(1,len(t))))
        assert np.isclose(actual_poisson_per_rate,increment_route,rtol=1e-13,atol=1e-15)
        reported_poisson_per_rate = 1/T
        result.append({"power":power,"weighted_time_mean_seconds":weighted_mean_t,
                       "Q_seconds_squared":Q,"source_denom":denom,
                       "h_per_second":h.tolist(),
                       "reported_read_variance_e_per_s_squared":reported_read,
                       "actual_read_variance_e_per_s_squared":actual_read,
                       "reported_to_actual_read_variance":reported_read/actual_read,
                       "reported_poisson_coefficient_per_second":reported_poisson_per_rate,
                       "actual_poisson_coefficient_per_second":actual_poisson_per_rate,
                       "reported_to_actual_poisson_variance":reported_poisson_per_rate/actual_poisson_per_rate,
                       "poisson_increment_route_gap":actual_poisson_per_rate-increment_route,
                       "read_common_zero_cancellation_gap":actual_read-independent_single_read})
    out={"status":"synthetic fixed-weight algebra complete; no real-pixel outcomes read",
         "protocol_sha256":sha(HERE/"protocol.json"),
         "source_ref":protocol["source_ref"],
         "formula_reported":"Var_code = reported_read_variance + total_electron_rate * reported_poisson_coefficient",
         "formula_actual":"Var_fixed_weight = actual_read_variance + total_electron_rate * actual_poisson_coefficient",
         "rate_domain":"nonnegative total source+dark electrons/s; branch power held fixed",
         "results":result}
    (HERE/"result.json").write_text(json.dumps(out,indent=2)+"\n")
    with (HERE/"coefficients.csv").open("w",newline="") as f:
        keys=[k for k in result[0] if k != "h_per_second"]
        w=csv.DictWriter(f,fieldnames=keys)
        w.writeheader()
        w.writerows({k:r[k] for k in keys} for r in result)
    print(json.dumps({"result_sha256":sha(HERE/"result.json"),
                      "coefficients_sha256":sha(HERE/"coefficients.csv"),
                      "branches":len(result),
                      "max_independent_poisson_gap":max(abs(x["poisson_increment_route_gap"]) for x in result)},
                     indent=2))


if __name__ == "__main__":
    main()
