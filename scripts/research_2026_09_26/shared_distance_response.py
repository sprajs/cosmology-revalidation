"""Local pre-BBC distance response; no physical correction is inferred."""
from pathlib import Path
import argparse
import hashlib
import json

import numpy as np
import pandas as pd
from scipy.linalg import solve_triangular

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "runs/research_2026_09_26/astra_design"
PROTOCOL = ROOT / "docs/research-2026-09-26/distance-response-protocol.md"
K = .4 * np.log(10)
STANDARDIZE = np.array([1., .16087, -3.1178, 0.])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--modes", type=int, choices=[9, 12], required=True)
    args = parser.parse_args()
    out = ROOT / f"runs/research_2026_09_26/shared_distance_response_{args.modes}"
    out.mkdir(exist_ok=False)
    inputs = {}

    def record(path):
        path = Path(path)
        key = str(path.relative_to(ROOT))
        if key not in inputs:
            inputs[key] = hashlib.sha256(path.read_bytes()).hexdigest()
        return path

    def read(path):
        return np.load(record(path), allow_pickle=False)

    record(__file__)
    record(PROTOCOL)
    discovery = read(BASE / "shared43/projected-modes.npz")
    validation = read(BASE / "shared1020/projected-modes.npz")
    matrices = read(BASE / "exact43/comparison/matched-matrices.npz")
    meta = pd.read_csv(record(BASE / "validation1020/cohort.csv"), dtype={"CID": str})
    zmap = meta.set_index("CID").zHEL.to_dict()
    units = None
    if args.modes == 12:
        expanded = json.loads(record(BASE / "expanded12/result.json").read_text())
        assert expanded["discovery_objects"] == 43 and expanded["validation_objects"] == 1020
        units = json.loads(record(BASE / "expanded12/calspec-unit-map.json").read_text())["bands"]
        posterior = read(BASE / "expanded12/posterior-models.npz")
        names = ["systematics_only", "systematics_plus_observer"]
        model_post = [(posterior[n + "_discovery_posterior_mean"],
                       posterior[n + "_discovery_posterior_covariance"]) for n in names]
    else:
        names = ["calibration_only", "calibration_plus_observer"]
        model_post = [(validation[n + "_discovery_posterior_mean"],
                       validation[n + "_discovery_posterior_covariance"]) for n in names]
    responses = {}
    ledger = []
    checks = []
    for cohort, data, directory in [("discovery", discovery, "shared43"),
                                     ("validation", validation, "shared1020")]:
        ids = data["CID"].astype(str)
        full = []
        for i, cid in enumerate(ids):
            prefix = cid + "__"
            nominal_dir = "exact43" if cohort == "discovery" else "validation1020"
            with read(BASE / nominal_dir / "objectives" / f"objective_{cid}.npz") as nominal:
                order = pd.DataFrame({"MJD": nominal["MJD"], "band": nominal["band"]}).sort_values(
                    ["MJD", "band"], kind="stable").index.to_numpy()
                flux = nominal["model_flux"][order]
                bands = nominal["band"][order]
                if cohort == "discovery":
                    jac = matrices[prefix + "jacobian_flux"].copy()
                    cov = matrices[prefix + "exact_covariance"]
                else:
                    with read(BASE / "validation1020/analysis/objects" / f"{cid}.npz") as cache:
                        jac = cache["jacobian_flux"].copy()
                        cov = cache["exact_covariance"]
                        assert np.array_equal(flux, cache["official_flux_model"])
                assert np.array_equal(cov, nominal["frozen_flux_covariance"][np.ix_(order, order)])
                jac[:, 0] = -K * flux
                chol = np.linalg.cholesky(cov)
                jw = solve_triangular(chol, jac, lower=True)
                u, singular, vt = np.linalg.svd(jw, full_matrices=True)
                assert np.sum(singular > singular[0] * 1e-10) == 4
                inverse = (vt.T / singular) @ u[:, :4].T
                columns = []
                for mode in range(1, 10):
                    with read(BASE / directory / f"model{mode:03d}/objectives/objective_{cid}.npz") as variant:
                        for key in ["MJD", "band", "data_flux", "data_fluxerr", "parameters_x0_x1_c_t0"]:
                            assert np.array_equal(variant[key], nominal[key]), (cid, mode, key)
                        columns.append(.3 * (variant["model_flux"][order] - flux))
                if args.modes == 12:
                    for name, scale in [("CALSPEC", 1.), ("MWEBV", 1.), ("COLORLAW", .3)]:
                        with read(BASE / f"expanded12/{cohort}/{name}/objectives/objective_{cid}.npz") as variant:
                            for key in ["MJD", "band", "parameters_x0_x1_c_t0"]:
                                assert np.array_equal(variant[key], nominal[key])
                            mapped = variant["model_flux"][order]
                            if name == "CALSPEC":
                                mapped = mapped / np.array([units[b]["flux_scale_float32"] for b in bands])
                            columns.append(scale * (mapped - flux))
                griz = np.column_stack([-K * flux * (bands == b) for b in "griz"])
                observer = np.column_stack([griz[:, 0] - griz[:, 1],
                                            griz[:, 2] - griz[:, 1],
                                            griz[:, 3] - griz[:, 1]])
                delta = np.column_stack([np.column_stack(columns), .02 * observer])
                whitened = solve_triangular(chol, delta, lower=True)
                dp = -inverse @ whitened
                independent = -np.linalg.solve(jw.T @ jw, jw.T @ whitened)
                assert np.allclose(dp, independent, rtol=1e-7, atol=1e-9), cid
                residual = whitened + jw @ dp
                orth = u[:, 4:] @ (u[:, 4:].T @ whitened)
                assert np.max(np.abs(residual - orth)) < 1e-8, cid
                gray = solve_triangular(chol, -K * flux, lower=True)
                gray_dp = -inverse @ gray
                assert np.allclose(gray_dp, [-1, 0, 0, 0], atol=1e-9, rtol=0), (cid, gray_dp)
                gray_projection = float(np.linalg.norm(u[:, 4:].T @ gray))
                assert gray_projection < 1e-8, cid
                response = STANDARDIZE @ dp
                full.append(response)
                for k in range(delta.shape[1]):
                    ledger.append({"cohort": cohort, "CID": cid, "mode": k + 1,
                                   "zHEL": float(nominal["zHEL"][0]),
                                   "delta_mB": float(dp[0, k]), "delta_x1": float(dp[1, k]),
                                   "delta_c": float(dp[2, k]), "delta_t0": float(dp[3, k]),
                                   "delta_standardized_mag": float(response[k]),
                                   "projected_flux_norm_squared": float(orth[:, k] @ orth[:, k]),
                                   "absorbed_flux_norm_squared": float((jw @ dp[:, k]) @ (jw @ dp[:, k]))})
                checks.append({"cohort": cohort, "CID": cid,
                               "normal_equation_max_error": float(np.max(np.abs(dp - independent))),
                               "decomposition_max_error": float(np.max(np.abs(residual - orth))),
                               "gray_coordinate_max_error": float(np.max(np.abs(gray_dp - [-1, 0, 0, 0]))),
                               "gray_projected_norm": gray_projection})
            if i % 200 == 0:
                print(cohort, i + 1, "of", len(ids), flush=True)
        responses[cohort] = np.array(full)
    ids = validation["CID"].astype(str)
    ordered = sorted(range(len(ids)), key=lambda i: (zmap[ids[i]], ids[i]))
    low, high = np.array(ordered[:255]), np.array(ordered[-255:])
    h = responses["validation"][high].mean(0) - responses["validation"][low].mean(0)
    contrasts = {}
    for name, (mean, covariance) in zip(names, model_post):
        vector = h[:len(mean)]
        prior_variance = float(vector @ vector)
        posterior_variance = float(vector @ covariance @ vector)
        object_sd = np.sqrt(np.einsum("ij,jk,ik->i", responses["validation"][:, :len(mean)],
                                    covariance, responses["validation"][:, :len(mean)]))
        contrasts[name] = {"mean_mag": float(vector @ mean),
                           "prior_sd_mag": float(np.sqrt(prior_variance)),
                           "discovery_conditioned_sd_mag": float(np.sqrt(posterior_variance)),
                           "remaining_variance_fraction": posterior_variance / prior_variance,
                           "validation_object_conditional_sd_quantiles": np.quantile(object_sd, [0, .5, .95, 1]).tolist(),
                           "contrast_response_per_unit_mode": vector.tolist()}
    check = pd.DataFrame(checks)
    result = {"scope": "Retrospective finite-mode local pre-BBC compensation diagnostic; not measured distance bias or a cosmology posterior.",
              "calibration_modes": args.modes, "discovery_objects": 43, "validation_objects": 1020,
              "coordinates": ["mB", "x1", "c", "t0"], "standardization_vector": STANDARDIZE.tolist(),
              "contrast": "Unweighted high255 minus low255 validation zHEL ranks; CID tie break; same target for every model.",
              "low_redshift_range": [zmap[ids[low[0]]], zmap[ids[low[-1]]]],
              "high_redshift_range": [zmap[ids[high[0]]], zmap[ids[high[-1]]]],
              "models": contrasts,
              "verification_maxima": {k: float(check[k].max()) for k in check.columns if k not in ["cohort", "CID"]},
              "limitations": ["Finite Gaussian mode prior is an approximation, not a bound on all calibration/dust/population errors.",
                              "Only discovery residuals condition modes; no validation coefficient fit.",
                              "Fixed nominal nuisance tangent and covariance; nonlinear transport, alpha/beta, selection and BBC not refitted.",
                              "Observer extension is an empirical mean family, not an identified physical calibration correction.",
                              "Gray luminosity/distance evolution lies outside these residual constraints."]}
    pd.DataFrame(ledger).to_csv(out / "object-mode-responses.csv", index=False)
    check.to_csv(out / "checks.csv", index=False)
    pd.DataFrame({"CID": ids, "zHEL": [zmap[c] for c in ids],
                  "low_quartile": np.isin(np.arange(len(ids)), low),
                  "high_quartile": np.isin(np.arange(len(ids)), high)}).to_csv(out / "contrast-membership.csv", index=False)
    np.savez_compressed(out / "responses.npz", discovery=responses["discovery"], validation=responses["validation"],
                        discovery_CID=discovery["CID"], validation_CID=validation["CID"], contrast=h)
    (out / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    outputs = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file()}
    (out / "manifest.json").write_text(json.dumps({"inputs_sha256": inputs, "outputs_sha256": outputs}, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
