"""Frozen synthetic one-factor amplitude-recovery screen; no observed flux or SNANA."""
from __future__ import annotations

from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import math
import time
import warnings

import numpy as np
from scipy.integrate import quad
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import brentq, minimize_scalar
from scipy.special import erfcx, log_ndtr

ROOT = Path(__file__).resolve().parents[2]
DESIGN = ROOT / "runs/research_2026_09_26/raisin_profile_solver_review/onefactor-review/v2"
OUT = ROOT / "runs/research_2026_09_26/onefactor_recovery"
SCRIPT = Path(__file__).resolve()
PROTOCOL = DESIGN / "recovery-protocol.json"
PREPARED = OUT / "execution-protocol.json"
NEW_DRAWS = OUT / "draws.npz"
RESULTS = OUT / "draw-results.jsonl"
KERNEL_PATH = DESIGN / "recovery-kernel.py"

spec = importlib.util.spec_from_file_location("frozen_recovery_kernel", KERNEL_PATH)
kernel = importlib.util.module_from_spec(spec)
spec.loader.exec_module(kernel)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path: Path, data: object) -> None:
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def finite_float(x: float) -> float | None:
    return float(x) if np.isfinite(x) else None


def source_checks() -> dict:
    frozen = json.loads(PROTOCOL.read_text())
    assert sha(PROTOCOL) == "39e07714c9c9de2b61d1a3ec3378335e14e9d3c6822202da1c8543d9fffb4bf0"
    for path, expected in frozen["inputs_sha256"].items():
        assert sha(ROOT / path) == expected, path
    with np.load(DESIGN / "recovery-design.npz", allow_pickle=False) as z:
        C, d, v, h = (z[k].copy() for k in ("covariance", "diagonal_variance", "loading", "template"))
        truth = float(z["true_amplitude"])
    assert C.shape == (74, 74) and d.shape == v.shape == h.shape == (74,)
    assert truth == 1 and np.all(d > 0)
    assert np.allclose(C, np.diag(d) + np.outer(v, v), atol=1e-9, rtol=1e-12)
    assert np.array_equal(h, .5 * np.sqrt(np.diag(C)))
    return frozen


def prepare() -> None:
    assert not PREPARED.exists() and not NEW_DRAWS.exists()
    frozen = source_checks()
    with np.load(DESIGN / "recovery-design.npz", allow_pickle=False) as z:
        C, d, v, h = (z[k].copy() for k in ("covariance", "diagonal_variance", "loading", "template"))
    rng_old = np.random.Generator(np.random.PCG64(frozen["randomness"]["independent_archive_seed"]))
    old_z = rng_old.standard_normal()
    old_eps = rng_old.standard_normal(74)
    old_y = h + v * old_z + np.sqrt(d) * old_eps
    archive_mask = old_y > 0
    rng = np.random.Generator(np.random.PCG64(frozen["randomness"]["new_noise_seed"]))
    z = np.empty(128)
    eps = np.empty((128, 74))
    y = np.empty((128, 74))
    for i in range(128):
        z[i] = rng.standard_normal()
        eps[i] = rng.standard_normal(74)
        y[i] = h + v*z[i] + np.sqrt(d)*eps[i]
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(NEW_DRAWS, mean=h, diagonal_variance=d, loading=v,
                        latent=z, independent_noise=eps, flux=y, positive=y > 0,
                        archived_latent=np.array(old_z), archived_independent_noise=old_eps,
                        archived_flux=old_y, archived_positive=archive_mask)
    prep = {
        "state": "Executable and all 128 complete paired draws frozen before estimator outcomes",
        "frozen_protocol_sha256": sha(PROTOCOL),
        "executable_sha256": sha(SCRIPT), "kernel_sha256": sha(KERNEL_PATH),
        "design_sha256": sha(DESIGN / "recovery-design.npz"),
        "draws_sha256": sha(NEW_DRAWS),
        "n": 128, "coordinates": 74, "archive_positive": int(archive_mask.sum()),
        "amplitude_domain": [0, 4], "grid_sizes": [129, 257],
        "tie_loglikelihood_tolerance": 1e-10,
        "candidate_dedup_amplitude_tolerance": 2e-8,
        "resource_limit_seconds": 600,
        "details": "One worker. Each draw uses one PCG64 normal latent and 74 normals. The independent archive uses its own seed once; the same 74-coordinate new y is used by all five estimators. Adaptive reference uses mode-centered scipy.quad, independent of Hermite nodes. The complete draw arrays precede any fit.",
    }
    dump(PREPARED, prep)
    (OUT / "executed-source.py").write_bytes(SCRIPT.read_bytes())
    print(json.dumps({"execution_protocol_sha256": sha(PREPARED),
                      "draws_sha256": sha(NEW_DRAWS), "archive_positive": int(archive_mask.sum())}))


def load_frozen() -> tuple[dict, dict, dict]:
    frozen = source_checks()
    prep = json.loads(PREPARED.read_text())
    assert prep["executable_sha256"] == sha(SCRIPT) == sha(OUT / "executed-source.py")
    assert prep["draws_sha256"] == sha(NEW_DRAWS)
    with np.load(NEW_DRAWS, allow_pickle=False) as z:
        draw = {k: z[k].copy() for k in z.files}
    h, d, v, y = (draw[k] for k in ("mean", "diagonal_variance", "loading", "flux"))
    assert y.shape == (128, 74)
    assert np.array_equal(draw["positive"], y > 0)
    assert np.allclose(y, h[None, :] + draw["latent"][:, None]*v[None, :]
                       + draw["independent_noise"]*np.sqrt(d)[None, :], atol=1e-12, rtol=1e-12)
    return frozen, prep, draw


def gauss_fit(y: np.ndarray, h: np.ndarray, d: np.ndarray, v: np.ndarray,
              mask: np.ndarray) -> dict:
    if not mask.any():
        return {"status": "unidentified_empty_mask", "amplitude": None,
                "unconstrained_amplitude": None, "score_at_truth": 0.0,
                "boundary": None}
    x, t, dd, vv = y[mask], h[mask], d[mask], v[mask]
    pinv = 1 / (1 + np.sum(vv*vv/dd))
    def bilinear(a, b):
        return float(np.sum(a*b/dd) - pinv*np.sum(a*vv/dd)*np.sum(b*vv/dd))
    numerator = bilinear(t, x)
    denominator = bilinear(t, t)
    assert denominator > 0
    raw = numerator/denominator
    fit = float(np.clip(raw, 0, 4))
    return {"status": "identified", "amplitude": fit,
            "unconstrained_amplitude": float(raw),
            "score_at_truth": numerator-denominator,
            "information": denominator,
            "boundary": "zero" if fit == 0 else "upper" if fit == 4 else "interior",
            "best_loglikelihood": kernel.log_gaussian(x, fit*t, dd, vv)}


def adaptive_log_orthant(mean: np.ndarray, d: np.ndarray, v: np.ndarray,
                         signs: np.ndarray) -> tuple[float, float | None, float]:
    n = len(mean)
    if n == 0:
        return 0.0, None, 0.0
    if n == 1:
        return float(log_ndtr(signs[0]*mean[0]/np.sqrt(d[0]+v[0]**2))), None, 0.0
    if np.all(v == 0):
        return float(log_ndtr(signs*mean/np.sqrt(d)).sum()), 0.0, 0.0
    aa = signs*mean/np.sqrt(d)
    bb = signs*v/np.sqrt(d)
    def mills(t):
        t = np.asarray(t)
        result = np.empty_like(t)
        negative = t < 0
        result[negative] = np.sqrt(2/np.pi)/erfcx(-t[negative]/np.sqrt(2))
        result[~negative] = np.exp(-.5*t[~negative]**2-.5*np.log(2*np.pi)-log_ndtr(t[~negative]))
        return result
    def score(z):
        return float(-z + np.sum(bb*mills(aa+bb*z)))
    radius = 4.0
    for _ in range(60):
        if score(-radius) >= 0 and score(radius) <= 0:
            break
        radius *= 2
    else:
        raise ValueError("adaptive mode could not be bracketed")
    mode = brentq(score, -radius, radius, xtol=1e-13)
    def ell(z):
        return float(-.5*z*z-.5*np.log(2*np.pi)+np.sum(log_ndtr(aa+bb*z)))
    top = ell(mode)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        integral, error = quad(lambda u: np.exp(ell(mode+u)-top), -12, 12,
                               epsabs=2e-13, epsrel=2e-13, limit=150, points=[0])
    assert integral > 0
    return float(top+np.log(integral)), float(mode), float(error/integral)


def reference_gaussian(y: np.ndarray, h: np.ndarray, d: np.ndarray,
                       v: np.ndarray, mask: np.ndarray):
    x, t, dd, vv = y[mask], h[mask], d[mask], v[mask]
    if len(x) == 0:
        return lambda amp: 0.0
    C = np.diag(dd) + np.outer(vv, vv)
    factor = cho_factor(C, lower=True, check_finite=True)
    logdet = 2*np.log(np.diag(factor[0])).sum()
    constant = -0.5*(len(x)*np.log(2*np.pi)+logdet)
    def loglike(amp):
        residual = x-amp*t
        return float(constant-.5*residual @ cho_solve(factor, residual))
    return loglike


def reference_proper(kind: str, amp: float, y: np.ndarray, h: np.ndarray,
                     d: np.ndarray, v: np.ndarray, positive: np.ndarray,
                     ref_gauss) -> tuple[float, dict]:
    neg = ~positive
    residual = y[positive] - amp*h[positive]
    vu = 1/(1 + np.sum(v[positive]**2/d[positive]))
    mu = vu*np.sum(v[positive]*residual/d[positive])
    lp_neg, mode_neg, err_neg = adaptive_log_orthant(
        amp*h[neg] + v[neg]*mu, d[neg], v[neg]*np.sqrt(vu), -np.ones(neg.sum()))
    joint = ref_gauss(amp) + lp_neg
    details = {"censor_latent_mode": mode_neg, "censor_adaptive_relative_error": err_neg}
    if kind == "joint_censored":
        return joint, details
    signs = np.where(positive, 1.0, -1.0)
    lp_pattern, mode_pattern, err_pattern = adaptive_log_orthant(amp*h, d, v, signs)
    details.update({"pattern_latent_mode": mode_pattern,
                    "pattern_adaptive_relative_error": err_pattern,
                    "full74_log_pattern": lp_pattern})
    return joint-lp_pattern, details


def proper_likelihood(kind: str, y: np.ndarray, h: np.ndarray, d: np.ndarray,
                      v: np.ndarray, positive: np.ndarray):
    neg_signs = -np.ones((~positive).sum())
    full_signs = np.where(positive, 1.0, -1.0)
    cache = {}
    def evaluate(amp: float, order: int = 64) -> float:
        key = (float(amp), order)
        if key in cache:
            return cache[key]
        mean = amp*h
        value = kernel.log_censored(y[positive], mean, d, v, positive, neg_signs, order)
        if kind == "full_pattern_conditioned":
            value -= kernel.log_orthant(mean, d, v, full_signs, order)
        cache[key] = float(value)
        return cache[key]
    return evaluate


def local_candidates(grid: np.ndarray, values: np.ndarray, evaluate, xtol: float) -> list[dict]:
    # Always retain endpoints and polish their neighboring intervals.
    n = len(grid)
    indices = [i for i in range(1, n-1) if values[i] >= values[i-1] and values[i] >= values[i+1]]
    intervals = [(grid[0], grid[1]), (grid[-2], grid[-1])]
    intervals.extend((grid[i-1], grid[i+1]) for i in indices)
    candidates = [dict(amplitude=float(grid[0]), loglikelihood=float(values[0]), origin="left_endpoint"),
                  dict(amplitude=float(grid[-1]), loglikelihood=float(values[-1]), origin="right_endpoint")]
    for left, right in intervals:
        solution = minimize_scalar(lambda a: -evaluate(float(a)), bounds=(left, right),
                                   method="bounded", options={"xatol": xtol, "maxiter": 500})
        if not solution.success or not np.isfinite(solution.fun):
            raise RuntimeError(f"polish failure in [{left},{right}]: {solution.message}")
        candidates.append(dict(amplitude=float(solution.x),
                               loglikelihood=float(-solution.fun), origin="polished_interval"))
    return candidates


def best_candidates(candidates: list[dict], tie_tol: float) -> tuple[dict, list[dict]]:
    distinct = []
    for candidate in sorted(candidates, key=lambda x: -x["loglikelihood"]):
        if all(abs(candidate["amplitude"]-old["amplitude"]) > 2e-8 for old in distinct):
            distinct.append(candidate)
    top = max(x["loglikelihood"] for x in distinct)
    ties = [x for x in distinct if top-x["loglikelihood"] <= tie_tol]
    # Reproducible convention only; distinct ties are explicitly reported.
    chosen = min(ties, key=lambda x: x["amplitude"])
    return chosen, ties


def numerical_fit(kind: str, y: np.ndarray, h: np.ndarray, d: np.ndarray,
                  v: np.ndarray, positive: np.ndarray, protocol: dict) -> dict:
    if kind == "full_pattern_conditioned" and not positive.any():
        return {"status": "unidentified_all_negative", "amplitude": None,
                "score_at_truth": 0.0, "boundary": None, "candidates": []}
    evaluate = proper_likelihood(kind, y, h, d, v, positive)
    xtol = protocol["optimization"]["xtol"]
    coarse_x = np.linspace(0, 4, 129)
    fine_x = np.linspace(0, 4, 257)
    try:
        fine_values = np.array([evaluate(float(a), 64) for a in fine_x])
        coarse_values = fine_values[::2]
        coarse_candidates = local_candidates(coarse_x, coarse_values, evaluate, xtol)
        fine_candidates = local_candidates(fine_x, fine_values, evaluate, xtol)
        coarse_best, coarse_ties = best_candidates(coarse_candidates, 1e-10)
        fine_best, fine_ties = best_candidates(fine_candidates, 1e-10)
        best_a = fine_best["amplitude"]
        coarse_log_gap = abs(fine_best["loglikelihood"]-coarse_best["loglikelihood"])
        coarse_a_gap = abs(best_a-coarse_best["amplitude"])
        step = protocol["score_step"]["central"]
        half = protocol["score_step"]["half"]
        score1 = (evaluate(1+step, 64)-evaluate(1-step, 64))/(2*step)
        score2 = (evaluate(1+half, 64)-evaluate(1-half, 64))/(2*half)
        points = [1.0] + [c["amplitude"] for c in coarse_candidates+fine_candidates]
        unique_points = []
        for amp in points:
            if all(abs(amp-other) > 2e-8 for other in unique_points):
                unique_points.append(amp)
        ref_gauss = reference_gaussian(y, h, d, v, positive)
        checked = []
        max_reference_gap = 0.0
        max_reference_error = 0.0
        unsupported = []
        for amp in unique_points:
            try:
                order_vals = {str(order): evaluate(amp, order) for order in (64, 128, 256)}
                reference, details = reference_proper(kind, amp, y, h, d, v, positive, ref_gauss)
                gap = max(abs(value-reference) for value in order_vals.values())
                max_reference_gap = max(max_reference_gap, gap)
                max_reference_error = max(max_reference_error,
                    details["censor_adaptive_relative_error"],
                    details.get("pattern_adaptive_relative_error", 0.0))
                checked.append({"amplitude": amp, "hermite": order_vals,
                                "independent_adaptive": reference,
                                "max_abs_gap": gap, **details})
            except (ValueError, AssertionError, RuntimeError, Warning) as exc:
                unsupported.append({"amplitude": amp, "reason": str(exc)})
        gates = {
            "coarse_fine_loglikelihood": coarse_log_gap <= protocol["optimization"]["fine_coarse_loglikelihood_tolerance"],
            "coarse_fine_amplitude": coarse_a_gap <= protocol["optimization"]["fine_coarse_amplitude_tolerance"],
            "score_half_step": abs(score1-score2) <= protocol["score_step"]["absolute_derivative_agreement"],
            "quadrature_and_adaptive": not unsupported and max_reference_gap <= protocol["optimization"]["max_loglikelihood_quadrature_error"],
        }
        return {
            "status": "identified" if all(gates.values()) else "gate_failure",
            "amplitude": best_a,
            "best_loglikelihood": fine_best["loglikelihood"],
            "boundary": "zero" if best_a <= 1e-8 else "upper" if best_a >= 4-1e-8 else "interior",
            "coarse_best": coarse_best,
            "coarse_fine_amplitude_gap": coarse_a_gap,
            "coarse_fine_loglikelihood_gap": coarse_log_gap,
            "coarse_candidate_count": len(coarse_candidates),
            "fine_candidate_count": len(fine_candidates),
            "coarse_ties": coarse_ties, "fine_ties": fine_ties,
            "score_at_truth": float(score1), "score_at_truth_half_step": float(score2),
            "score_step_gap": abs(score1-score2),
            "gates": gates,
            "candidate_checks": checked,
            "candidate_unsupported": unsupported,
            "max_hermite_vs_adaptive_gap": max_reference_gap,
            "max_adaptive_relative_error_estimate": max_reference_error,
        }
    except (ValueError, AssertionError, RuntimeError, Warning) as exc:
        return {"status": "optimization_failure", "amplitude": None,
                "failure": repr(exc), "boundary": None}


def score_draw(i: int, draw: dict, protocol: dict) -> dict:
    y = draw["flux"][i]
    h, d, v = (draw[k] for k in ("mean", "diagonal_variance", "loading"))
    positive = draw["positive"][i]
    archive = draw["archived_positive"]
    assert np.array_equal(positive, y > 0)
    estimators = {
        "signed_gaussian": gauss_fit(y, h, d, v, np.ones(74, bool)),
        "naive_positive_gaussian": gauss_fit(y, h, d, v, positive),
        "archive_mask_gaussian": gauss_fit(y, h, d, v, archive),
        "joint_censored": numerical_fit("joint_censored", y, h, d, v, positive, protocol),
        "full_pattern_conditioned": numerical_fit("full_pattern_conditioned", y, h, d, v, positive, protocol),
    }
    for result in estimators.values():
        a = result.get("amplitude")
        result["delta_log_distance"] = -2.5*math.log10(a) if a is not None and a > 0 else None
    return {"draw": i, "latent": float(draw["latent"][i]),
            "positive_count": int(positive.sum()),
            "positive_mask": positive.astype(int).tolist(),
            "archive_positive_count": int(archive.sum()),
            "estimators": estimators}


def read_existing() -> list[dict]:
    if not RESULTS.exists():
        return []
    rows = [json.loads(line) for line in RESULTS.read_text().splitlines()]
    assert [x["draw"] for x in rows] == list(range(len(rows)))
    return rows


def run(phase: str) -> None:
    protocol, prep, draw = load_frozen()
    rows = read_existing()
    if phase == "benchmark":
        assert not rows
        target = 8
    else:
        assert len(rows) == 8 and (OUT / "benchmark-gate.json").exists()
        gate = json.loads((OUT / "benchmark-gate.json").read_text())
        assert gate["pass"] and gate["forecast_total_seconds"] <= protocol["resource"]["maximum_elapsed_seconds"]
        target = 128
    start = time.monotonic()
    with RESULTS.open("a") as stream:
        for i in range(len(rows), target):
            result = score_draw(i, draw, protocol)
            result["elapsed_seconds_from_phase_start"] = time.monotonic()-start
            stream.write(json.dumps(result, allow_nan=False) + "\n")
            stream.flush()
            if phase == "full" and time.monotonic()-start + float(json.loads((OUT / "benchmark-gate.json").read_text())["elapsed_seconds"]) > 600:
                break
    all_rows = read_existing()
    if phase == "benchmark":
        elapsed = time.monotonic()-start
        forecast = elapsed * 128/8
        failures = [(r["draw"], name, value["status"]) for r in all_rows
                    for name, value in r["estimators"].items()
                    if value["status"] in ("gate_failure", "optimization_failure")]
        benchmark = {"pass": not failures and forecast <= 600,
                     "elapsed_seconds": elapsed, "forecast_total_seconds": forecast,
                     "draws_scored": len(all_rows), "estimator_failures": failures,
                     "protocol_sha256": sha(PROTOCOL),
                     "execution_protocol_sha256": sha(PREPARED),
                     "draws_sha256": sha(NEW_DRAWS),
                     "results_sha256": sha(RESULTS)}
        dump(OUT / "benchmark-gate.json", benchmark)
        print(json.dumps(benchmark, indent=2))
    else:
        state = {"draws_scored": len(all_rows), "complete": len(all_rows) == 128,
                 "elapsed_seconds_this_phase": time.monotonic()-start,
                 "results_sha256": sha(RESULTS),
                 "protocol_sha256": sha(PROTOCOL),
                 "execution_protocol_sha256": sha(PREPARED)}
        dump(OUT / "full-run-state.json", state)
        print(json.dumps(state, indent=2))


def mcse(values: list[float]) -> float | None:
    if len(values) <= 1:
        return None
    return float(np.std(values, ddof=1)/np.sqrt(len(values)))


def summarize() -> None:
    protocol, prep, draw = load_frozen()
    rows = read_existing()
    names = ("signed_gaussian", "naive_positive_gaussian", "archive_mask_gaussian",
             "joint_censored", "full_pattern_conditioned")
    report = {}
    for name in names:
        completed = [r["estimators"][name] for r in rows if r["estimators"][name]["status"] == "identified"]
        amplitudes = [x["amplitude"] for x in completed]
        errors = [a-1 for a in amplitudes]
        scores = [x["score_at_truth"] for x in completed]
        finite_d = [x["delta_log_distance"] for x in completed if x["delta_log_distance"] is not None]
        report[name] = {
            "identified": len(completed),
            "statuses": {status: sum(r["estimators"][name]["status"] == status for r in rows)
                         for status in sorted({r["estimators"][name]["status"] for r in rows})},
            "amplitude_bias": float(np.mean(errors)) if errors else None,
            "amplitude_bias_mcse": mcse(errors),
            "amplitude_rmse": float(np.sqrt(np.mean(np.square(errors)))) if errors else None,
            "amplitude_rmse_mcse_delta": (mcse([e*e for e in errors])/(2*np.sqrt(np.mean(np.square(errors))))) if errors and np.mean(np.square(errors)) > 0 else None,
            "score_mean_at_truth": float(np.mean(scores)) if scores else None,
            "score_mcse": mcse(scores),
            "score_mean_over_mcse": float(np.mean(scores)/mcse(scores)) if scores and mcse(scores) not in (None, 0) else None,
            "zero_boundary": sum(x["boundary"] == "zero" for x in completed),
            "upper_boundary": sum(x["boundary"] == "upper" for x in completed),
            "tied_fine_candidates": sum(len(x.get("fine_ties", [])) > 1 for x in completed),
            "log_distance_finite_count": len(finite_d),
            "log_distance_undefined_zero_count": sum(x["amplitude"] == 0 for x in completed),
            "log_distance_mean_finite_only": float(np.mean(finite_d)) if finite_d else None,
            "log_distance_mcse_finite_only": mcse(finite_d),
        }
        if name.endswith("gaussian"):
            raw = [x["unconstrained_amplitude"] for x in completed]
            report[name]["unconstrained_amplitude_bias"] = float(np.mean(np.array(raw)-1)) if raw else None
            report[name]["unconstrained_amplitude_bias_mcse"] = mcse([a-1 for a in raw])
            report[name]["unconstrained_outside_domain"] = sum(a < 0 or a > 4 for a in raw)
    paired = {}
    for name in names[1:]:
        pairs = [(r["estimators"][name]["amplitude"],
                  r["estimators"]["signed_gaussian"]["amplitude"])
                 for r in rows if r["estimators"][name]["status"] == "identified"
                 and r["estimators"]["signed_gaussian"]["status"] == "identified"]
        diffs = [a-b for a,b in pairs]
        paired[name+"_minus_signed"] = {"n": len(diffs),
                                        "mean": float(np.mean(diffs)) if diffs else None,
                                        "mcse": mcse(diffs)}
    max_kernel_gap = max((x.get("max_hermite_vs_adaptive_gap", 0)
                          for r in rows for x in r["estimators"].values()), default=0)
    max_score_gap = max((x.get("score_step_gap", 0)
                         for r in rows for x in r["estimators"].values()), default=0)
    all_failures = [{"draw": r["draw"], "estimator": name, "status": value["status"],
                     "failure": value.get("failure"), "gates": value.get("gates"),
                     "unsupported": value.get("candidate_unsupported")}
                    for r in rows for name, value in r["estimators"].items()
                    if value["status"] not in ("identified", "unidentified_all_negative", "unidentified_empty_mask")]
    result = {"scope": "Synthetic conditional one-factor recovery; no observed outcomes or physical correction",
              "draws_scored": len(rows), "complete": len(rows) == 128,
              "positive_count_range": [min(r["positive_count"] for r in rows),
                                       max(r["positive_count"] for r in rows)] if rows else None,
              "archive_positive_count": int(draw["archived_positive"].sum()),
              "estimators": report, "paired_amplitude_changes": paired,
              "max_hermite_vs_independent_adaptive_gap": max_kernel_gap,
              "max_score_halfstep_gap": max_score_gap,
              "failures": all_failures,
              "frozen_protocol_sha256": sha(PROTOCOL),
              "execution_protocol_sha256": sha(PREPARED),
              "draws_sha256": sha(NEW_DRAWS), "results_sha256": sha(RESULTS),
              "source_sha256": sha(SCRIPT)}
    dump(OUT / "summary.json", result)
    print(json.dumps({"draws_scored": len(rows), "complete": len(rows)==128,
                      "estimator_biases": {k:v["amplitude_bias"] for k,v in report.items()},
                      "failures": len(all_failures), "max_adaptive_gap": max_kernel_gap}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=["prepare", "benchmark", "full", "summarize"])
    args = parser.parse_args()
    {"prepare": prepare, "benchmark": lambda: run("benchmark"),
     "full": lambda: run("full"), "summarize": summarize}[args.action]()
