"""Verify late-time nested outputs, sensitivity arithmetic and derived kinematics."""

import json
import re
from pathlib import Path
import numpy as np
from scipy.special import logsumexp
from late_nested import ROOT, RESULTS, WORK, quantities, qj, sha


def main():
    sources = ["late_nested.py", "late_nested_refine.py", "late_tail.py"]
    paths = [
        RESULTS / name
        for name in [
            "late-nested-design.json",
            "late-nested-927085.json",
            "late-nested-927086.json",
            "late-nested-comparison.json",
            "late-nested-refinement-design.json",
            "late-nested-refined-927088.json",
            "late-nested-refined-927089.json",
            "late-nested-refined-comparison.json",
            "late-tail-design.json",
            "late-tail-profile.json",
            "late-tail-conditional.json",
        ]
    ]
    records = {p.name: json.loads(p.read_text()) for p in paths}
    checked = {}
    for record in records.values():
        for key in [
            "input_sha256",
            "dependencies_sha256",
            "outputs_sha256",
            "source_records_sha256",
        ]:
            for name, expected in record.get(key, {}).items():
                actual = sha(ROOT / name)
                assert actual == expected, name
                checked[name] = actual
    for seed, source, prefix in [
        (927085, sources[0], "late-nested-"),
        (927086, sources[0], "late-nested-"),
        (927088, sources[1], "late-nested-refined-"),
        (927089, sources[1], "late-nested-refined-"),
    ]:
        assert records[f"{prefix}{seed}.json"]["code_sha256"] == sha(
            Path(__file__).with_name(source)
        )
    for name in [
        "late-tail-profile.json",
        "late-tail-conditional.json",
        "late-tail-design.json",
    ]:
        assert records[name]["code_sha256"] == sha(Path(__file__).with_name(sources[2]))

    # Independent derivative definitions: q=(1+z)E'/E-1 and j=q^2+(1+z)^2E''/E.
    rng = np.random.default_rng(927091)
    errors = []
    for _ in range(200):
        om, w, wa, z = (
            rng.uniform(0.01, 0.99),
            rng.uniform(-3, 1),
            rng.uniform(-5, 3),
            rng.uniform(0, 2.33),
        )

        def expansion(t):
            return np.sqrt(
                om * (1 + t) ** 3
                + (1 - om)
                * (1 + t) ** (3 * (1 + w + wa))
                * np.exp(-3 * wa * t / (1 + t))
            )

        h = 1e-4
        f = np.array([expansion(z + k * h) for k in [-2, -1, 0, 1, 2]])
        first = (f[0] - 8 * f[1] + 8 * f[3] - f[4]) / (12 * h)
        second = (-f[0] + 16 * f[1] - 30 * f[2] + 16 * f[3] - f[4]) / (12 * h * h)
        q = (1 + z) * first / f[2] - 1
        j = q * q + (1 + z) ** 2 * second / f[2]
        production = qj(np.array([[om, 10000.0, w, wa]]), z)
        errors.append([abs(q - production[0][0]), abs(j - production[1][0])])
    maximum = np.max(errors, axis=0)
    assert maximum[0] < 1e-8 and maximum[1] < 2e-6

    arrays = []
    conditional_errors = []
    for seed in [927088, 927089]:
        path = WORK / f"refined-seed-{seed}" / "samples.npz"
        a = np.load(path)
        x, w = a["samples"], a["weights"]
        arrays.append((x, w))
        assert np.max(abs(w - np.exp(a["logwt"] - logsumexp(a["logwt"])))) < 1e-14
        saved = records[f"late-nested-refined-{seed}.json"]
        assert (
            saved["single_run_diagnostic_passed"] and saved["status"] == "completed_run"
        )
        masks = {
            "none": np.ones(len(x), bool),
            "Omega_m>0.1": x[:, 0] > 0.1,
            "w0+wa<=0": x[:, 2] + x[:, 3] <= 0,
            "Omega_m>0.1 and w0+wa<=0": (x[:, 0] > 0.1) & (x[:, 2] + x[:, 3] <= 0),
        }
        for label, keep in masks.items():
            row = records["late-tail-conditional.json"]["seeds"][str(seed)][label]
            total = float(w[keep].sum())
            cw = w[keep] / total
            assert abs(total - row["retained_baseline_probability"]) < 1e-14
            assert abs(1 / (cw @ cw) - row["weighted_ESS"]) < 1e-8
            for name, v in quantities(x[keep]).items():
                mean = float(cw @ v)
                sd = float(np.sqrt(cw @ ((v - mean) ** 2)))
                conditional_errors.extend(
                    [
                        abs(mean - row["posterior"][name]["mean"]),
                        abs(sd - row["posterior"][name]["sd"]),
                    ]
                )
    assert max(conditional_errors) < 1e-8

    # Rebuild the cross-seed CDF statistic from weighted histograms on the joint grid.
    cdfs = {}
    qa, qb = [quantities(x) for x, _ in arrays]
    for name in qa:
        grid, inverse = np.unique(np.r_[qa[name], qb[name]], return_inverse=True)
        n = len(qa[name])
        aa = np.bincount(inverse[:n], weights=arrays[0][1], minlength=len(grid))
        bb = np.bincount(inverse[n:], weights=arrays[1][1], minlength=len(grid))
        cdfs[name] = float(np.max(abs(np.cumsum(aa) - np.cumsum(bb))))
    comparison = records["late-nested-refined-comparison.json"]
    for name, value in cdfs.items():
        assert abs(value - comparison["maximum_marginal_CDF_differences"][name]) < 1e-12
    seeds = [records[f"late-nested-refined-{s}.json"] for s in [927088, 927089]]
    evidence_difference = abs(
        seeds[0]["logZ_unnormalized"] - seeds[1]["logZ_unnormalized"]
    )
    evidence_error = float(np.hypot(*[r["logZ_error"] for r in seeds]))
    passed = max(cdfs.values()) <= 0.05 and evidence_difference <= 3 * evidence_error
    assert comparison["status"] == (
        "passed_independent_seed_diagnostics" if passed else "diagnostic_failed"
    )
    # A completed failed convergence gate is retained as a scientific limitation.
    assert records["late-nested-comparison.json"]["status"] == "diagnostic_failed"

    note = ROOT / "studies/unified_cosmology/notes/late-time-tails.md"
    links = []
    for link in re.findall(r"\]\(([^)]+)\)", note.read_text()):
        if "://" not in link:
            destination = (note.parent / link.split("#")[0]).resolve()
            assert destination.exists(), link
            links.append(link)
    record = {
        "status": "passed",
        "refined_sampling_diagnostic_status": comparison["status"],
        "initial_sampling_diagnostic_status": records["late-nested-comparison.json"][
            "status"
        ],
        "verified_distinct_identities": len(checked),
        "kinematic_derivative_probes": len(errors),
        "maximum_q_and_jerk_derivative_errors": maximum.tolist(),
        "maximum_conditional_mean_or_sd_error": max(conditional_errors),
        "independent_CDF_differences": cdfs,
        "resolved_note_links": len(links),
        "code_sha256": sha(__file__),
        "input_sha256": {
            str(p.relative_to(ROOT)): sha(p)
            for p in paths
            + [note, Path(__file__).with_name("nested-requirements-lock.txt")]
            + [Path(__file__).with_name(n) for n in sources]
        },
        "qualification": "Pass means identity and arithmetic checks agree, including retained failed diagnostics. It is not a proof of global mode coverage or physical adequacy.",
    }
    (RESULTS / "late-tail-validation.json").write_text(
        json.dumps(record, indent=2) + "\n"
    )
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
