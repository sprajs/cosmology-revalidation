#!/usr/bin/env python3
"""Plot verified constructed responses; infer no physical correction or interval."""
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "runs/research_2026_09_26/sed_nonlinear_validation"
OUT = ROOT / "runs/research_2026_09_26/sed_ambiguity_figure"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    paths = [SOURCE / "resolved/paired-responses.csv", SOURCE / "resolved/result.json",
             SOURCE / "contrast-membership.csv"]
    result = json.loads(paths[1].read_text())
    assert result["all_gates_pass"] and result["objects"] == 1020
    rows = list(csv.DictReader(paths[0].open()))
    members = list(csv.DictReader(paths[2].open()))
    ids = [r["CID"] for r in sorted(members, key=lambda r: (float(r["zHEL"]), r["CID"]))]
    membership = {r["CID"]: r for r in members}
    assert len(ids) == len(set(ids)) == 1020
    z = np.array([float(membership[c]["zHEL"]) for c in ids])
    low = np.array([membership[c]["low_quartile"] == "True" for c in ids])
    high = np.array([membership[c]["high_quartile"] == "True" for c in ids])
    assert low.sum() == high.sum() == 255 and not (low & high).any()
    series = {}
    contrasts = {}
    for target in ("native_noiseless", "observed_flux"):
        for mode in ("observer", "sed"):
            selected = [r for r in rows if r["target"] == target and r["mode"] == mode]
            data = {r["CID"]: r for r in selected}
            assert len(selected) == len(data) == 1020 and set(data) == set(ids)
            assert all(r["paired_gate"] == "True" for r in selected)
            values = np.array([float(data[c]["nonlinear_standardized_mag"]) for c in ids])
            assert np.isfinite(values).all()
            series[target, mode] = values
            contrasts[target, mode] = float(values[high].mean() - values[low].mean())

    OUT.mkdir(parents=True, exist_ok=True)
    # Ten equal-count redshift groups are descriptive, not uncertainty intervals.
    groups = np.array_split(np.arange(1020), 10)
    with (OUT / "plotted-data.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["CID", "zHEL", "low_quartile", "high_quartile", "redshift_group",
                         "observer_observed_mag", "sed_observed_mag",
                         "observer_noiseless_mag", "sed_noiseless_mag"])
        bins = {int(i): b + 1 for b, group in enumerate(groups) for i in group}
        for i, c in enumerate(ids):
            writer.writerow([c, z[i], bool(low[i]), bool(high[i]), bins[i],
                             *[series[t, m][i] for t in ("observed_flux", "native_noiseless")
                               for m in ("observer", "sed")]])

    plt.rcParams.update({"font.size": 11, "axes.spines.top": False,
                         "axes.spines.right": False, "font.family": "DejaVu Sans"})
    fig, axes = plt.subplots(1, 2, figsize=(13, 6), gridspec_kw={"width_ratios": [1.5, 1]})
    fig.subplots_adjust(top=.77, bottom=.26, wspace=.32, left=.085, right=.965)
    fig.suptitle("Similar fitted flux changes can imply different distance responses",
                 x=.085, y=.97, ha="left", fontsize=17, fontweight="bold")
    fig.text(.085, .902, "Constructed alternatives • 1,020 DES supernovae • 39,606 accepted epochs",
             fontsize=12)
    colors = {"observer": "#147D92", "sed": "#B44732"}
    labels = {"observer": "Observer-filter change", "sed": "Smooth SED change with redshift"}
    markers = {"observer": "o", "sed": "^"}
    ax = axes[0]
    for mode in ("observer", "sed"):
        values = series["observed_flux", mode]
        ax.scatter(z, values, s=7, alpha=.17, color=colors[mode], marker=markers[mode],
                   rasterized=True, linewidths=0)
        ax.plot([z[g].mean() for g in groups], [values[g].mean() for g in groups],
                color=colors[mode], marker=markers[mode], lw=1.8, ms=6, label=labels[mode])
    ax.axhline(0, color=".55", lw=.7)
    ax.set(xlabel="Heliocentric redshift", ylabel="Change in fixed-reference distance coordinate (mag)",
           title="Responses when refitting observed fluxes")
    ax.legend(loc="lower left", frameon=False, fontsize=9)
    ax.grid(axis="y", alpha=.15)
    ax = axes[1]
    for mode in ("observer", "sed"):
        vals = [contrasts[t, mode] for t in ("native_noiseless", "observed_flux")]
        ax.plot([0, 1], vals, marker=markers[mode], color=colors[mode], lw=1.8, ms=8)
        for i, value in enumerate(vals):
            ax.annotate(f"{value:+.5f}", (i, value), xytext=(0, 11 if mode == "observer" else -20),
                        textcoords="offset points", ha="center", color=colors[mode], fontsize=11)
    ax.axhline(0, color=".55", lw=.7)
    ax.set(xticks=[0, 1], xticklabels=["Noiseless target", "Observed fluxes"], xlim=(-.4, 1.4),
           ylim=(-.235, .11), ylabel="High 255 minus low 255 response (mag)",
           title="Same frozen redshift groups")
    ax.text(.5, -.074, "Alternative separation\n0.2466 / 0.2467 mag", ha="center", va="center",
            fontsize=11, bbox={"facecolor": "white", "edgecolor": "none", "alpha": .95})
    ax.grid(axis="y", alpha=.15)
    fig.text(.085, .135,
             "Hypothetical responses, not measured biases or a cosmological correction. Covariance, training and selection are fixed.\n"
             "Left: all objects; lines join 10 equal-count redshift means. Points and lines are not uncertainty intervals.\n"
             "The refitted response vectors differ by 1.90% in squared covariance-weighted norm on observed fluxes.\n"
             "Source: resolved nonlinear fits and exact contrast membership; physical attribution remains unresolved.",
             fontsize=10, ha="left", va="top", linespacing=1.5)
    for extension in ("png", "pdf"):
        fig.savefig(OUT / f"constructed-ambiguity.{extension}", dpi=180, facecolor="white")
    plt.close(fig)
    output = {
        "purpose": "Descriptive figure of two constructed alternatives, not inferred bias or uncertainty",
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in paths},
        "script_sha256": sha(Path(__file__)),
        "objects": 1020, "groups": 10, "objects_per_group": 102,
        "high_minus_low_contrasts_mag": {f"{t}__{m}": value for (t, m), value in contrasts.items()},
        "outputs_sha256": {p.name: sha(p) for p in OUT.iterdir() if p.suffix in ("csv", "png", "pdf")},
    }
    (OUT / "manifest.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output["high_minus_low_contrasts_mag"], indent=2))


if __name__ == "__main__":
    main()
