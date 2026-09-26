"""Fixed descriptive two-panel view of already saved secondary pair ledgers."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np

HERE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    protocol = json.loads((HERE / "plot-protocol.json").read_text())
    assert sha(HERE / "signed_pair_ledger.csv") == protocol["pair_ledger_sha256"]
    assert sha(Path(__file__)) == protocol["script_sha256"]
    with (HERE / "signed_pair_ledger.csv").open(newline="") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 340 and {r["branch"] for r in rows} == {"secondary"}
    stage_a = json.loads((HERE / "stage_a.json").read_text())
    xlo, xhi, ylo, yhi = stage_a["bounds_east_north_arcsec"]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.9), sharex=True, sharey=True,
                             constrained_layout=True)
    artist = None
    for ax, visit in zip(axes, ("search", "template")):
        data = [r for r in rows if r["visit"] == visit]
        assert len(data) == 170
        east = np.array([float(r["east"]) for r in data])
        north = np.array([float(r["north"]) for r in data])
        v = np.array([float(r["V"]) for r in data])
        z = np.array([float(r["d"]) / np.sqrt(float(r["V"])) for r in data])
        # Fixed visual encodings. Marker area conveys quoted contrast V;
        # neither limits nor size mapping are fitted to the observed pattern.
        area = 35 + 50 * np.sqrt(v)
        artist = ax.scatter(east, north, c=z, s=area, cmap="RdBu_r", vmin=-3, vmax=3,
                            linewidths=.35, edgecolors="#1a1a1a", alpha=.91)
        ax.add_patch(Circle((0, 0), 5, fill=False, ec="#777777", lw=.85, ls="--"))
        ax.set_title(f"{visit.capitalize()} visit · 170 masked-secondary positions")
        ax.set_xlabel("East offset (arcsec)")
        ax.set_xlim(xlo, xhi)
        ax.set_ylim(ylo, yhi)
        ax.set_aspect("equal", adjustable="box")
        ax.grid(alpha=.16, lw=.5)
    axes[0].set_ylabel("North offset (arcsec)")
    cbar = fig.colorbar(artist, ax=axes, extend="both", shrink=.87, pad=.02)
    cbar.set_label("Signed repeat contrast z = (F₂−F₄)/√(V₂+V₄)")
    fig.suptitle("Current WFC3/IR F160W repeats · fixed mask and operator", fontsize=12)
    fig.text(.5, -.015, "Point area = 35 + 50√V points²; dashed circle is the frozen 5″ SN exclusion. "
             "Strict DQ0 branch stopped (n=0).", ha="center", fontsize=9)
    out = HERE / "secondary-spatial-signed-z.png"
    fig.savefig(out, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(json.dumps({"plot":str(out), "sha256":sha(out), "bytes":out.stat().st_size}))


if __name__ == "__main__":
    main()
