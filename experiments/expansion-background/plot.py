# %% Load a saved, accepted execution. Compatible with Jupytext percent format.
"""Plot Irreducible outputs; no scientific model is reimplemented here."""
import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def plot(run):
    run = Path(run).resolve()
    record = json.loads((run / "run.json").read_text())
    with (run / "stdout.json").open("rb") as stream:
        if hashlib.file_digest(stream, "sha256").hexdigest() != record["output_sha256"]["stdout.json"]:
            raise ValueError("Saved engine output changed after execution")
    reply = json.loads((run / "stdout.json").read_text())
    if record["execution"] != "completed" or reply["receipt"].get("accepted") is not True:
        raise ValueError("Plot requires a completed, accepted engine execution")
    if record["experiment"] != "expansion-background":
        raise ValueError("Wrong experiment for this plot")
    rows = reply["result"]["evaluations"]
    z, expansion, distance = [], [], []
    for row in rows:
        groups = row["groups"]
        for name in ("expansion", "physical"):
            if groups[name]["availability"] != "available" or groups[name]["status"] != "ok":
                raise ValueError("Unavailable engine result; preserve the failed row")
        z.append(row["source"]["z_expansion"])
        expansion.append(groups["expansion"]["value"]["H_km_s_mpc"]["value"])
        distance.append(groups["physical"]["value"]["luminosity_mpc"])
    # %% Render only reported values, retaining the request order.
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6), layout="constrained")
    axes[0].plot(z, expansion, "o-")
    axes[0].set_ylabel("H(z) [km/s/Mpc]")
    axes[1].plot(z, distance, "o-")
    axes[1].set_ylabel("Luminosity distance [Mpc]")
    for ax in axes:
        ax.set_xlabel("Redshift z")
        ax.grid(alpha=0.2)
    model = reply["result"]["models"][0]["source"]
    h0 = rows[0]["source"]["physical_scale"]["h0_km_s_mpc"]
    fig.suptitle(f"Flat matter + Λ control · Ωm = {model['omega_m']:g}, H₀ = {h0:g} km/s/Mpc")
    output = run / "figures"
    output.mkdir(exist_ok=True)
    fig.savefig(output / "background.png", dpi=160)
    fig.savefig(output / "background.svg")
    plt.close(fig)
    return output / "background.png"


# %% Command-line entry point. In a notebook, import plot and pass a saved run.
if __name__ == "__main__" and "__file__" in globals():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    print(plot(parser.parse_args().run))
