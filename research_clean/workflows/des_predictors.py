"""Fit selected-sample predictors; held-out scores require convergence."""

import numpy as np
import pandas as pd
import jax
import jax.numpy as jnp
from jax.scipy.special import logsumexp
from numpyro.infer import MCMC, NUTS, init_to_median
from lib.paths import DATA
from lib.predictors import model, predictive, convergence_diagnostics, MODELS
from lib.cosmology import mu
from lib.records import write_json
import json

DEFAULTS = {
    "model": "tripp",
    "noise": "gaussian",
    "warmup": 1200,
    "draws": 2000,
    "chains": 4,
    "seed": 2026092160,
    "reference": "coasting",
    "covariance_scale": 1.0,
}


def run(out, cfg):
    if cfg["model"] not in MODELS or cfg["noise"] not in {"gaussian", "student4"}:
        raise ValueError("Unknown predictor model/noise law")
    if (
        cfg["chains"] < 4
        or cfg["draws"] < 400
        or cfg["warmup"] < 100
        or cfg["covariance_scale"] <= 0
    ):
        raise ValueError(
            "At least four chains, 400 draws and 100 warmup steps are required"
        )
    refs = {"coasting": 0.0, "eds": 0.5, "desitter": -1.0}
    if cfg["reference"] not in refs:
        raise ValueError("Unknown reference distance")
    with np.load(DATA / "des/predictors/data.npz", allow_pickle=False) as source:
        arrays = dict(source)
    rows = pd.read_csv(DATA / "des/predictors/rows.csv", dtype={"CID": str})
    cohort = json.loads((DATA / "des/predictors/cohort.json").read_text())
    keep = (arrays["survey"] == 10) & (arrays["pIa"] > 0.999)
    actual = dict(zip(rows.loc[keep, "CID"], map(int, arrays["fold"][keep])))
    if actual != cohort["cid_to_fold"]:
        raise ValueError("Frozen selected cohort or fold membership changed")
    arrays["cov"] *= cfg["covariance_scale"]
    train, test = keep & (arrays["fold"] != 0), keep & (arrays["fold"] == 0)

    def subset(mask):
        data = {k: jnp.asarray(v[mask]) for k, v in arrays.items()}
        data["distance_reference"] = jnp.asarray(
            mu(
                np.asarray(data["z"]),
                [refs[cfg["reference"]], 0.0],
                "kinematic",
                np.asarray(data["zhel"]),
            )[0]
            + 5 * np.log10(299792.458 / 70)
            + 25
        )
        return data

    dtrain, dtest = subset(train), subset(test)
    sampler = MCMC(
        NUTS(
            model,
            target_accept_prob=0.9,
            max_tree_depth=10,
            dense_mass=True,
            init_strategy=init_to_median(num_samples=20),
        ),
        num_warmup=cfg["warmup"],
        num_samples=cfg["draws"],
        num_chains=cfg["chains"],
        chain_method="sequential",
        progress_bar=False,
    )
    sampler.run(
        jax.random.PRNGKey(cfg["seed"]),
        dtrain,
        cfg["model"],
        cfg["noise"],
        extra_fields=("diverging", "num_steps", "accept_prob"),
    )
    samples = sampler.get_samples(group_by_chain=True)
    np.savez_compressed(
        out / "chains.npz", **{k: np.asarray(v) for k, v in samples.items()}
    )
    gate = convergence_diagnostics(
        samples, sampler.get_extra_fields(group_by_chain=True)
    )
    write_json(out / "convergence.json", gate)
    if not gate["valid_for_scoring"]:
        raise RuntimeError(
            "Convergence failed before held-out scoring: " + "; ".join(gate["reasons"])
        )
    ll, means, variances = jax.jit(
        jax.vmap(lambda p: predictive(p, dtest, cfg["model"], cfg["noise"]))
    )(sampler.get_samples())
    scores = np.asarray(logsumexp(ll, axis=0) - jnp.log(ll.shape[0]))
    prediction = np.asarray(means.mean(axis=0))
    variance = np.asarray(variances.mean(axis=0) + means.var(axis=0))
    table = rows.loc[test, ["CID", "IDSURVEY", "field", "depth", "fold"]].copy()
    table["observed_mB"] = np.asarray(dtest["y"][:, 0])
    table["predictive_mean"] = prediction
    table["predictive_sd"] = np.sqrt(variance)
    table["log_predictive_density"] = scores
    table.to_csv(out / "heldout.csv", index=False)
    return {
        "n_train": int(train.sum()),
        "n_test": int(test.sum()),
        "heldout_log_score_sum": float(scores.sum()),
        "heldout_rmse": float(np.sqrt(np.mean((table.observed_mB - prediction) ** 2))),
        "convergence": gate,
        "scope": "Empirical apparent-magnitude prediction in a frozen high-purity selected DES cohort. Flexible distance offsets and model-specific predictors are conditional descriptions; no parent-population or selection-free cosmological inference.",
    }
