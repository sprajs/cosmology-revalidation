"""Selected-sample DES predictor comparison with enforced convergence gates."""

import numpy as np
import jax
import jax.numpy as jnp
import numpyro
import numpyro.distributions as dist
from numpyro.diagnostics import summary

jax.config.update("jax_enable_x64", True)
FEATURES = [
    "x",
    "c",
    "host",
    "host_colour",
    "red_hinge",
    "x_evolution",
    "c_evolution",
    "x_squared",
]
PRIORS = {
    "x": (-0.15, 0.2),
    "c": (3.0, 2.0),
    "host": (0.0, 0.3),
    "host_colour": (0.0, 2.0),
    "red_hinge": (0.0, 2.0),
    "x_evolution": (0.0, 0.5),
    "c_evolution": (0.0, 3.0),
    "x_squared": (0.0, 0.1),
}
MODELS = {
    "none": [],
    "stretch": ["x"],
    "colour": ["c"],
    "tripp": ["x", "c"],
    "host": ["x", "c", "host"],
    "host_colour": ["x", "c", "host", "host_colour"],
    "broken_colour": ["x", "c", "host", "red_hinge"],
    "evolution": ["x", "c", "host", "x_evolution", "c_evolution"],
    "flexible": FEATURES,
}
KNOTS = jnp.array([0.01, 0.10, 0.20, 0.35, 0.50, 0.70, 0.90, 1.20])


def features(data, h):
    x = data["y"][:, 1]
    c = data["y"][:, 2]
    hc = h - 0.5
    g = data["z"] / (1 + data["z"]) - 0.2
    one = jnp.ones_like(x)
    zero = jnp.zeros_like(x)
    f = jnp.stack(
        [x, c, hc * one, hc * c, jnp.maximum(c, 0), x * g, c * g, x * x], axis=-1
    )
    dx = jnp.stack([one, zero, zero, zero, zero, g, zero, 2 * x], axis=-1)
    dc = jnp.stack(
        [zero, one, zero, hc * one, (c > 0).astype(float), zero, g, zero], axis=-1
    )
    return f, dx, dc


def components(params, data, name):
    indices = jnp.asarray([FEATURES.index(x) for x in MODELS[name]], dtype=int)
    baseline = (
        data["distance_reference"]
        + params["M"]
        + jnp.interp(
            data["z"], KNOTS, jnp.concatenate([jnp.zeros(1), params["offsets"]])
        )
    )
    coef = params.get("coefficients", jnp.empty(0))
    means = []
    variances = []
    for h in [0, 1]:
        f, dx, dc = features(data, h)
        mu = baseline + f[:, indices] @ coef
        # With the pinned JAX 0.6.2, negating a closure-constant indexed array
        # before dot gives incorrect compiled values/gradients for 2+ features.
        # Negating the dot results preserves the intended projection algebra.
        dmean_dx = dx[:, indices] @ coef
        dmean_dc = dc[:, indices] @ coef
        v = jnp.stack([jnp.ones(len(mu)), -dmean_dx, -dmean_dc], axis=-1)
        variance = (
            jnp.einsum("ni,nij,nj->n", v, data["cov"], v) + params["scatter"][h] ** 2
        )
        means.append(mu)
        variances.append(variance)
    return jnp.stack(means, axis=-1), jnp.stack(variances, axis=-1)


def predictive(params, data, name, noise):
    means, variances = components(params, data, name)
    law = (
        dist.Normal(means, jnp.sqrt(variances))
        if noise == "gaussian"
        else dist.StudentT(4.0, means, jnp.sqrt(variances / 2.0))
    )
    lps = law.log_prob(data["y"][:, 0, None])
    ph = jnp.clip(data["host_prob"], 1e-10, 1 - 1e-10)
    lp = jnp.logaddexp(jnp.log1p(-ph) + lps[:, 0], jnp.log(ph) + lps[:, 1])
    mean = (1 - ph) * means[:, 0] + ph * means[:, 1]
    var = (
        (1 - ph) * variances[:, 0]
        + ph * variances[:, 1]
        + ph * (1 - ph) * (means[:, 1] - means[:, 0]) ** 2
    )
    return lp, mean, var


def model(data, name, noise):
    p = {
        "M": numpyro.sample("M", dist.Normal(-19.3, 1.0)),
        "offsets": numpyro.sample(
            "offsets", dist.Normal(0.0, 0.5).expand([7]).to_event(1)
        ),
        "scatter": numpyro.sample(
            "scatter", dist.HalfNormal(0.7).expand([2]).to_event(1)
        ),
    }
    if MODELS[name]:
        prior = jnp.asarray([PRIORS[x] for x in MODELS[name]])
        p["coefficients"] = numpyro.sample(
            "coefficients", dist.Normal(prior[:, 0], prior[:, 1]).to_event(1)
        )
    numpyro.factor(
        "conditional_selected_magnitude", predictive(p, data, name, noise)[0].sum()
    )


def serial(x):
    if isinstance(x, dict):
        return {k: serial(v) for k, v in x.items()}
    if isinstance(x, (np.ndarray, jax.Array)):
        return np.asarray(x).tolist()
    if isinstance(x, np.generic):
        return x.item()
    return x


def convergence_diagnostics(samples, extra, max_tree_depth=10):
    """Reject chains before any held-out prediction is evaluated."""
    stats = summary(samples, group_by_chain=True)
    rhats = np.concatenate([np.asarray(v["r_hat"]).ravel() for v in stats.values()])
    neffs = np.concatenate([np.asarray(v["n_eff"]).ravel() for v in stats.values()])
    acceptance = np.asarray(extra["accept_prob"])
    steps = np.asarray(extra["num_steps"])
    divergences = int(np.asarray(extra["diverging"]).sum())
    reasons = []
    n_chains = next(iter(samples.values())).shape[0]
    if n_chains < 4:
        reasons.append("fewer than four chains")
    if not np.all(np.isfinite(rhats)) or np.max(rhats) > 1.01:
        reasons.append("nonfinite or R-hat > 1.01")
    if not np.all(np.isfinite(neffs)) or np.min(neffs) < 400:
        reasons.append("nonfinite or effective sample size < 400")
    if divergences:
        reasons.append("divergent transitions")
    if not np.all(np.isfinite(acceptance)) or float(acceptance.mean()) < 0.6:
        reasons.append("nonfinite or mean acceptance < 0.6")
    if np.any(steps >= 2**max_tree_depth - 1):
        reasons.append("maximum tree depth reached")
    for name, values in samples.items():
        if not np.all(np.isfinite(values)) or np.any(np.std(values, axis=1) == 0):
            reasons.append(f"{name} has a nonfinite or immobile chain")
    report = {
        "valid_for_scoring": not reasons,
        "reasons": reasons,
        "n_chains": n_chains,
        "max_r_hat": float(np.max(rhats)),
        "min_n_eff": float(np.min(neffs)),
        "divergences": divergences,
        "max_steps": int(np.max(steps)),
        "mean_acceptance": float(np.mean(acceptance)),
        "parameters": serial(stats),
    }
    return report
