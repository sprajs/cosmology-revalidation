"""Cosmic star formation convolved with a specified delay-time law."""

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import PchipInterpolator

HUBBLE_GYR = 977.7922216807892
COSMOS = {
    "lcdm_H70": dict(H0=70.0, Om0=0.3, w0=-1.0, wa=0.0),
    "son_cpl_H70": dict(H0=70.0, Om0=0.353, w0=-0.42, wa=-1.75),
    "son_cpl_H63p6": dict(H0=63.6, Om0=0.353, w0=-0.42, wa=-1.75),
}
DTDS = {
    "C14_smooth": (0.3, -1.0, 20.0),
    "cut300": (0.3, -1.0, None),
    "W26_cut40": (0.04, -1.13, None),
}


def csfh(z, form):
    if form == "B13":
        return 0.180 * np.exp(
            -np.logaddexp(
                np.log(10) * (-0.997) * (z - 1.243), np.log(10) * 0.241 * (z - 1.243)
            )
        )
    if form == "MD14":
        return 0.015 * (1 + z) ** 2.7 / (1 + ((1 + z) / 2.9) ** 5.6)
    raise ValueError(form)


def clock(c, n=60001):
    # Integrate dt/d(ln a)=1/H; matter-era boundary is negligible at a=e^-16.
    la = np.linspace(-16.0, 0.0, n)
    a = np.exp(la)
    de = (
        (1 - c["Om0"])
        * a ** (-3 * (1 + c["w0"] + c["wa"]))
        * np.exp(-3 * c["wa"] * (1 - a))
    )
    einv = 1 / np.sqrt(c["Om0"] * a**-3 + de)
    age = (
        (
            2 * a[0] ** 1.5 / (3 * np.sqrt(c["Om0"]))
            + cumulative_trapezoid(einv, la, initial=0)
        )
        * HUBBLE_GYR
        / c["H0"]
    )
    return PchipInterpolator(la, age), PchipInterpolator(
        age, 1 / a - 1, extrapolate=False
    )


def calculate(c, sfh, dtd, n=12001, clock_n=60001):
    age_of_lna, z_of_age = clock(c, clock_n)
    tp, s, alpha = DTDS[dtd]
    rows = []
    for z in np.linspace(0, 2.5, 251):
        t = float(age_of_lna(-np.log1p(z)))
        # Separate continuous domain for truncated models avoids grid misplacement at cutoff.
        lo = 1e-7 if alpha is not None else tp
        delay = np.geomspace(lo, t * (1 - 1e-9), n)
        formed_time = t - delay
        sf = csfh(z_of_age(formed_time), sfh)
        sf = np.nan_to_num(
            sf
        )  # negligible tail earlier than first cosmic-time grid point
        if alpha is not None:
            x = np.log(delay / tp)
            phi = np.exp(alpha * x - np.logaddexp((alpha - s) * x, 0))
        else:
            phi = delay**s
        pdf = sf * phi
        cdf = cumulative_trapezoid(pdf, delay, initial=0)
        area = cdf[-1]
        mean = np.trapezoid(pdf * delay, delay) / area
        median = np.interp(area / 2, cdf, delay)
        # Galaxy integrated formed-mass and surviving-mass means, NOT SN host weighted.
        dg = np.geomspace(1e-7, t * (1 - 1e-9), n)
        sg = np.nan_to_num(csfh(z_of_age(t - dg), sfh))
        survive = 1 - 0.046 * np.log(dg / 0.000276 + 1)  # C14 A2: 0.276 Myr -> Gyr
        gm = np.trapezoid(sg * dg, dg) / np.trapezoid(sg, dg)
        gms = np.trapezoid(sg * survive * dg, dg) / np.trapezoid(sg * survive, dg)
        rows.append(
            dict(
                z=z,
                cosmic_age_gyr=t,
                mean_delay_gyr=mean,
                median_delay_gyr=median,
                formed_mass_mean_gyr=gm,
                surviving_mass_mean_gyr=gms,
            )
        )
    df = pd.DataFrame(rows)
    for col in [
        "mean_delay_gyr",
        "median_delay_gyr",
        "formed_mass_mean_gyr",
        "surviving_mass_mean_gyr",
    ]:
        df["delta_" + col] = df.loc[0, col] - df[col]
        df["correction_subtracted_mag_" + col] = 0.030 * df["delta_" + col]
    return df
