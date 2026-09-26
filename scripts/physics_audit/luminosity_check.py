"""Independent physical checks; this is not an explosion or SALT fitter.

Run: .venv/bin/python scripts/physics_audit/luminosity_check.py
Equations and approximation boundaries: docs/physics-audit/luminosity.md.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import scipy
from scipy.constants import Boltzmann, Stefan_Boltzmann, c, h, physical_constants
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "runs/physics_audit/luminosity-check.json"
DAY = 86400.0
MSUN = 1.98847e33  # g; specified convention, not fitted data
MU = physical_constants["atomic mass constant"][0] * 1e3  # g
MEV = 1.602176634e-6  # erg
TN, TC = 8.80, 111.3  # mean lives in days; historical Nadyozhin convention
QNG, QCG, QCP = np.array([1.75, 3.61, 0.12]) * MEV
QUNIT = 1e43  # erg/s


def chain(t):
    t = np.asarray(t)
    ni = np.exp(-t / TN)
    co = TC / (TC - TN) * (np.exp(-t / TC) - ni)
    return np.array([ni, co, 1 - ni - co])


def radioactive_power(t, mass=0.6):
    ni, co, _ = chain(t)
    number = mass * MSUN / (56 * MU)
    return number * (QNG * ni / TN + (QCG + QCP) * co / TC) / DAY


def deposited_power(t, mass=0.6, gamma_time=35.0):
    ni, co, _ = chain(t)
    number = mass * MSUN / (56 * MU)
    fgamma = 1.0 if t == 0 else -np.expm1(-(gamma_time / t) ** 2)
    return number / DAY * (fgamma * (QNG * ni / TN + QCG * co / TC)
                          + QCP * co / TC)


def main():
    results = {
        "scope": "Independent equation/limit checks, not validation of actual ejecta or SALT training.",
        "constants": {"day_seconds": DAY, "solar_mass_g": MSUN, "atomic_mass_unit_g": MU,
                      "mean_lives_days": {"Ni56": TN, "Co56": TC},
                      "energies_MeV": {"Ni_gamma": 1.75, "Co_gamma_including_annihilation": 3.61,
                                       "Co_positron_kinetic": 0.12}},
        "versions": {"numpy": np.__version__, "scipy": scipy.__version__},
    }
    times = np.linspace(0, 1500, 2001)
    sol = solve_ivp(lambda t, n: [-n[0]/TN, n[0]/TN-n[1]/TC, n[1]/TC],
                    (0, times[-1]), [1, 0, 0], t_eval=times, rtol=2e-12, atol=2e-14)
    assert sol.success
    chain_error = float(np.max(np.abs(sol.y-chain(times))))
    assert chain_error < 2e-11
    assert np.min(chain(times)) >= -1e-14
    assert np.max(np.abs(chain(times).sum(axis=0)-1)) < 1e-14
    number = MSUN/(56*MU)
    integrated = quad(lambda t: radioactive_power(t, 1)*DAY, 0, np.inf,
                      epsabs=0, epsrel=2e-11)[0]
    expected = number*(QNG+QCG+QCP)
    energy_relative_error = float(abs(integrated/expected-1))
    assert energy_relative_error < 2e-10
    co_coefficient = number*(QCG+QCP)/((TC-TN)*DAY)
    ni_coefficient = number*QNG/(TN*DAY)-co_coefficient
    results["decay_chain"] = {
        "analytic_vs_independent_ODE_max_number_fraction_error": chain_error,
        "energy_integral_relative_error": energy_relative_error,
        "energy_per_initial_solar_mass_erg": integrated,
        "two_exponential_coefficients_erg_s_per_solar_mass": {"Ni_term": ni_coefficient,
                                                            "Co_term": co_coefficient},
        "Co_power_fraction_underestimate_if_buildup_factor_omitted": TN/TC,
        "Co_population_peak_days": TN*TC/(TC-TN)*np.log(TC/TN),
    }

    # Luminosity ODE and independent quadrature use different numerical methods.
    # Time is days and q,L are scaled by QUNIT in the following equations.
    tm = 13.0
    q = lambda t: deposited_power(t)/QUNIT
    arnett = solve_ivp(lambda t, ell: [2*t/tm**2*(q(t)-ell[0])],
                       (0, 120), [0], dense_output=True, rtol=2e-11, atol=2e-13)
    assert arnett.success
    ell = lambda t: float(arnett.sol(t)[0])
    peak = brentq(lambda t: q(t)-ell(t), 1, 80, xtol=1e-10)
    evaluation_times = np.array([0.01, 0.1, 1, 5, 10, 18, 30, 60, 100])
    integral = np.array([quad(lambda s: 2*s/tm**2*np.exp(-(t*t-s*s)/tm**2)*q(s),
                              0, t, epsabs=2e-13, epsrel=1e-11)[0]
                         for t in evaluation_times])
    ode_values = np.array([ell(t) for t in evaluation_times])
    quadrature_error = float(np.max(np.abs(ode_values-integral)/np.maximum(integral, 1e-8)))
    assert quadrature_error < 3e-7
    conservation = []
    for t in [5, 18, 40, 100]:
        # E/(QUNIT*DAY) = tm^2 L/(2t); both sides have day^2 after scaling.
        stored_weighted = tm**2*ell(t)/2
        heating_weighted = quad(lambda s: s*q(s), 0, t, epsabs=1e-10, epsrel=1e-11)[0]
        light_weighted = quad(lambda s: s*ell(s), 0, t, epsabs=1e-9, epsrel=1e-9)[0]
        conservation.append(abs(stored_weighted+light_weighted-heating_weighted)/heating_weighted)
    assert max(conservation) < 3e-9
    results["one_zone_diffusion"] = {
        "parameters": {"initial_Ni_solar_mass": .6, "diffusion_time_days": tm,
                       "gamma_escape_time_days": 35, "positron_kinetic_trapping": 1},
        "peak_days_after_explosion": peak,
        "peak_luminosity_erg_s": ell(peak)*QUNIT,
        "peak_L_over_deposited_Q": ell(peak)/q(peak),
        "peak_L_over_undeposited_decay_power": ell(peak)*QUNIT/radioactive_power(peak),
        "quadrature_vs_ODE_max_relative_error": quadrature_error,
        "time_weighted_energy_identity_max_fractional_residual": max(conservation),
    }

    # Derive Stefan-Boltzmann coefficient by independent Planck integration (SI).
    planck_integral = quad(lambda x: x**3/np.expm1(x), 1e-10, 700, epsabs=1e-12)[0]
    sigma_from_planck = 2*np.pi*Boltzmann**4/(h**3*c**2)*planck_integral
    sigma_error = abs(sigma_from_planck/Stefan_Boltzmann-1)
    assert sigma_error < 1e-8
    results["blackbody_limit"] = {"dimensionless_Planck_integral": planck_integral,
        "derived_sigma_W_m2_K4": sigma_from_planck, "relative_error_vs_scipy_constant": sigma_error}

    assets = ROOT / "phase2/official/portable_pilot/assets"
    a = np.loadtxt(assets/"salt3_template_0.dat.gz")
    b = np.loadtxt(assets/"salt3_template_1.dat.gz")
    assert np.array_equal(a[:, :2], b[:, :2])
    region = (a[:, 0] >= -15)&(a[:, 0] <= 45)&(a[:, 1] >= 3500)&(a[:, 1] <= 8000)
    # A linear function reaches its minimum on [-3,3] at an endpoint.
    lower = a[:, 2]-3*np.abs(b[:, 2])
    results["SALT_empirical_support"] = {
        "scope": "Native grid and x1 in [-3,3]; no actual epoch/passband impact inferred.",
        "nodes_in_phase_minus15_to45_wavelength_3500_to8000": int(region.sum()),
        "nodes_allowing_negative_M0_plus_x1_M1_in_that_rectangle": int(np.sum(region&(lower<0))),
        "inputs_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in [assets/"salt3_template_0.dat.gz", assets/"salt3_template_1.dat.gz"]},
    }
    results["script_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    results["passed"] = True
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(results, indent=2)+"\n")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
