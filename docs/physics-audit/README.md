# Physical equations and implementation audit

Prepared 21 September 2026 from the current local workspace. This is the physical foundation and equation register for later tests of assumptions. It supplements the existing scientific reports; it does not recompute their cosmological posteriors.

**The implemented propagation and low-redshift cosmology equations are largely correct within their declared approximations. The project does not derive standardized Type Ia luminosity, intrinsic colour, dust populations or age corrections from fundamental physics. Those are empirical models whose assumptions remain testable. Four concrete implementation problems were corrected; an unphysical extrapolation of a dust law was independently reproduced and is explicitly unresolved as a population-model choice.**

The 123 numbered equation groups are written in LaTeX within Markdown, with stable labels, derivations, units, validity conditions, implementation pointers and primary references. Read in this order:

| Chapter | Physical content |
|---|---|
| [Luminosity and explosion physics](luminosity.md) | Nuclear binding energy, Ni–Co decay, deposition, energy conservation, expansion, diffusion, Arnett approximation, emergent spectrum and what SALT omits |
| [Extinction and attenuation](extinction.md) | Grain cross sections, radiative transfer, optical depth, scattering, magnitude extinction, CCM/O'Donnell/F99 equations, broad-band effects and dust geometry |
| [Observation chain](observation-chain.md) | Photon energy/redshift/time dilation, spectral and bolometric flux, band integration, detector charge/noise, calibration, K corrections, SALT coordinates and residual distance biases |
| [Cosmology and populations](cosmology-populations.md) | GR/FLRW, Friedmann and acceleration equations, distances, CPL, kinematic models, cosmic time, SFH/DTD/host ages, directional expansion, BAO and the CMB boundary |
| [Equation index](equation-index.json) | Machine-readable label-to-file/line index, generated from the chapter equations |
| [Finding register](findings.json) | Corrected errors, validated identities and assumptions that still require evidence |
| [Handoff for future agents](HANDOFF.md) | Completed work, verification, unresolved questions and workspace precautions |
| [Audit manifest](../../runs/physics_audit/manifest.json) | Code, chapter, result and baseline hashes plus runnable validation commands |

## What “established from physics” means here

Four different kinds of statement must remain separate:

1. **Physical law:** energy conservation, photon redshift, radioactive rate equations, radiative transfer, or GR field equations. Even these have a stated physical regime; a fixed laboratory decay rate, for example, assumes the appropriate ionization state.
2. **Derived approximation:** homologous expansion, diffusion, gray opacity, a passive dust screen, FLRW, a truncated distance series, or omission of radiation. We derive the resulting equation and identify the extra assumptions.
3. **Empirical model:** SALT surfaces/colour law, F99/CCM coefficients, a Tripp relation, a dust population, a star-formation history, DTD, or luminosity–age slope. Physics constrains these without uniquely determining their fitted coefficients or population distributions.
4. **Inference/calibration rule:** magnitudes, filter zero points, bias correction, sample selection, covariance and regression. These connect measurements to physical claims but are not additional fundamental laws.

Correct algebra under a chosen model does not verify that the model describes the selected supernova population. Conversely, the use of an empirical law does not itself invalidate a distance measurement. The issue is whether its calibration, support, errors and transport to other populations are tested.

## Corrections made

| Finding | Previous behavior | Correction | Effect on existing results |
|---|---|---|---|
| FIX-01: binned expansion rate | Calling `efunc(..., model='qbins')` fell through to a dark-energy formula and misinterpreted q-bin values as other parameters | Implemented the continuous, piecewise-integrated \(E(z)=\exp\int(1+q)d\ln(1+z)\) | Existing q-bin distance calculations bypassed that function; no claim of changed posterior |
| FIX-02: binned distance domain | The comoving integral saturated beyond the last configured bin, while q was extended or a plausible distance was still returned | NumPy functions reject redshifts outside 0–2.5; the JAX distance returns NaN outside its positive-redshift 0–1.3 domain and its runner rejects such data before sampling | Inspected saved Pantheon and hierarchy inputs are inside their respective domains |
| FIX-03: incomplete spectral band | The direct flux integrator silently dropped nonzero throughput outside the spectral model grid | Reject an evaluated band without complete wavelength support; unused invalid bands do not block valid bands | All 12 saved reference objects / 528 epochs give exactly unchanged flux predictions |
| FIX-04: transparent slab | The mixed dust-slab attenuation evaluated 0/0 at zero optical depth | Implemented zero attenuation at zero dust, a stable thin-limit expansion and rejection of invalid optical depth | Original positive optical-depth grid is bitwise unchanged |

These repairs follow directly from integral definitions and physical limits. Original source acquisitions and historical outputs were preserved. Changed implementations should have new manifests when used for future science runs.

## Important physical findings

**Luminosity is not predicted from an explosion here.** The new luminosity chapter supplies the energy-production/transport equations missing from the earlier mathematical overview. The working engine instead fits a trained spectral basis. It has no closed, validated mapping from progenitor age/metallicity/binary history to nickel production, ejecta structure, opacity and emergent light. A first-principles claim of universal post-standardization luminosity would exceed the implementation.

**Dust extinction can become unphysical under extrapolation.** A passive screen with nonnegative opacity must have \(A_\lambda\geq0\). Freshly compiled extinction routines and an independent equation implementation agree, yet sufficiently small extrapolated \(R_V\) produces negative red/near-IR extinction. At E(B−V)=0.1 and \(R_V=0.4\), the checked F99 curve gives \(A_{8000}=-0.01843\) mag, a flux multiplier about 1.0171. Re-reading 71,946 rows from the 25 pinned simulation HEAD files gives 4,601 such rows at 8000 Å under the historical prescription, and 4,791 under the exact F99 spline. These are **simulation truth rows, not observed supernovae, cosmology-selected fractions or measured distance biases**.

The scientifically justified immediate conclusion is that this extrapolation cannot be interpreted as the stated passive dust screen. Replacing it requires a nonnegative, validated extinction family and supported population distribution, followed by simulation, selection and bias recalibration. Clamping values to zero or arbitrarily truncating \(R_V\) would introduce a new physical/population assumption. This audit records the defect and the required discriminating test; it does not silently install such a replacement into vendor models.

**Colour coefficients and broad-band dust ratios are conditional.** SALT colour is not exactly E(B−V); fitted \(\beta\) is not generally \(R_V+1\). The identity \(R_B=R_V+1\) holds for extinction ratios defined from the *same* B/V attenuations and their difference. Identifying those quantities with trained colour coordinates or fixed monochromatic wavelengths is an additional approximation. A source comment equating the F99 input A_V with exact A(5495 Å) is also inaccurate: the verified ratio at \(R_V=3.1\) is about 0.97905. Source names do not establish physical normalization.

**Host age is not progenitor delay.** The SFH–DTD convolution is a rate-weighted population construction. An integrated host age is a different weighted functional. A local age slope does not derive the redshift evolution of standardized luminosity without a transport model, and a slope times a median age is not generally the mean luminosity shift.

**Acceleration remains conditional on luminosity evolution and geometry.** The relation \(q=-1+(1+z)H'/H\) and the tested flat distance integrals are correct. Arbitrary remaining luminosity evolution can still exactly compensate a different distance curve. Flat geometry, GR dynamics where used, population transport and early-universe ruler calibration are assumptions, not consequences of fitting a Hubble diagram.

## What was verified now

Fresh checks were run against the current code, rather than accepting the older audit reports as proof:

| Area | Independent check / result |
|---|---|
| Radioactive/thermal luminosity | Rate-equation ODE versus analytic decay-chain solution; deposited-energy normalization; diffusion energy balance and peak relation under the stated model; Planck/Stefan–Boltzmann integral |
| Extinction | Independent CCM/O'Donnell/F99 formulae versus freshly compiled pinned C routines; wavelength grid and low-R_V stress cases; broad-band derivative/curvature; host/MW frame contrast; slab integral and limits |
| Spectral propagation and photometry | Bolometric energy consistency to \(2.3\times10^{-16}\) relative; frequency/wavelength photon integrals to \(6.7\times10^{-16}\); K-correction sign; amplitude/inverse-square and gray-dimming checks |
| Cosmology | CPL distances against independent Astropy integrations; derivative and conservation identities; constant-q limits and JAX derivatives; bin-domain rejection; BAO combinations and cosmic-age calculations |
| Population clock | Cosmic SFH expression, normalized delay distribution and analytic constant-SFH / inverse-delay control; host-age and selection interpretation inspected separately |
| Regression boundaries | Saved in-domain fluxes unchanged; original positive slab grid unchanged; current saved redshift inputs fall inside bin support |

Exact metrics, tolerances, packages and source hashes are in the four JSON reports under `runs/physics_audit/`. Numerical roundoff agreement is evidence for the specified equation implementation; it is not evidence at that precision for astrophysical truth.

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 phase2/env-official/bin/python scripts/physics_audit/luminosity_check.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 phase2/env-official/bin/python scripts/physics_audit/extinction_check.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 phase2/env-official/bin/python scripts/physics_audit/observation_check.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/physics_audit/cosmology_check.py
```

The cosmology check needs the local hierarchy JAX environment; the extinction check needs a C compiler and the pinned source/data assets. This repository already required a local data bundle beyond Git. No new cosmology sampling, BBC production rerun or raw-image photometry reconstruction was performed for this equation audit.

## How to use the record for later assumption tests

Choose an equation label, then specify which physical quantity or approximation changes. Hold observation identities and units fixed. Propagate the change through the emitted spectrum, dust/propagation, instrumental response, selection and fitted coordinates as applicable. Compare observables and held-out predictions before translating the change into cosmology. A new physical dust or luminosity model must predict measured multiband/time-dependent fluxes, not only make a preferred Hubble residual disappear.

For each future test record the affected equation, parameter domain, observational discriminator, required recalibration and result. `findings.json` supplies an initial register. The chapters identify where a closed physical prediction is absent so that an assumed relation cannot later be mistaken for a derived law.

## Scope and remaining limits

The coverage is the physical equation chain underlying the local supernova, dust, population, directional and distance calculations, together with the BAO/CMB assumptions needed to interpret the imported probes. The audit directly inspects the first-party kernels and selected consumed vendor routines/assets. It does **not** certify all code in downloaded SNANA/Cobaya/Planck/ACT repositories, derive a complete nuclear reaction network or non-LTE atomic database, validate every instrument calibration, or reproduce the CMB-map likelihood. Those are distinct, much larger physical models; their absent implementation is stated explicitly rather than replaced by an invented derivation.
