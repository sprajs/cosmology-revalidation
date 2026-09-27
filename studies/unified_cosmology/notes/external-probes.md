# CMB and BAO information for the joint measurement

Two public CMB likelihood combinations are now executable with DESI DR2 BAO. Both evaluate predicted CMB spectra and lensing; neither replaces the CMB by a distance prior. This note establishes the inputs, numerical checks and limitations of those likelihoods. Evaluating them does not by itself establish a cosmological posterior.

## Measurements and shared physics

The initial reference combination uses Planck 2018 full Plik TTTEEE, Commander low-multipole TT, simall low-multipole EE and the public Planck lensing likelihood. A separate comparison replaces high-multipole Plik by its foreground-marginalized `Plik-lite` bandpowers. That replacement removes explicit foreground nuisance parameters by adopting the released marginalization; it does not make the CMB an uncalibrated distance measurement. The [Planck likelihood paper](https://www.aanda.org/articles/aa/full_html/2020/09/aa36386-19/aa36386-19.html) describes the measurement, foregrounds, calibrations and approximations.

The [official DESI DR2 release](https://data.desi.lbl.gov/public/papers/y3/bao-cosmo-params/README.html) supplies the 13-component BAO vector and covariance used here. It spans seven effective redshifts, 0.295–2.330, with the released combinations of transverse distance, radial distance and volume-averaged distance divided by the sound horizon. Within-bin covariance is retained. The released combined Gaussian matrix has zero covariance between different effective-redshift blocks; that is an input approximation, not a new measurement of their independence. These are released BAO summaries, rather than a new fit to raw galaxy clustering.

All likelihoods use one CAMB calculation for spectra, distances, expansion and the drag-epoch sound horizon. In particular, the joint CMB calculation does **not** give BAO an independently adjustable sound horizon. Supernova absolute magnitude can remain free in the joint analysis; the early-universe assumptions still calibrate BAO's absolute distance scale.

Declared theory choices are flat geometry, one massive neutrino with total mass 0.06 eV, effective neutrino number 3.044, standard CAMB recombination and BBN helium, and phenomenological PPF dark-energy perturbations. These choices are conditions of the measurement. PPF allows crossing of (w=-1); it is not a physical model of a particular dark-energy field. See the [CAMB dark-energy documentation](https://camb.readthedocs.io/en/latest/dark_energy.html).

## What the numerical checks establish

At the declared reference cosmology, which is **not** an optimized fit:

- Independent DESI distance construction and Cholesky whitening give χ² = 29.8689796119, compared with 29.8689796096 from Cobaya. The difference is (2.3\times10^{-9}), including scalar/vector background-rounding differences.
- The full Plik released check vector agrees within (4.1\times10^{-6}) in log likelihood. Native low-T, low-E and Plik-lite implementations agree with the independent `clipy` implementation at common spectra and three calibration values.
- The public native lensing text export differs from the original binary likelihood by about 0.00203 in reference log likelihood. Its observed bandpowers are rounded to a grid of (10^{-11}). After reconciling the different bin normalizations, the largest bandpower difference is (4.91\times10^{-12}). An exact quadratic decomposition attributes 0.00203521 to that rounding and −1.88×10⁻⁶ to rounded windows/covariance, inside the separately calculated 1.95×10⁻⁶ bound. Installed likelihood data are unchanged. This is a quantified representation difference, not exact binary identity.
- Repeated likelihood evaluations are deterministic. Cosmology changes recompute CAMB; changes to calibration/foreground parameters reuse spectra.

The compact records are [validation.json](../results/external_probes/validation.json), [acquisition.json](../results/external_probes/acquisition.json) and their executable sources. Reference log likelihoods include arbitrary normalization constants and should not be compared across different data combinations as if they were a model-selection statistic.

## The analytic CPL boundary is real

CAMB 1.6.6's analytic equation-of-state setter rejects (w_0+w_a>0). Thus the baseline CMB inference has an **effective implementation domain (w_0+w_a\leq0)**, even without an additional user-written prior. A rejected theory point is not evidence that observations excluded it. This corrects the initially incomplete “no extra prior” description in the design.

The documented tabulated-(w(a)) interface can evaluate some points beyond that boundary. The standalone [theory-domain audit](../results/external_probes/theory-domain.json) compares that interface with an explicitly diagnostic direct-field bypass, without changing the active adapter. At (w_0=-0.5,w_a=0.7), both produce finite spectra, (r_d\approx34.904\) Mpc and an early dark-energy fraction of 0.931 at (z=1100). A 4000-node table matches analytic background density evolution to (7.1\times10^{-9}); sampled spectra differ at about (6.5\times10^{-6}) relative. Numerical feasibility does not establish the validity of all perturbation/recombination assumptions in that extreme regime.

An additional [actual-likelihood scan](../results/external_probes/table-cpl-sensitivity.json) fixes the reference physical densities, Hubble constant, primordial parameters and nuisance parameters, then crosses the boundary at (w_0=-0.8,-1,-2). At (w_0+w_a=0.02), the respective changes in combined reference χ² are approximately 47,160, 17,857 and 20,350. The 2000-to-4000-node change in log likelihood is at most about 0.0064 for the checked components; the tabulated ΛCDM control reproduces the original likelihood to a few parts in (10^{-6}). One point just inside the boundary, (w_0=-1,w_a=0.98), emits 192 CAMB integrated-time warnings and is flagged explicitly. These are conditional evaluations, **not** a profile likelihood or a marginalized exclusion of the entire outside domain. The separate late-universe unrestricted-CPL analysis must retain that distinction.

## The contemporary Planck–ACT–SPT comparison

[The DES-Dovekie paper](https://lss.fnal.gov/archive/2025/pub/fermilab-pub-25-0842-csaid-ppd.pdf), arXiv:2511.07517v2, uses a different CMB combination from the initial Planck-only reference. The executable public-data comparison includes:

| Component | Input and treatment |
|---|---|
| Planck low multipoles | Commander TT and simall EE |
| Planck high multipoles | Plik-lite, complete bins below TT ℓ=1000 and TE/EE ℓ=600 |
| ACT DR6 primary CMB | Real NASA LAMBDA foreground-marginalized TTTEEE data; 135 bins |
| SPT 2025 primary CMB | Camphuis et al. D1 T&E lite data; 196 bins |
| ACT + Planck lensing | `actplanck_baseline` v1.2, 10 ACT and 9 Planck bins; released joint covariance |
| SPT 2023 lensing | Pan et al. likelihood with CMB-dependent response corrections; 12 bins |

Sources are the [ACT primary-CMB release](https://github.com/ACTCollaboration/DR6-ACT-lite), [ACT lensing release](https://github.com/ACTCollaboration/act_dr6_lenslike), [official SPT data repository](https://github.com/SouthPoleTelescope/spt_candl_data) and [candl data implementation](https://github.com/Lbalkenhol/candl_data). The newer SPT 2026 lensing data present in the installed library are **not** substituted for the paper's 2023 lensing data. Standalone Planck lensing is omitted when ACT+Planck joint lensing is present.

The joint ACT/Planck lensing covariance is materially non-diagonal: the largest cross-experiment correlation is 0.2063. The released inverse-covariance correction is retained, with Hartlap factor 0.9498746867. Independent Cholesky arithmetic reproduces the joint likelihood. Primary ACT, cropped Planck and both SPT likelihoods also pass independent selected-bin/covariance checks. SPT fixed-spectrum release tests agree to (3.7\times10^{-9}) in log likelihood and (8.6\times10^{-14}) in χ². ACT and joint-lensing reference tests pass the authors' declared floating-point tolerances; exact numbers and tolerances are preserved in [modern-release-checks.json](../results/external_probes/modern-release-checks.json).

There remain material convention and independence limits:

1. The paper literally says (A_{\rm ACT}=P_{\rm ACT}), while explaining that ACT is calibrated against Planck. The official combined SPT/ACT/Planck example instead fixes (A_{\rm ACT}=A_{\rm Planck}). Both mappings are exposed by the adapter. The default follows the official implementation; this is a declared independent configuration. The recovered author-chain header instead fixes A_ACT=1, as discussed below.
2. The paper does not state the ACT upper multipole cut. The ACT wrapper defaults to 6500 and a current combined example requests 8500. Both yield the same 135 retained bins and likelihood for the downloaded CMB-only SACC data.
3. The paper argues that residual Planck/ACT multipole overlap and ACT/SPT sky overlap do not materially affect its results. This analysis adopts the corresponding factorization; it does not recover a full cross-experiment sampling covariance. **SPT primary CMB and SPT lensing also share sky and underlying photons.** Their response corrections do not by themselves supply the missing primary/lensing sampling cross-covariance. The same distinction applies to shared primary-CMB/lensing information from Planck and ACT. The released joint ACT/Planck lensing covariance only addresses that joint lensing vector.
4. Calibration parameters are shared explicitly where required and their priors are applied once. Internal candl optical-depth/calibration priors are cleared before the declared external priors are applied, avoiding duplicated low-E/calibration information. Removing duplicate priors does not remove shared-observation covariance.
5. The cosmological priors in the baseline adapter differ from the paper's Bayesian-evidence priors. The latter include different parameterizations/ranges and an explicit (w_0+w_a<0) condition. This comparison is not an exact reproduction of the reported evidence or posterior.

The native candl Cobaya wrapper rounds returned log likelihoods to float32. At the reference, the errors relative to independent double-precision calculations are (3.4\times10^{-6}) for SPT T&E and −1.6×10⁻⁷ for SPT lensing. Increasing CAMB `AccuracyBoost`, `lAccuracyBoost` and `lSampleBoost` from 1 to 2 changes total reference log likelihood by +0.0874, mainly SPT T&E (+0.0813). This is a measured numerical sensitivity; parameter-dependent high-accuracy checks are still needed before precise posterior claims. See [modern-validation.json](../results/external_probes/modern-validation.json).

## What the released author chains actually specify

The [DES-Dovekie public release](https://github.com/des-science/DES-SN5YR/tree/c9a4fcafc4cbd19bd750dee47fc76194a45c181f/5_COSMOLOGY/chains) contains weighted Nautilus samples with embedded parameter, prior and module configurations. Their weighted CPL moments are w₀=−0.80288±0.05370 and wₐ=−0.72950±0.20131. These recover the published result from **the authors' samples**; they are not our independently inferred result. The chain has 52,224 stored rows and weight effective sample size 10,570. The corresponding ΛCDM chain has 51,200 rows and weight effective sample size 10,967.

The prior density reconstructed from the headers agrees with all 103,424 recorded prior values to 1.5×10⁻¹⁴, including normalization. Several differences are therefore verified configuration differences rather than interpretations of paper prose:

| Choice | Released chain header and public source | Our explicit modern comparison |
|---|---|---|
| Low multipoles | PlanckPy two Gaussian low-T bins; no active Commander or simall module | Commander TT and simall EE |
| Planck high-ℓ crop | Retain bins whose lower edge is ≤1000 TT or ≤600 TE/EE | Retain complete bins below these cuts |
| ACT calibration | A_ACT fixed to 1; P_ACT Gaussian σ=0.1, truncated to [0.9,1.1] | A_ACT=A_Planck; P_ACT Gaussian σ=0.003 |
| SPT T&E priors | Keep internal one-dimensional priors; historical YAML includes τ=0.051±0.006 and Tcal=1±0.0036, alongside the external Tcal prior | Clear internal priors; apply external calibration priors once and use low-E data |
| SPT lensing | `lens_only` variant | `lens_and_CMB` response variant |
| ACT/Planck lensing | Public header-interface default v1.1 and trim ℓ=2998 | v1.2 with declared modern settings |
| Neutrinos and helium | N_eff=3.046, three massive species totalling 0.06 eV, fixed YHe=0.245341 | N_eff=3.044, one massive species, BBN helium |
| Nonlinear spectrum | Mead2020; additional high-accuracy transfer settings | Mead2016 with separately declared accuracy settings |

The historical July 2025 SPT lite YAML is identical to the selected current YAML. The combination of the retained internal Tcal prior and the external prior would apply that Gaussian twice **if that public YAML was installed in the author's runtime**. This is distinct from the [official April 2026 SPT notice](https://lambda.gsfc.nasa.gov/product/spt/spt3g_d1_bandp_liklyhood_info.html) about repeated Planck priors in some cropped-Planck CMB-SPA chains; that notice alone does not establish that the Dovekie chains have the same bug.

Runtime identity remains unresolved. The chain records a working-directory Git revision, not a complete clean source tree or dependency lock. At that public revision the PlanckPy wrapper does not pass A_Planck to its calculator, while the released A_Planck posterior is informative. The named local Dovekie likelihood is also absent at that revision. We recovered the separately released Dovekie module and all public likelihood data needed for a conditional density reconstruction, but cannot silently identify these with unrecorded local files. [author-configuration.json](../results/external_probes/author-configuration.json) preserves the exact header options, prior check, sample summaries and proposal-only transformations; [author-acquisition.json](../results/external_probes/author-acquisition.json) pins the recovered source and historical lensing data.

The [five-row density check](../results/external_probes/author-density-check.json) resolves an important part of that gap numerically. We selected the maximum-posterior row and four weighted w₀ quantiles, then evaluated the recovered historical inputs. The literal public header/source differs from the recorded likelihood by a range of 17.015 in log likelihood. A diagnostic that applies the **same A_Planck calibration to both Planck and ACT** reduces the largest absolute discrepancy to 0.000109. There is no fitted zero point: the apparent common offset is exactly the independently calculated supernova magnitude-marginalization term, ½log(A/2π)=4.42633272, with A=1ᵀC⁻¹1=43938.3371. The small remaining differences are recorded, not rounded to exact equality.

This establishes a useful numerical reconstruction at five actual retained points, including the historical low-T treatment, retained SPT priors and shared calibration. It does not recover unrecorded source modifications or dependency versions, and cannot alone prove agreement throughout the posterior. The independent modern target remains unchanged. A configuration difference is not evidence of new cosmology; subsequent posterior comparisons must separate these input choices from changes in the observations.

## Reproduction and interfaces

Third-party source, data and generated spectra remain in ignored `.work/unified-cosmology/external-probes`. Public acquisition manifests hash the acquired likelihood data trees and selected likelihood source trees; core numerical package versions are pinned separately. Repository files contain first-party adapters, validators, compact results and pinned dependency lists.

The [executable reproduction guide](../code/external_probes/README.md) gives the full environment, acquisition and verification order, including the separate modern dependency lock and preserved acquisition manifest.


The modern lock deliberately selects the NumPy implementations of the 2025 primary-CMB and 2023 lensing datasets. The larger SPT package also declares JAX/GP dependencies for unused newer datasets; those are not installed. Use `--no-deps` with the complete recorded selected-runtime lock, rather than allowing those extra dependencies to change the numerical backend.

`adapter.external_info()` returns the initial full/lite CMB+BAO model; `modern_adapter.modern_info()` returns the modern comparison. Callers add one supernova likelihood using the same provider. `proposal.py` transforms released acoustic-angle proposal covariances to Hubble-constant coordinates. `sn_proposal.py` adds the local supernova distance Fisher information after projecting out a free absolute magnitude. Magnitudes are unused in that proposal construction. `modern_proposal.py` produces compatible 11/13-dimensional modern proposals. These operations change sampling efficiency, not scientific priors or data.

`spectral_training.py` completed a frozen local design of exact spectra for a possible computational approximation: 512 training and 96 independently seeded held-out attempts. Of these, 509 training points and all 96 held-out points have finite exact spectra and likelihoods. The three other training points have optical depths below the declared τ=0.01 prior boundary; they remain recorded and were not redrawn. There were no numerical theory errors. These are **not posterior samples**. The design uses acoustic-angle coordinates, stores the actual angle returned by CAMB, and records both the finalized requested multipoles and internal CAMB limits. Acquisition alone does not validate an emulator. Any accelerated inference needs held-out likelihood-error checks, exact-background and exact-spectrum closure, coverage/fallback diagnostics, and exact CAMB posterior reweighting with effective-sample and tail diagnostics.


The [optional lensing response compression](../results/external_probes/fast-lensing-validation.json) preserves the native joint covariance and all CMB/lensing response terms. The correction is linear in the input spectra about fixed fiducial spectra, so its large response matrices can be multiplied by the binning matrix once. Sixteen actual frozen spectra agree with the unmodified native likelihood within 2.9×10⁻¹⁴ in log likelihood and 6.3×10⁻¹⁶ in relative binned predictions; sixteen algebraic perturbations also pass. The one-thread median time for that likelihood fell from 0.496 s to 0.00166 s. This changes multiplication order, not observations or a theory approximation, and leaves the original adapter untouched. It does not accelerate CAMB itself.

An [independently authored response check](../results/external_probes/fast-lensing-independent.json) also compares 20 random response-matrix cases, including nonzero fiducial low multipoles and vanishing low-multipole normalization, against the original native implementation. Its largest absolute binned difference is 5.7×10⁻¹⁴. This separately checks orientation and low-multipole treatment; its artificial curves are not physical cosmologies.

A [lower-multipole CAMB check](../results/external_probes/theory-trim-check.json) was **not adopted**. At the reference and four frozen training points, reducing the request from ℓ=9001 to 6500 or 6325 changes total log likelihood by −0.0827 to +0.1691, without a consistent runtime gain. The exact target and all training spectra retain the original settings. The [aggregate verification](../results/external_probes/aggregate-validation.json) checks source/data identity and all 608 attempted spectra records while retaining these qualified outcomes.
