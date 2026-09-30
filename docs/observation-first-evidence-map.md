# Observation-first evidence for expansion, calibration and gravity

**Evidence-design snapshot: 27 September 2026.**

## Objective and completion boundary

Build an observational target for explanations that reconcile the local and early-universe expansion inferences and explain the evidence usually attributed to dark energy. Select no preferred mechanism. A successful explanation must predict the same observations with a consistent calibration, selection and uncertainty treatment.

This document is a source-linked observation register and proposed analysis protocol. Existing project records and the cited primary literature were inspected. New external data have not been downloaded, all cross-covariances have not been measured, and no new combined cosmological inference has been performed. Existing numerical fits are conditional benchmarks, not the starting observations or verdicts for a new theory.

“No separately postulated dark-energy component” and “no accelerating expansion” are distinct hypotheses. Modified dynamics could produce genuine acceleration without that component; a different spacetime interpretation could instead change what acceleration is inferred. Background distances alone do not uniquely separate a new gravitational law from an effective energy component. That interpretation ultimately needs a specified physical theory; the evidence target can be assembled first.

Likewise, resolving H0 does not require forcing two published central values to agree. It requires explaining the observations that produced them. In a sufficiently inhomogeneous spacetime, a local redshift–distance slope and a volume-averaged expansion rate need not be the same quantity.

## 1. Six observation groups

| Group | Observations to retain | Main question | Principal dependencies to expose |
|---|---|---|---|
| A. Absolute distances and local redshifts | Parallaxes; maser angles, Doppler velocities and accelerations; stellar fluxes and periods; GW strain and host redshifts; lens images, time delays and stellar spectra | Where does the absolute distance-scale disagreement enter? | Shared anchors, instrument calibration, local velocity field, source physics, lens mass models, GW generation/propagation and source selection |
| B. Expansion geometry and elapsed time | Galaxy angular positions and redshifts; radial/transverse clustering; SN flux histories; differential galaxy spectra; repeated spectral-line positions | What distance and time relations are supported across redshift, with no chosen dark-energy equation of state? | Ruler calibration, fiducial clustering reconstruction, stellar ageing, SN luminosity evolution, curvature, photon propagation; redshift drift is a future precision discriminator |
| C. Early-universe observations | CMB temperature/polarization maps or likelihood-supported spectra; primordial deuterium absorption; helium spectroscopy as a subsequent separately sourced addition | Which early physical conditions must an explanation preserve or change? | Beams, foregrounds, common sky modes, atomic/nuclear rates, recombination, assumed matter content and gravitational law |
| D. Our location and nearby structure | Galaxy counts and redshifts across the sky; redshift-independent distance indicators; line-of-sight lensing and cluster microwave signals | Does the local expansion pattern track an independently constrained density/velocity field? | Catalogue completeness, galaxy-to-matter mapping, distance calibration, overlapping distance objects, gas physics and observer location |
| E. Gravity acting on light and matter | Galaxy shapes, clustering anisotropy, stellar/galaxy velocities, CMB lensing and correlations with structure | Do matter motion and light deflection require the same gravitational response at different scales and epochs? | Shared sky/sample variance, galaxy bias, intrinsic alignments, redshift calibration, nonlinear/baryonic modelling and assumptions mapping observables to gravitational potentials |
| F. Source physics and independent gravity checks | Host spectra and optical–IR photometry; resolved stellar populations; binary-pulsar timing; lunar ranging; GW and electromagnetic timing | Could a proposed effect alter the measuring instruments supplied by nature, and does it satisfy gravity tests elsewhere? | Stellar evolution/atmospheres, dust, gravity-dependent source luminosities and clocks, orbital modelling, emission delays and environmental screening |

These are different measurement mechanisms, not six automatically independent likelihood factors. Individual observations can serve more than one question and must then be counted once.

## 2. Concrete observation register

Status labels: **project** = a local implementation or audit is documented; **candidate** = primary source inspected, acquisition and likelihood admission remain to be done; **future** = a prospective discriminator rather than an established precision detection. A project label does not certify every upstream assumption or current numerical fit.

| ID | Starting source and status | Retain or acquire | Inference boundary |
|---|---|---|---|
| A01 | [Released SH0ES reconstruction](../studies/unified_cosmology/notes/distance-ladder.md), project | The selected photometric/anchor system and covariance; separately the 37-host Cepheid-only factor after removal of all SN rows | Do not append the full ladder H0 posterior to its constituent SNe. Cepheid distances remain conditional on stellar/calibration physics. |
| A02 | [Megamaser Cosmology Project XIII](https://arxiv.org/abs/2001.09213), candidate | Galaxy-specific disk observables or distance likelihoods, redshifts and velocity treatment | SN-independent does not mean gravity-independent. NGC 4258 is also a ladder anchor, so its reuse must be explicit. |
| A03 | [GWTC-4 public products](https://gwosc.org/GWTC-4.0/) and [GWTC-5 cosmology manuscript](https://arxiv.org/abs/2605.27227), candidate | Event likelihood information or posterior samples with removable sampling priors, detector calibration, detection injections and redshift/host information | Choose one explicitly versioned cumulative catalogue. Source-mass features are statistical redshift information, not measured host redshifts. If gravity changes, both wave generation and propagation require admission checks. No GW-only verdict on acceleration is asserted here. |
| A04 | [TDCOSMO 2025](https://arxiv.org/abs/2506.03023), candidate | Lens images, measured time delays, kinematics and environment likelihoods with flexible mass distributions | Its headline H0 result incorporates a Pantheon+ matter-density constraint. It is not an untouched SN-free likelihood. Preserve mass-sheet and stellar-orbit uncertainties. |
| B01 | [DESI DR2 BAO measurement paper](https://arxiv.org/abs/2503.14738) and [local input audit](../studies/unified_cosmology/notes/external-probes.md), project | Released 13-component vector/covariance as an initial compressed branch; lower-level clustering for extensions outside its validated domain | Retain a free sound-horizon ruler initially. BAO ratios do not separately identify H0 and the ruler. Zero inter-redshift covariance blocks are a release approximation. |
| B02 | [Supernova measurement programme](experimental-plan.md), project | A single identified compilation at a time, signed fluxes where available, epochs, calibration/training lineage, selection and host information | Distance tables contain previous fits/corrections. Arbitrary achromatic luminosity evolution is degenerate with the inferred distance shape. Optical–IR information alone cannot eliminate that degeneracy. |
| B03 | [Chronometer observations](https://cluster.difa.unibo.it/astro/CC_data/) and [systematic-covariance method](https://arxiv.org/abs/2003.07362), candidate | Differential spectral-age information and the stellar-population covariance recipe; record data/model versions together | The source page says its tabulation does not itself provide the full covariance. Stellar-library uncertainties can link different redshifts and connect to host-age work. Avoid cosmological age ceilings when independently testing cosmic time. |
| B04 | [ESO redshift-drift programme](https://www.eso.org/sci/facilities/eelt/science/drm/C2/), future | Repeated line positions, observing times, instrumental drift and source-acceleration control | Forecast sensitivity is not a detected cosmic signal. Even this comparatively direct observable requires an observer/source and redshift–expansion mapping. |
| C01 | [CMB likelihood inputs](../studies/unified_cosmology/notes/external-probes.md) and [dependence audit](../studies/unified_cosmology/notes/probe-dependence.md), project | Temperature, polarization and lensing likelihood components, window functions, masks, calibrations and shared covariance | Use neither a Planck H0 prior nor a standard-model sound horizon as a theory-free measurement. Existing CAMB targets impose dynamics and have separately documented numerical qualifications. |
| C02 | [Primordial deuterium spectroscopy](https://arxiv.org/abs/1710.11129), candidate | Absorption-profile likelihoods and inferred abundance with astrophysical systematics, followed by explicit nuclear-rate inference | An abundance is closer to observation than a BBN-derived baryon-density prior. Changing expansion/gravity or nuclear physics requires recomputing that translation. |
| D01 | [CosmicFlows-4 direct-distance void test](https://arxiv.org/abs/2506.10518), candidate | Underlying galaxy measurements, selection and distance-indicator calibration, rather than the fitted void profile | A velocity reconstructed by assuming a background H0 cannot independently confirm that same background. Audit any SN-derived distances before claiming a SN-free branch. |
| E01 | [DESI full-shape analysis](https://arxiv.org/abs/2411.12022), [release description](https://data.desi.lbl.gov/doc/releases/dr1/vac/full-shape-cosmo-params/) and [DES galaxy/shear measurements](https://www.darkenergysurvey.org/des-year-3-cosmology-results-papers/), candidate | Clustering and shear observables with estimator windows and joint covariance; retain scale cuts and nuisance response | Published growth/gravity parameters are model-conditioned summaries. BAO and full shape from the same galaxies are not separate independent observations. |
| F01 | [Calibrated host observations and physical-age bounds](../studies/unified_cosmology/notes/calibrated-host-physics.md), project | Signed spectral-band vectors, their covariance and response, local/global aperture information and full spectra when available | 54/55 hosts admit very broad formed-mass ages under the tested flexible mixtures. These are compatibility sets, not precise ages or progenitor delays. |
| F02 | [Lunar-ranging gravity tests](https://research.uni-hannover.de/en/publications/relativistic-tests-with-lunar-laser-ranging/), candidate | Ranging/timing information or likelihood-supported gravitational-parameter constraints with orbital nuisances | Local bounds do not automatically apply to cosmological scales in an environment-dependent theory. Specify the mapping rather than discarding the bounds. |
| F03 | [GW170817 and GRB 170817A](https://arxiv.org/abs/1710.05834), candidate | Arrival times, sky localization and emission-delay uncertainty | Speed and propagation tests constrain some changes to gravity, not every possible gravity law. The same event used in A03 is counted once in a joint likelihood. |

The cited releases are starting points, not a claim of an exhaustive or fully acquired data inventory. Exact file hashes, observation IDs, masks and available covariance products are acquisition requirements before numerical use. Do not substitute posterior parameter tables for missing likelihood information.

## 3. Keep assumptions in separate layers

1. **Recorded signal:** detector counts/strain, timestamps, angles, wavelength calibration and observing conditions. These already need instrument models.
2. **Measurement response:** fluxes, spectral features, redshifts, angular correlations, time delays and calibrated distance likelihoods. Publish nuisance parameters and selection here.
3. **Kinematic interpretation:** distances and expansion versus redshift with explicit geometry and propagation assumptions, without choosing a dark-energy equation of state. A flexible fit still has resolution, smoothness and boundary assumptions.
4. **Dynamical explanation:** translate those measurements into matter content, gravity, spacetime evolution or a dark-energy component. Keep this stage replaceable.

Carry a homogeneous/isotropic kinematic branch and a separate direction/environment-dependent branch. An arbitrary H(z) reconstruction is not independent of homogeneity: a single global H(z) may be inadequate in the second branch. In that case retain observed redshift–distance relations by direction and distance range before introducing an averaged H0 or q.

Within an FLRW kinematic branch, q(z) = (1+z) H'(z)/H(z) − 1 does not require a dark-energy equation of state, but it does require that geometry and a stable derivative reconstruction. Report interval-average acceleration separately from instantaneous q(0). A free overall H0 normalization does not determine that derivative. A purely uniform distance rescaling therefore cannot by itself remove the relative-distance acceleration signal.

For nonstandard gravity, reassess source and propagation response before accepting any compressed likelihood. In particular, changing gravity can affect stellar pulsation, stellar lifetimes, supernova luminosity, lens dynamics, the early sound horizon and GW amplitude. A common gravity parameter is shared physics to propagate, not an extra independent correction to insert into each dataset.

## 4. Dependency ledger: retain covariance, remove duplicate information

Separate three things: shared observations; shared calibration/astrophysical uncertainty; and common physical predictions. Correlated cosmological signal is useful information. It is not contamination to erase.

| Link | Present evidence | Required handling |
|---|---|---|
| Calibrator SNe ↔ Hubble-flow SNe | Local ladder audit reports nonzero covariance and a maximum absolute correlation of 0.2783 | Retain the released joint covariance; do not use its H0 output as another independent prior. |
| Dovekie ↔ Pantheon+ | Local matching identifies 273 common events in 301 measurement rows under the stated rules | Alternative compilation branches, or a deduplicated observation model; never naive multiplication. |
| Cepheid/TRGB/maser calibration routes | Potential reuse of geometric anchors and photometric systems | One shared anchor likelihood; distinguish different stellar methods from independent absolute zero points. |
| ACT lensing ↔ Planck lensing | Current project joint covariance has maximum absolute cross correlation 0.206284 | Keep the joint block. Do not replace it with separate factors. |
| CMB primary ↔ lensing; different CMB experiments | Several cross blocks omitted under published approximations | Match maps, masks, estimators and noise before transferring validation. Track unquantified terms explicitly. |
| DESI BAO ↔ full shape; galaxy clustering ↔ shear | Reused galaxies and overlapping sky/large-scale modes | Use matched joint covariance or a documented nonoverlap/scale-selection branch. Disjoint sky still requires checking global shared systematics. |
| Dovekie ↔ DESI BAO; current CMB lensing ↔ DESI BAO | Exact-combination cross-covariance bounds are not established in the local audit | No automatic zero-covariance certification. Estimate with matched simulations/response or publish conditional sensitivity only. |
| Chronometers ↔ host ages ↔ lens kinematics | Potential common stellar libraries, spectral calibration and population assumptions | Explicit shared model alternatives or nuisance response; strength must be measured, not guessed. |
| Nearby SNe ↔ masers ↔ sirens ↔ local distance catalogues | Common local velocity field, some shared objects/anchors | Use source-level identity and a common velocity treatment with uncertain reconstruction. |
| GW distance ↔ gravity tests on the same event | Shared strain, calibration, source and timing information | Joint event-level analysis or explicit conditional factors; avoid independent multiplication of published summaries. |

The numerical examples above are read from existing [calibration](../studies/unified_cosmology/notes/distance-ladder.md) and [probe-dependence](../studies/unified_cosmology/notes/probe-dependence.md) records, not recomputed for this document.

For small, validated linear responses, a shared nuisance u with covariance S can contribute J S Jᵀ, including cross-probe terms. Use this only if that uncertainty is not already in the supplied covariance. Nonlinear calibration and selection may instead require explicit latent variables and forward simulations. Unknown covariance remains unknown: an arbitrary positive-definite matrix is a sensitivity scenario, not an empirical estimate.

Whitening can decorrelate Gaussian residual coordinates mathematically; it neither creates independent experiments nor removes shared physics. Removing one probe is useful, but is not itself a measurement of the omitted cross-covariance. Hold out entire sources or calibration/sky blocks where appropriate, rather than randomly splitting repeated observations of the same object.

## 5. The combined observations that would discriminate explanations

These are proposed tests and decision rules, not completed results.

| Test | Combine | Distinguishing outcome and limit |
|---|---|---|
| T1. Locate the absolute-scale discrepancy | A01–A04, with source/anchor overlap and a common local velocity treatment | If independently calibrated non-SN distances disagree with the stellar ladder, investigate calibration/source physics; if they agree, a SN-only explanation loses scope. Limited precision can leave both possibilities open. |
| T2. Separate geometry from brightness evolution | B01+B02+B03, adding A03/A04 where supported; retain a free BAO ruler | A brightness-specific drift with consistent geometric/clock evidence motivates luminosity or opacity work. Consistent drift across mechanisms points toward geometry or shared physics. Neither pattern alone identifies a unique cause. |
| T3. Test a local-environment explanation | D01 plus independently calibrated A measurements and sky-dependent B measurements | Predict the sign, size, direction and distance dependence of redshift/distance residuals using an independently constrained density field. Fit held-out regions/objects. Do not infer the void solely from the same residuals then call agreement independent. |
| T4. Test gravitational response | E01 plus A04 and F02/F03 | Compare light deflection, matter motion and wave propagation with one parameterized response. An expansion fit alone does not establish an altered gravitational law. |
| T5. Locate early-versus-late disagreement | Late-time A+B first; then C01+C02 | Infer which ruler/early-physics combinations are required by late observations, then test their early-universe predictions. Recompute recombination/BBN when altered physics requires it. Never insert a standard-model Planck H0 prior into this test. |
| T6. Test whether the measuring sources change | F01 plus A01/B02/B03, with independent environments and held-out flux/spectral predictions | Quantify residual source evolution after existing corrections. If gravity is varied, propagate it through clocks and luminosities too. Host-age trends alone do not measure a redshift-dependent distance correction. |

Two useful consistency relations are available without choosing a dark-energy equation of state, but each has conditions. Distance duality compares luminosity and angular distances assuming metric light propagation and photon conservation; a discrepancy can mix opacity, source evolution and calibration. A [recent SN–DESI study](https://arxiv.org/abs/2604.02433) demonstrates this observational comparison, but is not independent evidence beyond its reused inputs. The [Clarkson–Bassett–Lu relation](https://arxiv.org/abs/0712.3457) tests consistency of distances and expansion with FLRW geometry; derivative noise and observational calibration must be propagated. Neither is a universal test outside its assumptions.

## 6. Order of work and gates

1. **P0 — Assemble and identify.** Freeze source versions, observable units, redshift/sky ranges, raw-versus-derived status, source IDs, selection, covariance and calibration lineage. Record missing likelihoods and cross blocks. Gate: no observation or anchor counted twice; every claimed independence has a reason and scope.
2. **P1 — Absolute distances and late geometry.** Execute T1 and T2 without using early-universe ruler calibration. Reconstruct only redshift resolution supported by the measurements. Keep luminosity evolution, photon propagation, curvature and directional alternatives explicit. Gate: calibrated synthetic recovery, adequate precision and physical support; otherwise publish a bound or non-identification.
3. **P1 — Source and environmental explanations.** Execute T3 and T6 with held-out galaxies/regions and source properties. This connects directly to the existing host-age and local-calibration work. Gate: the proposed effect predicts more than the original distance residuals.
4. **P2 — Gravity and early-universe consistency.** Add T4 and T5 with supported joint covariance and forward physics. Gate: one explanation must account for multiple measurement mechanisms without inconsistent recalibration between them.
5. **P3 — Theory comparison.** Once a physical proposal is specified, predict all admitted observables, including the previously held-out ones, and apply a common likelihood with a declared complexity/prior treatment. Gate: distinguish a good fit, compatibility, predictive success and unique identification. Retain failures and indistinguishable alternatives.

The existing [probe-omission helper](../studies/unified_cosmology/notes/probe-omission.md) is useful only within its qualified parent model; the inspected note reports synthetic validation but no observational omission result. Reweighting a flat-GR parent cannot discover a theory excluded from that parent's physical support.

Success could be: an adequate standard model after better calibration; acceleration caused by different dynamics without a separate dark-energy component; a changed inhomogeneous interpretation of apparent acceleration; or an unresolved family of observationally equivalent explanations. None is selected in advance.

## 7. Present assessment

An observation-first comparison is feasible; complete theory independence is not. The project already has substantial components for stellar/SN measurements, calibrated distances, BAO and CMB, plus useful dependence audits. It does not yet have an admitted joint measurement spanning all six groups. The strongest next contribution is to expose exactly which combinations are identified by the measurements and which depend on calibration, shared uncertainty or a physical assumption.

The general question is therefore: **Can one physically coherent explanation reproduce the absolute distance measurements, relative distance/time history, local environment, gravity response and early-universe observations, while explaining why standard reductions produced both discrepancies?** A flexible empirical fit can identify the required changes; it cannot by itself establish that an underlying no-dark-energy theory exists or is preferred.
