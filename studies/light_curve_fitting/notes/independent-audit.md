# Independent falsification audit

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../README.md) for current execution status.

Independent Astra Ultra audit, 20 September 2026. This investigator received the other branches' initial written reports after they were recorded, then read the relevant source passages and code. Agreement between investigators is not counted as independent scientific evidence. The tests below change the mathematical implementation or interrogate a specific inference; they reuse the same observations.

The audit attempted to falsify (1) the numerical cosmology shift, (2) the evidence that fully corrected matched ages have a weak additional residual association, (3) the critique of the C26 mock and host-to-progenitor conversion, and (4) the directional Taylor-truncation result. It found **no consequential numerical implementation failure in the tested central calculations**. It did find a source-level sign discrepancy that needs to remain visible, and several boundaries beyond which the calculations do not support the claimed inference.

## A source discrepancy requiring explicit qualification

Hoyt et al. Union3.1 Appendix F states a mean distance-modulus shift of **−0.073 mag** for 114 low-redshift rows, −0.012 mag averaged over the low-redshift bin, and describes adding the shifts to `m_b_corr`. The later Roy Choudhury Eq. 2 and its public reconstructed vector use `new − old = f_low(c) − f_high(c)`, with mean **+0.07220 mag** on the affected rows. The standardization branch's machine-readable result preserves this discrepancy; the audit requested that it also be prominent in the final narrative. The cosmology report now explicitly records it.

The investigation's code follows the positive public-vector difference correctly; its BAO+SN sensitivity is consequently traceable to that vector. This is **not an independently verified reproduction of Hoyt's executed private correction or cosmological run**. A sign convention or a prose error is possible, but the absent executed Hoyt vector prevents establishing which. The mismatch is not evidence of a sign bug in the local code. The principal investigator was notified before the final synthesis, and the final record must preserve both signs and this boundary. Primary passages: `sources/updates/2026-09-20-standardization/2601.19424.txt`, Appendix F; `2607.24443.txt`, Eq. 2 and its following paragraph. [Hoyt et al.](https://arxiv.org/abs/2601.19424), [Roy Choudhury](https://arxiv.org/abs/2607.24443).

## Independent distance likelihood and correction-sign test

The audit implementation imports none of the investigators' cosmology functions. It reads the original magnitude table and covariance, subsets their original row indices, whitens with Cholesky, and explicitly fits the common magnitude offset. It solves the coupled differential equations

`d(rho_DE/rho_DE0)/dz = 31+w(z) (`docs/investigations/rho_DE/rho_DE0`; historical local reference)/(1+z)` and `dchi/dz = 1/E(z)`

with DOP853 instead of using the main implementation's explicit CPL density and Gauss–Legendre distance quadrature. The BAO predictions are independently assembled and whitened with the released covariance. It then reoptimizes all four cosmological/scale parameters after moving away from the reported mode. The correction is explicitly subtracted from the observed magnitudes.

| Data/correction, BAO+Pantheon+ | Direct-minus-main chi-square at main optimum | Independent optimum q0 | Difference from main optimum q0 |
|---|---:|---:|---:|
| Original distances | +1.08e−7 | −0.42816820 | −9.41e−8 |
| C14 fixed median template | +5.28e−7 | +0.06202534 | −2.06e−8 |
| C14 fixed mean template | +4.37e−7 | −0.03624928 | +1.38e−7 |

These are maximum-likelihood values, **not posterior means**. All optimizations converge. Matter-only and de Sitter distance limits agree within 9.5e−13 in dimensionless comoving distance. The independent calculation therefore does not support a correction-sign, covariance-order, intercept, distance-integral or interpolation explanation for the large posterior shift. It confirms the *conditional calculation*, not the astrophysical validity of the template or exact reproduction of Son's unavailable run bundle.

The contrast between the mean and median corrections is a serious interpretive sensitivity: it moves the optimum across zero using identical distance observations and cosmological priors. A linear luminosity law predicts a consistently weighted **mean** luminosity from the corresponding mean delay. A cosmic median-delay curve has no general right to replace that moment. The selection weights and response distribution could motivate another estimand, but they would have to be specified and tested. This is a model-identification issue, not numerical quadrature uncertainty.

## Posterior and prior checks

The auditor recomputed q0 from every retained CPL sample rather than substituting marginal parameter means. Results match the main summaries. Reweighting to a flat-q0 measure rather than flat-w0, on the **same original support**, changes the fixed-median result from q0=+.06320, P(q0<0)=.19816 to +.06219, .20207. A log-H0rd measure gives +.06357, .19676. For the shared-slope case, the corresponding probabilities are .33029, .33670 and .32700. Importance weight ESS fractions exceed .999 for these corrected cases. These two benign prior-measure changes do not remove the conditional shift.

This does not test arbitrary priors or enlarge the original support. BAO-only has 2.31% of draws in the lowest .05 interval of its wa range and its q0 is poorly localized; low-redshift BAO alone is not a model-independent observation of present deceleration. The explicit `w0+wa<0` cut is also part of the answer. Near-unit sampled probabilities in uncorrected cases must be reported as finite Monte Carlo tail resolution. The initial `1.0` probability and zero block-MC-error fields are empirical draw fractions, not literal posterior certainty. The principal investigator corrected the narrative/plot presentation after this warning.

Even perfect convergence cannot turn a fixed cosmological clock, fixed SFH/DTD, disputed transported slope and chosen mean/median into a uniquely measured correction. The broad-amplitude sensitivity and alternative kinematic fits in the main branch appropriately delimit the claim; this audit did not independently run every posterior ensemble or a full CMB likelihood.

The subsequent main-branch Gaussian posterior-predictive check is also consequential. Its observed quadratic discrepancy is about 1404–1406, versus 1589 expected for simulated residuals after profiling one intercept. The analytic lower-tail probabilities average .000352, .000355 and .000395 for the uncorrected SN, fixed-median SN and shared-slope joint runs. Inspection confirms the `chi-square(N−1)` replicate law is appropriate for this **conditional Gaussian fixed-covariance generator**: it projects out one intercept and does not refit cosmology for each replicate. This is an underdispersion warning about that generative interpretation, not evidence selecting an age correction. A total release covariance encodes systematic modelling choices as well as measurement variance, and a selected corrected-distance sample need not follow that simple Gaussian generator. No ad hoc error rescaling is justified by this test alone. The auditor checked the formula and source/result chain; it did not independently repeat the 4,000 posterior selections or the 200 recovery simulations.

## Independent age/mapping checks

The auditor independently parsed both original age tables and joined numeric SN IDs only within the SDSS survey. It recovers 199 and 102 input rows, **33 shared explosions**, and the 196-object union with cumulative zHD cut counts **100,144,176,192,196**. The release's tabulated-error versus covariance-diagonal RMS discrepancy is independently reproduced as **0.0603696 mag**. Thus the sample and uncertainty issues are real properties of the acquired inputs.

Direct full-covariance fixed-age GLS gives corrected slopes −.004950 and −.005082 mag/Gyr for the two overlap policies; reversing the complete bias correction gives −.012940 and −.012831. These differ from the latent-age slopes, as expected: treating noisy ages as exact attenuates the slope. They are an independent algebra/covariance check, not a competing physical estimate. The original LINMIX runs and profile intervals were inspected, not independently rerun in this bounded audit.

With full STAT+SYS weighting, the observed age variation projected onto one LCDM matter-density distance derivative is 1.23%/0.73%, or 1.55%/1.12% onto three CPL derivatives. The age branch's 4.1%/2.6% uses the separate **tabulated-error WLS** weighting. The difference is an explained weighting choice, not a contradiction. These calculations concern noisy measured ages. They do not bound how much an unobserved redshift-dependent true-age component can be absorbed by cosmology or standardization.

For the population evolution calculation, the auditor switched the integration variable from delay to **formation redshift**, used the analytic flat-LCDM age formula, adaptive quadrature and root finding for the median, and independently checked B13 Eq. F1/Table 6. Under the disclosed LCDM70/B13 convention:

| DTD | Mean-delay change to z=1 | Median-delay change | .030 times median change |
|---|---:|---:|---:|
| Smooth C14 | 4.47354156 Gyr | 5.29036810 Gyr | .158711043 mag |
| Hard 40-Myr cutoff, slope −1.13 | 3.33336988 Gyr | 2.42079534 Gyr | .072623860 mag |

Both mean and median computations agree with the mapping branch within 9e−7 Gyr. The published W26 cosmic-SFH-only ∼1.5-Gyr comparison therefore remains unreproduced by this explicit B13 empirical implementation, and cannot be explained by the mapping branch's numerical integration. The unidentified B13 representation/clock/implementation remains an exact-reproduction barrier. This is not a demonstrated arithmetic error by W26, whose executed inputs are unavailable.

The mapping derivations survive the attempted counterexamples: a stable affine conditional mean permits cancellation of slope compression and age evolution even with broad latent delay scatter; a nonlinear or changing conditional map does not generally permit replacement by one slope times one mean or median. The mapping report correctly avoids claiming that broad scatter alone destroys the cancellation. Conversely, algebraic cancellation for a relabelled proxy cannot validate a causal luminosity law or repair selection/confounding.

## Falsification of the claimed necessity in the C26 mock

The published C26 text explicitly converts assigned age-dependent magnitudes to residuals with the same mean in every redshift bin, and says this is required when LCDM defines the baseline. It is not required by `HR=mu_observed−mu_reference`: with the true fixed baseline and injected `Y=bA`, HR retains `bA` plus the intercept.

An independent noiseless realization with b=−.034 gives:

| Operation | Fitted slope |
|---|---:|
| Fixed correct distance baseline | −.03400000 |
| Centre residuals within redshift bins only | −.02069613 |
| Within/total age-variance identity | −.02069613 |
| Centre both ages and residuals | −.03400000 |
| Subtract a fitted linear-redshift luminosity baseline | −.02082245 |

The last row is the strongest counterevidence to an overly broad critique of the mock: **fitting a baseline can indeed hide a real redshift-correlated luminosity component**, and in this constructed population a smooth fitted baseline nearly reproduces bin-centering dilution. The correct conclusion is therefore narrow. C26's mock does not demonstrate that a fixed LCDM residual definition inevitably causes its dilution; it does not establish the size of that effect in the real jointly fitted and selected sample. The example also does not prove that dilution is absent in real cosmological fitting. That requires an end-to-end injection and recovery with the actual fitting, noise, selection and host likelihood.

Source: [published Chung reply](https://doi.org/10.1093/mnras/stag1513), Section 2 and Figure 2 methods. No exact replication of its private age-PDF arrays or noise configuration is claimed.

## Directional branch cross-check

The auditor independently rebuilt the C1 heliocentric reversed-bias likelihood from the raw row/covariance inputs using direct Cholesky, including its fitted additional scatter and Gaussian determinant. It matches the directional investigator's eigenbasis computation for isotropic and dipolar models to **1.1e−11 in −2 log L**. This verifies a delicate normalization step independently.

For 1,552 actual redshift/covariance rows below zHEL=.8, noiseless exact flat-LCDM data with true q0=−.55 fit by the cubic distance approximation give **q0=−.523632565**, a **+.026367435** bias, agreeing with the directional branch within 9.3e−9. The magnitude truncation error at z=.8 is **−.029467369 mag**. It is measurable and must be included in interpretation, but is much smaller than the ∼+.35 coefficient shift from the chosen age correction in these fits. It does not explain that shift by itself.

A dipole improvement is not directly an acceleration-sign test. The dipole decay scale is unidentified at zero amplitude, and some fitted scales reach their lower bound. Uncalibrated Wilks conversion is therefore unsafe; coefficient substitution into a cubic isotropic formula is also a phenomenological angular fit rather than a unique covariant deceleration-field reconstruction. These limits were communicated to the directional investigator. C2 covariance acquisition/order was resolved during the initial audit; the bounded follow-up below checks the resulting likelihood at fitted points rather than performing another optimization.

## Directional follow-up: C2 and the revised Ray paper

The independent C2 check uses a different observable basis: `u=(mB+alpha*x1−beta*c,x1,c)`. It applies this unit-determinant transformation to the **entire released measurement covariance**, then adds diagonal intrinsic population variances `(sigmaM²,sigmaX²,sigmaC²)` for each object. The directional implementation instead leaves observations in `(mB,x1,c)` and adds the induced off-diagonal population covariance. These formulations must be equivalent, including the Gaussian determinant. All three population means are independently profiled by whitened least squares. The raw table, explicit index array and interleaved covariance subset select the same 1,547 heliocentric rows.

The completed baseline and age fits, including their separately optimized dipole nulls, all pass. The independent −2 log likelihood agrees within **9.1e−12**, and independently profiled `(M,X,C)` within **4.3e−14**. The audited numerical snapshot is:

| C2 heliocentric fit | q_m | q_d | −2 log L | Dipole improvement over fitted isotropic null |
|---|---:|---:|---:|---:|
| Baseline dipole | .009517294 | −31.773406 | −184.912666612 | 36.890477741 |
| Baseline isotropic null | .028723577 | 0 | −148.022188872 | — |
| Age-corrected dipole | .353990601 | −32.133725 | −184.735623741 | 36.940917089 |
| Age-corrected isotropic null | .372890318 | 0 | −147.794706652 | — |

Both dipoles reach S=.00938 at the lower allowed boundary. These calculations reproduce the central numerical behavior of the released-input C2 model, including its distinct response from C1. They are **pointwise likelihood checks**, not independent proofs of global optimality or calibrated dipole significance. The input `c2-fits.json` SHA-256 is `c8eef00d771d2783bcbba6c67c26570fccec1d14b8d5243b56f7b6fd95a58771`. Full parameter values and this hash are retained in the audit output. After the root's authoritative 19-fit cosmology replay and final diagnostics completed, all three audit scripts were rerun against the final inputs. Every numerical audit output remained byte-identical to its pre-replay snapshot; all **58 recorded input/output hashes** passed. The final C2 snapshot also remained unchanged. Any subsequent revised fits require explicit comparison with these exact numbers and reexecution of `c2_crosscheck.py`.

The auditor independently read [Ray et al. 2607.20570v2](https://arxiv.org/abs/2607.20570v2), including its methods and discussion, and verified the archived arXiv metadata: version 2, submitted **3 August 2026 at 17:55:09 UTC**. Its title is now *Accelerating expansion and isotropic sky-hemisphere consistency in Pantheon+ supernovae: a revised analysis in the dark energy debate*. It reports q=−.490 before and −.267 after the supplied age correction, replacing its earlier near-zero/strongly-decelerating numerical claims. It explicitly acknowledges the Galactic-as-equatorial coordinate error, **and states that this cannot explain its failed reproduction of the original full-sample baseline**; another unidentified pipeline difference remains possible. The directional investigator's interpretation of this primary revision is correct.

Ray v2 fixes J=1, uses a uniform .15-mag error, and fits standardization coefficients. This is materially different from C2's population likelihood and full covariance. Its H0 and magnitude intercept are exactly degenerate at this data level. The source does not unambiguously specify which cosmological-redshift column generated its printed numbers. The directional branch's simplified reconstructions remain accelerating but do not exactly match them. Reading the revision verifies the changed claim; it does not certify the printed cosmological numbers or establish an independent observational confirmation of isotropy. Its several decompositions reuse the same observations, and similar hemispheric mean coefficients need not exclude a local, redshift-dependent dipole.

An independent fixed J2000 rotation, without calling Astropy, maps Galactic `(264°,48°)` to equatorial `(167.7866166°,−7.1453929°)`. Dot products computed wholly in Galactic or wholly in equatorial coordinates agree within 2.2e−16. On the 1,564-row hemisphere sample:

| Axis interpretation | Positive/negative hemispheres |
|---|---:|
| Mistaken Galactic numbers as equatorial `(264°,48°)` | 724 / 840 |
| Properly rotated direction | 538 / 1026 |
| Ray's rounded equatorial `(167.8°,−7.1°)` | 538 / 1026 |
| Sah main-analysis code's rounded `(168°,−7°)` | 539 / 1025 |

Three near-boundary rows change sides between the last two rounded axes: `2007is` (survey 65), `PS15cwx` (150), and `2007hu` (65), producing a net difference of one. The [Sah response](https://arxiv.org/abs/2608.02484) prints **539/1025 while labelling its axis `(167.8°,−7.1°)`**. The alternate rounded main-code axis reproduces that count, providing a concrete plausible convention explanation, but the response's printed direction/count pair remains inconsistent. This does not invalidate the much larger Galactic/equatorial error, nor establish which exact response code was executed.

## What survived, what did not, and decisive missing evidence

- The central conditional distance shifts survive an independent numerical implementation. They are not evidence that the particular astrophysical correction is uniquely identified.
- The fully corrected matched-sample residual association remains weak in the tested summary-data models. This neither rejects every slope near −.03 nor proves that cosmologically relevant population evolution has been removed.
- The C26 mock's claimed mathematical necessity fails, while the possibility of absorption by an actually fitted cosmology survives. An injection under the real pipeline is needed to distinguish them.
- The scalar host-to-progenitor conversion cannot by itself select either a large or negligible correction. Stable affine transport is a valid special case; population-dependent mappings, selection and inconsistent moments remain consequential alternatives.
- The Roy public correction is exactly identified; its equivalence to Hoyt's private executed correction is unresolved because the source signs differ. This should remain in the final discrepancy ledger.
- Dust/age causal decompositions and covariance-swapping diagnostics were read for model and matrix-subsetting logic. No new source-data dust fit or independent full DES-Dovekie numerical reproduction was performed by this auditor. The inverse-to-covariance-before-subsetting order in the standardization script is appropriate.

The most informative missing products remain author-specific matched-age/error crosswalks and host-property posterior samples; C26 mock age-PDF/run inputs and mapping tables; selection-aware W26/TITAN host-delay products; Son's executed correction/clock/uncertainty/likelihood bundle; and Hoyt's executed corrected magnitude vector. An independent, selection-calibrated age/dust relation predicting held-out distances across redshift would add empirical information. More fits to the same corrected distances alone cannot supply that information.

## Reproduction

Preregistered checks: `docs/experiments/independent-audit-plan.md`. Run in the pinned environment:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/audit/independent.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/audit/directional_crosscheck.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/audit/c2_crosscheck.py
```

Machine-readable calculations and input/code/output hashes are in `runs/audit/independent-results.json`, `runs/audit/manifest.json`, `runs/audit/directional-crosscheck.json`, `runs/audit/directional-manifest.json`, `runs/audit/c2-crosscheck.json` and `runs/audit/c2-manifest.json`. The scripts do not import the investigated cosmology/mapping/age/directional implementations. The source sign discrepancy is a primary-text audit, not an inferred numerical failure. Original acquisitions were preserved; no contact or publication occurred.
