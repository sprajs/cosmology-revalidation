# Preregistration: host age, progenitor delay, and population transport

Registered 2026-09-20 before executing calculations or inspecting their outputs. Investigator: independent mapping/physics Astra Ultra branch. Published numbers and source methods have been read; other investigators' reports have not been read. This is an analytic and controlled-simulation investigation, plus a source/data-availability audit. Synthetic experiments below are not reproductions of author data.

## M1: slope times evolution invariance

Question: When does consistently mapping host age H to progenitor delay T preserve a redshift-dependent mean magnitude correction?

Competing predictions: a common deterministic affine map preserves the slope–mean-evolution product exactly; nonlinear maps, latent scatter, measurement noise, changing host mixtures, and selection generally do not preserve it without further assumptions. Latent scatter alone need not break the result when the conditional mean map is affine, invariant between populations, and the magnitude response is linear.

Method: derive the population covariance identities, then use seeded synthetic low/high-redshift populations. Test affine mapping, independent classical age noise, stochastic conditional delay, deterministic nonlinear mapping, changing mapping intercept, and selection on brightness. Compare known latent population mean luminosity change with predictions using host-age and progenitor-age fits. Explicitly distinguish fitted conditional expectation from random progenitor-age draws. Fit low-redshift slopes only; transport to high-redshift population. Use sample sizes 200,000 per population (or exact finite-support quadrature when advantageous), seed 20260920, and independent Monte Carlo error diagnostics. Fix magnitude-age coefficient at -0.03 mag/Gyr for scale illustration, not an empirical prior.

Diagnostics/interpretation: exact affine invariance must agree to floating-point tolerance; stochastic cases judged against analytic expectation and Monte Carlo uncertainty. Failures demonstrate absence of a mathematical guarantee, not that the failing scenario describes actual SNe. Also test mean versus median mixtures and show that identical median evolution can coexist with different mean luminosity evolution.

## M2: CSFH–DTD equation reconstruction

Question: What ages follow from published CSFH–DTD prescriptions, which summary statistic is transported, and how much is convention-dependent?

Inputs: C14 equations and appendices, S25 published Section 3/Figure 2, W26 Appendix A, C26 published Section 3/Figure 3, and public formula implementations if identified. Evaluate the C14 smooth DTD, and truncated power laws (300 Myr/-1, 40 Myr/-1.13), with Behroozi 2013 CSFH only if its exact adopted formula/grid can be established. Evaluate a clearly labelled Madau–Dickinson sensitivity independently. Include mass-loss convention where relevant. Flat LCDM (H0=70, Om=0.3) and Son's BAO+CMB CPL parameters (Om=0.353, w0=-0.42, wa=-1.75), with H0 explicitly reported and any ambiguity treated as a limitation. Radiation treatment and cosmic-time conversion recorded.

Outputs: normalized progenitor PDFs, means, medians, galaxy formed-mass mean age, and slope-times-relative-age curves for z=0..2. Numerical convergence by doubling integration/grid resolution and independent LCDM analytic age check. Record whether 5.3 Gyr S25 evolution is reproduced, approximate, or blocked; do not tune undocumented parameters to force agreement. A cosmic volume is not a survey-selected population. No inferred cosmological probability comes from this test.

## M3: public W22/W26 simulation reproducibility audit

Question: Does the supplied public repository contain enough identified input/configuration to reproduce W26's selected host/progenitor curves?

Method: audit source, fixed repo commit, supplied library/schema, DTD implementation, host rate weighting, mass functions and selection implementation; run a genuine derived-summary calculation from author-supplied model products only if the weights and columns can be interpreted. Exact W26 reproduction requires its figure-specific config, seeds, selection functions and output correspondence. Any failed run or unavailable path is recorded, never replaced with synthetic input under the same label.

## M4: causal/source assessment

Trace host integrated attenuation versus line-of-sight extinction; latent metallicity/SFH and age–dust–metallicity SED degeneracy; marginal versus partial correlations; selection and calibration; correction stacking. Evaluate strongest contrary evidence for all interpretations. Record initial independent conclusions before reading any peer report. Full TITAN reproduction is attempted only if the cited SFH posterior data are released; a paper statistic is not a locally replicated dataset.

## Decision log

- Initial registration: controlled models above deliberately isolate mathematical transport assumptions. Their success/failure cannot adjudicate their empirical prevalence. Acquisition may refine exact source formulas before evaluation; any such choice must be recorded here before outcomes.
- Source refinement before calculation: B13 Appendix F equation F1/Table 6 specifies CSFR=C/[10^(A(z-z0))+10^(B(z-z0))], with z0=1.243, A=-0.997, B=0.241, C=0.180. Its systematic uncertainty (Table 7) is not negligible and will not be fit away. C14 equation 3 is a smooth DTD, not a hard delay cutoff. Park Appendix A3 explicitly adopts alpha=20, tp=0.3 Gyr, s=-1 following Kang; use this disclosed same-group convention for the S25-like reconstruction, with exact S25 code equivalence unproven. Treat the 300 Myr truncated power law as a separate W26-style approximation.
- M3 source refinement before summaries: supplied W22 host library has mean_age in Myr and SN_age in Gyr. The inspected generation script samples mass with rate×mass-function probabilities, then SN delay from its host SPAD. Summarize each redshift's supplied rows with equal weight; multiplying by pred_rate_total again would double-weight the SN rate. These summaries describe the supplied host library, without a recovered W26 selection function or figure linkage. Record redshift support instead of extrapolating it to zero.
- M2 H0 clarification before outputs: DESI Table 5 gives H0=63.6 for the BAO+CMB CPL parameter combination quoted by Son. Include CPL at both H0=63.6 and 70. A fitted SN intercept absorbs the distance zero point, but it does not remove the dependence of physical age-in-Gyr corrections on H0; exact S25 treatment remains to be established.
- 2026-09-20 amendment after initial M1–M3 outputs: extend M2 redshift grid to 2.5 to cover the parent's full Pantheon+ sample, keeping the same formulas. Quantify H0 dependence by central finite differences at 70 +/- 1 percent; no retuning to a published curve. Record the discovered mismatch to W26's ~1.5 Gyr cosmic-CSFH description; do not substitute a different CSFH silently.
- 2026-09-20 amendment before M4 numerical counterexample: test Park's claim that age-correction/removal of mass-step asymmetry establishes age causation. Generate an unobserved common driver U that causes host age, host mass, and luminosity, with no causal arrow age->luminosity. Fit a continuous age correction and a binary mass step on the same synthetic sample. Any asymmetry in this deliberately noncausal model disproves the uniqueness of the causal inference, not the possibility of age causation in data. Use independent seed 20260921, N=200000, H=5+1.5U+N(0,.2), M=10+.5U+N(0,.2), Y=-.05U+N(0,.05), U standard normal.
- M3 audit decision after first source-summary output: retain exact floating-point redshift groups in the descriptive CSV, flag near-duplicate redshifts and nonuniform group sizes. Do not automatically combine physically unexplained fractional low-z values. Exclude zero host-age variance groups from slope summaries with explicit NaN rather than a runtime division warning. Detailed summaries are acquisition/product diagnostics, not recovered W26 outputs.
