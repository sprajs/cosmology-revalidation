# Joint ΛCDM measurement and the absolute calibration

Under spatially flat ΛCDM, the released CMB, BAO and Dovekie supernova data give **H₀ = 68.164 ± 0.251 km s⁻¹ Mpc⁻¹**. Replacing Dovekie with the calibrated Pantheon+SH0ES release moves the joint result to **68.514 ± 0.242**. The calibrated and uncalibrated supernova information therefore do not simply produce one mutually consistent absolute scale. A separate prediction that excludes all 77 calibrator measurements places their observed common calibration contrast in an extreme predictive tail under this model.

These are **qualified results for the declared CAMB accuracy-1 target**, with the original importance weights retained. Qualification of the sampling and likelihood replacement is distinct from convergence when native numerical accuracy is increased. That numerical-accuracy qualification remains outstanding here; the separate 32-point screen and any required follow-up must be assessed before treating the quoted precision as final.

## Three separate supernova choices

All three rows below use the same CMB and DESI DR2 BAO factors. Supernova releases are alternatives, never multiplied together. Values are posterior means and standard deviations, not errors on the posterior mean or symmetric credible intervals.

| Supernova factor | H₀ [km s⁻¹ Mpc⁻¹] | Ωₘ | Raw weight effective sample size |
|---|---:|---:|---:|
| Dovekie, free common magnitude | 68.1638 ± 0.2506 | 0.303727 ± 0.003399 | 1,783.3 |
| Pantheon+SH0ES, all 1,657 selected rows | 68.5139 ± 0.2421 | 0.299096 ± 0.003293 | 432.3 |
| Pantheon+SH0ES, 1,580 noncalibrator rows only | 68.1882 ± 0.2522 | 0.303389 ± 0.003419 | 1,806.3 |

The [Dovekie measurement](../results/inference/modern-lcdm-independence-measurement.json), [calibrated replacement](../results/inference/anchored-bridge-lcdm.json) and [noncalibrator replacement](../results/inference/calibration-holdout-lcdm.json) retain full weighted covariance and asymmetric quantiles. Each uses 2,000 selected native-evaluated parent slots. Their Pareto-k diagnostics are 0.353, 0.296 and 0.314, respectively, and all prescribed weight, independent-chain and batch-stability gates pass. The replacements use the original, untrimmed exact/proposal weights multiplied by the normalized new-to-old supernova likelihood ratio.

The noncalibrator replacement changes the mean H₀ by only +0.0244 relative to Dovekie. Including the calibrators changes it by +0.3501. These are paired changes using overlapping observations and shared CMB/BAO information; they are not differences between independent measurements. The calibrated joint posterior also must not be described as agreement with the supernova-only calibration result, [H₀ = 73.5496 ± 1.0173](calibration-lcdm.md). That latter result conditions on a different data combination and the specified supernova-only background and priors.

For the Dovekie joint fit, q(0) = −0.54425 ± 0.00510, q(0.5) = −0.10647 ± 0.00581 and q(1) = +0.16624 ± 0.00417. Thus this fitted ΛCDM history accelerates now and decelerated at redshift one. The form w = −1 and spatial flatness are assumptions of this comparison. These narrow conditional intervals neither establish ΛCDM adequacy nor exclude more flexible expansion or luminosity-evolution models. No empirical progenitor-age correction has been inferred or applied.

## Predicting the withheld calibrators

The noncalibrator posterior uses CMB + BAO + the 1,580 noncalibrator rows only. All 77 calibrator values are excluded from its weights. Their prediction retains the released calibrator-to-noncalibrator covariance and integrates the common supernova magnitude inferred from the noncalibrators. The [conditional Gaussian derivation](calibration-holdout.md) specifies a covariance-defined common calibration contrast before evaluating the observations.

The predictive CDF is averaged over the noncalibrator posterior first. Its two-sided equal-tail area is **2.5614 × 10⁻⁸**, with log area −17.4801. This is a Bayesian predictive diagnostic under the released Gaussian likelihood and the assumed cosmological model. It is **not a calibrated frequentist p-value, a Gaussian “sigma” significance, a Bayes factor, or the probability that ΛCDM is true**. The common contrast tests one specific residual direction; it does not test every possible mismatch among the 77 rows.

The smaller-tail integral has contribution effective sample size 395.8 and relative Monte Carlo error 0.0698 estimated between the four independent chains. Its 10-, 20- and 40-batch estimates are 0.0443, 0.0454 and 0.0442. The proper 77-dimensional log predictive density is 38.2699, expressed with respect to magnitudes in all 77 coordinates; its contribution effective sample size is 432.3 and its independent-chain relative error is 0.0657. Both predictive quantities pass their predeclared precision checks. A positive log density in dimensional coordinates is not itself a goodness-of-fit probability.

Within these assumptions, the withheld absolute calibration is difficult to predict from the CMB, BAO and noncalibrating supernovae. The result identifies a conditional discrepancy, not its physical cause. Possible changes to calibration, the released covariance, luminosity standardization, early-universe assumptions or late-time expansion have not been distinguished by this test.

## Assumptions and unresolved limits

The background fixes zero curvature, w = −1, one massive neutrino with total mass 0.06 eV and N_eff = 3.044. The primordial and matter parameters retain the declared broad scalar priors, including a uniform H₀ prior from 40 to 100 and a uniform log-amplitude prior. The four Gaussian CMB calibration priors are retained once; their treatment and the primary/lensing cross-covariance approximations are recorded in the [external-probe specification](external-probes.md) and [dependence audit](probe-dependence.md). The sound horizon is calculated with the same early-universe model, not independently freed or supplied as an extra prior. No residual grey luminosity evolution is allowed in these three results.

Pantheon+SH0ES supplies the calibrator distance moduli and the complete statistical-plus-systematic covariance. Calibrators use their released Cepheid distances rather than redshift distances. No additional Cepheid factor or H₀ prior is added. This follows the released construction and [official likelihood](https://github.com/PantheonPlusSH0ES/DataRelease/blob/c447f0fea703fcd0fff57de5000947b5ca81286b/Pantheon%2B_Data/5_COSMOLOGY/cosmosis_likelihoods/Pantheon%2BSH0ES_cosmosis_likelihood.py), using the pinned [distance and covariance release](https://github.com/PantheonPlusSH0ES/DataRelease/tree/c447f0fea703fcd0fff57de5000947b5ca81286b/Pantheon%2B_Data/4_DISTANCES_AND_COVAR). The free-magnitude integral uses the same improper measure in all relevant factors; no absolute model-evidence normalization is asserted.

The [calibration-interface audit](calibration-interface.md) and [independent covariance analysis](calibration-covariance-estimands.md) leave the precise shared-host covariance construction unresolved. Several distinct-supernova sibling terms are close to twice independently recovered Cepheid-only host variances; that pattern does not establish a duplication error or license a factor-of-two repair. Upstream photometry, selection and covariance also already contain shared calibration information. Holding out likelihood rows therefore does not recreate a wholly independent blind reduction of raw observations. Those limitations are especially material to interpreting the very small predictive tail.

## Numerical evidence and retained failures

The new supernova calculations reproduce each stored source likelihood to at most 3.87 × 10⁻¹² in log likelihood; density replacement closes to 1.14 × 10⁻¹³. The [independent conditional-calibration review](../results/inference/calibration-holdout-review.json) checks the Gaussian factorization, the calibration-shift sign, full-covariance normalization and predictive integration on synthetic systems and fixed positive distance curves. These validate the equations and implementation, not the empirical covariance construction.

The original finite postprocessing stopped because a strict chronology test rejected two legitimate repeated selection slots. The [independent reconstruction](../results/inference/selected-slot-chronology-lcdm-review.json) recovers every selected index, source chain row, parameter point and stored density from the frozen random seed. Flooring draws from adjacent real-valued strata can select the same integer chain position. Their multiplicity must remain. Neither repeated pair crosses the chronological training boundary or the 10/20-block bootstrap boundaries in this cohort. The stopped receipt and proposal-freeze refusal are preserved; no new proposal or sampling run was created.

The explicit [postprocessing continuation review](../results/inference/postprocessing-continuation-review.json) confirms that the original failure and completed outputs remain hash-bound and that native precision points are not silently retried. Completion of that orchestration is not scientific qualification. The separate native-accuracy assessment remains required for these declared-target results.
