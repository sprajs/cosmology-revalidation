# Forward population and selection checks

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../README.md) for current execution status.

The first actual comparison of all nine newly generated, calibrated-flux-refitted variants is complete. It finds clear luminosity/noise dependence of selection, and some observable discrepancies, but **does not establish a physical dust/population winner**. These are public-assets variants, with the historical population prescriptions already tuned using overlapping DES observations. The heldout calculation is a retrospective predictive check.

The [preregistration](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/forward_discrimination/preregistration.json) was recorded before opening fitted-model distribution differences. [Numerical falsifiers](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/forward_discrimination/numerical-validation.json) pass: independent double-sum CRPS, global magnitude-offset invariance, no test-magnitude contamination of training alignment, unsupported-field rejection, and paired selection-flip variance.

## Generated-attempt selection response

Each arm has exactly **26,518 generated attempts**, including failures. The common redshift/host/cadence design stream is verified. CIDs are reused after early failure, so pairing uses the unique generated-attempt index. Measured basic cuts are applied after successful flux fitting; none of these counts is the complete historical BBC sample.

| Variant | Written | Fitted | Basic quality |
|---|---:|---:|---:|
| P21 | 2,678 | 1,825 | 1,314 |
| BS21 | 3,129 | 2,153 | 1,560 |
| G10 | 3,085 | 2,164 | 1,631 |
| P21 +0.10 mag | 2,492 | 1,686 | 1,193 |
| P21 −0.10 mag | 2,876 | 1,947 | 1,425 |
| P21 +0.20z mag | 2,433 | 1,648 | 1,166 |
| P21 added-noise correlation 0 | 2,696 | 1,837 | 1,346 |
| P21 added-noise correlation 0.9 | 2,669 | 1,806 | 1,282 |
| P21 true noise ×1.2 | 2,682 | 1,764 | 1,227 |

The +0.10 mag injection loses 126 baseline quality survivors and gains five: its acceptance change is **−0.456 percentage points**, paired 95% Monte Carlo interval [−0.541,−0.372]. The −0.10 mag injection changes acceptance by **+0.419 points** [0.340,0.497]. The +0.20z injection changes acceptance by **−0.558 points** [−0.650,−0.466]. These are physical flux injections, not unsupported reweighting of an effectively singular source-luminosity distribution.

True noise ×1.2 barely changes initial written counts but reduces quality acceptance by **−0.328 points** [−0.411,−0.245], or 87/1,314=6.62% of the baseline survivor count. Changing added-noise correlation from the nominal0.6 to0 or0.9 changes quality counts by +32 and−32. This parameter applies to the added error component, not total flux correlation. The LIBID-cluster intervals are nearly identical in this realization; both are sampling diagnostics, not uncertainty in the noise calibration or repeated-seed coverage.

Full stage/bin counts, gained/lost membership, Wilson intervals and paired/LIBID-cluster intervals are in [selection-response.json](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/forward_discrimination/selection-response.json), with [source hashes](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/forward_discrimination/selection-manifest.json). Larger accepted BS21/G10 counts cannot favor those physical models without a valid observed-rate/parent-population likelihood.

## Conditional observed-distribution check

Start with the original 1,635 DES members and published SNNV19>0.999. The additional common basic fitted cuts leave1,062, including213 original fold0 test objects. Matching is on exact field, measured host-mass category below/above10, and a Gaussian redshift kernel of width0.10 truncated at3 widths. Every object must have at least10 nonzero simulated neighbors and effective count8 in **each of nine arms**. This leaves **817 training and202 test objects across all10 fields**. [The support ledger](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/forward_discrimination/primary/support.csv) retains all exclusions. Compound fields absent from the simulation are not silently relabelled.

Field-provenance correction found during the final audit: the real field used here is the already frozen fold table's **first PHOT-epoch field**, whereas simulated field comes from FITRES. The preregistration's assumption that real compound FITRES labels would be rejected is therefore not implemented for real objects whose frozen first-epoch field is simple. This is a disclosed conditioning discrepancy, not a post-score relabelling. A rerun on consistent first-epoch fields or excluding multi-field real objects is required; the saved primary result must remain unchanged.

For magnitude, simulations use actual fitted mB minus generated distance modulus. Real mB is measured relative to a declared coasting reference; seven free linearly interpolated redshift offsets are estimated separately for each model on training objects only. All seven nuisance directions are identified in this finite design. Scores then compare the heldout empirical distributions of fitted colour/stretch and, secondarily, magnitude/colour/stretch. They use the complete finite weighted ensemble energy score and marginal CRPS, not a KDE likelihood or assumed source-colour-plus-dust proxy. Lower scores are better. The model errors are already present in noisy fitted simulations and are not convolved a second time.

| Comparison against P21 | Colour/stretch energy change | 95% bootstrap interval | Magnitude/colour/stretch change | 95% interval |
|---|---:|---:|---:|---:|
| BS21 | −0.0210 | [−0.0602,+0.0076] | −0.0245 | [−0.0733,+0.0051] |
| G10 | −0.0253 | [−0.0765,+0.0110] | −0.0377 | [−0.0978,+0.0004] |
| P21 correlation0 | −0.0194 | [−0.0413,+0.0006] | −0.0162 | [−0.0382,+0.0022] |

All intervals in this table cross zero. G10's marginal magnitude CRPS contrast is −0.00861 mag with an unadjusted interval [−0.02051,−0.00007]; that barely excludes zero among many reported diagnostics and is not decisive evidence. The finite empirical forecast V statistic has simulation-sample-size effects; resampling includes pilot simulation uncertainty, but does not create an infinite-population density or exact conditional law. Wider scientific claims require the registered sensitivity runs and matched classifier selection.

The equal-field t sensitivity favors G10's colour/stretch score, interval [−0.0490,−0.0049], and its magnitude/colour/stretch score, [−0.0606,−0.0187]. That calculation changes to equal field weighting and treats the simulated forecast as fixed; it does **not** include the finite-simulation uncertainty present in the primary bootstrap. The two uncertainty calculations therefore answer different questions, and neither should be silently substituted for the other.

One interpretable discrepancy is the red-colour fraction: the observed fraction with fitted c>0.1 is **8.03 percentage points below P21's matched forecast**, interval [−13.16,−2.43]. The corresponding BS21 and G10 deficits are3.99 and4.54 points with intervals including zero. Mean-colour contrasts alone are less decisive. This could reflect population differences, classifier conditioning, noise, or unreproduced final membership; it does not identify dust causally.

There is also direct counterevidence to calling the nominal simulation fully calibrated. Mean real FITPROB is0.146 below matched P21, interval [−0.205,−0.102]. Correlation0 reduces this discrepancy to−0.070 [−0.124,−0.023]; true-noise×1.2 overshoots to+0.101 [0.044,0.150]. Real fitted magnitude Hessian errors are on average0.0122 mag smaller than P21's conditional forecast. These are outcome diagnostics, not weights used to make agreement. Author-conditioned data masks versus algorithmic simulated clipping and classifier/BBC selection can affect both. The independent flux/noise audit must be considered before estimating a noise parameter from these quantities.

Primary uncertainties use300 reproducible replicates with real-training/test object resampling, shared generated-attempt Poisson weights across all simulation arms, and refitted training magnitude offsets. The [full summary](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/forward_discrimination/primary/summary.json) includes paired equal-field t sensitivities, marginal PIT/68%/95% coverage and all nine contrasts; [the manifest](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/forward_discrimination/manifest.json) pins inputs and output hashes. Intervals are unadjusted, support is held fixed, and shared calibration/model uncertainty and independent-seed coverage are absent.

## Limits and saved continuation

Truth-Ia simulated selection is not the real high-SNN classifier selection. A reconstruction workstream is testing the released classifier, but no validated matched classifier output was used here. Historical valid-bias-grid membership, chi2max16, contaminant likelihood and intersection across systematic variants remain missing. Data author masks and simulated iterative clipping differ. Historical DES tuning leaks into the physical models, and measured host mass is not SN progenitor age. W22 cannot be run exactly without the missing SN_age generator. Free redshift offsets remove a preferred distance law as a scoring target but cannot undo luminosity-dependent detection, as the injection counts explicitly demonstrate.

The requested pause occurs after the primary300-bootstrap comparison. Registered all-DES, recovered-mask, stricter-cut, bandwidth0.075/0.15 and reference-q sensitivities have **not yet run**. Do not promote the primary results to a robustness claim. [HANDOFF.md](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/forward_discrimination/HANDOFF.md) gives exact commands and remaining checks.
