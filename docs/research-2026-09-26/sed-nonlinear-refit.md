# Six-object nonlinear nuisance-refit gate

**The bounded gate passed: the constructed ambiguity survives nonlinear
nuisance fitting in all six selected objects.** The observer and SED
alternatives leave almost identical changes in the fitted flux predictions
(cosine .9969), while their highest-two minus lowest-two fixed-Tripp responses
are +.04557 and -.14799 mag in the noiseless test. All 108 optimizer runs
passed the prespecified gates. This warrants expanding the same frozen
experiment to the full validation cohort; it does not validate that cohort's
previous tangent-order .24-mag contrast or identify a physical correction.

## Protocol frozen before fitting

Select six validation objects at zero-based ranks 85,255,425,595,765,935 after
sorting by (zHEL,CID): the midpoints of six equal-population redshift strata.
The metadata-only selection and source hashes are in
`runs/research_2026_09_26/sed_nonlinear_refit/freeze.json`. No flux outcomes
inform the cohort, perturbation family, coefficients, starts or numerical gates.
This samples central redshift support rather than the extreme tails.

Use the unchanged 52-component broad rest-frame SED plus log-redshift drift
with the previously saved `broad_drift_discovery_rms0.05` coefficient column3.
The observer comparator is the original frozen discovery griz vector. Both
alter the predicted mean; observations and exact frozen SNANA covariance C
are held fixed. The SED multiplier is exp[-kappa Delta_m(phase,rest wavelength,z)]
inside the actual photon-weighted passband integral. Phase follows fitted t0.
The observer multiplier is exp[-kappa dmag_band] outside that integral. The
four nuisance coordinates are delta(mB,x1,c,t0), with x0 multiplied by
exp[-kappa delta_mB]. Keep all signed native spectral and integrated flux.
No per-epoch official/native rescaling is used as a physical spectrum.

Fit three distinct targets with nominal, observer and SED mean models:

1. **Native noiseless closure:** target is the identical independent sncosmo
   SALT3 implementation at the reference parameters. This isolates nonlinear
   parameter compensation and must recover zero shift/zero objective for the
   nominal model.
2. **Official-mean sensitivity:** target is the exported official mean. This
   quantifies the independent/native mean difference without assigning a
   physical SED to per-row normalization factors.
3. **Observed-flux diagnostic:** target is the actual accepted flux array.
   This is separate from the constructed noiseless experiment and does not
   estimate the truth or prior probability of either perturbation.

Retain the exact accepted masks and full covariance. Use two prespecified
starts, (0,0,0,0) and (.05,.5,.03,2), relative to the reference coordinates.
Bounds are delta_mB in [-1,1] mag, delta_x1 in [-3,3], delta_c in [-.3,.3]
and delta_t0 in [-10,10] observer days. A boundary hit fails the local gate.
Use least_squares, a central explicit Jacobian with steps
(1e-4,1e-3,1e-4,.01), analytic amplitude derivative, x_scale=(.1,1,.1,3),
xtol=ftol=1e-11, gtol=1e-9 and max_nfev=300. Verify the Jacobian at half steps.
For each fit require optimizer success, whitened gradient norm in Fisher
coordinates < 1e-4, retained rank 4, no near-boundary parameter, and agreement
between starts in objective within 1e-6 and fixed-Tripp response within 1e-4 mag.
The native noiseless nominal fit must have chi-square <1e-10 and each fitted
coordinate within 1e-6 of zero. These checks do not prove a global minimum.

Compare each alternative's actual nuisance shift relative to that target's
nominal optimum with two tangent predictions there: an infinitesimal mean
perturbation and the finite mean displacement. Both use the freshly evaluated
nominal Jacobian and exact fixed C. Standardize with (1,.16087,-3.1178,0).
Report individual responses and a descriptive highest-two minus lowest-two
six-object contrast, not the full 1020 high255/low255 contrast or a cosmology fit.

Quantify native-minus-official mean in quoted-error and exact-C units. Record
all nonpositive native or official baseline epochs. If the fixed six contain
any, repeat affected fits on the fixed positive-baseline subset using the
**marginal covariance submatrix**, preserving the full-mask primary result.
Never clamp spectra or model-band fluxes. The separate eight-object native-J
audit remains a distinct check and does not establish native equivalence on
all 1020 objects. No coefficient or cohort changes follow outcome inspection.

Results and a decision about expanding this gate follow below.

## Executed result

The six selected CIDs are1927242,1301936,1660018,1341459,1531033 and 1319898,
covering zHEL=.24666,.34351,.43550,.52341,.61002,.72925 and 226 accepted
epochs. Source parameters and exact objective arrays provide the unrounded
redshifts used in integration. There are 54 target/model/object fits, each
run from both starts: 108 optimization runs. No row was removed.

Nominal noiseless fits return exactly zero parameter shift and zero objective.
Every fit retains rank 4 and lies comfortably inside its bounds (minimum
fractional boundary margin .4257). The maximum whitened gradient norm in
Fisher coordinates is 1.53e-6; the largest Jacobian half-step discrepancy is
4.36e-7 relative. Both starts agree to 1.23e-7 mag in fixed-Tripp response and
1.32e-11 in objective. The maximum number of function evaluations is 15.
These checks establish local numerical reliability, not global optimization
or calibrated parameter uncertainty.

The actual nonlinear responses in the native noiseless experiment are:

| CID | zHEL | Observer delta(mu_ref), mag | SED delta(mu_ref), mag | SED nonlinear minus finite-mean tangent, mag |
|---|---:|---:|---:|---:|
| 1927242 | .24666 | -.013865 | +.002070 | +.000526 |
| 1301936 | .34351 | +.000723 | -.001744 | +.000439 |
| 1660018 | .43550 | +.023199 | -.002067 | +.000692 |
| 1341459 | .52341 | +.028975 | -.037466 | +.000008 |
| 1531033 | .61002 | +.035583 | -.101126 | -.002217 |
| 1319898 | .72925 | +.042418 | -.194522 | -.006011 |

Here delta(mu_ref) means compensation from changing the mean model while
holding the target and C fixed, using alpha=.16087 and beta=3.1178. It is
not a physical distance measurement. Both alternative models move their
rest phase with fitted t0. No official/native flux-ratio field is treated
as a physical spectral change.

For the descriptive six-object highest-two minus lowest-two contrast:

| Target branch | Observer nonlinear, mag | SED nonlinear, mag | SED finite-mean tangent, mag |
|---|---:|---:|---:|
| Native noiseless | +.045572 | -.147987 | -.143390 |
| Official-mean sensitivity | +.045573 | -.147985 | -.143390 |
| Observed-flux diagnostic | +.044321 | -.149192 | -.143403 |

The SED-minus-observer contrast is -.193558 mag noiseless and -.193513 mag
in the observed-flux diagnostic. Nonlinear-minus-finite-tangent differences
reach .00601 mag noiseless and .00627 mag for observed data. Infinitesimal
and finite-mean tangent predictions are both retained in the result: their
errors need not have the same sign because finite SED and nonlinear SALT
compensation can partially cancel.

The separate observed-flux branch optimizes nuisance parameters against
measured flux only after freezing the cohort and perturbations. It does not
fit SED coefficients, their redshift evolution or their physical probability,
and its objectives are not used to declare a winning physical mechanism.

## Do the flux patterns remain equivalent after fitting?

Replay each saved nonlinear mean and subtract the same target branch's
nominal fitted mean. Stack these changes using the exact C whitening. The
observer and SED changes have norms 1.14178 and 1.08617 in the native noiseless
test; their difference has norm .103883 and cosine .996896. The relative
squared mismatch is .008278, so the constructed SED response retains about
99.17% agreement by this response metric in the fixed six-object sample.
The observed-data branch gives difference norm .104287, cosine .996863 and
relative squared mismatch .008335. This is comparison of fitted mean
responses, not explained variance of measured residuals or a p-value.

The replay reproduces all stored fitted means and objective values exactly
at saved precision. Together with the large opposite-sign standardized
responses, it confirms that the local degeneracy did not disappear merely
because four SALT parameters were allowed to refit nonlinearly.

## Native mean, early-phase support and extension boundary

The explicit integration reproduces cached independent native means to
relative 1.02e-15. Native-minus-official mean norm ranges .00203–.01304 in
exact-C units; the maximum individual difference is .01328 quoted-error
units. Replacing the noiseless target with the official exported mean changes
individual nonlinear responses by at most 2.63e-6 mag for observer and
3.29e-6 mag for SED. This establishes small mean sensitivity for these six
objects only.

All 226 nominal integrated means are positive in both implementations, so
the prespecified early-negative row-removal sensitivity has no affected
rows in this cohort. Signed spectral flux is still used inside every
integral; positive band-integrated means do not establish positive spectral
support. The separate
[eight-object native-J audit](../../runs/research_2026_09_26/astra_design/negmean/result.json)
targets the ten previously flagged early-phase validation epochs: replacing
the affected Jacobians changes the global fixed observer gain 27.233214 to
27.229805; marginally removing those ten epochs with native Jacobians gives
27.219293. That specific edge does not drive the existing transfer score.
It neither substitutes for this nonlinear SED gate nor proves native
equivalence over all 1020 objects.

**Expanding the frozen nonlinear experiment is warranted.** Retain identical
coefficients, accepted masks, covariance, targets, two-start rules and
failure gates. Full 1020 refits are needed to verify the actual high255/low255
contrast, redshift tails and weak-tangent objects. Record all failures rather
than selecting successful responses. A larger gate would quantify this
constructed ambiguity; physical population support, shared training/calibration
uncertainty, selection, classifier/BBC and any cosmological interpretation
remain separate requirements.

Artifacts: [executed script](../../scripts/research_2026_09_26/sed_nonlinear_refit.py),
[frozen cohort and hashes](../../runs/research_2026_09_26/sed_nonlinear_refit/freeze.json),
[full result and all optimizer gates](../../runs/research_2026_09_26/sed_nonlinear_refit/result.json),
[paired responses](../../runs/research_2026_09_26/sed_nonlinear_refit/paired-responses.csv),
[saved-fit replay](../../runs/research_2026_09_26/sed_nonlinear_refit/saved-fit-verification.json),
and [manifest](../../runs/research_2026_09_26/sed_nonlinear_refit/manifest.json).
The script refuses to overwrite its completed result; a fresh run directory
or an explicitly revised output path is required for a new execution.
