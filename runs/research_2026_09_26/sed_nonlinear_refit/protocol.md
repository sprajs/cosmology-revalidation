# Six-object nonlinear nuisance-refit gate

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
coordinates <1e-4, retained rank4, no near-boundary parameter, and agreement
between starts in objective within1e-6 and fixed-Tripp response within1e-4 mag.
The native noiseless nominal fit must have chi-square <1e-10 and each fitted
coordinate within1e-6 of zero. These checks do not prove a global minimum.

Compare each alternative's actual nuisance shift relative to that target's
nominal optimum with two tangent predictions there: an infinitesimal mean
perturbation and the finite mean displacement. Both use the freshly evaluated
nominal Jacobian and exact fixed C. Standardize with (1,.16087,-3.1178,0).
Report individual responses and a descriptive highest-two minus lowest-two
six-object contrast, not the full1020 high255/low255 contrast or a cosmology fit.

Quantify native-minus-official mean in quoted-error and exact-C units. Record
all nonpositive native or official baseline epochs. If the fixed six contain
any, repeat affected fits on the fixed positive-baseline subset using the
**marginal covariance submatrix**, preserving the full-mask primary result.
Never clamp spectra or model-band fluxes. The separate eight-object native-J
audit remains a distinct check and does not establish native equivalence on
all1020 objects. No coefficient or cohort changes follow outcome inspection.

Results and a decision about expanding this gate follow below.
