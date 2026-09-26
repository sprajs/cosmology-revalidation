# Discovery-conditioned shared calibration/SALT prediction

This is a post-validation exploratory attribution analysis, motivated by the
successful frozen observer prediction and by release-provenance evidence. It
is not a new pristine test or a correction to the original validation tail.
Extend all9 actual coupled native variants from the verified43 cohort to the
unchanged1020 validation objects, holding each object's nominal exported
x0/x1/c/t0 fixed with LFIXPAR_ALL. Match every accepted epoch. Require exact
parameter, epoch, band, observed-flux and quoted-error equality between each
variant and nominal before differencing model means. If data units differ,
stop and establish their map explicitly; do not subtract unlike flux units.
Nominal frozen C and J define one common projection. Retain9 signed columns
with amplitude scale .3, no normalization, and do not add any distance Csys.

H_cal: y_s=A_s a+noise, A=.3*[projected coupled native differences], a~N(0,I9).
H_cal+observer: y_s=A_s a+.02*T_s b+noise, b~N(0,I3), independently of a before
discovery. The observer prior is the same declared .02-mag scale, with the
actual three-contrast basis, not an independently added observer-only posterior.
For each hypothesis learn the joint latent posterior solely from discovery43.
The calibration coordinates are common between discovery and validation.
Retain all Gram and cross-Gram terms, particularly A^T T. Do not fit or select
additional amplitudes, directions, scales or subsets on validation.

For prior-whitened X, F=X^T X and g=X^T y, define
R(F,g)=.5*[g^T(I+F)^(-1)g-logdet(I+F)].
The joint validation predictive gain relative to fixed zero residual is
R(F_D+F_V,g_D+g_V)-R(F_D,g_D).
Compare the two gains; integrate shared coefficient uncertainty once across
the validation sample. Independently verify this with the conditional Gaussian
mean X_V*(I+F_D)^(-1)g_D and covariance
I+X_V*(I+F_D)^(-1)*X_V^T using low-rank/orthogonal algebra.

Also report the old fixed observer direction as a descriptive conditional
shift after centering the calibration-null predictive residual and using its
shared predictive covariance. Do not reuse zero-mean independent-field sign
probabilities, since calibration creates cross-field dependence. A global
joint fit may be saved descriptively but cannot replace held-out prediction.

The nine finite variants and .3 distance-covariance convention yield a useful
Gaussian second-moment approximation, not the full upstream latent posterior.
CALSPEC and other existing modes are omitted here, realized draws imperfectly
cover the Fragilistic DES marginal, and the release's exact seeds/mapping,
wavelength covariance, nonlinear refits, selection and BBC are not reproduced.
No outcome alone identifies a physical new correction or rejects the complete
published systematic budget. All discovery/validation numerical ingredients
and native code/data gates will be retained for independent review.

## Prior-coordinate sensitivity, frozen before shared-validation scores

The main inherited prior is unchanged. Before reading any shared-validation
predictive scores, add one declared sensitivity: an isotropic physical zero-sum
griz prior with the same total prior squared magnitude. For the original gauge
matrix B=[(1,0,0),(-1,-1,-1),(0,1,0),(0,0,1)] and tau=.02, use
Cov(c)=2*tau^2*(B^T B)^(-1). Both priors have trace[B Cov(c) B^T]=6*tau^2;
the new prior has equal marginal variance for each physical band. Preserve
calibration modes and all other choices. Integrate this prior jointly with the
same discovery-trained calibration coordinates, report both alternatives and
do not select the larger predictive score. This is a sensitivity within the
already post-validation exploratory attribution analysis.
