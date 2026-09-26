# Independent review: conditioning shared calibration modes before validation

The correct comparison is **calibration-only versus calibration plus an extra
observer-band term, with both models trained on the same 43 discovery objects
and predicting the same 1,020 validation objects**. Discovery and validation
must share each calibration coordinate. A zero-mean validation covariance
inflation does not implement that conditional null: discovery changes both the
mean and uncertainty of the calibration forecast. The old observer-only
posterior cannot be imported as an independent prior into the joint model.

This is a methodological review, an independent synthetic calculation, and an
independent arithmetic check of the saved observed-data scores. It does not run
new light-curve fits, tune the observed-data models, or replace another agent's
analysis. The 1,020 residuals and their nominal-model
result were already inspected before this calibration attribution was designed.
The new calculation is therefore **retrospective attribution under an explicit
finite-mode approximation**, not a new untouched confirmation or a rejection of
the complete released systematic budget.

## 1. Define the observational quantity once

For object i, let d_i be the accepted observed epoch fluxes, f_i the nominal
model at its fixed nominal fitted coordinates, C_i=L_i L_i^T the frozen nominal
measurement-plus-model covariance, and J_i the same four-column nominal SALT
nuisance Jacobian used in the previous analysis. Let Q_i have orthonormal columns
spanning the complement of L_i^{-1} J_i. Define

\[
 y_i=Q_i^T L_i^{-1}(d_i-f_i),\qquad
 b_{ij}=Q_i^T L_i^{-1}(\bar f_{ij}-f_i).
\]

The bar is essential: every variant prediction must first be expressed in
**nominal observed-flux units**. If a native variation represents the same data
as d'_{ij}=S_{ij}d_i for a known positive diagonal scale, use
bar f_{ij}=S_{ij}^{-1}f'_{ij}; a pure coordinate conversion maps covariance as
C'=S C S^T. For a documented affine transformation d'=S d+t, use
bar f=S^{-1}(f'-t). A data-dependent estimated scale, a noninvertible transform,
or an unexplained flux mismatch does not pass this gate. Known factors or
independent error-scale checks are safer than flux ratios near zero flux.
The covariance actually produced by the changed physical model can also change;
that is a separate likelihood variant, not a unit conversion to insert here.

Use one common C_i, J_i, Q_i and accepted mask for all columns in this local
mean-response experiment. Reprojecting each alternative using a different
Jacobian/covariance would change the likelihood being compared and can remove
its coherent response. Fixed-coordinate differences followed by the nominal
projection account for first-order nuisance absorption. They do not supply a
nonlinear refit, a parameter-dependent-covariance likelihood, or selection
regeneration. A nominal fitted coordinate and Jacobian estimated from the full
object data also make the design an explicitly local, fixed-design
approximation rather than an exact raw-data generative model.

The convention y=d-f means b is **variant mean minus nominal mean**: a positive
coefficient adds b to the predicted residual. Reversing one entire column in
both discovery and validation is a harmless reparameterization under a
symmetric prior. Reversing it only in validation reverses the cross-sample
prediction while preserving both within-sample covariance matrices. That is a
scientifically different and erroneous model.

Stack objects separately into discovery D and validation V. With B_s containing
the nine b_j columns, define A_s=0.3 B_s. No centering, column standardization,
renormalization, per-field duplication, or independently chosen column rotation
is part of this definition. The finite-mode covariance is A A^T=0.09 sum_j
b_j b_j^T. Its rank is at most nine even when there are tens of thousands of
epochs.

## 2. Preserve the original observer prior exactly

The existing source defines T from [G_g-G_r, G_i-G_r, G_z-G_r], where each
G_b is the projected flux derivative for an extra dimming in band b. The
three coefficients c have the prior N(0,tau^2 I_3), tau=0.02 mag. Their griz
equivalent dimming is

\[
 \delta m=\begin{pmatrix}1&0&0\\-1&-1&-1\\0&1&0\\0&0&1\end{pmatrix}c.
\]

This is a zero-sum gauge, **not an isotropic prior in an orthonormal griz
contrast basis**. In physical band coordinates its prior standard deviations
are [0.02, sqrt(3)*0.02, 0.02, 0.02] mag, with the corresponding correlations.
That asymmetry is inherited empirical regularization. An orthonormalized T
requires transforming its coefficient prior too; simply retaining 0.02^2 I
after a nonorthogonal basis change would silently change the hypothesis.

This is verified in `exact43_compare.py` where the three columns are constructed
and `residual64.py` where `prior=I/sigma**2` is applied. Preserve that precise
definition for the requested comparison. The parent subsequently requested a
declared isotropic physical prior sensitivity, frozen before the new shared-mode
score was computed, without choosing between them by validation score. Write B
for the displayed band map. Since B^T B=I_3+11^T, its coefficient covariance is

\[
 C_c=2\tau^2(B^T B)^{-1}
 =2\tau^2(I_3-11^T/4),\qquad
 B C_c B^T=2\tau^2(I_4-11^T/4).
\]

Both the inherited and isotropic physical griz priors have total variance
6 tau^2=0.0024 mag^2. The isotropic prior gives each band standard deviation
sqrt(1.5) tau=0.024495 mag and every pair correlation -1/3. If L_c L_c^T=C_c,
the prior-whitened observer design is T L_c, not T L_c^T. The full 12-coordinate
transform is blockdiag(I_9,L_c); use it on both sides of F and once on g.
The pending implementation and then its saved output passed this algebra check.

## 3. The two models and their discovery posteriors

Introduce independent standard-normal prior coordinates a in R^9 and h in R^3:

\[
\begin{array}{ll}
H_0:&y_s=A_s a+\epsilon_s,\quad a\sim N(0,I_9),\\
H_1:&y_s=A_s a+\tau T_s h+\epsilon_s,\quad
      (a,h)\sim N(0,I_{12}),\qquad \epsilon_s\sim N(0,I).
\end{array}
\]

The physical observer coefficients are c=tau h. Let X_{0s}=A_s and
X_{1s}=[A_s,tau T_s]. Within each hypothesis, there is one shared latent vector
theta for both samples. The nuisance posterior under H_1 is joint; its
calibration-observer cross covariance is generally nonzero even though the
prior blocks are independent.

For either model k and sample s, retain F_{ks}=X_{ks}^T X_{ks} and
g_{ks}=X_{ks}^T y_s. The discovery posterior is

\[
 S_{kD}=I+F_{kD},\qquad
 V_{kD}=S_{kD}^{-1},\qquad
 m_{kD}=V_{kD}g_{kD}.
\]

Consequently the **joint validation** predictive distribution is

\[
 y_V\mid y_D,H_k\sim N\left(X_{kV}m_{kD},\;
 I+X_{kV}V_{kD}X_{kV}^T\right).
\]

In particular, H_0 usually predicts a nonzero mean. The prior block covariance
is C_DV=A_D A_V^T; Gaussian conditioning gives the same mean
A_V A_D^T(I+A_D A_D^T)^{-1}y_D and covariance
I+A_V A_V^T-A_V A_D^T(I+A_D A_D^T)^{-1}A_D A_V^T.
This supplies an independent implementation route and exposes the missing
cross-sample term in separate covariance inflation. These are standard Gaussian
linear-model identities; the weight-space derivation is given in
[Rasmussen and Williams, Chapter 2, equations 2.8–2.9](https://gaussianprocess.org/gpml/chapters/RW2.pdf).

Using discovery twice is avoided by recomputing H_1's joint posterior from its
original calibration and observer priors. Do not first fit the observer-only
posterior, treat that posterior as independent of calibration, and then use the
discovery residuals again. Nor is an independent sum of the observer-only and
calibration-only posterior predictive covariances the posterior of H_1.

## 4. Low-rank predictive score; only 9-by-9 and 12-by-12 solves

Define

\[
 R(F,g)=\frac12\{g^T(I+F)^{-1}g-\log\det(I+F)\}.
\]

Completing the square gives log p(y_s|H_k)=log phi(y_s;0,I)+R(F_{ks},g_{ks}).
Therefore the discovery-conditioned validation log density is exactly

\[
 \ell_k=\log\phi(y_V;0,I)
 +R(F_{kD}+F_{kV},g_{kD}+g_{kV})-R(F_{kD},g_{kD}),
\]

and the requested score is

\[
 \Delta=\ell_1-\ell_0
 =\big[R_{1,D+V}-R_{1,D}\big]
  -\big[R_{0,D+V}-R_{0,D}\big].
\]

Positive Delta means H_1 predicts these validation residuals better under the
declared priors. It is a conditional predictive log-density ratio, or
log conditional Bayes factor within the stated model, not a measured
calibration error or a Gaussian significance. The joint Bayes factor using
discovery and validation would additionally contain discovery's Bayes factor;
do not conflate them. Computing the combined sufficient statistics in this
identity integrates the validation likelihood under the frozen discovery
posterior. It is not permission to tune a prior, family, amplitude or mask
using validation.

Use Cholesky solves and Cholesky log determinants of I+F. Do not explicitly
invert the 35,526-dimensional validation covariance. Store per-object and
per-field sufficient statistics, but sum them before this single shared
Gaussian integral. Summing independently marginalized per-object or per-field
scores corresponds to redrawing calibration for each block. A sequential
factorization with the shared posterior updated after each block can reproduce
the joint score; the independent marginal sum cannot. Per-field matched
contributions remain descriptive, and sequential predictive contributions
depend on the chosen order.

For the saved unscaled joint design U=[A,T], apply the scale matrix
W=diag(1,...,1,tau,tau,tau) to its statistics: F_1=W^T F_U W and
g_1=W^T g_U. Scaling only the observer 3-by-3 block misses the necessary
calibration-observer cross terms, which get one factor of tau. H_0 uses the
upper nine-coordinate block. Store the complete 12-by-12 Gram matrices.

## 5. What conditional calibration can mean

A discovery-selected fixed direction w can still be used for a legitimate
*conditional* diagnostic. Its H_0 statistic must be

\[
 Z_{0\mid D}=\frac{w^T(y_V-A_Vm_{0D})}
 {\sqrt{w^T(I+A_VV_{0D}A_V^T)w}}.
\]

For fixed design and the exact H_0 Gaussian generative model this is standard
normal conditional on discovery; w can be any fixed function of discovery.
For the covariance-weighted matched filter of a fixed predicted pattern v,
set w=Sigma_{0,V|D}^{-1}v. Its numerator is then
v^T Sigma^{-1}(y_V-mu_0), with denominator sqrt(v^T Sigma^{-1}v).
Thus the selection of w on discovery is not itself fatal. The omissions in
the earlier nominal statistic are the conditioned null mean, common-mode
uncertainty, and the broader modeling/selection limitations. A w chosen using
validation does not have this calibration. Its direction and sign must be
declared before any new attribution result is examined, while retaining the
fact that the overall attribution question arose after the earlier validation.

For Delta, if a model-conditional tail is useful, hold y_D fixed and draw
a*~N(m_0D,V_0D), epsilon_V*~N(0,I), then
y_V*=A_V a*+epsilon_V*. Recompute the two frozen predictive densities for each
draw, without re-estimating either discovery posterior. Use exactly the same
joint mode across every validation object/field. Simulation can use only the
12-dimensional score vector: g_1V*=X_1V^T A_V a*+eta, with
eta~N(0,F_1V), because the common y_V^T y_V term cancels in Delta. Verify this
shortcut against full small-dimensional simulations, and factor semidefinite
F by a documented rank-aware eigendecomposition if necessary.

This calibrates the declared finite-mode, fixed-design Gaussian predictive
model. It is not a uniform frequentist test over every unknown calibration
coefficient, every calibration prior, all SALT uncertainty, or the adaptive
research history. A composite frequentist claim would require a separate
guarantee or an explicit nuisance/selection/adaptation calibration. Do not use
sqrt(2 Delta) as a sigma value: this penalized predictive ratio does not have a
generic chi-square law. Report the signed score and scope even if a simulated
tail is small.

The ten original fields are not independent under a global latent mode.
Flipping their residual signs independently generally changes the common-mode
covariance and does not preserve this null. Conditioning on discovery leaves
V_0D nonzero, so subtracting the posterior mean alone does not restore field
independence. Use the joint Gaussian forecast or simulations from it; do not
carry forward the 1/1,024 field-sign result as a shared-calibration p-value.

## 6. Identification: span coverage and admissible amplitude are different

An observer vector almost contained in the calibration span is geometrically
degenerate with those calibration predictions. It is not necessarily plausible
under their **scaled** prior covariance, and it need not transfer with the same
coefficients to a different epoch/redshift design. Report both geometric span
overlap and prior-scaled singular values, then the actual conditional prediction.

The newly saved [shared43 summary](../../runs/research_2026_09_26/astra_design/shared43/result.json)
reports 95.4334% span coverage of the discovery observer prediction but only
55.4146% retained information after its stated covariance inflation, with
largest scaled calibration singular value 0.96525. These are inspected outputs,
not an independent numerical reproduction here. They support substantial
geometric overlap; they do not establish that calibration explains 95% of the
observed pattern. The displayed discovery matched statistic under the added
covariance remains descriptive because its direction was estimated from those
same discovery residuals.

If an observer column is exactly a calibration combination in the combined
design [D;V], the observations identify only their sum. Adding the extra term
changes the prior variance along an existing direction; a predictive preference
can reflect that variance choice without identifying a new physical mechanism.
With near-collinearity, report singular values/condition numbers, canonical
angles of the relevant template spaces, and posterior cross correlations. Use
rank-aware bases: unqualified reduced QR on a rank-deficient matrix can invent
spurious canonical directions. An observer posterior interval excluding zero
does not eliminate calibration degeneracy or prior dependence.

## 7. Nine released variations are not the whole physical uncertainty budget

The fixed-column construction is a new Gaussian approximation to the second
moments of specific coupled perturbations. The local saved `fitopts.yml` gives
nine cal_1..cal_9 scales 0.3. The historical `get_cov_from_diff` multiplies a
distance difference by its scale **before** taking the outer product. That
supports 0.09 sum_j b_j b_j^T as the explicitly inherited convention for this
epoch experiment; it does not prove an exact prior on a continuous physical
calibration parameter or a full distance-level reproduction.

The methods paper describes ten realizations and its Table 6 uses 1/10, whereas
the archived model set is nominal MODEL000 plus nine perturbed models and the
saved amplitude scale is 0.3. Preserve and report this provenance discrepancy.
Do not silently replace 0.3 by 1/3 or 1/sqrt(10), invent a tenth perturbation,
center the draws, or promote the empirical mean variant to a signed correction.
As a purely mathematical illustration, **if** nine columns were independent
zero-mean draws with covariance C, then E[0.09 sum b_j b_j^T]=0.81 C and the
relative standard deviation of a fixed-direction second moment would be
sqrt(2/9)=47.14%. Those premises are not verified for the archived releases;
these numbers are not a measured 19% covariance deficit or an error bar for DES.
The calibration/training coupling is described in
[Vincenzi et al., section 6.1](https://arxiv.org/abs/2401.02945).

The available 102-by-102 Fragilistic **zero-point** covariance is valuable
external information, but its existence does not supply a signed epoch-level
operator for every zero-point direction, the wavelength-shift covariance or its
cross covariance, or the nonlinear retraining response outside the nine
realizations. A four-band marginal is not a replacement for training's
cross-survey coupling. Do not condition that marginal on the other surveys
without an actual separately observed conditioning data set. Regressing the
nine flux responses on their physical shift lists yields an extension outside
their span only by additional assumptions; a pseudoinverse is not recovered
physical information.

Avoid these specific double counts or substitutions:

- Each native coupled mode already includes its direct calibration leg and
  retrained SALT surface; do not add those again as independent uncertainties.
- The nine-mode covariance and Fragilistic covariance are overlapping accounts
  of calibration, not independent priors to add. Their difference need not be
  positive semidefinite, so it is not automatically a missing covariance term.
- Released distance covariance has different units and a different inference
  layer. Do not add it to flux C or add it alongside explicit replacement modes
  for the same systematic.
- The saved `CALIBplusSALT3` distance grouping also includes the separately
  listed CALSPEC direction according to the existing grouping audit. Adding
  both distance files double counts it. Conversely, this **nine-mode flux model
  excludes CALSPEC** unless a separate, correctly signed epoch response is
  explicitly introduced in both H_0 and H_1.
- The nominal per-object SALT model/error maps remain in C. This exercise does
  not derive a universal decomposition of their empirical scatter and training
  uncertainty; it adds the particular released calibration perturbations under
  the stated convention, with remaining possible modeling overlap explicit.

A finite rank-nine model can assign zero systematic variance to a direction
whose true calibration variance is nonzero. Good or bad prediction under these
nine modes therefore cannot by itself accept or reject all released systematics.
More complete physical response operators and uncertainty distributions are
required before making that claim; selection/BBC and population changes are
additional requirements for corrected distances or cosmology.

## 8. Executed independent checks

Run with:

```sh
.venv/bin/python runs/research_2026_09_26/shared_mode_review/synthetic_checks.py
```

The [script](../../runs/research_2026_09_26/shared_mode_review/synthetic_checks.py)
reads no SN residual data. The [JSON](../../runs/research_2026_09_26/shared_mode_review/synthetic-checks.json)
records its source hash, seed 2026092607, and Python/NumPy/SciPy versions.

| Check | Result | Meaning |
|---|---:|---|
| Low-rank score versus direct posterior predictive versus dense joint-minus-discovery density | max difference 1.066e-14 | Three implementations close on a random synthetic 7/11-row, 3/5-mode design. |
| Same orthogonal prior-coordinate rotation in both samples | score change 3.553e-15 | Common latent-basis choice is immaterial. |
| Rotate validation coordinates alone | score change +0.76465 | Correct within-sample covariance alone does not preserve the cross-sample model. |
| Dimensionless 43/1,020-row coherent toy: all residuals 0.5 with a shared unit-prior calibration column | nominal Z=15.9687; proper conditional Z=0.07380 | Discovery calibration posterior mean is 0.488636 with variance 1/44; the held-out mean is expected under that null. |
| Same toy, reverse the calibration sign only in validation | conditional Z=6.42085 | Wrong cross-sample mode identity manufactures tension. |
| Add an exactly identical unit-prior extra column | predictive gain -0.003445; posterior correlation -0.97727 | The data cannot assign the common signal to its two identical sources. |
| Normalized difference of those identical coordinates | prior variance 1; posterior variance 1 | An exact likelihood-null direction remains unconstrained. |
| Integrate toy validation modes once versus separately for each row | log gains 125.9045 versus 115.9744 | Independent marginal sums evaluate a different latent model. |
| 200,000 draws from the discovery-conditioned toy null | proper 5% exceedance 5.0535%, Monte Carlo SE 0.04873 percentage points | Correct statistic closes to the conditional Gaussian target. |
| Apply the nominal zero-calibration statistic to those same conditional-null draws | 99.770% exceedance at its nominal 5% threshold | The earlier null is not the conditional shared-calibration null. This is a conditional teaching construction, not a DES false-positive estimate. |
| Nine Gaussian synthetic response draws in ten true latent dimensions | empirical covariance rank 9; omitted unit-direction variance about 0 versus true 1 | Finite-span coverage cannot certify the full physical covariance. |

The toy vectors are deliberately constructed; none of their amplitudes,
significances, priors or failure rates are estimates for the observed DES data.

## 9. Acceptance gates for the forthcoming observed-data result

1. **Provenance and paired coordinates.** Hash nominal and all nine input models,
   shift lists, executable and exporter. Verify full fixed fitted coordinates,
   redshift, MW extinction, filters, epochs and duplicate multiplicities. Require
   the intended calibration-shift flags and matched filters, together with
   explained nominal-unit data/error identity. Never drop failed objects or
   epochs silently. Preserve exact CID disjointness: validation excludes all
   64 original discovery IDs, not only the retained 43.
2. **One local likelihood.** Establish baseline objective closure and use the
   same nominal C/J/projection for y, A and T. Retain all four nuisance
   directions with a numerical rank gate and an amplitude/gray null check.
   Explain model-covariance changes as excluded from this local response test.
3. **Prior and mode identity.** Keep all nine signed coupled columns in the
   same order/orientation between samples, scale each once by 0.3, and preserve
   the original three-coordinate 0.02-mag prior. Retain all cross Gram blocks.
   Verify common sign/orthogonal-rotation invariance and failure of a deliberately
   incorrect validation-only rotation on synthetic data.
4. **Training boundary.** Recompute both discovery posteriors from original
   priors and the 43 discovery residuals only. Save their means, covariances,
   cross blocks, input hashes and prediction definition before reading the new
   shared-mode validation score. Do not substitute the earlier observer-only
   posterior or reinterpret this timing as pristine preregistration.
5. **Independent arithmetic.** Compare R differences with direct Gaussian
   conditioning or an orthogonal low-dimensional validation basis. Compare a
   few per-object F/g blocks with direct residual vectors; check symmetry,
   positivity of I+F, total dimensions and complete sums. Integrate each shared
   posterior once across the full validation sample.
6. **Scientific interpretation.** Report signed Delta and both models' actual
   predictive means/uncertainties, geometric overlap, prior sensitivity and
   calibration-observer posterior dependence. Mark any simulated tail as
   finite-mode/fixed-design/retrospective. Retain earlier noise-model sensitivity,
   accepted-mask/classifier conditioning, physical prior limitations and the
   omitted CALSPEC/other systematics. No corrected distance or cosmological
   inference follows from this score alone.

## 10. Executed independent observed-score check

After the new exports and scores became available, the parent requested an
independent arithmetic check. The
[checker](../../runs/research_2026_09_26/shared_mode_review/check_observed_scores.py)
rebuilds every discovery and validation per-object 12-by-12 Gram matrix and
score vector from the stored projected responses and residuals: all agree
**exactly** with the saved blocks. It verifies the 43/1,020 cohort sizes, all
64 discovery-ID exclusions, finite arrays, source/protocol hashes and the prior
transformation above. The saved designs have 1,477 discovery and 35,526
validation residual coordinates.

The independent calculation constructs the discovery posterior from the
reassembled design, then evaluates the conditional validation density in a
rank-aware SVD observation basis with Cholesky solves. This differs from the
producer's joint-minus-discovery R formula and its QR check. The maximum
predictive-score discrepancy is **7.68e-13**; posterior means/covariances agree
within 6.4e-15. The
[verification JSON](../../runs/research_2026_09_26/shared_mode_review/observed-score-check.json)
hashes the inputs and executed checker.

| Discovery-trained predictive model | Validation log-density gain versus fixed zero residual | Difference versus calibration-only |
|---|---:|---:|
| Nine shared calibration modes | 49.8034410823 | 0 |
| Calibration plus inherited observer prior | 58.7470793547 | +8.9436382724 |
| Calibration plus equal-total-variance isotropic griz prior | 58.6779387725 | +8.8744976902 |

Thus the new nine-mode null itself predicts substantial coherent structure.
The observer extension still predicts better in this finite-mode comparison;
this conclusion changes little under the single declared prior-orientation
sensitivity. These log-score gains are not additive percentages of a physical
explanation: dividing 49.8 by 58.7 would not measure a fraction of the residual
explained by calibration. The agreement of these two observer priors does not
establish robustness to all prior scales or to omitted calibration modes.

An independent SVD-based inverse of the calibration predictive covariance also
reproduces the fixed-pattern check: centered inner product 22.5190629229,
information 46.9668752804, Z=3.2859010466 and one-sided Gaussian tail
0.0005082836 **only under this conditional nine-mode Gaussian model**. Adding
the entire original frozen observer vector as a deterministic mean shift gives
log-score change **-0.9643747173** relative to this calibration-conditioned null.
A positive remaining projection and an overlarge fixed amplitude coexist; the
original coefficient vector is not validated as a correction to apply.

The source now explicitly gates nuisance rank four, nominal model/data/error/
epoch/covariance cache identity, all 64 discovery exclusions, and exact
parameter/data/epoch identity for each native variation. This reviewer inspected
those source gates but did not independently rerun the native exports. The
independent score check starts from their saved projected arrays; its strong
arithmetic agreement must not be relabeled as independent raw-photometry or
calibration validation. None of these results supplies a complete released-
systematics significance, physical mechanism identification, selection/BBC
correction, or cosmological correction.

## 11. Expanded twelve-mode review: native transformations and shared scores

The bounded follow-up adds three explicitly released alternatives to the nine
coupled training/calibration columns: CALSPEC at amplitude 1, MWEBV_SCALE=0.95
at amplitude 1, and COLORLAW=CCM89 with R_V=3.0 at amplitude 0.3. These are the
saved fitopt conventions, not twelve independently identified physical causes.
The null has twelve prior coordinates and each observer extension has fifteen.
Discovery conditioning, observer priors, nominal C/J, masks and cohorts are
unchanged. The new
[protocol](../../runs/research_2026_09_26/astra_design/expanded12-protocol.md)
explicitly follows inspection of the nine-mode result.

The CALSPEC source derivation checks against the actual audit executable's
sources, rather than estimating an arbitrary scale from observed residuals:

- `sntools_calib.c`, lines 1627–1667, identifies transmission support above
  1e-6, retains one adjacent edge bin, and computes mean wavelength as
  sum(lambda*transmission)/sum(transmission).
- `snana.F90`, lines 3032–3040, converts that mean to single precision;
  `MAGOBS_SHIFT_USRFUN`, lines 26111–26164, evaluates the specified polynomial
  in wavelength measured in microns. Here delta_m=0.00714*lambda_mean/10000.
- `snana.F90`, lines 17389–17400, multiplies both observed flux and total quoted
  flux error by a=10^(-0.4 delta_m), in single precision. The intended comparison
  is therefore f_variant/a-f_nominal in nominal flux units, with nominal C kept
  fixed. Replacing C by variant C would be a different experiment.

These files are under
[`SNANA-audit-v3/src`](../../phase2/official/build/SNANA-audit-v3/src/).
Independent reconstruction from the hashed original KCOR table gives the
following values and exactly reproduces the saved unit map:

| Band | Transmission-weighted mean wavelength, Angstrom | Native single-precision flux/error scale |
|---|---:|---:|
| g | 4827.6585450 | 0.9968302845954895 |
| r | 6434.8324466 | 0.9957772493362427 |
| i | 7828.0982468 | 0.9948652982711792 |
| z | 9181.2204302 | 0.9939804673194885 |

The independent
[expanded checker](../../runs/research_2026_09_26/shared_mode_review/check_expanded_scores.py)
reads every added native export for all 1,063 objects: **3,189 object/variant
pairs and 123,765 epoch/variant rows**. The source-predicted forward CALSPEC
flux and error mappings match **exactly**; the native CALSPEC model itself is
unchanged before conversion. Inverse transformation cannot undo float32
rounding perfectly, but its maximum data discrepancy is only
2.23161e-5 times the quoted measurement error. The MWEBV float32 multiplication
by 0.95 matches exactly; data and quoted errors remain exactly nominal in
MWEBV and COLORLAW. Epoch, band, fitted-coordinate and redshift identity gates
pass throughout. No new native fits were run by this checker.

All three added columns are independently reconstructed with the fixed nominal
C/J projector, agreeing within 9.26e-17; the saved nominal projected residuals
agree within 1.23e-15. The first nine old mode columns are byte-for-byte
unchanged. All per-object fifteen-coordinate Gram/score blocks rebuild exactly.
Independent rank-aware SVD conditional predictive densities agree with the
producer's R and QR calculations within **9.24e-13**:

| Expanded discovery-trained model | Joint validation log gain over fixed zero | Difference versus expanded null |
|---|---:|---:|
| Twelve-mode systematics null | 50.0708679166 | 0 |
| Null plus inherited observer prior | 58.7521827187 | +8.6813148021 |
| Null plus isotropic observer prior | 58.6815802909 | +8.6107123742 |

The fixed original observer direction has centered statistic Z=3.2339188802
under this partial Gaussian null, while its entire deterministic mean shift
has log-score change -1.2245823811. The
[verification JSON](../../runs/research_2026_09_26/shared_mode_review/expanded-score-check.json)
records native, unit, projection, posterior and score checks with hashes.

There is no CALSPEC distance-group double count here. Independently executing
only the inspected `apply_filter` function from the archived covariance source
confirms that the `+cal` group includes cal_1 through cal_9 **and CALSPEC**.
The expanded experiment starts from the nine individual native modes and adds
one CALSPEC column once; it does not import `CALIBplusSALT3` or any distance
matrix. MWEBV and COLORLAW are separate named columns. The pre-existing
within-object MW variance nevertheless requires a separate physical accounting
qualification below.

## 12. An exact CALSPEC–observer degeneracy

The expanded alternative's fifteen latent coordinates have **likelihood rank
fourteen**, in both observer-prior conventions. This is expected from the
native CALSPEC response: since its unconverted model is unchanged,
delta_f_b=f_b(1/a_b-1). It is a fixed bandwise multiplicative perturbation.
Writing kappa=0.4 ln(10), its equivalent linear magnitude vector is
m_b=-(1/a_b-1)/kappa. Subtracting its griz mean removes the gray amplitude
tangent, leaving exactly a combination of the three observer contrasts.

The zero-sum equivalent vector is approximately
[+0.00160647,+0.00045465,-0.00054481,-0.00151631] mag in griz.
If c_C=(m_g,m_i,m_z), then after projection
A_CALSPEC=T c_C. In the inherited prior-whitened alternative, the latent vector
with CALSPEC entry 1 and observer entries -c_C/0.02 is consequently a likelihood
null direction. After normalizing that vector to prior variance one, the
independent [null check](../../runs/research_2026_09_26/shared_mode_review/calspec-observer-null-check.json)
finds predicted-vector norms 1.78e-14 on discovery and 8.78e-14 on validation.
Its discovery posterior variance remains 0.9999999999999999 and its mean is
consistent with zero to 2.5e-14.

The data cannot separately attribute this combination to CALSPEC and the
empirical observer term. A proper prior still makes the fifteen-dimensional
Gaussian integral well defined; dropping a column without transforming the
prior would change the model. The positive predictive extension score does
not resolve this exact physical-label ambiguity.

## 13. What the two Milky Way five-percent terms represent

The actual nominal local covariance already contains an MW uncertainty term.
In the executed pinned source, `MWgaldust.c` lines **302–307** implement
OPT_MWEBV=3 as E(B-V)=0.86 E(B-V)_SFD98 and set its per-object error to 5% of
that value. `snlc_fit.F90` line **6929** defaults OPT_COVAR_MWXTERR to 1; the
executed log explicitly records that value. Lines **10727–10745** form
ERR_i*ERR_j for epochs of one object and add it to that object's model/data
covariance. `snana.F90` lines **17201–17208** then apply the requested
MWEBV_SCALE to both the map value and its local error. The present sensitivity
retains nominal C, so it does not substitute the changed error covariance.

Assembled into a block-diagonal nominal likelihood, those local terms act like
separate per-SN errors. Adding one common 0.95 normalization response introduces
cross-SN covariance, under an additional independent global coordinate. It is
not a second application of the 0.86 mean correction. Nor does matching the
two amplitudes establish that the random quantities are physically independent.
To first order, if both five-percent responses were the same direction, their
independent addition would give twice that direction's within-object variance
and a new cross-object component; the two fractional errors would combine to
about 7.07%, not stay at 5%. This is a description of the adopted latent model,
not a claim about measured map uncertainty.

The pinned DES methods paper's [local text](../../papers/text/2401.02945v2.txt),
section 6.3, lines **1255–1259**, explicitly describes its systematic experiment
as a global 5% scaling. Reviewed sections 2.7, 3.6 and 6.3 do not supply a
physical variance decomposition that establishes independence from the local
5% term in the executable. This source evidence therefore warrants the label
**release-style added shared normalization sensitivity, with unresolved
physical variance decomposition**. It does not establish that nominal C should
be reduced, so no covariance term was silently subtracted or another fit arm
run. Retaining this qualification is part of distinguishing an already-applied
mean correction from empirical uncertainty and a shared systematic model.

## 14. COLORLAW paper/configuration weight discrepancy and frozen sensitivity

The paper's section 6.3, text lines **1270–1274**, specifies color-law covariance
weight W_S^2=1/3, with the square also appearing in its covariance definition
at lines **651–660**. The saved `fitopts.yml` instead supplies COLORLAW=0.3;
the inspected `get_cov_from_diff`, lines **743–757**, scales the response before
its outer product, producing covariance factor 0.09. The expanded12 result
faithfully uses that saved configuration/source convention. These records alone
do not establish which numerical convention produced the historical release.

At the parent's request, a
[separate protocol](../../runs/research_2026_09_26/shared_mode_review/colorlaw-weight-protocol.md)
was frozen before one explicit saved-array sensitivity. It changes only the
COLORLAW prior-scaled column from amplitude 0.3 to sqrt(1/3), multiplying
column 11 by 1.9245008973 in both cohorts; its covariance contribution therefore
increases by 3.7037037. The shared discovery posteriors are recomputed under
both declared observer priors. All native exports, original coefficients,
other mode weights, data, projection, masks, and historical results are retained.

| COLORLAW amplitude | Null log gain | Inherited-observer log gain | Observer-minus-null | Isotropic-observer-minus-null |
|---|---:|---:|---:|---:|
| Saved configuration: 0.3 | 50.070868 | 58.752183 | 8.681315 | 8.610712 |
| Paper covariance convention: sqrt(1/3) | 50.462646 | 58.763303 | 8.300657 | 8.225968 |

The inherited preference decreases by 0.380658 log-score units; the isotropic
preference decreases by 0.384744. The original direction's conditional Z changes
from 3.233919 to 3.172842, and its full deterministic shift score from -1.224582
to -1.572909. Independent SVD and R calculations agree within 5.55e-13.
This particular documented weighting inconsistency modestly weakens, but does
not remove, the finite-model predictive preference. It neither resolves the
historical execution provenance nor establishes robustness to every missing
mode or uncertainty prior. The executed
[script](../../runs/research_2026_09_26/shared_mode_review/colorlaw_weight_sensitivity.py),
[results](../../runs/research_2026_09_26/shared_mode_review/colorlaw-weight-sensitivity.json)
and saved posterior are separate from the original expanded12 artifacts.

## 15. Local distance-response sign, target and prior pushforward

The root's [distance-response script](../../scripts/research_2026_09_26/shared_distance_response.py)
uses delta_f=variant mean minus nominal mean. Holding observations fixed, its
first-order compensation has the correct sign:

\[
 \delta p=-(L^{-1}J)^+L^{-1}\delta f,\qquad
 \delta\mu_{\rm ref}=(1,\alpha,-\beta,0)\delta p,
\]

with p=(m_B,x_1,c,t_0), alpha=0.16087 and beta=3.1178. Since J_mB=-kappa f,
an extra +1-mag gray dimming of the model has compensation delta_mB=-1; its
projected residual response is zero, as the saved gray check requires. This
is a reference-coordinate, local pre-BBC compensation. Retrained SALT
coordinates, refitted alpha/beta, nonlinear covariance/parameter transport,
selection and BBC prevent interpreting it directly as a measured distance
correction.

For the saved high255-minus-low255 validation target, the contrast response
row h must be formed first and propagated as h m_D with variance h V_D h^T.
The same global coordinates affect all SNe, so treating object uncertainties
independently or dividing their common-mode uncertainty by sqrt(N) would be
wrong. The residual likelihood can constrain the part perpendicular to J while
leaving a substantial absorbed distance response uncertain. An unrestricted
redshift-dependent gray luminosity term remains outside this finite-mode bound.

The final bounded check freezes the exact existing 255/255 membership and
uses the root's completed twelve-mode response arrays. It independently
reconstructs that contrast exactly and transforms either its row or the full
1,020-object response matrix, agreeing within 6.60e-17. For the isotropic prior,
the observer response row is h_obs/0.02 times chol(C_c); for the paper-weight
COLORLAW sensitivity, h_11 is multiplied by sqrt(1/3)/0.3. The matching saved
discovery posterior is used in each case, with no new fitting or bin choices.

| COLORLAW convention | Discovery posterior family | Conditional contrast mean, mag | Conditional mode SD, mag |
|---|---|---:|---:|
| Configuration amplitude 0.3 | Twelve-mode null | +0.01134173 | 0.01279334 |
| Configuration amplitude 0.3 | Plus inherited observer | +0.05831766 | 0.02571699 |
| Configuration amplitude 0.3 | Plus isotropic observer | +0.06291010 | 0.02604077 |
| Paper amplitude sqrt(1/3) | Twelve-mode null | +0.01153640 | 0.01287203 |
| Paper amplitude sqrt(1/3) | Plus inherited observer | +0.05831759 | 0.02571709 |
| Paper amplitude sqrt(1/3) | Plus isotropic observer | +0.06291099 | 0.02604089 |

These are **conditional finite-mode contributions to the chosen compensation
contrast**, not its total statistical errors or measured biases. The roughly
0.00459-mag mean change under the observer-prior orientation sensitivity,
despite very similar flux scores, illustrates that reference/prior dependence.
The [frozen protocol](../../runs/research_2026_09_26/shared_mode_review/distance-prior-sensitivity-protocol.md),
[independent arithmetic](../../runs/research_2026_09_26/shared_mode_review/supplemental_checks.py)
and [result](../../runs/research_2026_09_26/shared_mode_review/distance-prior-sensitivity.json)
record this boundary and input hashes. No significance is assigned to these
contrast means.

## Source and review boundary

Reviewed the requested calibration investigation, shared-uncertainty audit,
scientific design, validation verification and frozen validation protocol;
checked the archived fitopt scaling code and original observer-template source;
and inspected the newly written shared43 result during this review. The
existing [shared-uncertainty audit](shared-calibration-uncertainty.md) contains
the detailed local release/source links and full-budget gaps. The original
1,020-object score remains an independently verified **nominal-model** transfer
result as documented in [validation verification](validation1020-verification.md).
The new shared-mode scores and conditional statistic were independently
recomputed as recorded above; no full-systematics significance is assigned.
