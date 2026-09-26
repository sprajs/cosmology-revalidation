# A fixed-mask SNooPy profile gate before bias simulations

The current native multi-start result does not establish a unique fitted
distance. Across 240 fits of 80 cases, all native return flags are zero and
reported covariances positive, yet 20 cases move DLMAG by more than .01 mag;
the maximum range is .224064 mag with the same accepted rows. This is a
numerical-readiness failure, not evidence for that size of physical distance
bias. Matching rows is insufficient: native covariance and peak-prior state
also depend on previous iterations. `FITCHI2` excludes priors and is not the
complete minimized function.

The source-level check below executes the unchanged native rest-grid
interpolation functions on predetermined coordinates. It finds sharp slope
changes at stretch knots, with tiny value discontinuities from float32 grid
arithmetic. A profile calculation must preserve those kinks and all relevant
branches. It must also distinguish a broad or multimodal objective from a
failed optimizer; replacing every such case by the smallest local reported
error would be incorrect.

No observed profile or native fit was launched in this review. Astra science
owns the native exporter and the first one-object proof on DES16C1cim,
author-signed B, fixed header peak, reference start 1. This note supplies
independent source checks, an executable amplitude-profiling kernel and the
solver acceptance criteria.

## The native objective is stateful

Sources are the pinned `phase2/official/build/SNANA-current/src` files.
`snlc_fit.F90:3942–3945` explicitly says the minimized `CHI2TOT` includes
`CHI2INI`, whereas `FCN_FITCHI2` does not. `FCNCHI2_PRIOR` at line 4847
centers the peak term on `INIVAL(IPAR_PEAKMJD)`, **not necessarily the original
header**. `FITINI_PARVAL`, lines 7919–7950, sets each later iteration's initial
values to the previous fitted values plus `FITERR*INIVAL_SHIFT_ERRFRAC`
(the default multiplier is zero), then returns before first-iteration
initialization. The peak prior thus recenters. Three iterations are an
algorithm with state, not automatically minimization of one fixed-header
Gaussian-prior objective. The native step sizes also inherit previous fitted
errors; a tiny knot-related error can restrict the next iteration's search.

The nominal log confirms `USE_MODEL_MAGERR=T`, `OPT_CHI2_SIGMA=0`,
`OPT_COVAR_MWXTERR=1`, `OPT_COVAR_FLUX=0`. There are two regimes:

* In iteration 1, the diagonal residual denominator is
  `sigma_data^2 + [f*(1-10^(-.4*magerr_model))]^2 + sigma_fudge^2`.
  The model flux and magnitude error are evaluated at each candidate
  (`snlc_fit.F90:4285–4344`), so this denominator can vary with parameters.
  With `OPT_CHI2_SIGMA=0` there is no log-variance contribution. This native
  estimating objective must not be relabelled a normalized Gaussian
  likelihood or silently replaced by one.
* In iterations 2/3, `FITINI_COV` enables full covariance through Galactic
  extinction even when `OPT_COVAR_FLUX=0` (around lines 10490–10500).
  The assembled matrix includes previous-state diagonal model-flux errors,
  the Galactic-error outer product, and data/fudge errors
  (10651–10730). `LOAD_EPALL` carries forward model-flux quantities and
  stores their preceding values in `R4EP_LAST`. `FCNSNLC` recomputes this
  matrix **inside** minimization only for `OPT_COVAR_FLUX=2`
  (4563–4580). Nominal option zero therefore uses a frozen matrix within
  that iteration, but different runs/iterations can freeze different matrices.

For one fixed later-iteration state S the comparison is well defined:

\[
 Q_S(p)=[y-f(p)]^T C_S^{-1}[y-f(p)]+P_S(p).
\]

Freeze row identities/order, data, covariance, prior center/table/bounds,
model/KCOR/MW settings and any cache/iteration-dependent model behavior.
An accepted-mask match alone does not imply the same S. Even full native
FCN values from different S cannot rank a single likelihood. A profile of
one S is a conditional diagnostic; establishing a stable iterative solution
later requires an explicit outer-state update and convergence definition.
Holding a common header-centered prior instead would be a deliberately
different objective, not a repair claimed to reproduce the original algorithm.

The existing audit-v3 exporter can provide `PHASE2_FLUX`, `COVINV`,
`OBJECTIVE`, `SETUP` and `COVSET`. `SETUP` already exports the final iteration's
initial values and previous fitted values. Its legacy parameter slot names
must be relabelled: for SNooPy, X0 aliases **DLMAG**, and X1 aliases stretch.
No SALT logarithmic conversion belongs in this parser. The available export
is final-iteration only; reconstructing the full three-iteration algorithm
requires separate per-iteration state exports. Original binaries and source
remain unchanged; any added export belongs in an isolated audit build with
nominal-output closure.

The native-engine proof also found that audit-v3's per-object list reader
accepts SALT column names x0/x1/c, while later SNooPy initialization requests
DLMAG/STRETCH/AV. Its attempted all-fixed list run aborted and is being
preserved. The one-object proof uses explicit native NML parameter values
instead, with their seed precision recorded. A fast batch oracle must not
reuse the SALT list reader by renaming physical coordinates; it needs a
separately audited generic-column extension or direct native function calls.
Steps smaller than input-coordinate precision cannot diagnose continuity.

## Executed rest-grid continuity check

The native `gridinterp_snoopy` in `genmag_snoopy.c` interpolates **magnitude
and magnitude error**, despite local variable names `gridFlux*`. It uses
bilinear interpolation in phase and stretch. `INDEX_GRIDGEN` in
`sntools_modelgrid_read.c:515` computes cell boundaries from the global
minimum plus a uniform step, while interpolation subtracts the stored
float32 knot. `fits_read_SNGRID` recomputes the step from float32 endpoints
and bin count. These boundaries need not coincide exactly with every stored
float32 knot. Packed magnitudes/errors have .001-mag table resolution.

[grid_continuity.py](../../runs/research_2026_09_26/raisin_profile_solver_review/grid_continuity.py)
extracts the two C function bodies unchanged, uses the original grid-structure
header and loads the released FITS arrays with the native one-based and
float32 conventions. The small adapter changes no interpolation arithmetic.
The C library is compiled locally in the owned run directory; no native
installation is modified. A separately written Python four-corner formula
provides an independent numerical comparison.

Before evaluating outcomes, the script freezes all 129 internal algorithm
stretch boundaries, all nine rest filters BVugriYJH, phases -7,0,10,30,45
days and one-sided offsets of .001 and .0001 cell. It extrapolates each
linear side to the boundary to separate a limiting jump from an ordinary
slope difference. The [protocol](../../runs/research_2026_09_26/raisin_profile_solver_review/grid-protocol.json)
hash is `927af208919bd8e916913557dc1498d135a82eb87b3693090d75702b7d96878c`.

The [result](../../runs/research_2026_09_26/raisin_profile_solver_review/grid-result.json)
is:

| Check | Result |
|---|---:|
| C versus independent bilinear formula | max 3.55e-15 mag |
| Native stretch spacing | .00461538415402174 |
| Stored knot versus algorithm boundary | max 1.24e-7 stretch |
| Limiting magnitude jump | max 5.06e-7 mag |
| Limiting magnitude-error jump | max 4.34e-6 mag |
| Change in boundary jump after reducing step tenfold | max 4.71e-14 |
| One-sided magnitude-slope difference | max 1.9500 mag per stretch |

The largest sampled slope switch is rest g at phase 30 and stretch
1.1892307084: left +.6500001, right -1.3000001 mag/stretch. All sampled
rest magnitudes remain inside the native rest-magnitude validity limits;
none is clamped. The scan took .27 seconds. It is a rest-grid component
check, **not** full observer-flux/KCOR continuity or objective closure. It
does not prove that half-micro-mag jumps cause the observed .224-mag
distance range.

Knot minima can be legitimate minima of this piecewise model. If the left
objective slope is negative and the right positive, a cusp is a true local
minimum even though no ordinary Hessian exists. A centered finite-difference
curvature at such a cusp can grow like 1/step and give an error shrinking
like sqrt(step). That is a reason to inspect profile widths and separated
branches, not to smooth the source grid without changing the model or to
declare every knot solution an optimizer bug. Observer K corrections,
colour warping, phase knots and native flux validity limits introduce
additional structure that the native oracle must retain.

## The bounded profile fallback

The first proof must hold the accepted mask, covariance and peak fixed at the
frozen reference S. Use the complete native observer-flux calculation for
every candidate; a monochromatic filter approximation, Python SNooPy variant,
or newly smoothed template would answer a different question.

For fixed stretch s, AV and peak t, source `USRFUN` adds DLMAG to the
observer magnitude before flux conversion. If all participating fluxes stay
inside its validity domain, native mean flux must satisfy
`f(D,s,AV,t)=a*h(s,AV,t)`, with `a=10^(-.4*(D-Dref))` and h the native
prediction at Dref. Establish this numerically at Dref±.5 and ±2 mag, using
the unchanged mask. Repeated point calls in shuffled order must close too,
so stale model caches cannot masquerade as derivatives. Export validity
flags: native `USRFUN` returns zero outside observer magnitude (5,40), and
crossing that switch can invalidate simple amplitude scaling.

With frozen W=C^-1 and an inactive distance prior,

\[
 b=h^TWy,\quad q=h^TWh>0,\quad \hat a=b/q,\quad
 \hat D=D_{ref}-2.5\log_{10}\hat a.
\]

Evaluate the residual quadratic directly at the solution to avoid
cancellation in `y^TWy-b^2/q`. The
[implemented kernel](../../runs/research_2026_09_26/raisin_profile_solver_review/fixed_covariance_profile.py)
requires finite inputs, matching dimensions and a positive-definite symmetric
inverse covariance. It refuses nonpositive amplitudes and candidate distances
outside the verified flat-prior interval (default [10,60]); those require a
bounded one-dimensional **native-prior** distance profile. The source native
hard DLMAG interval is [10,70], with additional guards above 60. The kernel
does not silently clip or pretend an active prior is flat. Any other active
prior is supplied from the same frozen native state.

An algebra-only test recovers amplitude 1.3 to 2.22e-16, with residual
quadratic 3.50e-29, and gives unchanged distance after simultaneous row
permutation or per-band flux-unit transformation of y,h,C. It also verifies
the two refusal paths. These vectors are a kernel unit check, not a
replacement SN model or native-physics validation.

After these gates, profile the remaining s/AV surface as follows:

1. Restrict stretch to the **actual released table support**, approximately
   [.7,1.3]. Native broad prior/hard bounds permit extrapolation, so this is
   a declared supported-model diagnostic. Include stored knots, algorithm
   boundaries and both one-sided limits; the tiny float32 discrepancies
   must not be hidden by choosing only one set. Include every cell midpoint.
2. Begin with the complete source AV hard interval [-2,6], a .1 grid and
   all active native-prior breakpoints. Evaluate every shape cell, not only
   cells near existing fitted solutions. At fixed shape bracket and refine
   **every** candidate AV minimum and retain endpoints; an AV profile need
   not be unimodal. Use derivative-free bounded refinement within cells,
   with all local branch results saved. Analytic distance profiling reduces
   this to two dimensions where its gates pass. Native clipping or an active
   distance prior invokes the bounded distance oracle instead.
3. Refine each stretch interval continuously, retaining all knot candidates
   and interior branches. Repeat with AV step .05, quarter-cell stretch
   seeds, reverse traversal and independent nuisance starts. A coarse global
   scan plus local minimization is not a mathematical global proof; require
   stability of the lower envelope and all competitive branches under this
   refinement. Any unresolved boundary/clipping region that could improve
   the minimum prevents a global-support claim.
4. Save profile minima and distance ranges at **declared** objective offsets
   Delta Q=1,4,9, plus all disconnected competitive branches. These are
   descriptive profile supports, not automatically calibrated confidence
   intervals or posterior probabilities. Do not use a single covariance
   error at a grid cusp to summarize them. To obtain a distance profile,
   minimize s/AV at fixed D or otherwise construct the joint profiled
   surface; the set of amplitudes minimizing Q for each shape alone is
   not a full distance-likelihood curve.

The proof gate is one object, one fixed mask/state/peak. Use a batch native
mean oracle inside one initialized process if possible. Launching a fresh
native fit for every surface point is unnecessary and can silently refresh
the state being held fixed. Time the first 100 native vector evaluations;
cap the first proof at 100,000 evaluations and 20 wall minutes on at most
two CPUs. Save the incomplete grid/checkpoint if the cap is reached, with
no convergence claim. A rough 261-shape by 81-AV initial grid already needs
about 21,000 evaluations before refinements, so timing the oracle is a real
feasibility gate. The native engine agent owns this execution and the exact
grid protocol; this review has not started it.

## Gates for a defensible result

At the exported reference coordinate, require reconstructed
`residual^T CINV residual + native prior` to reproduce the full native
FCN to 1e-7 absolute (or a separately justified printed-precision bound).
Require audit/current nominal output equality and exact row-identity gates.
Record covariance symmetry/positive definiteness without a hidden PSD repair.
Native multiplicative-distance and repeated-oracle checks must agree to
1e-8 relative on nonzero supported flux, with a 1e-7 quoted-error absolute
gate for vanishing flux. A tiny finite difference alone is insufficient;
compare explicit saved objective values.

For the profile refinement require minimum Q stability within 1e-4 and
interior minimizing DLMAG stability within 1e-4 mag. If a flat region or
multiple branches prevents a unique D, report its stable support instead
of forcing parameter agreement. Compare the native fitted point to the
profile using the **same** S; a lower Q then demonstrates incomplete
minimization of that frozen objective. If native multi-start points cease
to disagree materially after all are evaluated under one S, iteration-state
dependence is a more direct explanation. Both mechanisms can coexist.

Only after the fixed-peak proof should time become a continuous nuisance.
Hold the same peak-prior center and covariance throughout that diagnostic;
account for the phase-grid breakpoints at
`t0=MJD_j-(1+z)*phase_k`. Keep the fixed accepted mask and intersect the
time bounds with valid phase support for every row. A second diagnostic
can compare separate frozen states or iterate state updates, but each must
retain its own identity. Adding a parameter-dependent Gaussian covariance
likelihood, a log determinant, a fixed-header timing prior or a smoothed
template is a separate model/algorithm experiment, not an invisible solver fix.

The global/profile gate must pass before expensive sign-selection bias
simulations. Otherwise the apparent simulated bias mixes the epoch-selection
effect with optimizer or iterative-state branch selection. Even a repaired
conditional solver does not establish a physical extinction law, historical
SNANA equivalence, or a cosmology correction; the modern MW99 model and the
measurement/selection limits in the preceding design remain explicit.

## Completed independent native proof and first global-profile review

The subsequent saved-artifact checks are durable under
`runs/research_2026_09_26/raisin_profile_solver_review/native-proof-review/`
and `global-profile-review/`, each with its own report, executable checker
and hashes. No native fits were launched by these independent reviews.
The native proof passes original/instrumented output identity, exact accepted
74-row membership, covariance/quadratic closure, float32-coordinate distance
multiplicativity and order replay. The final pre-score global-profile
amendment corrects the native float32 cell-boundary calculation and adds
explicit native support and competitive-amplitude checks.

The completed first-object global numerical gate also passes. All 86,813
saved means/four metrics and all 93,159 raw native stream calls are checked;
maximum reconstructed Q disagreement is 9.095e−13 and DLMAG disagreement is
7.106e−15 mag. Both covariance anchors, the 65-positive/74-signed membership,
all coarse/fine scans, and all 1,548 boundary probes close. The finest
minima are within 1.012e−6 in Q and 1.666e−7 mag of the coarse minima.

This does **not** identify a unique distance. The primary frozen-B covariance
gives minimum-to-minimum B−A=−0.164949213 mag (alternate anchor
−0.165638773 mag), but the signed-data profile has ten sampled local branches
within ΔQ=0.551. Its high-stretch branch has DLMAG=42.767926739, only
ΔQ=0.115738988 above the global DLMAG=42.588805188: a 0.179121551 mag
distance separation. Adding the nine negative rows changes the preferred
branch through a small net balance of the positive and conditional-negative
quadratic terms. This is a conditional processing sensitivity, not an
identified correction.

The primary B saved-coordinate ΔQ≤1 distance envelope, including analytic
amplitude freedom, spans 42.471390–42.865636. It is not a confidence interval
or posterior; the ΔQ9 envelopes also reach the declared shape boundary.
All branch and boundary details are in `global-profile-review/README.md`.
The numerical result supports extending the unchanged diagnostic to the
predeclared cohort while retaining modes and failures. It does not yet
support physical bias/cosmology correction or interpreting A/B minimum
quadratics from different row sets as a likelihood ratio.

## Completed independent ten-object extension

The final all-ten numerical extension is independently verified under
`runs/research_2026_09_26/raisin_profile_solver_review/cohort10-review/`.
All560 closed artifacts and14 executed sources hash correctly; all670,988
saved mean vectors/four metrics reproduce (maxQ discrepancy5.457e−12,
maxDLMAG discrepancy7.106e−15 mag). All20,958 actual native amplitude
checks are independently reconstructed from replies (maxrelative4.633e−9).
Exact source-row multiplicity, reference/nominal/prefix/final identity,
marginal covariance restrictions, all grids, minima, branch geometry and
edge gates close. All ten pass the frozen numerical gates.

The two zero-change controls have exactly identical A/B metric arrays and
exactly zero global-minimum distance differences. A small nonzero offset
in the separate nearest-fine-grid-branch diagnostic is explicitly a
discretization artifact against a refined minimum. It must not be called
a processing effect. The largest conditional minimum shifts are−164.949,
−49.376 and−19.045 mmag, with competing branch and broad distance support
retained. These do not establish population bias or corrected distances.
At ΔQ4 DES16E1dcx reaches the shape boundary; at ΔQ9 five objects do.
The full report and `final-certification.json` retain those restrictions
and the exact closure of provisional runner-log prefix hashes.
