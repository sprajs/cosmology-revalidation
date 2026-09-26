# Independent expansion-shape checks

These are new analyses of the pinned DESI DR2 **compressed** BAO likelihood.
They do not remeasure galaxy clustering. The initial tests assume flat FLRW metric geometry,
positive expansion, one constant comoving BAO ruler and the released Gaussian
effective-redshift likelihood. They use no supernova luminosity calibration,
stellar-age correction, CPL dark-energy family, calibrated H0 or calibrated rd.
The curvature extension below removes flatness while retaining additional
branch conditions. The Astra investigator independently checked the principal inequalities and
null calibration. All runs preserve executed source and input hashes.

## Registered primary result

Writing t=ln(1+z), d=DM/rd and g=(1+z)DH/rd gives d'=g and g'=-qg.
Nonacceleration, q>=0, requires positive nonincreasing g. The six anisotropic
distance pairs obey a convex cone of integral and monotonicity inequalities.
The fit's chi-square distance to that cone is **11.4671**. A globally conservative
composite-null calibration gives **p=.4166** (10,000 simulations; binomial 95%
interval .4069–.4263). This test does **not** reject an always nonaccelerating
history; it also has limited power. Its power at the calibrated 5% threshold is
40.62% for the chosen flat LCDM Omega_m=.3 alternative, and only .33% at .4.
Failure to reject cannot establish nonacceleration.

Coasting-only and fitted-null simulation tails are .3057 and .1092 respectively;
neither is a global composite-null probability. An independent constrained
optimizer agrees in chi-square within 5.6e-9. Common scale changes cancel.

Source/results: [bao_shape.py](../../scripts/research_2026_09_26/bao_shape.py),
[calibrated results](../../runs/research_2026_09_26/bao_shape/calibrated/results.json).

## Exploratory BGS inequality

After inspecting the primary result, a new necessary inequality brought in the
previously held-out isotropic BGS measurement at z_b=.295. This is explicitly
post-pilot exploration, with its [amendment](bao-bgs-amendment.md) saved before
its numerical evaluation. Define

    s_b = (DV_b/rd) / [z_b/(1+z_b) * ln(1+z_b)^2]^(1/3).

For q>=0, the mean g below z_b and its endpoint g_b are both at least any later
g_j; therefore s_b=(mean(g)_b^2*g_b)^(1/3)>=g_j. This gives linear Gaussian
contrasts of released observables. Coasting makes all contrasts zero and is
least favourable for their maximum negative standardized contrast.

BGS gives s_b=32.04368±.30702; the radial point at z=.934 gives g=34.11865.
The contrast is **4.1883 marginal standard deviations**. Accounting for the
explicit family of 28 BGS, AP, radial and integral contrasts gives **p=.000341**
(341/1,000,000 draws; 95% interval .000306–.000379). This family correction does
not erase the adaptive decision to introduce the statistic after the pilot.
It is an interesting independent indication of **acceleration somewhere in
0<z<.934**, conditional on the assumptions above; the inequality does not
localize it solely between the BGS and radial endpoints. It is not a
new-physics discrepancy, a measurement of q(0), or a refutation of every recent
deceleration model.

Closing the central-value inequality would require BGS DV to increase by
6.475%, or DH(.934) to decrease by 6.082%. An overall ruler or H0 calibration
cannot change the contrast. A relative evolving ruler would need
rd(.934)/rd(.295)=.93918 at the coasting boundary. Allowing arbitrary unknown
cross-correlation while retaining the two marginal errors weakens the same
contrast to at least 2.982 standard deviations; a six-comparison Bonferroni
upper bound is .00859. This does not allow unreported mean biases or extra
variance. A necessary closed-curvature bound requires |Omega_k|>=8.291 to close
the central values under no acceleration; this is an analytic bound, not a
curvature posterior or a full curved-universe fit.

Sources/results: [BGS test](../../scripts/research_2026_09_26/bao_bgs.py),
[results](../../runs/research_2026_09_26/bao_shape/bgs_exploratory/results.json),
[assumption budget](../../runs/research_2026_09_26/bao_shape/assumption_budget/results.json).

## Removing the flatness assumption

A further exploratory extension allows arbitrary open/flat/closed FLRW
curvature, retaining positive expansion, one constant comoving ruler and the
positive first-antipode branch in closed geometry. The high-redshift transverse
measurement bounds a closed geometry's curvature radius: `R/rd>=DM(2.33)/rd`.
The BGS/radial measurements also establish that BGS must precede transverse
distance turnover if the history never accelerated. These statements yield
the necessary inequality in the [frozen protocol](bao-curvature-protocol.md).

At the central measurements, `DM(2.33)/rd=38.98897` limits curvature's reduction
of the z=.934 BGS bound to a factor **.994320**. Even at that limit,
nonacceleration requires `DV(.295)/rd>=8.40791`, versus the observed **7.94168**.
Curvature alone therefore cannot close that central-value gap within this
geometry/branch family. This is stronger than quoting an enormous required
Omega_k without checking the other distance measurements.

A conservative confidence test uses one-sided bounds for BGS DV, highest-z
DM and all six radial quantities. At k=2.708665 marginal standard deviations,
the limiting values are BGS<=8.14778, radial g(.934)>=33.06547 and
high-z DM>=37.54883. The inferred BGS transverse upper bound is only 8.47414,
so its pre-turnover gate is safely satisfied. The union bound for all eight
one-sided events is **8 Phi(-k)=.02702**, valid without specifying their
cross-correlation if the quoted Gaussian marginal errors are correct. This
bound accounts for the six tested radial choices, not the broader adaptive
research history. It remains an exploratory indication of past acceleration,
not evidence for present acceleration, LambdaCDM specifically, or new physics.

Five thousand random nonaccelerating histories satisfy 63,704 valid geometric
certificates, including 1,897 cases whose high-z transverse distance is after
turnover. An independent scalar bisection reproduces thresholds within
2.4e-13 and independently verifies the eight-event probability argument.
[Executed result](../../runs/research_2026_09_26/bao_shape/curvature_certificate/result.json),
[script](../../scripts/research_2026_09_26/bao_curvature_certificate.py), and
[Astra review](../../runs/research_2026_09_26/astra_design/curvature-certificate-independent.json)
preserve the derivation, bounds and checks. Ruler evolution, unmodelled mean
bias/variance, anisotropy and multiple-antipode geometry are not covered.

## Constructive recent-deceleration counterexample

Continuous-H, piecewise-q fits to all 13 BAO observables give chi-square 2.4064
with free recent q, and **2.5887 with q=+.5 throughout z<=.3**, allowing earlier
acceleration. A coarser binning gives 3.6311 and 3.8998 respectively. Thus these
compressed observations do not establish the sign of present acceleration.
These flexible constructive fits are not calibrated model comparisons or
posterior probabilities. Bounds, knots and every optimizer start are retained
in the [results](../../runs/research_2026_09_26/bao_shape/recent_counterexample/results.json).

An exploratory reconciliation with 1,590 Pantheon+ rows profiles a common
intercept and retains their full covariance. For the fixed BAO history with
q=+.5 through .3, the required standardized-luminosity shape spans .3241 mag.
A fixed age-template shape absorbs part of that difference with a conditional
amplitude corresponding to .0752 mag/Gyr, while the SN chi-square changes from
1858.1 to 1429.8. This is **not** an inferred age law: the BAO history is not
unique, stellar ages are not fitted, selection is not regenerated, and modern
SN corrections already contain correlated host information. These fixed
examples are not minimum correction budgets over all acceptable histories.
See [reconciliation](../../runs/research_2026_09_26/reconciliation/summary.json).

## Boundaries still requiring work

BAO extraction, reconstruction, tracer modelling, redshift kernels and
fiducial-distance compression remain upstream assumptions. Sharp cone histories
could invalidate effective-redshift compression, although a smooth constrained
fit reaches chi-square 11.535, close to the cone optimum. The signal is not
fully independent of observational calibration simply because it is independent
of supernova luminosity. An unrestricted gray SN luminosity evolution remains
exactly degenerate with SN distance and requires external information to break.

That degeneracy exists at the flux and selection level, not merely in a Hubble
diagram. Schematically, observed flux is proportional to
`L_lambda(t,z) * exp[-tau_lambda(z)] / d_L(z)^2`, with known frame factors.
Replacing distance by d'_L and luminosity by `L' = L*(d'_L/d_L)^2` leaves every
flux, colour and light-curve shape unchanged if the extra factor is gray.
A selection rule acting on those same measurements is also unchanged. A fully
reproduced BBC pipeline therefore cannot by itself rule out arbitrary gray
luminosity evolution. Counts help only with independently constrained rates
and volume/selection assumptions. Optical–NIR colours, spectra, host variables
and phase information constrain departures from gray evolution; absolute or
geometric information is needed for its unrestricted component.
