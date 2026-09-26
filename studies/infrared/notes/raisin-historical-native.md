# Source-version-matched RAISIN refits

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../README.md) for current execution status.

Rebuilding the source version named in the archived RAISIN optical fit table,
and retaining its original Milky Way extinction-law default, reproduces all ten
DES16 nominal distances within **0.00019 mag**. Seven agree exactly at the
table's printed distance precision. This substantially narrows the historical
reproduction gap. It does not establish that the original compiled executable,
libraries or complete execution environment have been reproduced.

The same source branch still fails the starting-point stability test:
**25 of 80 configurations vary in fitted distance by more than .01 mag**,
35 by more than .001 mag, with a maximum spread of **.39153 mag**. All 240
primary fits report `ERRFLAG_FIT=0`, and all 80 comparisons retain identical
accepted measurement sets across starts. Thus this particular instability is
not explained by the modern code version or by changed epoch membership.
The largest spread is in the fixed-peak sensitivity. In the released-data,
floating-peak branch that most closely matches the archived procedure, two of
ten objects exceed .01 mag and the maximum spread is **.02172 mag**. Neither
number is a measured supernova luminosity bias.

## Source, calibration and execution

The pinned author `FITOPT000` explicitly names `v11_04k`. The official
[SNANA tag](https://github.com/RickKessler/SNANA/tree/v11_04k) resolves to
commit `0583dd7ed547c56093c4e1191c4a1c9af205389e`, annotated tag object
`175c778aef1f12f5141c1bc66ef449de8b7074f5`. All 212 source files are hashed in
the [acquisition record](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_historical_native/acquisition.json).
The clean source is retained separately from the local build. No scientific
source edits or global installations were made. The build uses the existing
local GCC/gfortran, GSL and CFITSIO environment; legacy compiler warnings remain
in the build log.

The original RAISIN NML comments out `OPT_MWCOLORLAW=99`. The historical source
default is **94 (O'Donnell)**, whereas the earlier modern comparison explicitly
used **99 (Fitzpatrick)**. The primary historical branch therefore leaves that
setting unspecified, matching the original configuration and source default.
A separate old-source/law-99 control keeps the modern law for a controlled
version comparison. Neither law is declared correct by proximity to an
archived fit. The [source review](raisin-historical-source-review.md) traces
the grid, covariance, prior and mask mechanisms in both versions.

The [protocol and amendment](../specifications/experiments/raisin_historical_native/fit-protocol-shortprefix.json)
were frozen before accepted fit outcomes; final protocol SHA256 is
`1f7cdf1f16feafee068e36fc2633db75659515afffe3b88046ae54e83daf2a24`.
The initial metadata collection failed because the compiler lived outside the
inherited PATH; using its absolute path repaired metadata collection only.
The first two engineering attempts then exposed the historical C output-prefix
buffer's 100-character limit. Their long absolute prefixes produced a malformed
FITRES filename. Those runtime outputs are not accepted for inference and are
preserved with the original protocol, NMLs and logs. Repeating in fresh
directories with the relative output prefix `fit` avoids that limit without
changing data, source or scientific settings.

After this infrastructure correction, the first object's original and
byte-identical copied input produce identical printed FITRES fields and
byte-identical LCPLOT data. The
[engineering gate](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_historical_native/engineering-gate.json)
records 65 accepted epochs, fixed RV=1.518, and distance 42.77922 versus the
archived 42.77917. The full run then completes 240 primary fits: ten objects,
four frozen data arms, free/fixed peak, and initial stretch 1.0/.85/1.15.
Thirty additional old-source/law-99 free-peak release fits are retained as
controls. All 270 outputs are kept.

## What reproduces, and what does not

For the nominal old-source/default-94 released-data branch, every object has
the same NDOF as its archived fit. Distances agree to printed precision for
DES16C3cmy, E1dcx, E2clk, E2cqq, S1bno, X3cry and X3zd. Differences for
C1cim, S1agd and S2afz are +.00005, +.00010 and +.00019 mag. Their peak dates
and some reported errors also differ slightly; exact historical reproduction
is not claimed.

| Numerical comparison | Modern controlled source/law 99 | Historical source/default law 94 |
|---|---:|---:|
| Primary fits across three initial stretches | 240 | 240 |
| Configurations compared | 80 | 80 |
| Same accepted mask across starts | 76 | 80 |
| Distance spread >.001 mag | 28 | 35 |
| Distance spread >.01 mag | 20 | 25 |
| Largest distance spread | .22406 mag | .39153 mag |

The largest historical spread occurs for DES16E2cqq, fixed peak, released
photometry. Other large spreads occur in its author-positive/signed branches,
and in the restored signed free-peak DES16C1cim branch. These are algorithm
sensitivity measurements on the same input, not estimates of omitted-flux bias.
The entire cohort and all failed stability gates remain in the output.

Native `FITCHI2` excludes prior terms. More fundamentally, later iterations
construct covariance from the previous fit and recenter the peak penalty on
the preceding iteration. Two solutions with identical epoch masks may therefore
still have different covariance and penalty states. No preferred solution,
likelihood ratio, confidence interval or new correction is obtained by choosing
the smallest printed `FITCHI2`.

The summary also mechanically combines squared reported parameter errors with
the tabulated off-diagonal `COV_*` entries. Some such reconstructed matrices
are not positive definite. This alone does not establish that the native
Hessian covariance is invalid: reported one-parameter errors and Hessian
covariance need not be the same uncertainty representation. Their source mapping
has now been checked independently: the historical source reports averaged
positive/negative one-parameter errors, while the off-diagonal entries come
from MINUIT's separate `MNEMAT` matrix. It explicitly warns that the matrix
diagonal need not equal those squared errors. Thus the 105/270 nonpositive
reconstructions are a warning against combining incompatible uncertainty
representations, not a finding of 105 invalid native Hessians. No matrix is
clamped or repaired for inference.

The [modern independent stability review](raisin-native-stability-review.md),
[fixed-objective solver review](raisin-profile-solver-review.md), and
[selection-simulation design](raisin-selection-simulation-design.md) specify
the next gate: use common accepted rows and a stated, fixed covariance/prior
state, map competing template-shape solutions, and verify profile resolution.
Only a numerically supported signed-selection comparison can inform the
subsequent regenerated bias correction. Optical timing passed into NIR and
survey-level selection remain downstream requirements.

Runnable sources: [build](../code/build_raisin_historical.sh),
[frozen native execution](../code/raisin_historical_refit.py),
and [summary](../code/raisin_historical_summary.py).
The [result bundle](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_historical_native/summary.json)
links every fit's raw outputs and hashes; the
[baseline comparison](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_historical_native/author-and-version-baseline.csv)
and [start-spread table](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_historical_native/start-spreads.csv)
retain all objects and controls. The independent raw FITRES/LCPLOT and hash
review reproduces all 270 outcomes, 80 mask comparisons and the archived-distance
comparison; the [completed independent review](raisin-historical-native-review.md)
preserves the covariance qualification and all 346 verified hash entries.
