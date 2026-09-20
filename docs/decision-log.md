# Decision log

## 2026-09-20: scientific phase opened

The user explicitly authorized substantive analysis in the goal objective, superseding the acquisition-only restriction. Originals remain unchanged. Acquisition dossier committed as `472f13f`; scientific branch is `codex/independent-supernova-investigation`. Three independent Astra Ultra investigators assigned age-signal reproduction, standardization/dust/calibration, and host/progenitor mapping/causal structure. Initial reports must be recorded before sharing interpretations. Agent agreement is not data.

Pinned Python 3.12 environment and `uv.lock` created locally; no global configuration changed. Independent root work covers cosmology, kinematics and reproducibility. Published source versions, including the August Chung reply, take precedence over the earlier dossier. Missing exact author arrays/configurations do not authorize silent substitution.

## 2026-09-20: cosmology numerical validation before posterior fits

The released Pantheon+ covariance has maximum antisymmetry 3e-8 mag² from printed entries. Symmetrized as (C+Cᵀ)/2 for Cholesky operations; direct inversion of the unsymmetrized original changes the baseline chi-square by only 2e-9. No added diagonal errors. Compact Cobaya inputs match both full Pantheon snapshots byte-for-byte.

The initial 16-node-per-interval interpolation failed the predeclared 1e-5-mag maximum-error criterion on extreme piecewise-q draws (maximum 0.000139 mag), although chi-square error was below .01. Retained this failed check in `runs/cosmology/validation-initial-16/`. Increased to 40 nodes per interval and reran the same checks before inference; this is numerical refinement, not outcome-selected cosmological model adjustment.

## 2026-09-20: explicit correction sensitivity templates

Before corrected cosmology fitting, choose B13 cosmic-SFH equation reconstruction with C14 smooth DTD, comparing median delays under LCDM/H0=70 and Son's BAO+CMB CPL/H0=63.6; also compare mean delays and the W26 short-delay power law. The mapping investigator supplies curves covering z=0–2.5 and independent clock/grid checks. These are fixed-cosmology templates, not Son's self-consistently recomputed author correction. Mean vs median tests concern the statistic transported into a mean magnitude likelihood. Shared slope error uses one amplitude N(1,(4/30)^2); a uniform amplitude on [−2,4] is a broad identifiability sensitivity, not an empirically justified evolution prior.

The literature update found a September revision of arXiv:2607.24443 and a public 114-row low-z host-mass correction reconstruction. When independently checked, apply this as a separately named sensitivity with unchanged covariance, not a new official Pantheon+ release. It is distinct from adding a cosmic age-evolution curve.

## 2026-09-20: final inference diagnostics registered

Before these new outcomes, register 200 Gaussian simulations with the fixed released Pantheon+ covariance and true flat-LCDM Ωm=.33, seed 61020, fitting the same offset-marginal likelihood and measuring bias, empirical scatter and approximate curvature-interval coverage. This tests conditional numerical recovery, not whether the released covariance or astrophysics is true. Perform posterior predictive chi-square checks for the original CPL and fixed/shared-age-template models: simulate the Gaussian distance-product discrepancy after offset removal (N−1 degrees of freedom), compare against observed discrepancy across posterior draws, and report lower-than-expected scatter as well as excess scatter. Add replicate-ensemble block MC errors and freeze exact code snapshots referenced by manifests. No parameterization or correction will be selected because it passes these checks better.

## 2026-09-20: provenance repair before final release

During concurrent authoring, the first cosmology manifest implementation hashed source files at completion, so a file could change during an already-running fit. The scientific equations were unchanged after validation; later changes added revision-table support and reporting. Nevertheless those end-of-run hashes do not reliably identify the executed version. Capture exact code bytes/revision at process import, archive them by SHA-256, and rerun every retained fit from its saved configuration under frozen source before finalizing. Preserve the first-run summaries separately and check numerical agreement; do not retroactively relabel their manifests as exact execution provenance.

## Final freeze and independent review, 20 September 2026

All 19 cosmology configurations were replayed with exact code bytes captured at import. Scientific posterior fields and sign fractions are identical to the initial seeded runs; only wall time and an added zero-revision metadata field differ. Validation, template generation, reference-chain summaries, finite-bin likelihood limits, Gaussian diagnostics and figures were regenerated afterward. The independent audit was rerun against these outputs and retained byte-identical numerical results.

The public C2 covariance was recovered and its principal heliocentric likelihood reproduced: q_m=.009517→.353991 under the age shift, independently checked at four final points to 9.1e−12 in −2logL. No full multiframe/global-optimum or calibrated dipole-significance certification is inferred. Ray v2 and the Sah coordinate response were added to the primary-source baseline; the original Ray full-sample reversal cannot be explained by the coordinate error alone.

Final synthesis review requested two clarifications, both adopted: arbitrary SN luminosity evolution limits SN-only inference, not all independent geometric probes; the supplied W22 library is rate-weighted simulation output, not the exact recovered W26 survey selection. The reviewed report bytes are preserved in `runs/audit/provenance/`; final edits are these clarifications, the completed C2 numbers and record links. No conclusions were selected by agent agreement.
