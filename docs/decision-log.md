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
