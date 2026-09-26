# DES colour/shape cut sensitivity in the frozen 1,020-object cohort

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../README.md) for current execution status.

On the exact previously selected DES cohort, the retrospective strict `|x1|<3`, `|c|<0.3` check retains **1,019/1,020** using the independent native nominal coordinates, **1,016/1,020** under the frozen full nonlinear observer model, and **1,018/1,020** under the frozen broadSED model. Three native-pass objects newly fail the observer check; one of those also fails the SED check. No object enters either arm from a native failure. These are coordinate-based partial eligibility changes, not changes to the released DES sample or distances.

The [protocol](sed-quality-cut-protocol.md) was hashed **before computing cut-change outcomes** (`5a809ad6cdc03cb0275a2b19bea7e2aa1584d1386bab95dabcb9e0ef017b2617`; [freeze](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/sed_quality_cut_probe/freeze.json) SHA-256 `9b0608b9f5626a9775e4c23b47a0072ecd0f8fe5fb61729be156c7bbaf9f411b`). The published DES-SN5YR analysis (`papers/text/2401.02945v2.txt`; historical local reference), section 5.1/Table 4, specifies strict `-3<x1<3`, `-0.3<c<0.3`, `x1ERR<1`, `PKMJDERR<2`, and `FITPROB>0.001`, together with earlier light-curve, redshift and host requirements and later BBC/systematics membership. The release BBC input (`sources/repos/des-science__DES-SN5YR@1.3/7_PIPPIN_FILES/base_files/bbc/BBC_des5yr.input`; historical local reference) additionally lists `CUTWIN cERR 0 1.5`, redshift range and `chi2max=16`; the NOMINAL Pippin configuration (`sources/repos/des-science__DES-SN5YR@1.3/7_PIPPIN_FILES/D5yr_analysis.yml`; historical local reference) points at this base and applies `chi2max=16`. The exact historical BBC run is not replayed here.

## Same objects, cuts and rows

All 1,020 CIDs join uniquely among the frozen cohort, released HD metadata (`sources/repos/des-science__DES-SN5YR@1.3/4_DISTANCES_COVMAT/DES-SN5YR_HD+MetaData.csv`; historical local reference), original validation FITRES (`runs/research_2026_09_26/astra_design/validation1020/fit.FITRES.TEXT`; historical local reference), objective NPZ arrays and resolved [nonlinear responses](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/sed_nonlinear_validation/resolved/paired-responses.csv). The resolved observer and SED rows all have passing observed-target fit gates. The native baseline uses the unrounded objective `[x0,x1,c,t0]` with zero shift. Each alternative uses `native_x1+alternative_theta[1]`, `native_c+alternative_theta[2]`, as defined by the frozen nonlinear Engine. The observed-target nominal fit is reported only as a diagnostic, not substituted for the exact native baseline. The original 39,606 accepted epochs, photometry, errors and classifier outcomes are untouched. `x1ERR`, `PKMJDERR`, `cERR`, `FITPROB`, redshifts, detection, host criteria, BBC and common-systematics membership remain fixed at their released decisions.

| Inputs to the same strict cut | Shape pass | Colour pass | Joint pass |
|---|---:|---:|---:|
| Released HD metadata | 1,020 | 1,020 | 1,020 |
| Exact native nominal coordinates | 1,019 | 1,020 | 1,019 |
| Observer observed-flux nonlinear fit | 1,018 | 1,018 | 1,016 |
| BroadSED observed-flux nonlinear fit | 1,018 | 1,020 | 1,018 |
| Nominal observed-flux refit, diagnostic | 1,018 | 1,020 | 1,018 |

All released values also pass the fixed numerical `x1ERR`, `PKMJDERR`, `cERR`, `FITPROB`, and `zHD` proxy windows in the protocol. No x1 or c coordinate is exactly on a cut boundary. The objective-vs-FITRES mapping closes at the expected FITRES float32 half-ULP; the maximum absolute difference is 0.001951 day in t0. The resolved theta-shift identities and observer/SED nominal theta equality close exactly. For 24 objects FITRES lists compound field labels such as `C2+C3`; the frozen cohort's single canonical field is a member of each corresponding compound label and is used in the tables below.

## Transitions

| Comparison | Pass→pass | Pass→fail | Fail→pass | Fail→fail |
|---|---:|---:|---:|---:|
| Native → observer | 1,016 | 3 | 0 | 1 |
| Native → broadSED | 1,018 | 1 | 0 | 1 |
| Observer → broadSED | 1,016 | 0 | 2 | 2 |

The complete [per-object table](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/sed_quality_cut_probe/per-object.csv) retains every CID, released/native/arm coordinates, transition, frozen zHEL/zHD, and canonical field. The only CIDs in a native-failure or native-to-alternative transition are:

| CID | Field | zHEL | Released x1,c | Native x1,c | Observer x1,c | BroadSED x1,c | Interpretation |
|---|---|---:|---|---|---|---|---|
| 1307748 | C3 | .52589 | −.8105, −.0729 | +4.5237, −.0280 | +4.6789, −.0422 | +4.6712, −.0267 | Native and both alternatives fail shape. This is the previously documented native/released fit-branch discrepancy, not an observed loss from either perturbation. |
| 1936353 | X2 | .69290 | −2.6587, −.1078 | −2.6548, −.1079 | −3.2351, −.1385 | −3.2391, −.0873 | Both alternatives fail shape. The **unperturbed observed-target nominal refit** also reaches x1=−3.2432, so this shared crossing cannot be attributed specifically to either perturbation. |
| 1390343 | X3 | .89270 | −.5035, −.2659 | −.5035, −.2648 | −.5030, −.3228 | −.5229, −.1854 | Observer only fails colour. |
| 1559048 | X3 | .98157 | −.8734, −.2595 | −.8733, −.2585 | −.8675, −.3128 | −.8948, −.1486 | Observer only fails colour. |

The frozen field and redshift distributions show precisely where these four objects lie. Counts are joint-pass counts, with cohort size in `n`.

| Field | n | Native | Observer | BroadSED |
|---|---:|---:|---:|---:|
| C1 | 108 | 108 | 108 | 108 |
| C2 | 110 | 110 | 110 | 110 |
| C3 | 134 | 133 | 133 | 133 |
| E1 | 91 | 91 | 91 | 91 |
| E2 | 94 | 94 | 94 | 94 |
| S1 | 55 | 55 | 55 | 55 |
| S2 | 96 | 96 | 96 | 96 |
| X1 | 92 | 92 | 92 | 92 |
| X2 | 95 | 95 | 94 | 94 |
| X3 | 145 | 145 | 143 | 145 |

| zHEL bin | n | Native | Observer | BroadSED |
|---|---:|---:|---:|---:|
| [0.025,0.2) | 31 | 31 | 31 | 31 |
| [0.2,0.4) | 339 | 339 | 339 | 339 |
| [0.4,0.6) | 382 | 381 | 381 | 381 |
| [0.6,0.8) | 230 | 230 | 229 | 229 |
| [0.8,1.2] | 38 | 38 | 36 | 38 |

An [independent readback](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/sed_quality_cut_probe/independent-verification.json) recomputed the native and alternative coordinates from the source objective arrays and resolved theta strings, then reconstructed strict-cut counts and both baseline transition matrices without importing the probe script. It agrees exactly for all 1,020 coordinates and counts. The [result](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/sed_quality_cut_probe/result.json) and [manifest](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/sed_quality_cut_probe/manifest.json) preserve source/input hashes and the machine-readable field/bin details.

This selected cohort cannot reveal objects excluded by the original pipeline. The alternative fit errors, classifier, detection efficiency, bias-correction validity, common-systematics intersection and cosmological standardization were not recomputed. The counts therefore measure only this partial coordinate-cut sensitivity under the frozen hypothetical means; they do not establish a changed DES selection, corrected Hubble diagram or cosmology.


## Common-cut influence follow-up

A separately frozen post-audit diagnostic removes the union of all four flagged
objects from both alternatives and the original residual score. It retains the
original high/low quartile identities, renormalizing their means to 252/255
remaining objects rather than selecting new quartiles. The primary 1,020-object
analysis is unchanged. This is descriptive influence analysis after the cut
outcomes, not a new significance test or a replay of survey selection.

The original fixed observer residual gain changes **27.2332→27.5511 nats**;
its descriptive amplitude changes **.58244→.58348**. On the observed-flux
nonlinear fits, the observer high-minus-low response changes
**+.060585→+.059612 mag**, and the SED response **−.186125→−.184683 mag**.
Their separation is therefore **.244295 mag**, versus the original .246710.
The four flagged objects do not drive either the transfer or the constructed
distance-response ambiguity in this fixed-cohort calculation.

The [protocol](../specifications/experiments/sed_cut_influence/protocol.json),
SHA256 `f2bf58864250b82662ef22576f6ee4b7cc0491ff2189c49b21fb97d4c50c8e6c`,
predates this influence score and hashes all inputs and the executed code.
[Results](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/sed_cut_influence/result.json),
[removed-object scores](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/sed_cut_influence/removed-object-scores.csv)
and [script](../code/sed_cut_influence.py) are retained.
The untouched original totals and nonlinear quartile contrasts are reproduced
in the same run before computing the reduced-cohort values.
