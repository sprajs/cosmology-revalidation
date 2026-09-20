# Experimental preparation, without running the experiment

> Historical acquisition-stage record. Current analysis, resolved gaps and remaining barriers are in [the final report](final-report.md) and [current literature/data map](literature-and-data-status.md). The original scientific content below is preserved.

The current stage is acquisition, provenance and study design. No regression, residual recomputation, cosmological likelihood evaluation, chain reweighting or scientific hypothesis test has been run. Existing author chains are evidence artifacts. Counts, file hashes, archive CRCs, FITS headers and array shapes have been inspected to establish what is actually available.

## Separate the possible reproduction targets

| Target | Observable to reproduce later | Inputs already available | Principal missing specification |
|---|---|---|---|
| A. Published W26 matched-sample result | Age–HR relation before and after corrections | Revised ages; exact Pantheon+ snapshot; published figure and method | Author crosswalk, duplicate/cut rules and regression/error settings |
| B. Published Chung response | Four Figure 1 variants; cumulative-cut and mock Figure 2 | Published article, same broad public sample components | Exact matched table, uncertainty subtraction, mocks and cut configuration |
| C. Son cosmology | Corrected SN distances, combined constraints and `q0` | SN releases, BAO vector/covariance, specified CMB likelihood data and DESI reference configurations | Son's correction array/covariance, exact release pins, settings and corrected chains |
| D. Host/progenitor mapping | Redshift evolution and slope transformation under each DTD/SFH/selection | C14 paper, W22 code/library, S25/W26/C26 equations | Author run configurations and selection/posterior products |
| E. Dust-model comparison | Galaxy attenuation versus SN extinction under specified selections | GSWLC catalogue, dust papers, DES metadata and calibration/simulation products | Exact dust posterior sample and forward-model selection mapping |
| F. Dovekie calibration changes | Published old-versus-new distance/covariance changes | Both DES releases, SALT/calibration products and chains | Verified runnable historical pipelines; unavailable current LFS mocks for full simulation reproduction |
| G. Sah directional analysis | Isotropic/dipole cosmography with specified age correction | Public scripts, input tables and correction arrays | Environment audit, defaults and interpretation of incomplete upstream output files |

These targets answer different questions. Select one explicitly before interpreting a difference as confirmation or refutation of another paper.

## Proposed sequence for a later analysis stage

1. Freeze a target's exact input files and original paper version. Preserve original tables separately from normalized working copies. Resolve high-priority access gaps where exact replication depends on them.
2. Build a documented ID crosswalk. Retain survey IDs, aliases, repeated light-curve measurements, common G11/R19 objects, redshift definitions and cut masks. A table row is not automatically an independent SN.
3. Establish a correction ledger for every distance: light-curve model/version, stretch, colour, host term, simulation bias, peculiar velocity, intrinsic scatter and any age correction. Retain signs, units and zero points explicitly.
4. Reproduce a small published descriptive result on an identical sample before varying assumptions. Define how asymmetric/non-Gaussian age errors, outliers, redshift and host selection enter the regression.
5. Separate sample/cut changes from model changes. Compare age versus mass adjustments and scatter choices one at a time, with the same data whenever possible. The published C26 four-panel sequence changes multiple ingredients overall.
6. For a host-to-progenitor comparison, specify both the transformation of the age distribution and the transformation of the inferred luminosity slope. Propagate their dependence rather than assuming they are independent multipliers.
7. Only then reproduce unmodified likelihood baselines and introduce explicit correction variants. Do not combine overlapping SN compilations or CMB information as independent evidence. Record priors, nuisance marginalization, likelihood versions, seeds, convergence and numerical settings.
8. Distinguish an exact published reproduction, a controlled variant and a new analysis in separate outputs. Publish no conclusion until the data and method actually support it.

No future decision here is pre-committed to a particular value of `q0`, a preferred dust law, or either side's interpretation. The current source collection makes those choices inspectable; it does not decide them.

## Practical conventions already established

- Inputs live in `data/` and pinned `sources/repos/`; acquired archives remain intact. Future generated scientific outputs belong in a separate dated `runs/` tree, which has not been created.
- File dimensions/headers are in `catalog/table_inventory.json` and `catalog/npz_inventory.json`. These checks do not establish correct covariance modelling, convergence, or absence of selection bias.
- DES-Dovekie uses a precision matrix and a whitespace HD file despite its suffix. Preserve its row order and explicitly configure the loader. Existing release code/default names require inspection before use; no runtime compatibility is assumed.
- For uncalibrated SN distances, absolute magnitude and `H0` are degenerate. The DES likelihood marginalizes the offset. Importing SH0ES calibrators changes the experiment and must be explicit.
- At low redshift retain the distinction between heliocentric, CMB-frame and peculiar-velocity-corrected redshift. For BAO retain the ruler convention `r_d` and which observables are `D_M/r_d`, `D_H/r_d` or `D_V/r_d`.
- The data support a future study, but host ages and SN progenitor ages are inferred quantities. “Public data acquired” is not equivalent to “physical interpretation verified.”
