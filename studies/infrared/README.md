# Infrared distances, timing and signed flux

The released optical-minus-infrared high-minus-low-redshift contrast is +0.00124 mag; reversing exported mass and bias terms makes it +0.07541 mag on the same 79 objects. The published infrared systematic covariance is not fully reconstructed from the available raw shift vectors. All 29,995 matched infrared simulated peak times equal their initializer. Native-fit stability failures and missing selection information prevent claiming a new infrared distance correction.

## Code and scope

Correction ledgers, object pairing, source acquisition, historical native fits, timing profiles and population simulations, signed-flux likelihoods and censoring controls.

The [source directory](code/) retains the original calculations. See the [execution guide](../README.md#execution) before running them: historical imports and relative paths are assembled into an isolated workspace, and some studies require additional data or native software. Source preservation does not certify that an incomplete experiment succeeded.

## Supporting findings

- [Signed-optical BayeSN engineering gate](notes/bayesn-signed-optical-engineering.md)
- [Signed optical BayeSN fits with held-out NIR: a feasible conditional pilot](notes/bayesn-signed-optical-pilot.md)
- [Coherent 2021 RAISIN simulated-fit cut accounting](notes/raisin-2021-selection-accounting.md)
- [Coherent 2021 RAISIN simulation timing metadata](notes/raisin-2021-timing-population.md)
- [Independent RAISIN cadence and DES error-bridge review](notes/raisin-cadence-error-bridge-review.md)
- [Does the RAISIN covariance mismatch affect observable distance shape?](notes/raisin-covariance-quotient.md)
- [RAISIN WFC3/IR count-rate nonlinearity: inclusion audit](notes/raisin-crnl-inclusion.md)
- [RAISIN: paired optical–NIR distance audit](notes/raisin-differential.md)
- [RAISIN photometry flux-sign integrity audit](notes/raisin-flux-sign-audit.md)
- [Independent review of the historical-source RAISIN refits](notes/raisin-historical-native-review.md)
- [Source-version-matched RAISIN refits](notes/raisin-historical-native.md)
- [RAISIN fitter source-version review (read-only)](notes/raisin-historical-source-review.md)
- [Native recovery and archived-simulation timing: next computation](notes/raisin-native-recovery-timing-design.md)
- [Independent native RAISIN signed-refit stability review](notes/raisin-native-stability-review.md)
- [RAISIN NIR timing inheritance: source and ten-object audit](notes/raisin-nir-inheritance-review.md)
- [Ten DES16 NIR distances under fixed header-date changes](notes/raisin-nir-timing-sensitivity.md)
- [A same-object NIR phase contrast is useful, but the released cadence fails the frozen feasibility gate](notes/raisin-phase-identification.md)
- [Why the first signed-light-curve distance shift is not an identified correction](notes/raisin-profile-decomposition.md)
- [A fixed-mask SNooPy profile gate before bias simulations](notes/raisin-profile-solver-review.md)
- [Prospective native timing recovery: engineering record](notes/raisin-prospective-native-engineering.md)
- [Native noise and support review for the prospective timing experiment](notes/raisin-prospective-noise-review.md)
- [RAISIN signed-flux selection: a paired generative simulation design](notes/raisin-selection-simulation-design.md)
- [Signed pre-explosion RAISIN ancestor measurements](notes/raisin-signed-baseline.md)
- [Signed-epoch response in all ten DES16 objects](notes/raisin-signed-cohort.md)
- [Paired RAISIN optical SNooPy signed-epoch test](notes/raisin-signed-refit-design.md)
- [Source assets behind the RAISIN simulation corrections](notes/raisin-simulation-assets.md)
- [RAISIN DES simulation fit execution provenance](notes/raisin-simulation-execution-review.md)
- [Archived DES NIR simulation LCPLOT input feasibility](notes/raisin-simulation-lcplot-feasibility.md)
- [Independent archived RAISIN simulation coherence audit](notes/raisin-simulation-replay-coherence.md)
- [Are omitted negative epochs represented in RAISIN bias simulations?](notes/raisin-simulation-sign-accounting.md)
- [Archived DES simulation FITRES timing audit](notes/raisin-simulation-timing-audit.md)
- [Approximate archived DES NIR timing pilot: baseline gate failed](notes/raisin-simulation-timing-pilot.md)
- [RAISIN timing assets: source-selected historical recovery](notes/raisin-timing-asset-recovery.md)
- [A distinct quantization-aware timing experiment: design only](notes/raisin-timing-quantization-design.md)
