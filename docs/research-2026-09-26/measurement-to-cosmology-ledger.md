# Where calibration and bias corrections enter the inference

The present priority is to identify supported changes to the measurement and
population model before refitting cosmology. Published corrected distances
already contain several coupled operations. Removing one visible term does
not undo those operations; adding a second correction for the same physical
effect can count it twice.

The earlier [product correction ledger](../correction-ledger.md) retains exact
DES, Pantheon+, Dovekie and host-age conventions. This continuation adds the
observation-level tests and the RAISIN optical/NIR branch.

| Stage | Already included or assumed | What has been checked here | Remaining requirement |
|---|---|---|---|
| Extracted transient flux | Survey image/reference subtraction, exposure zero point, quoted error, epoch quality and retention | Native DES flux/error identity and covariance exports; repeated-row controls. RAISIN DES signed ancestors expose positive-flux-only product retention, including negative rows within the fit phase | Exact upstream retention/quality code and reference-image covariance. Signed pre-explosion measurements constrain a baseline regime, not all near-peak errors |
| HST detector response and pixel errors | Accumulated-count nonlinearity, dark and flat corrections, ramp slope/error construction, followed by downstream aperture/reference subtraction; separate count-rate nonlinearity may require a manual correction | Current files identify CALWF3 3.7.3 and exact calibration references. Native-pixel masked repeats have .708/.729 of the summed quoted diagonal variance. Supplied flat-error terms are too small to explain that difference; read-covariance accounting remains open. Released RAISIN code begins from a missing upstream aperture catalogue, so historical count-rate correction inclusion is unresolved | Establish count-rate correction and standard-star pivot upstream, measure temporal/spatial ramp covariance, and reproduce SN/host subtraction. Current FLT error diagnostics cannot be transferred automatically to historical empirical aperture errors or simulations |
| Photometric reference and passband | Survey standard stars, reference spectra, zero points, transmission curves and colour transformations | Source-defined calibration responses, CALSPEC unit mapping, stellar branch tests, foreground-map response; all 23 BayeSN pilot filter conventions bridged explicitly. CSP WIRC J rows traced unchanged into a label selecting RC1, conflicting with the paper's RC2 grouping; independent spectra verify shape-dependent instrument differences. Stable all-42 RC1-to-RC2 label comparison gives mean +0.483432 mmag, with ten exact controls | Correct physical WIRC operator and residual uncertainty; execution-linked calibration/training inputs and independent standards. The conditional label response is not a bound on training, extinction or population effects. Changing a reference convention is not automatically a measured calibration error |
| Foreground extinction | Map E(B−V), scale convention and wavelength law inside the flux model; associated uncertainty can be local or shared | Source sign/scale accounting and released MW variations. Historical RAISIN default law94 differs from modern explicit law99; source-matched nominal distances reproduce closely | Validate spatial/map/law uncertainties and their dependence before adding local and shared terms as independent |
| Empirical light-curve/SED model | SALT training, colour law and standardization; or SNooPy/BayeSN spectral and intrinsic-population assumptions | 1,020 nonlinear refits expose nearly matching flux responses with substantially different distance responses; training overlap audited; BayeSN distance-amplitude algebra and native forward gates verified | Independent wavelength prediction, justified intrinsic-SED priors, stable numerical inference and retraining where required. Good flux fit alone does not identify distance or physical dust |
| Event and epoch selection | Detection, host redshift, classifier, fit quality; simulated bias corrections condition on these processes | Modern classifier implementation reproduced but released probabilities remain unreproduced; matched simulation controls performed. RAISIN SIMLIB exactly inherits surviving positive epochs | Regenerate the same observation-dependent retention on every simulated realization. Fixed surviving times plus fresh signed noise do not reproduce positive-flux selection |
| Distance and host corrections | DES/Pantheon stretch/colour terms, mass dependence inside bias corrections, explicit residual step; RAISIN exported bias and mass terms | Exact signed arithmetic closes for released products. RAISIN paired optical–NIR contrast changes from +.00124 to +.07541mag when exported bias/mass terms are removed, at fixed fits | Infer which correction the data support. The change on removing existing terms is not the amount of true bias |
| Systematic uncertainty | Shared calibration/training variations, extinction, populations and covariance weighting | Twelve shared modes explain substantial residual structure but leave an additional predictive preference. RAISIN mass-threshold variant has a source-linked propagation defect; total covariance reconstruction fails even after removing the intercept | Reproduce the covariance recipe and propagate repaired response directions. Do not replace unresolved covariance or mix marginal one-parameter errors with Hessian off-diagonals |
| Cosmological inference | Standardized luminosity versus redshift, distance duality/propagation and geometry; external rulers or distance priors | BAO-only kinematic tests separate past from present acceleration; outcome-free SN–BAO designs assess drift sensitivity. Exact grey luminosity–distance degeneracy is retained | Admit correction models only after measurement, selection and recovery gates. Preserve shared objects/systematics when combining probes |

The **positive-epoch finding concerns RAISIN's released DES DIFFIMG photometry**,
used here to test optical/NIR consistency. The main DES-SN5YR cosmology uses
SMP photometry, a different extraction. The finding is not silently transferred
to that dataset. Shared physical SNe and calibration/training dependencies also
prevent treating every release or model as an independent experiment.

The immediate experiments follow this causal order:

1. Establish the signed measurement/error lineage and actual epoch selection.
2. Recover stable light-curve inference on those measurements, with a clearly
   specified covariance, epoch set and objective.
3. Compare the observed conditional response with forward simulations that
   apply the same measurement and event selection.
4. Test optical predictions against held-out NIR, accounting for inherited
   peak timing, training overlap and calibration.
5. Propagate supported population/calibration corrections and uncertainty to
   distance and then cosmology.

Current evidence links: [calibration investigation](calibration-investigation.md),
[CSP passband and magnitude lineage](csp-passband-lineage.md),
[independent NIR spectra](csp-spectral-identification.md),
[CSP native numerical and support gates](csp-native-baseline-review.md),
[complete CSP filter response](csp-stable-filter-response.md),
[generation and existing bias-correction accounting](csp-filter-bias-propagation.md),
[coherent simulation cuts and weighting revisions](raisin-2021-selection-accounting.md),
[physical WIRC reference feasibility](csp-wirc-operator-feasibility.md),
[Vega reference bridge](csp-wirc-operator-vega-bridge.md),
[CSP peak-header lineage](csp-peak-header-provenance.md),
[signed-flux lineage](raisin-flux-sign-audit.md),
[error-conversion and cadence review](raisin-cadence-error-bridge-review.md),
[native error-map applicability and image-to-noise construction](raisin-prospective-noise-review.md),
[HST detector count-rate correction inclusion](raisin-crnl-inclusion.md),
[public HST image and calibration-version lineage](raisin-hst-pixel-feasibility.md),
[completed current-image signed-aperture test](raisin-hst-signed-pixel-noise.md),
[ramp and shared-flat variance accounting](hst-ramp-variance-accounting.md),
[exact science/error component replay](raisin-calwf3-native-replay.md),
[RAW read-prefix metadata](raisin-dark-raw-pilot.md),
[dark-exposure quality and exact gain references](raisin-dark-quality-reference-audit.md),
[baseline noise](raisin-signed-baseline.md),
[historical native refits](raisin-historical-native.md),
[stable-profile identification](raisin-profile-decomposition.md),
[NIR timing inheritance](raisin-nir-inheritance-review.md),
[coherent historical timing recovery](raisin-timing-asset-recovery.md),
[same-photon timing and quantization analysis](raisin-timing-quantization-design.md),
[prospective native recovery and preserved failures](raisin-prospective-native-engineering.md),
[selection simulation](raisin-selection-simulation-design.md),
[BayeSN identification](bayesn-distance-identification.md),
[RAISIN corrections](raisin-differential.md), and
[all E00–E17 directions](programme-status.md).

No finite collection of null tests can confirm that every possible bias is zero.
Each null result must state the alternatives and effect sizes it had power to
detect. Likewise, an implementation or selection discrepancy must be quantified
before it can support a claim of different cosmology or new physical science.
