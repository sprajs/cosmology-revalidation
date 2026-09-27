# Survey selection, event identity and observational closure

The released Dovekie distance likelihood is ready for a joint expansion fit, conditional on its published calibration, light-curve, selection and contamination treatment. The observed SuperNNova classification can now be reproduced for **all 17,733 released objects**. New public host spectra substantially improve the observational basis for testing population evolution. These advances do not yet constitute a fully regenerated joint physical survey likelihood.

| Route | What is established | What remains conditional |
|---|---|---|
| Released Dovekie distances | 1,820 rows, exact source ordering, full released statistical and systematic covariance, common absolute-magnitude marginalization | Author light-curve model, population and selection correction, BEAMS contamination treatment and calibration systematics |
| Observed DES photons and classification | 19,706 SMP objects; exact candidate joins; epoch-pointer and flux-error checks; all 17,733 released SNN probabilities reproduced to release precision | Population calibration of probabilities, classifier retraining under changed physics, full CC-prior and multidimensional BBC regeneration |
| Historical spectroscopic DES3YR | 207 DES plus 122 low-redshift light curves; spectral classifications; two complete, well-populated bias-training grids; exact released 20-bin likelihood | Historical SALT2/calibration/scatter/selection assumptions; a new age/dust response requires new simulation and likelihood validation |

## Released distance interface

The actual Dovekie file contains **1,623 DES, 117 Foundation and 80 other low-redshift objects**. Its top-level README retains a 1,635-DES description; the data rows govern this analysis. The author metadata table has a different order and is explicitly joined to the distance table by CID and survey. Every distance row joins its original photometry header. Unknown IAU aliases retain a survey-qualified identity instead of an invented cross-survey name.

The `.npz` matrices contain packed upper-triangular **precision**, not covariance. Both total and statistical matrices are positive definite, and their reconstructed inverses close to numerical precision. Subsamples must be formed by taking the corresponding covariance block and then inverting it: taking a DES block of the full precision changes its entries by as much as 14.92. The total matrix already incorporates the released uncertainty model; adding `MUERR²` or multiplying classifier probabilities again would count effects twice. **122 accepted DES rows have published pIa below 0.5** and belong to the BEAMS-weighted likelihood, not to a pure-Ia threshold sample.

The supplied author likelihood has stale default filenames and reads these whitespace tables as CSV. The adapter repairs those packaging issues. Its projected quadratic form agrees with the extracted author likelihood function to $2.85\times10^{-14}$. The redshift convention is $d_L=(1+z_{\rm HEL})(1+z_{\rm HD})D_A(z_{\rm HD})$. A single unknown absolute magnitude is integrated out; supernovae alone therefore do not measure an absolute Hubble scale.

The nominal numerical SALT3 files in the current official bundle match the public Dovekie Git files byte for byte. SMP flux has the documented chromatic/DCR and uncertainty corrections; the matching calibration files supply remaining calibration offsets. This identity check is not an independent detector or calibration reanalysis.

## Candidate and classifier checks

All **19,706 SMP objects** join the **31,636 DIFFIMG candidates** by exact survey identifiers and consistent coordinates. The remaining **11,930 detected candidates** were not in SMP. The 1,779,030 SMP and 18,074,580 DIFFIMG science epochs have valid object pointers, finite fluxes and positive reported errors. Negative fluxes are retained. These products expose individual epoch errors, not a full cross-epoch covariance matrix. The authors explicitly prohibit using DIFFIMG flux for cosmology; it is used here only to establish the detected-candidate denominator.

The initial 256-object classifier diagnostic used public header `PEAKMJD` and failed to reproduce 19 probabilities. Full-catalogue comparison found 615 differences beyond the stated rounding tolerance and 35 disagreements at pIa=0.5. This was a reproducibility error in that diagnostic input, not evidence against the published classifier.

The documented Pippin operation uses separately generated **`PKMJDINI`, `OPT_SETPKMJD=16`, PHOTFLAG rejection mask 1016, and `REDSHIFT_FINAL`**. Applying that measured operator to the same photons reproduced **17,733/17,733** probabilities within $5\times10^{-5}$ rounding plus $10^{-6}$ float32 allowance, with maximum difference $5.0111\times10^{-5}$, and zero threshold disagreements. All 1,623 Dovekie objects pass. For the largest initial discrepancy, CID 1340917, the public peak was MJD 57319.289 but the native estimate is 57372.074; pIa changes from 0.00000374 to 0.998291, matching the released 0.9983. No peak or probability was tuned to obtain agreement.

The prior physical age-injection campaigns already used the correct measured native timing. A read-only audit verified **41 jobs and 55,431 classified occurrences**, matching every stored peak to `PKMJDINI`; option 20 is clump estimation plus a no-abort flag, not simulation-truth bit 2048. Their model weights, arguments and normalization match the newly checked model. Thus the newly discovered header-timing error does not invalidate those historical results. This closure verifies observed inference, not classifier truth accuracy, contaminant population priors or retraining.

## Host follow-up observations

[OzDES DR2](https://cdsarc.cds.unistra.fr/ftp/J/MNRAS/496/19/ReadMe) supplies 38,624 observed targets, including unsuccessful redshift attempts, plus stacked and individual spectra with variance and bad-pixel arrays. Matching **host coordinates**, requiring a unique target within 1 arcsecond, quality 3/4, a host-target tag and \(|\Delta z|\le0.003\), yields **1,089 accepted Dovekie SNe in 1,088 distinct hosts**. Their host redshifts span 0.074–1.145; 673 are above 0.5. Another 3,381 SMP objects outside Dovekie have clean host spectra. The 211 accepted objects with multiple nearby targets remain flagged and excluded from this clean set. One host shared by two SNe must contribute only one host likelihood.

Redshift agreement corroborates identity, not independent measurement: many released redshifts originate in OzDES. The delivered crosswalk preserves association distance, host DDLR/confusion, target identity and exact spectrum URL. Spectra are in counts per wavelength, so continuum response, aperture and normalization require explicit nuisance treatment before physical SED inference; they are not automatically absolutely calibrated host fluxes.

Among observed SN-host targets, secure-redshift success falls from 88.4% at aperture r=21–22 to 54.5% at r=23–24 and 6.48% at r=24–25. Wilson intervals and counts are retained. These are **success probabilities conditional on being observed**, not targeting probabilities or a complete detection-to-redshift selection function. OzDES aperture magnitudes are not interchangeable with Kron magnitudes in simulation efficiency tables.

## New selection assets and remaining production gates

The [official April 2026 SNDATA_ROOT release](https://zenodo.org/records/19503606) corrects 12 of 684 nominal host-efficiency nodes, replacing values above one with valid probabilities. The native `OPT_EXTRAP=1` implementation clamps coordinates to the tabulated boundary and then interpolates; bounded nodal values therefore remain bounded. All grids and 180,000 sampled interior/boundary evaluations pass. The separate `OBS` variant is unchanged and still reaches 1.008: the nominal repair must not be claimed for every variant.

The bundle contains cadence, host libraries, detection efficiency, calibration, SALT surfaces, trained classifiers, CC templates and Pippin/BBC configurations. Separate public [PLAsTiCC Iax and 91bg libraries](https://zenodo.org/records/6672739) were recovered and checksummed. Equality with the particular `model_libs_updates` and binary caches used in production is not established. A fresh Git LFS batch request still returns HTTP 403 for the original 250 Dovekie mock objects. Configurations referencing generated external CC-prior/BBC outputs do not supply those completed outputs automatically; several legacy path aliases also need explicit resolution. We therefore distinguish a recoverable simulation recipe from exact production closure.

The prior sparse-map experiments are not superseded by the host-efficiency repair. A complete new physical measurement still needs a selection-normalized age/dust/host/SN population model, validated contaminants and sufficiently supported correction training, calibration/observation response and coverage. Reproducing published probabilities alone does not establish any of those.

## Spectroscopic alternative and geometry checks

The [original DES3YR archive](https://desdr-server.ncsa.illinois.edu/despublic/sn_files/y3/tar_files/) contains 329 light curves and BBC rows, of which 207 are DES; all 207 join released spectral classifications. Of these, 185 have at least one `SNIa` label and 22 have only provisional `SNIa?` labels. Spectroscopic selection therefore does not establish perfectly known type for every object. Its G10 and C11 bias-training tables contain **1,361,829 and 1,340,209** events, respectively, with each of the four released alpha/beta cells populated by over 328,000 events. These are alternative intrinsic-scatter models, not independent observations or newly fitted physical age distributions. Matching DES3YR spectroscopic-efficiency tables apply to this historical programme, not all typed five-year objects.

The exact author cosmology interface contains **20 bins representing 329 SNe, with 18 nonempty bins**. DES-only has 20 rows and 13 nonempty bins. No 329×329 systematic covariance is supplied. The adapter adds statistical `dmb²` to the author systematic matrix once, following the pinned CosmoMC source with peculiar-velocity and intrinsic-dispersion additions disabled. The model is evaluated at each released bin redshift. A conventional +19.3 magnitude offset in the normalized interface cancels under the free absolute magnitude. At flat \(\Omega_m=0.3,w=-1\), the combined projected chi-square is 17.0022 with total covariance and 20.8182 with statistical errors; these are fixed-geometry checks, not an independently inferred cosmological result. **147 DES3YR objects overlap Dovekie**, so the likelihoods are alternatives.

An independent audit of the joint-inference geometry verified the full-covariance supernova quadratic form against direct FLRW integration, global magnitude/Hubble-scale degeneracy, signed luminosity drift, Gaussian spline marginalization and analytic CPL acceleration/jerk derivatives. It found no material mathematical discrepancy. The resulting cosmological inference remains conditional on the explicitly chosen expansion and luminosity-evolution models.

## Reproduction

Run from the repository root. Downloads, expanded arrays and generated tables remain under `.work/unified-cosmology/survey-selection`; compact findings and hashes remain in `results/survey_selection`.

```bash
.venv/bin/python studies/unified_cosmology/code/survey_selection/acquire.py --bundle
.venv/bin/python studies/unified_cosmology/code/survey_selection/acquire_followup.py
.venv/bin/python studies/unified_cosmology/code/survey_selection/acquire_processing.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python studies/unified_cosmology/code/survey_selection/dovekie.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python studies/unified_cosmology/code/survey_selection/observations.py
.venv/bin/python studies/unified_cosmology/code/survey_selection/ozdes.py
.venv/bin/python studies/unified_cosmology/code/survey_selection/assets.py --plasticc
.venv/bin/python studies/unified_cosmology/code/survey_selection/calibration.py
.venv/bin/python studies/unified_cosmology/code/survey_selection/des3yr.py
.venv/bin/python studies/unified_cosmology/code/survey_selection/des3yr_checks.py
.venv/bin/python studies/unified_cosmology/code/survey_selection/des3yr_likelihood.py
```

The native timing calculation requires a configured, built SNANA source tree and its dependencies. `dataprep.py` accepts `--snana-tree`, `--sndataroot` and `--library-path`, writes its executable and object copies only inside this new workspace, and preserves source/build identities. Defaults refer to the previously restored local environment; on a fresh machine, supply explicit paths following the [native recovery documentation](../../host_ages/code/survey_physics/README.md). The classifier requires Python 3.10, PyTorch 1.13.1 CPU and the pinned SuperNNova revision described there; set `SURVEY_CLASSIFIER_PYTHON` and `SURVEY_CLASSIFIER_SOURCE` to that environment. The initial 256-object diagnostic is intentionally retained before the corrected full comparison.

```bash
.venv/bin/python studies/unified_cosmology/code/survey_selection/dataprep.py
PYTHONPATH="$SURVEY_CLASSIFIER_SOURCE" "$SURVEY_CLASSIFIER_PYTHON" studies/unified_cosmology/code/survey_selection/classifier_closure.py
PYTHONPATH="$SURVEY_CLASSIFIER_SOURCE" "$SURVEY_CLASSIFIER_PYTHON" studies/unified_cosmology/code/survey_selection/classifier_full.py
.venv/bin/python studies/unified_cosmology/code/survey_selection/diagnose_classifier.py
.venv/bin/python studies/unified_cosmology/code/survey_selection/prior_classifier_audit.py
OPENBLAS_NUM_THREADS=1 .work/unified-cosmology/external-probes/.venv/bin/python studies/unified_cosmology/code/survey_selection/review_geometry.py
.venv/bin/python studies/unified_cosmology/code/survey_selection/validate.py
```

The prior-campaign audit additionally requires its already generated historical inputs, and the independent geometry audit requires the external-probe environment. Neither is substituted by an unrecorded mock or missing-data assumption. Archived output timestamps may change when regenerated; scientific values, row ordering and explicit identity checks determine agreement.
