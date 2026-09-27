# A separate accuracy-two posterior

The original native-CAMB posterior and its numerical screen are preserved. This executor can calculate a distinct posterior with `AccuracyBoost`, `lAccuracyBoost`, and `lSampleBoost` set to 2. All physical parameters, priors, data, calibration factors and other numerical settings remain those of the original native configuration. Implementation and synthetic validation do not establish an observational posterior at the new accuracy.

The prerequisites deliberately differ. An eight-request persistent-worker pilot requires the freshly qualified original posterior and the original 32-point screen, with the auxiliary thermal-background review applied, to identify genuine accuracy-one to accuracy-two density variation. It may precede the accuracy-two to accuracy-three screen. The full calculation additionally requires the latter screen to pass its unchanged numerical thresholds, its separate fresh accuracy-two replay to pass, and the persistent-worker pilot to pass.

The pilot uses exactly two original points: audit indices 0 and 31. Two independent processes request them in the orders `[0,31,0,31]` and `[31,0,31,0]`. Each process constructs one model; there is no hidden warmup or cache reset. Every requested density is compared against its original accuracy-two likelihood components, priors, total density, derived quantities, expansion rates and archived spectra. Repeated evaluations also agree within the registered tolerances. Construction and request times are reported separately. An invocation count is not an internal CAMB solver count.

The full calculation keeps all 2,000 selected slots, including repeated expanded-chain indices. It never recomputes the original proposal density. It reuses the original 32 accuracy-two calculations at their exact slot indices after checking all identities and converting their archived angular spectra under the validated provider convention. The remaining 1,968 slots are assigned to at most four persistent worker processes. The two equivalent untrimmed weights must agree:

$$\log w_2=\log p_2-\log q=\log w_1+\log p_2-\log p_1.$$

The old and new native densities, original proposal contributions, new derived quantities, raw spectra, slot identity and execution provenance are retained in each record. The unchanged importance-weight, independent-chain and chronological batch gates decide whether the new posterior is admitted. An unsuccessful report has no admitted posterior, even if diagnostic moments can be calculated.

The new `QualifiedNativeAccuracy2` boundary freshly verifies the original parent, both numerical screens, the persistent pilot, numerical configuration, scientific sources, assets, package versions, every record and spectrum, and all statistical gates. It performs no physical calculations. Its records explicitly identify accuracy two; they are never passed off as an original accuracy-one correction or accepted through a substituted filename.

Exclusive tickets, process claims, sealed records and a final file manifest preserve incomplete attempts. A claim without a result has an unknown invocation count; an unclaimed ticket has a known zero count. A failed worker stops its remaining sequence. The same execution folder cannot be retried. A numerical or support failure therefore remains inspectable rather than silently replaced by a later successful run.

The synthetic validator exercises the actual persistent-worker loop with an analytic Gaussian model, both pilot orderings and repeats, the full 2,000-slot reconstruction, the 32/1,968 reuse split, duplicate slots, original-density accounting, new derived values, immutable cache checks and failed statistical gates. Observational parent and numerical-screen boundaries are explicitly mocked in that fixture. Native model/background/spectrum calls are forbidden. These checks establish implementation properties, not the speed or numerical adequacy of real cosmological calculations.

A typical execution order is below. `PY` denotes the pinned modern Python environment; each command requires `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CLIPY_NOJAX=1`. Each work directory and report must be new.

```bash
$PY studies/unified_cosmology/code/inference/native_accuracy_execution_validate.py
$PY studies/unified_cosmology/code/inference/native_accuracy_pilot.py prepare \
  --screen ORIGINAL32.json --thermal-review THERMAL_REVIEW.json --work PILOT_WORK
$PY studies/unified_cosmology/code/inference/native_accuracy_pilot.py execute \
  --work PILOT_WORK --output PILOT_REPORT.json
$PY studies/unified_cosmology/code/inference/native_accuracy_execute.py prepare \
  --screen ORIGINAL32.json --thermal-review THERMAL_REVIEW.json \
  --refinement-work REFINEMENT_WORK --pilot-work PILOT_WORK --work FULL_WORK --workers 4
$PY studies/unified_cosmology/code/inference/native_accuracy_execute.py execute \
  --work FULL_WORK --output FULL_REPORT.json
```

The `report` action verifies an already completed report without evaluating likelihoods. The preparation-v1 contract files remain unchanged. Passing the fixed 32-point accuracy-two to accuracy-three screen cannot prove global numerical convergence; finite importance diagnostics cannot prove missing-mode coverage or validate the physical model.
