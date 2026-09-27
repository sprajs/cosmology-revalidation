# Reading the original and supplemental precision screens

The supplemental consumer distinguishes a repaired auxiliary-background check from a genuinely failed likelihood-precision screen. It preserves both decisions. A verified supplemental result can support the statement **no large numerical variation was detected on the fixed 32 points**; it cannot qualify a posterior, establish convergence with numerical accuracy, or erase the original failed record.

The original precision worker evaluated the higher-accuracy native spectra and likelihoods, then used a different CAMB thermal sampling path for one auxiliary background comparison. The separately validated [thermal review](../results/inference/native-precision-thermal-review-validation.json) recalculates those backgrounds while preserving the native spectra, likelihood values and all thresholds. This [new consumer](../code/inference/native_precision_review_consumer.py) performs no such recalculation. It reads the original screen and already completed thermal review, then:

- Freshly qualifies the unchanged parent through the existing measurement consumer, including its current scientific assets and numerical environment.
- Rebuilds the prescribed chronological 32-point selection from the qualified parent and checks the complete plan, native record seals, log/spectrum hashes, source identities and record manifest.
- Replays every thermal-review calculation using its **stored** backgrounds, preserving nonbackground failures and rejecting incomplete or nonfinite evidence.
- Recomputes both original and reviewed screen summaries and requires exact agreement with their saved results.

The importable interface is:

```python
from native_precision_review_consumer import verify
receipt = verify(original_screen_path, thermal_review_path)
```

The returned `original_screen_status`, `reviewed_screen_status`, diagnostic flags, failures and `point_statuses` retain the comparison. `reviewed_precision_screen_supported` is true only when the verified reviewed status is `no_large_variation_detected_on_fixed32`. `posterior_qualification` and `posterior_reweighting_performed` remain false. Integrity errors raise an exception instead of creating a scientific success or failure. All model, background and spectrum entry points are forbidden inside this consumer.

The [validation record](../results/inference/native-precision-review-consumer-validation.json) uses synthetic densities and backgrounds, with the fresh-parent qualifier and native configuration construction explicitly mocked. Actual files, hashes, seals, selection and consumer arithmetic execute normally. It verifies five outcomes: a background-only repair, genuine likelihood variation, a prior-density failure, an unknown pre-existing failure and a higher-accuracy background-closure failure. Fourteen changed-parent, changed-report, modified-point/log/spectrum, missing-binding, incomplete and NaN cases are rejected. The test makes no cosmological calls and does not represent an observational screen.

Run the synthetic validation with the pinned modern environment:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CLIPY_NOJAX=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/inference/native_precision_review_consumer_validate.py
```

After both actual screen files exist, write a separate, fresh receipt with `native_precision_review_consumer.py --screen ORIGINAL.json --review THERMAL.json --output NEW-RECEIPT.json`. This operation does not change either screen or an existing pipeline receipt. Its finite-screen support remains distinct from posterior qualification and from any subsequent conditional lensing analysis.
