# Native CMB precision sensitivity

Eight native calculations compare the declared CAMB accuracy with doubled
`AccuracyBoost`, `lAccuracyBoost` and `lSampleBoost` at four fixed physical
parameter points. All other physics, likelihood inputs and priors remain
unchanged. No spectral interpolation is used. This is a numerical sensitivity
check; it changes neither the active target nor its posterior weights.
Here the reported total includes the declared CMB and BAO likelihoods; these
eight external-probe evaluations do not evaluate a supernova likelihood.

| Fixed point | Doubled minus declared CMB + BAO log likelihood | Time warnings, declared / doubled |
|---|---:|---:|
| Existing reference | +0.087431 | 0 / 0 |
| Broad holdout 34 | −0.202506 | 0 / 0 |
| Broad holdout 37 | −0.779860 | 0 / 0 |
| Near-baseline holdout 90 | −0.223023 | 0 / 0 |

The three declared-accuracy holdout calculations reproduce every original
likelihood component exactly. The tested native accuracy changes are far smaller
than the original cubic interpolation errors of approximately +139.90 and
+110.24 at points 34 and 37. They therefore do not explain those large
interpolation discrepancies. The nonzero native likelihood shifts still matter
for precision assessment. In particular, these four selected points do not
establish convergence of a posterior or show that doubled accuracy is sufficient.
A separate check on qualified posterior points is needed before interpreting
accuracy changes statistically.

At fixed physical parameters, the five sampled values of H(z), CAMB's acoustic
scale approximation, and its H0 inversion at fixed requested acoustic scale are
identical between accuracy settings. The drag sound horizon changes by at most
0.000185 Mpc. CAMB's reported age changes by roughly 0.00021 Gyr, consistent with
the separate finding that its default physical-time integration tolerance is
looser than the independently converged age integral. These are numerical
properties of the assumed cosmological model, not independent age measurements.

The native warning comes from
[CAMB 1.6.6's `CalcScalarSources`](https://github.com/cmbant/CAMB/blob/1.6.6/fortran/cmbmain.f90#L1024-L1033).
At the last source step, it reports a difference greater than 5×10⁻⁵ conformal
Mpc between the perturbation integration time and the tabulated transfer time.
The routine runs for individual wavenumbers, so 269 warning lines in the earlier
interleaved broad-acquisition log do not imply 269 affected cosmologies. The
printed warning gives neither the mismatch size nor an observable error. None
of the eight present calculations emits it. The earlier log cannot securely
identify which original point emitted the warnings; these results do not resolve
that attribution or certify the original warning-producing calculation.

The preregistered points, settings and source hashes are in
[`native-precision-design.json`](../code/inference/native-precision-design.json).
The source audit is
[`native-precision-source-review.json`](../results/inference/native-precision-source-review.json),
and the full component differences, background quantities, spectrum comparisons,
input identities and timings are in
[`native-precision-audit.json`](../results/inference/native-precision-audit.json).
Spectrum differences there cover the entire requested multipole range, including
tails beyond some experiments' fitted ranges; they are not fractional errors on
the observed band powers. Total likelihood comparisons use the actual experiment
windows and covariance matrices.

The exact validated fast-lensing algebra is used in both arms. Final CAMB
settings retain `lmax=9001`, `lens_margin=1250` and
`lens_potential_accuracy=4`; the actual native `max_l` is 10251 and the provider
array length is 10152. The accuracy change also recomputes the drag horizon and
the BAO likelihood. Native calculations take 8.7–18.3 seconds at declared
accuracy and 135.7–186.0 seconds at doubled accuracy on the reserved single
worker in this run; these timings are not a hardware-independent benchmark.

After acquiring the pinned modern likelihoods and the recorded broad input
snapshots, replay or resume the audit with:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/inference/native_precision_audit.py
```

The driver preserves completed points and refuses to silently repeat an
incomplete native attempt. Each point runs in its own subprocess so that buffered
Fortran warnings can be attributed to that point. Original broad-input snapshots
are hash-checked; their complete physical parameter values are also recorded in
the public design. No additional native calls were made beyond the declared eight.
