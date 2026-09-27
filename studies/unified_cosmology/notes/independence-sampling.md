# A frozen global proposal for the same cosmological likelihood

A bounded LCDM pilot accepted **258 of 399 proposed moves (64.7%)** after its
initial state. Its 400 independently requested proposal points had a raw
importance-weight effective count of **289.2**, with largest normalized weight
**0.00859**. There were 259 occupied states, no hold longer than ten evaluations,
and 36 prior rejections, all retained. This is evidence for a more efficient
sampling method; it is **not a new cosmological measurement**.

The [local LCDM pilot record](../results/inference/independence-pilot.json)
preserves the initial efficiency experiment. The pilot learned an eleven-dimensional location and covariance from a frozen
snapshot of the ongoing chains. Each chain's integer frequency weights were
expanded; its first half was discarded, the following 40% used for learning,
and the final 10% withheld. The proposal is 90% Gaussian with covariance
$1.05^2 C$ and 10% Student-t with five degrees of freedom and scale matrix
$4(5-2)C/5$. Thus the Student component has covariance $4C$, rather than scale
matrix $4C$. A diagonal regularizer adds only $10^{-8}\operatorname{diag}(C)$.

The density is normalized and positive everywhere in the sampled real
coordinates. At a proposed state $y$ from current state $x$, the acceptance
probability is

$$
\alpha(x,y)=\min\left[1,
\exp\{\log\pi(y)-\log\pi(x)+\log q(x)-\log q(y)\}\right].
$$

This correction is essential because the proposal is not the target. Learning
is frozen before sampling; a training-chain bias affects efficiency rather than
changing the declared prior. Prior/domain failures are rejected without redrawing
the candidate. Independent production random-number streams are separate from
training and the pilot. Full support and finite convergence diagnostics do not
prove that distant posterior modes have been adequately visited.

## Separate CPL efficiency pilot

A separately frozen thirteen-dimensional CPL proposal was tested using 400
new candidates, seeds 273260/273261, and the validated native initialization.
It accepted **235/399 moves (58.9%)**; 367 targets were finite and all 33 invalid
candidates were prior rejections. There were 236 occupied states, a longest
hold of 14, and a raw candidate-weight effective count of **210.5/400**; the
largest normalized candidate weight was 0.0322. Five successful native-spectrum
fallbacks contributed to 148.4 seconds of evaluation time. Short-chain bulk ESS
ranged from 66.5 to 211.8 and does not qualify a posterior.

The [CPL pilot record](../results/inference/independence-cpl-pilot.json) and
[independent arithmetic audit](../results/inference/independence-cpl-pilot-validation.json)
preserve this efficiency test. SciPy mixture densities agree to
$1.29\times10^{-12}$, prior-plus-component accounting to
$1.14\times10^{-13}$, and all candidate decisions and holds reproduce.
The result supports a fresh independent CPL run under the same convergence and
native-correction requirements; it does not establish its scientific posterior.
Suggested fresh CPL production seed roots are 273262/273263, distinct from
LCDM, source-chain and pilot streams. Only the controlling task launches runs.

## Numerical checks and preserved failures

The original pilot did not set `CLIPY_NOJAX=1`; the baseline replay did. The
replay used exactly the same 400 candidates and random acceptance numbers.
**Every likelihood component and every decision was identical in that same
order.** Baseline evaluation time was 93.1 seconds, including four slow calls
that together took 71.8 seconds; median evaluation time was 0.059 seconds. These
measurements came from a shared machine and are not a production speed guarantee.
Single-chain, 400-step ESS estimates are too short to qualify a posterior.

An eight-point reordered comparison failed its predeclared $10^{-8}$ absolute
log-density tolerance. The largest difference, $3.66\times10^{-4}$, came solely
from DESI BAO. All CMB, supernova and expansion-diagnostic values were unchanged.
This was **not demonstrated to be a CLIPY discrepancy**. A controlled
`cached=False` sequence isolated a change in CAMB's background drag radius after
its first native-spectrum evaluation: at the fixed diagnostic point,
$r_d$ changed from 147.2555386193 to 147.2555968905 Mpc while spectra, $H(z)$ and
angular distances were unchanged. The failed closure is retained. Three distinct native-spectrum controls left
the subsequent warm state identical. The new sampler performs one explicit
native CAMB initialization at a fixed declared point before sampled densities.
All eight diagnostic points then matched the warm pilot exactly in both forward
and reverse orders. [Initialization validation](../results/inference/independence-runtime-validation.json).
This is a numerical initialization, not a changed physical prior. Runtime metadata
distinguish NumPy coordinate precision from third-party internal dtypes and
record actual JAX availability/settings without changing them.

A separate [full-native order check](../results/inference/native-order-review.json)
uses two fixed, previously evaluated controls in the order A, B, A in a fresh
process, with three forced uncached evaluations and no extra warmup spectrum.
The first and revisited A values, every likelihood component and every derived
quantity agree exactly at stored precision; their spectrum hashes are identical.
Both first occurrences also exactly reproduce the saved full-native controls.
Thus the earlier tested cold-background discrepancy is not reproduced in this
bounded full-native check. Two points do not establish universal order
independence or numerical convergence. The [fixed design](../code/inference/native-order-design.json)
and [reproducer](../code/inference/native_order_review.py) retain the complete comparison.

The reusable sampler is a narrow subclass of the installed Cobaya MCMC. It
replaces the proposal transition and includes the Hastings density ratio while
retaining Cobaya's original completed-hold collection, burn-in and convergence
equations. Proposal adaptation, dragging, thinning and non-unit temperature are
forbidden. Its [design](../code/inference/independence-sampling-design.json)
requires the original multivariate $R-1<0.005$ at two successive checks and the
original 95% bound agreement below 0.05. It additionally requires rank
$\widehat R\leq1.01$ and bulk and tail ESS at least 400 for every sampled and
existing diagnostic parameter, using the same segments as the unchanged
`diagnostics.py` consumer.

Checkpoint publication waits until both sets of checks pass. The additional
checks read completely flushed chain text, including its stored numerical
precision; a failed additional gate cannot leave a provisional convergence
flag on disk. Every candidate ledger contains the current and proposed target
and proposal densities, component, decision and holding weight. At normal
stopping, all completed post-burn holds have been written, and only the new
terminal state of weight one remains outside the standard Cobaya collection.
No rejected holding time is lost. Interrupted runs are preserved and cannot be
silently resumed or overwritten by this first implementation.

The actual four-chain synthetic validation used a target different from the
proposal. It independently reconstructed every candidate decision and every
recorded holding weight, including prior rejections. The final original
multivariate values were 0.00294 and 0.000802, bound agreement 0.0384; additional
rank diagnostics were below 1.002 with bulk/tail ESS above 1,993. The final
flushed-data diagnostics exactly matched the sampler's additional gate. Known
target mean errors were below 0.004 and standard-deviation errors below 0.007.
The final audit reconstructed all 16,406 candidate decisions, 6,948 completed
states and 717 prior rejections.
[Validation record](../results/inference/independence-sampler-validation.json).

A first synthetic invocation failed before sampling because Cobaya requires a
sampler's importable class name instead of a likelihood-style `external` option.
That interface was corrected. A subsequent audit exposed the difference between
full-precision in-memory diagnostics and rounded output text; the sampler now
checks the actual flushed text without weakening a threshold. These development
attempts and their logs remain in ignored working storage.

## Reproduction and output contract

Run the synthetic exercise and independent audit with separate fresh output
folders:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CLIPY_NOJAX=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/mpiexec -n 4 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/inference/independence_validate.py synthetic \
  .work/unified-cosmology/inference/independence-validation/FRESH

OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CLIPY_NOJAX=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/inference/independence_validate.py audit \
  .work/unified-cosmology/inference/independence-validation/FRESH
```

`independence_run.py freeze SOURCE_FOLDER PROPOSAL_FOLDER` copies complete
source rows and immutable rank manifests before learning. When the chain rows
already come from a separate frozen snapshot, `--manifest-folder` supplies the
matching original immutable manifests. This is proposal training, and does not
certify the source chains as a cosmological result.

The `sample PROPOSAL_FOLDER FRESH_OUTPUT_FOLDER` command requires exactly four
MPI ranks, the baseline one-thread NumPy environment, and an unchanged scientific
target. Its immutable `run-[0-3].json` records preserve the established
`target_identity` and `arguments` contract, while sampler/proposal/runtime
provenance is separate. Standard `chain.[1-4].txt` columns include integer
`weight`, `minuslogpost`, sampled/derived parameters and component contributions.
The unchanged diagnostic, exact-native-correction and qualified-summary
consumers remain mandatory. No production launch or qualified scientific result
is supplied by this note.
