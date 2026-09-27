# Synthetic optical-to-infrared prediction validation

This executor implements the existing two-cadence, three-truth, two-noise-seed screen. It generates synthetic signed fluxes using quoted source errors, and never fits observed supernova flux. Its working model and limitations are in [design.json](design.json) and the [feasibility note](../../notes/bayesn-heldout-feasibility.md). The earlier J/H forward gate is required and hash-verified.

The likelihood retains all 47 coordinates per supernova: D, AV, RV, shape, 42 whitened spectral residuals and optical-trigger peak time. D is sampled explicitly with the exact proper convolved prior; the other coordinates' marginal target is the same as integrating D out. `validate.py` independently compares adaptive scalar-D integrals, a dense grid and a full-support importance-sampling mixture. The latter checks the joint NIR-vector prediction, including the nonzero zero-amplitude likelihood tail. A product of separately marginalized epochs is deliberately shown to disagree.

Use the preserved isolated BayeSN environment (Python 3.12.14, NumPy 2.2.6, SciPy 1.15.3, JAX/JAXlib 0.6.2, NumPyro 0.19.0, ArviZ 0.21.0, Astropy 7.1.0, extinction 0.4.9, ruamel.yaml 0.18.10), or build an isolated copy from the [recorded complete package versions](../../environment/bayesn-heldout-requirements.txt). No global packages are installed. For a new replay choose fresh output/work paths; existing execution records are never overwritten.

The [input restoration helper](../restore_bayesn_heldout.py) obtains all 187 exact inputs from pinned author releases and historical execution fixtures. It verifies every byte hash against the original forward gate. Historical fixtures include preserved signed-photometry inputs as well as generated metadata/filter tables; restoring them is not an independent re-extraction. Downloaded vendor code, measurements and models stay ignored. The helper does not parse measured flux. A full Git clone can add `--prefer-local-git` to reuse available historical blobs. For a new checkout:

```bash
python studies/infrared/code/restore_bayesn_heldout.py \
  --destination .work/infrared/heldout-inputs
uv venv --python 3.12.14 .work/infrared/bayesn-environment
uv pip sync --python .work/infrared/bayesn-environment/bin/python \
  studies/infrared/environment/bayesn-heldout-requirements.txt
TASK_ARCHIVE="$PWD/.work/infrared/heldout-inputs"
TASK_PYTHON="$PWD/.work/infrared/bayesn-environment/bin/python"
```

The following commands can also use the existing historical archive and environment, as shown. Pass the fresh `TASK_ARCHIVE` and `TASK_PYTHON` above when restoring elsewhere.

```bash
TASK_ARCHIVE=/home/szymon/.local/share/cosmology-revalidation/archive/17487bf659fcbdeeea072221492bac14b04a0a85
TASK_PYTHON="$TASK_ARCHIVE/runs/research_2026_09_26/bayesn_distance_identification/.venv/bin/python"
"$TASK_PYTHON" studies/infrared/code/bayesn_heldout/validate.py \
  --output .work/infrared/heldout-synthetic-validation-replay.json
"$TASK_PYTHON" studies/infrared/code/bayesn_heldout/synthetic.py prepare \
  --archive "$TASK_ARCHIVE" --work .work/infrared/heldout-synthetic-replay \
  --validation .work/infrared/heldout-synthetic-validation-replay.json
timeout --signal=TERM --kill-after=30s 3600 \
  "$TASK_PYTHON" studies/infrared/code/bayesn_heldout/synthetic.py run-case \
  --work .work/infrared/heldout-synthetic-replay --case 0 --workers 2
```

Preparation creates all twelve datasets before fits, with separate optical and synthetic-NIR files and explicit random streams. It verifies the new single-object row layout against the author public simulation interface. Case zero was selected before outcomes: DES16E2clk, zero-dust truth, noise seed 4101. Two worker processes fit its two distance arms; each process uses one CPU and runs four independent chains, each with 1000 warmup and 1000 retained draws. No truncation is accepted as convergence. Workers cannot open synthetic NIR values, truth files or observed photometry. Their identities do not contain a held-out-payload hash that could change optical initialization or fitting when the held-out values change.

The old two-object runtime benchmark projected 237.8 and 153.8 seconds per 1000+1000 chain for the two arms. Multiplying by four gives 15.9 and 10.3 minutes; two simultaneous arm workers would therefore take approximately the longer duration plus setup if that short benchmark transferred. It is not a measured forecast for these noisy single-object posterior geometries. Allow roughly 20–35 minutes for the first complete case, supervised by the 60-minute process-group limit. The internal cap is checked between complete chains, so the external `timeout` is required for a strict in-flight limit. Preserve partial output after timeout; do not restart it in the same directory or call it a passed posterior. The original full-screen allowance is one hour on two CPUs: twelve independent cases may exceed it, so measure the first case and obtain a new resource allocation before running the rest. This is a resource decision, not a change to scientific gates.

Only if both optical arms pass all four-chain gates can the synthetic held-out score run:

```bash
"$TASK_PYTHON" studies/infrared/code/bayesn_heldout/score.py \
  --work .work/infrared/heldout-synthetic-replay --case 0
```

The score freezes every chain and diagnostic hash before reading synthetic NIR. Each posterior draw predicts all held-out bands and epochs with the same latent vector. The output retains full predictive covariance and the normalized joint-vector log density. Numerical error is estimated from four independent chain means and twenty contiguous blocks per chain, using the larger estimate; this is Monte Carlo precision, not extra astrophysical uncertainty or a calibrated population confidence interval. Whole-chain numerical underflow is preserved as failed support, not deleted. Timing/broad-distance boundary mass, all latent-coordinate diagnostics and failed cases remain visible.

The source supports one explicitly selected case per invocation. It does not automatically proceed to the remaining eleven, read observed outcomes, alter a prior after a failed gate, or substitute shorter chains. Synthetic recovery cannot by itself establish that the quoted observational noise, trained spectral population, calibration or selection model is correct.
