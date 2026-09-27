# Local gradients in the synthetic optical BayeSN experiment

No smooth-region gradient discrepancy was found in this bounded check. All **422** declared smooth coordinate comparisons passed, including every coordinate at two saved states associated with divergent transitions. A deliberate exact phase-boundary control showed a derivative jump. This does not identify the cause of the original sampling failures: saved transition states are not the failed leapfrog locations, and local checks cannot establish global gradient correctness or adequate posterior sampling.

The audit used only the synthetic optical data from case00 and completed chain files. The frozen M20 `Kernel` and `numpyro_model` were unchanged. The six retained states were declared in advance: LCDM-prior chain0/draw500, chain1/first flagged draw963, chain3/draw999; broad-distance-prior chain0/draw500, chain2/first flagged draw68, and chain2/draw999. Three controls kept the first state's other coordinates fixed and placed its closest optical phase exactly at integer9, or offset tau by ±0.0001day. No NIR outcome, generating-truth file or observed photometry was opened; an audit hook blocked those paths.

All47 coordinates were checked in the sampler's unconstrained parameterization, including the positive/interval transform Jacobians. Autodiff of the frozen NumPyro density was compared with central differences of an independently summed SciPy density using the same forward flux operator. Eleven relative steps from0.01 to10⁻⁷ were retained for each coordinate. A smooth comparison required two adjacent steps that did not cross an integer-phase boundary and had `abs(AD−FD)/(1+abs(AD)) < 10⁻⁴`. Every supported comparison passed; the largest, across coordinates, of the best supported errors was4.50×10⁻⁹. All4,653 step records are retained. The independent density agreed within2.84×10⁻¹⁴, and the actual NumPyro potential gradients were bitwise equal to the explicitly transformed density gradients at all nine states. A separate saved-record verifier checks this equality and input identities without reevaluating the model.

The exact phase control is intentionally nonsmooth. At the smallest physical tau step,10⁻⁷day, the left and right log-density derivatives were approximately **−1.681770** and **−1.649321** per day, a jump of0.03245. The pointwise autodiff value was−1.380122. The released Hsiao interpolation uses floor and ceil phase indices; exactly at an integer both indices select the same template, giving a special pointwise derivative. There is no unique ordinary derivative at that boundary. The two nearby controls passed their smooth-region checks. This establishes an interpolation feature, not a general smooth autodiff bug or proof that crossing this feature caused the recorded divergences.

The one-CPU calculation took102seconds. The verifier rehashed all134 opened upstream files:129 matched the frozen forward-gate registry, and five were recorded Python bytecode caches. Existing model sources, likelihoods, priors and sampler outcomes remain unchanged. The original failed sampling gates remain failed; no new chains, smoothing prescription, held-out score or cosmological conclusion is produced.

From the repository root, use the pinned BayeSN environment and a fresh output directory:

```bash
BAYESN_ARCHIVE="$HOME/.local/share/cosmology-revalidation/archive/17487bf659fcbdeeea072221492bac14b04a0a85"
BAYESN_PYTHON="$BAYESN_ARCHIVE/runs/research_2026_09_26/bayesn_distance_identification/.venv/bin/python"
"$BAYESN_PYTHON" studies/infrared/code/bayesn_synthetic_gradients.py \
  --work .work/infrared/heldout-synthetic --output NEW_GRADIENT_AUDIT_FOLDER
python studies/infrared/code/bayesn_synthetic_gradient_summary.py \
  --audit-folder NEW_GRADIENT_AUDIT_FOLDER --output NEW_SUMMARY_JSON
```

The exact selection and failure rules are in `code/bayesn-synthetic-gradient-design.json`; the compact verified record is `results/bayesn-synthetic-gradients.json`. Full step records and the pre-evaluation state manifest are under the ignored `.work/infrared/gradient-audit-v1` directory. Reproduction requires the finished optical chain inputs identified by SHA256 in that record. Neither script restarts a chain.
