# Joint lensing at native accuracy two

This is a prepared sensitivity calculation, with no new cosmological result yet. It replaces the existing ACT/Planck and Pan2018 SPT lensing terms by the released Qu ACT/Planck/MUSE joint likelihood, separately for its 35-bin baseline and 38-bin extended versions. The two versions are fixed before evaluation. The baseline changes the SPT estimator and joint covariance convention; the extended version also adds ACT multipoles. Neither comparison isolates covariance alone.

The consumer requires `native_accuracy_qualify.qualify_run` to freshly qualify the full accuracy-two posterior. It preserves all 2,000 selected slots, repeated selections, independent-chain groups and original proposal densities. Each prediction uses the corresponding accuracy-two raw spectrum, with the released response matrices and full joint precision. It performs no CAMB, background or model-construction calls. All remaining primary-CMB, BAO, SN and prior factors remain unchanged.

For each point the new log weight is

$$
\log w_{\rm joint}=\log w_{2}
+\log L_{\rm joint}-\log L_{\rm ACT+Planck,2}-\log L_{\rm PanSPT,2}.
$$

The original normalized uniform prior on the SPT nuisance parameter $A_{\rm fg}\in[0,2]$ remains in the extended target. This auxiliary coordinate is included in overlap and chain/batch diagnostics, then marginalized out of physical reporting. Removing that dimension before density accounting would give the wrong replacement.

Original accuracy-one spectral sidecars supply only the two fixed Gaussian normalization constants. Every sidecar is checked against its sealed original parent, index, physical point, density, derived values and capture-source hashes, and the constants must be identical across all slots. Their accuracy-one spectral arrays are never opened. Actual predictions use only the typed accuracy-two archives. The new constant is calculated from the released Hartlap-adjusted precision, so the effective covariance is $C/h$, not $C$. These fixed offsets cancel normalized posterior weights; comparisons across these different data vectors and units are not Bayes factors.

Both variants retain their own unchanged importance and chain/batch gates. Failure withholds the posterior, sign fractions and published covariance, while preserving diagnostic weights and the failed tests. The fresh cache contains an immutable plan, sealed point records and complete source/input hashes. No incomplete cache is silently resumed or overwritten.

After the parent qualifies, the reusable commands are:

```bash
.work/unified-cosmology/external-probes/.modern-venv/bin/python studies/unified_cosmology/code/inference/native_accuracy_lensing_validate.py
.work/unified-cosmology/external-probes/.modern-venv/bin/python studies/unified_cosmology/code/inference/native_accuracy_lensing.py --work PATH_TO_NATIVE2_WORK --summary PATH_TO_NATIVE2_SUMMARY --cache .work/unified-cosmology/native2-joint-lensing-FRESH --output PATH_TO_NEW_RESULT
```

The validation uses independent density algebra, synthetic qualified/failed parent fixtures, and existing archived spectra evaluated by both the compressed and unmodified released likelihood. It makes no new physical calculations. The observational run remains dependent on a genuinely qualified accuracy-two parent. Finite importance diagnostics cannot establish unvisited-mode coverage; the released Gaussian MUSE compression and omitted primary–lensing and other cross-probe covariance remain explicit limitations.
