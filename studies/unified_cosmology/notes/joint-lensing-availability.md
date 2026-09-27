# A released joint ACT–Planck–SPT lensing alternative

A reproducible alternative to the current lensing combination is available. It combines ACT DR6, Planck PR4 and **SPT 2019–2020 polarization MUSE lensing** in one Gaussian likelihood, including cross-experiment lensing covariance. It would replace both the current ACT–Planck lensing block and the separate Pan et al. SPT 2018 temperature-lensing block. Adding it to those blocks would repeat information. No measurement or active target has been changed.

The official [NASA LAMBDA release](https://lambda.gsfc.nasa.gov/product/spt/spt3g_act_spt_joint_prod_get.html) supplies a 98,518-byte joint-likelihood archive. Its likelihood source is byte-identical to [the author repository at commit cfd88b8](https://github.com/qujia7/spt_act_likelihood/tree/cfd88b8ebe24e969d407637fd2ac20a31220ac29). Bandpowers, windows, covariances and the SPT compression are included. ACT/Planck response matrices come from the separately released ACT v1.2 assets already acquired here. The [audit record](../results/inference/joint-lensing-availability.json) pins the official download URLs, archive hashes, source and input bytes. This verifies evaluation of the released likelihood, not independent regeneration of the reconstruction or covariance from sky simulations.

## What is jointly modelled

The ACT–Planck cross block uses 480 FFP10 common-sky simulations. The SPT cross blocks use analytic approximations, checked against noisier simulation estimates; they are not a fully simulation-estimated joint covariance. The paper uses the extended ACT range to L=1300. Its tests of covariance approximations concern its stated ΛCDM and neutrino-mass analyses, not a bound on CPL dark-energy uncertainty. [Qu et al., sections III, VII.1 and VIII](https://arxiv.org/html/2504.20038v2#S7.SS1)

Direct inspection and execution of the released data gives:

| Configuration with primary CMB | Baseline ACT range | Extended ACT range |
|---|---:|---:|
| Lensing vector order | 10 ACT + 9 Planck + 16 SPT | 13 ACT + 9 Planck + 16 SPT |
| ACT window support | L=40–763 | L=40–1300 |
| Planck / SPT window support | L=8–400 / 20–3073 | Same |
| Largest absolute ACT–Planck correlation | 0.206284 | 0.206284 |
| Largest absolute ACT–SPT correlation | 0.157377 | 0.157377 |
| Largest absolute Planck–SPT correlation | 0.109159 | 0.109159 |
| Released Hartlap precision factor | 0.909774 | 0.902256 |
| Fiducial χ² reproduced | 38.669010 | 41.730392 |
| Author's rounded test χ² | 38.67 | 41.73 |

The covariance is positive definite; its asymmetry is only 3.45×10⁻¹⁵ in correlation units. Independently evaluating the quadratic agrees with the native likelihood within 7.2×10⁻¹⁵ in χ². The two CMB-marginalized variants also reproduce their author tests. These are fixed-spectrum checks, not best fits or cosmological posterior results. Zero native CAMB calls were made.

The joint code applies ACT/Planck normalization and N1 theory-response corrections, while SPT's released 16-bin Gaussian compression is used directly. The SPT bandpower values and windows match the original MUSE release, but its covariance and likelihood representation differ from the original transformed-space product. The latter also contains delensed EE and calibration/other systematics, marginalizing components not explicitly sampled. Its full likelihood cannot be substituted silently for this joint Gaussian product. [Official MUSE release and implementation](https://lambda.gsfc.nasa.gov/data/suborbital/SPT/muse_3g_like_march_2025.zip)

## What remains factorized

No primary TT/TE/EE–lensing cross block is included in this release. SPT reports simulation support for neglecting correlation between its **lensed** primary spectra and MUSE lensing. It explicitly distinguishes the still insufficiently quantified dependence involving MUSE **unlensed EE**. Thus retaining the existing SPT primary likelihood and adding only the new lensing block has a published conditional justification; adding MUSE EE as another independent factor does not. This is not a quantitative validation of the exact current CPL combination. [Camphuis et al., section VII.2.2](https://arxiv.org/html/2506.20707v2#S7.SS2.SSS2)

The remaining ACT/Planck/SPT primary-covariance approximations, and the absent lensing–DESI BAO and lensing–Dovekie cross blocks, remain as described in the [probe-dependence review](probe-dependence.md). The new release improves the lensing-to-lensing treatment; it does not supply a complete covariance for all probes. The baseline alternative preserves the current ACT lensing range; the extended alternative additionally changes scale coverage and needs a separate label.

## Concrete interface and repeatability

The verified generic interface is `load_data('actplanckspt3g_baseline', lens_only=False, like_corrections=True)` followed by `generic_lnlike`. Supply raw convergence C_l and lensed TT/EE/TE/BB C_l in μK², without D_l factors; convert potential spectra using C_l(κκ)=[l(l+1)]² C_l(φφ)/4. Keep the released full covariance and precision correction. No additional joint-lensing nuisance parameter is required.

Two source-level interface pitfalls were executed or inspected: the default Cobaya requirement declaration requests only `pp` despite full-mode response corrections needing TT/EE/TE/BB; a future wrapper should request all spectra explicitly. Also, `indep=True` leaves the full-primary covariance unchanged: that flag's SPT-block deletion is implemented only for `lens_only=True`. Neither issue affects the generic fiducial audit. The archive does not include a reconstruction-level covariance-generation pipeline.

Run from the repository root, after the existing ACT v1.2 acquisition:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/external_probes/joint_lensing_availability.py
```

The script retrieves checksum-pinned small official archives, uses a separate ignored release directory and links existing response inputs read-only. It neither installs over the current likelihood nor changes frozen results.

## A separate posterior sensitivity calculation

A qualified current posterior could support an importance bridge with

`new log weight = original exact/proposal log weight + log L_joint − log L_ACTPlanck − log L_PanSPT`.

Each point needs newly evaluated native spectra under the parent's exact CAMB settings, plus source-component closure against its stored likelihoods. Existing records do not archive C_l arrays. Any new spectral cache must bind parameters, native settings/versions, source and data hashes and immutable payload checksums. Primary CMB, BAO, SN and priors must remain identical; normalized likelihood constants must be recorded consistently even when a common offset cancels normalized weights.

The Pan-only foreground amplitude `A_fg` becomes unused by the replacement lensing data. For a same-dimensional bridge, retain its original normalized uniform [0,2] prior as an auxiliary variable and integrate it out in reporting. Silently deleting or changing that density would invalidate prior cancellation. A newly sampled target could instead remove it explicitly.

This route does not intrinsically require new chains. It does require the unchanged importance-overlap, chain and batch-stability gates to pass on untrimmed weights. A failed bridge would call for independent chains for the new target. No bridge evaluations, chains, or new native spectra were launched in this availability audit.
