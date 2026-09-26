# Phase 2 literature and selection audit

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../README.md) for current execution status.

The [machine-readable map](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/literature/literature-map.json) contains **61 pinned sources**, including 20 supplementary primary-paper acquisitions. It separates 28 focused method checks, 20 carried phase-1 method audits, and 13 carried supporting summaries. This is a broad dependency map, **not a claim that every paper has received a new complete methodological review**. Each entry records exact source/version and SHA-256, data level, likelihood, population/dust/age/colour/stretch/selection assumptions, priors, sample reuse, critiques, counterarguments, reproduction limits and preregistered tests. [Schema/hash validation](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/literature/literature-map-validation.json) passes for 147 source references. The [source ledger](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/literature/source-hash-ledger.json) preserves original objects; phase 1 was not modified.

## Consequential methodological findings

- [Original DES methods](https://arxiv.org/abs/2401.02945v2), §§3–5, retain 1,635 DES events in a BEAMS mixture; 1,499 have nominal classification probability above 0.5. The latter count is not the complete likelihood sample. Bias corrections, residual host steps and BEAMS-adjusted errors are part of the released distance product, not raw observations.
- [Vincenzi host selection](https://arxiv.org/abs/2012.07180v1), §3, defines host-redshift success relative to **galaxies satisfying OzDES targeting criteria**, conditional on magnitude, observed host colour, field and discovery epoch. This is not the entire galaxy population or an unconditional SN detection rate.
- [PS1 photometric cosmology](https://arxiv.org/abs/1710.00846v3), §3/Table 4, explicitly fits Ia/CC mixtures and class-probability calibration. Its PS1 events and nearby anchors overlap other compilations. Different analysis pipelines or sample masks do not provide independent observations.
- [Nielsen et al.](https://arxiv.org/abs/1506.01354v3), Eqs 4–10, adopt global independent Gaussian latent magnitude/stretch/colour populations. [Rubin–Hayden](https://arxiv.org/abs/1610.08972v3) show that survey/redshift-dependent populations alter latent shrinkage and acceleration significance. A raw hierarchy must test that population assumption, while also recognizing that unconstrained luminosity drift can remain degenerate with cosmology.
- [Efstathiou v3](https://arxiv.org/abs/2408.07175v3) demonstrates sensitivity to small inter-compilation offsets. [DES's response](https://arxiv.org/abs/2501.06664v1) shows why differing selection functions, host estimates and scatter prescriptions generate different corrected distances for identical events. These are competing explanations and sensitivity tests; neither is an unconditional proof that all supernova evolution is absent or that a measured offset is entirely calibration error.
- [DESI DR2](https://arxiv.org/abs/2503.14738v3), §IV, distinguishes BAO with a free ruler, BBN calibration, acoustic-scale/early-universe compressed priors, and full CMB likelihoods. [Planck](https://arxiv.org/abs/1807.06209v4), §2, requires a perturbation/foreground model. A Gaussian Ωm prior is not a full CMB reproduction, and evolving dark energy is not synonymous with present deceleration.

The modern age exchange and its corrections, September ZTF/host updates, BayeSN, Dovekie, Union3.1, and revised directional responses are carried with explicit phase-1 limits. In particular, the current Ray result differs from its withdrawn deceleration version; host mean stellar age is not observed progenitor delay; and a small fitted residual mass step does not mean that dust/host corrections are absent from the distances.

## What the released mocks actually contain

The [computed denominator audit](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/literature/selection-audit.json) checks every released original v1.3 mock README, DUMP and HEAD table:

| Type | Realizations | Generated | Written/HEAD/DUMP | Overall acceptance |
|---|---:|---:|---:|---:|
| Ia | 25 | 662,950 | 71,946 | 0.108524 |
| non-Ia | 25 | 12,347,500 | 152,959 | 0.012388 |

Every DUMP declares `Pass Trigger + Cuts`; every HEAD search mask is 5 (pipeline + host-z). Counts match the written totals. These are **selected simulation light curves**, before the complete final fitted-SALT/analysis selection. The README generated totals permit an overall acceptance estimate under the original generated distribution, but not a multidimensional efficiency fit from selected rows alone.

Selected Ia truth (`phase2/literature/selected-ia-truth.csv.gz`; historical local reference) is exported with realization/CID, field, true and measured host properties, source SALT parameters, dust variables, and a candidate released-grid density. It contains 71,946 finite grid densities when conditioned on **true** host mass. Using `HOSTGAL_LOGMASS` instead incorrectly assigns missing/mismatched observed host properties to the generated prior.

The historical simulation version is `v11_05c-61-g2fe0f56`. Its pinned source shows `SIM_SALT2c` and `SIM_SALT2mB` describe the SALT source before separate host extinction is applied. Thus `c_intrinsic + AV/RV` and a corresponding additive magnitude approximation are not verified noiseless fitted SALT summaries. Across every released Ia row,

`SIM_mB − SIM_mu + alpha*SIM_x1 − beta*SIM_c = −19.365 ± 0.0000029 mag`,

and `SIM_MAGSMEAR_COH = 0`. At fixed full latent conditions the source luminosity relation is effectively singular. Reweighting to a freely shifted luminosity/cosmology distribution therefore requires an explicit support proof or new forward injections.

## Available empirical components and remaining likelihood work

[selection_assets.py](../code/literature/selection_assets.py) implements the released conditional GENPDF grids, empirical host-z grid and per-epoch detection table plus historical two-period trigger combinatorics. [Validation](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/literature/selection-validation.json) checks density normalization, trigger probability against complete enumeration of 32 detection states, and host-map probabilities across every field/season and extrapolation boundary.

The trigger uses **expected `SNR_CALC`**, not measured `FLUX/FLUXERR`, in historical SNANA. The host grid uses `r_obs_auto`, `obs_gr_auto`, field and `PEAKMJD`; it does not expose per-cell target/success counts needed for a binomial uncertainty model. Its evaluated distribution on selected mocks has median 0.874, which is descriptive rather than an unbiased population mean.

The [selection contract](../specifications/literature/selection-contract.json) specifies the complete conditioning and likelihood normalization. For identical fixed selection and exact proposal density with adequate support, selected-simulation importance sampling can identify **relative** normalization through `E_selected,p0[pθ/p0] = Zθ/Z0`. That is a legitimate route, but source-file equivalence, conditioning, measured-fit selection and support must all be established. Colour-only stress shifts of +0.02/+0.05/+0.10 yield ESS fractions 0.866/0.414/0.0476; these quantify overlap loss, not a cosmological result or completed selection correction.

The full likelihood still needs quality cuts evaluated on **fitted noisy light curves**, target/host association, class calibration, and any final membership requirements (outlier cut, valid bias correction and intersection across systematic variants). The hierarchy must either reproduce these or explicitly define a new sample. Treating post-fit `PBEAMS` as an independent prior double-counts distance information; photometric class probabilities from the same flux also require conditional treatment.

## Preregistered comparison and execution

The map fixes tests `SEL-01`–`SEL-04`, `POP-01`, `CONT-01`, `CAL-01`, `PRIOR-01`, `COS-01`, and `OVL-01` before new cosmological comparisons. Models are assessed using joint colour/stretch/host distributions and held-out within-redshift magnitude patterns with free distance-bin offsets, rather than accepting a population model because it yields a preferred cosmology. Nominal-mock recovery alone is insufficient; held-out population/luminosity alternatives and coverage are required.

Rebuild with:

```bash
python scripts/phase2/literature/build_map.py
uv run --no-project --with jsonschema==4.25.1 python scripts/phase2/literature/validate_map.py
.venv/bin/python scripts/phase2/literature/audit_selection.py
.venv/bin/python scripts/phase2/literature/check_selection_assets.py
```

The supplementary acquisition script is separately reproducible and refuses to overwrite existing paper objects. All simulation/toolchain work lives under phase 2. This report contains no new cosmological fit.

## Forward injections and tests of the simulator itself

The [predeclared forward pilot](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/literature/simulations/inputs/preregistration-pilot02.json) prepares P21, BS21 and G10 physical flux simulations plus P21 luminosity offsets of ±0.10 mag and +0.20z mag. Each requests 26,518 generated attempts, retaining failures through `SIMGEN_DUMPALL`. Released 2024 survey/population assets are used explicitly; reconstructed BS21/G10 choices are listed and are not asserted to be byte-identical historical systematic jobs. A common random seed does not by itself guarantee paired host/cadence draws across different generators. The complete generated-versus-recorded denominator must be audited after execution.

W22's released stretch distribution requires `SN_age`. The nominal DES host library does not contain it. The public Wiseman host-library builder needs an external HDF5 host model and sibling `SN_ages/*.dat` distributions; the exact inputs and the crosswalk to the DES cadence host library are not in the acquired release. An exact W22 run is therefore pending a defined age generator; replacing age with host mass would change the model.

The expanded [selection contract](../specifications/literature/selection-contract.json) and [hashed assumption sources](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/literature/selection-assumption-source-ledger.json) separate empirical calibration from imposed assumptions:

- The DES SIMLIB records 10,000 observing-position libraries, with exposure/MJD/filter/gain/sky/PSF/zero-point information derived from DiffImg. Sampling and compressing those conditions remains a model choice.
- The released flux-error map reports 421,634 fake and 390,257 simulated observations and calibrates errors versus field, band and host surface brightness. It sets same-band correlation parameter `REDCOV=0.6` for the **added error component**, not the complete flux error; some edge bins have no fake observations. These features need sensitivity checks, not an assumption of exact known noise.
- Dust and colour distributions are fitted through selected DES observables; they are not a directly observed parent population. A simulation that recovers its own inputs establishes code closure, while a test of physical robustness must inject a different population, noise model or luminosity drift.
- The fiducial distance law changes flux and detection probability. Conditioning on redshift can omit event-rate information, but cannot remove luminosity-dependent selection. Flexible cosmological alternatives require adequate new injection support.

Independent falsifiers include held-out real flux residual tails and correlations, reported-error distributions at matched observing conditions, wavelength/phase residual patterns, and joint colour/stretch/host distributions within redshift bins. SMP fixes reference-epoch transient flux to zero, so those fitted reference values are not independent evidence for that assumption. Model comparison uses free redshift-bin magnitude intercepts before any cosmological interpretation.

The [first released-mock mapping audit](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/literature/mock-fit-mapping-audit.json) joins an independent calibrated-flux SALT fit to physical truth: 1,936 fitted events, 1,384 passing the explicitly listed measured quality subset. It reports both naive source-plus-dust proxies and the known −0.12 mag generation-offset adjustment. The residuals combine noise, selection and spectral projection; they do not certify an exact noiseless SALT-to-dust mapping or establish cosmological bias.

The [DES photometry paper](https://arxiv.org/abs/2406.05046) describes calibrated light-curve products, separate DiffImg and SMP pipelines, and atmospheric corrections. [NOIRLab's image-service documentation](https://datalab.noirlab.edu/docs/manual/UsingAstroDataLab/DataAccessInterfaces/SimpleImageAccessSIA/SimpleImageAccessSIA.html) lists public DES DR2 single-epoch and raw/calibrated archive services. Availability of those services is not an exhaustive verification that every original SN search/reference image, PSF, mask and calibration state is recoverable. This task has not performed that crossmatch or a new pixel-level SMP/image-injection analysis.

The revised nominal forward run has now completed: **26,518 generated attempts → 2,678 written light curves**. [Run manifest](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/literature/simulations/manifests/PH2_pilot02_P21.json) and [denominator audit](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/literature/simulations/PH2_pilot02_P21-denominator.json) record the result. The historical-code initialization crash was resolved by relocating an exact constant alpha from the population file to the simulation input, with original and derived hashes preserved. Early rejected attempts reuse CIDs, so the unique generated denominator key is the row index; written membership matches only the unique pass rows. [The handoff](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/literature/simulations/HANDOFF.md) includes executable commands for the other eight prepared variants and outstanding scientific limits. Their flux fitting and model ranking are not yet complete.
