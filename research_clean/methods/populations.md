# Host ages and progenitor populations

## Published age extraction

`ages` extracts both literal tables from the bundled Chung supplementary ZIP, joins physical IDs to the Pantheon+ SDSS rows and retains unmatched records. The tables contain 199 G11 and 102 R19 rows. Their overlap is resolved explicitly using `G11_first` by default or `R19_first`; duplicates are not treated as independent supernovae.

With `0.06 < zHD < 0.42` and G11-first, expect 196 unique selected objects, 226 matched age rows before deduplication, and 75 unmatched age rows. `all_age_crosswalk.csv` contains the complete join; `selected.csv` is the analysis sample.

Three outcomes make already-applied corrections explicit:

- `hr_corrected = MU_SH0ES - mu_reference`;
- `hr_no_bias = hr_corrected + biasCor_m_b`;
- `hr_tripp = mB + 0.148*x1 - 3.112*c + 19.253 - mu_reference`.

The reference uses flat LCDM, `H0=70`, `Omega_m=0.3`, `zHD` for the radial integral and `zHEL` for the photon redshift factor. Reversing an exported bias term is not equivalent to reconstructing a complete raw Tripp pipeline.

WLS is a diagnostic with quoted magnitude errors and no age uncertainty. The default corrected-distance slope is about `-0.00417 ± 0.00642 mag/Gyr`; the uncorrected-bias diagnostic is about `-0.01327 ± 0.00642`. These standard errors condition on the chosen sample and model.

The optional Gaussian latent-age fit instead models `a_true ~ Normal(m,tau²)` and `a_observed | a_true ~ Normal(a_true,sigma_a²)`. Conditioning on the age observation gives an analytic mean and variance; the magnitude likelihood includes the full released covariance, latent-age variance times slope squared, and free intrinsic scatter. Three bounded optimizer starts are retained. A boundary solution is flagged. This is a Gaussian-summary sensitivity, not original per-object age PDFs, original LINMIX, or a population/selection reconstruction. Only MLEs are supplied; no posterior probability or profile interval is inferred from them.

```bash
uv run --frozen python research.py run ages --name ages-g11
uv run --frozen python research.py run ages --name ages-r19 --config configs/ages-r19.json
```

## SFH and delay-time convolution

`populations` calculates `p(delay|z) proportional to SFH[t(z)-delay] * DTD(delay)` under a chosen flat radiation-free clock. It supports the B13 and MD14 cosmic star-formation laws and three explicit delay laws: smooth C14, a 0.3-Gyr cutoff power law, and a 0.04-Gyr cutoff with slope -1.13.

The calculation separates cosmic age, mean/median SN progenitor delay, formed-mass mean age, and surviving-mass mean age. These quantities are not interchangeable with fitted host-galaxy age or the distribution selected by a survey.

The default grid uses the specified Son-like CPL clock (`H0=63.6, Omega_m=0.353, w0=-0.42, wa=-1.75`), B13, 12,001 delay points and 60,001 clock points. It evaluates 251 redshifts from zero to 2.5. Astropy clock agreement must be below `2e-6 Gyr`; doubled integration resolution must change retained age summaries by less than `1e-4 Gyr`.

```bash
uv run --frozen python research.py run populations --name populations
```

`delay_curves.csv` contains all three delay laws. `mean_correction.csv` and `median_correction.csv` apply an imposed `0.030 mag/Gyr` slope to the decrease from the z=0 age. They are normalized to zero at z=0. This reproduces a mathematical correction template under stated assumptions; it does not establish the survey's age evolution, the empirical coefficient or an identified correction to distance.
