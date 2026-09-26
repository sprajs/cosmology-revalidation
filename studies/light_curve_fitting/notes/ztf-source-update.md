# ZTF BayeSN environmental analysis: source update (20 September 2026)

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../README.md) for current execution status.

## Scope and status

This is an acquisition and source-backed extraction for Ginolin et al., arXiv:2605.06799. It records what the authors state, what the linked public products contain, and what they do not contain. It does not assess whether the model or conclusions are correct, and no regression, cosmology fit, or BayeSN fit was run.

The latest arXiv record at the literature freeze is **v3, submitted 17 September 2026**, titled *On the origin of the environmental dependence of SN Ia magnitudes: A BayeSN view of the ZTF SN Ia DR2*. The record says “Accepted for publication in MNRAS,” 14 pages, 13 figures, 8 tables. Submission history is v1 (7 May), v2 (11 May), and v3 (17 September). Originals, previous-version source, exact hashes, URLs, retrieval times, public data, and pinned code are under `sources/updates/2026-09-20-ztf/` (`sources/updates/2026-09-20-ztf/PROVENANCE.md`; historical local reference).

## What the paper analyses

The measured survey sample is ZTF SN Ia DR2: spectroscopically typed SNe Ia observed from March 2018 through December 2020. The authors state that most classifications came from BTS spectra. The released SN table contains redshifts, light-curve fit parameters, sky positions, classifications, and Milky Way reddening information. The host products contain four environmental proxies: global stellar mass, global rest-frame `(g-z)` colour, local stellar mass, and local rest-frame `(g-z)` colour. “Local” is a 2 kpc-radius aperture around the SN. Host matching uses `d_DLR`; host properties are SED fits to Pan-STARRS photometry. These are the quantities described in **Section 3.1, “ZTF SN Ia DR2.”**

The analysis does **not** use a measured host age or progenitor age. Age appears as a possible physical explanation and in literature comparisons. The fitted environmental axes are mass and colour, local or global. Dust is treated through BayeSN latent `A_V` and `R_V` distributions; intrinsic SED terms and an achromatic absolute-magnitude offset are also model components. Consequently, this paper is evidence about an environmental magnitude dependence under its model, not a released measured-age/Hubble-residual dataset.

The BayeSN training is not trained on the ZTF sample. The authors say the G26 training uses 1,024 SNe from Foundation, CSP, CfA3-4, PS1, SDSS, SNLS, and DES (**Section 2.2, “BayeSN-Dovekie model”**). Thus ZTF supplies distinct photometry for the paper’s application, while the fitted SED model has training connections to surveys used by Pantheon+/DES-era cosmology. The companion Grayling et al. paper, arXiv:2606.19429v1, describes the joint BayeSN×Dovekie training and applies it to DES-SN5YR; its arXiv record had no revision after 17 June by the freeze.

## Selection and transformations stated by the authors

The authors call the selected sample volume limited and state that this design minimises observational bias such as Malmquist bias. **Section 3.2, “Volume-limited sample,”** specifies:

- `z <= 0.06`;
- at least seven 5-sigma detections between phases -10 and +40 days, in at least two bands, with at least two detections before and two after peak;
- SALT fit probability greater than `1e-7`, `-3 <= x1 <= 3`, `sigma_x1 < 1`, `-0.2 <= c <= 0.8`, `sigma_c < 0.1`, and peak-time error below one day;
- removal of peculiar Ia subtypes such as 91bg and Ia-CSM, while retaining 91T;
- removal of objects missing either local or global host measurements.

The paper states a final sample of **932 objects**. It also says selection follows the earlier Ginolin SALT analyses to make the samples identical, and notes that some exclusions depend on bad or extrapolated SALT fits even if BayeSN might fit them acceptably.

The input light curves are released photometric measurements, not raw detector counts. **Section 3.3, “Preprocessing of the light curves,”** says the authors convert them to zero point 27.5, apply the `ztfcosmo` `lightcurve.get_lcdata()` quality procedure, average in 12-hour rolling bins, apply per-band error floors, and compute CMB-frame redshifts. The paper’s model quantities are farther downstream:

- SALT `x0`, `x1`, and `c` are fitted light-curve parameters.
- BayeSN `theta_1`, `A_V`, `R_V`, `mu`, intrinsic residual SED terms, and population hyperparameters are fitted/latent quantities.
- **Section 4.4, “Hubble residuals,”** defines `Delta mu = mu_fit - mu_cosmo`, using flat LambdaCDM with `H0=73.24 km/s/Mpc` and `Omega_m=0.28`, then subtracts the sample mean.
- BayeSN fitting mode keeps distance free and does not apply a cosmology prior, but the population analysis in **Section 5.1, “Method,”** conditions distances on redshift and the stated cosmology to constrain individual `R_V` and population dust parameters (see also **Section 2.1**).

The authors also state that ZTF DR2 is not cosmology ready because of known non-linearities and insufficient photometric accuracy (**Sections 3.1 and 6.6, “Full retraining”**). Section 6.6 quotes a 90 mmag magnitude calibration uncertainty and a 10 mmag colour shift, and says 52% of this analysis sample has only `g` and `r` detections. The paper therefore does not present a ZTF cosmological acceleration fit.

## Reported author claims, with exact locations

The following are transcriptions/condensations of claims made by the authors, not independent confirmations:

- **Section 4.4 and Table 1:** for the paper’s residual sample, BayeSN has `sigma_STD=0.225`, `sigma_nMAD=0.177`; SALT has 0.230 and 0.181. For the host-redshift subset, the respective pairs are 0.204/0.165 and 0.213/0.169. The table caption labels these samples 929 total and 729 host-redshift SNe.
- **Section 4.5, “A posteriori steps,” and Table 2:** BayeSN steps are `0.103 +/- 0.010` mag for global mass and `0.086 +/- 0.010` mag for local colour. The authors say agreement with their SALT values shows that the large ZTF steps are not an artefact of SALT or its standardisation.
- **Section 5.1:** the split model allows the achromatic offset, `R_V` distribution, and `A_V` exponential scale to differ across low/high-mass or blue/red environments. It uses a probability-weighted smooth split to propagate uncertainty in each environmental proxy. The authors fit ten population hyperparameters and state that four MCMC chains of 500 steps each, including 250 warm-up steps, were checked visually and with `R-hat`.
- **Section 5.2 and Table 3:** the authors report no significant `R_V`-mean difference for any of the four proxies. They report a roughly 7-sigma difference in the `A_V` scale `tau_A` for local mass.
- **Section 5.3 and Table 4:** after allowing environment-dependent dust distributions, the reported intrinsic achromatic steps are `0.103 +/- 0.018` mag for global mass and `0.085 +/- 0.019` mag for local colour. The authors describe this as evidence that the environmental step has an intrinsic component.
- **Section 6.1 and Tables 5-6:** repeating the analysis with the T21 BayeSN training changes residual scatter but yields step estimates described as compatible with the fiducial G26 estimates.
- **Section 6.2 and Table 7:** using a sharp rather than smooth split changes the reported intrinsic-step values within their quoted uncertainties.
- **Section 6.3:** an ad hoc broken-`W1` refit gives a global-mass a posteriori step of `0.097 +/- 0.009` mag, which the authors compare with the fiducial `0.103 +/- 0.010` mag.
- **Section 6.4 and Table 8:** allowing the intrinsic-colour covariance `Sigma_epsilon` to vary with global mass gives an intrinsic step of `0.111 +/- 0.020` mag; the authors state that this additional freedom does not change their main result.
- **Section 7, “Conclusion”:** the authors conclude that the environment dependence is at least partly intrinsic under their BayeSN model and argue that it should be accounted for in cosmological analyses.

## Public inputs and reproducibility barriers

The paper’s **Data Availability** section links the ZTF DR2 archive, official BayeSN, the authors’ BayeSN fork, and `standax`. Those products were acquired and pinned. The 1.27 GiB `ztfsniadr2_lite.zip` passes ZIP integrity. Its principal tables each contain 3,628 rows:

- `snia_data.csv`: redshift and uncertainty/source; SALT `t0,x0,x1,c` and errors; the six within-object covariance terms among `t0,x0,x1,c`; Milky Way reddening; fit probability; coordinates; type/subtype; flags; IAU name.
- `globalhost_data.csv`: host coordinates, stellar mass and error, rest-frame `(g-z)` and error, and `d_DLR`.
- `localhost_data.csv`: local stellar mass and error and local rest-frame `(g-z)` and error.

The archive contains calibrated light-curve CSVs and spectra but no author-produced table of per-SN BayeSN `mu`, `theta_1`, `A_V`, `R_V`, residuals, posterior chains, step posterior samples, or cross-object covariance. The SALT table’s covariance terms are within-object light-curve-fit covariances; they are not a covariance matrix for the paper’s Hubble residuals or environmental steps.

The exact analysis mask is not reconstructible from the linked list alone:

- Section 3.2 says 932 objects.
- Table 1 labels the residual sample 929 objects.
- The linked `Ginolin25ab_masterlist.csv` has 944 unique ZTF names and only that one data column.
- Of those 944 names, 933 have a corresponding light-curve CSV in the public lite archive; 11 do not. This still differs from both 932 and 929.
- No released final 932/929 mask or explanation of these count transitions was located in the paper source or linked repositories.

The archived code supplies the G26 model in official BayeSN commit `08eef9e5` (12 August 2026). The cited `smooth_step` analysis branch is pinned at `0b3cac09` (9 February 2026) and does not itself contain the later official `G26_model` directory. Neither code archive contains the article’s fitted ZTF posterior products. An exact result reproduction would therefore require resolving the mask/count transitions and obtaining the per-analysis configuration, trained-model/branch combination, seeds, and posterior outputs rather than substituting the earlier 944-name list.

## Does this enable an independent measured-age residual check?

**No.** The public ZTF tables expose mass and rest-frame colour, locally and globally, but no measured stellar-population age, progenitor age, age uncertainty/posterior, or age-residual covariance. The paper’s BayeSN residual vector is also not released. It would be possible in principle to create a new age estimator from other data and to refit light curves, but that would be a new analysis with new modelling assumptions, not an independent check using the paper’s released measured-age inputs. No age proxy was silently substituted here.

## Overlap with Pantheon+ and DES

The paper does not publish an overlap table. A deterministic identifier audit was performed only to describe the public sample relationship:

- Comparing the 944 linked ZTF names and their ZTF-release IAU aliases with the acquired Pantheon+ table gives three normalised IAU-name matches: SN 2018fop (`ZTF18abtfvsk`), SN 2019np (`ZTF19aacgslb`), and SN 2019ein (`ZTF19aatlmbo`). Only SN 2018fop has a light-curve CSV in the lite archive. Because the final 932/929 mask is unavailable, membership of these objects in each paper result is not established by this audit.
- The same literal/normalised-name check finds no match to the acquired DES-Dovekie distance table.
- This is a lower-bound identifier check, not a sky-coordinate/epoch cross-match. Alias incompleteness can hide physical overlap. Shared model training and calibration also remain a separate dependence: G26 includes DES plus several low- and high-redshift surveys.

The ZTF application therefore uses different survey photometry from the Pantheon+/DES distance releases, while not being statistically isolated from all of their training/calibration ingredients.

## Revision audit through 20 September 2026

All three source bundles are preserved. The v1-to-v2 TeX change is confined to shortening the abstract. The v2-to-v3 source diff is larger. Source-visible changes include:

- the new title and accepted-paper metadata;
- an explicit 932-object sample count in the abstract;
- replacement of the in-preparation G26 description with citation to arXiv:2606.19429;
- added methodological text and the broken-`W1` step result;
- a new intrinsic-colour-variation subsection and Table 8;
- sign changes for all four `Delta tau_A` entries in the sharp-step Table 7 (magnitudes and quoted significances unchanged);
- a correction in the Appendix A post-processing expression from applying the high-mass shift to both populations to using the corresponding low/high-mass shift.

The latest v3 source and PDF are the baseline. Earlier numerical statements should not be assumed unchanged where the source diff shows a correction.

## Targeted post-27-August primary-literature search

A date-bounded search was frozen at 20 September 2026. Saved arXiv API queries covered exact response phrases, core authors, `progenitor age`, `host age`, and a broader `Type Ia supernova` + `age` query. Saved OpenAlex queries looked for works dated 27 August through 20 September that cite Son (DOI `staf1685`), Wiseman (`stag797`), or the published Chung reply (`stag1513`). OpenAlex returned zero indexed citing works in that interval, and the exact-title/author arXiv queries returned no new explicit response.

Two relevant primary papers appeared in the broader arXiv query, but neither labels itself as a response to the W26/C26 exchange:

- **arXiv:2609.12083v1, 10 September**, Kim et al., *Which Type Ia supernova observables best indicate the ages of their progenitor stars?* It explicitly cites Son et al. 2025 and uses a volume-limited SN sample, a modelled progenitor-age distribution, and local age as a proxy. Its stated focus is whether `x1` or colour indicates progenitor age, not a re-analysis of the W26/C26 residual dispute.
- **arXiv:2609.16972v1, 15 September**, Kelsey, *The UV-to-FIR environments of nearby Type Ia supernovae with the DustPedia Galaxy Catalogue -- I. Local versus global host properties*. It explicitly cites Wiseman et al. 2026 and studies local/global properties for 90 nearby SNe with UV-to-FIR host SED fitting. Its stated focus is environment measurement, not an explicit response or cosmology fit.

Separately, **arXiv:2608.02484v1 (3 August)** is a primary one-page Sah, Rameez & Sarkar response to Ray et al. about an Equatorial/Galactic coordinate mix-up in a directional Pantheon+ argument. It predates 27 August and is not a response to the host-age exchange, but its metadata is preserved because it was absent from the original dossier and is relevant to the separate directional branch.

Search limits: arXiv queries search deposited metadata, not every citation in every full text; author strings can vary; OpenAlex indexing and citation links can lag; publisher-only accepted manuscripts may not yet be indexed; and the interval ended on 20 September. Therefore “no new explicit response located” is a bounded search result, not proof that no such manuscript exists.

## Full-source extraction of the two September papers

Both v1 PDFs, arXiv metadata records, landing pages, and source bundles were preserved. The arXiv records had no later revision at the 20 September freeze. Kelsey also reached MNRAS Advance Access on 18 September as accepted manuscript DOI `10.1093/mnras/stag1765`; its publisher supplement was acquired and passes ZIP integrity.

### Kim et al., arXiv:2609.12083v1

This paper does not measure a distance-residual–age relation. Its new observational comparison is between a host-environment age estimate and SALT2 light-curve parameters, while its progenitor-age distribution is model constructed:

- **Section 2.1, “Constructing the SN Ia progenitor star age distribution,” and Figure 1:** SPAD is the product of an assumed SN Ia delay-time distribution and the cosmic star-formation history. The adopted DTD has late-time slope `s=-1`, `t_prompt=0.3 Gyr`, and polynomial parameter `alpha=5`; the SFH is the Madau & Fragos (2017) cosmic SFR-density form. It is not an object-by-object measured progenitor-age distribution.
- **Section 2.2, “Sample,” and Table 1:** the 1,456-object observable sample contains 890 ZTF DR2 SNe after the Rigault et al. cuts, `z<0.06`, and the paper's additional `c<0.3` cut; and 566 objects from the Nicolas et al. compilation: SNfactory 111 (`z_lim=0.08`), SDSS 167 (`0.20`), PS1 160 (`0.31`), SNLS 102 (`0.60`), and HST 26 with no new volume-limit cut. The authors explicitly say the 26 HST objects are not volume limited. SDSS, PS1, SNLS, and HST came from the original Pantheon catalogue; this is sample reuse, not an independent Pantheon-family data set.
- **Section 2.3, “Empirical descriptions for the light-curve shape and the colour”:** the comparison uses SALT2.4 `x1` and `c`. The ZTF fits used the Taylor et al. retrained SALT2.4 surface. The `x1(z)` description is adopted from Nicolas et al.; the colour distribution is a Gaussian intrinsic-colour term convolved with an exponential dust term, with parameters adopted from the earlier ZTF and Popovic analyses.
- **Section 3.2, “Local age versus stretch and colour,” and Figures 4–6:** the local-age test is a separate 103-SN SDSS-II sample adopted from Rose et al. (2019), selected there using cosmology cuts. Rose et al. measured local `ugriz` photometry in a 1.5 kpc-radius aperture and inferred mass-weighted stellar-population age by Bayesian FSPS SED matching. Kim et al. call local age a proxy and state explicitly that no progenitor-star age has been determined. Their 10,000 LINMIX regressions give an age–`x1` slope of `-0.46 +/- 0.08` and linear correlation coefficient `-0.71`, and an age–`c` slope of `0.00 +/- 0.01` and coefficient `0.08`.
- **Section 3.2 and Figures 5–6:** four local-age bins contain 25, 26, 26, and 26 objects. The authors report that the youngest-bin `x1` distribution differs from each other bin with KS `p<=0.026`, including `5.41e-6` versus the oldest bin, whereas the `c` comparisons have `p>=0.05`.
- **Section 4 and Figure 7:** the empirical mapping draws 100,000 ages from the model SPAD at each redshift, then samples `x1` from KDEs built from the four Rose local-age bins and applies `-3<x1<3`. The reported KS p-values between predicted and observed `x1` are 0.09, 0.35, and 0.28 for the first three bins through `z=0.3`. The source therefore combines modeled ages, previously measured local host ages, and observed SALT2 parameters; those are distinct data levels.

No Data Availability section, machine-readable object list, local-age table, covariance product, posterior samples, or analysis code is included in the v1 paper/source bundle. The source links only the generic LINMIX package. The paper's Figures 4–6 encode the Rose-sample result graphically, but the underlying 103 rows are not supplied with this paper. Its Hubble-residual discussion in **Section 5.2** cites scatter values from earlier studies and proposes a high-`x1` selection; it is not a new residual analysis. Consequently this paper cannot supply an independent measured-age residual check without reacquiring and harmonising external object-level samples and performing a new analysis.

### Kelsey, arXiv:2609.16972v1 / MNRAS `stag1765`

This is a host-environment measurement paper. It does not fit SN light curves or Hubble residuals:

- **Section 2.1, “DustPedia Galaxy Sample”:** DustPedia contains 875 extended galaxies with heliocentric velocity below 3,000 km/s and Herschel observations. The homogeneous imaging spans GALEX, SDSS, 2MASS, WISE, Spitzer, Herschel, and Planck, up to 42 bands and an average of 25 bands per galaxy.
- **Section 2.2, “SN Ia Sample Selection,” and Table 1:** TNS (bulk snapshot last updated 21 May 2026) and Open Supernova Catalogue transients were matched to DustPedia by host name and a deliberately broad 10 arcmin cone, then manually reviewed. All classifications containing “Ia,” including peculiar subtypes, were retained. The sequence is 107 SNe/92 hosts after association review; removal of nine contamination-flagged hosts (12 SNe) and three hosts/SNe with insufficient SED coverage leaves 92/80; the 3 kpc-aperture-size requirement removes two more, leaving **90 SNe in 78 hosts**.
- **Section 3.1, “Global Photometry”:** major contamination flags remove entire galaxies; major artefact and null-coverage flags remove affected bands. Surviving photometry must cover `0.35–3.6 micrometres` and the `60–500 micrometres` dust-emission region.
- **Sections 3.2.1–3.2.2:** images are background subtracted, assigned variance maps, PSF matched to the coarsest retained band, and reprojected. The fiducial local measurement is a fixed **3 kpc radius** aperture. The paper also tests 1 kpc where resolution permits.
- **Section 3.3, “SED Fitting,” and Table 2:** local and global SEDs use CIGALE 2025.1 with a flexible delayed SFH, Bruzual–Charlot stellar populations, Salpeter IMF, fixed solar stellar/gas metallicity, a modified Calzetti attenuation law, and THEMIS dust emission. Redshift is set to zero in CIGALE and redshift-independent DustPedia/HyperLeda distances set the luminosity scale. The grid contains 80,041,500 models; published values are Bayesian SED estimates.
- **Section 4 and Figures 4–5:** studied outputs are stellar mass, `u-r`, sSFR, mass-weighted age `AgeMW`, host-population attenuation `A_V`, and dust mass. These are modeled host/aperture properties. The paper explicitly distinguishes its aperture-integrated attenuation `A_V` from point-source line-of-sight SN extinction (**Section 6**).
- **Section 4 and Figure 4:** the authors report median global-minus-local offsets of `-0.063 mag` in `u-r`, `+0.10 dex` in log sSFR, `-0.013 mag` in `A_V`, and `-312 Myr` in `AgeMW`, each with signed-rank `p<0.01`; the corresponding `AgeMW` object-to-object NMAD is 1,087 Myr.
- **Section 4 and Figure 5:** for mass-weighted age, only 24.4% of objects exceed a one-sigma local/global discrepancy versus 31.7% for an ideal half-normal comparison. The authors note the large per-object age uncertainties and report a formal KS `p=0.007` whose direction differs from the other properties; **Section 6** cautions that broadband-SED ages remain SFH sensitive.
- **Section 5 and Table 3:** the 22 sibling SNe in ten hosts have an in-source table of local/global environment values. This table concerns environment differences within shared hosts, not standardized luminosity residuals.
- **Section 7:** the planned comparison with SN light-curve standardisation properties is deferred to the second paper in the series.

The published supplement materially improves reuse. `cigale_local_3kpc.txt` contains 90 named SNe with local `AgeMW` and its uncertainty plus the other CIGALE outputs; `cigale_global.txt` contains 78 hosts; the archive also contains 3 kpc photometry and 1 kpc photometry/CIGALE tables for the 41-object resolvable subset. These tables permit independent checks of the paper's local-versus-global host-property comparisons. They do **not** contain SALT parameters, standardized magnitudes, distance moduli, residuals, or a residual covariance, so they do not alone permit a measured-age residual check.

A deterministic name-only comparison of the 90-object supplement found 18 matches in the acquired Pantheon+ distance table: SN 1981B, 1990N, 1992G, 1994ae, 1997bq, 1998aq, 2007af, 2007gi, 2007le, 2007on, 2009Y, 2011B, 2011fe, 2011iv, 2012fr, 2016coj, 2019np, and 2021hpr. It found no name match in DES-Dovekie and one match to the linked 944-name ZTF list after joining through the public ZTF IAU aliases: SN 2019np (`ZTF19aacgslb`). This is an identifier audit, not a coordinate/epoch cross-match, and does not establish which objects survive any external analysis mask.

## Directional-paper metadata addendum

The primary Sah, Rameez & Sarkar version of record (DOI `10.1093/mnras/stag844`) was recovered from Oxford's ORA repository and preserved. The currently indexed Ray et al. arXiv record is **arXiv:2607.20570v2, revised 3 August 2026**, titled *Accelerating expansion and isotropic sky-hemisphere consistency in Pantheon+ supernovae: a revised analysis in the dark energy debate*; v1 was submitted 21 July. The v2 abstract explicitly says the earlier sky-hemisphere direction contained a coordinate error, reports accelerating full-sample fits with and without the Son correction, and does not retain the original strong-deceleration claim. This metadata supersedes any unqualified description of v1 in the literature freeze; detailed directional interpretation remains outside this ZTF acquisition note.
