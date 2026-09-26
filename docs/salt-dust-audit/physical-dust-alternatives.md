# Physical alternatives to the extrapolated dust tail

This is a follow-up design, based on primary literature, the pinned public BayeSN implementation, and the local SNANA source. No new population dust fit, BayeSN fit, or cosmology fit was run for this note. Its purpose is to replace an inadmissible **instantaneous passive-screen interpretation** with testable alternatives, without declaring a convenient cutoff to be the true population.

The [implementation audit](snana-implementation.md) identifies negative red-wavelength extinction in some **generated, written mock objects**, and the [independent review](synthesis-review.md) checks that the constructed broadband example is not caused by a negative SALT SED. Those counts are neither observed-SN counts nor final BBC-selected counts. The change in the mean observed distance or cosmology remains unknown. Correcting the old F99 approximation and replacing extrapolation outside a law's physically admissible domain are distinct operations.

## What low extinction ratios do and do not establish

For an unresolved source behind a passive foreground screen with scattered light excluded,

\[
F_\lambda^{\rm obs}(t)=F_\lambda^{\rm int}(t)e^{-\tau_\lambda},\quad
A_\lambda=(2.5/\ln10)\tau_\lambda\ge0.
\]

Nonnegative optical depth is necessary for this interpretation. It does **not** imply a universal lower bound of 2, 1.2, or any other chosen number on \(R_V=A_V/(A_B-A_V)\). It excludes a particular negative-extinction curve, not every other curve sharing its nominal \(R_V\).

There is direct evidence for steep SN reddening: the 14-band UV–NIR comparison of SN 2014J with SN 2011fe favored a Galactic-family \(R_V=1.4\pm0.1\), and also a positive power law \(A_\lambda/A_V=(\lambda/\lambda_V)^{-2.1\pm0.1}\). These are conditional on the intrinsic comparison and extinction parameterization; the different families need not assign identical numerical \(R_V\). They rule out treating \(R_V\ge2\) as a universal observational fact. [Amanullah et al. 2014](https://arxiv.org/abs/1404.2595).

The broader UV–NIR sample found different extinction laws even for two highly reddened SNe, with estimates near 1.4 and 2.8. UV data improved discrimination, but intrinsic UV diversity and broad-filter wavelength shifts still require modeling. This evidence does not measure the frequency of an \(R_V\simeq0.3\) tail in a magnitude-selected cosmology sample. [Amanullah et al. 2015](https://arxiv.org/abs/1504.02101).

Keep four quantities separate:

| Quantity | What it estimates | What would be an unjustified substitution |
|---|---|---|
| Screen \(R_V\) | Ratio of extinction quantities for a specified sightline, SED/filter convention and law | SALT \(\beta-1\), or an integrated-galaxy attenuation slope |
| SALT \(c\), \(\beta\) | Empirical color coordinate and conditional brightness regression | A directly observed nonnegative dust column and grain law |
| Hierarchical dust \(\mu_R,\sigma_R\) | Latent population parameters conditional on intrinsic-color model, selection and priors | The sample mean/variance of noisy independent \(R_V\) estimates |
| Effective attenuation | Transmission after aperture, geometry, scattering and source-population averaging | A unique point-source microscopic extinction law |

For example, with \(c=u+E\), \(m_B-\mu=M+b\,u+R_B E\), constant \(R_B\), and independent \(u,E\), the population least-squares slope is

\[
\beta_{\rm eff}=
\frac{b\,\mathrm{Var}(u)+R_B\,\mathrm{Var}(E)}
{\mathrm{Var}(u)+\mathrm{Var}(E)}.
\]

Intrinsic–dust correlation adds \((b+R_B)\mathrm{Cov}(u,E)\) to the numerator and \(2\mathrm{Cov}(u,E)\) to the denominator. Variable \(R_B\), selection and measurement errors require the more general \(\mathrm{Cov}(c,m_B-\mu)/\mathrm{Var}(c)\). Thus a low effective coefficient can be real without measuring a steep physical screen. Equating it to \(R_V+1\) additionally assumes the same physical B,V definitions and that all relevant color variation is dust. See [the model-assumptions derivation](salt-model-assumptions.md) for the remaining intrinsic/dust degeneracies.

Broadband extinction must be calculated from the extinguished SED,

\[
A_b(t)=-2.5\log_{10}\frac{\int \lambda T_b(\lambda)F_\lambda(t)10^{-0.4A_\lambda}\,d\lambda}
{\int \lambda T_b(\lambda)F_\lambda(t)\,d\lambda}.
\]

For host dust, evaluate \(A\) at rest wavelength, with the observer-frame SED and passband in this integral. Effective-wavelength substitution misses SED, phase, reddening and redshift dependence. Even a time-independent monochromatic screen can give phase-dependent **broadband** color excess; phase dependence alone is not evidence for circumstellar scattering.

## Alternative families worth testing

**A nonnegative power-law screen is the smallest useful diagnostic.** On a declared wavelength interval, use

\[
A_\lambda=A_V(\lambda/\lambda_V)^{-p},\qquad A_V\ge0,\ p>0.
\]

It is an admissible attenuating function, empirically motivated over UV–NIR by the SN 2014J comparison. It is not automatically a validated grain model or proof of a scattering geometry. Fixing monochromatic B,V anchors gives

\[
p=\frac{\ln(1+1/R_V)}{\ln(\lambda_V/\lambda_B)}.
\]

At anchors 4400 and 5495 Å, \(R_V=1.4\) implies \(p=2.4254\), whereas \(R_V=0.4\) implies \(p=5.6371\). The latter can be made positive mathematically but lies far beyond the cited SN 2014J slope and is an **extreme stress function**, not an empirically justified replacement population. In this construction \(A_{8000}/A_V=0.40213\) and 0.12035 respectively. These numbers are algebraic calculations, not fits, and the numerical \(R_V\) is not automatically equal to an F99 parameter or a broadband ratio.

The inspected SNANA already implements the Goobar-shaped approximation
\(A_\lambda/A_V=1-a+a(\lambda/5495\,{\rm Å})^P\)
as option 208 in `src/MWgaldust.c:485–518`. It enforces \(0<a\le1\), \(-2.5\le P\le-0.5\), and 2000–22000 Å. This gives a positive curve and includes the pure power law at \(a=1\). The allowed range has minimum *monochromatic* B,V ratio about 1.346 at those anchors; it cannot preserve the original 0.3–0.4 tail. A runtime rejection outside that range is not evidence that the data rule it out. Also, the inspected host and MW SED paths obtain the same color-law option and parameter list (`genmag_SEDtools.c:2687–2698,2784–2795`): blindly changing a global option would alter MW extinction too. A host-only experiment must verify this separation. [Pinned SNANA source](https://github.com/RickKessler/SNANA/blob/886408a4e171896db5eaa97e735a655f50cec2db/src/MWgaldust.c).

**A nonnegative grain-opacity mixture is a physical second branch.** Use

\[
\tau_\lambda=\sum_s\int n_s(a)C_{{\rm ext},s}(a,\lambda)\,da,
\quad n_s(a)\ge0,\quad C_{{\rm ext},s}\ge0.
\]

Infer extinction ratios from the resulting curve. Grain composition, size distribution and abundance are model assumptions to vary. Silicate/graphite fits to SN 2014J found a small-grain-rich solution near \(R_V\simeq1.7\). [Gao et al. 2015](https://arxiv.org/abs/1507.00417). Nozawa showed that changing size distributions can reproduce steep target curves with nominal ratios 1–2; that calculation fits assumed CCM curves and is not an independent survey measurement of their prevalence. [Nozawa 2016](https://arxiv.org/abs/1608.06689). Public tabulated absorption/scattering efficiencies permit a small quadrature calculation without installing a radiative-transfer framework; [Draine's optical-property tables](https://www.astro.princeton.edu/~draine/dust/dust.diel.html) specify materials, grain sizes and approximations. This audit identified those tables but did not acquire or validate a new grain grid.

**Circumstellar scattering is a different forward model.** A schematic causal model is

\[
F_\lambda^{\rm obs}(t)=T_\lambda F_\lambda^{\rm int}(t)
+\int_0^\infty K_\lambda(\Delta t)F_\lambda^{\rm int}(t-\Delta t)d\Delta t,
\]

with nonnegative transmission/kernel and the radiative-transfer energy budget for the chosen geometry. Delayed photons can increase late-time flux relative to the unscattered source at that same epoch; that is not an instantaneous negative optical depth. Goobar's simulations produce effective ratios around 1.5–2.5 for ordinary grain types in a surrounding geometry. [Goobar 2008](https://arxiv.org/abs/0809.1094). A static power-law fit alone does not establish this mechanism: Brown et al. found SN 2014J time evolution inconsistent with the tested circumstellar contribution and favored interstellar dust. [Brown et al. 2015](https://discovery.ucl.ac.uk/id/eprint/1449756/). Use multi-epoch color/spectral shapes, late light, polarization or echo constraints to distinguish geometries, after modeling broadband evolution from a static screen. This is a later branch, not the minimal repair to the current simulation.

**Changing the label of a Milky Way law is insufficient.** Gordon et al. 2023 calibrate their relation over \(2.3\le R_V\le5.6\); it is a valuable overlap-domain comparison, not validation of tiny \(R_V\). [Gordon et al. 2023](https://arxiv.org/abs/2304.01991). The inspected CCM IR expression is \(A_\lambda/A_V=(0.574-0.527/R_V)x^{1.61}\), also negative for \(R_V<0.91812\). Merely imposing \(R_V>0\), swapping F99 for CCM, or using a lognormal \(R_V\) distribution does not guarantee admissible extinction. Pointwise clipping \(A_\lambda\) to zero or flooring \(R_V\) is an explicitly labeled sensitivity branch, not a derived correction.

## Likelihoods and data that can discriminate

The most direct next comparison is a joint intrinsic-SED/dust model using optical **and** NIR fluxes. In 86 CSP SNe, hierarchical BayeSN analysis explicitly demonstrated that the variance of independently estimated \(R_V\) values can overstate population scatter, especially when residual intrinsic color variation is ignored. Its low/moderate-reddening subset gave \(\mu_R=2.59\pm0.14\), \(\sigma_R=0.62\pm0.16\), conditional on that model/sample. [Thorp & Mandel 2022](https://academic.oup.com/mnras/article/517/2/2360/6713961). These estimates are useful comparison likelihoods, not ready-made priors for a differently selected DES population.

The public [CSP data page](https://csp.obs.carnegiescience.edu/data) provides DR3 photometry, filter/zero-point information, optical spectra and newer NIR spectra. The [RAISIN release](https://github.com/djones1040/RAISIN_DataRelease) includes SNANA-format CSP/RAISIN photometry, calibration inputs, SNooPy models, and published distance/systematic products. It is a concrete small-data starting point rather than an unavailable ideal observation. Refit fluxes for law testing; previously corrected distance tables cannot distinguish arbitrary replacement extinction laws.

The optical–NIR RAISIN comparison uses 42 low-z CSP and 37 higher-z objects and finds a 95% interval \(-1.16<\Delta\mu(R_V)<1.38\). That is a weak bound on population evolution, and the authors explicitly allow remaining cosmological impact. [Thorp et al. 2024](https://arxiv.org/html/2402.18624v2). Do not translate “no significant evolution” into zero dust drift, or transplant its conditional posterior without selection and sample-overlap checks.

The already preserved BayeSN archive is
`sources/updates/2026-09-20-ztf/originals/official-bayesn--08eef9e54188f601506ef9f7f47fc65f9f52a946.tar.gz`.
Inspection of its `bayesn/bayesn_model.py` establishes material limitations:

- It still uses F99 (`615–651`). It changes the intrinsic/dust hierarchy, not the extinction-law family by itself.
- Population fitting truncates \(R_V\) at 1.2 (`263,1032–1043`); the uniform option instead uses 1–6 (`863`). An absent posterior tail below a prior boundary is not evidence excluding that tail.
- G26's `BAYESN.INFO` declares ugriz training and 2800–10800 Å support. M20 declares BVRIYJH training, 3000–18500 Å, and explicitly warns against u/U. Use an appropriate trained model and full passband support; do not extrapolate G26 into J/H to manufacture an NIR check.

The [G26 paper](https://arxiv.org/abs/2606.19429) and [local source audit](../investigations/standardization-dust.md) motivate calibration-aware optical comparisons, but its improved selected DES scatter is not a completed selection-corrected cosmology likelihood. The preserved ZTF lite archive has optical photometry and spectra and can test phase dependence, intrinsic spectral diversity and environment dependence. Its exact published sample mask and BayeSN posterior products remain unresolved; [the source audit](../investigations/ztf-source-update.md) identifies the 932/929/944/933 count differences. A new documented subset is runnable but is not a reproduction of the published fit. Distinct survey observations do not remove shared training/calibration dependence.

## Prioritized minimal experiments

1. **Separate law shape from population replacement using existing inputs.** Extend the existing positive-SED broadband example and response pilot with a positive power-law control. At a fixed reference SED/phase, match **integrated** \(A_B,A_V\), solving for its amplitude and slope; record failures and parameter excursions. Keep the same seeds, intrinsic SED, observer epochs, MW treatment and noise. Integrate the entire used passband, check every supported wavelength for positivity, and inspect early/late phase predictions. A branch holding the old numeric \(A_V,R_V\) is a separate convention test, not matched physical reddening. Test the empirical slope range first; an extreme slope needed to retain \(R_V\simeq0.4\) must retain the “unsupported stress” label. Do not pool alternative-law columns as independent random perturbations.

2. **Fit a small optical subset and predict its held-out red/NIR measurements.** Before a costly full hierarchy, use a fixed, explicitly documented intrinsic template and a small grid over \(A_V,p,t_0,\) shape, with free normalization; retain full likelihoods rather than best-fit dust points. On CSP, compare optical fits' Y/J/H predictions for F99, positive power law and the selected grain basis. Propagate intrinsic/template and shared calibration uncertainty; train or cross-validate intrinsic components without the held-out SN. Repeat with a NIR-trained BayeSN/SNooPy model as a model-dependence test. Optical-only DES/ZTF held-out bands are useful, but weaker, preliminary versions. Report predictive log likelihood, residuals versus wavelength/phase/color, and prior sensitivity rather than only minimized Hubble scatter. Do not use the SALT3.DES5YR wavelength extension as a substitute for validated J/H training.

3. **Refit the population before promoting any law to the bias simulation.** Jointly vary intrinsic color/brightness structure, nonnegative dust-column distribution, law-shape distribution, host dependencies, and their correlations. Allow a residual intrinsic host term; forbidding it can force dust to explain the entire step. Compare a truncated-F99 control against positive-law families, including different tail priors. For a detected sample at fixed host/redshift, the likelihood must account for detection normalization:

   \[
   \mathcal L(\Phi,\kappa)\propto\prod_i
   \frac{\int p(y_i\mid\theta_i,z_i,h_i,\kappa)
   p(\theta_i\mid z_i,h_i,\Phi)d\theta_i}
   {P(\mathrm{det}\mid z_i,h_i,\Phi,\kappa)}.
   \]

   Here \(\kappa\) represents shared calibration and \(\Phi\) population hyperparameters; the detection numerator must also match the survey's actual conditioning if selection involves additional data. Rate/count inference adds its own likelihood. Model actual detection, spectroscopic/classification and analysis cuts rather than labeling a fitted selected distribution the parent population. Red-tail exclusions deserve their own sensitivity branch because they remove the objects most informative about dust.

4. **Only then measure the distance consequence.** Regenerate matched-seed surveys for each empirically competitive population/law, refit light curves, repeat classification/quality selection and BBC, and propagate paired distance differences with shared calibration/training uncertainty. Compare the true/recovered distance bias, its redshift and host dependence, and cosmology with the intercept marginalized. Distinguish a frozen-population law swap from a newly fitted population. Small changes in a bounded stress test validate only those objects/branches; neither broad physical positivity nor a lower scatter proves end-to-end bias closure.

The immediate deliverable should be a comparison of admissible models and falsifiable held-out flux predictions. A justified replacement for the mock tail requires those predictions and selection-aware population inference; the present evidence does not select a unique floor or establish a numerical global bias bound.
