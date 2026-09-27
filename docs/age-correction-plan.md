# Testing whether an additional supernova age correction is needed

**Study design, 27 September 2026.** The [execution results](age-correction-results.md) now report the feasible observational, recovery and transport tests. The full selection-aware physical inference remains incomplete; no new cosmological result is claimed. Availability updates below are findings from execution, not retrospectively blinded design choices. The [evidence review](age-correction-evidence.md) gives the case on each side and the existing calculations. This plan makes the age question a focused part of the broader [unified cosmology programme](experimental-plan.md).

The primary objective is to measure **the redshift-dependent distance bias left after standardization**, and determine whether an age-aware model reduces it without correcting the same effect twice. Identifying the ultimate explosion physics is a separate, harder objective. We may resolve the first question within a bounded population domain without uniquely resolving the second.

The later [physical extension](physical-program-results.md) executes independent spectroscopic/aperture validation, high-redshift and infrared data recovery, flexible physical stellar fits, and photon-level survey interventions with refitted corrections. Its final section identifies the remaining measured and computational gates. It supersedes earlier availability limits where new public observations were recovered; it does not retrospectively alter this design or claim complete survey-wide inference.

## 1. Define the alternatives and the quantity to measure

| Hypothesis | Testable prediction |
|---|---|
| Existing treatment is sufficient | Age adds no material improvement in held-out distance prediction, and the standard analysis recovers distances across admissible age/dust population changes. Physical age dependence may still exist. |
| Existing treatment is incomplete | A reproducible additional age-related term improves held-out predictions and removes a redshift-dependent recovery bias under independently supported population models. |
| The proposed extra template double counts | Adding the full historical-age template creates overcorrection when compared with a jointly refitted model and known injected truths. |
| Current observations do not identify the correction | Different latent age/dust/population models fit the observations but predict materially different distance evolution. Report the remaining range and the observations needed to narrow it. |

Let `S=1` denote inclusion in the analysis, and let η describe the true population and observing process. Define

`B_std(z; η) = E[mu_hat_std − mu_true | z, S=1, η] − reference_mean`.

Use a declared weighted low-redshift reference sample for the additive zero point. This removes an arbitrary absolute-magnitude offset; it does not remove redshift evolution. Positive B means inferred distances are too large. A validated residual correction would subtract B. In a joint generative analysis, the corresponding effect belongs in the likelihood instead of being subtracted again.

Report the whole bias function, its covariance and its projection onto the cosmological quantities of interest. For paired simulated populations with an age mechanism enabled/disabled, the difference in B isolates that mechanism **within the declared generative model**. It is not a model-free causal measurement.

Before inspecting new outcomes, freeze the sample, redshift range, correction definitions, candidate models, holdouts, power targets and decision rules. Existing inspected data cannot retrospectively become a blind test.

## 2. Data required

| Input | Present starting point | What is still needed and why |
|---|---|---|
| Exact disputed SN sample | C25 G11/R19 summary tables; frozen Pantheon+ distance/covariance release; our 196-object crosswalk | Author age choice for duplicate IDs, precise error definitions, nuisance coefficients and residual conventions. Execution recovered the original Gupta quality columns and the 175/70 Park sample counts. This separates changed samples from changed estimators. |
| Original age inference | Published age summaries; original R19 code; public host photometry | C25 full joint age/mass/dust/metallicity/SFH posterior samples **and original priors**, likelihood definition, synthesis-library versions, apertures, calibration and photometric covariance. Posterior samples alone are insufficient to change the population prior safely. |
| Redshift adjustments and mocks | C25 adjusted/original G11 residuals; published C26 mock description | Young-host fit coefficients and full covariance with the age regression; exact C26 age-PDF arrays, randomization/noise choices and conversion code; exact Park regression settings. The original quality crosswalk is now recovered. These are exact-reproduction blockers, not reasons to substitute newer inputs silently. |
| Independent low-redshift check | DES low-z photometry/distances; ZTF DR2 photometry and host mass/colour; Foundation/CSP optical–NIR assets | W26 morphology IDs/labels and weights; independent age likelihoods; final ZTF per-analysis masks and configurations/posteriors. The public master list does not uniquely specify the published samples. |
| Host-to-progenitor mapping | Childress/Wiseman model code; TITAN paper; DustPedia environmental supplement | Object-level SFH likelihoods, selection and host associations, DTD uncertainty, local/global mapping and held-out age validation. Execution located an older public 8,610-host TITAN summary catalogue and used it. The final 6,983-object sample and joint SFH files remain unreproduced; summaries are not replacements for their likelihoods. |
| Complete survey treatment | Released SALT/photometry, covariances, simulations and first-party research code | Pinned training sets and overlaps, selection/follow-up/classifier models, rejected events/epochs, signed host/bias correction components, correction-training covariance and complete simulation configurations. A covariance matrix alone cannot reconstruct a signed correction. |
| Current comparison releases | Pantheon+, DES-SN5YR, Dovekie and public host-mass reconstruction | Audit Union3.1 v2 and its linked products against the earlier sensitivity inputs. Each release must keep its own calibration, light-curve model and covariance. |
| Cosmological cross-checks | DESI BAO likelihood and existing geometry analysis | Separately validated CMB/anchor likelihoods only for the combinations actually reported; explicit shared calibration and sample ancestry. |

The [source manifest](../provenance/age-correction-literature.json) records acquisition evidence. The [original study guide](../studies/README.md) explains how to restore historical layouts and pinned inputs. First-party analysis code stays versioned; downloaded papers, author dependencies, input arrays and regenerable tables remain outside the tracked manuscript. No author contact has been made for this plan.

Public entry points are the [C25 journal supplement](https://academic.oup.com/mnras/article/538/4/3340/8098234#supplementary-data), [Pantheon+ release used in the exchange](https://github.com/PantheonPlusSH0ES/DataRelease/tree/7fc680548d1ea9fe6ee983a89d8b5635bb5784a8), [DES-SN5YR release](https://github.com/des-science/DES-SN5YR), [ZTF DR2 archive](https://ztfcosmo.in2p3.fr/download/data), [Kelsey/DustPedia supplement](https://doi.org/10.1093/mnras/stag1765), [BayeSN model](https://github.com/bayesn/bayesn), and [galaxy-population simulation code](https://github.com/wisemanp/des_sn_hosts). Pin actual files and commits before analysis; these entry points do not imply that every required posterior or selection product is public.

## 3. Experiment A: put both sides on exactly the same footing

Build one object-level correction ledger with immutable source row IDs. Distinguish original author residuals, literal released corrected distances, a reconstructed Tripp-only estimate, and each separately reversed exported component. Reversing a column while retaining fitted coefficients and selected objects is an accounting diagnostic, not a new uncorrected survey measurement.

Cross the following choices without selecting only the favourable combinations: original versus updated host measurements; common versus paper-specific samples; full versus cumulative redshift ranges; tabulated versus covariance-based uncertainties; Gaussian age summaries versus joint host likelihoods. Estimate differences with paired resampling of physical SNe, respecting shared calibration and peculiar-velocity modes. Catalogue entries for the same object are not independent observations.

Reproduce the published summaries first, then fit all standardization coefficients jointly with an extra age term. For C25/Park's young-host redshift adjustment, propagate its shared fitted uncertainty and compare with a joint model and a cross-fitted estimate trained on other objects. Do not treat the adjusted residuals as independent measurements with unchanged error bars.

Separate age variation **within redshift and host-property strata** from changes in the population mean **between redshifts**. Require sufficient overlapping age support; publish where matching or weighting fails. Fit cosmology/redshift dependence jointly, rather than removing an arbitrary brightness trend first and then interpreting the remaining slope as the total effect. Repeat with a free-distance-bin diagnostic, explicitly acknowledging the evolution modes it makes unidentifiable.

**Deliverable:** a reconciled table explaining how each reported slope changes when one choice changes, with shared uncertainties and posterior/likelihood diagnostics. This should determine which disagreements survive harmonization. It cannot alone certify high-redshift transport.

## 4. Experiment B: measure what the real analysis absorbs

This is the central missing experiment. Simulate the same eligible events and observing noise in paired realizations, introduce age effects **before detection and fitting**, and rerun every downstream step affected by the intervention. Recompute selection after injection; keep excluded objects and fit failures in the accounting.

| Simulated truth | Standard treatment | Standard treatment + imposed historical template | Joint age/dust treatment |
|---|---|---|---|
| Dust variation; no extra intrinsic age effect | Tests false age detection and baseline distance recovery | Tests overcorrection | Tests whether the extra model invents an age term |
| Intrinsic age effect; dust correctly specified | Measures how much width/colour/host/BBC fitting absorbs | Tests whether the full template adds too much or too little | Tests recovery of the injected residual contribution |
| Age and dust correlated, both varying | Tests partial absorption | Tests double counting of correlated effects | Tests separation and uncertainty coverage |
| Misspecified SFH, metallicity/channel mixture, dust law or selection | Tests robustness outside the training assumptions | Tests accidental agreement | Tests failure detection and identification limits |

Within each row vary age-effect amplitude, nonlinearity, within/between-redshift dependence and colour/phase dependence. Include the proposed approximately −0.030 mag/Gyr scale as one scenario, zero as another, and smaller effects; do not assume every physical age effect is a grey magnitude shift. Inject both grey and supported spectral/temporal changes. Hold the generation model separate from the inference model so that successful recovery is not guaranteed by identical assumptions.

Run two explicitly different branches: a frozen-trained-model branch to locate where information is absorbed, and a full refit/retraining branch to assess the actual measurement. Refit α, β, host terms and population parameters; regenerate BBC-style bias corrections and selection normalization. Retrain the light-curve surface where the hypothesized population change affects its training. A joint flux likelihood that already models selection must not also apply an independent correction for that same selection.

Measure `B_std(z)`, corrected B, coverage, false positive rate, and failures as functions of age, mass, colour, redshift and survey. Compare the actual fitted-distance projection with the existing bin-centering surrogate. A simulation showing that bin centering attenuates a slope is not enough.

**Deliverable:** an empirically calibrated response from injected population changes to recovered distances, including off-diagonal uncertainty. The standard pipeline must be allowed to fail; the proposed age model must face the same tests.

## 5. Experiment C: use independent observations to separate mechanisms

Construct a low-redshift sample spanning young/old environments **at overlapping mass, metallicity, colour and redshift**, with coverage of passive and star-forming hosts. Combine local spectroscopy, UV–optical–IR host measurements and multi-epoch optical–near-infrared SN photometry. Use physical apertures and PSF/deblending models consistently; central-galaxy spectra cannot automatically stand in for the explosion site. Infer local and global SFHs separately, including dust/metallicity uncertainty.

Fit dust-only, intrinsic-population-only and mixed generative models to a training subset; test their predictions in a held-out survey or observing programme. Hold out whole physical SNe and track shared hosts, fields, model training and calibration. Evaluate both colour/phase structure and brightness contrasts. Optical–IR contrasts suppress the distance contribution and help constrain chromatic effects; they do not identify arbitrary grey luminosity evolution.

Use the following complementary comparisons:

- At similar mass and redshift, test whether age adds predictive information after jointly fitting dust and light-curve properties.
- In low-reddening objects, test whether environmental differences persist, while modelling the selection introduced by this restriction.
- For sibling SNe in the same host, use paired distances to remove much of the distance/host nuisance; require local-age contrast and account for correlated calibration and dust. Siblings without age contrast cannot estimate an age slope.
- Test environment-dependent spectral or light-curve-shape predictions on data not used to infer the environmental correction.
- Use field/instrument/redshift-shuffled controls designed to preserve the relevant selection and covariance. A naïve permutation across an evolving population is not a valid null.

ZTF provides a useful starting population. Execution joined 401 SNe to older global TITAN age summaries; independent local-age likelihoods, documented mask/calibration gaps and full host posterior covariance remain unresolved. Foundation/CSP and the DustPedia hosts supply complementary observables, subject to overlap and follow-up selection. Neither an alternative fitter nor a larger compilation creates independent data by itself.

**Deliverable:** held-out predictive comparisons that constrain which mixtures of intrinsic population and dust effects can explain the observations. An “intrinsic” component should not be renamed “age” unless age actually discriminates it from the alternatives.

## 6. Experiment D: test the host–progenitor mapping and high-redshift transfer

Infer progenitor-delay distributions from local/global SFHs and an uncertain DTD, keeping their covariance with the brightness slope. Compare mappings on the *same observed hosts* before applying them to redshift evolution. Validate synthetic host observations through the age-inference procedure, including priors, to measure age bias and regression dilution rather than assuming posterior medians are unbiased ages.

The slope–age-contrast compensation argument is exact only under restricted mappings. If `tau = a + k A_host` everywhere, changing units gives `b_tau = b_host/k`, so `b_tau × Delta tau = b_host × Delta A_host`. A nonlinear, scattered, selected or redshift-dependent mapping need not preserve that identity. Test both the simple invariant case and departures from it; never change the DTD in the age contrast while silently fixing an incompatible slope.

Estimate the age distribution for the selected population of each survey. Cosmic star-formation history convolved with a DTD describes a parent population; it cannot automatically replace a survey's detected hosts. Infer or validate detection, follow-up, host-measurement and classification selection. Integrate the brightness model over the resulting distribution: the mean effect need not equal the effect evaluated at a median age.

A second route is to match the distribution of measured environmental likelihoods across redshift and construct an age-restricted Hubble diagram. It must use observed high-redshift host information, an unconstrained brightness intercept and an explicit selection model. Matching only low-z ages to an assumed high-z mean would not supply the required independent check. Report loss of age support and the bias–variance cost of restriction.

**Deliverable:** a redshift-dependent residual correction or bound with uncertainty from host inference, DTD, slope, selection and calibration. Shared parameter errors produce correlated distance shifts, including cross-covariance with existing corrections.

## 7. Precision, model comparison and decision rules

Choose a target distance-bias tolerance from the intended uncertainty in q(z), w0 and wa using the existing likelihood response. State both a magnitude tolerance over the redshift interval and a permitted cosmological shift. Do not choose the tolerance after seeing the inferred correction.

Plan sample size by injection/recovery with the actual age spread, covariance, selection and calibration floor. As a rough design diagnostic only, a linear Gaussian model has slope information proportional to the covariance-weighted age variation left after projecting out nuisance predictors. If ages and mass/redshift have almost no independent variation, more objects of the same kind will not resolve the question. Independent local-age accuracy and optical–IR coverage may be more valuable than raw sample size.

Report calibrated interval coverage and detection power for both zero and scientifically relevant injected effects. Use robust residual models, alternative age likelihoods and justified population priors. Compare joint likelihoods/predictive scores, not raw χ² values from different samples or different covariance normalizations. Repeated redshift cuts form a sensitivity family, not multiple independent confirmations.

Cross-fitting nuisance relationships can reduce reuse of the same data to construct and test a correction. A partial-linear or flexible regression is a useful diagnostic only where overlap and sample size support it. Conditioning on width/colour/host properties targets the **incremental** predictive age effect; it does not recover a total causal age effect when those quantities mediate age. No machine-learning method removes missing physical observables or selection ignorance.

| Outcome | Evidence required |
|---|---|
| Existing treatment adequate within the tested domain | Upper bound on residual distance evolution lies below the declared tolerance, simulations show sensitivity to larger effects, and held-out observations support the population model. Failure to reject zero alone is insufficient. |
| Additional correction supported | Incremental predictions replicate, recovery tests distinguish the competing explanations, and the selected-population correction improves distances with calibrated uncertainty. Its size may differ from the proposed template. |
| Full proposed correction overcounts | Paired simulations and independent observations favour a smaller residual adjustment or a replacement model; appending the full template creates systematic offsets. |
| Unresolved | Competing models remain observationally equivalent or essential input/selection likelihoods are absent. Publish conditional bounds and the observable that would separate them. |

Only then propagate admitted models into SN-only, SN+BAO and separately SN+BAO+CMB likelihoods. Keep cosmology out of the model-selection target: choosing an age correction because it agrees with a desired BAO/CMB contour would circularly validate it. Update the cosmic time–redshift relation inside age/SFH inference wherever relevant, rather than imposing a clock from the preferred final cosmology.

## 8. Execution order and stopping points

1. **Public-data reconciliation first:** restore the existing age calculations; freeze current releases and the new Union3.1 version; complete Experiment A's common-sample accounting and list exact missing author products. Existing results provide a strong starting point; reimplementing every historical fit from scratch is unnecessary.
2. **Small paired recovery test next:** implement Experiment B using a supported survey/model subset, including null and nonzero injections and the real cosmology fit. This directly tests suppression versus genuine correction. A simplified pilot must remain labelled until complete selection and retraining are available.
3. **Independent population validation:** obtain missing public/author data where available, then run Experiment C and the mapping tests in D. If observations cannot distinguish the models, design new spectroscopy/IR follow-up from forecast information gain, not an invented fixed sample-size target.
4. **Survey-wide inference last:** regenerate the complete corrections or fit the joint selected-data model; validate transport and coverage; propagate the remaining uncertainty to cosmology.

Every stage should retain inputs, code identity, source/release version, parameter definitions and validation results. Original research remains intact. Exact author reproduction, approximate reconstruction, new experiment and blocked comparison must remain visibly distinct. This plan authorizes no claim that a correction, its absence, or a new cosmology has already been established.
