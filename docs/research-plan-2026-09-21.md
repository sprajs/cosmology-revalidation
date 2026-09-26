# A prioritized research programme to test supernova acceleration and its corrections

Prepared 21 September 2026. This is a research design grounded in the existing workspace and a targeted refresh of primary literature. It does not report a newly executed cosmology analysis. Existing observations and many results have already been inspected, so this is a transparent registered reanalysis plan, not a claim of pristine preregistration or blinding.

**The decisive task is to measure which corrections the observations support, propagate those corrections through the measurement and selection process, and then test the expansion history. Agreement between fitted cosmologies is a result to explain, not a target used to select corrections.**

The programme covers both meanings of “original supernova results”: the 1998/1999 discovery analyses and the modern host-age challenge to supernova cosmology. It keeps historical reproducibility, correction validity, present acceleration, earlier acceleration, and the cosmological-constant hypothesis as distinct questions.

Companion files: [structured experiment register](research-plan-2026-09-21/experiment-register.json) and [local evidence manifest](research-plan-2026-09-21/source-manifest.json). The existing [phase-1 report](final-report.md), [phase-2 plan](phase2/plan.md), and [SALT/dust audit](salt-dust-audit/plan.md) remain their own scientific records.

## 1. The questions and what would answer them

The historical studies inferred acceleration from standardized supernova brightness versus redshift, under specified luminosity, dust, population, and cosmological assumptions. Reproducing their calculations answers whether those published inferences follow from those inputs. Establishing robustness to newly measured systematics is a further test. See [Riess et al. 1998](https://arxiv.org/abs/astro-ph/9805201) and [Perlmutter et al. 1999](https://arxiv.org/abs/astro-ph/9812133).

For a homogeneous, isotropic expansion history,

$$q(z)=-1+(1+z)\frac{d\log H(z)}{dz}.$$

Negative q means acceleration. Supernovae measure an integrated distance relation, so resolving q requires explicit geometry, photon-propagation, and smoothness assumptions. A point estimate at z=0 is particularly sensitive to extrapolation and local velocities.

Use five separate outcomes:

| Question | Primary quantity | What a defensible result would mean |
|---|---|---|
| Do the original numerical analyses reproduce? | Original likelihoods, nuisance fits, contours and sample membership | Their conclusions follow, or fail to follow, from their documented implementation and assumptions. |
| Does a correction improve the description of observations? | Held-out flux and conditional luminosity predictions, with calibrated uncertainty | Evidence for the specified correction model; not automatically proof of its physical interpretation. |
| Is q negative on average over a resolved recent interval? | Primary \(\bar q_{0.05:0.30}=(0.25)^{-1}\int_{0.05}^{0.30}q(z)\,dz\) | An average recent-acceleration statement, with an explicitly reported resolution kernel; it does not establish the same sign at every redshift. |
| Is it accelerating at present, or did acceleration occur anywhere in the observed range? | q(0); separately the constrained family q(z)≥0 over a prespecified range | Two different claims. A positive q(0) can coexist with earlier acceleration. |
| Is a cosmological constant adequate? | Comparison of ΛCDM with specified alternatives | A dark-energy-model test, which does not by itself determine the sign of q(0). |

The interval 0.05–0.30 is a proposed design choice to reduce reliance on the very nearest objects; use 0.02–0.20 and 0.10–0.40 as declared sensitivities. Simulations must establish the effective resolution. If the data cannot resolve the proposed primary interval, report that limitation and publish an amended estimand before examining the corresponding real-data result.

In a spatially flat matter+CPL dark-energy model, neglecting radiation at the present epoch,

$$w(z)=w_0+w_a z/(1+z),\qquad q_0=\tfrac12+\tfrac32w_0(1-\Omega_m).$$

Thus departure from \((w_0,w_a)=(-1,0)\) and evidence for present deceleration are different statistics. The [DESI DR2 analysis](https://arxiv.org/abs/2503.14738) provides an important external comparison, with its early-universe assumptions retained explicitly.

## 2. What is already done, and what is actually unresolved

The local project is well beyond a literature collection. Its phase-2 literature map contains 61 pinned sources, but explicitly distinguishes focused method checks from carried summaries. This plan uses that map, the completed local investigations, and current files; it does not claim to have exhaustively reviewed every supernova paper.

| Existing work | Evidence available now | Consequence for priorities |
|---|---|---|
| Released-distance cosmology | The local Pantheon+ plus free-ruler DESI BAO CPL baseline gives q0≈−0.436±0.081. A fixed, specified age/median template moves this to about +0.063±0.074. | The ability of an imposed correction to move q is already demonstrated. Repeating that result has low information value. |
| Modern matched host ages | On the same modern corrected objects the additional slope is much weaker than historical estimates; the simple latent-age profile interval is approximately [−0.032,+0.010] mag/Gyr. | Neither a universal extra correction nor absence of a consequential effect has been established. Joint age errors and population transport matter. |
| Host-to-progenitor mapping | Mean/median choice, delay-time law, clock, nonlinear mappings and sample weighting change the transported correction. | Estimate the selected population distribution and its uncertainty, not just a universal slope times a cosmic median. |
| Original DES versus Dovekie | Released distances/covariances give different conditional cosmological fits. Calibration, dust implementation and covariance changes interact. | Treat versions as branches of the same observations. Attribute changes with controlled swaps and full reruns. |
| Observation inventory | Public DES SMP, low-z and DIFFIMG tables have been audited below fitted distances. Original Hubble membership is 1,829 objects. | Use existing integrity work. The remaining uncertainty is not principally missing basic table parsing. |
| Calibrated-flux refits | Latest inspected recovered run fits 1,828/1,829 objects; 66 have \(|\Delta m_B|>0.01\) mag. The high nominal-Ia-probability subset still has 22 such cases. | Most fits are close, but the reproduction discrepancy cannot be dismissed using only the median or the largest low-weight outlier. |
| Distance algebra | Published distance algebra closes to ≈0.000266 mag maximum using unrounded host mass. | This verifies algebra; it does not regenerate BBC corrections or their uncertainty. |
| Measurement covariance | A local output-only SNANA patch exports coherent full-precision Hessian covariance. | This is useful for summary likelihoods. MINOS errors must remain separate; joint survey systematics still need modelling. |
| Selection assets | Trigger, host-z and population assets exist. A new nominal forward pilot records 26,518 generated attempts and 2,678 written light curves. | Generated failures are available for that pilot. The final measured-fit, classification, BBC and membership chain is incomplete. |
| Independent hierarchy | Synthetic likelihood integration is independently checked; 30 fixed-distance MLE recoveries and four seeds per free-q stress exist. | Numerical closure is encouraging. It is not sufficient power analysis, posterior SBC, or validation of real DES selection. |
| SALT/dust response audit | The initial 12-object pilot required fixes to its own magnitude-offset, primary-calibration and MW-reddening conventions. A corrected pilot appeared during preparation; its 20 input/output manifest hashes were checked and match. | Use the corrected, identified pilot as a conditional response calculation. The broader planned experiment and full scientific audit remain separate. These diagnostic errors are not discoveries of corresponding DES errors. |
| Noise and repeated records | Public flux records expose repeated groups and off-signal diagnostics. SMP and DIFFIMG behave differently. | Test covariance and extraction at the appropriate data layer; do not turn a DIFFIMG anomaly into an SMP distance correction. |
| Independent numerical audit | Separate distance and likelihood implementations agree at tested points. | Useful implementation verification on shared data, not an independent astrophysical confirmation. |
| Current external challenge | The Son–Wiseman–Chung exchange, newer ZTF environment results, calibration revisions and DESI comparisons are already mapped. | The next gain comes from tests distinguishing the explanations, not counting papers on each side. |

Sources for the local numbers: [final report](final-report.md), [DES refit report](phase2/official-reproduction.md), [newer recovered diagnostics](../phase2/official/diagnostics/recovered_refit_summary.json), [weighted diagnostic](../phase2/official/diagnostics/recovered_weighted_gate.json), [data audit](phase2/data-audit.md), [selection audit](phase2/literature-selection.md), [synthetic validation](phase2/independent-engine-validation.md), and [dust-pilot amendment](salt-dust-audit/plan.md).

**Status caution:** these workstreams were changing while this plan was prepared. Some reports describe earlier runs. File existence alone is not completion: require a successful log, complete output, manifest, and an interpretation matching the executed source. The companion manifest records inspected bytes; E00 below is the gate for later execution.

The latest weighted refit diagnostic holds published standardization/BBC quantities fixed and uses diagonal weights. Its linear Ωm shifts are not new cosmological fits. Likewise, current monochromatic F99 implementation tests verify an extinction formula; they do not determine the broadband, refitted, selected-sample cosmological correction.

The targeted literature refresh supports the following map of the live dispute. “Current” means the cited release/revision checked for this plan, not an assertion of exhaustive coverage.

| Research direction | Primary source and checked version | Claim or issue to test here |
|---|---|---|
| Age-dependent luminosity and transported correction | [Chung I](https://arxiv.org/abs/2411.05299), v2; [Son II](https://arxiv.org/abs/2510.13121), v1 | Separate the local empirical age relation, its extrapolation to selected progenitors, and the resulting cosmology. E06/E08/E13 test those links separately. |
| Modern corrections and population-based critique | [Wiseman et al.](https://arxiv.org/abs/2601.13785v2), 8 May 2026 | Test identical-object residuals and selected delay distributions, with correction state and uncertainty explicit. |
| Counterargument about residual dilution and age mapping | [Published Chung reply](https://doi.org/10.1093/mnras/stag1513); [arXiv v1](https://arxiv.org/abs/2605.21586) | The published version contains additional mock material. Test flexible-baseline absorption through the actual selection/fitting pipeline; bin centering is not automatically a physical correction. |
| Error accounting in earlier age claims | [Rose et al.](https://arxiv.org/abs/2002.12382v2) | Error models, influential observations and sample extrapolation must be tested, without assuming an older critique settles every later estimator. |
| Official DES measurement and selection methodology | [Vincenzi et al.](https://arxiv.org/abs/2401.02945v2); [DES photometry](https://arxiv.org/abs/2406.05046) | E01–E04 reconstruct the appropriate data layer and distinguish calibrated flux from detector-image remeasurement. |
| Calibration, training and extinction revision | [DES-Dovekie](https://arxiv.org/abs/2511.07517v3), 27 March 2026 | Attribute interacting changes instead of assigning the total shift to one mechanism. |
| Environmental signal under flexible dust modelling | [ZTF BayeSN](https://arxiv.org/abs/2605.06799v3), 17 September 2026 | Test whether an intrinsic environmental component predicts held-out data; then test its population evolution rather than assuming it. |
| Better age indicators and host environments | [Kim et al.](https://arxiv.org/abs/2609.12083); [Kelsey et al.](https://arxiv.org/abs/2609.16972) | Improve host/population measurement. These are not by themselves new age-corrected acceleration likelihoods. |
| External geometry and dark-energy inference | [DESI DR2](https://arxiv.org/abs/2503.14738v3), 9 October 2025 | Separate free-ruler geometry, early-universe calibration and full CMB information. |
| Orthogonal wavelength and redshift leverage | [RAISIN](https://arxiv.org/abs/2201.07801v2); [high-z HST test](https://arxiv.org/abs/astro-ph/0402512v2) | These extend the local core map. Their executable data products, covariance and overlap still need an acquisition audit. |

The local [literature/selection review](phase2/literature-selection.md) also covers global versus evolving latent populations, inter-compilation offsets, photometric contamination, host-redshift selection and their counterarguments. Those methodological disputes are represented in Gpop/Gcoef, E03, E04 and E10 rather than treated as votes for a cosmological conclusion.

## 3. Priority order and dependencies

P0 closes errors capable of invalidating every later comparison. P1 measures the corrections and their identifiability. P2 performs statistical discrimination and cosmological inference. P3 addresses broader historical or physical alternatives when the necessary input is ready. Within a tier, prioritize expected information per cost; several cards can run in parallel.

| ID | Priority | Experiment | What it decides | Readiness |
|---|---|---|---|---|
| E00 | P0 | Freeze evidence, aliases and validation partitions | Which inputs and objects the comparison actually uses | Local; immediate |
| E01 | P0 | Close DES fitting and covariance discrepancies | Whether a changed result is a physical model or an implementation difference | Local; immediate, with an exact-runtime gap |
| E02 | P0 | Repeated-record and noise/extraction tests | Whether effective flux information and uncertainty are correctly represented | Local diagnostics; images may be needed |
| E03 | P0 | Complete generated-to-selected simulation chain | Whether selection corrections remain valid under alternatives | Partly implemented; final cuts incomplete |
| E04 | P0/P1 | Signed calibration, dust and covariance audit | Which known implementation/calibration changes can account for differences | Local; corrected 12-object pilot verified, broader audit pending |
| E05 | P1 | Constrained reconciliation and minimum required correction | Whether independently allowed corrections can bridge each cosmological distance gap | Local screening now; inference after P0 |
| E06 | P1 | Matched-age residual and confounding experiment | Whether age predicts residual luminosity beyond competing measured properties | Summary data now; joint PDFs limited |
| E07 | P1 | Dust versus intrinsic-colour/host luminosity experiment | Which wavelength/phase predictions discriminate the models | DES local; NIR extension requires acquisition |
| E08 | P1 | Selected host-to-progenitor transport experiment | Whether a local relation predicts the needed redshift evolution | Partial; exact W22 and richer ages missing |
| E09 | P1 | Adversarial recovery, SBC and power | Whether the planned decision rule works under correct and incorrect models | Existing engine; survey extension required |
| E10 | P2 | Common-observation hierarchical model tournament | Which corrections predict data best without cosmology choosing them | Depends on E01–E04 and E09 |
| E11 | P2 | Survey, wavelength, field and redshift transfer | Whether predictive gains survive changes in observing conditions | DES available; external products must be checked |
| E12 | P2 | Acceleration reconstruction and no-acceleration test | Which acceleration statements survive measured uncertainty | Depends on admitted models and validation |
| E13 | P2 | Separate-probe prediction, then joint inference | Whether one model can explain SN and non-SN observations | BAO local; exact CMB branch has gaps |
| E14 | P2, parallel | Original 1998/1999 and high-z historical reconstruction | Whether original claims reproduce and survive justified modern changes | Published papers available; exact legacy assets unverified |
| E15 | P3 | Anisotropy, velocities, curvature and propagation alternatives | Whether a specific broader assumption explains residual structure | Bounded local groundwork; calibrated tests needed |
| E16 | P3, escalates on E02 | Targeted pixel and training-data reconstruction | Whether unresolved effects originate below released photometry or in shared SED training | Archive and training completeness must be established |
| E17 | P3 | Design the smallest missing measurement | Which additional data would actually break the remaining degeneracy | Follows E05/E09/E11 information analysis |

The critical path is E00 → E01/E02/E03/E04 → E09 → E10/E11 → E12/E13. E05, E06, E08 and E14 have useful work before that path completes. E16 moves to P0 if E02 demonstrates an unresolved effect large enough to threaten the primary inference. Literature and public-data checks continue only where they change an experiment.

## 4. Model comparison must have two axes

An age or dust prescription is a model of the supernova population and measurement process. ΛCDM or a no-acceleration history is a model of the expansion. Cross the two axes fairly; do not give one expansion history much more flexible corrections than another.

**Observation/population axis**

| Family | Specification | Key discriminator |
|---|---|---|
| G0 | Tripp relation with latent stretch/total colour, host term, intrinsic scatter, and measured covariance | Baseline conditional magnitude–colour–stretch patterns |
| Gpop | Redshift/host-dependent stretch and colour distributions; luminosity coefficients initially fixed | Population shrinkage and selection versus luminosity evolution |
| Gcoef | Prespecified drift of standardization coefficients, with Gpop retained | Whether the colour/shape–brightness relation changes |
| D | Intrinsic colour plus nonnegative reddening, separate intrinsic and extinction coefficients | Multiband/phase prediction and colour-dependent scatter |
| DH | Host-dependent reddening scale and/or extinction law | Environmental pattern explained through dust |
| DHM | DH plus an intrinsic environmental luminosity component | Whether dust alone predicts held-out host contrasts |
| A | An additional latent-age contribution, conditional on dust, metallicity, host information and selection | Within-redshift age contrasts and cross-population transport |
| AD | Joint age, intrinsic colour and dust model with regularization | Whether an age effect persists when competing mechanisms coexist |
| F | Prespecified flexible empirical colour/host relation or small mixture | Tests whether physical labels are overinterpreting an inflexible functional form |
| L | Smooth luminosity-evolution sensitivity, anchored at a reference redshift | Identifies how much acceleration information depends on restricting evolution |

These are families, not a license for an unrestricted search through hundreds of models. Freeze a small main roster, fixed basis dimensions, prior predictive envelopes and limited sensitivities. Treat later inventions as exploratory and test them in an outer validation loop.

Start from the [existing hierarchical specification](phase2/generative-model.md), but do not silently treat its effective SALT dust coefficients as physical extinction-law parameters. A physically interpreted dust branch must apply reddening to the time-dependent SED through actual bandpasses and then refit.

**Expansion axis**

1. Flat ΛCDM as a reference.
2. Flat constant-w and CPL as specified parametric alternatives.
3. Positive-H flexible kinematics with a small, fixed number of q or log-H coefficients.
4. The same flexible family constrained to q(z)≥0 over 0–1, with an explicitly separate extension to the full supported redshift range.
5. Curvature sensitivity; additional propagation or anisotropy branches only with their own measurement equations.

For the population-screening stage use common free redshift-bin distance offsets or a shared flexible distance relation. This prevents an imposed ΛCDM curve from selecting the correction. It also deliberately removes sensitivity to a purely redshift-only luminosity term: that term must remain unidentified until external information is introduced.

## 5. How to ask whether corrections can bring the models into agreement

Suppose model A predicts distance modulus μA and model B predicts μB. Write the measured standardized apparent magnitude as

$$m=M+\mu(z;\theta)+\delta M(z,\mathrm{host},\ldots).$$

With this sign convention, B reproduces A's apparent-magnitude curve if

$$\delta M_{\mathrm{required}}(z)=\mu_A(z)-\mu_B(z)-\mathrm{constant}.$$

The constant is absorbed into M. A correction to an already inferred distance has the opposite sign when the luminosity term is subtracted. Record the convention in every table.

This equality demonstrates an identification limit: unrestricted luminosity evolution can exactly mimic a different distance history. It is not evidence that the required evolution exists.

E05 should produce a **reconciliation budget** for each pair:

1. Calculate the required magnitude-shape difference at the actual redshifts, profiling the common intercept and preserving the same sample.
2. Construct signed responses J to physically named perturbations: zero points, passbands, Galactic reddening, extinction implementation, intrinsic colour, host properties, selection, and a specified age law.
3. Identify the uncertainty of each perturbation using independent calibration, host, or other data. A stress-test amplitude is not an empirical uncertainty.
4. Solve the locally linear diagnostic

$$\min_a \left[(\Delta\mu-Ja)^T W(\Delta\mu-Ja)
 +(a-a_0)^T S^{-1}(a-a_0)\right],$$

where W is the inverse remaining covariance projected off the intercept, after removing only nuisance components explicitly represented by J and S. Retain other statistical and systematic uncertainties. Statistical-only weighting is appropriate only when all material nonstatistical terms are represented explicitly. S describes the independently constrained shared corrections. If S is not known, show an amplitude-versus-mismatch sensitivity curve and label the missing constraint.
5. Report residual mismatch, individual signed corrections, correlated prior displacement, observable predictions, and sensitivity to priors. Do not convert the penalty to a universal “sigma” without calibrating the statistic and dimensions.
6. Follow promising cases through nonlinear flux fitting, population inference, selection, and cosmology. Linear response is a prioritization tool.

Covariance eigenvectors represent allowed uncertainty directions, with arbitrary sign; they do not reveal the actual signed correction. When an explicit nuisance parameter replaces an already marginalized covariance component, remove the duplicate representation. Compare total marginalized uncertainty with a single coherent generative model.

Classify each reconciliation as:

- **Empirically supported:** correction measured independently, improves held-out observables, and bridges the gap with uncertainty propagated.
- **Allowed but unresolved:** the data permit it, but equally predictive alternatives remain.
- **Disfavoured within tested assumptions:** bridging the gap requires corrections inconsistent with independent measurements or loses predictive performance.
- **Nonidentified:** the observations constrain only a combination of distance and evolution.
- **Not testable with current products:** a named likelihood component or measurement is absent.

“There exists a flexible correction that fits” is not an empirical reconciliation. Conversely, failure of one age template is not disproof of every evolving population model.

The existing local fits provide a useful scale for this experiment:

| Treatment in the existing Pantheon+ + BAO CPL sensitivity study | Conditional q0 mean ± SD | Planning implication |
|---|---:|---|
| No added correction | −0.436 ± 0.081 | Reference for that specified released-distance likelihood |
| Public 114-row host-mass revision | −0.402 ± 0.081 | This particular public correction is insufficient by itself to reverse the sign in that fit |
| Short-delay DTD, cosmic mean, fixed 0.030 mag/Gyr coefficient | −0.132 ± 0.076 | Population transport can substantially weaken the conditional acceleration inference |
| C14 DTD, cosmic mean, same fixed coefficient | −0.034 ± 0.075 | The moment and delay distribution can move the result near the boundary |
| C14 DTD, cosmic median, same fixed coefficient | +0.063 ± 0.074 | A specified correction can reverse the mean; uncertainty still spans zero |
| C14 median with one shared 0.030 ± 0.004 coefficient | +0.042 ± 0.097 | Shared slope uncertainty materially changes the sign probability and precision |

These are existing controlled sensitivities from the [local cosmology report](investigations/cosmology.md), using a declared frozen clock and the same observations. They are not independent empirical measurements of the corrections, and changing the DTD while holding the coefficient fixed is not automatically a consistently transformed host-age law. The full corrected author CMB likelihood was not reproduced. The table explains why transport and common uncertainty deserve higher priority than another unconstrained cosmological fit.

## 6. Detailed experiment protocols

### E00 — Freeze inputs, aliases, prior knowledge and partitions

**Question.** Can another analyst reconstruct exactly which evidence entered each comparison, without counting reused objects twice?

**Work.** Freeze hashes for raw acquisitions, executable binaries, source code, calibration/SED assets, effective runtime configuration, fitting masks, selection assets and derived tables. Preserve literal, recovered and amended branches separately. Reconcile current outputs with older narrative reports; a successful result requires source/configuration/output agreement. Snapshot active-work inputs before launching a new dependent run.

Build a physical-supernova alias table across Pantheon+, DES, Foundation, SDSS, ZTF, historical samples and any optical/NIR extension. Use names plus sky position/redshift and manual review of ambiguities. Preserve repeated observations and their covariance; grouping them does not mean deleting them. Record overlap in light-curve training, host inference, calibration standards, nearby anchors and external-probe inputs.

Replace any validation partition based only on survey+CID where the same explosion can cross surveys. Freeze physical-object groups, survey/field blocks, redshift bins, model roster, score targets and prior families. Record which outcomes were already inspected. A deterministic new split of already studied data is locked reanalysis validation, not untouched prospective evidence.

**Output/gate.** One manifest, alias/dependency graph, stable mask ledger and fold file; every row has a source and disposition. Undefined aliases or mutable inputs block claims of independent replication, not unrelated development work.

### E01 — Close the calibrated-flux reproduction gap

**Question.** Which remaining differences arise from initialization, clipping, quality cuts, numerical fitting or covariance conventions?

**Inputs.** Original DES v1.3, current and historical local SNANA builds, published accepted-epoch LCPLOT products, recovered refits and full Hessian outputs.

**Protocol.** Retain the original all-object denominator. For each discrepant or failed object compare available and accepted epochs, bit masks, phase/wavelength support, fluxes/errors, initial parameters, iteration history, fit bounds, objective value and covariance. Run prespecified starts and a fixed published-epoch diagnostic separately from the operational selection algorithm. Public parameters may initialize a reproduction diagnostic but cannot be treated as an independent fit result.

Separate genuine same-likelihood multiple optima from different effective likelihoods. Quantify discrepancies in mB, x1, c, standardized brightness and full-covariance cosmological projection. Regenerate standardization/BBC for any proposed production branch. Never reuse old bias corrections with changed fits as though the pipeline were closed.

**Gate.** Exact reproduction requires matching effective likelihood, membership and documented tolerances; a small median is insufficient. If exact settings remain unavailable, permit a separately named independent pipeline only after its own selection, synthetic coverage and observation-level adequacy pass. Its conclusion cannot be advertised as an exact official reanalysis.

Proposed numerical budget for an accepted independent implementation: residual numerical effects below 0.1 posterior SD in the primary acceleration estimand and below 10% of the smallest correction effect the experiment is powered to discriminate. Define object-level tolerances from serialization precision and actual error scales. These future validation budgets do not retroactively turn the failed historical gate into a pass.

### E02 — Test repeated measurements, flux noise and extraction

**Question.** Is the model assigning the correct amount of independent information to the flux measurements?

**Inputs.** The audited SMP/DIFFIMG tables, repeated-record ledgers, accepted-epoch map, off-signal tables and calibration/exposure metadata.

**Protocol.** Resolve repeated records using exposure, detector, source position and extraction provenance. The observed A,B,A,B pattern is not automatically four independent exposures or four identical values. Compare only justified alternatives: original treatment, removal of demonstrated serialized duplicates, or explicit same-exposure covariance. Ambiguous provenance stays a sensitivity bracket.

For off-signal data prespecify phase windows, bands and quality masks; account for late SN light, reference subtraction, host background and selection on measured flux. Estimate bias, central scatter, tails, time/band correlations and shared reference-image terms. Resample by physical SN and relevant exposure/calibration blocks, not by millions of purportedly independent rows.

The inspected SMP >200 observer-day diagnostic includes 26,109 rows from 832 objects and has a near-zero mean but heavy tails; no SMP >365-day result is present in that diagnostic. The strong DIFFIMG offset is a different reduction and cannot simply be subtracted from SMP. See [saved diagnostics](../phase2/assumptions/off_signal_statistics.json).

Fit Gaussian, robust-tail and empirically constrained covariance alternatives on training noise data. Carry them into SN amplitude fitting and rerun measured selection. Test held-out flux residuals and injected recovery, not merely reduced chi-square.

**Output/gate.** A noise/extraction error budget and its propagated luminosity/acceleration response. A material unresolved effect escalates E16. A flux-tail anomaly alone does not establish a cosmological bias.

### E03 — Rebuild the full selection denominator

**Question.** Would each candidate universe/population produce the observed selected sample under the real observing and analysis process?

**Protocol.** Track generated attempt → sky/cadence/host → trigger → measured photometry → host-z success → fitted light curve → quality/classification → bias-grid eligibility → final membership. Preserve failed attempts and reasons.

Use the generated-attempt index as the key where early rejections reuse CIDs. Verify every simulated output and manifest before crediting a variant as completed. Nominal P21, alternative BS21/G10 populations, luminosity offsets/drift and noise perturbations should pass through the same measured-flux fitting and cut chain. A common RNG seed is not proof of identical underlying hosts/cadences; compare streams before paired estimation.

The released selected mocks supply useful conditional information, but at fixed latent conditions their source luminosity relation can have inadequate support for arbitrary luminosity changes. Importance reweighting is admitted only after a support test, effective-weight analysis by field/redshift/host, and agreement with newly generated alternative populations. Validate relative normalization separately from absolute detection efficiency.

Model uncertainty in trigger and host-z efficiencies using underlying fake/target counts if available; otherwise publish sensitivity families. Do not invent binomial calibration precision from a grid without its denominators. Reproduce the historical expected-SNR trigger, and distinguish it from cuts on measured SNR.

**Gate.** Successful recovery of selection probabilities and primary estimands for known generated populations, including alternatives. Failed support or an unknown final cut blocks decisive evidence ratios for that region. An explicitly new scientifically defined sample may be valid if its own complete selector is modelled and applied equally to data and simulations.

### E04 — Attribute calibration, extinction and covariance changes

**Question.** Which measured changes explain differences between releases, and which remain model uncertainty?

**Protocol.** Construct a signed correction ledger at the stage each effect enters: flux, bandpass, trained SED, fitted parameters, population, selected distance or covariance. Audit units, redshift frame, observer/rest wavelength and whether a correction is already incorporated.

Use the amended SALT response pilot only with matching executed source, inputs and outputs before expanding to its prespecified redshift/coverage sample. At the final planning check, all 20 input/output hashes in the corrected 12-object manifest matched; its summary records the corrected conventions and successful nonlinear checks. This verifies that saved run's provenance, not its full empirical validity. Check direct flux normalization, finite differences, nonlinear recovery, model covariance and the MW map convention. Existing public HEAD MWEBV values and fresh SFD-map values require different bookkeeping; apply the scale once to the appropriate representation.

Separate historical polynomial F99 approximation from the current spline implementation. Propagate differences through broadband SEDs and fits; do not report a monochromatic maximum as a distance bias.

For DES→Dovekie compare fixed common-object branches and native release samples separately. Use a small factorial design for calibration, extinction implementation, SED training, population/selection and covariance changes; report interactions or explicit unresolvable attribution. Update selection after flux changes.

Maintain the covariance/precision distinction: reconstruct a covariance before selecting a subset from a precision release. A subblock of precision is generally a conditional precision, not the marginal precision for the selected data. Preserve signed physical modes and cross-survey calibration correlations.

**Output/gate.** A reproducible stage-by-stage attribution table and uncertainty model. The conflicting signs in the Hoyt description and the public low-z mass-revision vector remain separate hypotheses until the executed vector is available.

### E05 — Measure the correction required by each alternative

**Question.** Can empirically constrained nuisance changes actually reconcile competing expansion histories?

**Protocol.** Apply the reconciliation calculation in section 5 to the same data, intercept treatment and geometry. Begin with standard-versus-age-template histories, ΛCDM versus best-fitting flexible no-acceleration, and native DES versus Dovekie at matched membership. Include correction interactions only where the linear approximation fails or the physics requires them.

Plot required ΔM(z), the measured correction envelope and unexplained remainder. Along each important degeneracy, report a profile likelihood and the calibration/host data responsible for any finite bound. Give conditional results both with and without imported slope/DTD information.

**Primary output.** A table naming how much correction is needed, how much is measured, which observable tests its sign/shape, and what remains unidentified. This is the highest-value early diagnostic for deciding where to spend later compute and seek missing data.

**Gate.** No physical reconciliation claim until the candidate survives E07/E08/E10/E11. If data support only a joint distance–evolution combination, publish that allowed region instead of choosing a point on it.

### E06 — Test age on matched objects and within redshift

**Question.** Does age predict additional standardized luminosity after accounting for measured dust, colour, stretch, metallicity, host mass and selection?

**Inputs.** Matched G11/R19/modern distance objects, exact cumulative cuts, any available host photometry/spectroscopy and age posterior products. Retain shared objects and correlated corrections.

**Protocol.** Reproduce historical estimators as historical estimators, then compare modern models on identical objects. Use a joint latent measurement model for age, metallicity, mass and attenuation; inspect prior dependence when only asymmetric summaries exist. Age estimates sharing photometric calibration or population-synthesis assumptions require common nuisance terms. If reusing host/SFH posterior samples as likelihood information, account for their original prior and selection: reweight by the new prior divided by the old prior only where support is adequate, or refit the host observations. Treating a posterior as an independent likelihood would count its prior again.

Estimate age–brightness association within overlapping redshift/host/colour support, with common flexible distance offsets. Do not extrapolate across nonoverlapping age populations. Check residual dependence on SNR, dust, mass, fit quality and redshift; distinguish measurement-error attenuation from genuine lack of association.

Specify adjustment sets using the causal graph. Conditioning on variables affected by age or selection can remove a mediated effect or create a collider; present both the total predictive association and the explicitly defined conditional/direct association. A small conditional coefficient is not a proof that age plays no physical role.

Test no added age term, a freely estimated term, and the fixed claimed term. Treat externally measured slopes as one shared quantity with common uncertainty. Reanalyse the cut sequence jointly; overlapping redshift cuts are not independent detections.

**Output.** Coefficient/profile/posterior, incremental held-out luminosity score, support/missingness diagnostics, and joint age–dust–metallicity contours. If true progenitor ages are absent, label the result as a host-age association.

### E07 — Distinguish dust from intrinsic colour and luminosity

**Question.** Do competing physical explanations predict distinguishable wavelength, phase and environmental patterns?

**Protocol.** First compare effective SALT-coordinate models on identical observables. Then confront physically interpreted dust models with signed multiband fluxes and actual bandpasses, accounting for rest-frame wavelength coverage at each redshift.

Hold out filters or phase ranges to test predictions from the remaining light curve. Retain the actual discovery and final-selection conditioning when a withheld band/epoch helped select the object: withholding it from fitting does not remove its role in selection. Separate within-object held-band prediction, new-object prediction and new-survey transfer: they answer different questions. Low-reddening/high-SNR strata are diagnostic, but must keep their selection function and may have little power to estimate a dust law.

Compare universal dust, host-dependent dust, dust plus an intrinsic host term, age plus dust, and a small flexible empirical alternative. Use spectra or independently constrained reddening indicators where available. Galaxy attenuation is not automatically SN sightline extinction; SALT c is not a direct E(B−V) measurement.

Add an optical–NIR test after acquiring its actual data and uncertainty model. [RAISIN](https://arxiv.org/abs/2201.07801) supplies a relevant existing-data route, with higher-redshift HST NIR observations and a low-z CSP comparison. These include PS1/DES discoveries: the advantage is different wavelength sensitivity, not independent explosions. Selection and NIR host dependence still matter.

**Output/gate.** A wavelength/phase-resolved predictive comparison. If age and dust explain the same measured subspace, report observational equivalence and quantify which extra bands or spectra would separate them.

### E08 — Test transport from hosts to the selected progenitor population

**Question.** Does a relation measured in one host sample imply the redshift-dependent luminosity change required for another?

**Protocol.** Model progenitor delay τ using host star-formation histories convolved with a delay-time distribution, then apply the survey's detection, host-redshift and quality selection. Distinguish p(τ|z) for the cosmic rate from p(τ|host,z,selected) for the analyzed sample. Preserve the prior-removal rule in E06 when importing SFH posteriors. Record whether the host sample is galaxy-selected or already SN-event-selected; do not apply an additional event-rate weight to an already event-weighted sample without a generative justification.

Propagate uncertain star-formation histories, DTD shape/cutoff, metallicity, host–progenitor mapping, age-estimator bias and luminosity law jointly. H0 remains relevant to the clock in Gyr even when the SN absolute-distance intercept is free. Keep that clock consistent with the chosen cosmology or state a fixed-clock approximation and test it.

For a luminosity law linear in τ, use the appropriate selected weighted mean in the magnitude likelihood, not an automatically substituted median. For a nonlinear law, average the predicted observable under the full population distribution; do not evaluate the function only at an average age.

Compare alternative DTDs with coefficients transformed consistently when representing the same empirical host relation. Also retain fixed-physical-slope changes as a separately labelled sensitivity. Affine rescaling can preserve a slope×age-change product; nonlinear mappings, measurement error and changing populations need not.

**Validation.** Predict age/colour/stretch/host distributions and within-redshift brightness in a held-out population. Corrected cosmology is a later consequence, not the transport model's training target.

**Gate.** An exact W22 or TITAN-based inference requires the missing age-bearing host inputs/posteriors and executed configuration. Approximate transport can proceed with explicit assumptions, but cannot be relabelled as that reproduction.

### E09 — Establish recovery, uncertainty calibration and power

**Question.** Can this experiment tell the candidate explanations apart, and does it give honest acceleration uncertainty?

**Protocol.** Use a generator independent of the fitted engine where feasible. Separate four checks: numerical integration/derivatives; fixed-truth repeated coverage; prior-predictive posterior simulation-based calibration; and deliberately misspecified physical/noise/selection generators.

Inject both accelerating and nonaccelerating histories, including q near zero. Cross them with no age drift, drift near the reconciliation threshold, the claimed age template, dust drift, intrinsic-host effects, heavy tails, calibration tilts, correlated noise, contamination, redshift/velocity errors and combinations in the weakly identified directions. Do not validate only the easiest nominal generator.

For a transparent stress grid use anchored luminosity changes of 0, 0.02, 0.05, 0.10 and 0.15 mag at z=1, plus the E05 threshold. These values are experimental amplitudes, not measured priors. Use actual observing/redshift/host distributions and explicitly state where they are conditioned upon.

Run a cheap pilot first. For near-final cheap-engine coverage, target approximately 1,000 independent repetitions per key truth region or enough that the binomial 95% Monte Carlo interval has half-width ≤0.02 for nominal 95% coverage. At 1,000 trials the simple approximate half-width near 95% is 0.014. Use an exact/binomial interval in reporting. Full flux simulations can be a smaller strategically chosen validation set, with their weaker precision disclosed.

Measure bias, coverage, interval width, wrong-model selection, false acceleration/deceleration decisions and power for the scientifically relevant correction threshold. Proposed progression criterion: ≥80% power at that threshold with a declared ≤5% false-preference rate under the tested null generators. If unattainable, the experiment remains useful for bounds but cannot be sold as decisive discrimination.

For Bayesian computations inspect parameter and data-dependent SBC quantities, effective sample sizes, rank diagnostics, modes and prior-boundary behaviour. Existing fixed-truth MLE/Hessian checks are not SBC. See [Talts et al.](https://arxiv.org/abs/1804.06788) and [Modrák et al.](https://arxiv.org/abs/2211.02383).

**Output/gate.** A power map and calibrated decision rule, including scenarios where it fails. No nominal high-significance claim beyond the simulations' tail resolution; zero exceedances in B trials supports a finite upper bound, approximately 3/B at 95%, not zero probability.

### E10 — Fit the common-observation model tournament

**Question.** Which candidate correction models make the best adequate predictions without receiving a preferred cosmology in advance?

**Protocol.** Fit the frozen roster to identical physical objects using a common distance treatment, calibration inputs, selection, contamination and scoring representation. Use summary observables only after validating their likelihood approximation against flux-level evaluations across representative SNR, colour, cadence and clipping cases.

Before launching real-data comparisons, export a numeric configuration for every family: parameter definitions and units, prior distributions and bounds, drift bases, host thresholds, covariance components, latent mixture structure, distance knots and integration tolerances. Start from the existing implementation's explicit priors, then freeze justified amendments after prior-predictive checks. Separate empirical calibration constraints from regularizing population priors; broaden/narrow the latter in declared sensitivities. Check predictions over the full redshift/host support and expose any posterior that presses against an arbitrary bound. New A/AD/F/physical-dust branches in this plan require implementation; their inclusion in the table is not a claim that the present engine already fits them.

The selected likelihood must normalize over objects that could have been observed. With y the measurements, h host data, u latent variables and s observing conditions, a schematic form is

$$p(y,h|z,s,S=1,\theta)=
\frac{\int p(y,h,u|z,s,\theta)\,p(S=1|y,h,u,z,s,\theta)\,du}
{\Pr(S=1|z,s,\theta)}.$$

Actual selection may involve further detector/classification variables; include them rather than pretending selection depends only on final magnitude. Shared calibration and training parameters couple supernovae, so integrate those at the joint level. This is a conditional-redshift/fixed-count likelihood. If event counts are added as information, introduce a separate rate/exposure likelihood, including the expected count; do not import rate information implicitly.

Model Ia/non-Ia mixtures using declared class information. A classifier output derived from the same light curve is not independent evidence; BEAMS-updated Ia probability already contains distance information. Either use a validated conditional approximation or a coherent generative/classification treatment. Repeat a high-purity analysis with its changed selection and information loss.

**Score.** Use outer-fold paired log predictive density on common observations, integrated over training uncertainty. Publish the luminosity prediction conditional on colour/stretch/host as well as the joint score; a win from colour-distribution prediction alone does not settle distance calibration. Models using different internal summaries must be scored on a common observable or a demonstrably equivalent likelihood measure.

**Adequacy.** Inspect flux residuals by phase/band, brightness, host surface brightness, redshift, field, SNR, colour, stretch and classification. All models can fail. The least bad model is not therefore physically correct.

**Output.** A ranked set with predictive uncertainty, failures, computational diagnostics and a record of which observables support each preference. Only adequate models advance to strong cosmological interpretation.

### E11 — Test transfer across the dimensions that could expose a hidden correction

**Question.** Does an explanation work outside the data conditions that selected it?

**Protocol.** Use four separate transfer tests: new physical objects under similar conditions; held-out fields/depths/calibration groups; redshift-range extrapolation; and survey/wavelength transfer. Prespecify overlap/support requirements for each.

Train population relations on one appropriate low-z sample and predict within-redshift brightness/colour relations in another; reverse the direction where data permit. DES and ZTF offer different observing programmes, but training/calibration dependencies remain. Pantheon+, DES revisions and legacy compilations overlap heavily and cannot be multiplied as independent confirmations.

The [September ZTF BayeSN revision](https://arxiv.org/abs/2605.06799v3) is an important external challenge to a universal dust-only explanation. Reproduce its actual masks, preprocessing and posterior products if available. Its environmental proxies are not a direct progenitor-age measurement, and a residual environmental step alone does not measure redshift evolution.

For common SNe observed in optical and NIR, evaluate paired residual differences with their covariance. For disjoint surveys, integrate shared calibration and permit independently constrained survey-specific noise/population features. Do not fit an unrestricted survey offset that silently erases the correction being tested.

**Gate.** A model preferred only within one field or one wavelength range is a local explanation. A robust correction requires consistent transfer or an explicitly justified domain of validity.

### E12 — Infer acceleration under the admitted correction models

**Question.** Which claims about acceleration survive the models that passed measurement and predictive tests?

**Protocol.** Run SN-only first with a free absolute-magnitude/H0 intercept. Fit the finite-resolution kinematic and parametric families using the same admitted nuisance models. Infer positive H and integrate to distance, rather than differentiating a noisy unconstrained distance curve.

Report the primary recent-interval quantity, q0, and whether an accelerating interval is required. Repeat fixed, prespecified smoothness/binning and curvature sensitivities; give effective resolution and prior-to-posterior information. Include peculiar velocities/coherent flows at low z, redshift uncertainties, lensing and relevant common calibration.

Compare the flexible unconstrained family with the same family restricted to q≥0. Calibrate the likelihood-ratio statistic under constrained null simulations, nuisance uncertainty and boundary cases. A bootstrap only at one convenient fitted null is not proof of uniform control over a composite null; test difficult/least-favourable regions identified by simulation.

Report posterior \(P(\bar q<0)\) and \(P(q_0<0)\) within each stated prior/model, plus calibrated frequentist results where available. A probability from one conditional fit is not a model-free probability that the Universe accelerates.

An affirmative claim of positive average q needs its own prespecified one-sided test of \(\bar q\leq0\), with false-deceleration control under accelerating and boundary generators. Provisionally require a calibrated 95% lower confidence bound above zero for this interval average, with adequate predictive fit; report the corresponding Bayesian sign probability separately. This would not establish positive q at every redshift in the interval. Failure to reject the no-acceleration family is not that test. A claim of approximate coasting instead requires an equivalence margin fixed from the scientifically relevant distance-shape scale and E09 power before examining its real-data result.

**Output.** A robustness matrix of acceleration estimands versus correction family, geometry, selection and prior. If plausible predictive peers yield different signs, report unresolved acceleration under those assumptions rather than average away the disagreement.

### E13 — Predict other probes before combining them

**Question.** Can the same inferred expansion history explain non-SN measurements without having used those measurements to tune the correction?

**Protocol.** Freeze the SN correction candidates, then compare their inferred expansion shapes with DESI BAO geometry. If the free ruler scale is fitted on those BAO data, label this a conditional/profiled cross-probe shape test, not an untouched SN-only prediction of the BAO distances. For a literal predictive test, constrain the ruler only with declared training/external data and score held-out BAO. Conversely fit BAO alone and calculate the SN luminosity evolution it would require, before imposing that evolution. These directional checks expose tension that a joint best fit can conceal.

Next run joint SN+BAO, retaining any relevant shared nuisance or selection dependence. Add CMB in a separately labelled branch using a documented full likelihood or clearly named compression with its validity tested for the proposed model. Keep sound-horizon, BBN, early-dark-energy, gravity, perturbation and foreground assumptions explicit.

Do not use an Ωm Gaussian as a substitute for reproducing the author's complete corrected CMB result. Do not count reference chains as a new rerun. An exact Son-style CMB claim remains an access/reproduction target until its configuration, correction covariance and posterior products can be reconstructed.

**Output.** SN-only, BAO-only, predictive cross-check, SN+BAO and full/validated-compressed CMB results in separate columns. A combined preference can be driven by one probe; identify that contribution.

### E14 — Reconstruct the historical result directly

**Question.** What survives when the actual discovery analyses are repeated and then updated transparently?

**Protocol.** Acquire the original machine-readable tables, photometry and calibration where available; record unavailable assets. Reproduce the original sample, light-curve method, nuisance/intercept treatment, peculiar-velocity prescription, colour/dust cuts, cosmology and statistical assumptions separately for the two teams. Compare published likelihood contours and tables.

Then make one justified change at a time: modern calibration, revised redshifts/flows, extinction, population/host information, covariance and selection. A correction learned from modern data must predict its effect on the legacy sample under an explicit transport model. Do not import a mean modern shift blindly.

Keep analysis of original objects separate from a modern much larger sample. The historical calculation may be reproducible while its quoted uncertainty needs revision; modern independent evidence may remain even if a historical implementation is corrected.

Add the high-z distance-turnover test represented by [Riess et al. 2004](https://arxiv.org/abs/astro-ph/0402512). Its value is redshift-shape discrimination against specified simple dust/evolution models; it cannot exclude every arbitrary evolving luminosity function. Account for object overlap with later compilations.

**Output/gate.** “Original reproduced,” “original modified by measured correction,” or “exact legacy reconstruction unavailable,” alongside the separate modern robustness conclusion. Published-table reproduction and raw-photometry remeasurement must not share the same label.

### E15 — Audit broader geometry and propagation alternatives

**Question.** Could a specific velocity, anisotropy, curvature or photon-propagation effect account for the remaining signal?

**Protocol.** Extend the existing directional work only with correct coordinates/redshift frames, survey angular selection, calibration fields, flexible population drift, exact distance relations where required and full covariance. Check low-z velocity fields and sky/depth coupling before interpreting a preferred direction.

Calibrate scans over directions and decay scales using null simulations with the actual geometry. Parameters absent under the null invalidate an automatic Wilks conversion. Revised or withdrawn results cannot be treated as current independent evidence.

For opacity/grey-dust or distance-duality alternatives write the additional observable predictions and external constraints. For modified gravity or inhomogeneous expansion specify luminosity distance, rate/population clock and perturbation treatment; do not append a label to a fitted q curve. Curvature is a prespecified sensitivity, while wholly new physical frameworks require their own controlled protocol.

**Output.** A bounded test of each specified alternative and its distinguishing predictions. This branch should not delay the measurement/selection closure unless existing evidence points to a material effect.

### E16 — Escalate to images and SED training where justified

**Question.** Are the remaining uncertainties caused by common assumptions below released fluxes or in trained spectral surfaces?

**Protocol.** Start with representative discrepant/control objects. Inventory search/reference images, PSFs, masks, zero points and processing versions; public archive availability does not prove recovery of the exact original inputs. Independently remeasure flux or inject artificial sources only where those ingredients are adequate.

Test host-background and template-noise effects at the pixel level and compare their predictions to E02. Account for the fact that reference-epoch SN flux constrained to zero cannot independently verify its own zero-flux assumption.

For shared SED training, identify training/test overlap, calibration links, phase/wavelength support and regularization. Compare independent training/model families where available, or retrain targeted variants with complete calibration propagation. A second fitter using the identical trained SALT surface is computational independence, not independence of the standardization model.

**Gate.** Produce a measured explanation or explicit remaining floor. Do not require an entire survey pixel reconstruction before exploiting valid lower-cost tests, but do not label calibrated-flux results detector-level validation.

### E17 — Design the smallest decisive additional measurement

**Question.** Which missing measurement most reduces uncertainty in the contested correction and acceleration sign?

**Protocol.** Use E05's weak directions and E09's realistic simulations to compare: host spectroscopy separating age/metallicity; UV–MIR host SEDs; optical–NIR SN light curves; repeated calibrated measurements; additional high-z objects with characterized selection; or improved host-z/fake-injection denominators.

Prefer existing archival measurements first. Rank proposed data by expected reduction in the correction–distance degeneracy, model-discrimination power, calibration dependence and cost. Adding many optical SNe under the same unresolved systematic can be less useful than a smaller sample with a new observable.

Estimate required numbers using the full correlated model and selection. Simple \(N^{-1/2}\) projections are only sanity checks. Specify redshift/host coverage and observing precision, not just a sample count. Include publication/availability of calibrated inputs, failed detections and selection denominators in any proposed new programme.

**Output.** A quantitative design with an explicit go/no-go criterion. New telescope time is justified only after existing-data limits and the expected information gain are demonstrated.

## 7. Statistical rules for deciding “more likely to be correct”

**Separate prediction, parameter inference and physical truth.** A higher predictive score supports a model's ability to describe the tested observations. It is not automatically the posterior probability of a physical mechanism. A p-value is not the probability a model is false.

**Use an outer evaluation loop.** Freeze the main roster before outer-fold scores are opened; tune nuisance complexity only inside training/inner folds. If the same observations already informed the roster, disclose that history and reserve stronger confirmation for genuinely new data or a separately frozen cross-survey test.

**Respect dependence in scoring.** Aggregate by physical SN and calibration/field blocks. For global shared parameters calculate appropriate joint-block predictive probabilities; do not assert that summing independent pointwise scores captures all dependence. Use paired score differences on identical held-out data, with block uncertainty or a generative replication distribution. A small number of fields limits the precision of field-bootstrap uncertainty.

**Require relevance, not only a large total score.** Publish joint flux/summary scores and conditional luminosity scores. A model may fit the colour distribution better without reducing uncertainty in distance evolution. Give per-redshift and per-environment score contributions and influence diagnostics.

**Treat predictive thresholds as an operational rule to calibrate.** A paired score advantage greater than roughly two estimated standard errors is a useful screening criterion, not universal evidence of a true mechanism. Before a strong claim, it must survive E09's false-preference calibration, observation-level adequacy and E11 transfer. Multiple model comparisons and exploratory cuts belong in the simulated decision procedure.

**Use Bayesian evidence only when it is meaningful.** If reporting Bayes factors, use proper, scientifically stated priors, identical observations, normalized selection likelihoods and stable evidence estimates. Report sensitivity to prior widths and measures. Convert to model probabilities only conditional on an explicit model set and prior model odds, while acknowledging model incompleteness. Do not call AIC/BIC, normalized likelihoods or predictive stacking weights probabilities of physical truth.

**Use likelihood tests with the right null distribution.** Non-nested families, mixtures, boundaries, unidentified directions and optimized sky scans often require simulation rather than a standard chi-square reference. Report both effect sizes and uncertainty.

**Retain absolute adequacy.** If every family fails residual, selection or calibration checks, the outcome is model inadequacy, not a winner. Broaden only the failure implicated by the data and record the amendment.

**Check numerical error separately.** Require convergence across chains/starts, suitable ESS for reported tail probabilities, likelihood/quadrature agreement, stable gradients and positive well-understood covariance. A nominal R-hat threshold such as 1.01 is necessary monitoring, not proof of mode exploration. Monte Carlo uncertainty in a score should be small relative to its sampling uncertainty, provisionally below 10%.

**Do not hide model uncertainty.** Give within-model and across-model acceleration results separately. If reporting a model-averaged quantity, identify the weighting rule and also display disagreement among well-predicting models. Predictive weights need not be posterior model probabilities.

The rationale for out-of-sample predictive comparisons is developed in [Vehtari, Gelman and Gabry](https://arxiv.org/abs/1507.04544); the block structure, luminosity-specific target and decision gates here are proposed for this project.

## 8. Missing products and the experiments they unlock

Use the existing [author-data request ledger](author-data-requests.md), updated against current public releases before contacting anyone. No external messages are part of preparing this plan.

| Missing or incomplete product | Why it matters | Useful work without it |
|---|---|---|
| Exact effective DES fit/rejection runtime, full final cut/BBC chain | Exact official flux-to-cosmology reproduction | Independent, separately validated pipeline and detailed mismatch attribution |
| Host-z efficiency target/success counts and fake-injection uncertainty products | Empirical uncertainty of the selection model | Transparent sensitivity envelopes; no invented grid precision |
| W22/W26 executed age-bearing host population and crosswalk | Exact selected-progenitor transport and simulation | Explicit alternative SFH×DTD transport, labelled approximate |
| Joint host age/mass/metallicity/attenuation posterior or likelihood products | Separating confounded host effects and errors | Bounded sensitivity using summaries; no replacement by zero-error proxies |
| Son correction configuration, executed vectors/covariance and corrected CMB chains | Exact claimed corrected cosmology and significance | Existing fixed-template sensitivity and independently specified models |
| Hoyt's executed low-z correction vector | Resolve sign discrepancy against the positive public revision | Treat published text and public vector separately |
| Exact ZTF analysis masks, preprocessing and posterior draws | Reproduce environmental result and perform shared-observable comparison | Public photometry reanalysis under a separately declared pipeline |
| Complete legacy photometry/calibration and raw-image provenance | Exact historical/pixel reconstruction | Published-table reproduction with explicit boundary |
| Optical/NIR data, covariance and selection assets for extensions | A wavelength discriminator | Acquisition/overlap feasibility audit; no assumption that an abstract is a usable likelihood |

A missing product should block only the claim it is necessary for. It should not stop unrelated tests that can be completed with current public data.

## 9. Execution programme and practical stopping rules

These are suggested effort allocations, not measured runtime promises. Benchmark one full chain before budgeting the large simulation campaign.

**First work package: establish a reliable baseline.** Run E00 and audit E01–E04 in parallel. Deliver the refreshed status, exact unresolved refit causes, corrected dust pilot, exposure/noise evidence and a verified generated-to-selected chain. Allocate approximately one third of initial effort here. Stop broad new cosmological sampling while consequential implementation differences remain unexplained.

**Second work package: quantify discriminating corrections.** Run E05–E08 while E09 expands from cheap pilots to survey-like recovery. Deliver the required-versus-measured correction budget, within-redshift associations, wavelength diagnostics and transport uncertainties. Allocate approximately another third here. If the primary age/dust directions are not identifiable, state that early and redirect effort toward E11/E17.

**Third work package: compare and interpret.** After gates pass, run E10–E13 on the frozen roster. Deliver common-observation predictive comparisons, external transfer and the acceleration robustness matrix. Use the remaining core effort here, preserving an independent audit budget. E14 proceeds as a bounded parallel historical thread; E15/E16 expand only when evidence or the requested historical scope justifies their cost.

With four concurrent investigators, a useful allocation is: (1) measurement/calibration/refits; (2) selection/forward simulation/noise; (3) host-age/dust/transport evidence; (4) likelihood/statistics/integration. Reassign one investigator to an independent challenge of the final ranking before declaring completion. Reusing the same data remains explicit regardless of how many implementations agree.

Provisional computational stages: analytic/unit checks; 12–64 representative flux objects; all selected objects; cheap repeated simulated catalogues; targeted full flux simulations in the dangerous parameter directions. Increase fidelity and sample count only when it can change the conclusion or validate a needed approximation.

**Stop or narrow a branch when:**

- Its uncertainty interval or calibrated equivalence test excludes the scientifically relevant correction scale. A small or nonsignificant point estimate, even in a well-powered design, is insufficient.
- It is mathematically nonidentified by the available observables and no new constraint is available.
- Importance weights have inadequate support; regenerate rather than increase sampler iterations.
- Every model fails the same diagnostic; repair that likelihood layer before ranking.
- A prior, hard bound or smoothness choice determines the sign; report the sensitivity.
- Legacy/private products prevent exact reproduction; finish the bounded public-data result and label it.

## 10. Deliverables and possible conclusions

The completed programme should provide:

1. A versioned evidence and overlap map, with source hashes and explicit data levels.
2. Reproduction records separating exact, approximate and independently validated pipelines.
3. A signed correction ledger with units, physical source, already-applied status, common uncertainty, interactions and nonlinear checks.
4. A required-versus-measured correction budget for every central alternative.
5. A selected-sample likelihood with tested normalization, contamination and calibration treatment.
6. Recovery, power and failure maps for correction and acceleration inference.
7. Paired held-out prediction tables on common observations, including absolute model inadequacy.
8. Separate SN-only and external-probe acceleration results, with prior/geometry/selection/model robustness.
9. The historical reconstruction, clearly separated from modern validation.
10. A quantified list of remaining identification limits and the smallest observations that would address them.

Use one of the following conclusions, or a combination when they answer different questions:

| Finding | Defensible conclusion |
|---|---|
| Measured corrections improve observations but recent acceleration survives across adequate models | The acceleration inference is robust to the tested measured corrections; numerical values or uncertainties may still change. |
| A correction is reproducibly measured, transfers, and substantially changes cosmological parameters while acceleration persists | Revise the standardization/cosmological estimate, without claiming the discovery is disproved. |
| Adequate models spanning both signs remain predictive peers | Existing data do not identify the present/recent acceleration sign under the tested correction freedom. This is loss of support, not positive evidence for deceleration. |
| A calibrated interval-average sign test supports positive average q after independent correction validation, and the model passes predictive checks | Evidence for positive average q over the defined interval and assumptions, not deceleration at every redshift. A weaker preference for a no-acceleration family alone does not establish this. State separately what is known about q0 and earlier epochs. |
| ΛCDM loses to an evolving model while q remains negative | Evidence against the cosmological-constant model in that comparison, with acceleration retained. |
| The required correction conflicts with independent measurements and fails transfer | The specified alternative correction is disfavoured; this does not establish that every possible correction is absent. |
| All candidates fail, or required inputs are missing | The tested model set or reproduction boundary is inadequate; publish the failed discriminator and required next measurement. |

The most useful first result is the E05 reconciliation budget grounded in E01–E04, followed by E06/E07's observable discriminators. Another plot showing that a sufficiently large imposed age correction changes q would add little. The real advance is showing whether that correction is independently measured, survives selection and transport, and predicts observations that its competitors do not.
