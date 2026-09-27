# What native supernova corrections learn from an imposed age signal

In this conditional simulation, a correction trained to treat an imposed age term as an unknown luminosity mismatch substantially reduces its redshift-dependent mean distance shift, while a strong age association remains within redshift bins. Those are different quantities. A residual age slope alone does not show that a cosmology-relevant population-mean bias remains uncorrected. These calculations establish that distinction for a declared mock population; they do not establish the true supernova age law or a correction to observed cosmological distances.

## Physical intervention and training target

We generated 240,000 attempts per arm, in sixteen independent paired shards, through DES cadence, noise, detection, native light-curve fitting and recovered classifier selection. One arm receives

\[
\Delta m_{\rm age}=-0.030\,(t_{\rm mock}-3)\ {\rm mag}
\]

before photon generation. Counterparts share cadence, host and dust draws, latent parameters and random seeds. The training seeds are distinct from the existing 3,000-attempt paired evaluation sample. The final quality/classifier-selected training sets contain 14,043 nominal and 13,724 age-injected objects. The independent photon audit confirms the sign and amplitude of the injected change, without changing true distance or the separate native intrinsic-amplitude field.

Ages are W22 **simulated progenitor delays**, not observed or independently measured ages. Both arms use positive-support F99 dust screens with fixed \(R_V=3.1\), and noiseless mock host masses. The SALT3.DOVEKIE surface and recovered SNN classifier remain fixed. These pure-Ia calculations measure signal acceptance; they omit core-collapse contamination, classifier retraining, SALT-surface retraining and host SED measurement uncertainty.

The native `opt_biascor=4336` calculation retains its field groups, redshift/width/colour/mass grid, SNR>60 intrinsic-scatter sample, minimum counts, interpolation and MUCOV rules. We compare two age-generated correction-training targets:

- **Literal generated truth:** the direct-MU bias target subtracts the known `SIM_gammaDM` grey shift. It therefore does not ask the mean correction to learn the entire imposed age signal as an unknown bias.
- **Retained age term:** set only training `SIM_gammaDM=0`, while preserving the original shift in `TRUE_AGE_MAGSHIFT`. All photons, measured fits, generated distances and evaluation rows remain identical. This lets the mean correction learn the imposed luminosity mismatch. It also changes the native gamma zero point and may change MUCOV; it is not exclusively a mean-map switch.

Native dust mode builds a single alpha/beta correction cell at \(\alpha=0.15,\beta=2.87\), even when final nuisance coefficients float. Thus the following refits do not continuously regenerate the correction surface at each fitted alpha and beta.

## Support and nuisance identification

Training support grows materially with computation. At 15,000 attempts per arm no evaluation object passes the native multidimensional correction. At 60,000, the four native cases retain 19/16/15/15 objects; at 120,000 they retain 51/50/47/47. Only 26 objects in the latter four-case intersection lie in bins meeting the prespecified ten-object minimum, below the thirty-object slope rule. At 180,000, the common-input, fixed-nuisance control has 63 common objects, of which 47 meet the bin rule. All these intermediate records remain available.

At 240,000 attempts per arm, the common-input control retains 94/92/90/90 objects, with **84 in the four-case intersection**, all in supported bins. No scientific support threshold was relaxed. Nevertheless, the full free mass-step fit reaches its lower bound, \(\gamma\simeq-0.5\), and lacks two-sided retained host-mass support. Its age inference is withheld. A native success exit and positive fitted-bin covariance do not make that mass step identified.

A separately declared, exploratory stratum uses pre-BBC host mass \(\log_{10}(M_*/M_\odot)\ge10.04\), reducing the 161 common evaluation candidates to 141. The native host-step factor is then exactly −0.5 at double precision throughout the stratum. Fixing gamma to zero removes a global-offset direction while allowing alpha, beta and the supported native distance bins to refit. All four final designs pass their rank, covariance and alpha/beta boundary checks. The cases retain 92/90/88/88 objects, with **82 common objects**. This measures width/colour adjustment in selected high-mass hosts; it does not measure a mass step or survey-wide standardization.

## Within-redshift age response

The primary quantity is the injected-minus-nominal change in `MU−SIM_DLMAG` on exactly common occurrence IDs. Each comparison removes one weighted global offset. Age slopes control for the prespecified redshift-bin intercepts.

| Correction and nuisance treatment | Fixed alpha/beta/gamma, 84 objects | High-mass stratum, 82 objects |
|---|---:|---:|
| Nominal training and frozen nominal nuisance values | −0.02983 | −0.02980 |
| Nominal training, allowed nuisance refit | −0.02983 | −0.02737 |
| Age training, literal generated-truth target | −0.03120 | −0.02895 |
| Age training, retained-age target | −0.02695 | −0.02791 |

Entries are mag/Gyr. In the second column all nuisance coefficients remain fixed. In the third column only width and colour are refitted; gamma remains the declared zero-point gauge. The imposed slope is −0.030 mag/Gyr.

Each supported variant has 200 joint training/evaluation attempt bootstraps and 100 training-only bootstraps, with coupled draws for nominal/injected counterparts and independent draws between training and evaluation. The all-case age slopes are supported in 139/200 joint replicates and 65/100 training-only replicates. In the high-mass stratum, joint bootstrap standard deviations are 0.00105, 0.00371, 0.00421 and 0.00476 mag/Gyr for the four rows above. The paired width/colour refit change is +0.00243 mag/Gyr; its supported-bootstrap percentile range is [−0.00010,+0.01198]. The retained-target versus literal-target slope change is +0.00104 mag/Gyr, with range [−0.00239,+0.00501]. The present simulation does not precisely identify those small within-bin changes.

These ranges summarize supported resamples, not coverage-validated confidence intervals. Sparse interpolation gates change the retained sample between replicates, and the bootstrap distribution need not be centered on the original point estimate. Slope ratios must not be interpreted as the surviving fraction of a cosmological bias.

## Redshift-mean response

The high-mass common sample contains 12/44/15/11 objects in the four populated bins. The following means are the paired distance changes after a single global normalization, in magnitudes:

| Redshift | Frozen nominal correction | Nominal map, width/colour refit | Literal age target | Retained age target |
|---|---:|---:|---:|---:|
| 0.05–0.30 | −0.02535 | −0.01560 | −0.00389 | +0.00647 |
| 0.30–0.50 | +0.00515 | −0.00109 | −0.00736 | +0.00002 |
| 0.50–0.70 | −0.01459 | −0.00032 | +0.00401 | −0.00887 |
| 0.70–0.90 | +0.05741 | +0.04883 | +0.05408 | −0.00027 |

There is no supported bin above redshift 0.9. Retained-target joint marginal bootstrap standard deviations are 0.0276/0.0159/0.0272/0.0324 mag, with respectively 97/139/96/60 supported replicates. Their full four-bin joint covariance has only **25 complete replicates** and is recorded rather than approximated by independent errors.

An exploratory summary, chosen after seeing these point estimates, subtracts the lowest-bin change from the highest-bin change. This cancels the arbitrary global offset. It is +0.08276 mag under frozen nominal standardization and correction and −0.00674 mag under retained-target training. The paired change is **−0.08949 mag**, with supported-bootstrap standard deviation 0.03596 mag and percentile range [−0.16775,−0.04574]. Only **39/200** joint replicates support both extreme bins. The retained-target high-minus-low range itself is [−0.08823,+0.04953] mag. These sparse, conditional, post-pattern summaries demonstrate the mechanism, not a calibrated cosmological significance.

The central implication is therefore narrower than either “age is already completely corrected” or “the residual age slope proves an extra cosmological correction.” A correction target that retains the imposed age mismatch can remove much of its population-mean redshift distortion while leaving the within-bin age slope. Establishing which target and population correspond to observations still requires empirical progenitor/environment constraints and a survey likelihood with uncertainty in dust, populations and selection.

## Numerical and uncertainty limits

BBC derives field-group redshift map bounds from the evaluation input. We retain the original selected-input pipeline comparison, then use identical pre-BBC occurrence IDs and identical observed redshift/field lists for the common-input alternatives. Actual native grid bounds, cell counts and training counts are checked equal between nominal and injected evaluations using nominal training. Freezing nuisance coefficients externally on the unchanged injected nominal-map output isolates nuisance refitting without silently rebuilding the map.

`MURES` removes fitted redshift-bin offsets and is unsuitable for this drift test. The checked identity is `MU−MUMODEL=MURES+M0DIF`; native MU is independently reconstructed from x0, coefficients, host term, mean correction and one global offset. The observed-redshift truth convention gives identical paired results because counterpart redshifts agree. Native `.COV` is the raw fitted-bin Hesse block with a free global normalization and placeholders for unfitted bins. It differs from the bootstrap covariance of centered means and from survey-systematic uncertainty.

A twenty-replicate intermediate pilot had one native scan-storage overflow and five zero-MAD pull failures. A separately built executable extends only storage through the existing scatter scan; its step, lower bound, estimator, interpolation and scientific gates are unchanged. Four already successful science tables and covariance files are byte-identical in replay. The repaired five-object cell requires 121 steps from a starting value 0.90995 and reaches the unchanged −0.3 lower limit. That exposes sparse-cell variance instability. The five zero-MAD failures persist and are not repaired by changing the estimator.

In the final joint bootstraps, the fixed case has 38 nominal native failures; the high-mass case has 37, plus one replicate without an analyzable slope. These are preserved zero-MAD failures. All-case support is lower still because other arms can fail. Each training-only series has 20 nominal failures. Expanded-storage diagnostics and all unsupported replicates remain in the results. No missing distance, covariance or bootstrap draw is replaced by zero. Independent source and photon checks are summarized in the [method review](survey-bbc-method-review.md).

## Reproduction and records

First restore the public assets, configured native build and pinned classifier environment described in the [survey-physics instructions](../code/survey_physics/README.md), including the positive-RV training/evaluation inputs. Downloaded sources, models and generated arrays remain under ignored `.work/survey-physics/`; they are not vendored.

From the repository root, the enlarged campaign is:

```bash
.venv/bin/python studies/host_ages/code/survey_physics/scaled_native.py prepare --shards 16 --attempts 15000
.venv/bin/python studies/host_ages/code/survey_physics/scaled_native.py run --workers 8 --classifier-python "${SURVEY_CLASSIFIER_PYTHON:?Set classifier Python path}" --classifier-pythonpath "${SURVEY_CLASSIFIER_SOURCE:?Set pinned classifier source path}"
.venv/bin/python studies/host_ages/code/survey_physics/scaled_build.py
.venv/bin/python studies/host_ages/code/survey_physics/scaled_bbc.py merge --shards 16
```

`scaled_build.py` copies the separately configured legacy-sensitivity build, applies `scaled_scatter_capacity.patch` and builds only BBC in a new directory; it refuses to overwrite an existing build identity. It does not change historical simulation binaries. Existing dependencies can be selected with `SURVEY_PHYSICS_LIBRARY_PATH`; no system installation is performed.

Run `scaled_bbc.py bbc --shards 16` and `scaled_response.py --shards 16` for each of four variants: no extra flag, `--common`, `--common --fixed`, and `--highmass-gauge`. Use the same variant flags with `scaled_bootstrap.py --shards 16 --mode joint --replicates 200 --workers 8`, then `--mode training_only --replicates 100`. Unsupported baseline variants write explicit nonexecution records. `scaled_redshift_summary.py --shards 16` derives the exploratory joint-bin contrasts; `scaled_validate.py` verifies inputs, identities, failures and manifests. On a fresh regeneration without the ignored historical pilot, use its explicit `--skip-historical-replay` option: the final scientific inputs and outputs are still checked, while the older regression is recorded as archival evidence that was not rerun. `scaled_capacity_replay.py` checks the preserved local intermediate pilot when those generated files are available. A fresh checkout can reconstruct the bounded regression after the first twelve native shards and fixed evaluation pool exist:

```bash
.venv/bin/python studies/host_ages/code/survey_physics/scaled_capacity_replay.py --restore-fixtures
.venv/bin/python studies/host_ages/code/survey_physics/scaled_validate.py --restored-fixtures
```

This uses the **original** BBC executable for four baseline cases and twenty nominal pilot cases at seeds 972000–972019, then the **separate repaired** executable for the four successful baselines and six failures. The new sidecar writes ignored `capacity-fixtures/` and a separate `scaled-native-capacity-restoration.json`; it leaves the original pilot records intact. It restores the regression fixtures, not the entire historical bootstrap experiment, and generates no new photons. Regenerated timestamps and provenance hashes can differ; agreement is assessed from the actual scientific rows, covariance and failure semantics. Existing restored fixtures are not overwritten.

Compact records are in [survey-physics results](../results/survey_physics/). `scaled-validation.json` binds source code, final inputs, native executables and generated data; `scaled-redshift-summary-k016.json` includes full cross-bin/cross-estimand covariance and support counts. Original lower-volume campaigns remain unchanged.
