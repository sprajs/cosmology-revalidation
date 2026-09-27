# Independent environmental measurements and the extra age correction

The new observed-host comparisons do **not** provide positive evidence that the proposed full age correction should be appended to standardized supernova distances. They also do **not** prove that dust and standardization have removed every physical age effect. The distinction matters: the best available new age catalogue contains model-dependent posterior summaries with substantial uncertainty, while the small nearby sample with independently measured host environments has insufficient statistical power.

These are new comparisons, rather than an exact reproduction of either side's cosmological result. The [code and execution instructions](../code/environment_validation/README.md), [initial design](../code/environment_validation/design.json), [TITAN extension design](../code/environment_validation/titan-design.json), and compact [numerical records](../results/environment_validation/) retain their scope and source identities.

## 1. Local and global host properties are measurably different

The [Kelsey DustPedia supplement](https://doi.org/10.1093/mnras/stag1765) provides 90 local 3-kpc apertures in 78 hosts, with 41 local 1-kpc counterparts. All joins are unique by supernova or host as appropriate. We independently recover the published global-minus-local medians. Resampling whole hosts 2,000 times gives:

| Property | Global minus local 3-kpc median | Host-cluster bootstrap 95% interval |
|---|---:|---:|
| Mass-weighted stellar age | −0.312 Gyr | [−0.480, −0.042] |
| Host-population attenuation, A_V | −0.0134 mag | [−0.0326, +0.0013] |
| Rest-frame u−r colour | −0.0630 mag | [−0.0807, −0.0296] |
| log specific star-formation rate | +0.104 dex | [+0.006, +0.196] |

The object-to-object normalized median absolute deviation of the age difference is **1.087 Gyr**. Thus a global age cannot simply be substituted for the explosion-site age. The median 1-minus-3-kpc age difference is only −0.004 Gyr, but its object-to-object NMAD is 0.359 Gyr; absence of a mean offset does not make individual apertures equivalent. Their overlapping photometry also means that adding local/global marginal errors in quadrature is not a validated joint uncertainty.

Among the 90 local measurements, age correlates with u−r (Pearson 0.832), specific star formation (−0.844), and mass (0.548). Its simple linear correlation with attenuation is only −0.121. These are correlations **between SED posterior summaries**, not the unavailable age/dust covariance inside each host's posterior. In particular, neither strong colour correlation nor weak attenuation correlation identifies progenitor age or excludes dust. Kelsey's A_V describes an extended stellar population, not point-source extinction along the supernova line of sight.

## 2. Nearby distance anchors give a useful but underpowered age check

We joined DustPedia to Pantheon+ using normalized supernova names and an independent 3-arcsec sky-position check. Coordinates identify two additional legitimate release aliases: `2008fv_comb` and `1994DRichmond`. Their observing epochs match the named event years. SN 1998aq has a 4.59-arcsec catalogue discrepancy and remains excluded under the declared threshold; it was not silently repaired. Nineteen distinct SNe pass, of which **11 SNe in nine Cepheid hosts have 19 released light-curve rows**.

For these calibrators we compare `m_b_corr − CEPH_DIST` with the local CIGALE age, fitting a free brightness intercept. The pinned [Pantheon+ distance README](https://github.com/PantheonPlusSH0ES/DataRelease/blob/7fc680548d1ea9fe6ee983a89d8b5635bb5784a8/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/README) explicitly states that the released covariance includes Cepheid-host uncertainty and correlations between repeated measurements of the same SN. We retain the corresponding full submatrix, accept that covariance as released, and do not add those errors again. We have not independently rebuilt the SH0ES covariance from the underlying Cepheid measurements.

| Regression | Age coefficient, mag/Gyr | Nominal fixed-age 95% interval |
|---|---:|---:|
| Free intercept plus local age | −0.0063 ± 0.0245 | [−0.0543, +0.0417] |
| Also local attenuation and global host mass | +0.0034 ± 0.0252 | [−0.0459, +0.0528] |

A Gaussian marginal-age-error sensitivity gives an age-only nominal profile interval [−0.0656, +0.0415] mag/Gyr. That working model lacks the joint host posterior and its interval coverage is not certified. Leaving an entire host out for prediction does not improve the error or predictive score when age is added. Even with ages treated as exact, this sample has only **23% power** to detect a −0.030 mag/Gyr coefficient at a two-sided 5% threshold. Its null result cannot rule out that value.

NGC 3147 supplies three observed sibling SNe, but only two independent contrasts. For SN 1997bq minus SN 2021hpr, the local-age difference is 2.99 Gyr and the corrected-magnitude difference is −0.0186 ± 0.0893 mag. The proposed coefficient predicts −0.0897 mag. The marginal-age errors are themselves large, and the three pairwise contrasts share objects; this is not three independent age tests.

## 3. ZTF confirms environmental prediction, then supplies a conditional age test

We fit a declared subset of the released ZTF SALT parameter tables, keeping reliable catalogue host redshifts, 0.01≤z≤0.06, public quality flags, normal/91T subtype labels, and complete local/global host properties. The result is **537 SNe**, rather than the paper's unavailable exact 932/929-object BayeSN selection. No fit is claimed to recreate that paper's result.

On these fitted observables, a jointly fitted high-minus-low global mass step is **−0.110 ± 0.020 mag** after width, colour, and a redshift nuisance term. It improves grouped held-out log predictive score by about 14.3. The sign means that the higher-mass group is brighter after width/colour standardization. A local-colour term adds a smaller and uncertain held-out improvement; allowing nonlinear width/colour terms also matters. This supports an environmental brightness association under the declared working model. It does not identify its cause as age rather than dust, metallicity, light-curve inadequacy, calibration, or a mixture.

A new primary input was located while executing the plan: the [public TITAN author catalogue at commit 9b9ed9e](https://github.com/SterlingYM/age-of-titans/tree/9b9ed9e5faf8b2f6267c0dce8cea97704d9b95ea). Its 8,610 transient rows are a **prepaper snapshot**, not the final published 6,983-host sample. Exact unique IAU aliases match 475 ZTF SNe; requiring `d_DLR<4` in both host catalogues and finite ordered summary intervals leaves **401**. The snapshot has no host coordinates, so the TITAN association itself cannot be independently sky-verified here. Its global masses correlate with the ZTF estimates at Spearman 0.980; the median scale offset is +0.155 dex, with NMAD 0.080 dex. This consistency is useful, but is not proof of identical host assignment or independent host photometry.

The comparison fits all coefficients jointly on the same 401 objects: width, colour, quadratic width/colour terms, redshift, global mass step, local colour, and TITAN continuous global mass, attenuation, and metallicity. It then adds either global mass-weighted stellar age or the catalogue's inferred mean-delay summary:

| Incremental predictor | Coefficient, mag/Gyr | Change in held-out log predictive score |
|---|---:|---:|
| Global mass-weighted host age | +0.0103 ± 0.0110 | −0.23; conditional bootstrap 95% [−2.10, +1.58] |
| SFH-derived mean progenitor delay | +0.0036 ± 0.0070 | −0.29; conditional bootstrap 95% [−1.64, +1.06] |

Five deterministic folds keep each of the 400 host groups together. The bootstrap intervals condition on the fitted overlapping training folds and omit unavailable shared calibration modes. Independent marginal host-error propagation changes these coefficients little, but it is not a substitute for the full correlated host likelihood. With c<0.1, the age coefficient is −0.0006 ± 0.0100; this is a SALT-colour restriction, **not an independently measured low-extinction sample**.

The key identification limit is quantitative: the host-age variation remaining after the controls is only **0.784 Gyr**, while the median marginal 68% half-width of the inferred ages is **0.929 Gyr**. Posterior medians cannot be assumed to be unbiased, independent noisy ages in this regime. The inferred delay is a posterior summary of the SFH-specific mean delay under an adopted delay-time model, not a measurement of the individual exploded star's age. The global TITAN age scale also differs from the local R19 and C25 age definitions, so transplanting a slope between them is an additional hypothesis.

## 4. Blinding, calibration, and numerical checks

The [ZTF release description](https://arxiv.org/html/2409.04346v2#S6.T3) explicitly marks x0 as **“blinded flux zeropoint near 30.”** A [primary companion paper](https://inspirehep.net/files/91e12bac894170dcfd70a703f09aac5c), A&A 697 A125, page 12, describes a blinding factor added to the fitted ZTF absolute-magnitude parameter M. This supports an unknown common brightness zero point, which our free intercept removes. The precise public-release transformation code was not recovered. These results therefore remain conditional internal contrasts, not fully certified unblinded luminosity tests.

We explicitly shifted all magnitudes by +0.2 mag: the intercept changed by +0.2 and the age coefficient by less than 10⁻⁸ mag/Gyr. Replacing the linear redshift trend with free 0.005-wide redshift-bin intercepts gives **+0.0125 ± 0.0110 mag/Gyr**. Applying arbitrary bin-constant magnitude shifts leaves that within-bin age coefficient unchanged below 10⁻⁸ mag/Gyr. These checks protect against the stated constant and bin-constant shifts; unrestricted object-dependent or within-bin calibration changes remain unconstrained. The release itself says DR2 is not yet calibrated for cosmological inference.

For the fixed observed design and conditional Gaussian working likelihood, 300 paired null and −0.030 mag/Gyr simulations all converge. Nominal 95% intervals cover the input in 96.0% of cases, null rejection is 4.0% (Monte Carlo standard error 1.1 percentage points), and detection of the nonzero coefficient is 80.0% (2.3 percentage points). This tests optimizer/interval behaviour under that model. It does **not** test recovery of latent physical ages, missing host posterior correlations, selection, light-curve training, calibration, or transport to high redshift.

The observations have therefore narrowed the practical issue: a blanket additional linear age term is not supported by these conditional predictions, but the assumptions required to identify an unabsorbed physical age effect remain unvalidated. The next discriminating data are joint host age/dust/SFH likelihoods with local age information, independently constrained SN extinction or optical–infrared contrasts, and a verified selected-survey recovery calculation. Increasing the number of posterior-median ages alone would not resolve those gaps.

## 5. Nearby Dovekie and optical–infrared coverage

A final descriptive coverage audit checks the available named events and sky/peak-epoch bridges, rather than treating numerical survey IDs as physical names. TITAN matches 15 uniquely named events in the Dovekie distance compilation; 12 retain finite host summaries and `d_DLR<4`. All 12 are **nearby Foundation SNe**, with 0.0256≤z≤0.0870, not the main DES survey. There are no TITAN matches to the 31,636 DES HEAD aliases, nor through the 2,691 TITAN events with usable ZTF coordinate/epoch bridges. TITAN itself supplies no coordinates, so this does not prove that every possible external alias has been resolved.

The [Dovekie release](https://github.com/des-science/DES-SN5YR/tree/main/4_DISTANCES_COVMAT) provides already standardized and bias-corrected distances, with BEAMS uncertainty treatment. Its `STAT+SYS.npz` is the **total inverse covariance**. We unpack all 1,820 rows in Hubble-diagram order, solve for the required covariance columns and then select the 12 events; using the differently ordered metadata table or selecting the precision matrix first would be incorrect. We neither add statistical errors again nor refit the published corrections. The full solve closes below 3×10⁻¹⁵, and independent whitened least squares reproduces both fits.

| Conditional Foundation regression | Age coefficient, mag/Gyr | Nominal fixed-age 95% interval | Power for −0.030 mag/Gyr |
|---|---:|---:|---:|
| Intercept plus global age | −0.0416 ± 0.0326 | [−0.1055, +0.0222] | 15.1% |
| Also global attenuation and mass | −0.0091 ± 0.0384 | [−0.0844, +0.0661] | 12.2% |

The weak negative age-only coefficient moves towards zero after host controls, but neither fit distinguishes zero from −0.030 mag/Gyr. These are fixed posterior-summary regressions on a selected 12-event sample; joint host errors, host-assignment verification and high-redshift transport are not established. These distances provide a useful independent measurement route, although the same host inference assumptions recur.

The 79 RAISIN events with common optical, infrared and combined distances contain **zero identifiable TITAN or DustPedia age matches** under direct names and the available 3-arcsec/30-day bridges. The bridges do recover 61 RAISIN events in Pantheon+, so the distance-event matching itself works. Thus the available optical–infrared distances cannot yet perform this age comparison. In addition, their separate branch covariance matrices do not provide the optical–infrared cross-covariance needed for the uncertainty of a difference. Zero cross-covariance is not assumed. The [coverage ledger](../results/environment_validation/coverage.json), [nearby distance fit](../results/environment_validation/dovekie-summary.json), and [extension validation](../results/environment_validation/coverage-validation.json) preserve these measured coverage and precision limits.
