# Dependence between the cosmological probes

The contemporary CMB, BAO and supernova combination is a **joint measurement conditional on a factorized likelihood**. Some shared information is explicitly accounted for; other cross-covariances are neglected using published approximations. The latter have not been validated quantitatively for this exact CPL dark-energy combination. This does not establish a large bias, but it prevents describing its uncertainty or significance as fully validated against all probe dependence.

This review covers the default `official_planck`, ACT maximum multipole 6500 configuration in [modern_adapter.py](../code/external_probes/modern_adapter.py). It examines measurement dependence, not convergence, the adequacy of the cosmological model, or supernova correction identification. The [machine-readable audit](../results/inference/probe-dependence.json) records source snapshots, executable runtime checks, and their exact scope.

## What the current likelihood includes

The CMB blocks retain their released internal TT/TE/EE or lensing covariance. ACT and Planck lensing are one joint block: its 19-by-19 covariance contains ten ACT and nine Planck band powers. Independent inspection gives maximum absolute cross-experiment correlation **0.206284**. There is no second standalone Planck lensing likelihood. This joint covariance must not be replaced with independent factors. The configured ACT lensing likelihood also applies its released theory-response corrections.

The Planck primary likelihood discards high multipoles; 252 bins remain active and 361 receive effectively infinite variance. ACT begins at multipole 600 in TT, TE and EE. The crop removes complete Planck bins above TT 1000 and TE/EE 600; it does **not** eliminate the TT 600–1000 overlap. Louis et al. explicitly neglect this residual covariance because Planck has smaller errors there. Their example of less than 5% improvement in the uncertainty of the effective number of neutrino species from a more optimal combination is not a bound on CPL uncertainty. [Louis et al., section VI.3](https://arxiv.org/html/2503.14452v2#S6.S3)

`A_act = A_planck` shares the ACT and Planck temperature-amplitude parameter. SPT temperature and polarization calibration retain separate scalar priors. This is partial representation of common calibration information, not a complete cross-experiment calibration covariance. Sharing a theoretical spectrum also does not supply the sampling covariance of two measurements of the same sky.

## Which remaining CMB correlations are approximated

| Measurements | Treatment in this adapter | Evidence and scope |
|---|---|---|
| ACT primary and cropped Planck primary | Cross-covariance omitted after scale cuts | Residual TT overlap explicitly accepted by Louis et al.; the multipole ranges are not completely disjoint. |
| SPT D1 primary and ACT/Planck primary | Separate likelihood factors | Camphuis et al. assume independence; the SPT footprint overlaps about 10% of ACT's mask. This is a documented approximation. |
| Pan 2023 SPT lensing and primary CMB | Separate likelihood factors, with lensing response corrections | Pan et al. neglect the two-point/four-point covariance, citing small cosmological effects at then-current noise. Response corrections alter the predicted mean, not the omitted sampling cross-covariance. |
| Pan 2023 SPT lensing and ACT/Planck lensing | Separate likelihood factors | Pan et al. justify approximate independence from Planck using differing sensitivity/noise and sky fractions, 3.6% versus 67%. No matching Pan–ACT cross-covariance is included here. |

Sources: [Camphuis et al., section VII.1](https://arxiv.org/html/2506.20707v2#S7.SS1); [Pan et al., sections VI.A.1–2, printed pages 17–18](https://arxiv.org/pdf/2308.11608).

The SPT release name `SPT3G_2018_Lens_and_CMB` is easy to misread. Its data vector here contains **12 lensing band powers**, with primary-CMB response corrections. It does not add a second primary TT/TE/EE measurement.

Newer SPT results are useful evidence but are not interchangeable with this lensing block. Camphuis et al.'s simulated primary/lensing independence check uses **MUSE polarization lensing**. The ACT/SPT/Planck joint-lensing analysis of Qu et al. likewise uses the newer SPT polarization estimator: its neighboring ACT–SPT lensing bins have correlations around 10–20%, with small posterior effects even in a pessimistic covariance test. This does not measure the omitted covariance for the **Pan 2018 temperature-based** reconstruction used here. [Camphuis et al., section VII.2.2](https://arxiv.org/html/2506.20707v2#S7.SS2.SSS2), [Qu et al., section VIII](https://arxiv.org/html/2504.20038v2#S8)

The [released joint-lensing alternative](joint-lensing-availability.md) has now been recovered and independently evaluated at the supplied fiducial spectra. It includes ACT–Planck–SPT MUSE cross blocks and is available for a separately identified posterior comparison. Using it replaces both current lensing factors and changes the SPT estimator; it is not just a covariance switch for the current data. It still supplies no primary–lensing, lensing–BAO or lensing–SN cross block.

The model extension matters. Trendafilova's idealized primary/lensing forecast finds changes of at most about 3% in six-parameter ΛCDM uncertainties, but an effect around 5% for free constant `w0` even at 10 μK-arcmin noise with minimum lensing multipole 50. The forecast assumes half-sky coverage, a 1.4-arcminute beam and simplified noise/foreground conditions; it is not this survey combination. Covariance can increase or decrease an uncertainty. These numbers therefore neither supply a correction factor nor bound the CPL `w0, wa` significance here. [Trendafilova, sections III–IV and figure 3](https://arxiv.org/pdf/2308.11588)

## A separate published prior-counting problem

The official SPT release notice dated April 1, 2026 reports that some CMB-SPA and CMB-SPA+DESI chains counted Gaussian Planck nuisance priors repeatedly when cropping Planck's high-multipole likelihood. Corrected chains change parameter values by approximately 0.2–0.3 standard deviations. This is a prior-construction error, not a cross-covariance effect. It does not by itself establish that a different author's chain is affected. [Official LAMBDA notice](https://lambda.gsfc.nasa.gov/product/spt/spt3g_d1_bandp_liklyhood_info.html)

Independent initialization of **our current adapter** found one cropped Planck instance and no remaining internal priors in either SPT candl likelihood. Moving each of `A_planck`, `P_act`, `Tcal` and `Ecal` from its prior mean by one declared standard deviation changed the global log prior by −0.5, within 3.6×10⁻¹⁴. This rules out multiplicity in those configured scalar prior factors; it is not a validation of every historical chain or all information embedded in released likelihoods.

## Supernovae and BAO

The adapter uses DESI DR2 BAO and a single supernova compilation at a time. It does not multiply overlapping supernova compilations or add DES Y6 BAO to DESI. DESI's published cosmological analysis also presents separate SN-compilation combinations. Its test correlating BAO modelling systematic contributions by 0.5 across redshift bins concerns the **internal BAO systematic budget**, not SN–BAO independence. [DESI DR2, section IV.3 and appendix B](https://arxiv.org/html/2503.14738v2#S4.SS3)

The DES supernova-plus-BAO analysis explicitly assumes independent dataset likelihoods. It supplies precedent for this factorization, not a measured cross-covariance bound. This review located no quantitative bound for the exact Dovekie–DESI DR2 or CMB-lensing–DESI BAO cross-covariance in the primary sources examined. Their omission remains an unquantified assumption at the precision of this fit. [DES SN+BAO analysis, section III.2](https://arxiv.org/html/2503.06712v2#S3.SS2)

A guardrail for future additions: the DES multiprobe analysis removes roughly 1000 square degrees of DES–DESI BAO overlap and excludes CMB lensing when including its galaxy/shear three-by-two-point likelihood, because the latter would require cross-covariance. Those extra likelihoods are not included here. Their treatment must not be used as evidence that the current SN–BAO combination has the same overlap problem. [DES multiprobe analysis, section II.2](https://arxiv.org/html/2605.27221v1#S2.SS2)

## Consequence for interpretation

Report the resulting posterior as conditional on the declared covariance factorization. A numerical preference for evolving dark energy is not yet a cross-covariance-validated significance for a complete unified measurement. No source reviewed here establishes that neglected correlations reverse that preference, and no defensible universal uncertainty inflation follows from these sources.

An exact assessment requires matched cross-covariance or common-sky simulations for these estimators, masks and noise levels, propagated through the same CPL fit. Omitting individual overlapping blocks can provide a useful sensitivity test, but a difference between such fits is not itself an estimate of missing covariance. This review makes no likelihood changes and performs no additional cosmological fit.
