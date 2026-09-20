# Standardization, calibration and dust experiment register

Investigator: independent standardization/dust branch. Created 2026-09-20 22:03 UTC, before numerical outcomes. No other investigator conclusions have been read. Primary-source claims are hypotheses, not truth labels.

## STD-01: release identity, correction accounting and controlled DES comparison

Question: how much does replacing the original DES-SN5YR distance product with DES-Dovekie change identical objects, and can the observed shift be assigned uniquely to calibration or dust?

Inputs: original DES repository tag 1.3 (e3493cb3…), Dovekie c9a4fca…, Pantheon+ W26 snapshot 7fc6805. Use HD row order for matrices and explicit (CID, IDSURVEY) joins for metadata and matching. Preserve inputs. Record hashes, duplicate counts, sample overlap, redshift agreement and metadata differences.

Method: verify original systematic-only covariance plus MUERR_FINAL squared diagonal, and reconstruct total Dovekie precision from packed upper triangle. Test symmetry, Cholesky positivity, diagonal agreement for statistical matrices. Calculate delta mu (Dovekie minus original), and decomposition into x0 amplitude, alpha*x1, -beta*c, host step, -biasCor_mu and intercept from release metadata/global parameters. Report residual if documentation rounding or hidden normalization prevents exact closure. Do not attribute jointly changed terms to a single physical cause. Summarize matched shifts by survey and fixed redshift bins (0,.1,.3,.5,.7,1,2); common offset is arbitrary.

Controlled model implication: fit uncalibrated flat LCDM Omega_m with an analytic free magnitude intercept, using each native total covariance and optionally the identical matched sample with covariance obtained by subsetting the covariance (not subsetting precision). Cross native distance vectors and covariances on the matched set to separate distance and weight effects. These are corrected-distance analyses, not raw photometry reconstructions or reproductions of combined CMB+BAO w0wa results. No significance for release differences will be assigned without cross-release covariance. Examine identical IDs shared with Pantheon+ as an independence audit (unresolved aliases mean lower-bound overlap).

Prediction/criterion: an overall magnitude shift has no SN-only cosmological effect; a redshift/survey-dependent shift can affect inferred curvature of distance-redshift relation. Different standardization/training/calibration/selection changes mean total release difference cannot identify dust error alone. Significant covariance mishandling invalidates any affected fit before interpretation.

## STD-02: attenuation versus extinction and selection identifiability

Question: does a different host-mass trend of galaxy-wide effective attenuation RV logically falsify an opposite trend of SN sightline extinction RV? What additional data are required?

Method: a physically explicit absorption-only mixed slab with homogeneous emissivity and dust has transmission T(tau)=(1-exp(-tau))/tau. Compute integrated B,V attenuation and effective RV=A_V/(A_B-A_V), fixing microscopic extinction tau_B/tau_V=1+1/RV. Compare to foreground-screen extinction. Calculate grid tau_V=.01,.1,.3,1,3,10 with microscopic RV=2,3.1,4. Plot monotonic optical-depth effects and an explicitly illustrative pair of host populations with inverse microscopic RV and different optical depths. State excluded scattering and age-dependent emissivity; this is a counterexample/identifiability exercise, not fitted galaxy astrophysics.

For SNe uniformly distributed in slab depth, each sightline keeps the microscopic extinction RV. Impose a deterministic extinction detection limit A_B<.5 mag and compare selected mean reddening with unselected. Show selection can change observable dust-column distribution and suppress information about dustiest objects without equating SN extinction and integrated galaxy attenuation.

Then derive colour/brightness mixture c=c_int+E and r=beta_int*c_int+R_B*E+b_A*A+noise. With independent latent components, beta_eff=[beta_int Var(c_int)+R_B Var(E)]/[Var(c_int)+Var(E)] and demonstrate two physically distinct decompositions yielding the same (co)variances. Age-correlated dust and progenitor luminosity effects are not identified from mean residual versus age plus c alone without justified latent priors or additional spectral/IR/environment measurements. No synthetic fit is evidence for either real-population hypothesis.

Diagnostics: screen and optically thin slab limits recover microscopic RV; mixed-slab result should approach grey attenuation with increasing optical depth; numerical integration of uniform-depth point sources matches analytic slab transmission. Selection and intrinsic colour assumptions are exposed in outputs. Seed any simulations and store configs/manifests.

## Decision log

- 2026-09-20: use released corrected distances rather than claiming full light-curve reconstruction. Preregister exact/approximate closure diagnostic because Dovekie README nuisance coefficients are rounded; classify any failed reconstruction explicitly.
- 2026-09-20 22:13 UTC: STD-01 schema inspection confirms that the Dovekie `.csv` files are actually whitespace SNANA tables. The original release has host masses rounded to two decimals; reconstructing a sigmoid at a boundary from these cannot close to rounding precision. Retain the failed boundary reconstructions and quantify closure away from a 0.02-dex boundary neighborhood. The released distances, not these reconstructions, remain inference inputs. Report both profiled chi-square and release-style chi-square with the magnitude-marginalization normalization.

## STD-03: newly discovered low-redshift Pantheon+ mass revision

Registered 2026-09-20 22:15 UTC, before inspecting the numerical reconstructed table. Live revision search identified Roy Choudhury 2607.24443v2 (15 September 2026), which implements Hoyt et al. 2601.19424 Appendix F. This was absent the initial dossier. Preserve the original baseline and create a named sensitivity branch only.

Question: can the public reconstruction of the low-z host-mass correction be reproduced, and precisely which measurements and uncertainties change?

Method: pin the public reconstruction repository at c5583379eb06f9a4b045353bb39a0991f8c15e4c, inspect code before running, compare supplied original data with W26 snapshot and current Pantheon release. Independently reconstruct two unweighted degree-3 colour-to-bias fits on .01<=zCMB<.15 split at logM=10; apply f_low(c)-f_high(c) to m_b_corr for 9.4<=logM<10 in that interval. Match by (CID,IDSURVEY), preserve exact row order; compare published and reconstructed values, counts, redshifts and unchanged columns. Record maximum discrepancies. Root investigator will handle SN-only cosmological sensitivity using the original covariance as a conditional approximation, with no claim that the unchanged covariance captures uncertainty in the +0.6-dex mass shift or branch transfer.

Criterion: success establishes reproduction of the stated deterministic branch-transfer prescription only. It does not reproduce raw photometry/host masses, prove correctness of all shifted masses, or independently validate the publication's cosmological conclusions. A material mismatch of input release invalidates direct substitution.
