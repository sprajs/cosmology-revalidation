# Preregistered cosmology investigations

Recorded 2026-09-20 before evaluating new likelihoods. Existing papers and reported parameter values were read in advance. The acquisition dossier is baseline commit `472f13f`. These are public distance-product analyses, not light-curve refits.

## COS-01: baseline and conditional age-correction response

Question: can the public Pantheon+ likelihood recover its published flat-LambdaCDM constraint and the unmodified Son wCDM/CPL SN-only constraints, and how does a specified evolution correction alter them?

Inputs: exact public Pantheon+ distance vector and STAT+SYS covariance; preserve original ordering and repeated light-curve rows. Use zHD>0.01, zHEL in the luminosity-distance prefactor, and analytically marginalize a common magnitude offset. Do not add the tabulated diagonal errors to an already total covariance. No Cepheid distance information. Hash/check the compact Cobaya copy against the full source release. Independently implement FLRW distances and compare against Astropy and numerical quadrature. Report release differences.

Models: flat LambdaCDM, flat constant-w, flat CPL. Use DESI-style broad priors identified in the archived YAML where appropriate, disclose all actual ranges. First validate uncorrected likelihood. Corrections are explicitly labelled: an equation-level CSFH/DTD reconstruction (if defensible), and independent sensitivity templates. Never call the public Sah correction Son's exact array. Shared slope uncertainty is correlated across redshift, not independent per SN; compare fixed slope with Gaussian amplitude uncertainty and a broader sensitivity prior.

Predictions: a positive dimming/evolution term subtracted from high-z SN distances will generally weaken the inferred acceleration. Recovery of that shift alone does not validate its astrophysical premise. A residual after fitting population or calibration nuisance terms is conditional on those terms. Reproduction tolerance: published posterior centres within 0.2 quoted standard deviations for same inputs/priors; failures investigated and labelled approximate/blocked, never hidden.

Diagnostics: positive covariance, sample/ID counts, full versus accelerated likelihood differences, analytic distance limits, offset invariance, optimizer starts, effective sample size/autocorrelation, independent ensemble runs where consequential, posterior predictions and prior-boundary occupancy. Derived q0 must be evaluated per posterior sample, not only at marginal parameter means. Separate profile likelihood ratios from posterior probabilities and from significance against LambdaCDM.

## COS-02: kinematics and identifiability

Question: which evidence concerns an accelerating epoch, which concerns present acceleration, and what survives freedom in luminosity evolution?

Fit a smooth kinematic history q(z)=q0+q1 z/(1+z) and a piecewise-constant q history over predeclared z intervals [0,0.1,0.3,0.6,1.0,2.5], with continuity of H and distance. Test adjacent-bin smoothness scales 0.5 and 1.5 and broad bounds [-3,2]; the first bin is a finite-redshift average, not a point measurement of q(0). Compare constrained all-nonnegative-q maximum likelihood with unconstrained histories, without using an uncalibrated chi-square asymptotic law for the composite boundary null. These models assume flat homogeneous/isotropic metric expansion, transparent propagation/distance duality and a standardizable luminosity population, but do not require GR dark energy dynamics.

Before fitting, derive the exact additive degeneracy between an arbitrary redshift-dependent standardized absolute magnitude and luminosity distance. Numerically demonstrate two different expansion histories giving identical predicted apparent magnitudes under the compensating evolution. This is an identifiability result, not evidence that the compensating function is physically realized. Compare its scale to age/dust/calibration sensitivities without treating visual similarity as model selection.

## COS-03: BAO and original chains

Independently evaluate DESI DR2's 13-component BAO Gaussian likelihood with a free H0*r_d scale. Fit BAO alone and BAO+Pantheon+ within CPL, retaining correlations and reporting probe-specific contributions at the joint optimum. No compressed CMB substitute. If exact Son full-CMB reproduction is blocked by author-specific correction, priors or likelihood configuration, state that barrier. Derive q(z) and P(q0<0) from the public DESI reference chains as a reference calculation only, not as a reproduction of the original inference.

## COS-04: directional branch numerical audit

Inspect public Sah C1/C2 scripts, row/covariance ordering, redshift frames and likelihood normalization. Test cubic luminosity-distance Taylor truncation against exact known cosmologies through z=0.8, including mock fits and low-z restrictions if material. Reproduce a tractable published fit only with correct input/covariance. Missing C2 products/configuration must be resolved publicly or named precisely. Directional and isotropic analyses share Pantheon+ observations and age assumptions; do not combine as independent measurements.

All changes to this plan will be dated in the decision log before their new outcomes are inspected. Synthetic experiments will be labelled as such. A separate Astra Ultra falsification audit is required before the final scientific report.
