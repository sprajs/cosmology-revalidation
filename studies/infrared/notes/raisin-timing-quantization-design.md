# A distinct quantization-aware timing experiment: design only

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../README.md) for current execution status.

The current archived timing branch remains stopped. Its prescribed precision screen fails six chi-square caps even though every tested distance response is below 0.000132 mag. That result must not be relabelled as successful historical reproduction. A different, useful target is possible: an interval for the **change in a declared distance estimator when the same photons are fitted at two fixed peak times**. This is a conditional measurement-processing response, not a timing-bias correction, calibrated likelihood ratio or cosmological result.

No JOINT-peak outcome, native timing fit or new observed-flux score was evaluated for this note. The companion executable reads synthetic vectors only.

## What the new target could identify

Let x contain the original retained fluxes, uncertainties, epoch times and required metadata. The rounded archive specifies a compatible set X, with shared x in both arms. Let t_H be the archived NIR header peak and t_J the coherent same-CID joint-fit peak. For a fully stated estimator T, the target is

    Δ(x) = T(x, t_J) − T(x, t_H),    x ∈ X.

The identifiable object is a set of responses, not the value at whichever reconstruction gives the nicest chi-square. Bounds must preserve the same flux/error/time realization in both arms; independent perturbation envelopes for the two fits would throw away this pairing. The eight original CIDs, all masks and all failed historical gates remain in the record. No new tolerance or pass/fail threshold on Δ is needed: report the response set, its numerical certification and its conditional scope.

The historical native three-iteration algorithm is one possible T, but certified global bounds for that algorithm are not yet available. Its changing covariance, clipping and parameter initialization make local derivative screens insufficient as a mathematical proof. The completed endpoint/corner screen is evidence that the sampled response is smooth, not certification over the entire input box.

A simpler, explicitly distinct T uses fixed positive covariance and the native SNooPy mean with shape=1, AV=0 and RV=1.518. This changes the statistical model relative to the iterated historical estimator; it cannot be called its repaired reproduction. A common covariance across the two peak arms gives the clearest controlled comparison. If an observed-data-derived covariance is frozen, the result is conditional on that estimated covariance, rather than a new generative validation of it.

## A conditional likelihood and an exactly solvable part

For a fixed mean shape h(t), distance amplitude a>0 and covariance C independent of a,

    y ~ N(a h(t), C),    W=C⁻¹,
    â(t) = h(t)ᵀ W y / [h(t)ᵀ W h(t)],
    D̂(t) = Dref − (2.5/ln 10) ln â(t).

The interior formula requires a positive numerator. A nonpositive numerator is a boundary case that must be retained and handled by a stated proper amplitude/distance prior; it must not be clipped into a fabricated finite distance. Covariance determinants and normalization are constant in this particular amplitude fit. They are not optional if the model changes C with a or compares probabilistic models.

For flux quantization alone, with means and covariances fixed, define v_k=W_k h_k and q_k=h_kᵀW_kh_k. The paired amplitude ratio is

    â_J/â_H = (q_H/q_J) (v_Jᵀy)/(v_Hᵀy).

If both numerators are positive throughout the flux box, the extrema of this linear-fractional function occur at box vertices. With only 2–6 retained epochs in the current eight-object pilot, all flux vertices can be enumerated exactly. A direct positivity check is the minimum linear form on the box. Each coordinate derivative of a linear-fractional function has a sign independent of that same coordinate; hence endpoints contain an extremum. This is an analytical bound, not a Monte Carlo confidence interval. The decreasing log map reverses the ratio extrema when converting to ΔD.

Time, redshift and error rounding also change h and W. The flux-only vertex proof does **not** certify those additional dimensions. A complete interval claim therefore needs either a source-based global derivative/interval bound for those inputs, or a plainly restricted conditional target that holds them fixed. There is no justification for ignoring their uncertainty merely because the flux calculation is exact. Native template interpolation-cell boundaries, mean positivity, model support and active parameter boundaries are mandatory gates. A finite grid or derivative screen may be reported as a screen, not relabelled as an enclosure.

A quantized Gaussian likelihood is another possible model:

    L(a,t) = ∫[flux bin] N(u; a h(t), C) du.

This integrates the measurement distribution over an observed bin. It is not the same as randomly imputing a uniformly distributed “true flux” and fitting each draw. Times and quoted-error precision have no automatically justified uniform prior; their uncertainty can instead remain a set of nuisance values. Moreover, accepted-epoch and object selection may be informative. Without the original selection likelihood, this is at most a conditional working likelihood for retained measurements, not a population-normalized distance likelihood. Joint peak estimates extracted from the same photons cannot be added as independent timing data in that likelihood.

## Small timing errors: amplitude, curvature and correlated noise

The following expansion assumes a correct mean family, fixed C, positive interior amplitude, a fixed epoch mask, smooth native mean within the perturbation domain, and small timing and relative flux errors. It excludes covariance updates, clipping and nuisance refits. Peak error η is in observer-frame days; native derivatives must include the (1+z) conversion and phase-dependent K correction, not only an isolated rest-frame light-curve derivative.

Write y=A h₀+n and h(η)=h₀+ηh₁+η²h₂/2. Define

    S=h₀ᵀWh₀,    a₁=(h₁ᵀWh₀)/S,
    b=(h₁ᵀWh₁)/S,    c=(h₂ᵀWh₀)/S,
    e₀=(h₀ᵀWn)/(A S),    e₁=(h₁ᵀWn)/(A S),
    K=2.5/ln 10.

Here a₁ has units day⁻¹, b and c day⁻², e₀ is dimensionless and e₁ day⁻¹. Keeping total order two in η and n gives

    D̂(η)−Dtrue ≈ K[−e₀+a₁η + ½e₀²
                         −ηe₁+a₁ηe₀
                         +(b+c/2−3a₁²/2)η²].

For the **same-photon** difference from a fit at the true peak, the ordinary amplitude-noise/Jensen terms cancel:

    ΔD ≈ K[a₁η − ηe₁+a₁ηe₀
                 +(b+c/2−3a₁²/2)η²].

Consequently its expectation includes timing mean, timing curvature and noise/timing cross moments:

    E[ΔD] ≈ K{a₁E[η]
             +(b+c/2−3a₁²/2)E[η²]
             −E[ηe₁]+a₁E[ηe₀]}.

Cross moments include covariance and products of nonzero selected-sample means. These coefficients differ by object, passband, redshift and observed phase; a single pooled timing standard deviation cannot be inserted into one mean coefficient without further assumptions. Curvature can have either sign. The leading full-distance variance also contains

    K²[Var(e₀)+a₁²Var(η)−2a₁Cov(e₀,η)].

An independently injected zero-mean Gaussian with the fitted peak standard deviation sets the timing mean and noise cross moments to zero. It therefore omits effects of using the same NIR photons to estimate the joint peak, phase/selection-dependent errors, model mismatch, non-Gaussian tails and shared calibration or intrinsic-scatter fluctuations. Matching the standard deviation alone does not match the estimator's bias or scatter. Adding that extra timing dispersion to an existing intrinsic-scatter simulation may double-count a mechanism already represented by the joint fit; the decomposition must be measured with paired draws.

For arbitrary fixed y, without assuming the correct mean or Gaussian noise, a useful exact local identity is available. Put B=hᵀWy, S=hᵀWh, T=h₁ᵀWh, U=h₁ᵀWh₁ and V=h₂ᵀWh. Then

    D′ = K[2T/S − (h₁ᵀWy)/B],
    D″ = K[2(U+V)/S − 4T²/S²
            − (h₂ᵀWy)/B + ((h₁ᵀWy)/B)²].

These are derivatives of the specified fixed-C estimator, not an average bias formula. The archived simulation actually contains varying truth stretch/dust while the NIR estimator fixes its nuisance parameters, so the well-specified y=A h₀+n expansion must not silently be applied to that population. Prospective source-consistent simulations can separately vary this misspecification. For the historical finite-iteration algorithm, derivatives must follow the actual sequence of covariance updates; the fixed-C formula cannot substitute for that sequence.

## Recommended next bounded decision

1. Keep the failed historical precision branch closed and retain its artifacts. Do not replace its inputs with the individually successful SIMLIB branch.
2. Prioritize the separately proposed fresh source-consistent paired native simulation with known full-precision photons, same latent/noise draw in both fits, and a peak estimated from those same photons. Its eight engineering and 64 conditional draws can measure the joint cross moments that independent timing smearing discards. It still does not identify a real-population or cosmological correction.
3. If extracting more from the coherent archive is worthwhile, first freeze a **new conditional interval-response design**, with a declared estimator, shared uncertainty set and a numerical enclosure objective rather than an arbitrary relaxed chi-square cap. A fixed-C flux-vertex proof is immediately feasible. A joint time/error/redshift enclosure remains the main technical prerequisite for a full quantization-aware interval. No JOINT-peak score is authorized by this design note.
4. Only a subsequent execution-linked full-precision dataset and selection/noise validation could support an actual bias-corrected population distance likelihood. Current archive responses and synthetic algebra cannot supply those missing measurements.

## Synthetic verification and provenance

`runs/research_2026_09_26/astra_design/raisin_timing_assets/quantization_design/algebra_check.py` uses invented vectors only. Finite-difference checks agree with the exact first/second derivatives within 4.15×10⁻¹¹ and 7.87×10⁻¹⁰; halving η reduces the noiseless second-order remainder by approximately eight. Enumerated flux-box extrema contain all 10,000 independently sampled interior points, alongside the analytical endpoint proof above. An illustrative correlated-noise toy changes the sign of the expected timing increment relative to independent smearing; it is not calibrated to any SN.

The source/version and precision claims are grounded in [the acquisition and failed precision report](raisin-timing-asset-recovery.md), its immutable source excerpts and the exact official v11_04d commit. The statistical derivations here are new algebra under explicitly stated assumptions. No new observed outcomes or native fits were evaluated.
