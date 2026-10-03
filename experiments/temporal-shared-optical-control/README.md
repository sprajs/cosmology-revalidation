# Shared optical state across evolving exposures

Does one supplied optical state shared across two exposures and two bands change
the joint likelihood compared with separately mixing each exposure? This packet
uses a fixed synthetic source and two-state law to test the compiled temporal
photon and detector interfaces. It does not reconstruct measured supernova
photometry or qualify a released likelihood.

The original engine-owned control came from component revision
`bd2dff2d7d84fec45e384326d8e0e9d2aa1e9d09`, ignored source directory
`evidence/next-wave-temporal-optical-sdk/` in Irreducible. That commit alone does
not identify the ignored files. Their exact source bytes and history are retained
in `.work/temporal-optical-admission-20261003/`, admission SHA256
`1ee77abd937be8342f3be946c17453a162f75e9cd315a9a53b05f8fe2db401cc`.
The original v2 contract has SHA256
`e3a1ceb57d96f599cf82d938462d03feef5d424ad9d2c70f82316d63c34ba6d7`.
All original source and attempt identities, including the missing-mpmath and
missing-time-wrapper failures, remain unchanged. The earlier engine-owned
attempt reported 169 passing comparisons and a signed -3.8038735 percent
zero-read-noise all-detected numerator difference (magnitude 3.8039 percent).
That is prior handoff evidence; this public controller
needs a fresh committed Reproducible attempt under its own request/source identity.

The source is the original bilinear spectral/time law
`L_lambda = 2^-80 * (1 + u/4 + x/2 + u*x/16)` W/m,
where `u = tau/(1 second)` and `x = (lambda_rest - 1 metre)/(1 metre)`.
The rest-frame time support is -1 to 2 seconds and wavelength support 1 to 2
metres. Redshift is 1, so `tau = (t_observer - 10 seconds)/2`; luminosity distance
is 1 metre and collecting area is 1 square metre.
These are declared controls, with no measured source, fitted calibration offset
or SALT2 template substituted. The source is exactly zero outside its declared
support. The same native immutable source and four state/band curves are
prepared once; only transmitted photons are requested.

There are two observer exposures with full durations of 2 and 3 seconds.
Both have 2 seconds of source overlap, so the second has partial coverage 2/3.
Photon integration uses the supported part of the source; normalization and dark
current retain the full exposure duration. No extra coverage multiplier or
covered-duration renormalization is introduced. The nonseparable `u*x/16` term
changes the spectral shape across exposures.

One latent state has masses 1/4 and 3/4. Its two band transmissions are shared
across all four readouts. Conditional Poisson arrivals and read errors are
independent across the readouts given that state; this is a synthetic assumption.
The state mixture occurs once after multiplying the four conditional likelihoods.
The selected all-detected likelihood uses the same common-state ALL-four event
probability as its denominator. Count, continuous ADU density, selected and
censored records are alternatives, not extra events to multiply together.
Nondetections, structural zeros and required channel failures remain explicit.

Independent references use the original polynomial antiderivatives, observer
time/frequency quadrature and Poisson/characteristic-function density/CDF routes.
Their refinement, arithmetic reservations and tail errors stay inside the fixed
per-quantity allocations. Those routes share the declared synthetic source and
law; method agreement does not establish observational independence. The
empirical photon endpoint allowance is not a universal temporal certificate.
These exact 2/3-second durations do not qualify arbitrary-duration dark-current
sensitivity.

The dedicated controller admits clean committed source, the original contract,
new request and complete pinned SDK/toolchain/runtime identities, then compiles
one native consumer and runs it and the reference sequentially. Typed complete
and failed-control groups are admitted before comparison. Captured log bytes are
bound to parsed objects; before/after source, SDK, runtime and output identities
are retained on failures as well as success. Fresh bounded result directories
cannot overwrite previous attempts and are sealed read-only at completion.
The generic one-request runner does not execute this blocked packet.

Run the dedicated route from a clean committed checkout with the read-only SDK
and reference runtime at their reviewed pinned paths:

```sh
uv run python -B experiments/temporal-shared-optical-control/controller.py \
  --attempt fresh-shared-optical-control
```

The attempt name is one new bounded ASCII component. The controller admits no
numerical overrides or arbitrary engine commands. A fresh attempt lives under
`results/temporal-shared-optical-control/`; complete local records, original
consumed buffers and comparisons stay there. A changed live output revokes
acceptance while retaining the originally consumed value and earned checks.

The supplied optical state law is not a measured transmission distribution.
Measured template/support/training, optical calibration, raw electron/extraction,
parent population/selection and cross-observation covariance remain open. Fitted
magnitude offsets and covariance cannot supply those missing laws. No physical
calibration/systematic uncertainty or full supernova likelihood is qualified.
Production temporal/detector mathematics belongs to Irreducible; a convenience
native composition remains optional and separately reviewed.

Full local receipts and generated products stay under ignored results storage.
A durable archive of exact sources, lawful input routes and full receipts is
required before scientific publication.
