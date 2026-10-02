# Conditional SN observer and historical optical control

The question is whether native supplied-redshift observer geometry and sampled
optical photometry preserve their source identities and agree with independent
mathematical controls. This is a bounded deterministic experiment, with a
synthetic coasting background and steady flat spectrum. It does not reproduce
the Pantheon+ fit or provide a measured instrument calibration distribution.

The pinned Prospector candidates retain the broader scientific blockers. The
selected Pantheon+ sample is `IDSURVEY=1`: 321 light-curve rows in original table
order. Every row has a unique source SDSS header/update join and valid photometry
pointer interval, covering 35,181 calibrated points. The header-named and renamed
PHOT bodies have distinct byte identities; audited MJD/FLT/FIELD/FLUXCAL arrays
agree. Ten RA comparisons exceed final-table decimal rounding, updated velocity
map/group ancestry remains unresolved, and the PHOT asset alone does not define
physical flux units, time scale or an error law. The full original covariance is
mandatory for a future qualified likelihood; no covariance or magnitude fit is
performed here. Distinct CIDs are not a physical-event deduplication or proof of
independence from other surveys.

The released CosmoSIS code uses `(1+zHD)(1+zHEL) D_A(zHD)`, equivalent in the
declared flat background to `D_L=(1+zHEL)(c/H0) integral_0^zHD dz/E(z)`. The native
control chooses `ConstantQ(q=0)`, `E(z)=1+z`, `H0=70 km/s/Mpc`; it compares all
321 released pairs and 321 same-redshift twins, zero and two low-redshift limits.
This analytic model is a numerical control, not the fitted SN cosmology or a
reconstruction of the frame/velocity corrections.

The response is STScI `wfc3_uvis_f555w_001_syn.fits`: a historical ground optical
filter component with 910 binary32 samples on 1938–11002 Angstrom, zero endpoints
and peak transmission 0.9519490003585815. Amplitudes are unaltered. Wavelengths
use the native exact `1e-10 metre/Angstrom` conversion with final binary64 rounding;
the admitted function is linear in wavelength and zero outside its finite support.
It excludes telescope optics, detector quantum efficiency and electronic gain.
The later `_004` asset is a distinct calibrated curve: its history includes a
1.0059 normalization correction and its equal time columns are not random draws.
Neither is substituted for the SDSS observing system.

The source spectrum is chosen constant `L_lambda=1 W/m` on `[1e-7,1.2e-6] m`.
The operator uses `F_lambda=L_lambda/[4 pi D_L^2(1+z)]`, distance `1e20 m`, area
`1 m²` and observer exposure `1 s`. Incident band flux has no transmission
weight; energy uses transmission once; photons additionally use observed
`lambda/(hc)` once. Three redshifts, 0, 0.1 and the first selected SDSS row's
heliocentric 0.06707, exercise these conventions. Exact wavelength-linear moments
and frequency-coordinate Gauss–Legendre 8/16 refinement are independent controls.
The allocations are frozen in [config.json](config.json) before code/execution.

The actual released Fragilistic archive supplies a 102-coordinate covariance of
fitted magnitude offsets with literal source label order, including SDSS griz
at indices 34–37 and no HST coordinates. The paper describes 105 fitted filters;
the missing-coordinate mapping remains unresolved. Tiny source asymmetry is
retained rather than silently symmetrized. Nine published SALT2 calibration
variants carry systematic scale 0.3, not normalized optical-state probability
masses. A joint optical transmission law requires source-supported states,
weights and calibration mapping; this experiment fabricates none.

The dedicated [controller](controller.py) verifies the candidates, frozen source
bytes, complete SDK inventory, pinned source blobs, compiler and library before
and after the bounded run. It reuses the reviewed released-ladder SDK identity
checker. Production geometry/photometry stay in installed C++ kernels; Python
only transports structural inputs and runs the separately identified reference.
Full-paper readings remain pending, with zero new full reads. Raw sources and
attempts remain ignored; source redistribution terms are unverified.

For an admitted local source store and an isolated reference environment:

```sh
python experiments/sn-observer-passband/controller.py \
  --engine-source /path/to/irreducible --sdk /path/to/preserved/sdk \
  --sources /path/to/pinned/source/store --reference-python /path/to/reference/python \
  --output results/sn-observer-passband/new-immutable-attempt
```

The preserved SDK is commit `c9b7febc01ecc85d289056faba8722a8f24a6ad8`, build ID
`a08f62cb32a76097e68beff6a540ccf4da903b5b898ae75921ee53818c913903`.
This packet does not update previous experiment SDK pins or receipts. The general
runner remains blocked for the faithful candidate; the dedicated controller
executes only the stated conditional numerical slice.

No accepted run is claimed until its immutable receipt has been inspected.
