# Independent cross-review of the low-RV finding

The low-RV finding in [snana-implementation.md](snana-implementation.md) survives an independent, bounded numerical check. The claim supported is that the **extrapolated host-extinction component in these simulations can exceed unit transmission at red wavelengths**. It is not a measurement of negative extinction in observed supernovae, nor a measured residual distance or cosmological bias.

[matrix_mock_crosscheck.py](../../scripts/salt_dust_audit/matrix_mock_crosscheck.py) opens the original released FITS truth tables directly. It does not import the source agent's extinction wrapper, compiled library, mock-counting script, or flux-response model constructor. For the historical law it independently transcribes the optical O'Donnell polynomial and historical F99/O'Donnell multiplier from `phase2/official/build/SNANA-2fe0f56/src/MWgaldust.c`. For exact F99 it calls the separately implemented `extinction.fitzpatrick99` package. Source files and package bytes are hashed in [the result](../../runs/salt_dust_audit/matrix_mock_crossreview/results.json).

The first, middle, and last of the 25 released mock realizations were checked, covering 8,531 written simulated objects:

| Realization | Written objects | RV < 1 | RV < 2 | Historical A(8000 Å) < 0 | Exact-F99 A(8000 Å) < 0 |
|---|---:|---:|---:|---:|---:|
| 0001 | 2,829 | 419 | 1,564 | 159 | 161 |
| 0013 | 2,843 | 474 | 1,615 | 184 | 193 |
| 0025 | 2,859 | 443 | 1,575 | 161 | 170 |

Every count agrees with the corresponding entry in the source agent's independently produced per-realization table. This cross-review is a three-realization spot-check, **not an independent recount of all 71,946 objects**. The denominator is written simulated objects after generation/write selection; it does not represent all attempted transients, the final BBC sample, or a frequency of real interstellar dust properties.

For example, realization 0025, simulated SNID `226984`, stores `SIM_RV=0.3056366444` and positive `SIM_AV=0.1012773439`. Direct evaluation gives historical `A8000=−0.0814688514` mag and exact-F99 `−0.0828103216` mag. The historical host component has multiplier `10^(−0.4 A)=1.07792251`. This is one actual stored simulation draw, not an invented RV/AV pair. Other source components and noise still determine the final simulated observations.

The independent zero-extinction thresholds at 8000 Å are RV=0.65323401588 for the historical approximation and 0.66338573930 for exact F99, agreeing with the source audit. The historical `genmag_SEDtools.c::fill_TABLE_HOSTXT_SEDMODEL` evaluates `GALextinct` at `LAMOBS/(1+z)` and stores `pow(10,−0.4*XT_MAG)`. The historical SALT integrator reads this host multiplier. The source path does not transform a negative monochromatic extinction into a positive attenuation. This confirms the narrow interpretation of the component; it does not reproduce broadband fitted distances.

A separate constructed broadband check independently reads the released phase-zero SALT3 M0 surface and DES i transmission curve and directly integrates the photon-weighted spectrum, using cubic SED interpolation and linear filter interpolation at 0.25 Å spacing. It does not use sncosmo. For x1=c=0, z=0.1, E(B−V)=0.1 and RV=0.4, the result is `A_i=−0.0215412274` mag, corresponding to a flux ratio of 1.02003833. This agrees with [the root's constructed example](../../runs/salt_dust_audit/f99_support/summary.json), `−0.0215411921` mag, within 3.54×10⁻⁸ mag. The monochromatic finding therefore can survive integration in at least this particular passband/configuration. This example is not asserted to be one of the actual stored mock events, and its grid is not weighted by the selected DES population.

The implementation finding is robust enough to justify a targeted closure test and to reject a literal passive-screen dust interpretation of these particular extrapolated values. It does not establish the sign or size of a cosmological bias: the model may act phenomenologically, and SALT fitting, intrinsic colour, selection and BBC may absorb or transform the effect. Physically admissible alternatives need matched population generation, light-curve fitting, selection and bias-correction reruns before assigning a distance correction or uncertainty. Simply clipping RV would change the population and is not a validated correction.

Reproduce with:

```sh
OPENBLAS_NUM_THREADS=2 phase2/env-official/bin/python scripts/salt_dust_audit/matrix_mock_crosscheck.py
```
