# Signed-epoch response in all ten DES16 objects

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../README.md) for current execution status.

All ten source-defined objects pass the fixed-covariance native-profile
numerical gates and the separate independent review. Adding 33 measured
negative epochs changes the accepted sample from 556 to 589 measurements.
The point responses are conditional on fixed header timing, an empirical
SNooPy model and two stated covariance anchors. They are not estimates of
population bias or corrections to cosmology.

| Object | Added negative epochs | Signed-minus-positive distance minimum (mag) |
|---|---:|---:|
| DES16C1cim | 9 | −0.164949 |
| DES16C3cmy | 1 | −0.002569 |
| DES16E1dcx | 1 | −0.000015 |
| DES16E2clk | 0 | 0 exactly |
| DES16E2cqq | 0 | 0 exactly |
| DES16S1agd | 5 | −0.049376 |
| DES16S1bno | 6 | −0.019045 |
| DES16S2afz | 2 | −0.001431 |
| DES16X3cry | 4 | −0.004753 |
| DES16X3zd | 5 | −0.003156 |

The largest response switches between almost equally good shape/colour
solutions. For DES16C1cim, the nearest high-shape signed branch instead has
a +0.01417mag response and costs only 0.11574 in the same signed quadratic.
Its descriptive ΔQ≤1 distance envelope spans 42.471–42.866mag after allowing
amplitude freedom. These envelopes are not confidence intervals. The
[exact decomposition](raisin-profile-decomposition.md) explains why the
minimum's change cannot be interpreted as a unique measurement correction.

The second covariance anchor changes paired point responses by at most
0.000690mag. This checks one state-freezing choice; it does not validate
the observational covariance or full parameter-dependent likelihood.
No ΔQ≤1 set reaches a parameter-box edge, but wider sets do: DES16E1dcx
reaches a shape edge by ΔQ4, and five objects have a signed ΔQ9 shape-edge
flag. Extrapolation beyond the tested model domain is not supplied.

The [executed cohort report and figure](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/astra_design/raisin_signed_refit/fixed-c-profile/cohort10/report.md)
retain both anchors, all branches, boundaries, failed initialization and
immutable amendments. The [independent certification](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_profile_solver_review/cohort10-review/final-certification.json)
verifies 560 final files, 14 executed sources, 670,988 saved mean vectors,
712,924 native calls and 20,958 actual-amplitude checks. The largest
independent Q discrepancy is 5.46e−12 and distance discrepancy 7.11e−15mag;
the independently verified coarse/fine minima meet the frozen 0.001 gates.
Both identical-row controls have exactly identical paired metric arrays.

The source now establishes a sign-dependent omission and a measurable
conditional fitting response. The remaining scientific task is to model
the actual joint measurement, timing, population and selection process,
retaining the weak distance branches. A [correlated synthetic recovery
test](../code/onefactor_sign_validation_v2.py) checks the sign-likelihood machinery under a
declared generator, separately from this observed-data response.
