# CAMB release and numerical-support audit

This source-only audit found an upstream nonlinear-lensing spline defect in the **tagged CAMB 1.6.6 source**, as well as in 2.0.4. Its effect on our installed binary and likelihood is not measured here. A v2 upgrade would additionally change several physical approximations and cannot substitute for an isolated repair test. No CAMB calculation or package installation was performed for this audit.

The [source ledger](../results/inference/camb-v2-source-audit.json) pins CAMB 1.6.6 (`3ef0272d…`), 2.0.4 (`a6de8cc5…`), the September 25 fix (`56a95f78…`), the July paper, and technical documentation. Mutable issue/documentation responses are hash checked and refuse changed bytes on replay.

## A concrete legacy source defect

In 1.6.6, `fortran/cmbmain.f90` allocates `Transfer_Times(0:n+1)` at line 841, fills the physical times starting at element 1 at line 843, and passes the whole array to `spline_def` at line 1188. The latter declares explicit-shape `x(n)` in `fortran/subroutines.f90:12`. Fortran sequence association therefore sends element 0 as the first knot while the ordinates start at `scaling(1)`. Element 0 has not been initialized; the remaining knots are displaced by one. The same defect persists under the renamed spline routine in 2.0.4. [Legacy caller](https://github.com/cmbant/CAMB/blob/3ef0272d6f7ba1231128872e56e6d4c12af8267b/fortran/cmbmain.f90#L1188), [legacy callee](https://github.com/cmbant/CAMB/blob/3ef0272d6f7ba1231128872e56e6d4c12af8267b/fortran/subroutines.f90#L9).

The upstream fix supplies `transfer_times(1:n)`. Its author reports poisoned-memory failures on Linux and finite high-multipole lens-potential changes of about $5\times10^{-5}$ in their test. Those are upstream measurements, not estimates of our error. Our declared nonlinear-lensing calculation reaches this source branch, but the version string alone does not prove the installed wheel's compiled ancestry. Installed-wheel identification and matched unpatched/patched builds must precede any attribution of our observed accuracy dependence to this defect. [Exact fix](https://github.com/cmbant/CAMB/commit/56a95f78fdd72709c5af6b668141ec0c617777c3), [upstream report](https://github.com/cmbant/CAMB/issues/210).

## Why v2 is a separate comparison

| Item | Our 1.6.6 declaration/default | Explicit comparison requirement |
|---|---|---|
| Recombination | RECFAST Planck approximation | Select `recfast_approx_model='planck'`; the v2 default is a recalibrated CosmoRec fit. |
| Helium | PRIMAT 2021 BBN relation | Explicitly select `PRIMAT_Yp_DH_ErrorMC_2021.dat`, or pass each old model's exact helium fraction. |
| Neutrinos | One massive species, `mnu=.06`, `nnu=3.044` | Match the old `omnuh2`, species and degeneracies with `mnu=None, omnuh2_active=old_density`; v2 interprets the same mass through a thermal-density conversion. |
| Dark energy | Flat CPL, PPF, unit sound speed | Retain these assumptions and the same early-time domain. V2 replaces the hard high-k Gamma cutoff with capped relaxation; no public legacy PPF switch was found. |
| Nonlinear power | `mead2016` | Explicitly retain that model and its parameters, not the v2 default. |
| Lensing support | Final native request `lmax=9001`, `lens_margin=1250`, `max_l=10251`, `lens_potential_accuracy=4` | Inspect actual internal support. V2 renames the margin to `lens_output_margin` and also changes its Fortran meaning; copying the number is not an exact operator match. |

The RECFAST parameters are directly exposed in [v2 recombination code](https://github.com/cmbant/CAMB/blob/a6de8cc59124c8bbe3924f1c546f96cef73dcabb/camb/recombination.py). `use_rosenbrock=False` does not restore the old DVERK integrator: v2's alternative is RK45. Keeping the old recombination fit thus separates fit changes from a still-different numerical integrator. [Integrator description](https://github.com/cmbant/CAMB/blob/a6de8cc59124c8bbe3924f1c546f96cef73dcabb/docs/changelog/recfast_rosenbrock.md).

At our standard neutrino temperature, the old density proxy is approximately $6.448666\times10^{-4}$; v2's fixed-mass density is roughly $5.1\times10^{-5}$ larger fractionally. Supplying the old density avoids that parameterization change, but does not certify identical neutrino numerics. The inverted-hierarchy and very-light-neutrino fixes are not our one-species 0.06 eV configuration. [Neutrino mapping and limitations](https://github.com/cmbant/CAMB/blob/a6de8cc59124c8bbe3924f1c546f96cef73dcabb/docs/changelog/2026-07-29-neutrino-mass-density-mapping.md).

`AccuracyTarget=0` reproduces some older numerical choices, **not all older physics or source behavior**. PPF's new cap is compiled into its stress-energy routine independently of this setting. The source describes it as a numerical regularization of the same quasi-static limit, but its finite operator differs from the old hard cutoff and needs a distinct comparison label. It does not establish a physical dark-energy perturbation model. [PPF implementation](https://github.com/cmbant/CAMB/blob/a6de8cc59124c8bbe3924f1c546f96cef73dcabb/fortran/DarkEnergyPPF.f90#L134), [PPF experiments](https://github.com/cmbant/CAMB/blob/a6de8cc59124c8bbe3924f1c546f96cef73dcabb/docs/changelog/ppf_gamma_evolution_limit.md).

## Convergence axes still requiring attention

Increasing only `AccuracyBoost`, `lSampleBoost` and `lAccuracyBoost` does not independently test every relevant numerical choice. A bounded convergence study should separate differential-equation tolerance (`IntTolBoost`), recombination/background sampling, dense high-l interpolation (`min_l_logl_sampling=100000` reference), lensing angular sampling, convolution support, and source k support. Preserve physical parameters and output multipoles, then compare the actual likelihood components and derived distances as well as spectra. The reference must itself stabilize; generic pointwise spectrum tolerances do not imply a bound on a many-bin likelihood shift. [Official diagnostic workflow](https://camb.readthedocs.io/en/latest/check_accuracy.html).

Our explicit lens-potential setting is already 4; replacing it with the generic legacy recommendation 0 would reduce support. Conversely, v2's automatic rule at 9001 would choose about 15, so leaving it unspecified materially changes cost and k coverage. Nonlinear `AccuracyBoost * NonlinSourceBoost` can also enlarge transfer kmax and changes the redshift grid; this must be recorded when interpreting a boost comparison. [Multipole/support semantics](https://github.com/cmbant/CAMB/blob/a6de8cc59124c8bbe3924f1c546f96cef73dcabb/docs/changelog/lmax_settings.md), [nonlinear grid](https://github.com/cmbant/CAMB/blob/a6de8cc59124c8bbe3924f1c546f96cef73dcabb/docs/changelog/2026-07-29-nonlinear-lensing-redshift-grid.md).

High-l sparse interpolation is directly relevant to our request above 5000. The documented v2 repair increases the density of that sampling, and its reference disables the sparse logarithmic region. Several early tuning notes were later superseded, so their temporary floors must not be copied as final settings. [Final follow-up description](https://github.com/cmbant/CAMB/blob/a6de8cc59124c8bbe3924f1c546f96cef73dcabb/docs/changelog/HighL_accuracy_stability_followup.md). The separately documented Wigner-d fix concerns the full-sky method/new correlation implementation; it is not evidence of that defect in our legacy default method 1. [Scope of that fix](https://github.com/cmbant/CAMB/blob/a6de8cc59124c8bbe3924f1c546f96cef73dcabb/docs/changelog/2026-06-12-corrfuncfullsky-optimization-d4m4-fix.md).

## Performance implication

The July paper's four-thread flat, lmax=4000, matched-high-k benchmark is 4.2 versus 5.0 CPU seconds, or 1.1 versus 1.3 wall seconds: approximately 16% less time. Its much larger nonflat benefit is not applicable to our flat models. Our lmax, priors, component likelihoods, processor contention and single-thread execution differ, so no local speedup is established. Likewise, a faster recombination or HMCode subroutine is not a full-likelihood timing. The useful next step is the isolated legacy spline repair, followed only if needed by an explicitly controlled v2 comparison. [CAMB v2 paper, §XI/Table 1](https://arxiv.org/html/2607.14854v1).

Replay the source audit without CAMB:

```bash
python studies/unified_cosmology/code/external_probes/camb_v2_source_audit.py --acquire
```
