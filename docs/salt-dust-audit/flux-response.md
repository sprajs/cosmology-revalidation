# Dust can hide in a good SALT light-curve fit

This independent calculation uses the original public DES SALT3 surface, calibrated DES passbands, actual accepted observing times and quoted flux errors for 64 objects. It establishes a **conditional sensitivity and collinearity**, not a measured dust bias, new physical dust population, or corrected cosmology. It does not fit the dust parameters to the observed fluxes. Published SALT coordinates specify each reference SED, and deliberately altered model fluxes test which changes a SALT fit could absorb.

## What was calculated

Objects are equally spaced ranks in redshift among original DES Hubble members with published accepted light-curve points: zHEL=0.05999–1.12212; 2,581 accepted epochs. This deterministic diagnostic sample is not a probability sample. All 64 completed. Membership is in [selected_objects.csv](../../runs/salt_dust_audit/flux_response/selected_objects.csv), and every epoch's time, band, measured flux/error, and reference model flux is retained in [matrices.npz](../../runs/salt_dust_audit/flux_response/matrices.npz).

Reference parameters are published x0,x1,c,t0,zHEL. MW reddening is the full-precision HEAD value, already scaled: multiplying it by .86 again would be incorrect for these inputs. SNANA's SALT `MAG_OFFSET=0.27` and the released KCOR AB-primary magnitudes are applied in the flux normalization. This avoids mixing the amplitude conventions of SNANA and sncosmo. The script uses sncosmo 2.12.1 and extinction 0.4.9; it does not fetch a differently trained model from the online registry. F99 foreground dust acts at observer wavelength, and additional host dust acts at rest wavelength. The nominal fitted colour is not reinterpreted as the intrinsic colour or subtracted as physical dust.

Let theta=(dmB,dx1,dc,dt0), where increasing dmB dims the amplitude. At reference flux f, compute J=df/dtheta, and G=df/deta for nine nuisance modes: added host E(B−V) at each of R_V=(1.5,2,3.1,4), added foreground E(B−V) at R_V=3.1, and individual g,r,i,z dimmings in magnitudes. The four host columns are alternative laws, not four independently established dust populations. The R_V=1.5 column extrapolates the empirical F99 law below its documented recommended range and is a stress case.

For fixed flux covariance C, whiten with C=LLᵀ, Jw=L⁻¹J and Gw=L⁻¹G. SVD gives

```
R = pinv(Jw) Gw                       # fitted-parameter response to added flux effects
G_res = Gw - Jw R                     # part distinguishable from fitted SALT parameters
D = (1, alpha, -beta, 0) R            # standardized brightness response
V_theta = (Jw^T Jw)^(-1)              # conditional local parameter covariance
K_res = G_res^T G_res                 # information remaining about the nuisance modes
```

Here alpha=.16087 and beta=3.11780 are held fixed. The response excludes BBC changes, host assignment changes, nuisance re-estimation, retraining, clipping, selection and contamination. **D is not the derivative of the complete DES distance pipeline.** It cannot be pasted onto published corrected distances as an empirically measured correction.

Two weight choices are retained: quoted measurement variance alone, or measurement plus the sncosmo SALT3 model covariance. The latter includes the model's same-passband colour-dispersion correlations. Neither represents full cross-object training/calibration uncertainty or is asserted identical to SNANA's floors, interpolation, MW covariance and rejection algorithm. Covariance and epoch selection remain frozen under perturbations. The baseline treats accepted repeated rows as independent measurement noise; the other thread's duplicate audit remains relevant.

## Numerical results

For measurement-plus-model weights, the median fractions of **squared whitened perturbation norm** absorbed into the four fitted parameters exceed 99.98% for every host law and 99.99% for MW extinction. This percentage is neither a fraction of physical dust removed nor a posterior probability. It quantifies geometric alignment of flux responses.

| Added extinction hypothesis | Median D (mag per mag E(B−V)) | Range across these objects | Median remaining SNR for delta E=.01 |
|---|---:|---:|---:|
| Host, R_V=1.5, extrapolation | −0.4409 | −0.7033 to +0.4437 | 0.0108 |
| Host, R_V=2.0 | +0.0412 | −0.2025 to +1.0908 | 0.0082 |
| Host, R_V=3.1 | +1.1165 | +0.8383 to +2.3539 | 0.0089 |
| Host, R_V=4.0 | +1.9886 | +1.6022 to +3.3263 | 0.0093 |
| MW, R_V=3.1 | +0.1688 | −0.1699 to +0.6936 | 0.0055 |

The response varies with SED, redshift, phase coverage, available filters and weights. In particular c does not respond exactly as E(B−V), and beta is not automatically R_V+1. The smallest residual component is not necessarily the smallest standardized-brightness shift. At E=0, changing host R_V alone has zero flux response: an R_V estimate requires nonzero reddening and population or other information.

The numerical signs mean additional extinction in the *observed* flux. An error in an adopted foreground model changes the inferred fit in the opposite direction to the corresponding observed-flux perturbation. A passband's absent observations produce an exactly zero response, not an empirical calibration constraint. Summing the four equal passband dimmings gives exactly dmB=1 with no residual; it is degenerate with the amplitude/intercept.

An independent nonlinear check refits the injected model fluxes at every eighth object, holding the same covariance. All 24 fits converged. For R_V=3.1 and delta E=.01 the largest linear/nonlinear standardized-brightness difference is 0.0000878 mag. At delta E=.05 it grows to 0.00218 mag, so unit derivatives must not be extrapolated arbitrarily. In those eight delta E=.05 cases, the nonlinear brightness shifts span +.0416 to +.0803 mag, with maximum additional residual chi-square .110. These finite injections are constructed alternatives, not measured reddenings or errors.

Halving finite-difference steps changes the total Jacobian norm by at most 1.13e−6 fraction. Exact gray-response, projection orthogonality and analytic amplitude identities are asserted for every object/weight choice. The 64 reference median observed/model flux ratios have median .9944 and range .8716–1.0799; this is a normalization sanity check only, not proof of an exact independent baseline reproduction.

## Using the matrices without inventing uncertainty

[response.csv](../../runs/salt_dust_audit/flux_response/response.csv) gives signed R and D, object IDs, weights and residual norms. The NPZ additionally records C, J, G, V_theta, full and residual nuisance Gram matrices, and singular values. Keys use `CID__quantity`; both weight choices are explicit. Parameter and nuisance orders and units are in [summary.json](../../runs/salt_dust_audit/flux_response/summary.json).

For a justified nuisance covariance S, propagation is `C_theta=R S Rᵀ` and `C_mu=D S Dᵀ`. Across SNe, a shared calibration or MW mode must use a common nuisance column: it creates cross-object covariance, not just diagonal errors. Independent sightline perturbations require separate columns; a population change can correlate them. No S is inferred or fabricated here. Combining this with an existing released systematic covariance requires checking whether the same uncertainty is already included. The covariance represents uncertainty about a specified response, not a guarantee that its mean bias is zero.

Reproduce from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 phase2/env-official/bin/python scripts/salt_dust_audit/flux_response.py
```

The [manifest](../../runs/salt_dust_audit/flux_response/manifest.json) hashes source inputs, the script, package implementation and outputs. The [nonlinear checks](../../runs/salt_dust_audit/flux_response/nonlinear_checks.csv) preserve all fitted perturbations. The implementation follows the primary [sncosmo SALT3 interface](https://sncosmo.readthedocs.io/en/stable/api/sncosmo.SALT3Source.html), [observer/rest-frame propagation contract](https://sncosmo.readthedocs.io/en/stable/api/sncosmo.Model.html) and [F99 dust interface](https://sncosmo.readthedocs.io/en/stable/api/sncosmo.F99Dust.html), with the locally pinned implementation authoritative for this run.
