# Independent scientific review of the audit synthesis

Reviewed 2026-09-21: `README.md`, `uncertainty-register.json`, the SNANA implementation report, and the saved extinction/population/response/covariance summaries. No phase-2 results or release assets were changed by this review.

**The synthesis supports its present conclusion: a physical-support problem is demonstrated, but its final observed-sample distance or cosmological bias is not measured.** It does not justify a replacement correction vector or a global numerical bias bound. The uncertainty register correctly uses `global_bias_bound=null`, distinguishes model alternatives from stochastic error modes, and prohibits adding everything independently in quadrature.

## Independent checks of the negative-extinction finding

I reopened all 25 original P21 Ia mock HEAD files independently and evaluated exact F99 with `extinction.fitzpatrick99` on their stored `SIM_AV` and `SIM_RV`. Results reproduce **71,946 written simulated objects**, 40,369 with RV<2, 11,965 with RV<1, and **4,791 with A(8000 Angstrom)<0**. Every stored AV is positive; the minimum is `4.01226e-6` mag. An independent root solve gives RV=`0.663385739303` at A8000=0. The threshold and direct-law signs agree for every object. Historical-law counts remain supported by the separately checked SNANA-source calculation; this reviewer did not substitute exact-law values for the historical count.

The broadband example is not caused by negative SALT flux. At `z=.1`, phase zero, `x1=c=0`, the DES i filter spans 6750–8700 Angstrom observed (6136.36–7909.09 rest). The dust-free model flux is strictly positive at all 20,001 independently sampled wavelengths with nonzero transmission, minimum `4.141e-18` in the model's observer flux-density units. With E=.1 and RV=.4, F99 attenuation is negative throughout this supported interval: approximately −.023424 to −.008236 mag.

Independent dense photon-weighted trapezoidal integration gives reddened/unreddened flux ratio **1.0200383283**, or `A_i=−.02154122710` mag. The saved sncosmo integration gives `−.02154119210` mag, a difference of **3.50e-8 mag**. Zeropoint and SALT magnitude-offset constants cancel in this ratio. Thus the broadband counterexample follows from the extrapolated extinction law even with a positive SED.

These checks establish behavior of a simulated component and a constructed light curve. They do **not** establish how many final BBC objects are affected, how fitted c and mB respond in those mocks, whether selection compensates or amplifies the effect, or the resulting mean distance/cosmology shift. Replacing the law, truncating RV, or flooring extinction also changes the inferred population; none is automatically the true correction.

## Scope and completeness review

The synthesis appropriately keeps separate parent-PDF moments, written-mock frequencies, monochromatic extinction, broadband constructed changes, local pre-BBC tangent responses, and released covariance modes. In particular, its 71,946-object denominator is labeled after generation/write selection, not as all generated transients, actual SNe, or the final cosmology sample.

The twelve-channel register covers the central SALT/dust concerns: intrinsic/dust decomposition, extinction implementation/support, latent population transport, MW foregrounds, calibration/training, spectral dimensionality, likelihood/covariance assumptions, covariance serialization, BBC selection, host/environment degeneracy, and unrestricted luminosity evolution. Required follow-ups explicitly include mismatched-generator recovery and coverage, selection normalization, held-out flux/color/host/redshift tests, and shared uncertainty propagation. It is a scoped uncertainty inventory, not a claim that every possible cosmological systematic has been enumerated.

I found one factual synthesis error and reported it: the component count was initially 46. The authoritative matrix summary contains **23 original plus 24 Dovekie components = 47**. The README now records 47 correctly. The broadband response numbers refer to the measurement-plus-model weighting and fixed published alpha/beta; their labels and linked report must retain those conditions. Near-complete flux-direction absorption includes the free amplitude, so it does not itself quantify a physical dust posterior or bound a cosmological shift.

No remaining scientific overclaim was identified in the reviewed conclusion. This review is not end-to-end closure; it specifically endorses leaving that requirement open.

Reviewed evidence identifiers:

| File | SHA-256 |
|---|---|
| `uncertainty-register.json` | `cc179a3a7a8e9267f0ee7e44b4d7de78a4584bb51c9574aeacd6e183d095e1fe` |
| `../../scripts/salt_dust_audit/f99_support.py` | `359ac32dd739787ae1bd33ed6f66c39a5c4ef614d52aaa78db8deebbd8dc3c62` |
| `../../runs/salt_dust_audit/f99_support/summary.json` | `028180786552b7b0a3e9db560814cdbed80702fa7d378371dbac4fe8cf7208dc` |
| `../../runs/salt_dust_audit/snana_mock_support/results.json` | `4c27c0bb8da84ae85157c5ee7e75d41596fde18b1b748553dda1439a92aa2e11` |

The README is a live synthesis; the numerical evidence above is pinned independently of subsequent editorial edits.

## Follow-up review: resolved grouping and bounded SNANA stress

The earlier snapshot above is historical. The revised synthesis explicitly names measurement-plus-model weighting for the response example, resolves Dovekie's apparent covariance-group discrepancy, and adds a bounded actual-SNANA refit. I reviewed the updated text, `matrix_grouping.py`, its results, the source filter rule, and the stress contract/results and per-object response table. This follow-up did not rerun every covariance inversion or the SNANA executable.

The grouping explanation is supported by two complementary checks: the saved `[CAL_SALT2] [+cal,=DEFAULT]` rule and uppercase substring matching in `create_covariance.py:1473` include CALSPEC; the numerical subtraction leaves nine resolved positive calibration/surface modes in each release. The Dovekie component sum after removing the separately listed duplicate closes to `6.4624033257e-7` relative statistically whitened Frobenius norm. The remaining Dovekie residual is numerical-scale compared with the full systematic covariance; its many tiny signed eigenmodes must not be promoted into astrophysical components. Original DES still has one resolved positive unassigned remainder, with whitened eigenvalue `0.5294990891` and relative norm `0.0034652942`. Its possible identification with the unreleased separate `SIGINT_MODEL` component remains a provenance lead. The supplied total already includes the remainder; adding it again would double count. This resolves the particular Dovekie grouping concern without certifying completeness of its physical uncertainty model.

The actual-SNANA branch is appropriately described as a **fixed selected-sample, hybrid noise-preserving intervention**, not a regenerated physical population. It modifies only negative host extinction, leaves measurement errors fixed, and fits the changed fluxes with the actual SNANA fitter. Its effect is small in this fixture: four of twelve deliberately selected low-RV objects have nonzero fitted standardized shifts, the largest `+0.0006548857674` mag; the other eight and twelve comparison controls are unchanged. I independently recomputed `delta_mB + 0.16087*delta_x1 - 3.1178*delta_c` from all 24 saved rows, agreeing with the reported column to `2.1e-16` mag. The stored baseline/no-op/reference checks are exact and the epoch audit reports no accepted-mask changes or unsupported injected accepted epochs. Those latter checks were inspected as recorded outputs, not independently rerun here.

This is useful balancing evidence: the physical-support failure does not imply a large direct fitted-distance effect in every affected mock. The small fixture response also does not bound changes from other laws, intrinsic/dust population refitting, detection/classification, sample migration, or regenerated BBC. Its ad hoc floor is not promoted to the true dust law. The earlier statement that fitted response was wholly unknown is superseded only for this specific intervention and fixture; the mean population/cosmology bias remains unknown.

The new [physical-alternatives report](physical-dust-alternatives.md) identifies observationally motivated discriminators and documents an additional prior limitation: public BayeSN still uses F99 and ordinarily truncates its RV population at 1.2. It is an alternative intrinsic/dust likelihood, not independent empirical exclusion of a tail its prior forbids.

Updated reviewed evidence identifiers:

| File | SHA-256 |
|---|---|
| `README.md` | `c19ec4b1a9f9615af834a4c9b8bd4efa2888272a0569d6ba602fd16a8236d1a3` |
| `uncertainty-register.json` | `f7a42eeff4d8b22ec3bc3754c4fdb01193ee8001b8dcd6c237501ef7aee01c1b` |
| `../../runs/salt_dust_audit/matrix_grouping/results.json` | `8a26f7b7b3c545766e17a2b2fb83526fcc1941dff746b3bc0d673f1e9731d850` |
| `../../runs/salt_dust_audit/snana_stress/results.json` | `8c9d8edf7b277c517d61cf978f494f3023debb8746bdf02490486fcb2027cd28` |
| `../../runs/salt_dust_audit/snana_stress/contract.json` | `207854bd389c0e9f3e822aa4a0f54ea94117b3d0d01f514887c9b2a33bb1bae5` |
| `../../runs/salt_dust_audit/snana_stress/fit_response.csv` | `f160e62675d34359a61bc45728b4d266adba3d159b8b70b92de5d3f31b40e1ca` |
