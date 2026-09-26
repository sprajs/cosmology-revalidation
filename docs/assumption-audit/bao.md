# BAO, joint-probe and CMB-reference assumption audit

Audit date: 21 September 2026. This examines the currently executed released-distance analysis, its original DESI reference chains, and primary DESI validation evidence. It does not change the SALT/SNM work or any `phase2` implementation. Reproducible additions are `scripts/assumption_audit/bao_*.py` and `runs/assumption_audit/bao/`.

**Main finding:** the local BAO vector and covariance implementation passes independent checks. A consequential, previously unquantified assumption is the CPL prior used to extrapolate the BAO-only result to today: extending just the lower bound on `wa` shifts mean q0 by about +0.27. The specific published-scale cross-bin BAO systematic correlations tested here have effects hundreds of times smaller. Dust enters BAO through target selection and clustering estimation, so a claim that BAO is intrinsically free of dust assumptions would be incorrect; the inspected evidence does not establish a dust-induced BAO peak bias.

## 1. What the code actually uses

`scripts/cosmology/core.py:216–234` reads the official DR2 `ALL_GCcomb` mean and total covariance. A fresh retrieval of both files from the [official Cobaya BAO release](https://github.com/CobayaSampler/bao_data/tree/master/desi_bao_dr2) is byte-identical to the local inputs; URLs, times and hashes are in `runs/assumption_audit/bao/primary/retrieval.json`.

The 13-element vector is one BGS DV/rd point, five galaxy/quasar DM/rd–DH/rd pairs, and a Ly-alpha **DH/rd–DM/rd** pair. The last pair's reversed ordering is intentional and correctly handled by observable labels. The LRG3+ELG1 measurement occurs only once. Comparing all seven individual release files reproduces every mean and every covariance block, apart from 2.2e-10 rounding in the more precisely printed Ly-alpha covariance. There is no covariance/precision inversion confusion, extra statistical diagonal, lost within-bin correlation, or tracer double-counting in this path.

Independent adaptive quadrature plus Cholesky evaluation agrees with the executed prediction to 8.6e-14 and chi-square to 1.2e-10 at 18 test cosmologies, including extremes of the extended wa prior. The covariance is positive definite (minimum eigenvalue .00579); six within-bin correlation coefficients range from −.347 to −.494. These checks and the row-level mapping are in `release_check.json`, produced by `bao_release_check.py`.

## 2. Covariance that is retained, assumed zero, and tested

Every cross-redshift entry in the released 13×13 covariance is zero. This is inherited from the product, not introduced by the local inverse. QSO redshifts overlap other tracers, so block diagonality is an approximation; the DESI release describes the QSO sample as shot-noise dominated. Its LRG3/ELG1 overlap is handled by the combined sample. The official baseline also assumes no cross-bin systematic covariance. [DESI DR2, Table III and Appendix B](https://arxiv.org/html/2503.14738v3).

I added only cross-bin covariance for the published theoretical error scales: 0.1% in alpha_iso and 0.2% in alpha_AP, preserving all diagonal and within-bin terms. The Jacobian is `delta ln DM = delta ln alpha_iso − delta ln alpha_AP/3`, `delta ln DH = delta ln alpha_iso + 2 delta ln alpha_AP/3`. The six galaxy/QSO blocks are correlated; Ly-alpha retains its separate error budget. This follows the class of tests in [DESI validation §VII.1](https://arxiv.org/html/2503.14742v3#S7.SS1), with correlation coefficients .5 and 1. Existing chains are importance-reweighted, so no SN correction is recalculated.

| Existing fit | Baseline mean q0 | Shift with rho=.5 | Shift with rho=1 |
|---|---:|---:|---:|
| BAO only | +.03517 | +.00044 | +.00089 |
| BAO + original Pantheon+ | −.43608 | +.00019 | +.00038 |
| BAO + fixed median age template | +.06325 | +.00007 | +.00014 |
| BAO + shared age-slope uncertainty | +.04219 | +.00007 | +.00015 |

Weight effective-sample fractions exceed .996 for rho=1; paired block Monte Carlo errors on joint-fit shifts are about .00003–.00004. These results constrain this particular small theory covariance, not an arbitrary hidden systematic. A deliberately incorrect diagonal-only control changes the BAO-only mean q0 by +.060; retaining the existing within-bin anticorrelations is consequential. That control is not a proposed alternative likelihood.

**Unresolved rather than disproven:** `scripts/cosmology/run.py:61–66` adds SN and BAO chi-squares. This sets their measurement-error cross-covariance to zero and has no shared dust-map or calibration nuisance. Shared sky structure, Galactic foreground errors and imaging calibration can in principle correlate their errors; shared cosmological parameters alone are not a reason to insert measurement covariance. The compressed files contain no derivative of the BAO estimator with respect to the SN dust/calibration nuisance, so this audit cannot estimate its size or sign. The correct next test would perturb the same foreground/calibration field in both pipelines and remeasure both observables. Adding an arbitrary off-diagonal matrix would not be evidence-based.

## 3. A substantial prior sensitivity, newly measured

`scripts/cosmology/run.py:44–56,68–71` uses uniform wa∈[−3,2] and w0+wa<0. BAO's lowest effective redshift is .295; its q0 is an extrapolation through this finite CPL family. I changed only the lower wa bound, preserving the 13 observations, total covariance, flat geometry, remaining priors and free H0rd. These are sensitivity experiments, not claims that a broader prior is physically preferable.

| BAO-only wa support | Mean q0 ± posterior SD | P(q0<0) | Posterior mass below wa=−3 |
|---|---:|---:|---:|
| [−3,2], existing run | +.035 ± .284 | .388 | 0 by definition |
| [−5,2], new run | +.248 ± .375 | .249 | 38.0% |
| [−10,2], new ensemble 1 | +.296 ± .409 | .232 | 42.1% |
| [−10,2], new ensemble 2 | +.307 ± .400 | .222 | 42.6% |

Both [−10,2] ensembles retain over 101 autocorrelation times for every parameter. Their q0 means agree within 1.2 combined block Monte Carlo errors; probabilities agree within 1.3. Neither accumulates near wa=−10. Their best fit is unchanged (chi-square 5.6186), demonstrating a change in integrated posterior volume rather than a numerical failure or improved fit. Other priors, including w0≤1 and w0+wa<0, remain. Agreement with published BAO-only centres therefore verifies the conditional calculation; it does not establish a prior-robust present-acceleration probability.

**The large effect does not survive adding these SN distances.** Matched 30,000-step original-joint ensembles give q0=−.4360±.0801 for wa≥−3 and −.4350±.0806 for wa≥−10. Their .00105 mean difference is 1.1 combined block Monte Carlo errors; both retain over 143 autocorrelation times per parameter, have no wa<−3 draws, and no non-accelerating draws. The fixed-age joint extension gives +.0647±.0745 with P(q0<0)=.192, versus the existing +.0633±.0742 and .198; only .066% of its draws are below −3. These are checks of the existing frozen age template, not a new age correction. One shorter original-joint extended-prior ensemble gave −.4310±.0795 despite no excluded-support draws. That run is retained as a sampling discrepancy; it motivated the longer paired ensembles and is not attributed to a physical prior shift. Existing original-joint repeats and identical source/input hashes were also checked.

A smaller coordinate-prior effect is also measurable. Replacing a uniform H0rd prior by uniform log(H0rd), within the same bounds, moves the BAO-only mean q0 from .0352 to .0447; a uniform inverse-H0rd prior gives .0540. For original joint SN+BAO the respective shifts are only .0004 and .0008. These are explicit prior changes, not missing statistical error terms.

## 4. How dust, host age and fiducial cosmology enter BAO

DESI targeting uses SFD Galactic extinction corrections. The sample-validation study compares SFD with stellar-spectrum and CIB-cleaned maps. Some residual dust-map correlations trace extragalactic cosmic-infrared-background structure; blindly removing them can erase real clustering. It documents very large extinction-related changes to ELG broadband clustering while the BAO peak remains stable. Thus the inference distinction is between foreground-sensitive sample density and the fitted acoustic scale. [DESI sample definitions §6 and §10.1](https://arxiv.org/html/2411.12020v1).

DR2 additionally tests removing imaging weights, sky/imaging-region splits, reconstruction bias, alternative fiducial cosmologies and Gaussian compression. These tests support the released peak measurement within the examined alternatives; they are not proof that all possible foreground patterns are harmless. Published validation is inherited here, not independently reproduced from catalogues. [DESI DR2 validation §§VI.5–VII](https://arxiv.org/html/2503.14742v3).

No stellar age or luminosity-versus-age slope is inserted into this BAO likelihood. Host populations can still affect the galaxy selection, bias and reconstruction used upstream. A supernova age correction must not automatically be applied to the acoustic ruler. Conversely, “no age parameter in the 13-row likelihood” does not mean upstream selection is population independent.

`core.py:73–76,225–230` assumes flat FLRW CPL geometry and a common, redshift-independent ruler. A free H0rd removes the need for a numerical CMB sound-horizon calibration, but does not remove those assumptions or the extraction model behind the compressed data. A coherent 1.5% scaling of every BAO distance and its covariance is exactly absorbed by H0rd→H0rd/1.015 (chi-square difference 5.2e-14 in the check). Such a pure scale shift cannot explain a change in this uncalibrated distance-shape inference; a redshift/tracer-dependent bias can. The finite H0rd prior can qualify exact posterior invariance near its boundaries, which are unoccupied in the inspected runs.

For joint inference, `core.py:102–105` maps metric distance to luminosity distance with the standard redshift prefactor. There is no independent opacity/distance-duality nuisance. Residual extinction or gray attenuation varying with redshift can therefore enter as SN distance-shape error. Constant extinction would instead be absorbed into the marginalized SN magnitude offset. The existence of this degeneracy does not determine its physical amplitude.

## 5. CMB-reference assumptions and overlap

The actual `chain.updated.yaml` files specify CAMB/PPF, flatness, one massive neutrino with total mass .06 eV, Neff=3.044, BBN helium assumptions, Planck low-ell TT/EE, NPIPE CamSpec high-ell TTTEEE and ACT+Planck lensing v1.2. CamSpec has foreground amplitudes/slopes and calibration nuisance parameters; those chains cannot be described as foreground-uncertainty-free observations. The audit JSON records the exact configurations and hashes, rather than inferring them from directory labels.

Those reference chains also use wa∈[−3,2]. The BAO+CMB-only posterior puts 3.05% of its weight below −2.8 and .63% below −2.95; neither CMB+SN reference has retained weight below −2.8. This identifies where a wider-prior CMB run would be relevant, but cannot quantify the unseen posterior outside the original support. Reweighting existing chains cannot create that missing support, so the BAO-only prior result must not be numerically transferred to CMB.

The lensing component uses `variant: actplanck_baseline`, `lens_only: false`, likelihood corrections enabled and a Hartlap correction (`chain.updated.yaml:132–149` in the BAO+CMB-only directory). The inspected ACT code loads `covmat_actplanck.txt` for this branch (`sources/repos/ACTCollaboration__act_dr6_lenslike/act_dr6_lenslike/act_dr6_lenslike.py:354–430`); it does not add independent ACT and Planck lensing likelihoods. The primary ACT analysis describes a joint covariance from common sky simulations. [ACT DR6 lensing cosmological analysis](https://ntrs.nasa.gov/api/citations/20240002550/downloads/Madhavacheril_ACT%20DR6%20Lensing%20Cos%20Parameters%20final.pdf). This does not establish that every cross-covariance with SN or BAO is included.

`reference_chains.py:24–44` correctly uses sample weights and reports alternative burn fractions. Its radiation-free low-redshift q formula is an acknowledged approximation. Restoring today's radiation density inferred from `1−omegam−omegal` shifts mean q0 by only +.00010 to +.00014, with massive neutrinos still treated as nonrelativistic matter. This cannot explain the reported sign sensitivity. The output is derived from original posteriors; it is not a corrected CMB likelihood, nor an independent dust reanalysis.

## Reproduction and limits

Run from the workspace root with `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/assumption_audit/bao_release_check.py`. Prior experiments use `bao_prior.py --wa-min -10 --seed 210921` (repeat seed 210922), or `--wa-min -5 --seed 210923`. The long joint comparisons use `--data joint --wa-min -10 --seed 210926 --steps 30000 --burn 5000` and `--data joint --wa-min -3 --seed 210927 --steps 30000 --burn 5000`; the fixed-template check uses `--data joint-c14fixed --wa-min -10 --seed 210925 --steps 10000 --burn 2500`. Full run settings and source/input/output hashes accompany the outputs. Frozen source snapshots preserve the initially executed BAO-only script as well as the current extension.

**Concurrent-edit provenance:** another task changed only `qbins` handling in `core.py` during the two long joint runs. Their original manifests hash source at completion, so the recorded whole-file version cannot alone establish the imported version. Both exact snapshots are preserved under `runs/assumption_audit/bao/provenance/<sha256>.py`; the old hash starts `18cc923c`, the newer hash `83f852a3`. `bao_verify.py` confirms byte hashes and whole-module AST identity after removing only the `model == 'qbins'` branches and their added domain helper; the Pantheon, BAO, magnitude and interpolation implementations are exactly identical. These CPL results are therefore unaffected by that edit. Historical manifests remain unchanged, and future `bao_prior.py` runs capture source bytes at import. `verification.json` records the timing limitation, exact snapshots and report hash. Core line references above use the newer snapshot.

No confirmed BAO implementation bug was found. The material new limitation is BAO-only prior sensitivity; zero cross-probe covariance, restricted geometry/opacity, upstream selection and fixed likelihood products are additional explicit model boundaries. A complete foreground-field perturbation of both SN and DESI estimators requires raw catalogues and pipeline assets beyond the compressed products tested here.
