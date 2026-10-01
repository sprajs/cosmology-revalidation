# Released compact SH0ES ladder

This packet executes a bounded numerical experiment on the **unchanged 3492-row,
47-column compact release**. It also audits the primary source identities. The
paper reproduction packet remains blocked because the complete calibration,
selection and parameter-measure lineage is unresolved. Its [experiment.json](experiment.json)
has no generic runner route; the reviewed experiment-specific
[controller](controller.py) below runs the narrower deterministic full-design
calculation. The existing LCDM packet and all its original pins/receipts are unchanged.

[Riess et al. 2112.04510v3](https://arxiv.org/abs/2112.04510v3), printed pages
10–12, defines the empirical Cepheid/SN ladder and Equation (6). The release at
[commit c447f0f](https://github.com/PantheonPlusSH0ES/DataRelease/tree/c447f0fea703fcd0fff57de5000947b5ca81286b/SH0ES_Data)
provides y, L and full C. The source equation and its Python orientation establish
X=L.T and q=(y−Xβ)ᵀC⁻¹(y−Xβ). This experiment minimizes that quadratic in all
47 original coordinates, reports its relative score −q/2 and the conditional
fixed-design Gaussian estimator variance for original contrast e46. Production
factorization, fitting and variance are entirely in Irreducible's installed
native C++ library; [consumer.cpp](consumer.cpp) supplies bounded transport and
status handling. It retains one covariance/QR preparation across evaluations.
No determinant, normalized parameter density, evidence or posterior is reported.

## Primary source dictionary and unresolved identities

[The lineage manifest](lineage.json) pins exact URLs, byte sizes, SHA-256 hashes,
original source-axis indices and equation/script locators. The complete 37-host
order is established by an executable join, not by sorting names or matching
approximate fitted moduli. For **every** initial SN-host Cepheid row (0–2149),
its original column41 log-period coordinate and column43 metallicity coordinate
are matched to the pinned primary `table2.tex`. Matching uses the interval of each
printed decimal value (half its last decimal unit), the FITS binary32 half-ulp,
and a separately declared 1e-14 allowance in log-period for `log10` evaluation.
The intersection of possible source hosts across all rows in each host-column
support must be exactly one host. No unsuccessful row is dropped and no tolerance
is adapted. The result records original row indices and primary table line
locators. This verifies group identities, not unique individual event identities,
complete selection reconstruction or observational qualification.

Source table names `N105A` and `N976A` are preserved verbatim; an alias to paper
names N0105/N0976 is context and is not silently imposed. The compact N1365
column contains **46** initial Cepheid rows, all matching its source photometry;
paper Table 3 prints **45**. That release/paper selection discrepancy remains open.
The join corroborates the period/metallicity coordinate assignment; the release
README requires adding the fitted slope correction back to the initial −3.285.
Column46 is independently identified by the paper's vector before Equation (6)
and `read_chains_example.py`'s selected-index/conversion code as
5log10(H0/[1km/s/Mpc]). Its additive reference-coordinate shift would leave
contrast variance unchanged. This adapter reports the source logarithmic
coordinate, without introducing another H0 projection implementation.

The paper schematic suggests the identities of the remaining anchor/luminosity,
SN and ground zero-point coordinates. These correspondences are retained as
context, with **null physical identities** until the complete release mapping is
verified. Column44 has no established physical identity and no schematic label.
No coordinate is renamed, eliminated or regularized. Original external-constraint
rows3207–3214, their exact y/Cii values and nonzero design support are recorded by
the audit separately from physical labels. The full C retains dependence between
standardized calibrator/Hubble-flow SNe and all other encoded covariance; unknown
external probe overlap is not independence.

## Different target measures

The source MCMC example is not the same parameter-measure contract as this full47
relative profile. `run_mcmc.py` sets `PRIOR_WIDTH_RATIO=10`; `MCMC_utils.py` sets
each box half-width to ten times the corresponding `lstsq_results.txt` width,
initializes walkers uniformly inside that box, and returns a constant log-prior
inside it. Ordered row 44 is exactly `(0,0)`, so this code fixes β44=0 with
zero-width support. That is **not a proper density in a 47-dimensional Lebesgue
parameter measure**. A normalized posterior requires an explicit lower-dimensional
measure and treatment of that coordinate; this packet does not supply one. The
helper fit table supplies initialization/prior limits, not a final numerical
acceptance target. Numerical acceptance of the full47 profile does not qualify the
released chain or reproduce observational calibration uncertainty.

Paper Section 5 reports its fitted logarithmic coordinate as 9.318, its fitted
H0 as 73.04±1.01km/s/Mpc, and an increased uncertainty±1.04 after analysis-variant
systematics. The controller checks only agreement of β46 with the **rounded
9.318 coordinate**, within half the final printed digit. That coarse source
consistency check is separate from numerical QR/SVD comparisons and does not
assign the systematic uncertainty to the supplied fixed-C estimator variance.
The source's277 Hubble-flow SNe are selected at 0.0233<z<0.15 with matching quality
cuts and late-type hosts; standardized/selected products and constraints are
already encoded. Full raw selection and source-effect reconstruction is blocked.

The optional [reference](reference.py) is reference-only SciPy/LAPACK SVD and
pivoted QR, with full-C Cholesky whitening. Both reference algorithms share
covariance-whitening ancestry; native portable wide covariance factorization and
scaled pivoted QR have separate implementation ancestry. Agreement tests those
named inputs/methods, not universal independence or physical correctness. Frozen
allocations remain coefficient 1e-8 absolute+1e-9 relative, q 1e-7 absolute+1e-10
relative, variance 2e-12 absolute+2e-12 relative. No row removal, jitter or budget
relaxation is allowed. External refinement/high-precision evidence for the
permanent engine kernels remains in Irreducible; this reference is a separate
cross-algorithm check, not a high-precision refinement claim.

A declared synthetic sensitivity control adds/subtracts0.01 in turn to y of each
original constraint row 3207–3214, at unchanged X,C. The native profile reuses
retained factors and compares each perturbed β46/q with full47 SVD. This answers
how this conditional numerical experiment responds to those exact input-row
perturbations. The perturbations are synthetic controls, not estimated calibration
systematics, physical prior changes or an observational uncertainty budget.

## Execution and provenance

Acquire the exact named sources from the URLs in `lineage.json` under an ignored
local directory, without overwriting source originals. The covariance URL uses
GitHub media because the raw repository member is an LFS pointer. Papers, source
assets and FITS inputs are not redistributed: no upstream redistribution license
has been verified. The adapter/reference sources are original repository code
under BSD-3-Clause. The historical archives are read-only acquisition evidence.

The reviewed local SDK is clean source commit
`c9b7febc01ecc85d289056faba8722a8f24a6ad8`, build ID
`a08f62cb32a76097e68beff6a540ccf4da903b5b898ae75921ee53818c913903`.
Exact manifest/archive/CLI identities live in the controller. It verifies every
manifest source hash against **immutable pinned Git blobs**, reconstructs the
build ID, checks the complete installed header inventory, compiler/standard-library
bytes, archive and CLI. The primary checkout may have moved forward; none of its
mutable source files are compiled or used as this historical build's evidence.
This local SDK is not a distributed release. A new build needs a deliberately
reviewed new build identity; source commit alone cannot reproduce archive bytes.

From the Reproducible root, with the exact sources and preserved SDK available:

```sh
uv run python experiments/released-ladder/controller.py \
  --sources .work/released-ladder-lineage-20261001/sources \
  --engine-source ../irreducible \
  --sdk ../irreducible/evidence/project-review/science/prerequisites-20261001/final-primary/sdk \
  --name fresh-full47-attempt \
  --reference-python /absolute/path/to/python-with-numpy-and-scipy
```

The reference environment is optional and separate from production dependencies.
One compiler/compute job is used; reference threads are pinned to one. Native
preparation has a 2 GiB payload ceiling and 1e-10 sensitivity screen; child address
space is bounded to 3GiB, output to 1MiB (compiler output 8MiB), compile time 120s
and each scientific child 900s. Payload bounds are distinct from process memory.
The controller admits only exact unscaled binary32 primary FITS images and
unchanged canonical binary64 hashes; checks full-C symmetry and final contrast
support; and verifies source/SDK/packet identity again after execution, including failed
children once admission identities are available. It keeps
source adapter snapshots, actual compiler command/executable hash, SDK/input
identities, all host joins, native/reference outputs and terminal gates under a
fresh ignored `results/released-ladder/<name>/` directory. Existing attempts are
never overwritten; terminal files are made read-only, including failures.
Local receipts still require a durable external archive before scientific publication.

## History

2026-10-01: primary-source audit closed the complete 37-host column order using all
2150 initial Cepheid rows, corroborated period/metallicity coordinates, and pinned
the source-defined column46 contrast. It exposed the zero-width released MCMC
coordinate 44 and the N1365 selection-count discrepancy. Complete anchor/nuisance
identities, raw selection/calibration lineage, a normalized posterior measure and
analysis-variant systematics remain blocked.


2026-10-01: clean committed adapter `9ff0af96ddbf97444c604fea27a5b909ed761353`
completed the full source join and native/SVD/pivoted-QR experiment in 69.93s.
All 130 named comparisons passed at the unchanged allocations; maximum budget
fraction 1.52e-5. Native β46=9.317887632654045, q=3552.759330295523 and
conditional Var(β46)=0.0008964788320013900; SVD variance differs by 2.17e-18.
The source's rounded 9.318 coordinate agrees. All eight synthetic constraint-row
perturbations agree with the full47 SVD law; row 3210 leaves β46/q unchanged in
this named control. These results establish a bounded fixed-design numerical
comparison, not a parameter posterior or the paper's analysis-variant systematics.
The complete ignored receipt is
`results/released-ladder/committed-full47-lineage-20261001/record.json`, SHA-256
`3cc6f9529dfd4b78ba2bdff474b6b3b7ce0f68aa44582f022e25eadfabbc1fd8`.
The earlier diagnostic attempt is preserved as failed because the reference
Python executable was resolved outside its supplied environment; native outputs
remain in that receipt. The executable path handling was corrected before the
clean committed run. Subsequent integrity hardening records post-run SDK/source/
packet identities even after failed children; the original receipts are unchanged.
