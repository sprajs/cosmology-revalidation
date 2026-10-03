# Conditional SN symmetric working profile

Can Irreducible's Gaussian/profile kernel and an independent full-mode spectral calculation agree for the complete released SN+SH0ES selection on a specified derived covariance? The private v5 attempt passed this conditional numerical comparison. It does **not** reproduce the released likelihood, infer H0, or supply a normalized posterior or joint CMB/BAO/SN result.

This is a manual development reference with eight source files. The six code files and [request.json](request.json) are literal copies of the reviewed private v5 sources and proposal. The proposal deliberately has `execution_admitted=false` and null reference Python/runtime pins. Public paths/import identities are new; the successful private attempt does not qualify a cold replay of this copy.

## Inputs and derived target

The original PantheonPlusSH0ES/DataRelease revision is `c447f0fea703fcd0fff57de5000947b5ca81286b`. Select `(zHD > 0.01) OR calibrator`: retain 1657 of 1701 original occurrences, including all 77 calibrators; exclude 44. Preserve increasing original indices, repeated events, redshift frames, calibration fields and the complete STAT+SYS principal covariance. Occurrence IDs include the table hash and original index; names never collapse rows. SN+SH0ES calibration/host dependence remains.

The selected binary64 matrix B fails exact symmetry: 361 unequal unordered selected pairs, 389 in the original full matrix. Preserve the first refusal and complete diagnostic/pair scan. Printed-decimal rounding overlap is diagnostic, not permission to repair B. The released Cosmosis caller passes the full covariance unchanged; its inherited GaussianLikelihood implementation is unavailable here. No raw-B numerical score is evaluated.

The distinct target is `released-sn-symmetric-working-profile/v1`; its derivation is `selected-binary64-rational-pair-projection-RN-even/v1`. For only the 361 unequal stored binary64 pairs, form the exact rational average S0_ij=(B_ij+B_ji)/2, cast **once** with verified nearest-even binary64 rounding, and write identical chosen bits to both cells. Equal pairs and diagonals retain their original bits, including unchanged signed zeros. S_hat is this stored matrix, without a two-operation float half-sum, triangle choice, jitter, clipping or row removal.

Assembly retains separate row-major big-endian binary64 B/S_hat artifacts, each 21,965,192 bytes. The ledger records both original hex values, exact averages, R=S_hat−S0, A=(B−B^T)/2, changed cells and exact row norms. Physical/source uncertainty remains unknown. If the unavailable producer used r^T B^-1 r, its symmetric precision would be Sym(B^-1), with effective covariance S0−A S0^-1 A when applicable; S_hat is a different target.

## Means, profile and checks

Reuse the retained CLASS cases in order: anchor, precision, ns-minus, ns-plus. Pin their actual configuration/input, executable, argv, background, thermodynamics and receipt. The original theory/BAO/SN mean helpers execute from consumed pinned buffers in dependency order, with explicit code provenance. They are not replaced by later disk/cache imports or new Python physics.

Calibrators use CEPH_DIST. Other rows use `mu = 5 log10((1+zHD)(1+zHEL) DA_CLASS(zHD)) + 25`, with DA in Mpc and the source frame convention. The four retained backgrounds and assembled mean vectors are bitwise identical: this is one shared corpus, not four independent changes or an SN sensitivity test of n_s. Retained grid interpolation has no invented certified error bound.

For r=m_b_corr−mu, profile one free, unbounded common M:

`g = 1^T S_hat^-1 1`, `M = (1^T S_hat^-1 r)/g`,

`q = (r−M1)^T S_hat^-1 (r−M1)`, `relative profile score = −q/2`.

This differs from the supplied fixed-M route. There is no M prior, determinant normalization, evidence, posterior, independent SN/H0 information or combined calibration uncertainty.

The native consumer uses the existing historical f3b19 Gaussian/ProfileOperator API, build 6d495 and F02/binary64-legacy arithmetic. Prepare S_hat once, cache the ones response once, then evaluate four ordered residuals. Require finite/ok statuses and cached-response, coefficient and adjusted-solve forward-sensitivity diagnostics <=1e-8. Retain actual adjusted residuals; the internal wider coefficient and reporting casts forbid reconstructing them from printed M. Refused preparation retains statuses and null unearned cache fields. The API exposes cached response/Gram but no adjusted solution vector: native adjusted stationarity is null. The separate cached-response original-matrix residual is diagnostic only and never feeds the reference. Compare emitted 17-digit native decimals as Decimal values.

The reference uses actual NumPy 2.2.6 `eigh(S_hat, UPLO='L')`/LAPACK_syevd after full symmetry/finite admission. Keep all 1657 finite positive modes and eigenpair artifacts; no truncation, clipping, fallback or native output. Reconstruction/Q^TQ use binary64 BLAS; projection, positive weighted sums, solves and original-matrix residual/norm reductions use observed 64-significand-bit longdouble/RN. These are engineering checks, not interval certificates.

With b=Q^T1 and a=Q^Tr, use M=sum(ab/lambda)/sum(b²/lambda) and q=sum((a−Mb)²/lambda). Define rho=||S_hat−Q diag(lambda)Q^T||inf/||S_hat||inf, delta=||Q^TQ−I||inf, kappa=lambda_max/lambda_min and eta as the maximum response/adjusted original-matrix residual normalized by `||S_hat||inf ||x||inf + ||rhs||inf`. Require finite positive modes/Gram, delta<1 and

`kappa*(rho + delta + eta + 1657*2^-52) <= 1e-8`.

Each case also requires `abs(sum(x))/(sqrt(1657)*||x||2) <= 1e-8`; a zero denominator permits only the exact zero vector/numerator. Each same-S_hat native/reference relative-profile difference must be <=1e-6. Reference scalars retain 36-digit strings; solution reportcasts never feed gates. R/A and optional triangle sensitivity estimates have unvalidated eigenvalue bounds and are separate engineering diagnostics, not physical source-error allowances. Triangle variants and raw-B inversion are not run.

Controls retain full 1657 dimensions: diagonal-plus-rank-one/exact Sherman–Morrison truth, every adjusted lane/constant shifts; conditioning, positive-diagonal non-SPD projection, nonsymmetry, nonfinite, order and missing-row refusals. Original B is preparation-only. Exact skew adversaries cover near-singularity and sign/zero limits; RN-even ties, equal cells/diagonals and consumed-stream failure evidence are checked. Stop first failure and retain earned prefixes/raw outputs without retries or gate changes.

## Dated evidence and limits

2026-10-03: the exact-B structural route refused (a253e53f); the f094ef45 diagnostic retained all 389 asymmetric pairs without repair. Metadata v2 refused for missing threadpoolctl; v3 refused because its raw complete runtime exceeded 256 KiB. Neither performed numerical work. V4 introduced lossless canonical JSON/zlib9/base64: stored <=256 KiB, decoded <=the existing 4 MiB allocation. Strict envelope/type/base64/zlib EOF/tail/size/SHA/canonical finite duplicate-free JSON and roundtrip checks retain every identity. Controller checks still rehash full inventory/configuration after decoding.

V4 attempt01 failed at the controller's output `mkdir` because its parent directory was missing, before request admission; no controller record was earned (outer a1be00d9). Attempt02 passed required controls and emitted four accepted native payloads, then refused before actual spectral execution because the native map gate omitted transient `/etc/ld.so.cache`. Its controller pre-record timer was 32.61345997895114s; outer launch elapsed was 32.96761010901537s. Inner4505d6a3/outer d2dfb466 and the uncompared raw native prefix remain unchanged.

V5 adds one separately pinned `native_loader_cache`, checked before first child and by a fresh terminal read whether mapped or not. All five immutable SDK libraries plus the new binary must be observed; the exact lexical mapping set may equal those critical paths or those plus this cache. Unknown/aliased mappings refuse; actual paths/file identities retain the observed cache. No data, model, math, criterion, SDK or resource policy changed.

**V5 attempt03 passed:** controller status `accepted-conditional-working-target-engineering`, all controls, full modes and four comparisons; no terminal errors. All four score differences were `5.582442679157156817382E-13` (<1e-6). Reference q=`1.495917324311447356488535831431363476e+03`, M=`-1.943765800475993320642342787607503851e+01`, score=`-7.479586621557236782442679157156817382e+02`; these are identical across the shared four mean vectors. Recorded rho=`3.956354692537067615051392576576527994e-15`, delta=`1.045046046913556910852582948479039260e-13`, kappa=`3.052333740609954298328432287235045806e+03`; response/case scaled gate=`1.454222381088254416856756150542683661e-09`, stationarity=`4.829408374940031020242685467808957280e-18`.

The controller pre-record timer was 35.8845919110463s (before retained-file scan/record write/final seal/stdout), outer launch elapsed 36.28028687502956s, systemd service runtime36.228s/CPU36.513s/peak465.7M/swap0B. All 16 direct children were reaped with observed leader absence; descendant cleanup is not inferred from those checks. The original service cgroup was absent at terminal observation. All 24 controller source identities and 843 outer source/runtime witnesses matched before/after. Native adjusted stationarity and physical source uncertainty remain unavailable.

## Replay requires a new admission

The frozen public proposal cannot run as-is. In a **fresh separate executable request**, retain all science/SDK/input/criteria/resource fields, replace all six source_ports with the actual public-copy absolute paths/bytes/hashes, and fill positively admitted reference_python/reference_runtime and execution_admitted=true only after review. Do not rewrite historical source/receipt pins or edit the public proposal into a past acceptance. Changed original inputs, SDK/cache identity or model needs a separately reviewed successor. Private helper/data paths require lawful acquisition and exact role/axis admission, not path guessing or a replacement covariance.

Use the observed resolved regular Python and explicit primary site-packages/private official threadpoolctl 3.6.0 context; no venv alias relaxation. Before any run, bind reviewed source/bounds and obtain the shared one-job grant with external control-group supervision. A concrete metadata CLI is:

```sh
export PYTHONPATH=/home/szymon/Projects/reproducible/.venv/lib/python3.12/site-packages:/home/szymon/.codex/worktrees/a4f7/reproducible/.work/sn-threadpoolctl-3.6.0-20261003
export PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
/home/szymon/.local/share/uv/python/cpython-3.12.14-linux-x86_64-gnu/bin/python3.12 -B /ABS/PUBLIC-COPY/reference.py --runtime-fingerprint --output /ABS/FRESH-GETTER
```

Independently admit the full decoded inventory/configuration and exact before/after equality, then pin a fresh file containing its packed fingerprint as reference_runtime. Actual source/cache/module/namespace/lexical maps, Python/backend, longdouble/fenv, one-thread getters/environment and sys.path must remain complete. The loaded-process `ctypes.CDLL(None)` libc metadata binding avoids discovery children; packages/pools stay unchanged. A new public-path fingerprint is not the historical private fingerprint. After final executable-request review/hash and a new grant:

```sh
/home/szymon/.local/share/uv/python/cpython-3.12.14-linux-x86_64-gnu/bin/python3.12 -B /ABS/PUBLIC-COPY/controller.py --request /ABS/request.execution.json --request-sha EXACT_NEW_REQUEST_SHA256 --attempt /ABS/FRESH-ATTEMPT
```

These are manual CLI arguments, not a generic JSON command framework. The absolute output parent must already exist as a nonsymlink directory; the attempt itself must be absent. The fixed family is 1800s, one job/thread, AS2 GiB per child, file128 MiB, owned store256 MiB and logs16 MiB. Phase maxima: assembly120s, compile120s, two metadata30s each, controls840s, actual native300s/reference300s, terminal60s; <=16 direct children. Seven Gaussian preparation attempts reserve three condition estimates of at most 1657 inverse columns each, a source-derived upper bound rather than instrumented primitive counts; reserve three eigensolves/six dense checks. Fresh EXCL outputs, pre-child/during-run store checks and failure evidence remain mandatory. The source/output/public-tree caps are unchanged.

## Exact identities and preservation

The request and local full receipts hold complete roles/pins. Key SHA256 identities:

- Table (579283 B): `1cb0fc379ef066afdc2ffd1857681cc478024570d8a3eba284fb645775198cf8`.
- Original covariance (33284960 B): `abf806d966485e64afdb359c87bffc0ecc00d05eff0a31ced66f247385df0fdc`.
- Selection: `e6a3afa4789a6fa77bd977e653b6cdcaac4060d13a5918eec7c9c4ab459a5510`.
- CLASS four-case receipt: `8eb4f479207f53f7889f29ecb395dbdc4304a19d3404623f030b98f3a3340eee`.
- Original selected B / new S_hat: `25beb54b94df0d420eb9117ea3148d15ccb45b54e431ec2d879e0ada8c9f7182` / `07e942f5063e70523206975308ea724b96351b4c0d0b50ed66d53678ca5f040c`.
- Private v5 execution request: `67284ff6217e1e78116e908d8e5415390a7da793ff035a365f44bd3b03a4784c`.
- Private native binary: `f2b123370036719484560c5d35e3bb4bb40746a5706dc76c27047a36bd81397e`; native raw record (identical to v4 prefix): `ef52cde796006ad26850b047089c7729b32023c89ccca213e7b0f8b3a7589cce`.
- Actual spectral receipt: `1f9bc36bff8499ea97bedbb4b364db81df40ae4a9b6b3a6a227f0d1503946d63`.
- Full v5 controller receipt (74504 B): `3eab250c963820512deced760b0c3ab0e692256b3752748f0f63cab333aa447b`.
- Outer v5 seal: `06f5909555afed7970e57a1f5e594b382511d2c5b5efd592f2e8e9cf5e4b76a6`.

The pinned SDK revision is f3b19b6539c72fda5a8b489a0546ca4d0266ace3, build6d495efb166006c6ce651359a366af5d87686ce516ecd7f0f49eed4e8d8be07e. Archive2090606 B/be08158b and admission294226 B/406675b9 are separate identities; observations.hpp remains SHA2514472c with corrected size3622 B. Source commits do not identify new executables.

Full receipts are local under `results/sn-working-profile/full1657-v5-attempt03-20261003/`; private source/runtime/request/failure custody is under `.work/`. No durable archive exists yet. Inputs, third-party dependencies, matrices and full receipts are not redistributed here: review publication permission and provide a lawful exact-input acquisition route and durable research archive before relying on long-term preservation. Private threadpoolctl has official wheel/hash and BSD-3-Clause acquisition/license metadata. A CI pass or successful private engineering comparison supplies no faithful-source, normalized-posterior, H0 or joint-likelihood qualification.
