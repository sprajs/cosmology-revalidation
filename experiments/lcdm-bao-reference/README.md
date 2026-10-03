# Separate full-LCDM DESI BAO reference

This packet scores the released thirteen-coordinate DESI DR2 Gaussian BAO
compression against the same explicit full ΛCDM CLASS state used by the
[external reference](../lcdm-reference/README.md). It retains the full covariance
and keeps this score separate from the primary CMB and SN targets. It performs
no cosmological fit and defines no combined posterior.

[request.json](request.json) binds the input/order/state/source and resource
contract. [bao.py](bao.py) admits the exact release bytes, interpolates the full
emitted CLASS background table and forms DM/r_drag, DH/r_drag and DV/r_drag in
the released order. CLASS emits H/c in 1/Mpc, so DH=1/(H/c); DM is the radial
comoving distance only for the explicitly flat state. The final z=2.33 pair is
DH then DM. The ruler is the same state's predicted drag horizon, printed by
CLASS, never a supplied replacement. There are no Python cosmological or
Gaussian kernels here.

The source is [CobayaSampler/bao_data at bb0c1c9](https://github.com/CobayaSampler/bao_data/tree/bb0c1c9009dc76d1391300e169e8df38fd1096db/desi_bao_dr2).
The mean and covariance are released fitted summaries, not original galaxy
observations. [experiment.json](experiment.json) retains their byte identities
and scientific role. Acquisition and redistribution rights remain separate;
no third-party arrays or run products are tracked in this packet.

[consumer.cpp](consumer.cpp) uses the existing Irreducible Gaussian owner:
one preparation of the full 13×13 covariance, then four whole residual vectors
`observed - predicted`. [transport.py](transport.py) checks exact axes, source
and vector echoes, finite domains, work/status/partial availability and the
binary64 log-density identity. The native policy retains the existing `1e-8`
estimated-forward-sensitivity screen. This empirical screen is not a certified
score error or an observational calibration uncertainty.

2026-10-03: a fresh historical attempt C completed all four separate scores in
3.51 seconds. For anchor, changed precision, `n_s−0.005` and `n_s+0.005`,
each quadratic was **31.433594404478523** and each normalized log density was
**−9.15910744477944**. The full-covariance log determinant was
−37.007781378241134 and the positive normalizer was 13 log(2π) =
23.892401863321492. No prior term is included. The empirical forward sensitivity
was 3.8946463578645627e−13; preparation and all four evaluations returned finite
status with numerical status `ok`.

The four full background tables had identical SHA256
`51f02c4e9152bf2f9f009df9541d5e54def9fb9dd18f5a28ff7387f857f2e4c9`
(21,080,729 bytes); their thermodynamics tables likewise shared
`0e088a78921a82fb8b2171e7e8a1008497b50ad757c95e662bdbd3bd07b1dd51`
(8,557,352 bytes). Their printed predicted ruler was 147.054261 Mpc at
z=1059.928342. Spectral-index changes affect perturbation spectra while holding
this background fixed. The precision variant did not change these emitted
background/drag bytes; identical BAO scores therefore earn no refinement bound.
Linear background interpolation and six-decimal drag output still have no
propagated likelihood-accuracy certificate. The original parameter point was
conditioned on a joint Planck primary-plus-lensing result; it is not an asserted
primary-only best fit or independent posterior draw.

The historical SDK was honestly identified as Irreducible
`f3b19b6539c72fda5a8b489a0546ca4d0266ace3`, build ID
`6d495efb166006c6ce651359a366af5d87686ce516ecd7f0f49eed4e8d8be07e`.
The source admission checked all 327 manifest source blobs; the consumer used
46 exact committed headers and the pinned archive. The actual compiled consumer
SHA256 was `7de53db291ac96f8611124d3611828337e8e8339f2bd6ef0685a94bf644217fd`.
This is not relabeled as a later engine build. SDK discovery, source/artifact
and Python runtime identities agreed before/after the attempt.

Historical record SHA256 is
`97dd8c95ddd5021ef387e34628b24e8b2580d6fae788f8cee91da967fef31b33`;
the immutable courier is
`e1367af247403dbe8d2ccf8a6e9672d9c0e4cd351eb98dba63f1e080d6a07e0c`.
They remain under ignored local `results/lcdm-reference/` and `.work/` stores;
a durable publication archive has not been created. Failed attempt A admitted
an unnecessary oversized historical SDK test log; B compiled successfully but
refused the exact release header. Both failures remain immutable. C used the
narrow exact optional leading-header repair; five strict header/order/full
covariance controls passed. Arbitrary, repeated or nonleading comments refuse.
A separate root-owned arithmetic check on 2026-10-03 passed for all four
states, using independent pivoted Gaussian elimination, all thirteen positive
leading principal minors, and Machin π in Decimal at 80 and 120 digits. The
absolute native discrepancies were 2.2967e−15 in q, 8.4395e−17 in logdet,
7.8631e−16 in normalization and 1.5837e−15 in log density, below the fixed
`1e-8` comparison gate. The largest 80→120 drift was 1.66e−78, below `1e-10`;
the independent analytic det=35, q=12/5 controls also passed. This checks the
represented covariance and rounded residual arithmetic, not CLASS predictions
or interpolation. The four states share one identical numerical corpus.
Proposed native negative controls remain unexecuted. Independent checker source
SHA256 is `53fc00c52b1fd7ba23b5956a717a337cf7517882f99fe2f41f773ee7ea841ddd`,
and result SHA256 is `7f74e08ea555af3f1fe284b590dd4372f7385a6627ec777ad46fa83582626420`.
Shared input ancestry is retained; no physical or joint inference gate closes.

The adapter, transport and C++ consumer here preserve C's executed source bytes.
[prepare_probe.py](prepare_probe.py) and [controller.py](controller.py) are
distinct successors: they replace workstation literals with repository
relative sources and an explicitly admitted local path/identity configuration.
They reuse the exact existing [CLASS table adapter](../lcdm-reference/theory.py)
and [process capture](../lcdm-reference/run.py), loaded from hash-checked RAM
bytes. The old full CLASS tables and exact INI/receipt law remain the target;
rebasing a historical path is not automatically admissible. The pinned SDK
admission likewise retains its original paths and source identity. A fresh clone
needs those lawful local artifacts or a separately reviewed new identity.
Historical execution does not qualify this new controller source. The immutable
request records its initially unexecuted design; the later execution below
retains that request identity.

2026-10-03: the committed public route at
`17d6b0213996fcb8577c4dc5464e13928a8850db` completed a fresh bounded confirmation.
Its reported controller wall before final sealing was **3.691418 seconds**;
describe, compile and native evaluation all exited
0. Controller SHA256 was
`9d9573290e98684d7d5fdf5e95c7420c359c4cf10e6d599c833e39ea7f73fb96`;
preparer SHA256 was
`d084193d7568cb465e42573e39430acd09941e3fbc2511541baede727e3c7a31`.
The receipt, retained locally at
`results/lcdm-bao-reference/public-source-confirmation-v1-20261003/record.json`,
is 325,957 bytes with SHA256
`ab2c077d8b84123bffba774deaa660e55c24247d744024b27771ce52443502bc`.
The compiled consumer, native wire and raw native output were byte-identical to
C. All four complete prediction and admitted native-result objects agreed
exactly with C, including means, residuals and scores. All 105 consumed-file
identities agreed on common hash/stat fields before and after; actual Python
runtime inventories matched. No failures or sealing errors occurred. The
portable 14 synthetic header/authority controls and the full 247-test suite
passed; frozen dependency sync and packet/storage/link checks passed. This
confirms this configured source route, while the arithmetic, CLASS accuracy and
inference limits above remain distinct.

A separate final source review found that refused changed files lost their
observed identities. The correction retains consumed-stream byte counts, hashes,
EOF state and descriptor/link observations in initial and terminal failure
records, including failed stat capture. It preserves expected authorities and
never labels a changed stream as the current complete file. Controller SHA256
is now `76d53f71e53cbae91ff5bacde26029c0dc5f846ee3ba7773ec8f4f5accda3fb5`.
The 17 synthetic controls pass, including changed-byte/read-drift, missing-link,
terminal-stat and serialized-failure evidence.

2026-10-03: this correction at committed source
`cfd24e19d1a71ba8b4ab5bba73e264ae3ec61cc3` passed its separate bounded
confirmation. The reported pre-seal controller wall was **3.536021 seconds**;
describe, compile and native evaluation exited 0, with no failures or sealing
errors. Receipt
`results/lcdm-bao-reference/public-source-confirmation-v2-20261003/record.json`
is 325,927 bytes with SHA256
`93081f47596d11320bdb57a0b5ba3825434679c3b21fc003bb7a8b8ce67bf760`.
The closed local configuration was 981 bytes with SHA256
`9e042abdb1612a64802c416561be85979f43ea4db56281cd7660dc0af1fb7959`.
All four complete prediction and native-result objects matched both C and the
first public confirmation exactly. Compiled consumer, native wire and raw
native output likewise retained their original byte identities. All 105
consumed-file identities agreed on common hash/stat fields before and after,
and actual Python runtime inventories matched. The full 250-test suite passed.
The earlier receipts and source-review findings remain immutable; this new
execution closes the corrected launcher gate while scientific limits remain.

After complete source/runtime review and a one-job grant, the dedicated route is:

```sh
env -i PATH=/usr/bin:/bin LC_ALL=C PYTHONDONTWRITEBYTECODE=1 \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 \
  VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  timeout --signal=TERM --kill-after=2s 298s \
  /path/to/admitted/python -B experiments/lcdm-bao-reference/controller.py \
  --config /absolute/reviewed-local-runtime.json --attempt fresh-label
```

The closed local configuration has exactly `schema`, `controller`, `source_pins`,
`sdk_directory`, and `native_archive`; schema is
`lcdm-bao-reference-local-runtime/v1`. File authorities contain absolute `path`,
integer `bytes` and SHA256. The source-pins file preserves the request's 26
ordered source roles and raw byte identities; its concrete paths, the controller
pin and archive pin require explicit review. It supplies no commands, model,
numeric policy or seeds. The current Python/tool identities are the historical
reviewed ones; no installation or environment creation occurs here.

Execution is serial: describe and compile within 120 seconds, then retained-table
preparation, four native scores and terminal checks within 180 seconds, at most
300 seconds total. The process uses one job/thread, 2 GiB address space, 32 MiB
fresh attempt storage, bounded raw logs, kernel phase timers and exclusive
readonly outputs. Failed and partial objects remain available; terminal source
or runtime drift withdraws acceptance. Per-process CPU limits are not an
aggregate compiler-descendant bound, and descendant cleanup is not proved.
The generic packet runner remains blocked. Calibration, end-to-end numerical
accuracy, inference and cross-probe covariance/prior gates remain unqualified.
