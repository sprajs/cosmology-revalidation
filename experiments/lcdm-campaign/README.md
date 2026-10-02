# Standard LambdaCDM reproduction campaign

The tested question is which of the nineteen real-input/source and compiled
physics contracts support a faithful standard LambdaCDM reproduction. The full
Planck model remains **blocked**. [campaign.json](campaign.json) freezes all
nineteen questions, dependencies, references and missing contracts from the
coordinator audit; [candidate.json](candidate.json) preserves the Prospector
source design. This status snapshot serves the
[sole active Irreducible roadmap](https://github.com/sprajs/irreducible/blob/main/docs/roadmap.md).
It does not establish another work queue or promote proposed questions to findings.

Reproducible owns exact inputs, experiment transport, references and attempt
records. Irreducible owns production physics. Prospector owns source review.
The original BAO, ladder and SN packets keep their original requests/SDK pins.
No Planck posterior parameters enter as additional observations. Unknown
cross-covariances remain unknown and no combined likelihood is executed here.

## Current bounded consumer

[request.json](request.json) defines a new experiment/build identity, separate from
the historical massless BAO packet. [controller.py](controller.py) admits the
verified clean pinned engine source, its actual build manifest/archive/CLI,
all 302 manifest source hashes and compiler/standard-library identities. It reads
original inputs without changing them, snapshots the source and transport,
compiles one thin native consumer and compares all thirteen DESI DR2 predictions
and normalized density components with newly evaluated references.

The two variants are explicitly chosen physical-density massless and massive
FD controls with supplied drag. They do not complete Planck's species mapping,
recombination, predicted drag or perturbations. Data are released fitted
compressions under their selection/reconstruction/fiducial assumptions, not raw
galaxy observations. Full ordered covariance remains intact. Native distances
and ruler use one retained thermal state; no production equation is implemented
in this workspace.

[reference.py](reference.py) is reference-only mpmath1.3.0 at 60/90 digits with
independent direct momentum quadrature, direct-redshift distance integration,
sqrt(a) ruler integration and scalar LDLT. Its strategy shares ancestry with
Irreducible's original thermal peer references. Shared physical equations and
constants are explicit ancestry, not independent observations. The fresh
reference compares its decimal strings against exact native binary64 values.
Two additional 90-digit cross refinements isolate momentum and outer integration;
an analytic exponential-tail envelope is propagated through ratios and the full
covariance density separately. These refinement comparisons are empirical
discretization checks. Runtime fingerprints bind Python and the complete mpmath
source inventory before and after execution. Reference refinement and the summed
axis/tail envelope must consume at most 5% of fixed allocations. Ratios use
1e-10+5e-10*abs(reference), density components absolute1e-8, projection <=1e-8.
The producer uses separately frozen tighter settings and finite work limits.
No failed row is removed, covariance inflated or acceptance tolerance relaxed.
The native consumer returns ratios and density only. Its reference ruler output
does not exercise the request's native distance/ruler comparison allocation.

One compiler/scientific job and one test thread run at a time. This campaign
uses the verified primary native archive read-only; it does not rebuild or
install an engine behind its owner. Existing installed SDKs predate this build.
The operation name is local to this controller and is not an Irreducible CLI
operation. The generic packet runner refuses the blocked full experiment.

```sh
uv sync --frozen
uv run python experiments/lcdm-campaign/controller.py \
  --engine-source .work/irreducible-source-e9a9e6bd \
  --engine-artifacts /home/szymon/Projects/irreducible \
  --input-root /home/szymon/Projects/reproducible \
  --reference-python /home/szymon/Projects/irreducible/evidence/project-review/science/recovery-abundance-growth-20261002/reference-venv/bin/python \
  --name fresh-thermal-desi-control
```

The paths describe local verified evidence, not a downloadable SDK release.
A new host/build needs a separately reviewed identity. Runs create fresh ignored
`results/lcdm-campaign/` directories, preserve failures and make records read-only.
Original third-party bytes remain local: release licenses/publication permission
have not been established, so no input assets are redistributed. A durable
archive is required before a scientific result can be published.

## Findings

2026-10-02: admitted the full nineteen-scope campaign as blocked and froze
source-fixed44 box normalization and thermal DESI contracts. A preserved
uncommitted diagnostic compared all thirteen rows for both supplied-drag models
at engine source `e9a9e6bd5af3404c7efc66916dc6ed6c68d72b4f`, build
`e92a0e64ee8e7101864a2d221b67d02abe434ad6f8288f84475b15a4a2445999`.
Its native projection estimates were below1e-8, but the controller converted
60/90-digit references to binary64 before comparing refinement. Its recorded
numerical pass is provisional: small differences were erased. The immutable
diagnostic receipt SHA256
`24c4be5c0f1b29ef988d4d56b9583e13b081500ff7c14bd1e4c8eac184ab9d6e`
remains local and unchanged. A corrected committed comparison and deliberate
review are required before acceptance. Every remaining closure stays visible;
a passing control will not promote the full experiment.
