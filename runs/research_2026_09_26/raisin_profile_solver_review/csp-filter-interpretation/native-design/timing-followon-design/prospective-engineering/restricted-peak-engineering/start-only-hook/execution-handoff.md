# Start-only build and staged executor handoff

The reviewed patch compiled once in a separate private tree in15.795086334seconds, within the authorized120seconds. No native fit or photon generation was executed by this preparation task. The prior restricted binary/source and every failed/partial run remain unchanged.

New binary: `phase2/pte/restricted-peak-engineering/start-only-hook/build/bin/snlc_fit.exe`, SHA256 `5f8fb0add496c5c10400138ccb97d0c2d33d656fc9482c3a5c7b05ebb1badfcc`. The old restricted binary remains `4ba9c40660fd1d1df42dd300b64f56c0aaa42e1d6e5c891c6d4b56c60737cf14`. All actual source hashes match the prebuild freeze.

The executor is frozen as:

- `execution-protocol.json`: `98e4453f08bd32410d74aac0435c15ad4c3e241a3e0521a5827c1f6cd05af599`
- `execution-freeze.json`: `27af53c9c3b908c7f4e351565a5761b533ce78bad74b9a738576af515905feff`
- `run_start_only.py`: `8c7eba951c14a6edde64c84bbb1d100160bec0579dd4ada3d4993c24b2cfba12`

The freeze covers3568files, including all original inputs, exact baseline joint12/joint9 saved states, private compiled source/binary, inherited pure checkers and the carried resource ledger. Eleven standalone MINUIT-definition/readback tests and six read-only/synthetic executor tests pass. The generalized full-identity checker accepts complete nine-iteration coverage; its only change is an explicit expected-iteration argument.

Root releases the following stages separately:

1. `identity`: four new disabled controls, absent/zero offset for joint12 and joint9. Full CSP entry/model/error/W/objective/FITRES/HEAD/PHOT/support/domain identity against the corresponding saved baseline is mandatory. This establishes whether those baselines can be reused with the new private build.
2. `starts`: two new actual MINUIT entry−2/+2 checks with **byte-identical nominal NML initialization**. Common INIVAL and the first prior center remain unchanged. Actual MNPOUT stored values, all iteration records, first CSP_ENTRY and pre-grid bound records are checked before accepting the original numerical/support gates. Source/prehook state is shared; W evaluated at different parameter points may differ.
3. `nir`: nine remaining NIR true/fitted-peak12/9/postinit-D-start/null jobs, using the original nominal joint12 peak. No shifted-start winner is selected. The offset environment is cleared. The full null model/weight/objective/FITRES comparison and all original gates precede the engineering paired-effect table.

The activity ledger carries **11.42092445801245seconds** already used into the existing120-second restricted-estimator allowance. Each new native process is capped at40seconds and uses one worker. Failures remain charged and preserved. No source support, distance, peak, covariance, stationarity, fixed-mask, final±4day or actual-start threshold is relaxed. No new photons, reseeds, deleted objects, automatic64-draw study, survey correction or cosmology calculation follows from this handoff.

Example, only with the corresponding root release file:

```
phase2/env-official/bin/python phase2/pte/restricted-peak-engineering/start-only-hook/run_start_only.py identity --release <root-release.json>
```

The release binds both frozen protocol and freeze hashes plus the permitted stage. All new fits are under root-level `P/fits-start-only/`; new HEAD adapters are under `P/derived-start-only/`, preserving the validated relative NML path depth. The executor explicitly removes inherited `PROSP_MINUIT_PEAK_SHIFT` for every default/NIR call, setting it only for declared identity-zero or±2 jobs. Its output files are opened exclusively and no attempted case is silently overwritten.
