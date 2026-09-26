# Independent generation-stage failure review

**STOP: the first original-generator process failed before writing event photometry.** It returned native fatal-error status 11 after 17.047889647 seconds. The two FITS output files are both zero bytes, the DUMP contains only its header, and the ledger/noiseless branches were not run. The 120-second budget therefore has 102.952110353 seconds remaining; the failed work must not be discarded from that accounting.

The only prepared cadence is LIBID11, but the frozen input specifies `SIMLIB_IDLOCK: 1`. Exact v11_04d `SIMLIB_findStart` ultimately copies that literal input to `GENLC.SIMLIB_IDLOCK`. `keep_SIMLIB_HEADER` then rejects LIBID11 because it differs from lock1. The later main-loop special case that would lock the first accepted LIBID cannot be reached. The log explicitly prints `IDLOCK: 1 (will use only this LIBID)` and aborts after `NTRY=100000`.

The earlier static review missed this source-control interaction. Its design/hash checks remain reproducible but did not establish native readiness. This failure is not a photometric noise discrepancy, historical-simulation defect claim, or timing result. Generated-versus-reported error algebra and original/instrumented equality are currently untestable because there are no emitted epoch measurements.

An explicit lock on the already selected LIBID11 would implement the intended metadata-defined cadence without changing the physical template, observing conditions, population, selection or seed. It nevertheless requires a separately preserved source-control amendment and release; no input or executor was modified and no retry was invoked by this reviewer. All original failed outputs and hashes are retained.

Evidence: `failure-review.json` and the independent source/output checker `audit_failure.py`. No native invocation or executor import occurs in this review.
