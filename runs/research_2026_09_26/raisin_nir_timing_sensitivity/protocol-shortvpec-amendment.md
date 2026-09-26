# Prescore technical amendment

Original frozen protocol SHA-256: `3cb9ed02a920ddda6dd9742c607d2a49f672729c9a8b5e2d36c1d0fffcb14aeb`.

The first baseline launch exited −6 before any FITRES row, with a stack-smashing backtrace in `check_file_docana` during `INIT_HEADER_OVERRIDE`. The old NML, log and execution record remain in `fits/baseline/`. The exact buffer fault is not proven; v11_03c receives an absolute `HEADER_OVERRIDE_FILE` in this attempt.

The amended seven NMLs use the short work-local `HEADER_OVERRIDE_FILE='vpec.list'`, each a byte-identical copy of the released `vpec_baseline_raisin.list`. Source, model, KCOR, photometry, cohort, fit settings, timing cases and preregistered gates are unchanged. New isolated work/output paths prevent overwriting the failed run. No successful fit or timing outcome existed before this amendment.

Amended protocol SHA-256: `e42ceec5781fa2ae699ba51ac73d1a6cfe2b0c5d72b9e2bbf6ea5cebc7e489dd`.
