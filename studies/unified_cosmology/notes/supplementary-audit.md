# Supplemental result verification

The original finite calculation stopped when its quantile consumer rejected two legitimate repeated selection slots. Its failure record remains unchanged. The separately recorded continuation retains those slots and their weights, and completes the requested downstream calculations. A completed invocation does not establish a scientific result: each output keeps its own qualification or failed gates.

The supplemental auditor checks the original plan and failure hashes, the explicitly named continuation plans, each original command and producer, and every completed attempt, log, result and completion receipt. Only substitution of the independently validated ordered-quantile consumer is allowed. A missing completion receipt leaves its outputs pending. Duplicate receipts or outputs, altered commands, unbound parents, changed sources and unknown scientific statuses are rejected. Logs lack a producer-supplied digest; their bytes are bound at audit time and checked again before return.

The auditor calls the unchanged parent qualifier and recursively checks child source, data, cache and parent hashes. It reports the original execution and stage states separately from effective supplemental states. A completed omission calculation with insufficient importance overlap remains failed. Histories qualified at the declared accuracy appear in `verified_conditional_histories`; only histories with a supported finite native-precision diagnostic also appear in the plotting interface `verified_children`.

Scientific target identities contain logical asset groups such as `primary/planck_2018`. These are hashes of declared file inventories, not filesystem paths. The audit checks the complete group mapping against both inventories, then verifies each underlying file's hash and size. Ordinary source and data path manifests retain their file checks. Synthetic controls reject altered or missing groups, asset bytes and ordinary files; the earlier snapshot that failed on a logical asset name is preserved.

The native 32-point summary embeds each sealed original record with five additional file/log annotations. The original seal applies to the original file, not that annotated view. The audit recognizes only this exact projection: it verifies the original file and seal, unchanged payload, selected point and producer lineage, and the same-index log and warning count. Extra or altered fields, resealed altered views, and changed files are rejected. All other payload seals retain their ordinary strict checks, and original failed scientific checks remain failed. The earlier audit snapshot that treated the annotated view as a new sealed record is preserved.

An optional thermal-path review must be explicitly named and hashed in a new design manifest. The auditor then calls `native_precision_review_consumer.verify` afresh, retaining the original numerical flags and the reviewed decisions. Only its supported diagnostic permits native-accuracy-supported plots or lensing follow-ups. Genuine likelihood variation or other surviving failures continue to block that interpretation. The 32-point screen remains a numerical diagnostic, not posterior qualification or a guarantee of accuracy everywhere.

The default design names only the two LCDM continuation plans. There is no trusted directory discovery. A new cohort or thermal review requires an explicit versioned manifest supplied with `--design`; an optional saved pure-consumer receipt must equal the freshly verified result.

From the repository root, validation and a subsequent read-only audit use:

```bash
.work/unified-cosmology/external-probes/.modern-venv/bin/python studies/unified_cosmology/code/supplementary_audit_validate.py
.work/unified-cosmology/external-probes/.modern-venv/bin/python studies/unified_cosmology/code/supplementary_audit.py --output .work/unified-cosmology/supplementary-audit-snapshot.json
```

The output path must be new. Synthetic validation covers qualified, failed and incomplete receipts; exact producer substitution; duplicate and tampered evidence; parent rejection; and repaired versus genuinely failed precision screens. It uses no physical model or background calculations and does not run an observational audit.
