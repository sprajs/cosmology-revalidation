# Frozen ten-object NIR timing sensitivity

Frozen before this source clone's first fit: `protocol.json` SHA-256 `3cb9ed02a920ddda6dd9742c607d2a49f672729c9a8b5e2d36c1d0fffcb14aeb`.

Use SNANA v11_03c commit `06f2ccfdbf99d62c23dec66d0bf7b7e9452604d2`, built in an isolated clone with existing local compiler/libraries. Its only tracked build-generated source change disables unavailable HBOOK/ROOT flags; the patch, source hashes, build log and binary hash are frozen. Inputs are the ten released DES16 photometry files, released NIR JH-only NML settings, released KCOR/model and VPEC override. The seven immutable condition directories change only both `PEAKMJD` header lines per object. `author_t0` uses the pinned author file's unrounded values; shifts are observer-frame days from the released header. The baseline copy uses byte-identical photometry but a separate data/work directory.

Run baseline and its copy first. Require ten exact FITRES rows, byte-identical LCPLOT, ERRFLAG0, fixed peak/shape/AV/RV and accepted-epoch NDOF closure. Then compare all ten baseline DLMAG to author raw NIR FITRES at ≤0.001 mag and require exact NDOF. If any condition fails, preserve the outputs and diagnose; do not run timing variants. Only after both gates, run −1, −0.5, +0.5, +1 day and author-t0 headers, retaining any failure and accepted-mask changes.

This measures the conditional fitted-distance response to an imposed header date. The original header estimator, any true offset, selection/bias-correction recomputation and historical binary identity remain unknown. `INIVAL_DLMAG` is not separately varied: source initialization may replace it, so no broader optimizer-start or common-objective claim follows.
