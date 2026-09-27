# Released absolute-distance ladder

The [scientific note](../../notes/distance-ladder.md) reports the independent SH0ES matrix reconstruction, SN-free correlated host-distance factor, overlap inventory and integration limits. Its final section gives ordered reproduction commands.

- `acquire.py` verifies pinned official data, including the actual Git LFS covariance.
- `solve.py` solves the complete released Gaussian linear system.
- `cepheid_factor.py` removes all SN rows, identifies 37 hosts against Table 2 and compresses the non-SN likelihood with full covariance.
- `overlap.py` screens current Dovekie/Pantheon identities conservatively; it does not construct a cross-release covariance.
- `independent_review.py` and `validate.py` check independent algebra, source identities and evidence.

No script modifies or runs the cosmological targets. The H₀ summary must not be multiplied into overlapping corrected-distance likelihoods as an independent prior.
