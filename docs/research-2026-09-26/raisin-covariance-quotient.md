# Does the RAISIN covariance mismatch affect observable distance shape?

Yes. Removing the free common magnitude does not close the mismatch between
the supplied total systematic covariance and the available source-generated
sum. This is a reproduction failure; it does not establish that the supplied
matrix is itself incorrect.

The earlier [RAISIN audit](raisin-differential.md) found a source-linked
mass-threshold propagation defect and quantified a conditional repair of its
systematic variant. Its full-covariance reconstruction failed. Since a common
supernova magnitude offset is normally free, raw matrix equality is stronger
than the equality needed for distance-shape inference. This follow-up checks
that distinction before accepting the remaining reproduction barrier.

For an orthonormal Helmert contrast matrix B with B1=0, a Gaussian likelihood
after eliminating the common magnitude depends on By and BCBᵀ. Its precision
in the original coordinates is

    Bᵀ(BCBᵀ)⁻¹B = C⁻¹ − C⁻¹1(1ᵀC⁻¹1)⁻¹1ᵀC⁻¹.

We compare the published and source-reconstructed systematic covariances in
this 78-dimensional contrast space for the same 79 objects. The source sum
uses unit weights, the actual variant redshifts, and the source's original
inverse-variance centering. Both documented option ranges, 1–28 and 1–29, are
reported. No coefficient, option list or rounding tolerance is fitted.

The [protocol](../../runs/research_2026_09_26/raisin_covariance_quotient/protocol.json)
was frozen before these projected outcomes, SHA256
`a627dddbe72e408f12e3f06fc41606e598b4df9f5b91525e41739f090506d11f`.
It hashes all inputs and the executed script. The known raw-matrix failures
were already visible; this is an explicit follow-up, not a blinded discovery.
The literal covariance/error export rounding bounds T propagate conservatively
as |B|T|Bᵀ|, with a separate 1e−12 numerical allowance.

| Branch | Relative centered systematic matrix mismatch, options 1–28 | Largest mismatch / propagated rounding bound | Published high−low systematic variance | Source high−low systematic variance |
|---|---:|---:|---:|---:|
| NIR | 2.51844 | 1,327,364 | .002626624 | .003096551 |
| Optical | .00472796 | 28,709 | .001773799 | .001781372 |
| Combined | .00217277 | 23,717 | .002011310 | .002017556 |

Variances are mag², for the previously frozen high-37-minus-low-42 contrast.
All six branch/range combinations fail the propagated rounding gate. Adding
option 29 leaves NIR unchanged and increases the optical/combined mismatch.
In particular, the NIR discrepancy cannot be discarded as a common zero point.
The largest variance-direction difference relative to published total covariance
is 9.80 for NIR, while the optical and combined extrema through option 28 are
approximately −.05 to +.05. These describe differences between two artifacts;
they are not uncertainty estimates for the actual supernova measurement.

The independent contrast and centered-projector representations agree within
3.4e−16 mag². The displayed precision identity agrees within 1.14e−13, and an
injected pure-offset covariance contribution cancels within 1.12e−16 mag².
These checks validate the algebraic diagnostic. They do not identify which
historical covariance-generation inputs or operations are absent.

An independent [Astra check](../../runs/research_2026_09_26/astra_design/next_experiment_review/quotient-check/result.json)
verifies all 102 input hashes and reconstructs the NIR 1–28 sum using adaptive
cosmological quadrature and a separate QR contrast basis. The centered matrices
agree within 8.43e−15 mag². It also checks the flat-intercept determinant
identity within 2.04e−14: when comparing changing covariance models, the
likelihood normalization and prior measure must accompany the quadratic form.

The result supports retaining the current boundary: the mass-group repair is
reviewable, but a full corrected covariance and cosmology require the missing
generation provenance. No source or released input was modified.

Reproduce with [raisin_covariance_quotient.py](../../scripts/research_2026_09_26/raisin_covariance_quotient.py),
running `freeze` then `score` in the project environment. Existing outputs are
protected against overwriting. [Result and numerical checks](../../runs/research_2026_09_26/raisin_covariance_quotient/result.json)
and [contrast arrays](../../runs/research_2026_09_26/raisin_covariance_quotient/arrays.npz)
are preserved.
