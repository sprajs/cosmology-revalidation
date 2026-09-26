# Direct SN–transverse-BAO comparison: outcome-free design

A direct comparison can constrain the **total relative brightness drift** of
standardized supernovae against a geometric distance tracer without fitting
ΛCDM. It cannot by itself attribute that drift to dust, intrinsic luminosity,
photometric calibration or selection. This design reads redshifts, identities
and covariance, deliberately omitting the SN magnitude column. No observed
drift, opacity or corrected cosmology is reported.

For this diagnostic define

    r_SN(z) = m_b_corr − 5 log10[(1+zHEL) zHD],
    r_BAO(z) = 5 log10[(DM/rd)/z].

Under standard distance duality and a constant comoving BAO ruler,
r_SN−r_BAO is a common magnitude/ruler intercept plus any total relative
brightness drift. No absolute M, H0 or rd is fixed. Heliocentric redshift enters
the photon redshift factor and the corrected Hubble-diagram redshift enters
the distance coordinate, as in the existing SN branch. Constant luminosity or
opacity offsets remain unidentifiable.

This uses the published corrected Pantheon+ magnitudes as the eventual
estimand, retaining their standardization/BBC assumptions. The design uses
their full total covariance, including correlations between repeated rows of
the same physical SN. It marginally selects transverse BAO entries; radial
distances are not used to condition or tighten them. The existing BAO
[reconciliation](bao-and-identification.md) instead compared SN against a
particular expansion history. This is a distinct, more direct comparison.

The [DESI DR2 primary release, section III.3](https://arxiv.org/html/2503.14738v3#S3.SS3)
reports template-derived distances at effective redshifts. Those compressions
are supported by its validation tests, but are not arbitrary-redshift delta
function measurements. The overlapping QSO bin spans .8–2.1. Its response to
rapid distance evolution cannot be inferred from the compressed covariance.
The [distance-duality literature](https://arxiv.org/abs/1304.7317) illustrates
the additional assumptions required to interpret a distance difference
specifically as photon loss. We retain the broader total-drift interpretation.

The [frozen protocol](../../runs/research_2026_09_26/sn_bao_design/protocol.json)
specifies local GLS fits to r_SN, quadratic in ln(1+z), with halfwidth .10.
Quadratic halfwidths .075/.125 and a cubic at .10 are fixed sensitivities.
Their cross-node covariance is W C_SN Wᵀ. Each node must have at least 20
distinct SN identities, at least five strictly on either side, and stable
polynomial rank in **every** design. Eighteen prespecified positive-expansion
curves test interpolation recovery to .005 mag. This finite curve set is an
approximation check, not a smoothness theorem or a cosmological prior.

The data contain 1,590 selected observation rows for 1,473 distinct SN
identities. Five transverse BAO nodes lie within the SN redshift range:
.510, .706, .934, 1.321 and 1.484. Only **.510 and .706** pass every fixed
support gate. The .934 node has 29 distinct SNe in the primary window but only
12 in the narrower window. Higher nodes have insufficient support. We do not
widen a window or relax admission after seeing this result.

For the admitted pair, a Gaussian forecast gives a **.0451-mag standard error**
on the relative drift between the two redshifts; fixed interpolation
sensitivities span .0423–.0469 mag. Under the stated covariance assumptions,
the two-sided 5% test would detect endpoint drifts with the following power:

| True drift across .510–.706 | Primary forecast power |
|---|---:|
| .02 mag | 7.3% |
| .05 mag | 19.8% |
| .10 mag | 60.2% |

Thus agreement in this particular test would have little ability to certify
the absence of a 20–50-millimag bias. Nor does this redshift span constrain
drift below .51, where nearby-to-distant survey calibration is especially
relevant. The forecast is not a measured interval or physical-null calibration.
It assumes zero cross-probe covariance and the supplied Gaussian errors.

The 96/192-node quadratures agree within 1.96e−14 mag. Finite-family
interpolation errors on admitted nodes are far below the .005-mag gate.
Fifty thousand Gaussian BAO realizations check the logarithmic delta method:
its covariance differs by .50% in Frobenius norm, consistent with the scale
of Monte Carlo variation; second-order log-mean shifts are separately retained.
This does not validate galaxy extraction or the redshift compression.

Protocol SHA256 is
`0f0e3ee4a6c75734183ae3c8617409de7a3197e17008a03e31e6eaa8e5666a7d`.
An explicitly archived pre-execution amendment fixes a dictionary-iteration
error caught in code review, without changing the scientific specification.
Original source bytes are preserved. [Design table](../../runs/research_2026_09_26/sn_bao_design/design-table.csv),
[forecast and numerical checks](../../runs/research_2026_09_26/sn_bao_design/result.json),
and [script](../../scripts/research_2026_09_26/sn_bao_design.py) are reproducible.

An independent [raw-input verifier](../../runs/research_2026_09_26/sn_bao_design/independent-check.json)
reconstructed all 20 local GLS operators with whitened QR (maximum weight
difference 1.39e−15), checked quadrature with adaptive integration, and
reproduced admission and endpoint precision through a separate contrast
calculation. It also verified that the script reads no SN magnitude column.

This remains a secondary avenue. A valid observed comparison must specify the
domain in which BAO compression is adequate and retain interpolation
uncertainty. The optical–infrared amplitude-free flux experiment remains more
directly useful for separating chromatic extinction from inherited distance
information. Neither route can identify unrestricted grey luminosity evolution
from SN colours alone.
