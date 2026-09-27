# Shorter-wavelength infrared measurements of the DES hosts

**234 of the 265 DES C3/X3 hosts have a close catalogue counterpart with reported
3.6 and 4.5 µm fluxes.** The corresponding counts are 50 at 5.8 µm, 47 at 8 µm
and 24 at 24 µm. These are unions across catalogues, with repeated measurements
of the same host counted once. They provide additional observed constraints and
source positions for a future common-aperture analysis. They are not deblended
SPIRE fluxes, dust luminosities, stellar-age measurements or an age correction.

The sample is the unchanged 265-host subset with the earlier eight-band
photometry and valid measurements in all three SPIRE maps: 144 C3 and 121 X3
hosts. The official [SWIRE archive](https://irsa.ipac.caltech.edu/data/SPITZER/SWIRE/overview.html)
(DOI10.26131/IRSA406) and [SERVS archive](https://irsa.ipac.caltech.edu/data/SPITZER/SERVS/overview.html)
(DOI10.26131/IRSA407) supplied 31,820 catalogue sources in the union of
75-arcsec neighbourhoods around these hosts. Every host and all local candidates
remain in the ledgers; catalogue source identifiers remain strings.

## Positional evidence

The primary radius is 1 arcsec for the IRAC-positioned catalogues and 2 arcsec
for the separate MIPS24 catalogue. An association must be unique within that
radius, have the target as the nearest unflagged DES deep-catalogue galaxy, and
have no other such galaxy within the primary radius of its infrared centroid.
The primary associations all satisfy these geometric checks.

| Catalogue | C3 associations /144 | X3 associations /121 | Radius |
|---|---:|---:|---:|
| SWIRE IRAC-positioned | 104 | 84 | 1 arcsec |
| SERVS two-band | 118 | 96 | 1 arcsec |
| SWIRE MIPS24 | 16 | 8 | 2 arcsec |

At 2 arcsec, the IRAC counterpart counts increase to 108/86 for SWIRE and
123/98 for SERVS; one X3 SWIRE host then has multiple infrared candidates.
The complete 1/2/3/6-arcsec comparison is retained, rather than selecting a radius
after seeing its yield. No primary infrared source is assigned to two selected
hosts.

Independent spherical calculations agree to **2.9×10⁻¹¹ arcsec**. Thirty-six
shifted positions per host, at 20/40/60 arcsec, give primary-radius coincidence
fractions of 0.87%/0.76% for SWIRE, 1.25%/1.31% for SERVS and 0.27%/0.30%
for MIPS24, in C3/X3 respectively. The large positional excess is evidence for
real counterparts. Those control fractions are not posterior false-association
probabilities: footprint edges, clustered galaxies and shared sources remain.
The host-bootstrap intervals in the [control result](../results/infrared_photometry/matching-controls.json)
condition on the catalogue and offsets, without field or shared-source covariance.

Association does not imply isolated photometry. **22 of the 24 MIPS counterparts
have another optical galaxy within 6 arcsec.** The combined SWIRE position is
usually IRAC dominated. Its 24 µm flux therefore additionally requires the same
`detid_24` in the independently passing separate 24 µm association; those
duplicate records are not independent measurements. Flagged and optically
undetected possible infrared emitters remain outside the optical competitor
screen.

## Flux consistency and measurement limits

For counterparts passing both catalogue association checks and agreeing within
1 arcsec between infrared centroids, native aperture-2 fluxes give:

| Field | Hosts | SERVS/SWIRE at 3.6 µm: median [5th,95th percentile] | At 4.5 µm |
|---|---:|---:|---:|
| C3 | 91 | 0.998 [0.905,1.061] | 1.005 [0.912,1.243] |
| X3 | 75 | 0.984 [0.909,1.043] | 1.012 [0.877,1.246] |

These are observed ratio distributions, **not confidence intervals**. Two
additional X3 pairs fail the infrared-centroid agreement criterion; their rows
are retained. No valid common pair has a nonpositive flux. The near-unity
medians are useful release consistency, but do not independently validate
calibration: SERVS was explicitly calibrated to SWIRE, including an already
applied 1.02 factor at 3.6 µm in these fields. The surveys use the same IRAC
instrument; individual exposure overlap has not been audited. SERVS incorporates
earlier SIMPLE/GOODS and SpUDS observations, as described by
[Mauduit et al.](https://arxiv.org/abs/1206.4060), rather than establishing that
these particular SWIRE and SERVS entries share photons.

All fluxes and errors retain native microJy values and flags. No additional
aperture, foreground-extinction or K-correction is applied. The
[SERVS documentation](https://irsa.ipac.caltech.edu/data/SPITZER/SERVS/docs/SERVS_DR1_v1.4.pdf)
describes a detection-selected two-band catalogue, already aperture corrected
and with adjusted extraction errors. The SWIRE aperture-2 convention is 1.9
arcsec for IRAC. Its CDFS TAP error descriptions instead list 3.2 arcsec:
**269 associated CDFS IRAC catalogue-band measurements carry this unresolved
metadata conflict**, and no chi-square agreement claim uses those errors.

The [SWIRE delivery document](https://irsa.ipac.caltech.edu/data/SPITZER/SWIRE/docs/delivery_doc_r2_v2.pdf)
also contradicts itself about extension flag −1. This affects one matched host,
SN1292805, in five catalogue-band rows, including duplicated 24 µm measurements.
Those rows are retained and excluded from the strict flags-clear diagnostic.
Flag0 is indeterminate, not certified pointlike; aperture corrections derived
from stars need not recover extended-galaxy total light. For MIPS, the
[release-specific column definitions](https://irsa.ipac.caltech.edu/data/SPITZER/SWIRE/SWIRE_24only_columns.html)
explicitly start at aperture2, radius5.25 arcsec. The delivery document's
unlabelled radius sequence starts with the rounded value5.3; it does not map
aperture2 to7.5 arcsec.

There are 443 matched catalogue-band rows with coverage but missing flux, and
728 host/catalogue-band rows without a primary association. These are neither
zero measurements nor numerical upper limits. No joint interband calibration
covariance was recovered. IRAC contains stellar continuum as well as possible
dust or active-nucleus contributions; these points cannot be interpreted as
dust-only luminosity. Infrared nondetection also does not justify removing a
possible contributor from a future SPIRE image model.

The [compact analysis](../results/infrared_photometry/summary.json) records all
counts, radius sensitivities, flux comparisons and source/output hashes.
The [reproduction instructions](../code/infrared_photometry/README.md) retain
first-party acquisition and analysis code; downloadable catalogues and full
ledgers stay under the ignored `.work/infrared-photometry/` directory.
