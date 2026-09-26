# A same-object NIR phase contrast is useful, but the released cadence fails the frozen feasibility gate

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../README.md) for current execution status.

A same-supernova, same-physical-filter early-versus-late amplitude comparison can isolate a **phase-dependent discrepancy** without using a cosmological distance. It removes each object's constant distance/gray-luminosity factor and a constant filter zero point. It cannot establish a gray luminosity correction or distinguish every source of a phase-dependent discrepancy. The currently listed CSP light curves fail the predefined cadence threshold, and their independent timing, epoch-selection and training-provenance gates remain open. **No observed flux residuals, fitted phase contrasts or native fits were evaluated here.**

This is distinct from reproducing the published between-supernova Hubble-residual split. RAISIN section V.3 reports a .164±.051mag late-versus-early sample difference, a stretch imbalance, and smaller template-phase offsets; its illustrated rest-frame magnitudes use cosmological distances and SNooPy K-corrections. The authors already note shared training data between the two NIR models. None of those findings is new here. The proposed observable instead compares phases of the same object before any Hubble diagram. [RAISIN paper, V.3](https://arxiv.org/html/2201.07801v2#S5.SS3).

**Filter-provenance refinement:** the archive-level test below finds that the supplied SNpy J/H labels already combine Swope and duPont observations. The original8/6 result is therefore a release-token cadence count, not a certified count of distinct physical-filter pairs. Keeping raw physical filters separate gives7Y/5second-band objects among the71 fully linked DR3 objects;2012fr has noDR3 lineage here. The primary cadence gate still fails.

## Frozen metadata-only test and outcome

Before counting cadence support, [protocol.json](../specifications/experiments/raisin_phase_identification/protocol.json), SHA256 `552c1cbe8358961db26fdc19d101eb55a83a4d57c325baa8866dc5aef6dd445b`, defined:

- Primary observer filter **CSP-Y**, early rest phases **[−7,+7] days**, late **[+10,+20] days**, relative to the listed header peak only for this inventory.
- Each window must contain at least three distinct observing nights spanning at least two rest-frame days. Multiple measurements on one night do not count as independent nights.
- Require at least 15 eligible Y objects, including at least ten with the same coverage in a second exact J, j, H or h filter.
- Do not merge Y/y or J/j: their instrument/filter identities differ. Released KCOR input (`sources/repos/djones1040__RAISIN_DataRelease@a383c4b/kcor/kcor_CSPDR3_BD17.input`; historical local reference).
- J/H and the predeclared late[15,30] support ledger are diagnostics, never automatic replacements for a failed primary window.

The parser accesses only SN identifiers, zHEL, header peak, MJD, filter, field and column availability. It does not convert FLUXCAL, MAG or their errors. It uses the authoritative 76-file SNANA LIST, preserving 12,684 listed epoch records. Two other physical `.DAT` files are unlisted and are not silently added or double-counted.

| Metadata quantity | Result |
|---|---:|
| Listed CSP objects | 76 |
| Exact Y objects passing both windows | **8**, required15 |
| Those also passing a second NIR exact-filter pair | **6**, required10 |
| Exact J / H eligible objects, separate diagnostics | 10 / 11 |
| Exact y / j / h eligible objects | 0 / 0 / 0 |
| Y-eligible objects with the optical cadence proxy | 8 |
| Y-eligible names present in the M20 training roster | 6 |

The Y-eligible IDs are **2005el, 2005kc, 2006mr, 2007A, 2008gl, 2008hv, 2008ia, 2012fr**. A separately labelled membership-only check finds that six are in the published nominal w sample; all six are among the M20 name matches. The two nonmatches, 2006mr and2012fr, are not in that published sample. This does not imply they are outside SNooPy training or within its intended population support. No membership cut was imposed after counting support.

The eight optical cadence proxies require three optical nights in[−10,10], with observations on both sides of the header peak, plus three nights in(10,40]. **This is not a timing-precision or independent-timing result.** No thresholds were relaxed after the failed cadence gate. The [inventory result](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/astra_design/raisin_phase_identification/inventory-result.json), exact-filter and epoch ledgers, source hashes and membership labels retain the complete accounting.

## Identification: what cancels and what does not

For one object and exact observer band, write a conditional flux model as

\[
 y_b(t)=a_b\,h_b((t-t_0)/(1+z);s,A_V,R_V,\ldots)+b_b+\epsilon_b(t).
\]

Here `a_b` is free: it includes distance and any time-constant multiplicative calibration/gray factor. `b_b` is an additive template-subtraction/background offset. A late-phase magnitude perturbation θ can be parameterized by multiplying h by `exp(−k θ q(t))`, with `k=ln(10)/2.5`, q=0 in the early window and q=1 in the late window. This is a precisely stated relative-amplitude observable; a smooth transition across the unused gap can be fixed before scoring.

With a_b profiled or marginalized jointly, θ does not require an absolute luminosity, H0, expansion model or peculiar-velocity distance correction. Separate early and late amplitudes yield `θ=−2.5 log10(a_late/a_early)` only in an interior positive-amplitude regime; low-SNR cases require the flux likelihood and proper priors, not logging noisy ratios. Both amplitudes share nuisance parameters and covariance and must not be treated as independent fitted distances.

| Possible cause | Independent information or control | Remaining limitation |
|---|---|---|
| Distance, gray luminosity and constant same-filter zero point | Same-object, same-filter contrast cancels their common multiplier | They are unmeasured by this test; a null cannot rule out gray redshift evolution |
| Peak-time error | Optical-only timing likelihood, excluding both tested NIR windows; compare its derivative pattern across phases/bands | A timing shift can imitate a phase offset when the light-curve slopes lack sufficient variation |
| Stretch/secondary-maximum variation | Optical-only shape/color evolution and independently reserved NIR phases, with a joint nuisance likelihood | Intrinsic phase diversity and a mean-template error are not separately physical without a population model |
| Static dust | Constant monochromatic attenuation cancels; model/project the phase-dependent broadband dust derivatives | Effective broad-band extinction changes with the evolving SED; free intrinsic colors can remain degenerate with dust |
| Additive host/template subtraction | Independent reference-image/fake-source measurements, or explicit per-band additive nuisance with sufficient phase leverage | A constant flux offset does not cancel as a constant magnitude; absent constraints it may absorb the phase contrast |
| Calendar-dependent calibration | Local standard-star observations and night/instrument zero-point covariance; compare phases at different calendar dates | Same-filter cancellation covers constant calibration only |
| Missing faint/negative epochs or follow-up scheduling | Full attempted-epoch/limit/flag record and a selection/censoring likelihood | Later phases can be preferentially lost; same-object matching does not remove epoch selection |
| Population transport to high redshift | Common support in independently measured shape, color, host and spectral properties, followed by matching selection simulations | An internal low-z discrepancy does not establish the high-z distribution or a universal correction |

A useful design gate is the nuisance-projected information. Whiten by a fixed measurement covariance; let v be the θ derivative and Z contain free per-band amplitudes, timing, shape, dust and additive offsets. Then `Iθ = vᵀ[I−Z Z⁺]v`. If Iθ vanishes, better optimization cannot identify θ. An optical nuisance posterior can instead be integrated, preserving its cross-covariances; report both unrestricted projection and the explicit-prior result rather than silently borrowing prior precision.

The protocol freezes future design-only thresholds: optical timing standard deviation≤.5day, at least20% of raw contrast information remaining after nuisance projection, pooled conditional σθ≤.01mag, effective object count≥10, and no object supplying more than20% of information. These have **not** been evaluated or passed. Quoted errors alone do not establish the physical covariance. Model-only injections must separately recover a pure timing perturbation, stretch perturbation, additive background and a .05mag phase perturbation without mislabelling one as another. A failure is nonidentification, not evidence of a phase error.

## Local inputs and gates that prevent an observed likelihood now

The lowest local inputs are the original listed epoch photometry under `sources/repos/djones1040__RAISIN_DataRelease@a383c4b/photometry/RAISIN/CSPDR3_RAISIN/`, together with the exact CSP KCOR file, native SNooPy grid and their hashes. The released CSP columns contain magnitudes/derived FLUXCAL and quoted errors, but no per-epoch PHOTFLAG or complete attempted-exposure indicator. The prior sign audit found positive-only published values, but unlike DES it did not establish the CSP signed ancestor or retention rule. Positivity alone does not prove the DES mechanism also occurred in CSP. Ordinary Gaussian flux scoring cannot be justified merely by this metadata inventory. [Existing sign audit](raisin-flux-sign-audit.md).

The CSP README explicitly says **PEAKMJD was obtained from optical+NIR SNooPy fitting**. Using it as an independent fixed clock would leak the tested NIR observations into the nuisance estimate. It is used here only to inspect provisional cadence. Before any outcome score, fit timing/shape from a declared optical-only data stream with free amplitude and no redshift-distance prior, or obtain an execution-linked independent timing likelihood. Preserve its full uncertainty and validate the mapping between an observed-band peak and the template's rest-frame B peak. Optical calibration/erratum lineage matters for this step too.

The release's SNooPy grid header identifies `SNooPy_v2.5.2`; its `.info` file only names the grid. The historical native reader describes the grid as based on Burns2018. These files do not provide a per-object executed training roster. M20 has an existing 79-name training roster, with29 names matching the listed CSP data and six of the eight cadence-eligible objects. These are exact-name overlap proxies, not a proof of independent training for the other names. Even leaving the tested phases out of a new fit does not remove their information from an already trained template. A same-data closure check remains useful, but it is not external model validation.

## Precise next experiment

The immediate step is **an expanded provenance and cadence inventory, keeping the frozen windows and thresholds unchanged**, followed by a model-only identifiability calculation if sufficient data exist. The official CSP data page supplies the full CSP-I DR3 archive, filters and zero points, with134 objects in the release; it is larger than the RAISIN LIST. The parent is acquiring that archive and the2017/2020 errata separately. Apply this same metadata filter to its exact physical passbands, preserve the original eight-object result, and independently establish complete epoch/limit semantics, optical calibration lineage and a per-object training map. Additional objects are not automatically independent or ordinary SNe Ia. [Official CSP data](https://csp.obs.carnegiescience.edu/data), [DR3 release notice](https://csp.obs.carnegiescience.edu/news-items/csp-dr3-photometry-released).

At minimum the current primary cadence needs seven additional supported Y objects, including four additional second-band pairs, without weakening the phase/night thresholds. The full archive may or may not supply them. If it does not, this design remains a small internal engineering check; switching bands or widening phases requires a new explicitly motivated target rather than a retroactive primary.

If the measurement, timing, training and information gates close, freeze object IDs and a single primary signed phase coefficient before observed scores. Retain all predeclared objects, field/night calibration structure, and shared template uncertainty. Use each object's full early/late joint likelihood; do not impose a residual-dependent clipping or re-estimate the template on the validation objects. A leave-object-out retraining control is only meaningful if the original training procedure and input roster are executable.

Only after independently validated phase behavior is measured should it be transported: apply the fixed low-z phase modification to the actual high-z passbands and epoch distribution, retaining full spectral/K-correction and nuisance uncertainty, and rerun the selection-normalized distance inference. Shared low/high-z population and measurement assumptions must be explicit. Do not add the published .164mag split, or any mean contrast, directly to high-z distances. This route can test a **phase-dependent part** of the observational model. It cannot by itself break the remaining distance–gray-luminosity degeneracy or support new cosmology.

## Reproduction

`runs/research_2026_09_26/astra_design/raisin_phase_identification/` contains the frozen protocol, metadata-only parser, all epoch/filter/object ledgers, original-source hashes, published-membership labels, and captured primary-source HTML with SHA256 records in `web-sources.json`. No observed residual file or native-fit output is produced by this task.

## Separately frozen full-DR3 metadata extension

The official `CSP_Photometry_DR3.tgz` archive is now available with SHA256 `ea337b375a7da6223b11f2dd6b50e8a4546045c3133719cc66cbf1bdf4795026`. The original76-object result above is preserved. A separate [extension protocol](../specifications/experiments/raisin_phase_identification/dr3_extension/protocol.json) freezes the historical converter's candidate time rule (`source date+53000`) and filter map. In particular Jrc2→j and Ydw→y are explicit; Jdw/Hdw would remain unmapped. The archive README distinguishes these physical systems. This is a candidate source map, not proof that the historical analysis ran that converter.

The134 SNpy files contain19,376 epoch records. Seventy-three names link to the listed RAISIN CSP sample. A conservative predeclared identity gate also requires coordinates within5arcsec and zHEL within.0002. Two names have identical coordinates but fail that redshift check:2005al (source-minus-release−.00298),2009ab (−.00097). The first assertion failure is preserved; a versioned inventory only changes abort behavior to retain failed objects as unevaluated, without loosening the threshold. Their redshift provenance is not inferred.

For the71 passing objects, all11,115 mapped source epoch records have a released same-filter/time candidate within.00055day:10,529 have one candidate and586 have multiple candidates. Duplicate cardinality is retained and no flux-based assignment is made. No source Jdw/Hdw records occur in these SNpy files. This verifies the candidate time/filter linkage on these objects, **not** signed measurement completeness, error equivalence or independent repeated measurements.

With unchanged windows and provisional linked header peaks, seven exactY objects qualify: the original eight except2012fr, which is not in this archive. Sixty-one archive names have no linked peak; two additional names fail the conservative redshift bridge. The source SNpy headers contain name,z,RA,Dec, not PEAKMJD, so their phase support remains unevaluated. The historical converter's fallback peak is the brightest observation; it is neither an independent optical timing likelihood nor used here. Consequently the full134-object archive does not yet close the cadence feasibility gate. The minimum next metadata input is a provenance-linked optical-only peak estimate/likelihood for these remaining objects, with population and training labels; it must be obtained without inspecting NIR contrast outcomes.

The extension outputs and immutable protocol versions are under `raisin_phase_identification/dr3_extension/`. The parent is acquiring/reading the2017 and2020 photometric errata separately; this report does not claim to have read or applied them. No magnitude, error or flux values were parsed by the extension, and no observed phase score was generated.


## Raw table versus SNpy: instrument labels were combined upstream

A separately frozen [metadata multiset test](../specifications/experiments/raisin_phase_identification/filter_lineage/protocol.json) compares all19,376 `SN_photo.dat` records with all134 SNpy files, using only exact name, decimal time and filter labels. It does not parse magnitudes/errors. The raw table has **520 Jdw and496 Hdw** observations, whereas the SNpy files have none. Mapping Jdw→J and Hdw→H reproduces **all5,491 NIR records exactly, including duplicate multiplicities**. The only remaining global differences are1,768 raw V0 labels appearing as V in SNpy. Those optical labels need their own calendar/physical-filter provenance check before an optical timing fit.

Thus a released SNpy J/H token does not prove that every epoch used the same physical instrument bandpass. The exact metadata collapse does not tell us whether an S-correction was applied, whether magnitudes are unchanged, or whether stellar-equivalent bands are equivalent for an evolving SN SED. The original README/converter mapping cannot establish those properties. No photometric bias or faulty author procedure is inferred.

Using the **raw physical labels**, with the unchanged early/late/night/span criteria and the71 original identity passes, yields Y7, J4, H3, and zero eligible Jdw/Hdw/Jrc2/Ydw pairs. Five of the seven Y objects have an eligible second physical band:2005el,2005kc,2006mr,2008hv,2008ia. The seven Y IDs remain unchanged. This is a provenance refinement, not a new primary band/window or outcome-selected cut. It leaves the local design below the frozen15Y/10second-band thresholds and adds a clear prerequisite: either retain the original physical filters with verified throughput/zero points, or demonstrate the actual synthetic-photometry transformation for SNe across both phase windows before treating combined J/H as one band.
