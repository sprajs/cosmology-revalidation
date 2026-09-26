# Direct cross-instrument calibration: metadata feasibility

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../README.md) for current execution status.

A same-supernova, nearly simultaneous comparison can cancel distance and most
temporal evolution before imposing a fitted cosmology. It still measures a
combination of passband response, source spectrum and instrument calibration;
different natural passbands cannot simply be assigned equal supernova
magnitudes. Shared standards and duplicated measurements also prevent treating
every pair as independent evidence.

The [metadata-only protocol](../specifications/experiments/csp_direct_instrument_overlap/protocol.json)
and [complete pair inventory](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/csp_direct_instrument_overlap/result.json)
inspect object, physical filter and date from the original CSP DR3
`SN_photo.dat`, without reading magnitude or uncertainty columns into the
calculation. Pairing is every Cartesian same-object pair within the declared
time gap. These counts are feasibility measures, not independent sample sizes.

| Physical pair | Pairs / objects within 0.5 day | Pairs / objects within 1 day |
|---|---:|---:|
| WIRC J / RetroCam RC1 J | 4 / 1 | 28 / 7 |
| WIRC J / RetroCam RC2 J | 0 / 0 | 0 / 0 |
| WIRC H / RetroCam H | 0 / 0 | 14 / 4 |
| WIRC Y / RetroCam Y | 4 / 1 | 39 / 8 |

The 0.05- and 0.25-day windows have the same counts as 0.5 day. Every close pair
belongs to SN 2006kf, at the two repeated dates 1037.81 and 1037.82 in the
archive's date convention. This is the object already retained as ambiguous
in the [physical-label lineage](csp-passband-lineage.md). Metadata alone does
not establish that these are distinct exposures or independent instrument
measurements. No magnitude difference, calibration offset or likelihood score
has been calculated from this inventory.

Before a direct comparison, the close-pair measurement lineage must establish
independent exposures. Longer-gap comparisons require a declared interpolation
model and its uncertainty; those extra assumptions cannot be hidden inside a
claim of direct calibration. The current source inventory therefore provides
no immediate model-independent WIRC-versus-RC2 calibration constraint.
