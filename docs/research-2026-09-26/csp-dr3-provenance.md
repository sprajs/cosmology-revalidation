# Original CSP DR3 products and calibration errata

The [official CSP data page](https://csp.obs.carnegiescience.edu/data)
links the original DR3 archive. The acquired 289,496-byte archive has SHA-256
`ea337b375a7da6223b11f2dd6b50e8a4546045c3133719cc66cbf1bdf4795026`.
[Acquisition](../../runs/research_2026_09_26/csp_dr3_provenance/acquisition.json)
and [140 archive-member records](../../runs/research_2026_09_26/csp_dr3_provenance/archive-members.json)
preserve the source URL and individual uncompressed hashes. No archive code
was executed. The files are published photometry, not a newly acquired raw
exposure or nondetection ledger.

The archive README distinguishes Swope/RetroCam J before 2009 January 15
from Jrc2 afterward, and distinguishes du Pont/WIRC Ydw/Jdw/Hdw. These
names must not be pooled merely because their broad wavelength names agree.
The tagged [SNANA converter](../../sources/repos/RickKessler__SNANA@v11_04k/util/translate_CSPDR3.py:73)
maps J→J, Jrc2→j, Y→Y, Ydw→y and H→H, and adds 53000 to the published
times. This source is an implementation witness, not proof of the exact
historical conversion run. The subsequent [magnitude and passband lineage](csp-passband-lineage.md)
now verifies all 5,491 NIR archive rows and 73 matching RAISIN objects, and
identifies the WIRC J to RC1 assignment for a paired distance test.

A separate [foreground-scale check](../../runs/research_2026_09_26/csp_dr3_provenance/mw-scale-check.json)
finds the `SFD98 * 0.86` convention labelled in all 78 CSP DAT files,
including the duplicate calibrator file. All nine nominal release NMLs use
`OPT_MWEBV=1`, scale1 and shift0. Tagged `MWgaldust.c:256` leaves the
input E(B−V) and error unchanged for that option. Thus this source/configuration
path does not apply the 0.86 factor twice. The 0.95 multiplier is a separate
systematic variant. This validates the bookkeeping, not the physical map.

The converter also maps magnitude uncertainty to a one-sided positive flux
increment, `sigma_f = f [10^(0.4 sigma_m) − 1]`, rather than to the linearized
derivative. That convention must be reproduced in a baseline; it does not
by itself establish a symmetric Gaussian measurement distribution. No new
error rescaling has been applied here.

The two published errata are independently listed by the
[Caltech author repository](https://authors.library.caltech.edu/records/0f6hh-4pr18).
Direct PDF acquisition from Caltech failed with HTTP403; a second institution
returned an access-challenge page. These failures are retained, and no local
PDF is claimed. The accessible text of the original 2020 paper identifies
the BV problem as the conversion of field-star magnitudes back to the
standard system for publication. It explicitly states that the natural-system
supernova photometry was not affected. Thus this erratum cannot justify
applying the field-star change directly to the SN fluxes. Independent
cross-calibration must use the appropriate corrected stellar catalogue and
conversion. The same erratum corrects the published i colour coefficient
and one coordinate. [Krisciunas et al. 2020, AJ160,289, original-paper text](https://www.researchgate.net/publication/347185134_Erratum_The_Carnegie_Supernova_Project_I_Third_Photometry_Data_Release_of_Low-redshift_Type_Ia_Supernovae_and_Other_White_Dwarf_Explosions_2017_AJ_154_211)

The phase question also has prior literature. RAISIN's authors reported a
0.164±0.051mag difference between objects observed earlier and later in
the NIR, discussed stretch/selection confounding, and noted shared training
between their model alternatives. That published between-object statistic
is not a new result here. A same-object, same-filter early/late comparison
could remove distance and constant calibration offsets, but would still
need independent timing, adequate phase coverage and measurement-selection
provenance. [Jones et al. 2022, section V.3](https://arxiv.org/html/2201.07801v2#S5.SS3)
