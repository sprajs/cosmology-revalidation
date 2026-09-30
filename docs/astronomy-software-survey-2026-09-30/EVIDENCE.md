# Evidence, coverage and reproducibility

Collected 30 September 2026 in Europe/London. UTC timestamps in the logs begin on 29 September because the collection crossed the local midnight boundary.

## What was checked

- Curated entries: **528**, in **26** groups; every entry has a nominated source and a retrieval attempt.
- HTTP 200 responses: **526**. This is a transport-level count, not a count of verified scientific implementations.
- Extracted text of at least 450 characters: **459**. This threshold is only a content triage aid; it is not a scientific or identity-verification score.
- Short pages, metadata, redirects, access challenges or JavaScript shells: **67**.
- Unresolved nominated-source retrievals: **2**. Retrieval failure does not imply that software is unavailable elsewhere.
- ASCL breadth index: **4,260** records from **43** browse pages: **4,105 registered** and **155 submitted**. Unique record URLs and registry totals agree.

The review used project repositories, author pages, observatory documentation, software papers and registry records. Searches covered all 26 catalogue groups, including explicit legacy and alternative-language sweeps. Capability and language summaries were curated; repository contents were not cloned or exhaustively audited. Individual packages were not installed, tested or benchmarked. Release dates and activity metrics were not systematically collected, so this is not a maintenance ranking.

## Unresolved direct retrievals

| Entry | Nominated source | Retrieval result |
|---|---|---|
| SSE and BSE | [Source](https://astronomy.swin.edu.au/~jhurley/) | URLError: <urlopen error timed out> |
| Wilson Devinney code | [Source](https://www.astro.ufl.edu/~wilson/) | URLError: <urlopen error [Errno -2] Name or service not known> |

The Wilson–Devinney lineage also has [ASCL record 2004.004](https://ascl.net/2004.004). SSE/BSE are also referenced by coupling frameworks such as [AMUSE](https://github.com/amusecode/amuse); this does not resolve the old author-site retrieval.

## Short, redirected or limited source pages

These entries remain useful, but the automatic source-text snapshot is limited. Several have substantial documentation elsewhere, or were supported by additional web-search results. This table avoids treating a repository title or an HTTP 200 access challenge as full documentation.

| Entry | Extracted characters | Landing-page title |
|---|---:|---|
| [synphot](https://github.com/spacetelescope/synphot_refactor) | 167 | GitHub - spacetelescope/synphot_refactor: Synthetic photometry using Astropy · GitHub |
| [ESO MIDAS](https://www.eso.org/sci/software/esomidas/) | 0 | (no usable title) |
| [IDL Astronomy Library](https://github.com/wlandsman/IDLAstro) | 69 | GitHub - wlandsman/IDLAstro: Astronomy related procedures in the commercial IDL language · GitHub |
| [PDL](https://pdl.perl.org/) | 260 | Perl Data Language |
| [Source Extractor](https://github.com/astromatic/sextractor) | 214 | GitHub - astromatic/sextractor: Extract catalogs of sources from astronomical images · GitHub |
| [PSFEx](https://github.com/astromatic/psfex) | 156 | GitHub - astromatic/psfex: Generate PSF super-tabulated models · GitHub |
| [SCAMP](https://github.com/astromatic/scamp) | 159 | GitHub - astromatic/scamp: Compute astrometric solutions from SExtractor catalogs · GitHub |
| [SWarp](https://github.com/astromatic/swarp) | 191 | GitHub - astromatic/swarp: resample and coadd FITS images to an arbitrary astrometric projection · GitHub |
| [statmorph](https://github.com/vrodgom/statmorph) | 426 | GitHub - vrodgom/statmorph: Python code for calculating non-parametric morphological diagnostics of galaxy images. · GitHub |
| [ESO instrument recipes](https://www.eso.org/sci/software/pipelines/) | 0 | (no usable title) |
| [linetools](https://github.com/linetools/linetools) | 212 | GitHub - linetools/linetools: A package for the analysis of 1d astronomical spectra, especially quasar and galaxy spectra. · GitHub |
| [Fermipy](https://github.com/fermiPy/fermipy) | 259 | GitHub - fermiPy/fermipy: Fermi-LAT Python Analysis Framework · GitHub |
| [CubiCal](https://github.com/ratt-ru/CubiCal) | 236 | GitHub - ratt-ru/CubiCal: A fast radio interferometric calibration suite. · GitHub |
| [SoFiA 2](https://github.com/SoFiA-Admin/SoFiA-2) | 84 | GitHub - SoFiA-Admin/SoFiA-2: Version 2 of the Source Finding Application (SoFiA) · GitHub |
| [TiRiFiC](https://github.com/gigjozsa/tirific) | 7 | GitHub - gigjozsa/tirific: Tilted-Ring-Fitting-Code · GitHub |
| [PSRCHIVE](https://psrchive.sourceforge.net/) | 421 | PSRCHIVE |
| [DSPSR](https://dspsr.sourceforge.net/) | 409 | DSPSR Project |
| [TEMPO2](https://bitbucket.org/psrsoft/tempo2/) | 9 | Bitbucket |
| [gravlens and lensmodel](https://www.physics.rutgers.edu/~keeton/gravlens/) | 173 | Index of /~keeton/gravlens |
| [metadetect](https://github.com/esheldon/metadetect) | 79 | GitHub - esheldon/metadetect: Library for meta-detection, combining detection and metacalibration · GitHub |
| [LensTools](https://github.com/apetri/LensTools) | 434 | GitHub - apetri/LensTools: Useful computing tools for Weak Lensing analyses http://lenstools.readthedocs.io · GitHub |
| [CosmoSIS](https://github.com/joezuntz/cosmosis) | 238 | GitHub - cosmosis-developers/cosmosis · GitHub |
| [Colossus](https://bitbucket.org/bdiemer/colossus/) | 9 | Bitbucket |
| [sotodlib](https://github.com/simonsobs/sotodlib) | 146 | GitHub - simonsobs/sotodlib: Simons Observatory: Time-Ordered Data processing library. · GitHub |
| [SMICA and Planck likelihood code](https://pla.esac.esa.int/) | 157 | Planck Legacy Archive |
| [SNCosmo](https://github.com/sncosmo/sncosmo) | 110 | GitHub - sncosmo/sncosmo: Python library for supernova cosmology · GitHub |
| [ADIPLS](https://users-phys.au.dk/jcd/adipack.n/) | 168 | Obsolete version of adipack |
| [BPASS](https://bpass.auckland.ac.nz/) | 397 | Home \| BPASS Auckland |
| [kcorrect](https://github.com/blanton144/kcorrect) | 113 | GitHub - blanton144/kcorrect · GitHub |
| [ATLAS and SYNTHE](https://wwwuser.oats.inaf.it/castelli/sources.html) | 102 | Source Codes |
| [TLUSTY and SYNSPEC](https://tlusty.oca.eu/) | 16 | Tlusty Home Page |
| [Cloudy](https://gitlab.nublado.org/cloudy/cloudy) | 211 | cloudy / cloudy · GitLab |
| [KROME](https://bitbucket.org/tgrassi/krome/) | 9 | Bitbucket |
| [DESPOTIC](https://bitbucket.org/krumholz/despotic/) | 9 | Bitbucket |
| [GADGET-2](https://wwwmpa.mpa-garching.mpg.de/gadget/) | 36 | Cosmological simulations with GADGET |
| [MPI-AMRVAC](https://github.com/amrvac/amrvac) | 100 | GitHub - amrvac/amrvac: MPI-AMRVAC: A Parallel Adaptive Mesh Refinement Framework · GitHub |
| [SpEC](https://www.black-holes.org/code/SpEC.html) | 65 | Redirecting… |
| [Teukolsky and KerrGeodesics](https://github.com/BlackHolePerturbationToolkit/Teukolsky) | 380 | GitHub - BlackHolePerturbationToolkit/Teukolsky: A Mathematica package for computing solutions to the Teukolsky equation. · GitHub |
| [LALSuite](https://git.ligo.org/lscsoft/lalsuite) | 327 | lscsoft / lalsuite · GitLab |
| [GstLAL](https://git.ligo.org/lscsoft/gstlal) | 202 | lscsoft / GstLAL · GitLab |
| [BayesWave](https://git.ligo.org/lscsoft/bayeswave) | 222 | lscsoft / bayeswave · GitLab |
| [RIFT](https://git.ligo.org/rapidpe-rift/rift) | 221 | RapidPE-RIFT / RIFT · GitLab |
| [PESummary](https://github.com/pesummary/pesummary) | 217 | GitHub - pesummary/pesummary: Code to generate summary pages for Parameter Estimation results. Mirror of https://git.ligo.org/lscsoft/pesummary · GitHub |
| [ligo.skymap](https://git.ligo.org/lscsoft/ligo.skymap) | 262 | lscsoft / ligo.skymap · GitLab |
| [PyRing](https://git.ligo.org/lscsoft/pyring) | 172 | lscsoft / pyRing · GitLab |
| [batman](https://github.com/lkreidberg/batman) | 182 | GitHub - lkreidberg/batman: Fast transit light curves models in Python. · GitHub |
| [starry](https://github.com/rodluger/starry) | 39 | GitHub - rodluger/starry: Tools for mapping stars and planets. · GitHub |
| [petitRADTRANS](https://gitlab.com/mauricemolli/petitRADTRANS) | 197 | Paul Mollière / petitRADTRANS · GitLab |
| [Tudat](https://github.com/tudat-team/tudat) | 143 | GitHub - tudat-team/tudat: A C++ platform to perform astrodynamics and space research. · GitHub |
| [GMAT](https://gmat.atlassian.net/wiki/) | 0 | (no usable title) |
| [LISIRD solar model and data archive](https://lasp.colorado.edu/lisird/) | 53 | LISIRD — LASP Interactive Solar Irradiance Datacenter |
| [JuliaAstro ecosystem](https://juliaastro.org/) | 0 | (no usable title) |
| [Cosmology.jl](https://github.com/JuliaAstro/Cosmology.jl) | 164 | GitHub - JuliaAstro/Cosmology.jl: Cosmology library for Julia · GitHub |
| [VLBISkyModels.jl](https://github.com/EHTJulia/VLBISkyModels.jl) | 337 | GitHub - EHTJulia/VLBISkyModels.jl: Simple models for on-sky emission of radio data · GitHub |
| [Gradus.jl](https://github.com/astro-group-bristol/Gradus.jl) | 18 | GitHub - astro-group-bristol/Gradus.jl: (Moved to Codeberg) Extensible spacetime agnostic general relativistic ray-tracing (GRRT). · GitHub |
| [GROMACS](https://gitlab.com/gromacs/gromacs) | 251 | GROMACS / GROMACS · GitLab |
| [SAOImage DS9](https://ds9.si.edu/site/Home.html) | 24 | (no usable title) |
| [nalgebra](https://github.com/dimforge/nalgebra) | 86 | GitHub - dimforge/nalgebra: Linear algebra library for Rust. · GitHub |
| [faer](https://github.com/sarah-ek/faer-rs) | 352 | GitHub - sarah-quinones/faer-rs: Linear algebra foundation for the Rust programming language · GitHub |
| [COMPASS](https://github.com/ANR-COMPASS/shesha) | 151 | GitHub - ANR-COMPASS/shesha: ATTENTION: This repository is no longer up to date. The project is now available on GitLab https://gitlab.obspm.fr/cosmic-rtc/compass · GitHub |
| [pyKLIP](https://bitbucket.org/pyKLIP/pyklip/) | 9 | Bitbucket |
| [LITpro](https://www.jmmc.fr/english/tools/proposal-preparation/litpro/) | 270 | Making sure you're not a bot! |
| [OIFITS and OIFITSlib](https://www.jmmc.fr/oifits/) | 270 | Making sure you're not a bot! |
| [VARTOOLS](https://www.astro.princeton.edu/~jhartman/vartools.html) | 49 | VARTOOLS |
| [ROCKSTAR](https://bitbucket.org/gfcstanford/rockstar/) | 9 | Bitbucket |
| [Dask](https://github.com/dask/dask) | 143 | GitHub - dask/dask: Parallel computing with task scheduling · GitHub |
| [SHTns](https://bitbucket.org/nschaeff/shtns/) | 9 | Bitbucket |

## Registry fallbacks

Where an original site was inaccessible, a stable ASCL record was used for identity and scope. Both URLs are retained. A registry record does not guarantee that source code or a working build remains available.

| Entry | Registry source | Original nominated URL |
|---|---|---|
| STSDAS and TABLES | [ASCL](https://ascl.net/1206.003) | [Original](https://www.stsci.edu/institute/software_hardware/stsdas) |
| DAOPHOT and ALLSTAR | [ASCL](https://ascl.net/1104.011) | [Original](https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/community/STETSON/daophot/) |
| DOLPHOT | [ASCL](https://ascl.net/1608.013) | [Original](https://americano.dolphinsim.com/dolphot/) |
| Montage | [ASCL](https://ascl.net/1010.036) | [Original](https://montage.ipac.caltech.edu/) |
| SpecPro | [ASCL](https://ascl.net/1404.014) | [Original](https://specpro.caltech.edu/) |
| L.A.Cosmic original | [ASCL](https://ascl.net/1207.005) | [Original](https://www.astro.yale.edu/dokkum/lacosmic/) |
| ParselTongue | [ASCL](https://ascl.net/1208.020) | [Original](https://www.jive.eu/jivewiki/doku.php?id=parseltongue:parseltongue) |
| Lenstool | [ASCL](https://ascl.net/1102.004) | [Original](https://projets.lam.fr/projects/lenstool/wiki) |
| binary_c | [ASCL](https://ascl.net/2307.035) | [Original](https://binary_c.gitlab.io/) |
| PHOENIX | [ASCL](https://ascl.net/1010.056) | [Original](https://www.hs.uni-hamburg.de/EN/For/ThA/phoenix/) |
| RADEX | [ASCL](https://ascl.net/1010.075) | [Original](https://personal.sron.nl/~vdtak/radex/index.shtml) |
| TORUS | [ASCL](https://ascl.net/1404.006) | [Original](https://www.astro.ex.ac.uk/people/th2/torus_html/) |
| GIZMO | [ASCL](https://ascl.net/1410.003) | [Original](https://www.tapir.caltech.edu/~phopkins/Site/GIZMO.html) |
| GALPROP | [ASCL](https://ascl.net/1010.028) | [Original](https://galprop.stanford.edu/) |
| MPFIT | [ASCL](https://ascl.net/1208.019) | [Original](https://pages.physics.wisc.edu/~craigm/idl/idl.html) |
| AHF | [ASCL](https://ascl.net/1102.009) | [Original](https://popia.ft.uam.es/AHF/) |

## Identity and lineage findings

- GLEE and WSLAP+ are documented using author papers. Their implementation language and public source-distribution status are explicitly unconfirmed here.
- The proprietary-language identification is an inference: IDL is the strongest match. Mark Sullivan’s personal use was not established by a sufficiently direct primary source.
- IRAF community maintenance, ESO MIDAS releases, the IDL-library archive migration and AIPS maintenance show why age alone cannot determine lifecycle status.
- GRChombo/GRTeclyn, Tudat/tudatpy, COMPASS, SoFiA 2, Gradus and faer have successor, archive or migration information in their landing pages. Those relationships appear in the catalogue notes.
- Name collisions require domain identity: ISIS optical subtraction, ISIS X-ray spectroscopy and USGS ISIS are different systems; CLASS cosmology and GILDAS CLASS are different; pulsar PINT is not the Python units package.
- Wrappers, forks, suites and shared engines must not be counted as independent implementations. In particular, Python-facing software often depends on native C/C++/Fortran kernels.

## Files and reproducibility

- `catalogue.json`: curated entries, source links, language labels, capability notes and source-retrieval metadata.
- `source-checks.json`: retrieval attempts, final URLs, UTC timestamps, response metadata and extracted-text checksums.
- `source-extracts/`: local research snapshots of extracted source-page text, for traceability; these are not software source distributions or an independently publishable source collection.
- `ascl-index.json`: registry titles and record links. Abstracts and source code are not reproduced in this index.
- `ascl-index-provenance.json`: browse-page URLs, timestamps, page checksums, counts and collection scope.
- `research-searches.json`: retained thematic search results. Search snippets are leads; primary project sources are the main catalogue references.
- `REPORT.md`, `CATALOGUE.md`, `OVERLAP-MAP.md`, `overlap-map.json`, `EVIDENCE.md`, `catalogue.html`: reading and exploration views.

Rebuild the curated catalogue and presentation from this directory:

```sh
python build_catalogue.py
python verify_sources.py
python build_overlap_map.py
python render_catalogue.py
```

`verify_sources.py` reuses cached checks when the nominated URL is unchanged; use `--retry-unresolved` to retry failures. It does not refresh every cached page. `collect_ascl_index.py` performs a fresh network collection and can change the index. The rendered catalogue intentionally validates the snapshot counts; update its assertions deliberately when making a later survey edition.

## Limits of the sweep

The ASCL index was validated for pagination completeness, record-URL uniqueness and registry counts, not for the correctness of every title or scientific claim. The ASCL and curated lists overlap and their totals must not be added. Submitted registry entries are explicitly labelled. Registry search operates on titles and names, not full abstracts.

Coverage is broad but cannot be a census of private collaboration code, unpublished scripts, every instrument recipe or every proprietary extension. General physics is represented by major frameworks and capability families rather than every computational-physics package. License compatibility, source availability, build reproduction, algorithm ancestry and scientific validation remain per-component follow-up work. The Rust module map is an architectural proposal, not an upstream consensus or an implemented product.
