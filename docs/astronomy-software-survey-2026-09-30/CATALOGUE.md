# Curated astronomy and physics software catalogue

528 entries across 26 groups. Compiled 30 September 2026, Europe/London.

This is a capability inventory, including legacy and unmaintained software. An entry can be a package, a family, a suite, a standard or a data service. Entries are not independent numerical engines. Language labels identify principal implementation or interface families; they are not source-percentage audits. A retrieved source is not evidence of present maintenance, buildability, adoption or numerical correctness.

Read [the analysis](REPORT.md), use [the searchable catalogue](catalogue.html), or inspect [retrieval evidence and limitations](EVIDENCE.md). The separate ASCL index has 4,260 title-and-link records, including 155 submitted records. The two collections overlap.

| Domain | Entries |
|---|---:|
| 01 Foundations and interoperability | 24 |
| 02 Integrated suites and legacy environments | 13 |
| 03 Image reduction and measurement | 20 |
| 04 Spectroscopy and observatory pipelines | 24 |
| 05 X ray and gamma ray astronomy | 24 |
| 06 Radio interferometry and spectral cubes | 30 |
| 07 Pulsars fast transients and timing | 11 |
| 08 Strong weak and microlensing | 22 |
| 09 Cosmology CMB and large scale structure | 26 |
| 10 Supernovae and transient ecosystems | 18 |
| 11 Stars populations and spectral synthesis | 25 |
| 12 Radiation chemistry and atomic physics | 17 |
| 13 Dynamics hydrodynamics and cosmological simulations | 28 |
| 14 Relativity spacetime and black holes | 23 |
| 15 Gravitational waves and multi messenger inference | 14 |
| 16 Exoplanets planetary and solar system science | 19 |
| 17 Solar plasma and space physics | 14 |
| 18 Julia scientific and astronomy ecosystems | 18 |
| 19 General numerical inference and physics foundations | 25 |
| 20 General physics quantum and particle software | 17 |
| 21 Visualization archives and operations | 14 |
| 22 Existing Rust scientific building blocks | 14 |
| 23 Optical instruments adaptive optics and interferometry | 14 |
| 24 Additional legacy R and domain specialist software | 15 |
| 25 Additional simulation and emerging modelling capabilities | 24 |
| 26 Shared data storage execution and numerical infrastructure | 35 |

## 01 Foundations and interoperability

| Software / source | Language family | Capability | Lineage, integration and evidence notes |
|---|---|---|---|
| [Astropy](https://github.com/astropy/astropy) · ASTRO-0001 | Python; C extensions | Units, quantities, constants, time scales, coordinates, tables, FITS, WCS, cosmology, statistics and models | Core capability reference; separate from its affiliated ecosystem. Source text retrieved. |
| [Astroquery](https://github.com/astropy/astroquery) · ASTRO-0002 | Python | Query astronomical archives and catalogue services | Service adapters require their own schema and authentication handling. Source text retrieved. |
| [PyVO](https://github.com/astropy/pyvo) · ASTRO-0003 | Python | Virtual Observatory discovery and access protocols | TAP, SIA, SSA and registry interoperability. Source text retrieved. |
| [ASDF](https://github.com/asdf-format/asdf) · ASTRO-0004 | Python; format specification | Serialize scientific trees, arrays and metadata | A format and extension ecosystem as well as a Python implementation. Source text retrieved. |
| [GWCS](https://github.com/spacetelescope/gwcs) · ASTRO-0005 | Python | Composable generalized world coordinate transformations | Separate from FITS WCS conventions. Source text retrieved. |
| [WCSLIB](https://www.atnf.csiro.au/people/mcalabre/WCS/) · ASTRO-0006 | C | FITS world coordinate transformations and projections | Foundational native library used through many wrappers. Source text retrieved. |
| [CFITSIO](https://heasarc.gsfc.nasa.gov/fitsio/) · ASTRO-0007 | C; Fortran interfaces | Read and write FITS images and binary tables | Format implementation reference. Source text retrieved. |
| [ERFA and PyERFA](https://github.com/liberfa/erfa) · ASTRO-0008 | C; Python bindings | Astronomical positional algorithms derived from SOFA | Track SOFA lineage; PyERFA is a wrapper, not an independent model. Source text retrieved. |
| [IAU SOFA](https://www.iausofa.org/) · ASTRO-0009 | C; Fortran | Reference algorithms for time, Earth rotation and celestial coordinates | Standards and validation reference. Source text retrieved. |
| [Starlink AST](https://starlink.eao.hawaii.edu/starlink/AST) · ASTRO-0010 | C; Fortran interfaces | Coordinate frames, mappings and world coordinate systems | A substantial alternative coordinate architecture. Source text retrieved. |
| [Starlink HDS and NDF](https://starlink.eao.hawaii.edu/starlink) · ASTRO-0011 | C; Fortran | Hierarchical and N-dimensional scientific data structures | Retain variance, quality and history semantics in format conversion. Source text retrieved. |
| [HEALPix and healpy](https://healpix.sourceforge.io/) · ASTRO-0012 | Fortran; C++; Python | Equal-area spherical pixelization, maps and spherical harmonics | An algorithm and multi-language ecosystem; not just healpy. Source text retrieved. |
| [astropy-healpix](https://github.com/astropy/astropy-healpix) · ASTRO-0013 | Python; C | HEALPix indexing and coordinate operations | Different implementation scope from the full HEALPix suite. Source text retrieved. |
| [MOCPy and CDS MOC](https://github.com/cds-astro/mocpy) · ASTRO-0014 | Python; Rust | Multi-order sky coverage and set operations | Existing Rust core is relevant to reuse. Source text retrieved. |
| [Regions](https://github.com/astropy/regions) · ASTRO-0015 | Python | Pixel and sky regions with astronomy file interfaces | A shared selection primitive. Source text retrieved. |
| [Reproject](https://github.com/astropy/reproject) · ASTRO-0016 | Python | Reprojection and mosaicking of astronomical images | Flux conservation and surface brightness semantics matter. Source text retrieved. |
| [ndcube](https://github.com/sunpy/ndcube) · ASTRO-0017 | Python | Multidimensional data with coordinates, masks and uncertainties | Useful shared cube object model. Source text retrieved. |
| [specutils](https://github.com/astropy/specutils) · ASTRO-0018 | Python | Spectral data containers, manipulation and analysis | Shared spectral abstraction. Source text retrieved. |
| [synphot](https://github.com/spacetelescope/synphot_refactor) · ASTRO-0019 | Python | Synthetic photometry and passband integration | Calibration data are separate dependencies. Short page, metadata or access shell. |
| [dust_extinction](https://github.com/karllark/dust_extinction) · ASTRO-0020 | Python | Extinction laws for astronomical sources | Physical assumptions and valid wavelength domains differ by law. Source text retrieved. |
| [dustmaps](https://github.com/gregreen/dustmaps) · ASTRO-0021 | Python | Interfaces to Galactic dust maps | Code does not include every map dataset. Source text retrieved. |
| [Skyfield](https://github.com/skyfielders/python-skyfield) · ASTRO-0022 | Python | Ephemerides, positions and observation geometry | Depends on ephemeris and time data. Source text retrieved. |
| [SPICE Toolkit](https://naif.jpl.nasa.gov/naif/toolkit.html) · ASTRO-0023 | Fortran; C; IDL and MATLAB interfaces | Spacecraft and planetary geometry, time and reference frames | Kernels are essential data products; SpiceyPy is a Python binding. Source text retrieved. |
| [IVOA standards](https://www.ivoa.net/documents/) · ASTRO-0024 | Specifications; multiple implementations | VOTable, TAP, ADQL, SAMP, UWS, SIAv2, SSA, ObsCore, MOC and HiPS | Standards inventory; not counted as an executable package. Source text retrieved. |

## 02 Integrated suites and legacy environments

| Software / source | Language family | Capability | Lineage, integration and evidence notes |
|---|---|---|---|
| [AMUSE](https://github.com/amusecode/amuse) · ASTRO-0025 | Python with C, C++ and Fortran community codes | Couple gravitational dynamics, stellar evolution, hydrodynamics and radiative transfer | Closest architectural analogue for multiphysics integration; its code list warns it can be stale. Source text retrieved. |
| [Starlink Software Collection](https://github.com/Starlink/starlink) · ASTRO-0026 | Fortran; C; C++; Java; Perl; Tcl/Tk | Image and spectral reduction, coordinates, data models and pipeline tools | Longstanding suite with current community and observatory stewardship. Source text retrieved. |
| [IRAF Community Distribution](https://iraf-community.github.io/) · ASTRO-0027 | SPP; C; Fortran; CL | General astronomical image and spectral reduction | Institutional development ended; community maintenance continues. Do not label simply dead. Source text retrieved. |
| [PyRAF](https://github.com/iraf-community/pyraf) · ASTRO-0028 | Python | Python access to IRAF tasks and command workflows | Community continuation of a legacy interface. Source text retrieved. |
| [STSDAS and TABLES](https://ascl.net/1206.003) · ASTRO-0029 | SPP; CL; C; Python | Historical HST calibration and analysis tasks under IRAF | Legacy lineage; modern HST and Python successors must be considered per task. Source text retrieved. [Original nominated site](https://www.stsci.edu/institute/software_hardware/stsdas); source link uses ASCL. |
| [ESO MIDAS](https://www.eso.org/sci/software/esomidas/) · ASTRO-0030 | Fortran; C; MIDAS command language | Munich Image Data Analysis System for images and spectra | ESO lists it as legacy software while still advertising a 25MAY release; legacy does not mean abandoned. Short page, metadata or access shell. |
| [IDL](https://www.nv5geospatialsoftware.com/products/IDL) · ASTRO-0031 | Proprietary IDL runtime | Array-oriented numerical analysis and visualization environment | Likely the language recalled by the user; Sullivan attribution not established. Source text retrieved. |
| [IDL Astronomy Library](https://github.com/wlandsman/IDLAstro) · ASTRO-0032 | IDL | FITS I/O, astrometry, photometry and numerical utilities | NASA website frozen in 2022 points to this repository; runtime and library licenses differ. Short page, metadata or access shell. |
| [GNU Data Language](https://github.com/gnudatalanguage/gdl) · ASTRO-0033 | C++; IDL-compatible language | Open implementation of much of the IDL language and numerical environment | Compatibility reference; do not assume complete IDL equivalence. Source text retrieved. |
| [SolarSoft](https://www.lmsal.com/solarsoft/) · ASTRO-0034 | Primarily IDL; other languages | Solar mission analysis libraries, calibration and common utilities | Large distributed ecosystem; not one homogeneous package. Source text retrieved. |
| [PDL](https://pdl.perl.org/) · ASTRO-0035 | Perl; C | N-dimensional numerical arrays and scientific computing in Perl | Relevant numerical foundation behind Perl scientific workflows. Short page, metadata or access shell. |
| [ORAC-DR and PICARD](https://github.com/Starlink/ORAC-DR) · ASTRO-0036 | Perl; external Starlink tasks | Recipe-based instrument reduction and post-processing | Current and historical pipeline knowledge; include calibration selection rules. Source text retrieved. |
| [Gnuastro](https://www.gnu.org/software/gnuastro/) · ASTRO-0037 | C; command line | Astronomical arithmetic, detection, measurement and image processing | Integrated command-line suite with reusable libraries. Source text retrieved. |

## 03 Image reduction and measurement

| Software / source | Language family | Capability | Lineage, integration and evidence notes |
|---|---|---|---|
| [Photutils](https://github.com/astropy/photutils) · ASTRO-0038 | Python | Source detection, segmentation, aperture and PSF photometry, backgrounds | Major Astropy ecosystem capability. Source text retrieved. |
| [ccdproc](https://github.com/astropy/ccdproc) · ASTRO-0039 | Python | CCD calibration, combination, uncertainty and mask processing | Detector correction primitives. Source text retrieved. |
| [Astro-SCRAPPY](https://github.com/astropy/astroscrappy) · ASTRO-0040 | Cython; Python | Cosmic-ray rejection based on L.A.Cosmic | Preserve noise and saturation assumptions. Source text retrieved. |
| [Source Extractor](https://github.com/astromatic/sextractor) · ASTRO-0041 | C | Detect sources, deblend and measure image catalogues | Also known as SExtractor. Short page, metadata or access shell. |
| [SEP](https://github.com/sep-developers/sep) · ASTRO-0042 | C; Python | Source extraction library derived from Source Extractor | Related implementation lineage; do not count algorithm twice. Source text retrieved. |
| [PSFEx](https://github.com/astromatic/psfex) · ASTRO-0043 | C | Model spatially varying point-spread functions | Pairs with Source Extractor catalogues. Short page, metadata or access shell. |
| [SCAMP](https://github.com/astromatic/scamp) · ASTRO-0044 | C | Astrometric and photometric calibration of imaging surveys | Global catalogue and exposure calibration. Short page, metadata or access shell. |
| [SWarp](https://github.com/astromatic/swarp) · ASTRO-0045 | C | Resample and coadd astronomical images | Separate coaddition covariance from per-pixel variance. Short page, metadata or access shell. |
| [Astrometry.net](https://github.com/dstndstn/astrometry.net) · ASTRO-0046 | C; Python | Blind astrometric image solving | Index catalogues are separate resources. Source text retrieved. |
| [DAOPHOT and ALLSTAR](https://ascl.net/1104.011) · ASTRO-0047 | Fortran | Crowded-field stellar PSF photometry | Historical and specialist capability; distribution and rights need checking. Source text retrieved. [Original nominated site](https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/community/STETSON/daophot/); source link uses ASCL. |
| [DOLPHOT](https://ascl.net/1608.013) · ASTRO-0048 | C | Stellar photometry with instrument-specific PSFs and calibrations | Important HST and other instrument lineage. Source text retrieved. [Original nominated site](https://americano.dolphinsim.com/dolphot/); source link uses ASCL. |
| [DoPHOT](https://github.com/AbhijitSaha/DoPhot) · ASTRO-0049 | Fortran | Automated PSF stellar photometry | Legacy algorithm family with later implementations. Source text retrieved. |
| [GALFIT](https://users.obs.carnegiescience.edu/peng/work/galfit/galfit.html) · ASTRO-0050 | C | Parametric two-dimensional galaxy image decomposition | Executable availability is not the same as reusable source availability. Source text retrieved. |
| [Imfit](https://github.com/perwin/imfit) · ASTRO-0051 | C++ | Galaxy image fitting and structural models | Alternative to GALFIT with public implementation. Source text retrieved. |
| [The Tractor](https://github.com/dstndstn/tractor) · ASTRO-0052 | Python; native extensions | Forward modelling and forced photometry across images | Models sources in the individual observation space. Source text retrieved. |
| [scarlet](https://github.com/pmelchior/scarlet) · ASTRO-0053 | Python | Constrained multiband source separation and deblending | Research deblending model family. Source text retrieved. |
| [statmorph](https://github.com/vrodgom/statmorph) · ASTRO-0054 | Python | Nonparametric galaxy morphology measurements | Gini, asymmetry and other morphology statistics. Short page, metadata or access shell. |
| [HOTPANTS](https://github.com/acbecker/hotpants) · ASTRO-0055 | C | PSF-matched image subtraction | Alard-Lupton difference imaging lineage. Source text retrieved. |
| [ISIS image subtraction](https://www.iap.fr/useriap/alard/package.html) · ASTRO-0056 | C; scripts | Image registration and difference photometry | Distinct from the X-ray package named ISIS. Source text retrieved. |
| [Montage](https://ascl.net/1010.036) · ASTRO-0057 | C; Python interfaces | Astronomical image mosaics and background matching | Mosaic workflow reference. Source text retrieved. [Original nominated site](https://montage.ipac.caltech.edu/); source link uses ASCL. |

## 04 Spectroscopy and observatory pipelines

| Software / source | Language family | Capability | Lineage, integration and evidence notes |
|---|---|---|---|
| [PypeIt](https://github.com/pypeit/PypeIt) · ASTRO-0058 | Python | Optical and near-infrared spectroscopic reduction | Broad multi-instrument pipeline. Source text retrieved. |
| [Gemini DRAGONS](https://github.com/GeminiDRSoftware/DRAGONS) · ASTRO-0059 | Python | Recipe-driven Gemini data reduction | Includes instrument data models and calibration workflows. Source text retrieved. |
| [AstroData](https://github.com/GeminiDRSoftware/astrodata) · ASTRO-0060 | Python | Common interface for astronomical instrument data products | Separate reusable component of Gemini ecosystem. Source text retrieved. |
| [JWST pipeline](https://github.com/spacetelescope/jwst) · ASTRO-0061 | Python; compiled dependencies | Detector correction, calibration, imaging and spectroscopy | Requires instrument reference files and CRDS context. Source text retrieved. |
| [HSTCAL](https://github.com/spacetelescope/hstcal) · ASTRO-0062 | C | Hubble instrument calibration tasks | Includes ACS, WFC3 and STIS calibration lineage. Source text retrieved. |
| [DrizzlePac](https://github.com/spacetelescope/drizzlepac) · ASTRO-0063 | Python; C | HST image alignment, distortion correction and drizzle combination | Sampling and correlated-noise semantics matter. Source text retrieved. |
| [romancal](https://github.com/spacetelescope/romancal) · ASTRO-0064 | Python | Roman WFI imaging calibration and processing | Roman pipeline development; distinguish from completed survey use. Source text retrieved. |
| [stpipe and stcal](https://github.com/spacetelescope/stcal) · ASTRO-0065 | Python | Shared pipeline framework and detector calibration algorithms | Component relationship to JWST and Roman. Source text retrieved. |
| [CRDS](https://github.com/spacetelescope/crds) · ASTRO-0066 | Python | Select and manage calibration reference data | Scientific reproducibility requires pinning reference context. Source text retrieved. |
| [Rubin Science Pipelines](https://pipelines.lsst.io/) · ASTRO-0067 | Python; C++ | Instrument signature removal, calibration, coadds, detection and measurements | Includes task graph and Butler data management architecture. Source text retrieved. |
| [Rubin Butler](https://github.com/lsst/daf_butler) · ASTRO-0068 | Python | Dataset registry, storage and provenance for processing pipelines | Reusable data-management design reference. Source text retrieved. |
| [ESO CPL and EsoRex](https://www.eso.org/sci/software/cpl/) · ASTRO-0069 | C | Common pipeline library and recipe execution | Instrument recipes are additional packages. Source text retrieved. |
| [ESO Reflex](https://www.eso.org/sci/software/esoreflex/) · ASTRO-0070 | Java; Python; external recipes | Graphical reduction workflows for ESO instruments | Workflow and recipe orchestration layer. Source text retrieved. |
| [ESO EDPS and PyCPL](https://www.eso.org/sci/software/pipe_aem_main.html) · ASTRO-0071 | Python; C bindings | Data organisation and pipeline scheduling with Python recipe access | Newer ESO orchestration ecosystem. Source text retrieved. |
| [ESO instrument recipes](https://www.eso.org/sci/software/pipelines/) · ASTRO-0072 | Mostly C via CPL; varies | MUSE, X-shooter, UVES, ESPRESSO, CRIRES, FORS, KMOS, SPHERE and VLTI reduction | Family entry; each instrument needs its own future inventory. Short page, metadata or access shell. |
| [DESI spectroscopic pipeline](https://github.com/desihub/desispec) · ASTRO-0073 | Python; native dependencies | Spectral extraction, calibration, redshifts and survey processing | Includes relationships to specter and redrock. Source text retrieved. |
| [SDSS idlspec2d](https://github.com/sdss/idlspec2d) · ASTRO-0074 | IDL; supporting scripts | SDSS and BOSS optical spectroscopy reduction | Substantial IDL survey pipeline lineage. Source text retrieved. |
| [MaNGA DAP](https://github.com/sdss/mangadap) · ASTRO-0075 | Python | Integral-field spectral analysis and derived galaxy maps | Separate data analysis from upstream reduction. Source text retrieved. |
| [MPDAF](https://github.com/musevlt/mpdaf) · ASTRO-0076 | Python | MUSE integral-field data analysis | Cube operations, spectra and images. Source text retrieved. |
| [pPXF](https://users.physics.ox.ac.uk/~cappellari/software/) · ASTRO-0077 | Python; historical IDL | Stellar and gas kinematics and population fitting from spectra | Retain template resolution, regularization and noise conventions. Source text retrieved. |
| [GANDALF](https://ascl.net/1708.012) · ASTRO-0078 | IDL; Python ports | Joint stellar continuum and emission-line analysis | Legacy lineage; verify chosen port separately. Source text retrieved. |
| [linetools](https://github.com/linetools/linetools) · ASTRO-0079 | Python | Spectroscopic analysis of absorption lines and systems | Atomic line and redshift utilities. Short page, metadata or access shell. |
| [SpecPro](https://ascl.net/1404.014) · ASTRO-0080 | IDL | Interactive spectrum inspection and redshift identification | Proprietary runtime ecosystem. Source text retrieved. [Original nominated site](https://specpro.caltech.edu/); source link uses ASCL. |
| [L.A.Cosmic original](https://ascl.net/1207.005) · ASTRO-0081 | IDL | Laplacian cosmic-ray identification | Legacy algorithm reference alongside Astro-SCRAPPY. Source text retrieved. [Original nominated site](https://www.astro.yale.edu/dokkum/lacosmic/); source link uses ASCL. |

## 05 X ray and gamma ray astronomy

| Software / source | Language family | Capability | Lineage, integration and evidence notes |
|---|---|---|---|
| [HEASoft and FTOOLS](https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/) · ASTRO-0082 | C; C++; Fortran; Perl; Tcl; Python | High-energy mission reduction and general FITS analysis | Large suite; component entries below are not independent suites. Source text retrieved. |
| [XSPEC and PyXspec](https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/) · ASTRO-0083 | C++; Fortran models; Tcl; Python | Forward-folded spectral fitting and physical emission models | Responses, backgrounds and model tables are essential. Source text retrieved. |
| [XSELECT](https://heasarc.gsfc.nasa.gov/docs/software/ftools/xselect/) · ASTRO-0084 | Fortran; command interface | Filter event lists and extract images, spectra and light curves | HEASoft component. Source text retrieved. |
| [XRONOS](https://heasarc.gsfc.nasa.gov/docs/xanadu/xronos/) · ASTRO-0085 | Fortran; C | High-energy time-series analysis | Legacy and maintained-suite component lineage. Source text retrieved. |
| [XIMAGE](https://heasarc.gsfc.nasa.gov/docs/xanadu/ximage/) · ASTRO-0086 | Fortran; C; Tcl | High-energy image analysis and source detection | HEASoft component. Source text retrieved. |
| [HEASoft mission tasks](https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/versions.html) · ASTRO-0087 | Mixed compiled languages; Perl; Python | NICERDAS, NuSTARDAS, Swift, Suzaku, IXPE and XRISM instrument processing | Family entry with explicit component/version evidence. Source text retrieved. |
| [CIAO](https://cxc.cfa.harvard.edu/ciao/) · ASTRO-0088 | C; C++; Python; other compiled components | Chandra calibration, event processing, detection and analysis | Includes mission-specific CALDB dependencies. Source text retrieved. |
| [Sherpa](https://github.com/sherpa/sherpa) · ASTRO-0089 | Python; C++; Fortran | Statistical modelling and fitting of spectra and images | Mission-independent library also distributed with CIAO. Source text retrieved. |
| [XMM Newton SAS](https://www.cosmos.esa.int/web/xmm-newton/sas) · ASTRO-0090 | C++; Fortran; Perl; Python interfaces | XMM event calibration and science products | Instrument calibrations and observation metadata are separate. Source text retrieved. |
| [SPEX](https://www.sron.nl/en/pillars/science/astrophysics/) · ASTRO-0091 | Fortran; Python interface | High-resolution plasma spectral modelling and fitting | SRON describes SPEX as GPLv3 open source with source distributed through Zenodo. Source text retrieved. |
| [ISIS spectral analysis](https://space.mit.edu/cxc/isis/) · ASTRO-0092 | C; S-Lang | Interactive high-resolution X-ray spectroscopy | Different ISIS from optical difference imaging. Source text retrieved. |
| [Stingray](https://github.com/StingraySoftware/stingray) · ASTRO-0093 | Python | Spectral timing, variability, cross spectra and light curves | Event timing and noise-statistics capability. Source text retrieved. |
| [HENDRICS](https://github.com/StingraySoftware/HENDRICS) · ASTRO-0094 | Python | High-energy timing workflows built on Stingray | Pipeline layer; avoid double-counting its numerical base. Source text retrieved. |
| [Gammapy](https://github.com/gammapy/gammapy) · ASTRO-0095 | Python | Gamma-ray event, map, spectrum and likelihood analysis | Observation and instrument-response data model reference. Source text retrieved. |
| [ctools and GammaLib](https://cta.irap.omp.eu/ctools/) · ASTRO-0096 | C++; Python | Cherenkov gamma-ray analysis and general high-energy likelihood library | Related library and tool suite. Source text retrieved. |
| [Fermitools](https://github.com/fermi-lat/Fermitools-conda) · ASTRO-0097 | C++; Python | Fermi LAT event selection, response and likelihood analysis | Mission suite distribution repository. Source text retrieved. |
| [Fermipy](https://github.com/fermiPy/fermipy) · ASTRO-0098 | Python | Higher-level Fermi LAT analysis workflows | Uses Fermitools; not an independent response implementation. Short page, metadata or access shell. |
| [threeML](https://github.com/threeML/threeML) · ASTRO-0099 | Python | Multi-instrument likelihood modelling | Plugin architecture for multi-messenger analysis. Source text retrieved. |
| [Naima](https://github.com/zblz/naima) · ASTRO-0100 | Python | Nonthermal radiation modelling and inference | Particle populations to observed spectra. Source text retrieved. |
| [ctapipe](https://github.com/cta-observatory/ctapipe) · ASTRO-0101 | Python | Cherenkov telescope low-level event reconstruction | Distinct from high-level Gammapy analysis. Source text retrieved. |
| [eSASS](https://erosita.mpe.mpg.de/edr/DataAnalysis/) · ASTRO-0102 | Mixed; consult distribution | eROSITA event and survey analysis | Release-specific accessibility and calibration dependencies. Source text retrieved. |
| [ixpeobssim](https://github.com/lucabaldini/ixpeobssim) · ASTRO-0103 | Python | X-ray polarimetry simulation and analysis | Include Stokes statistics and modulation response. Source text retrieved. |
| [SOXS and pyXSIM](https://github.com/lynx-x-ray-observatory/soxs) · ASTRO-0104 | Python | Synthetic X-ray observations and photons from simulations | Couples plasma emission, responses and simulation outputs. Source text retrieved. |
| [MARX](https://space.mit.edu/cxc/marx/) · ASTRO-0105 | C | Chandra observation and ray-trace simulation | Instrument simulator and validation reference. Source text retrieved. |

## 06 Radio interferometry and spectral cubes

| Software / source | Language family | Capability | Lineage, integration and evidence notes |
|---|---|---|---|
| [CASA](https://casadocs.readthedocs.io/en/stable/) · ASTRO-0106 | C++; Python | Radio calibration, imaging, deconvolution and single-dish analysis | Full suite; preserve Measurement Set and measurement-equation semantics. Source text retrieved. |
| [casacore](https://github.com/casacore/casacore) · ASTRO-0107 | C++ | Tables, measures, coordinates and radio astronomy data primitives | Foundational library distinct from CASA tasks. Source text retrieved. |
| [AIPS](https://www.aips.nrao.edu/) · ASTRO-0108 | Fortran; C; POPS | Radio interferometry calibration and imaging | Longstanding NRAO system; age alone does not imply abandoned. Source text retrieved. |
| [ParselTongue](https://ascl.net/1208.020) · ASTRO-0109 | Python | Python scripting of AIPS | Wrapper lineage and workflow compatibility. Source text retrieved. [Original nominated site](https://www.jive.eu/jivewiki/doku.php?id=parseltongue:parseltongue); source link uses ASCL. |
| [MIRIAD](https://www.atnf.csiro.au/computing/software/miriad/) · ASTRO-0110 | Fortran; C | Radio interferometric reduction and analysis | Legacy ecosystem with observatory-specific branches. Source text retrieved. |
| [GILDAS](https://www.iram.fr/IRAMFR/GILDAS/) · ASTRO-0111 | Fortran; C; SIC | CLASS spectroscopy, CLIC calibration and MAPPING imaging | CLASS here is unrelated to the cosmology Boltzmann code. Source text retrieved. |
| [WSClean](https://wsclean.readthedocs.io/) · ASTRO-0112 | C++ | Wide-field radio synthesis imaging and deconvolution | Modern imager; dependent on shared radio data libraries. Source text retrieved. |
| [AOFlagger](https://aoflagger.readthedocs.io/) · ASTRO-0113 | C++; Lua | Radio-frequency interference detection and flagging | Selection and flag provenance matter. Source text retrieved. |
| [DP3](https://dp3.readthedocs.io/) · ASTRO-0114 | C++; Python plugins | Visibility preprocessing, averaging and calibration | LOFAR and SKA lineage. Source text retrieved. |
| [DDFacet](https://github.com/saopicc/DDFacet) · ASTRO-0115 | Python; C/C++ | Direction-dependent wide-field radio imaging | Facet and beam-aware imaging. Source text retrieved. |
| [killMS](https://github.com/saopicc/killMS) · ASTRO-0116 | Python; native kernels | Direction-dependent radio calibration | Companion to DDFacet. Source text retrieved. |
| [CubiCal](https://github.com/ratt-ru/CubiCal) · ASTRO-0117 | Python; accelerated kernels | Radio interferometric calibration | Gain solutions and calibration equations. Short page, metadata or access shell. |
| [QuartiCal](https://github.com/ratt-ru/QuartiCal) · ASTRO-0118 | Python; Numba; Dask | Scalable radio interferometric calibration | Related modern calibration ecosystem. Source text retrieved. |
| [MeqTrees](https://github.com/ratt-ru/meqtrees-timba) · ASTRO-0119 | C++; Python | Numerical solution of measurement equations | Flexible equation-tree architecture. Source text retrieved. |
| [RASCIL and SKA SDP components](https://developer.skao.int/projects/rascil/en/latest/) · ASTRO-0120 | Python; C/C++ components | Simulation, calibration, imaging and radio data models | Parts migrated to ska-sdp-datamodels and ska-sdp-func-python. Source text retrieved. |
| [OSKAR](https://github.com/OxfordSKA/OSKAR) · ASTRO-0121 | C++; CUDA; Python | Radio interferometer and beam simulation | Instrument forward modelling. Source text retrieved. |
| [Stimela and CARACal](https://github.com/caracal-pipeline/caracal) · ASTRO-0122 | Python; containerized tools | Compose radio calibration and imaging workflows | Orchestration, not a distinct physical model. Source text retrieved. |
| [PyBDSF](https://github.com/lofar-astron/PyBDSF) · ASTRO-0123 | Python; Fortran | Radio image source detection and characterization | Former PyBDSM naming lineage. Source text retrieved. |
| [Aegean](https://github.com/PaulHancock/Aegean) · ASTRO-0124 | Python | Radio source finding and catalogue measurement | Includes background estimation and related tools. Source text retrieved. |
| [SoFiA 2](https://github.com/SoFiA-Admin/SoFiA-2) · ASTRO-0125 | C | Neutral-hydrogen spectral-cube source finding | The GitHub landing page points to a GitLab successor repository; acquisition must follow the migration. Short page, metadata or access shell. |
| [BBarolo](https://github.com/editeodoro/Bbarolo) · ASTRO-0126 | C++ | Three-dimensional tilted-ring galaxy kinematic modelling | HI and other resolved spectral cubes. Source text retrieved. |
| [TiRiFiC](https://github.com/gigjozsa/tirific) · ASTRO-0127 | C | Tilted-ring fitting directly to spectral cubes | Alternative forward model and fitting strategy. Short page, metadata or access shell. |
| [spectral-cube and radio-beam](https://github.com/radio-astro-tools/spectral-cube) · ASTRO-0128 | Python | WCS-aware cube analysis and radio beam operations | Beam and brightness-temperature equivalencies are essential. Source text retrieved. |
| [GBTIDL](https://gbtidl.nrao.edu/) · ASTRO-0129 | IDL | Green Bank Telescope single-dish spectral analysis | IDL observatory workflow. Source text retrieved. |
| [dysh](https://github.com/GreenBankObservatory/dysh) · ASTRO-0130 | Python | Single-dish radio spectral analysis | Python ecosystem for GBT and single-dish workflows. Source text retrieved. |
| [eht-imaging](https://github.com/achael/eht-imaging) · ASTRO-0131 | Python | VLBI imaging and modelling for black-hole-scale sources | Visibility-domain likelihoods and closure quantities. Source text retrieved. |
| [SMILI](https://github.com/astrosmili/smili) · ASTRO-0132 | Python; C; Fortran | Sparse-model VLBI image reconstruction | Alternative EHT imaging approach. Source text retrieved. |
| [Difmap](https://sites.astro.caltech.edu/~mcs/difmap/) · ASTRO-0133 | C | VLBI visibility editing, imaging and model fitting | Longstanding specialist tool. Source text retrieved. |
| [DiFX](https://github.com/difx/difx) · ASTRO-0134 | C++; Python; other components | Software VLBI correlation | Correlation is a distinct upstream capability. Source text retrieved. |
| [pyuvdata and pyuvsim](https://github.com/RadioAstronomySoftwareGroup/pyuvdata) · ASTRO-0135 | Python | Radio visibility, beam and calibration formats; simulations | Shared objects and interoperability for low-frequency arrays. Source text retrieved. |

## 07 Pulsars fast transients and timing

| Software / source | Language family | Capability | Lineage, integration and evidence notes |
|---|---|---|---|
| [PRESTO](https://github.com/scottransom/presto) · ASTRO-0136 | C; Python | Pulsar searches, acceleration searches and analysis | Filterbank and time-domain signal-processing chain. Source text retrieved. |
| [PSRCHIVE](https://psrchive.sourceforge.net/) · ASTRO-0137 | C++; Python | Pulsar archives, calibration and profile analysis | Polarization and profile data model. Short page, metadata or access shell. |
| [DSPSR](https://dspsr.sourceforge.net/) · ASTRO-0138 | C++ | Pulsar coherent dedispersion and signal processing | Streaming baseband and folded data. Short page, metadata or access shell. |
| [TEMPO](https://tempo.sourceforge.net/) · ASTRO-0139 | Fortran | Pulsar timing model fitting | Legacy reference with scientific importance. Source text retrieved. |
| [TEMPO2](https://bitbucket.org/psrsoft/tempo2/) · ASTRO-0140 | C++; C; plugins | Precision pulsar timing | Different numerical and parameter conventions require validation. Short page, metadata or access shell. |
| [PINT pulsar timing](https://github.com/nanograv/PINT) · ASTRO-0141 | Python | Pulsar timing models and fitting | Not the unrelated Python units package named Pint. Source text retrieved. |
| [enterprise](https://github.com/nanograv/enterprise) · ASTRO-0142 | Python | Pulsar timing array noise and gravitational-wave inference | Correlated timing-residual likelihoods. Source text retrieved. |
| [PSRSIGSIM](https://github.com/PsrSigSim/PsrSigSim) · ASTRO-0143 | Python | Pulsar signal and observation simulation | Synthetic-data validation. Source text retrieved. |
| [SIGPROC](https://sigproc.sourceforge.net/) · ASTRO-0144 | C; Fortran | Filterbank processing and pulsar search utilities | Legacy file formats and algorithms remain useful. Source text retrieved. |
| [Heimdall](https://sourceforge.net/projects/heimdall-astro/) · ASTRO-0145 | C++; CUDA | GPU single-pulse transient searches | Fast-radio-burst search lineage. Source text retrieved. |
| [baseband](https://github.com/mhvk/baseband) · ASTRO-0146 | Python | Baseband data access and processing | Precise timestamps and sample formats. Source text retrieved. |

## 08 Strong weak and microlensing

| Software / source | Language family | Capability | Lineage, integration and evidence notes |
|---|---|---|---|
| [lenstronomy](https://github.com/lenstronomy/lenstronomy) · ASTRO-0147 | Python | Strong-lens imaging, mass models, time delays and cosmography | Broad reference for reusable lens primitives. Source text retrieved. |
| [PyAutoLens](https://github.com/Jammy2211/PyAutoLens) · ASTRO-0148 | Python | Automated strong-lens modelling and source reconstruction | Related PyAutoGalaxy and PyAutoFit ecosystem. Source text retrieved. |
| [Herculens](https://github.com/Herculens/herculens) · ASTRO-0149 | Python; JAX | Differentiable strong-lens modelling | Gradient-based pixel and parametric models. Source text retrieved. |
| [caustics](https://github.com/Ciela-Institute/caustics) · ASTRO-0150 | Python; PyTorch | Differentiable gravitational-lensing simulations | Accelerator-oriented lens forward models. Source text retrieved. |
| [jaxtronomy](https://github.com/lenstronomy/JAXtronomy) · ASTRO-0151 | Python; JAX | JAX implementations of lenstronomy-style lens calculations | Related implementation lineage. Source text retrieved. |
| [Lenstool](https://ascl.net/1102.004) · ASTRO-0152 | C; supporting interfaces | Parametric strong and weak cluster lens modelling | Cluster-scale mass reconstruction. Source text retrieved. [Original nominated site](https://projets.lam.fr/projects/lenstool/wiki); source link uses ASCL. |
| [glafic](https://github.com/oguri/glafic2) · ASTRO-0153 | C | Strong-lensing calculations and mass modelling | Point and extended sources; cluster and galaxy applications. Source text retrieved. |
| [gravlens and lensmodel](https://www.physics.rutgers.edu/~keeton/gravlens/) · ASTRO-0154 | Compiled code; language requires inspection | Lens calculations, image solving and model optimization | Legacy specialist code; distribution rights and source access need review. Short page, metadata or access shell. |
| [GLEE](https://arxiv.org/abs/2209.03094) · ASTRO-0155 | Implementation language not confirmed in this survey | Strong-lens mass and source reconstruction | Author paper documents strong-lens modelling and comparisons; public source-distribution status not established. Source text retrieved. |
| [GLASS](https://github.com/jpcoles/glass) · ASTRO-0156 | Python; native components | Free-form gravitational-lens mass modelling | Distinct from the cosmological simulation package also called GLASS. Source text retrieved. |
| [WSLAP plus](https://arxiv.org/abs/1304.2393) · ASTRO-0157 | Implementation language not confirmed in this survey | Free-form cluster lens reconstruction | Primary author paper documents WSLAP+ reconstruction; surviving source distribution requires follow-up. Source text retrieved. |
| [PyCS3](https://github.com/COSMOGRAIL/PyCS) · ASTRO-0158 | Python | Time-delay estimation from lensed quasar light curves | Time series and microlensing variability treatment. Source text retrieved. |
| [hierArc](https://github.com/TDCOSMO/hierArc) · ASTRO-0159 | Python | Hierarchical strong-lens cosmography inference | Population and nuisance assumptions sit above image modelling. Source text retrieved. |
| [SNTD](https://github.com/jpierel14/sntd) · ASTRO-0160 | Python | Analysis of strongly lensed supernova light curves | Connects transient and lensing modules. Source text retrieved. |
| [MulensModel](https://github.com/rpoleski/MulensModel) · ASTRO-0161 | Python; C/C++ kernels | Microlensing light curves and event modelling | Finite-source and binary-lens methods. Source text retrieved. |
| [pyLIMA](https://github.com/ebachelet/pyLIMA) · ASTRO-0162 | Python | Microlensing event modelling and fitting | Independent fitting and workflow reference. Source text retrieved. |
| [VBMicrolensing](https://github.com/valboz/VBMicrolensing) · ASTRO-0163 | C++; Python bindings | Fast binary and multiple-lens magnification calculations | Successor family to VBBinaryLensing. Source text retrieved. |
| [GalSim](https://github.com/GalSim-developers/GalSim) · ASTRO-0164 | C++; Python | Realistic astronomical image and weak-lensing simulations | PSFs, shear, noise and detector effects; not an N-body code. Source text retrieved. |
| [ngmix](https://github.com/esheldon/ngmix) · ASTRO-0165 | Python; compiled acceleration | Galaxy and PSF model fitting and shear estimation | Weak-lensing measurement algorithms. Source text retrieved. |
| [metadetect](https://github.com/esheldon/metadetect) · ASTRO-0166 | Python | Detection-aware metacalibration for weak lensing | Selection effects are part of the estimator. Short page, metadata or access shell. |
| [TreeCorr](https://github.com/rmjarvis/TreeCorr) · ASTRO-0167 | C++; Python | Two-point and three-point correlation functions | Shared large-scale structure and shear statistics. Source text retrieved. |
| [LensTools](https://github.com/apetri/LensTools) · ASTRO-0168 | Python; native components | Weak-lensing maps, simulations and analysis | Map-level statistics and ray-tracing workflows. Short page, metadata or access shell. |

## 09 Cosmology CMB and large scale structure

| Software / source | Language family | Capability | Lineage, integration and evidence notes |
|---|---|---|---|
| [CAMB](https://github.com/cmbant/CAMB) · ASTRO-0169 | Fortran; Python | Linear cosmological perturbations and CMB/matter predictions | Boltzmann solver; distinct from inference packages. Source text retrieved. |
| [CLASS cosmology](https://github.com/lesgourg/class_public) · ASTRO-0170 | C; Python classy wrapper | Background and perturbation evolution and spectra | Distinct from IRAM CLASS spectroscopy. Source text retrieved. |
| [Cobaya](https://github.com/CobayaSampler/cobaya) · ASTRO-0171 | Python | Theory/likelihood composition and Bayesian sampling | Framework can share CAMB or CLASS engines. Source text retrieved. |
| [CosmoMC](https://github.com/cmbant/CosmoMC) · ASTRO-0172 | Fortran; Python tooling | Cosmological Monte Carlo parameter inference | Historical and specialist inference lineage. Source text retrieved. |
| [MontePython](https://github.com/brinckmann/montepython_public) · ASTRO-0173 | Python | Cosmological likelihood and parameter inference | Typically coupled to CLASS. Source text retrieved. |
| [CosmoSIS](https://github.com/joezuntz/cosmosis) · ASTRO-0174 | Python; C; C++; Fortran modules | Modular cosmological parameter-estimation pipelines | Useful multi-language component architecture. Short page, metadata or access shell. |
| [CCL and pyccl](https://github.com/LSSTDESC/CCL) · ASTRO-0175 | C; Python | Cosmological predictions for large-scale structure | Core Cosmology Library. Source text retrieved. |
| [Colossus](https://bitbucket.org/bdiemer/colossus/) · ASTRO-0176 | Python | Cosmology, halo profiles, masses and large-scale structure | Useful halo-physics primitives. Short page, metadata or access shell. |
| [jax-cosmo](https://github.com/DifferentiableUniverseInitiative/jax_cosmo) · ASTRO-0177 | Python; JAX | Differentiable background and projected large-scale-structure spectra | Approximations and probe scope differ from full Boltzmann solvers. Source text retrieved. |
| [CosmoPower and CosmoPower-JAX](https://github.com/alessiospuriomancini/cosmopower) · ASTRO-0178 | Python; TensorFlow or JAX | Neural emulators of cosmological predictions | Learned model domain and training artefacts are part of the scientific product. Source text retrieved. |
| [CLASS-PT](https://github.com/Michalychforever/CLASS-PT) · ASTRO-0179 | C; Python | Perturbation-theory galaxy clustering predictions | Extends CLASS; shared lineage. Source text retrieved. |
| [PyBird](https://github.com/pierrexyz/pybird) · ASTRO-0180 | Python | Effective-field-theory large-scale-structure predictions | Nonlinear galaxy clustering analysis. Source text retrieved. |
| [velocileptors](https://github.com/sfschen/velocileptors) · ASTRO-0181 | Python | Perturbative large-scale-structure and redshift-space models | Alternative perturbative implementation. Source text retrieved. |
| [Corrfunc](https://github.com/manodeep/Corrfunc) · ASTRO-0182 | C; Python | Fast pair counting for correlation functions | Reusable performance-critical statistics kernel. Source text retrieved. |
| [nbodykit](https://github.com/bccp/nbodykit) · ASTRO-0183 | Python; MPI; native dependencies | Large-scale structure catalogues, meshes and statistics | Historical dependency stack; inspect maintenance before adoption. Source text retrieved. |
| [pycorr and pypower](https://github.com/cosmodesi/pycorr) · ASTRO-0184 | Python; native dependencies | Survey correlation functions and power-spectrum estimation | DESI ecosystem; pypower is a separate companion repository. Source text retrieved. |
| [NaMaster](https://github.com/LSSTDESC/NaMaster) · ASTRO-0185 | C; Python pymaster | Masked-sky angular power-spectrum estimation | Mode-coupling and spin-field treatment. Source text retrieved. |
| [pixell](https://github.com/simonsobs/pixell) · ASTRO-0186 | Python; native libraries | CMB maps, curved-sky transforms and analysis | Map geometry and Fourier conventions. Source text retrieved. |
| [TOAST](https://github.com/hpc4cmb/toast) · ASTRO-0187 | C++; Python | Time-ordered CMB data simulation and mapmaking | Distributed observation-level workflows. Source text retrieved. |
| [sotodlib](https://github.com/simonsobs/sotodlib) · ASTRO-0188 | Python | Simons Observatory time-stream analysis | Instrument-specific analysis atop common numerical infrastructure. Short page, metadata or access shell. |
| [Commander](https://github.com/Cosmoglobe/Commander) · ASTRO-0189 | Fortran; supporting languages | Bayesian CMB component separation | Foreground and instrument modelling. Source text retrieved. |
| [SMICA and Planck likelihood code](https://pla.esac.esa.int/) · ASTRO-0190 | Mixed compiled code; Python interfaces | CMB component separation and likelihood evaluation | Family and archive entry; individual source accessibility varies. Short page, metadata or access shell. |
| [lenspyx](https://github.com/carronj/lenspyx) · ASTRO-0191 | Python; native transforms | Curved-sky lensing deflections and map remapping | CMB lensing simulation. Source text retrieved. |
| [21cmFAST](https://github.com/21cmfast/21cmFAST) · ASTRO-0192 | C; Python | Semi-numerical reionization and 21 cm simulations | Distinct early-universe modelling assumptions. Source text retrieved. |
| [21CMMC](https://github.com/21cmfast/21CMMC) · ASTRO-0193 | Python | Bayesian inference around 21cmFAST | Shared simulation engine. Source text retrieved. |
| [GetDist](https://github.com/cmbant/getdist) · ASTRO-0194 | Python | Posterior-chain summaries and visualization | Canonical project name is GetDist; works beyond Cobaya. Source text retrieved. |

## 10 Supernovae and transient ecosystems

| Software / source | Language family | Capability | Lineage, integration and evidence notes |
|---|---|---|---|
| [SNANA](https://github.com/RickKessler/SNANA) · ASTRO-0195 | C; Fortran; Python; historical Perl | Supernova simulation, fitting and cosmology workflows | Survey selection, calibration and bias correction are major capabilities. Source text retrieved. |
| [SNCosmo](https://github.com/sncosmo/sncosmo) · ASTRO-0196 | Python | Transient light-curve models, synthetic photometry and fitting | Shared passband and model data dependencies. Short page, metadata or access shell. |
| [BayeSN](https://github.com/bayesn/bayesn) · ASTRO-0197 | Python; JAX; NumPyro | Hierarchical Type Ia optical/NIR emission and distance modelling | Different population and dust assumptions from SALT workflows. Source text retrieved. |
| [SALTShaker](https://github.com/djones1040/SALTShaker) · ASTRO-0198 | Python | Train SALT-family supernova spectral models | Training and fitting are distinct stages. Source text retrieved. |
| [SNooPy](https://github.com/obscode/snpy) · ASTRO-0199 | Python; native components | Supernova light-curve fitting and distance methods | Carnegie Supernova Project lineage. Source text retrieved. |
| [SNID](https://people.lam.fr/blondin.stephane/software/snid/) · ASTRO-0200 | Fortran | Supernova spectral classification and redshift matching | Template library is essential scientific input. Source text retrieved. |
| [Superfit](https://github.com/dahowell/superfit) · ASTRO-0201 | IDL | Supernova spectral-template fitting | Legacy proprietary-language workflow. Source text retrieved. |
| [SuperNNova](https://github.com/supernnova/SuperNNova) · ASTRO-0202 | Python; PyTorch | Neural supernova light-curve classification | Training populations and calibration determine applicability. Source text retrieved. |
| [MOSFiT](https://github.com/guillochon/MOSFiT) · ASTRO-0203 | Python | Physical transient light-curve modelling and inference | Modular model architecture. Source text retrieved. |
| [Redback](https://github.com/nikhil-sarin/redback) · ASTRO-0204 | Python | Bayesian transient and multi-messenger modelling | Links physical models and inference engines. Source text retrieved. |
| [TARDIS](https://github.com/tardis-sn/tardis) · ASTRO-0205 | Python; Numba; native components | Monte Carlo supernova spectral synthesis | Atomic data and radiative-transfer assumptions are separate dependencies. Source text retrieved. |
| [SEDONA](https://github.com/dnkasen/pubsed) · ASTRO-0206 | C++ | Time-dependent radiation transport for transients | Specialist simulation reference. Source text retrieved. |
| [SNEC](https://stellarcollapse.org/index.php/SNEC.html) · ASTRO-0207 | Fortran | Supernova radiation-hydrodynamics and light curves | Explosion modelling and microphysics inputs. Source text retrieved. |
| [SkyPortal](https://github.com/skyportal/skyportal) · ASTRO-0208 | Python; JavaScript | Collaborative transient data and follow-up management | Application/workflow capability, not a physical solver. Source text retrieved. |
| [TOM Toolkit](https://github.com/TOMToolkit/tom_base) · ASTRO-0209 | Python; JavaScript | Target observation management and telescope integration | Operational scheduling and provenance. Source text retrieved. |
| [Fink](https://github.com/astrolabsoftware/fink-broker) · ASTRO-0210 | Python; Spark | Transient alert processing and scientific enrichment | Distributed broker; external services and trained models. Source text retrieved. |
| [Lasair](https://github.com/lsst-uk/lasair-lsst) · ASTRO-0211 | Python; service stack | Transient alert ingestion, filtering and discovery | Service architecture rather than standalone numerical library. Source text retrieved. |
| [ALeRCE](https://github.com/alercebroker) · ASTRO-0212 | Python; service stack | Alert classification and time-domain discovery | Organization-level family entry. Source text retrieved. |

## 11 Stars populations and spectral synthesis

| Software / source | Language family | Capability | Lineage, integration and evidence notes |
|---|---|---|---|
| [MESA](https://github.com/MESAHub/mesa) · ASTRO-0213 | Fortran | Detailed stellar structure and evolution with microphysics | Equation-of-state, opacity and nuclear-network layers matter independently. Source text retrieved. |
| [GYRE](https://github.com/rhdtownsend/gyre) · ASTRO-0214 | Fortran | Stellar oscillation and asteroseismic calculations | Can consume MESA stellar models. Source text retrieved. |
| [ADIPLS](https://users-phys.au.dk/jcd/adipack.n/) · ASTRO-0215 | Fortran | Adiabatic stellar pulsation calculations | The cited landing page explicitly describes an obsolete adipack version and links to its replacement; retain the lineage. Short page, metadata or access shell. |
| [SSE and BSE](https://astronomy.swin.edu.au/~jhurley/) · ASTRO-0216 | Fortran | Rapid single and binary stellar evolution | Legacy fitting-formula lineage used by later population codes. Retrieval unresolved. |
| [COSMIC](https://github.com/COSMIC-PopSynth/COSMIC) · ASTRO-0217 | Python; Fortran | Compact-binary population synthesis | Builds on SSE/BSE lineage. Source text retrieved. |
| [COMPAS](https://github.com/TeamCOMPAS/COMPAS) · ASTRO-0218 | C++; Python analysis | Rapid binary population synthesis | Alternative binary-evolution assumptions. Source text retrieved. |
| [POSYDON](https://github.com/POSYDON-code/POSYDON) · ASTRO-0219 | Python; MESA grids | Population synthesis with detailed stellar/binary evolution grids | Grid generation and interpolation are distinct from inference. Source text retrieved. |
| [binary_c](https://ascl.net/2307.035) · ASTRO-0220 | C | Binary population synthesis and nucleosynthesis | Specialist evolution framework. Source text retrieved. [Original nominated site](https://binary_c.gitlab.io/); source link uses ASCL. |
| [BPASS](https://bpass.auckland.ac.nz/) · ASTRO-0221 | Mixed; stellar model and data distribution | Binary population synthesis and spectral predictions | Separate downloadable grids from code access. Short page, metadata or access shell. |
| [FSPS and python-fsps](https://github.com/cconroy20/fsps) · ASTRO-0222 | Fortran; Python | Flexible stellar population synthesis | Stellar libraries and isochrones affect predictions. Source text retrieved. |
| [Starburst99](https://www.stsci.edu/science/starburst99/docs/default.htm) · ASTRO-0223 | Fortran; supporting tools | Population synthesis for young stellar systems | Legacy and model-grid reference. Source text retrieved. |
| [Prospector](https://github.com/bd-j/prospector) · ASTRO-0224 | Python | Bayesian stellar-population and galaxy SED inference | Commonly uses FSPS; shared models imply dependence. Source text retrieved. |
| [Bagpipes](https://github.com/ACCarnall/bagpipes) · ASTRO-0225 | Python | Galaxy spectral-energy-distribution and star-formation-history fitting | Model assumptions and priors explicit. Source text retrieved. |
| [CIGALE](https://cigale.lam.fr/) · ASTRO-0226 | Python | Multiwavelength galaxy SED modelling | Energy-balance model family. Source text retrieved. |
| [MAGPHYS](https://www.iap.fr/magphys/) · ASTRO-0227 | Fortran; model libraries | Galaxy SED fitting with energy balance | Distribution and model-library access need separate review. Source text retrieved. |
| [LePhare](https://github.com/lephare-photoz/lephare) · ASTRO-0228 | C++; Python; historical Fortran | Photometric redshifts and SED fitting | Track current rewrite and historical algorithms separately. Source text retrieved. |
| [EAZY](https://github.com/gbrammer/eazy-py) · ASTRO-0229 | Python; historical C | Photometric redshift fitting | Template and prior dependence. Source text retrieved. |
| [kcorrect](https://github.com/blanton144/kcorrect) · ASTRO-0230 | Python; C; historical IDL interfaces | K corrections and galaxy SED reconstruction | Long-lived multi-language lineage. Short page, metadata or access shell. |
| [MOOG](https://www.as.utexas.edu/~chris/moog.html) · ASTRO-0231 | Fortran | LTE stellar line analysis and spectral synthesis | Classic stellar abundance software. Source text retrieved. |
| [Turbospectrum](https://github.com/bertrandplez/Turbospectrum2019) · ASTRO-0232 | Fortran | Stellar spectral synthesis | Atmosphere and line-list resources required. Source text retrieved. |
| [SME and PySME](https://github.com/AWehrhahn/SME) · ASTRO-0233 | C/C++; IDL or Python interfaces | Stellar spectral fitting and parameter estimation | Spectroscopy Made Easy lineage. Source text retrieved. |
| [iSpec](https://github.com/marblestation/iSpec) · ASTRO-0234 | Python; external synthesis engines | Stellar spectral analysis with multiple synthesis backends | Integration example and algorithm comparison environment. Source text retrieved. |
| [ATLAS and SYNTHE](https://wwwuser.oats.inaf.it/castelli/sources.html) · ASTRO-0235 | Fortran | Stellar atmospheres and line-blanketed synthetic spectra | Legacy numerical and line-list lineage. Short page, metadata or access shell. |
| [TLUSTY and SYNSPEC](https://tlusty.oca.eu/) · ASTRO-0236 | Fortran | Non-LTE stellar atmospheres and spectra | Distinct physical approximation regime. Short page, metadata or access shell. |
| [PHOENIX](https://ascl.net/1010.056) · ASTRO-0237 | Fortran; mixed infrastructure | Stellar and planetary atmosphere radiative transfer | Public grids do not imply full code is freely available. Source text retrieved. [Original nominated site](https://www.hs.uni-hamburg.de/EN/For/ThA/phoenix/); source link uses ASCL. |

## 12 Radiation chemistry and atomic physics

| Software / source | Language family | Capability | Lineage, integration and evidence notes |
|---|---|---|---|
| [Cloudy](https://gitlab.nublado.org/cloudy/cloudy) · ASTRO-0238 | C++ | Photoionization, plasma microphysics and emitted spectra | Code and atomic/molecular data both required. Short page, metadata or access shell. |
| [XSTAR](https://heasarc.gsfc.nasa.gov/lheasoft/xstar/xstar.html) · ASTRO-0239 | Fortran | Photoionized gas and X-ray spectral models | Also distributed with HEASoft. Source text retrieved. |
| [CHIANTI and ChiantiPy](https://www.chiantidatabase.org/) · ASTRO-0240 | IDL; Python; atomic database | Optically thin plasma emission and diagnostics | Model database and two language ecosystems. Source text retrieved. |
| [AtomDB and PyAtomDB](https://github.com/AtomDB/pyatomdb) · ASTRO-0241 | Atomic database; Python; compiled emissivity tools | X-ray plasma emissivities and spectral modelling | Data versions are part of the scientific model. Source text retrieved. |
| [PyNeb](https://github.com/Morisset/PyNeb_devel) · ASTRO-0242 | Python | Emission-line diagnostics and nebular abundances | Atomic coefficients and uncertainty propagation. Source text retrieved. |
| [RADMC-3D](https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/) · ASTRO-0243 | Fortran; Python tools | Dust and line radiative transfer | Monte Carlo and ray tracing. Source text retrieved. |
| [MCFOST](https://github.com/cpinte/mcfost) · ASTRO-0244 | Fortran | Three-dimensional dust and line radiative transfer | Disk and circumstellar applications. Source text retrieved. |
| [SKIRT](https://github.com/SKIRT/SKIRT9) · ASTRO-0245 | C++ | Three-dimensional Monte Carlo radiative transfer | Dust, gas and synthetic observations. Source text retrieved. |
| [Hyperion](https://github.com/hyperion-rt/hyperion) · ASTRO-0246 | Fortran; Python | Dust continuum radiative transfer | Radiation/density grid abstractions. Source text retrieved. |
| [LIME](https://github.com/lime-rt/lime) · ASTRO-0247 | C | Non-LTE molecular line transfer | Irregular-grid radiative transfer. Source text retrieved. |
| [RADEX](https://ascl.net/1010.075) · ASTRO-0248 | Fortran | Non-LTE molecular excitation and line intensities | Escape-probability approximation; not full spatial transport. Source text retrieved. [Original nominated site](https://personal.sron.nl/~vdtak/radex/index.shtml); source link uses ASCL. |
| [MOCASSIN](https://mocassin.nebulousresearch.org/) · ASTRO-0249 | Fortran | Three-dimensional photoionization and dust transfer | Specialist Monte Carlo lineage. Source text retrieved. |
| [TORUS](https://ascl.net/1404.006) · ASTRO-0250 | Fortran | Radiation transport and radiation hydrodynamics | Model and distribution scope require inspection. Source text retrieved. [Original nominated site](https://www.astro.ex.ac.uk/people/th2/torus_html/); source link uses ASCL. |
| [KROME](https://bitbucket.org/tgrassi/krome/) · ASTRO-0251 | Python generator; Fortran output | Chemical reaction networks and thermal microphysics | Generated solver architecture. Short page, metadata or access shell. |
| [Grackle](https://github.com/grackle-project/grackle) · ASTRO-0252 | C; Fortran; Python | Cooling, chemistry and radiative heating for simulations | Reusable microphysics library. Source text retrieved. |
| [DESPOTIC](https://bitbucket.org/krumholz/despotic/) · ASTRO-0253 | Python | Interstellar cloud thermal and chemical modelling | Reduced physical cloud models. Short page, metadata or access shell. |
| [pynucastro](https://github.com/pynucastro/pynucastro) · ASTRO-0254 | Python; generated C++ and Fortran | Nuclear reaction networks and rate handling | Rates and network generation are separate capabilities. Source text retrieved. |

## 13 Dynamics hydrodynamics and cosmological simulations

| Software / source | Language family | Capability | Lineage, integration and evidence notes |
|---|---|---|---|
| [GADGET-4](https://wwwmpa.mpa-garching.mpg.de/gadget4/) · ASTRO-0255 | C++ | Cosmological N-body and hydrodynamical simulations | Separate historical GADGET-2/3 branches and descendants. Source text retrieved. |
| [GADGET-2](https://wwwmpa.mpa-garching.mpg.de/gadget/) · ASTRO-0256 | C | Tree/SPH cosmological simulations | Legacy public reference with many descendants. Short page, metadata or access shell. |
| [AREPO](https://arepo-code.org/) · ASTRO-0257 | C | Moving-mesh hydrodynamics and gravity | Public version and collaboration extensions may differ. Source text retrieved. |
| [GIZMO](https://ascl.net/1410.003) · ASTRO-0258 | C | Gravity and multiple hydrodynamic methods | Public branch and optional physics vary. Source text retrieved. [Original nominated site](https://www.tapir.caltech.edu/~phopkins/Site/GIZMO.html); source link uses ASCL. |
| [SWIFT](https://github.com/SWIFTSIM/SWIFT) · ASTRO-0259 | C | Task-based cosmological gravity and hydrodynamics | Parallel execution and subgrid physics framework. Source text retrieved. |
| [RAMSES](https://github.com/ramses-organisation/ramses) · ASTRO-0260 | Fortran | Adaptive-mesh cosmological hydrodynamics and gravity | AMR and particle-mesh methods. Source text retrieved. |
| [Enzo](https://github.com/enzo-project/enzo-dev) · ASTRO-0261 | C++; Fortran; Python | Adaptive-mesh astrophysical fluid simulations | Longstanding AMR code. Source text retrieved. |
| [Enzo-E](https://github.com/enzo-project/enzo-e) · ASTRO-0262 | C++ | Scalable adaptive-mesh astrophysics | Distinct implementation using Cello infrastructure. Source text retrieved. |
| [FLASH](https://flash.rochester.edu/site/flashcode/) · ASTRO-0263 | Fortran; C | Adaptive-mesh multiphysics and astrophysical explosions | Access and redistribution terms require review. Source text retrieved. |
| [Athena and Athena++](https://github.com/PrincetonUniversity/athena) · ASTRO-0264 | C; C++ | Astrophysical hydrodynamics and MHD | Athena++ is distinct from the original C implementation. Source text retrieved. |
| [AthenaK](https://github.com/IAS-Astrophysics/athenak) · ASTRO-0265 | C++; Kokkos | Portable accelerated astrophysical fluid and relativity simulations | GPU-capable successor family. Source text retrieved. |
| [PLUTO](https://plutocode.ph.unito.it/) · ASTRO-0266 | C | Astrophysical fluid, MHD and relativistic flows | Multiple physical/numerical modules. Source text retrieved. |
| [MPI-AMRVAC](https://github.com/amrvac/amrvac) · ASTRO-0267 | Fortran | Adaptive-mesh hydrodynamics, MHD and related physics | Extensible PDE framework. Short page, metadata or access shell. |
| [Pencil Code](https://github.com/pencil-code/pencil-code) · ASTRO-0268 | Fortran; Python | High-order compressible MHD and multiphysics | Turbulence and dynamo applications. Source text retrieved. |
| [Phantom](https://github.com/danieljprice/phantom) · ASTRO-0269 | Fortran | Smoothed-particle hydrodynamics and astrophysical dynamics | Disks and other Lagrangian-fluid applications. Source text retrieved. |
| [ChaNGa](https://github.com/N-BodyShop/changa) · ASTRO-0270 | C++; Charm++ | Parallel cosmological N-body and SPH simulations | Task/runtime design reference. Source text retrieved. |
| [PKDGRAV3](https://github.com/CLS-project/PKDGRAV3) · ASTRO-0271 | C; C++ | Large cosmological N-body gravity simulations | Performance-critical gravity implementation. Source text retrieved. |
| [REBOUND](https://github.com/hannorein/rebound) · ASTRO-0272 | C; Python | N-body orbital integration | Multiple integrators and collision treatments. Source text retrieved. |
| [REBOUNDx](https://github.com/dtamayo/reboundx) · ASTRO-0273 | C; Python | Additional forces and physical effects for REBOUND | Extension mechanism reference. Source text retrieved. |
| [NBODY6 and NBODY6++GPU](https://github.com/nbodyx/Nbody6ppGPU) · ASTRO-0274 | Fortran; C/C++; CUDA | Direct collisional star-cluster dynamics | Legacy and accelerated family; regularization and close encounters. Source text retrieved. |
| [PeTar](https://github.com/lwang-astro/PeTar) · ASTRO-0275 | C++ | Star-cluster dynamics with hybrid gravitational methods | Coupling dynamics and stellar evolution. Source text retrieved. |
| [CMC](https://github.com/ClusterMonteCarlo/CMC-COSMIC) · ASTRO-0276 | C; Python interfaces | Monte Carlo star-cluster evolution | Related COSMIC evolution dependencies. Source text retrieved. |
| [galpy](https://github.com/jobovy/galpy) · ASTRO-0277 | Python; C | Galactic orbits, potentials and distribution functions | Useful reusable dynamics core. Source text retrieved. |
| [gala](https://github.com/adrn/gala) · ASTRO-0278 | Python; C/Cython | Galactic dynamics and orbit analysis | Astropy-compatible phase-space objects. Source text retrieved. |
| [AGAMA](https://github.com/GalacticDynamics-Oxford/Agama) · ASTRO-0279 | C++; Python | Galaxy dynamics, potentials, actions and distribution functions | Alternative dynamics implementation. Source text retrieved. |
| [NEMO](https://github.com/teuben/nemo) · ASTRO-0280 | C; C++; Fortran tools | Stellar dynamics toolbox and simulation analysis | Legacy integrated toolbox with ongoing stewardship. Source text retrieved. |
| [yt](https://github.com/yt-project/yt) · ASTRO-0281 | Python; Cython | Analysis and visualization of volumetric simulation data | Reader/analysis layer rather than a simulation engine. Source text retrieved. |
| [pynbody](https://github.com/pynbody/pynbody) · ASTRO-0282 | Python; C/Cython | N-body and hydrodynamic simulation analysis | Particle data, units, halos and derived quantities. Source text retrieved. |

## 14 Relativity spacetime and black holes

| Software / source | Language family | Capability | Lineage, integration and evidence notes |
|---|---|---|---|
| [Einstein Toolkit](https://einsteintoolkit.org/) · ASTRO-0283 | C; C++; Fortran; Python | Numerical relativity, initial data, spacetime evolution and diagnostics | A suite of separately maintained thorns; Cactus and Carpet/CarpetX infrastructure. Source text retrieved. |
| [Cactus and CarpetX](https://github.com/eschnett/CarpetX) · ASTRO-0284 | C++; C; Fortran components | Modular scientific simulation framework and AMR driver | Infrastructure related to Einstein Toolkit; not an independent gravity theory. Source text retrieved. |
| [SpECTRE](https://github.com/sxs-collaboration/spectre) · ASTRO-0285 | C++; Python | Parallel relativistic astrophysics using high-order numerical methods | Distinct from the collaboration's SpEC code. Source text retrieved. |
| [SpEC](https://www.black-holes.org/code/SpEC.html) · ASTRO-0286 | C++; other components | Spectral numerical relativity and compact-binary simulations | Code accessibility differs from public waveform availability. Short page, metadata or access shell. |
| [GRChombo](https://github.com/GRChombo/GRChombo) · ASTRO-0287 | C++ | Adaptive-mesh numerical relativity | Primary repository identifies GRChombo as predecessor to GRTeclyn; retain both implementations. Source text retrieved. |
| [NRPy](https://github.com/nrpy/nrpy) · ASTRO-0288 | Python; generated C/CUDA | Symbolic generation of numerical-relativity solvers | Useful symbolic-to-numerical architecture reference. Source text retrieved. |
| [EinsteinPy](https://github.com/einsteinpy/einsteinpy) · ASTRO-0289 | Python | Symbolic relativity, geodesics and spacetime calculations | Not a general replacement for numerical-relativity evolution codes. Source text retrieved. |
| [xAct](https://www.xact.es/) · ASTRO-0290 | Wolfram Language | Tensor algebra, perturbations and differential geometry | Open package ecosystem with proprietary Mathematica runtime. Source text retrieved. |
| [Cadabra](https://cadabra.science/) · ASTRO-0291 | C++; Python | Symbolic tensor algebra and field-theory calculations | Open symbolic environment. Source text retrieved. |
| [SageManifolds](https://sagemanifolds.obspm.fr/) · ASTRO-0292 | Python; SageMath | Differential geometry, tensor fields and relativity | Open symbolic alternative with numerical interfaces. Source text retrieved. |
| [GRTensorIII](https://github.com/grtensor/grtensor) · ASTRO-0293 | Maple | Tensor calculations in general relativity | Legacy GRTensorII lineage and proprietary Maple runtime. Source text retrieved. |
| [Black Hole Perturbation Toolkit](https://bhptoolkit.org/) · ASTRO-0294 | Wolfram Language; Python; C/C++; others | Perturbations, Kerr orbits, waveforms and special functions | Umbrella; inspect individual packages and data dependencies. Source text retrieved. |
| [Teukolsky and KerrGeodesics](https://github.com/BlackHolePerturbationToolkit/Teukolsky) · ASTRO-0295 | Wolfram Language | Black-hole perturbation solutions and Kerr orbital quantities | Components of BHPToolkit; not full nonlinear spacetime evolution. Short page, metadata or access shell. |
| [FastEMRIWaveforms](https://github.com/BlackHolePerturbationToolkit/FastEMRIWaveforms) · ASTRO-0296 | Python; C++; CUDA | Extreme-mass-ratio-inspiral waveforms | Specialized waveform approximations and model domains. Source text retrieved. |
| [qnm](https://github.com/duetosymmetry/qnm) · ASTRO-0297 | Python | Black-hole quasinormal modes | Perturbative ringdown spectrum. Source text retrieved. |
| [GYOTO](https://github.com/gyoto/Gyoto) · ASTRO-0298 | C++; Python; Yorick | General-relativistic ray tracing | Metric plugins and radiative models. Source text retrieved. |
| [ipole](https://github.com/moscibrodzka/ipole) · ASTRO-0299 | C | Polarized general-relativistic radiative transfer | Post-processing of relativistic fluid simulations. Source text retrieved. |
| [grtrans](https://github.com/jadexter/grtrans) · ASTRO-0300 | Fortran; Python | General-relativistic polarized radiation transport | Metric/radiation conventions need explicit comparison. Source text retrieved. |
| [BHAC](https://bhac.science/) · ASTRO-0301 | Fortran | General-relativistic magnetohydrodynamics | Compact-object plasma simulations. Source text retrieved. |
| [HARM and iharm3D](https://github.com/AFD-Illinois/iharm3d) · ASTRO-0302 | C | General-relativistic magnetohydrodynamics | HARM lineage; scope varies by implementation. Source text retrieved. |
| [KHARMA](https://github.com/AFD-Illinois/kharma) · ASTRO-0303 | C++; Kokkos; Parthenon | Portable general-relativistic MHD | Modern implementation lineage. Source text retrieved. |
| [gevolution](https://github.com/gevolution-code/gevolution-1.2) · ASTRO-0304 | C++ | Relativistic cosmological N-body simulations | Weak-field relativistic cosmology; not interchangeable with full NR. Source text retrieved. |
| [CosmoLattice](https://cosmolattice.net/) · ASTRO-0305 | C++ | Expanding-universe classical field simulations | Early-universe nonlinear fields. Source text retrieved. |

## 15 Gravitational waves and multi messenger inference

| Software / source | Language family | Capability | Lineage, integration and evidence notes |
|---|---|---|---|
| [LALSuite](https://git.ligo.org/lscsoft/lalsuite) · ASTRO-0306 | C; Python bindings | Gravitational-wave waveforms and analysis algorithms | Umbrella for LAL, LALSimulation, inference and other components. Short page, metadata or access shell. |
| [PyCBC](https://github.com/gwastro/pycbc) · ASTRO-0307 | Python; C/CUDA components | Compact-binary searches and inference | Often uses LAL waveform models. Source text retrieved. |
| [Bilby](https://github.com/bilby-dev/bilby) · ASTRO-0308 | Python | Bayesian inference for gravitational waves and other data | Likelihood and sampler interface layer. Source text retrieved. |
| [GWpy](https://github.com/gwpy/gwpy) · ASTRO-0309 | Python | Gravitational-wave detector time-series analysis | Data quality and time-frequency analysis. Source text retrieved. |
| [GWOSC client](https://github.com/gwpy/gwosc) · ASTRO-0310 | Python | Access public gravitational-wave strain and event data | Service/data interface. Source text retrieved. |
| [GstLAL](https://git.ligo.org/lscsoft/gstlal) · ASTRO-0311 | C; Python; GStreamer | Streaming gravitational-wave analysis | Low-latency pipeline architecture. Short page, metadata or access shell. |
| [cWB](https://gwburst.gitlab.io/) · ASTRO-0312 | C++; ROOT | Coherent burst gravitational-wave searches | Waveform-independent search family. Source text retrieved. |
| [BayesWave](https://git.ligo.org/lscsoft/bayeswave) · ASTRO-0313 | C; supporting scripts | Bayesian signal and glitch reconstruction | Transient noise and signal separation. Short page, metadata or access shell. |
| [RIFT](https://git.ligo.org/rapidpe-rift/rift) · ASTRO-0314 | Python; native dependencies | Rapid gravitational-wave parameter inference | Alternative inference architecture. Short page, metadata or access shell. |
| [PESummary](https://github.com/pesummary/pesummary) · ASTRO-0315 | Python | Posterior summaries and comparison reports | Provenance and result presentation layer. Short page, metadata or access shell. |
| [ligo.skymap](https://git.ligo.org/lscsoft/ligo.skymap) · ASTRO-0316 | Python; C | Sky localization maps and multi-messenger utilities | Probability on the sphere and distance information. Short page, metadata or access shell. |
| [gwpopulation](https://github.com/ColmTalbot/gwpopulation) · ASTRO-0317 | Python; NumPy/CuPy/JAX options | Hierarchical compact-object population inference | Selection and detection efficiency are central. Source text retrieved. |
| [gwsurrogate](https://github.com/sxs-collaboration/gwsurrogate) · ASTRO-0318 | Python; native components | Numerical-relativity waveform surrogate evaluation | Training model files have separate provenance. Source text retrieved. |
| [PyRing](https://git.ligo.org/lscsoft/pyring) · ASTRO-0319 | Python | Black-hole ringdown inference | Tests of remnant and perturbation models. Short page, metadata or access shell. |

## 16 Exoplanets planetary and solar system science

| Software / source | Language family | Capability | Lineage, integration and evidence notes |
|---|---|---|---|
| [Lightkurve](https://github.com/lightkurve/lightkurve) · ASTRO-0320 | Python | Kepler and TESS light-curve processing | Time-series extraction and systematics correction. Source text retrieved. |
| [batman](https://github.com/lkreidberg/batman) · ASTRO-0321 | C; Python | Exoplanet transit light curves | Finite-source occultation and limb-darkening models. Short page, metadata or access shell. |
| [starry](https://github.com/rodluger/starry) · ASTRO-0322 | Python; C++ | Analytic occultation and surface-map light curves | Differentiable map/geometry modelling. Short page, metadata or access shell. |
| [exoplanet](https://github.com/exoplanet-dev/exoplanet) · ASTRO-0323 | Python; probabilistic-programming dependencies | Exoplanet transit and orbital inference | Inspect current backend versions before integration. Source text retrieved. |
| [juliet](https://github.com/nespinoza/juliet) · ASTRO-0324 | Python | Joint transit and radial-velocity inference | Combines existing models and samplers. Source text retrieved. |
| [RadVel](https://github.com/California-Planet-Search/radvel) · ASTRO-0325 | Python | Radial-velocity orbit fitting | Instrument offsets and stellar variability. Source text retrieved. |
| [allesfitter](https://github.com/MNGuenther/allesfitter) · ASTRO-0326 | Python | Joint exoplanet and stellar-system light-curve fitting | Integrates several numerical/model engines. Source text retrieved. |
| [PHOEBE](https://github.com/phoebe-project/phoebe2) · ASTRO-0327 | Python; C/C++ | Eclipsing-binary light curves and stellar geometry | Wilson-Devinney lineage with modern model architecture. Source text retrieved. |
| [Wilson Devinney code](https://www.astro.ufl.edu/~wilson/) · ASTRO-0328 | Fortran | Eclipsing binary light and radial-velocity modelling | Legacy physical-model reference. Retrieval unresolved. |
| [petitRADTRANS](https://gitlab.com/mauricemolli/petitRADTRANS) · ASTRO-0329 | Python; Fortran components | Exoplanet atmospheric spectra and retrieval | Opacity tables can dominate scientific data requirements. Short page, metadata or access shell. |
| [TauREx](https://taurex3.readthedocs.io/en/latest/user/) · ASTRO-0330 | Python | Exoplanet atmospheric retrieval | Modular atmospheric forward models. Source text retrieved. |
| [Exo-Transmit](https://github.com/elizakempton/Exo_Transmit) · ASTRO-0331 | C | Exoplanet transmission spectra | Opacity and chemistry tables. Source text retrieved. |
| [VPLanet](https://github.com/VirtualPlanetaryLaboratory/vplanet) · ASTRO-0332 | C; Python tooling | Coupled planetary evolution modules | Multiphysics coupling reference. Source text retrieved. |
| [Tudat](https://github.com/tudat-team/tudat) · ASTRO-0333 | C++; Python | Astrodynamics, propagation and orbit determination | The original C++ repository is archived and points to tudatpy, which now contains the C++ core; follow the current distribution. Short page, metadata or access shell. |
| [Orekit](https://www.orekit.org/) · ASTRO-0334 | Java | Orbit propagation, frames, time and estimation | Important non-Python astrodynamics ecosystem. Source text retrieved. |
| [GMAT](https://gmat.atlassian.net/wiki/) · ASTRO-0335 | C++; scripting language | Mission analysis and trajectory design | NASA mission-analysis application. Short page, metadata or access shell. |
| [OpenOrb](https://github.com/oorb/oorb) · ASTRO-0336 | Fortran; Python bindings | Small-body orbit computation and uncertainty propagation | Asteroid orbit workflows. Source text retrieved. |
| [sbpy](https://github.com/NASA-Planetary-Science/sbpy) · ASTRO-0337 | Python | Small-body ephemerides, photometry and activity analysis | Astropy ecosystem for planetary astronomy. Source text retrieved. |
| [USGS ISIS](https://github.com/DOI-USGS/ISIS3) · ASTRO-0338 | C++; Python tooling | Planetary image calibration, mapping and geometry | Third distinct ISIS name in this inventory. Source text retrieved. |

## 17 Solar plasma and space physics

| Software / source | Language family | Capability | Lineage, integration and evidence notes |
|---|---|---|---|
| [SunPy](https://github.com/sunpy/sunpy) · ASTRO-0339 | Python | Solar data discovery, maps, coordinates and time-series analysis | Major domain ecosystem analogous to Astropy. Source text retrieved. |
| [aiapy](https://github.com/LM-SAL/aiapy) · ASTRO-0340 | Python | SDO AIA calibration and analysis | Instrument-specific layer. Source text retrieved. |
| [sunkit-image](https://github.com/sunpy/sunkit-image) · ASTRO-0341 | Python | Solar image processing methods | Shared solar image analysis. Source text retrieved. |
| [sunkit-spex](https://github.com/sunpy/sunkit-spex) · ASTRO-0342 | Python | Solar X-ray spectral analysis | Relationship to legacy SolarSoft spectroscopy. Source text retrieved. |
| [PlasmaPy](https://github.com/PlasmaPy/PlasmaPy) · ASTRO-0343 | Python | Plasma parameters, particles and analysis | General plasma foundation beyond astronomy. Source text retrieved. |
| [SpacePy](https://github.com/spacepy/spacepy) · ASTRO-0344 | Python; C; Fortran | Space physics coordinates, radiation belts and analysis | External empirical models and datasets. Source text retrieved. |
| [PySPEDAS](https://github.com/spedas/pyspedas) · ASTRO-0345 | Python | Space-physics mission data and analysis | Python counterpart to IDL SPEDAS. Source text retrieved. |
| [SPEDAS](https://spedas.org/) · ASTRO-0346 | IDL | Space Physics Environment Data Analysis System | Proprietary-runtime mission ecosystem. Source text retrieved. |
| [LISIRD solar model and data archive](https://lasp.colorado.edu/lisird/) · ASTRO-0347 | Mixed; models and data | Solar irradiance and space-weather model products | Data/model service entry; not an independent numerical software package. Short page, metadata or access shell. |
| [BOUT++](https://github.com/boutproject/BOUT-dev) · ASTRO-0348 | C++ | Plasma-fluid PDE simulations | General physics and fusion overlap. Source text retrieved. |
| [Smilei](https://github.com/SmileiPIC/Smilei) · ASTRO-0349 | C++; Python | Particle-in-cell plasma simulation | Kinetic plasma rather than fluid MHD. Source text retrieved. |
| [WarpX](https://github.com/BLAST-WarpX/warpx) · ASTRO-0350 | C++; Python; AMReX | Accelerated particle-in-cell simulations | Electromagnetic particles and fields. Source text retrieved. |
| [EPOCH](https://github.com/epochpic/epoch) · ASTRO-0351 | Fortran | Particle-in-cell plasma simulations | Alternative kinetic implementation. Source text retrieved. |
| [Vlasiator](https://github.com/fmihpc/vlasiator) · ASTRO-0352 | C++ | Hybrid-Vlasov space-plasma simulations | Distribution-function physics and large-scale simulation. Source text retrieved. |

## 18 Julia scientific and astronomy ecosystems

| Software / source | Language family | Capability | Lineage, integration and evidence notes |
|---|---|---|---|
| [JuliaAstro ecosystem](https://juliaastro.org/) · ASTRO-0353 | Julia | Astronomical types, FITS, coordinates, cosmology, time and images | Umbrella; individual packages have independent status. Short page, metadata or access shell. |
| [AstroImages.jl](https://github.com/JuliaAstro/AstroImages.jl) · ASTRO-0354 | Julia | Astronomical image arrays, FITS and WCS-aware display | Image/data model reference. Source text retrieved. |
| [FITSIO.jl](https://github.com/JuliaAstro/FITSIO.jl) · ASTRO-0355 | Julia; CFITSIO | FITS reading and writing | Wrapper, not an independent FITS engine. Source text retrieved. |
| [SkyCoords.jl](https://github.com/JuliaAstro/SkyCoords.jl) · ASTRO-0356 | Julia | Astronomical coordinate systems and transformations | Validate frames and precision scope. Source text retrieved. |
| [AstroTime.jl](https://github.com/JuliaAstro/AstroTime.jl) · ASTRO-0357 | Julia | Astronomical time scales and epochs | Precision time representations. Source text retrieved. |
| [Cosmology.jl](https://github.com/JuliaAstro/Cosmology.jl) · ASTRO-0358 | Julia | Cosmological background quantities and distances | Do not confuse similarly named repositories or assume full Boltzmann support. Short page, metadata or access shell. |
| [UnitfulAstro.jl](https://github.com/JuliaAstro/UnitfulAstro.jl) · ASTRO-0359 | Julia | Astronomical units for Unitful quantities | Shared units extension. Source text retrieved. |
| [PSFModels.jl](https://github.com/JuliaAstro/PSFModels.jl) · ASTRO-0360 | Julia | Point-spread-function models and fitting | Imaging primitive. Source text retrieved. |
| [Photometry.jl](https://github.com/JuliaAstro/Photometry.jl) · ASTRO-0361 | Julia | Astronomical photometry methods | Scope and maintenance require package-level review. Source text retrieved. |
| [Transits.jl](https://github.com/JuliaAstro/Transits.jl) · ASTRO-0362 | Julia | Exoplanet transit and occultation calculations | Geometry and differentiability reference. Source text retrieved. |
| [Octofitter.jl](https://github.com/sefffal/Octofitter.jl) · ASTRO-0363 | Julia | Joint orbital fitting and inference | Direct imaging, astrometry and related orbital observations. Source text retrieved. |
| [VLBISkyModels.jl](https://github.com/EHTJulia/VLBISkyModels.jl) · ASTRO-0364 | Julia | Sky intensity models and interferometric observables | EHT and VLBI forward modelling. Short page, metadata or access shell. |
| [Comrade.jl](https://github.com/ptiede/Comrade.jl) · ASTRO-0365 | Julia | Bayesian VLBI imaging and modelling | Probabilistic radio-imaging ecosystem. Source text retrieved. |
| [Gradus.jl](https://github.com/astro-group-bristol/Gradus.jl) · ASTRO-0366 | Julia | General-relativistic ray tracing and spectral models | GitHub repository states that development moved to Codeberg; follow its link for current source. Short page, metadata or access shell. |
| [Krang.jl](https://github.com/dchang10/Krang.jl) · ASTRO-0367 | Julia | Kerr ray tracing and black-hole image models | Specialist relativistic modelling. Source text retrieved. |
| [Skylight.jl](https://github.com/joaquinpelle/Skylight.jl) · ASTRO-0368 | Julia | Relativistic ray tracing and radiation transfer | Primary repository identity confirmed for general-relativistic ray tracing and radiative transfer. Source text retrieved. |
| [Bolt.jl](https://github.com/xzackli/Bolt.jl) · ASTRO-0369 | Julia | Differentiable cosmological Boltzmann integration | Research implementation; scope differs from CAMB/CLASS. Source text retrieved. |
| [SymBoltz.jl](https://github.com/hersle/SymBoltz.jl) · ASTRO-0370 | Julia | Symbolic-numeric Einstein-Boltzmann solver | Direct model equations and automatic differentiation. Source text retrieved. |

## 19 General numerical inference and physics foundations

| Software / source | Language family | Capability | Lineage, integration and evidence notes |
|---|---|---|---|
| [NumPy and SciPy](https://scipy.org/) · ASTRO-0371 | Python; C; C++; Fortran | Arrays, linear algebra, integration, optimization, statistics and special functions | Foundational capability families, not astronomy-specific. Source text retrieved. |
| [SymPy](https://github.com/sympy/sympy) · ASTRO-0372 | Python | Symbolic mathematics and code generation | Algebraic foundation for physics tooling. Source text retrieved. |
| [JAX](https://github.com/jax-ml/jax) · ASTRO-0373 | Python; compiled XLA backend | Automatic differentiation, array transformations and accelerators | Execution/model layer rather than a science model. Source text retrieved. |
| [Diffrax](https://github.com/patrick-kidger/diffrax) · ASTRO-0374 | Python; JAX | Differentiable ordinary and stochastic differential equations | Useful solver architecture. Source text retrieved. |
| [SciML DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) · ASTRO-0375 | Julia | ODE, SDE, DAE, delay and related numerical solvers | Large solver ecosystem with individually versioned components. Source text retrieved. |
| [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) · ASTRO-0376 | Julia | Symbolic system modelling and numerical code generation | Equation-to-solver architecture reference. Source text retrieved. |
| [Dedalus](https://github.com/DedalusProject/dedalus) · ASTRO-0377 | Python; compiled numerical libraries | Spectral PDE solvers and physical simulations | Equations and boundary conditions specified at high level. Source text retrieved. |
| [FEniCSx](https://fenicsproject.org/) · ASTRO-0378 | C++; Python | Finite-element PDE discretization and solution | Code generation and variational formulations. Source text retrieved. |
| [Firedrake](https://www.firedrakeproject.org/) · ASTRO-0379 | Python; generated/native kernels | Finite-element PDE solution | Alternative variational solver architecture. Source text retrieved. |
| [deal.II](https://www.dealii.org/) · ASTRO-0380 | C++ | Finite-element numerical methods | General multiphysics infrastructure. Source text retrieved. |
| [MFEM](https://mfem.org/) · ASTRO-0381 | C++ | Modular finite-element methods | High-order methods and accelerators. Source text retrieved. |
| [PETSc](https://petsc.org/) · ASTRO-0382 | C; Fortran; Python interfaces | Parallel linear and nonlinear systems and time stepping | Numerical engine used by many domain applications. Source text retrieved. |
| [SUNDIALS](https://sundials.readthedocs.io/) · ASTRO-0383 | C; C++; Fortran interfaces | ODE, DAE and nonlinear solvers with sensitivities | Stiff systems and sensitivity analysis. Source text retrieved. |
| [AMReX](https://github.com/AMReX-Codes/amrex) · ASTRO-0384 | C++; Fortran; Python interfaces | Block-structured adaptive mesh refinement | Infrastructure for astrophysics and plasma codes. Source text retrieved. |
| [Kokkos](https://github.com/kokkos/kokkos) · ASTRO-0385 | C++ | Portable parallel numerical kernels | Execution portability reference. Source text retrieved. |
| [emcee](https://github.com/dfm/emcee) · ASTRO-0386 | Python | Ensemble Markov chain Monte Carlo | Sampler, not a likelihood or physical model. Source text retrieved. |
| [dynesty](https://github.com/joshspeagle/dynesty) · ASTRO-0387 | Python | Nested sampling and evidence estimation | Posterior and evidence workflow. Source text retrieved. |
| [UltraNest](https://github.com/JohannesBuchner/UltraNest) · ASTRO-0388 | Python; Cython | Nested sampling and diagnostics | Alternative sampler and validation tooling. Source text retrieved. |
| [MultiNest and PyMultiNest](https://github.com/farhanferoz/MultiNest) · ASTRO-0389 | Fortran; Python | Nested sampling and multimodal evidence calculation | Legacy native sampler with wrappers. Source text retrieved. |
| [PolyChord](https://github.com/PolyChord/PolyChordLite) · ASTRO-0390 | Fortran; C++; Python | High-dimensional nested sampling | Sampler implementation family. Source text retrieved. |
| [NumPyro](https://github.com/pyro-ppl/numpyro) · ASTRO-0391 | Python; JAX | Probabilistic models and gradient-based Bayesian inference | Autodifferentiation and model compilation dependency. Source text retrieved. |
| [PyMC](https://github.com/pymc-devs/pymc) · ASTRO-0392 | Python; PyTensor | Probabilistic programming and Bayesian inference | Alternative modelling ecosystem. Source text retrieved. |
| [Stan](https://mc-stan.org/) · ASTRO-0393 | Stan language; C++ | Compiled probabilistic models and inference | Independent scientific modelling language. Source text retrieved. |
| [iminuit and Minuit2](https://github.com/scikit-hep/iminuit) · ASTRO-0394 | Python; C++ | Numerical minimization and uncertainty estimation | Optimizer diagnostics and likelihood assumptions must remain visible. Source text retrieved. |
| [celerite and tinygp](https://github.com/dfm/celerite2) · ASTRO-0395 | Python; C++; JAX variants | Gaussian-process time-series and covariance modelling | Family entry; separate algorithms and supported kernels. Source text retrieved. |

## 20 General physics quantum and particle software

| Software / source | Language family | Capability | Lineage, integration and evidence notes |
|---|---|---|---|
| [Geant4](https://geant4.web.cern.ch/) · ASTRO-0396 | C++ | Particle transport through matter and detector simulation | Relevant to instrument response and particle astrophysics. Source text retrieved. |
| [ROOT and RooFit](https://root.cern/) · ASTRO-0397 | C++; Python | Particle-physics data processing, statistics and likelihood models | Large integrated analysis ecosystem. Source text retrieved. |
| [CRPropa](https://github.com/CRPropa/CRPropa3) · ASTRO-0398 | C++; Python | High-energy cosmic-ray propagation | Astroparticle and magnetic-field modelling. Source text retrieved. |
| [GALPROP](https://ascl.net/1010.028) · ASTRO-0399 | C++; other components | Galactic cosmic-ray transport and diffuse emission | Data, source populations and propagation assumptions matter. Source text retrieved. [Original nominated site](https://galprop.stanford.edu/); source link uses ASCL. |
| [CORSIKA](https://www.iap.kit.edu/corsika/) · ASTRO-0400 | Fortran; C++ successor | Extensive air-shower simulation | Distinguish CORSIKA 7 and CORSIKA 8 implementations. Source text retrieved. |
| [QuTiP](https://github.com/qutip/qutip) · ASTRO-0401 | Python; Cython | Open quantum-system dynamics | General quantum physics capability. Source text retrieved. |
| [QuantumOptics.jl](https://github.com/qojulia/QuantumOptics.jl) · ASTRO-0402 | Julia | Quantum optics and open-system simulation | Julia counterpart in a related capability area. Source text retrieved. |
| [Kwant](https://kwant-project.org/) · ASTRO-0403 | Python; C/C++ | Quantum transport in tight-binding systems | General condensed-matter physics. Source text retrieved. |
| [Meep](https://github.com/NanoComp/meep) · ASTRO-0404 | C++; Python; Scheme | Finite-difference time-domain electromagnetics | Useful for optics and instrument modelling. Source text retrieved. |
| [MPB](https://github.com/NanoComp/mpb) · ASTRO-0405 | C; Scheme; Python interfaces | Photonic band-structure calculation | Frequency-domain electromagnetic eigenproblems. Source text retrieved. |
| [LAMMPS](https://github.com/lammps/lammps) · ASTRO-0406 | C++ | Classical molecular dynamics and particle simulation | General physics scope extension. Source text retrieved. |
| [GROMACS](https://gitlab.com/gromacs/gromacs) · ASTRO-0407 | C++; C | Molecular dynamics | Primarily molecular/biophysical domain. Short page, metadata or access shell. |
| [OpenMM](https://github.com/openmm/openmm) · ASTRO-0408 | C++; Python | Molecular simulation and programmable particle forces | Reusable GPU simulation architecture. Source text retrieved. |
| [Quantum ESPRESSO](https://www.quantum-espresso.org/) · ASTRO-0409 | Fortran; C | Electronic structure and density-functional theory | Materials physics, beyond core astronomy. Source text retrieved. |
| [GPAW and ASE](https://gpaw.readthedocs.io/) · ASTRO-0410 | Python; C; numerical libraries | Electronic structure and atomistic simulation workflows | Related but distinct calculator and environment. Source text retrieved. |
| [FeynCalc](https://github.com/FeynCalc/feyncalc) · ASTRO-0411 | Wolfram Language | Symbolic quantum-field-theory calculations | Open package with proprietary runtime. Source text retrieved. |
| [FeynRules](https://feynrules.irmp.ucl.ac.be/) · ASTRO-0412 | Wolfram Language | Lagrangian models and particle-physics interaction rules | Model-generation capability. Source text retrieved. |

## 21 Visualization archives and operations

| Software / source | Language family | Capability | Lineage, integration and evidence notes |
|---|---|---|---|
| [SAOImage DS9](https://ds9.si.edu/site/Home.html) · ASTRO-0413 | C++; Tcl/Tk | FITS image display, regions and interactive analysis | Widely integrated viewer; do not reduce scope to plotting. Short page, metadata or access shell. |
| [TOPCAT and STILTS](https://www.star.bris.ac.uk/~mbt/topcat/) · ASTRO-0414 | Java | Catalogue exploration, crossmatching and table processing | GUI and command-line pair on STIL library. Source text retrieved. |
| [Aladin Desktop and Aladin Lite](https://aladin.cds.unistra.fr/) · ASTRO-0415 | Java; JavaScript; Rust/WebAssembly | Sky atlas, catalogue overlays and survey visualization | Desktop and browser implementations differ. Source text retrieved. |
| [Glue](https://github.com/glue-viz/glue) · ASTRO-0416 | Python | Linked multidimensional scientific data exploration | Cross-dataset selection and visualization. Source text retrieved. |
| [Ginga](https://github.com/ejeschke/ginga) · ASTRO-0417 | Python | Astronomical image viewer and plugin framework | Embedded or standalone user interface. Source text retrieved. |
| [Jdaviz](https://github.com/spacetelescope/jdaviz) · ASTRO-0418 | Python; web widgets | Interactive imaging, spectral and cube analysis | STScI analysis application ecosystem. Source text retrieved. |
| [CARTA](https://cartavis.org/) · ASTRO-0419 | C++; TypeScript | Interactive large radio-image and cube visualization | Client/server architecture for large data. Source text retrieved. |
| [VisIt](https://visit-dav.github.io/visit-website/) · ASTRO-0420 | C++; Python | Large scientific-simulation visualization | General simulation visualization. Source text retrieved. |
| [ParaView and VTK](https://www.paraview.org/) · ASTRO-0421 | C++; Python | Scientific visualization, meshes and volume rendering | General infrastructure and application pair. Source text retrieved. |
| [OpenSpace](https://www.openspaceproject.com/) · ASTRO-0422 | Mixed | Interactive visualization of astronomical datasets and space missions | Scientific visualization application; distinct from yt and WorldWide Telescope. Source text retrieved. |
| [Astroplan](https://github.com/astropy/astroplan) · ASTRO-0423 | Python | Observation planning, observability and scheduling | Constraints depend on sites, ephemerides and time conventions. Source text retrieved. |
| [INDI](https://github.com/indilib/indi) · ASTRO-0424 | C++; XML protocol | Astronomical instrument control | Hardware/operations boundary. Source text retrieved. |
| [ASCOM and Alpaca](https://ascom-standards.org/) · ASTRO-0425 | Windows/.NET; network protocol | Telescope and instrument interoperability | Protocol compatibility matters more than implementation language. Source text retrieved. |
| [RTS2](https://github.com/RTS2/rts2) · ASTRO-0426 | C++; Python | Robotic telescope control and scheduling | Operations rather than science inference. Source text retrieved. |

## 22 Existing Rust scientific building blocks

| Software / source | Language family | Capability | Lineage, integration and evidence notes |
|---|---|---|---|
| [ndarray](https://github.com/rust-ndarray/ndarray) · ASTRO-0427 | Rust | N-dimensional arrays and numerical operations | Candidate foundation; benchmark needed, not chosen automatically. Source text retrieved. |
| [nalgebra](https://github.com/dimforge/nalgebra) · ASTRO-0428 | Rust | Linear algebra and geometric transformations | Candidate foundation. Short page, metadata or access shell. |
| [faer](https://github.com/sarah-ek/faer-rs) · ASTRO-0429 | Rust | Dense and sparse linear algebra | GitHub landing page states that development migrated to Codeberg and GitHub remains a mirror. Short page, metadata or access shell. |
| [uom](https://github.com/iliekturtles/uom) · ASTRO-0430 | Rust | Type-safe units of measurement | Astronomical equivalencies, logarithmic units and runtime units need separate design. Source text retrieved. |
| [hifitime](https://github.com/nyx-space/hifitime) · ASTRO-0431 | Rust | High-precision epochs, durations and time scales | Potential astronomy/astrodynamics foundation. Source text retrieved. |
| [ANISE](https://github.com/nyx-space/anise) · ASTRO-0432 | Rust; Python bindings | Spacecraft geometry, ephemerides and reference frames | Existing Rust astrodynamics capability. Source text retrieved. |
| [Nyx](https://github.com/nyx-space/nyx) · ASTRO-0433 | Rust | Astrodynamics and orbit determination | Existing domain engine. Source text retrieved. |
| [rust-fitsio](https://github.com/simonrw/rust-fitsio) · ASTRO-0434 | Rust; CFITSIO binding | FITS data access | FFI binding; not pure Rust FITS implementation. Source text retrieved. |
| [CDS HEALPix Rust](https://github.com/cds-astro/cds-healpix-rust) · ASTRO-0435 | Rust | HEALPix geometry and spatial indexing | Existing astronomy kernel. Source text retrieved. |
| [CDS MOC Rust](https://github.com/cds-astro/cds-moc-rust) · ASTRO-0436 | Rust | Multi-order spatial and space-time coverage | Core related to MOCPy. Source text retrieved. |
| [RustFFT](https://github.com/ejmahler/RustFFT) · ASTRO-0437 | Rust | Fast Fourier transforms | Shared signal-processing kernel. Source text retrieved. |
| [argmin](https://github.com/argmin-rs/argmin) · ASTRO-0438 | Rust | Numerical optimization | Candidate optimizer framework. Source text retrieved. |
| [diffsol](https://github.com/martinjrobins/diffsol) · ASTRO-0439 | Rust | Differential-equation solvers | Assess stiff solvers and sensitivity support against required physics. Source text retrieved. |
| [Burn](https://github.com/tracel-ai/burn) · ASTRO-0440 | Rust | Tensor computation and machine learning | Potential accelerator/model backend, not a physics package. Source text retrieved. |

## 23 Optical instruments adaptive optics and interferometry

| Software / source | Language family | Capability | Lineage, integration and evidence notes |
|---|---|---|---|
| [HCIPy](https://github.com/ehpor/hcipy) · ASTRO-0441 | Python | Wave optics, atmospheric turbulence, adaptive optics and coronagraphs | End-to-end high-contrast instrument models. Source text retrieved. |
| [POPPY](https://github.com/spacetelescope/poppy) · ASTRO-0442 | Python | Physical-optics propagation and point-spread functions | Reusable optical propagation engine. Source text retrieved. |
| [STPSF and WebbPSF](https://github.com/spacetelescope/stpsf) · ASTRO-0443 | Python | Space-telescope instrument PSF simulation | STPSF is the successor name; optical calibration data are required. Source text retrieved. |
| [dLux](https://github.com/LouisDesdoigts/dLux) · ASTRO-0444 | Python; JAX | Differentiable optical-system modelling | Wavefront and telescope inference. Source text retrieved. |
| [COMPASS](https://github.com/ANR-COMPASS/shesha) · ASTRO-0445 | C++; CUDA; Python | Adaptive-optics system simulation | This GitHub repository is marked obsolete and points to an Observatoire de Paris GitLab successor; preserve lineage when acquiring source. Short page, metadata or access shell. |
| [Soapy](https://github.com/AOtools/soapy) · ASTRO-0446 | Python | Adaptive-optics simulation | Alternative modular AO implementation. Source text retrieved. |
| [AOtools](https://github.com/AOtools/aotools) · ASTRO-0447 | Python | Adaptive-optics numerical utilities | Atmospheric screens, wavefronts and related primitives. Source text retrieved. |
| [VIP](https://github.com/vortex-exoplanet/VIP) · ASTRO-0448 | Python | High-contrast imaging and faint-companion analysis | Post-processing, detection and characterization. Source text retrieved. |
| [pyKLIP](https://bitbucket.org/pyKLIP/pyklip/) · ASTRO-0449 | Python | PSF subtraction and forward modelling for direct imaging | Karhunen-Loeve image-projection algorithms. Short page, metadata or access shell. |
| [PMOIRED](https://github.com/amerand/PMOIRED) · ASTRO-0450 | Python | Optical-interferometric data modelling | OIFITS observables and parametric models. Source text retrieved. |
| [MiRA](https://github.com/emmt/MiRA) · ASTRO-0451 | Yorick; C components | Image reconstruction from optical interferometry | Optimization under incomplete Fourier information. Source text retrieved. |
| [SQUEEZE](https://github.com/fabienbaron/squeeze) · ASTRO-0452 | C++ | Optical-interferometric image reconstruction | Alternative reconstruction and regularization strategies. Source text retrieved. |
| [LITpro](https://www.jmmc.fr/english/tools/proposal-preparation/litpro/) · ASTRO-0453 | Java client; numerical services | Parametric optical-interferometry model fitting | Distributed service access and model implementation differ. Short page, metadata or access shell. |
| [OIFITS and OIFITSlib](https://www.jmmc.fr/oifits/) · ASTRO-0454 | Standard; C implementations | Exchange optical-interferometric observations and uncertainties | Data-standard entry; implementations need separate audit. Short page, metadata or access shell. |

## 24 Additional legacy R and domain specialist software

| Software / source | Language family | Capability | Lineage, integration and evidence notes |
|---|---|---|---|
| [ProFound](https://github.com/asgr/ProFound) · ASTRO-0455 | R; C++ | Source detection, segmentation and photometry | Important astronomy ecosystem outside Python. Source text retrieved. |
| [ProFit](https://github.com/ICRAR/ProFit) · ASTRO-0456 | R; C++ | Bayesian galaxy image modelling | Galaxy fitting and structural inference. Source text retrieved. |
| [ProSpect](https://github.com/asgr/ProSpect) · ASTRO-0457 | R; C++ | Galaxy spectral-energy-distribution modelling | Stellar populations, dust and fitting. Source text retrieved. |
| [celestial](https://github.com/asgr/celestial) · ASTRO-0458 | R | Astronomical coordinate and cosmology utilities | Alternative-language foundation. Source text retrieved. |
| [KAPPA](https://starlink.eao.hawaii.edu/docs/sun95.htx/sun95.html) · ASTRO-0459 | Fortran; C | General image and NDF manipulation and analysis | Starlink component; useful mature numerical task catalogue. Source text retrieved. |
| [SMURF](https://starlink.eao.hawaii.edu/docs/sun258.htx/sun258.html) · ASTRO-0460 | C; Fortran | Submillimetre time-stream mapmaking and reduction | JCMT detector and mapmaking knowledge. Source text retrieved. |
| [FIGARO and DIPSO](https://starlink.eao.hawaii.edu/docs/sun86.htx/sun86.html) · ASTRO-0461 | Fortran; C components | Spectral reduction, fitting and interactive analysis | Legacy Starlink spectral-analysis lineage; FIGARO documentation cited. Source text retrieved. |
| [CCDPACK](https://starlink.eao.hawaii.edu/docs/sun139.htx/sun139.html) · ASTRO-0462 | Fortran; C | CCD calibration, registration and combination | Legacy pipeline algorithms and parameter semantics. Source text retrieved. |
| [CUPID](https://starlink.eao.hawaii.edu/docs/sun255.htx/sun255.html) · ASTRO-0463 | C; Fortran | Identify and characterize emission clumps in images and cubes | Multiple detection algorithms with different assumptions. Source text retrieved. |
| [PGPLOT](https://sites.astro.caltech.edu/~tjp/pgplot/) · ASTRO-0464 | Fortran; C bindings | Scientific plotting used by legacy astronomy software | Preserve diagnostic plots when recovering old codes. Source text retrieved. |
| [SM SuperMongo](https://www.astro.princeton.edu/~rhl/sm/) · ASTRO-0465 | C; command language | Interactive scientific plotting and scripting | Historical astronomy analysis environment. Source text retrieved. |
| [MPFIT](https://ascl.net/1208.019) · ASTRO-0466 | IDL; C ports | Nonlinear least-squares fitting based on MINPACK | Shared optimizer lineage behind many older analyses. Source text retrieved. [Original nominated site](https://pages.physics.wisc.edu/~craigm/idl/idl.html); source link uses ASCL. |
| [idlutils](https://github.com/sdss/idlutils) · ASTRO-0467 | IDL; C support | SDSS analysis utilities, spectra, masks and fitting | Utilities distinct from the idlspec2d pipeline entry. Source text retrieved. |
| [VARTOOLS](https://www.astro.princeton.edu/~jhartman/vartools.html) · ASTRO-0468 | C | Astronomical light-curve analysis | Broad time-series task collection. Short page, metadata or access shell. |
| [VaST](https://github.com/kirxkirx/vast) · ASTRO-0469 | C; shell; other components | Variability search from imaging time series | Legacy-friendly photometric workflow. Source text retrieved. |

## 25 Additional simulation and emerging modelling capabilities

| Software / source | Language family | Capability | Lineage, integration and evidence notes |
|---|---|---|---|
| [GRTeclyn](https://github.com/GRTLCollaboration/GRTeclyn) · ASTRO-0470 | C++; AMReX | Accelerated numerical relativity | Successor to GRChombo; preserve historical code as a comparison target. Source text retrieved. |
| [Castro](https://github.com/AMReX-Astro/Castro) · ASTRO-0471 | C++; Fortran; AMReX | Compressible astrophysical radiation hydrodynamics | Stellar explosions and microphysics coupling. Source text retrieved. |
| [MAESTROeX](https://github.com/AMReX-Astro/MAESTROeX) · ASTRO-0472 | C++; Fortran; AMReX | Low-Mach-number stellar hydrodynamics | Filters acoustic dynamics; different validity regime from compressible codes. Source text retrieved. |
| [MUSIC](https://www.oca.eu/en/olivier-hahn/1914-music) · ASTRO-0473 | C++ | Multiscale cosmological initial conditions | Initial-condition generation is separate from evolution. Source text retrieved. |
| [ROCKSTAR](https://bitbucket.org/gfcstanford/rockstar/) · ASTRO-0474 | C | Phase-space halo identification and merger analysis | Halo definition changes downstream catalogues. Short page, metadata or access shell. |
| [AHF](https://ascl.net/1102.009) · ASTRO-0475 | C | Adaptive Mesh Investigations of Galaxy Assembly halo finder | Alternative halo definitions and substructure. Source text retrieved. [Original nominated site](https://popia.ft.uam.es/AHF/); source link uses ASCL. |
| [VELOCIraptor](https://github.com/pelahi/VELOCIraptor-STF) · ASTRO-0476 | C++ | Phase-space structure and halo identification | Particle catalogues and merger workflows. Source text retrieved. |
| [NIFTy](https://github.com/NIFTy-PPL/NIFTy) · ASTRO-0477 | Python; JAX components | Information field theory and Bayesian field reconstruction | Astronomy and general inverse problems. Source text retrieved. |
| [JaxPM](https://github.com/DifferentiableUniverseInitiative/JaxPM) · ASTRO-0478 | Python; JAX | Differentiable particle-mesh cosmological simulations | Gradient-based initial-condition and cosmology inference. Source text retrieved. |
| [JAX-GalSim](https://github.com/GalSim-developers/JAX-GalSim) · ASTRO-0479 | Python; JAX | Differentiable astronomical image simulation | Project documentation explicitly labels early development and directs scientific use to reference GalSim. Source text retrieved. |
| [AstroML](https://github.com/astroML/astroML) · ASTRO-0480 | Python | Astronomical statistics, machine learning and examples | Useful algorithm and example collection. Source text retrieved. |
| [sbi](https://github.com/sbi-dev/sbi) · ASTRO-0481 | Python; PyTorch | Simulation-based Bayesian inference | Validation must include simulator support and coverage. Source text retrieved. |
| [nuSQuIDS](https://github.com/arguelles/nuSQuIDS) · ASTRO-0482 | C++; Python bindings | Neutrino evolution and oscillations in matter | Astroparticle and general physics capability. Source text retrieved. |
| [SNEWPY](https://github.com/SNEWS2/snewpy) · ASTRO-0483 | Python | Supernova-neutrino models, flavor transformations and detector predictions | Links source simulations to SNOwGLoBES detector calculations. Source text retrieved. |
| [SNOwGLoBES](https://github.com/SNOwGLoBES/snowglobes) · ASTRO-0484 | C; data and scripts | Supernova-neutrino detector event-rate calculations | Cross sections and detector tables are essential inputs. Source text retrieved. |
| [IceTray](https://github.com/icecube/icetray-public) · ASTRO-0485 | C++; Python | IceCube event processing framework | Public framework does not imply every collaboration module is public. Source text retrieved. |
| [Bifrost](https://ita-solar.github.io/Bifrost/Start_here/) · ASTRO-0486 | Fortran; IDL/Python analysis | Radiation-MHD solar atmosphere simulations | Official documentation is public; source-repository access and distribution require separate confirmation. Source text retrieved. |
| [MURaM](https://www2.mps.mpg.de/projects/solar-mhd/muram_site/) · ASTRO-0487 | C++; Fortran lineage | Solar and stellar radiative magnetohydrodynamics | Distribution conditions and branches require checking. Source text retrieved. |
| [RH](https://github.com/ITA-Solar/rh) · ASTRO-0488 | C; Python tools | Non-LTE radiative transfer and solar spectral modelling | Atomic models and formal-solver assumptions. Source text retrieved. |
| [STiC](https://github.com/jaimedelacruz/stic) · ASTRO-0489 | C++; Python | Non-LTE spectropolarimetric inversion | Infers atmospheric properties from polarized spectra. Source text retrieved. |
| [NICOLE](https://github.com/hsocasnavarro/NICOLE) · ASTRO-0490 | Fortran; Python tools | Non-LTE inversion and polarized spectral synthesis | Alternative solar-atmosphere inversion lineage. Source text retrieved. |
| [AstroLib.jl](https://github.com/JuliaAstro/AstroLib.jl) · ASTRO-0491 | Julia | Astronomical and astrophysical utility algorithms | Explicit alternative-language analogue of historical utilities. Source text retrieved. |
| [Trixi.jl](https://github.com/trixi-framework/Trixi.jl) · ASTRO-0492 | Julia | High-order hyperbolic PDE solvers with adaptive meshes | General conservation-law physics engine. Source text retrieved. |
| [GLoW](https://github.com/glow-astro/GLoW) · ASTRO-0493 | Python; C; Cython | Wave-optics gravitational lensing and amplification factors | Diffraction and interference extend geometric-optics lensing; connects directly to gravitational-wave inference. Source text retrieved. |

## 26 Shared data storage execution and numerical infrastructure

| Software / source | Language family | Capability | Lineage, integration and evidence notes |
|---|---|---|---|
| [HDF5](https://github.com/HDFGroup/hdf5) · ASTRO-0494 | C; Fortran; C++ interfaces | Hierarchical scientific arrays, metadata, chunking and compression | Shared storage infrastructure beneath many simulation and observational workflows. Source text retrieved. |
| [h5py](https://github.com/h5py/h5py) · ASTRO-0495 | Python; Cython; HDF5 binding | NumPy-oriented interface to HDF5 datasets | Wrapper around HDF5; not an independent storage engine. Source text retrieved. |
| [netCDF and netCDF4-Python](https://github.com/Unidata/netcdf4-python) · ASTRO-0496 | C; Fortran; Python | Self-describing multidimensional scientific datasets | Some netCDF variants use HDF5; preserve schema and dimension semantics. Source text retrieved. |
| [Zarr](https://github.com/zarr-developers/zarr-python) · ASTRO-0497 | Specification; Python and other implementations | Chunked compressed N-dimensional array storage | Relevant to object storage and large cubes; not a scientific observation model by itself. Source text retrieved. |
| [xarray](https://github.com/pydata/xarray) · ASTRO-0498 | Python | Named dimensions, coordinates and labelled multidimensional arrays | Useful shared data abstraction; can use Dask and several storage backends. Source text retrieved. |
| [Dask](https://github.com/dask/dask) · ASTRO-0499 | Python | Task graphs, parallel arrays and distributed data processing | Execution and chunk scheduling overlap with bespoke scientific pipeline frameworks. Short page, metadata or access shell. |
| [Apache Arrow and Parquet](https://arrow.apache.org/) · ASTRO-0500 | C++; Rust; Java; Python; format specifications | Columnar memory interchange and analytical table storage | Data interchange is distinct from astronomy metadata semantics; Arrow and Parquet serve different roles. Source text retrieved. |
| [Polars](https://github.com/pola-rs/polars) · ASTRO-0501 | Rust; Python | Columnar table processing, lazy queries and joins | Existing Rust engine relevant to catalogue processing; astronomy types still need integration. Source text retrieved. |
| [pandas](https://github.com/pandas-dev/pandas) · ASTRO-0502 | Python; Cython; native dependencies | In-memory table cleaning, joins, grouping and statistics | Overlaps with Astropy tables at generic operations, with different units and metadata behaviour. Source text retrieved. |
| [Vaex](https://github.com/vaexio/vaex) · ASTRO-0503 | Python; C++ | Large tabular data exploration and aggregation | Catalogue-scale processing; evaluate data model and execution differences. Source text retrieved. |
| [DuckDB](https://github.com/duckdb/duckdb) · ASTRO-0504 | C++; Python and other interfaces | Embedded analytical SQL over files and tables | Candidate common catalogue-query service, not an astronomy-specific package. Source text retrieved. |
| [fsspec](https://github.com/fsspec/filesystem_spec) · ASTRO-0505 | Python | Common interfaces to local and remote filesystems | Centralize access, caching and range requests; filesystem access is separate from scientific format decoding. Source text retrieved. |
| [MPI and mpi4py](https://github.com/mpi4py/mpi4py) · ASTRO-0506 | Standard; C/Fortran implementations; Python bindings | Distributed communication and collective numerical execution | Bindings and MPI implementations are separate layers, not independent physics engines. Source text retrieved. |
| [FFTW](https://www.fftw.org/) · ASTRO-0507 | C; Fortran interfaces | Fast Fourier transforms and reusable execution plans | Shared numerical kernel across signal, map and simulation workflows. Source text retrieved. |
| [ducc](https://github.com/mreineck/ducc) · ASTRO-0508 | C++; Python bindings | FFT, spherical-harmonic transforms, gridding and related numerical primitives | Includes development lineage related to libsharp; useful radio and CMB overlap. Source text retrieved. |
| [libsharp](https://github.com/Libsharp/libsharp) · ASTRO-0509 | C | Spherical-harmonic transforms for scientific maps | Shared transform engine; compare normalization and spin conventions with alternatives. Source text retrieved. |
| [SHTns](https://bitbucket.org/nschaeff/shtns/) · ASTRO-0510 | C; Python interfaces | Spherical-harmonic transforms | Alternative implementation and execution architecture. Short page, metadata or access shell. |
| [BLAS and LAPACK](https://www.netlib.org/lapack/) · ASTRO-0511 | Fortran; C interfaces; multiple implementations | Dense linear algebra, matrix factorizations and eigenproblems | Specifications, reference implementation and optimized implementations must be distinguished. Source text retrieved. |
| [OpenBLAS](https://github.com/OpenMathLib/OpenBLAS) · ASTRO-0512 | C; Fortran; assembly | Optimized BLAS and LAPACK-related numerical kernels | Multiple frontend packages may call the same underlying library. Source text retrieved. |
| [SuiteSparse](https://github.com/DrTimothyAldenDavis/SuiteSparse) · ASTRO-0513 | C; C++; MATLAB interfaces | Sparse matrix factorization and graph-based numerical kernels | Shared core for inverse problems and sparse solvers. Source text retrieved. |
| [GNU Scientific Library](https://www.gnu.org/software/gsl/) · ASTRO-0514 | C | Special functions, integration, interpolation, random numbers, optimization and ODE solvers | Broad overlapping numerical foundation used by many compiled scientific packages. Source text retrieved. |
| [MINPACK](https://www.netlib.org/minpack/) · ASTRO-0515 | Fortran; multiple ports and wrappers | Nonlinear least squares and nonlinear equations | Legacy algorithm lineage behind several fitting interfaces; ports are not independent methods. Source text retrieved. |
| [ODEPACK](https://www.netlib.org/odepack/) · ASTRO-0516 | Fortran | Ordinary differential equations including stiff and nonstiff systems | Legacy LSODE/LSODA family; modern wrapper availability does not change ancestry. Source text retrieved. |
| [QUADPACK](https://www.netlib.org/quadpack/) · ASTRO-0517 | Fortran | Adaptive one-dimensional numerical integration | Legacy quadrature algorithms behind multiple scientific interfaces. Source text retrieved. |
| [FITPACK and DIERCKX](https://www.netlib.org/dierckx/) · ASTRO-0518 | Fortran | Spline approximation, fitting and interpolation | Important interpolation lineage; basis, smoothing and extrapolation conventions matter. Source text retrieved. |
| [FFTPACK](https://www.netlib.org/fftpack/) · ASTRO-0519 | Fortran; later translations | Fourier, sine and cosine transforms | Legacy transform reference, not identical to modern FFTW or pocketfft implementations. Source text retrieved. |
| [Cuba](https://feynarts.de/cuba/) · ASTRO-0520 | C; Fortran and other interfaces | Multidimensional numerical integration | Algorithms have different stochastic and deterministic assumptions. Source text retrieved. |
| [FFTLog](https://jila.colorado.edu/~ajsh/FFTLog/) · ASTRO-0521 | Fortran; later ports | Fast logarithmic-grid Hankel and Fourier-Bessel transforms | Cosmological power-spectrum and correlation-function transformations. Source text retrieved. |
| [mcfit](https://github.com/eelregit/mcfit) · ASTRO-0522 | Python | Cosmological integral transforms based on FFTLog | Related algorithm family with reusable transform kernels. Source text retrieved. |
| [FAST-PT](https://github.com/JoeMcEwen/FAST-PT) · ASTRO-0523 | Python | Fast perturbation-theory convolution integrals | Overlaps with cosmological perturbation calculations at kernel level, not every physical approximation. Source text retrieved. |
| [fitsio Python](https://github.com/esheldon/fitsio) · ASTRO-0524 | Python; CFITSIO binding | Read and write FITS data through CFITSIO | Shares an underlying engine with Julia FITSIO and Rust fitsio bindings. Source text retrieved. |
| [ArviZ](https://github.com/arviz-devs/arviz) · ASTRO-0525 | Python | Posterior diagnostics, summaries and inference-data interchange | Unify sampler outputs without hiding sampler-specific weights or diagnostics. Source text retrieved. |
| [ChainConsumer](https://github.com/Samreay/ChainConsumer) · ASTRO-0526 | Python | Posterior-chain summaries and comparative plotting | Overlaps with GetDist, ArviZ and domain-specific posterior summaries. Source text retrieved. |
| [Snakemake](https://github.com/snakemake/snakemake) · ASTRO-0527 | Python | Dependency-aware scientific workflow execution | General orchestration reference; scientific calibration rules remain domain-specific. Source text retrieved. |
| [Nextflow](https://github.com/nextflow-io/nextflow) · ASTRO-0528 | Groovy; Java | Portable workflow and task execution | General infrastructure comparison; inclusion does not assert widespread cosmology adoption. Source text retrieved. |
