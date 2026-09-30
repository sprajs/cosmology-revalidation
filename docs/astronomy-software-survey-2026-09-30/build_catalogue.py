"""Rebuild the astronomy software catalogue from the curated records below.

This is an inventory, not a software installation or a claim of measured adoption.
Primary-source retrieval evidence is added by verify_sources.py.
"""
from pathlib import Path
import json
from collections import Counter

ROOT = Path(__file__).resolve().parent
records = []

def group(domain, text):
    for line in text.strip().splitlines():
        name, language, capability, url, note = line.split('|')
        records.append(dict(name=name.strip(), domain=domain, language=language.strip(), capability=capability.strip(), source_url=url.strip(), note=note.strip(), assessed_date='2026-09-30'))

group('01 Foundations and interoperability', '''
Astropy|Python; C extensions|Units, quantities, constants, time scales, coordinates, tables, FITS, WCS, cosmology, statistics and models|https://github.com/astropy/astropy|Core capability reference; separate from its affiliated ecosystem.
Astroquery|Python|Query astronomical archives and catalogue services|https://github.com/astropy/astroquery|Service adapters require their own schema and authentication handling.
PyVO|Python|Virtual Observatory discovery and access protocols|https://github.com/astropy/pyvo|TAP, SIA, SSA and registry interoperability.
ASDF|Python; format specification|Serialize scientific trees, arrays and metadata|https://github.com/asdf-format/asdf|A format and extension ecosystem as well as a Python implementation.
GWCS|Python|Composable generalized world coordinate transformations|https://github.com/spacetelescope/gwcs|Separate from FITS WCS conventions.
WCSLIB|C|FITS world coordinate transformations and projections|https://www.atnf.csiro.au/people/mcalabre/WCS/|Foundational native library used through many wrappers.
CFITSIO|C; Fortran interfaces|Read and write FITS images and binary tables|https://heasarc.gsfc.nasa.gov/fitsio/|Format implementation reference.
ERFA and PyERFA|C; Python bindings|Astronomical positional algorithms derived from SOFA|https://github.com/liberfa/erfa|Track SOFA lineage; PyERFA is a wrapper, not an independent model.
IAU SOFA|C; Fortran|Reference algorithms for time, Earth rotation and celestial coordinates|https://www.iausofa.org/|Standards and validation reference.
Starlink AST|C; Fortran interfaces|Coordinate frames, mappings and world coordinate systems|https://starlink.eao.hawaii.edu/starlink/AST|A substantial alternative coordinate architecture.
Starlink HDS and NDF|C; Fortran|Hierarchical and N-dimensional scientific data structures|https://starlink.eao.hawaii.edu/starlink|Retain variance, quality and history semantics in format conversion.
HEALPix and healpy|Fortran; C++; Python|Equal-area spherical pixelization, maps and spherical harmonics|https://healpix.sourceforge.io/|An algorithm and multi-language ecosystem; not just healpy.
astropy-healpix|Python; C|HEALPix indexing and coordinate operations|https://github.com/astropy/astropy-healpix|Different implementation scope from the full HEALPix suite.
MOCPy and CDS MOC|Python; Rust|Multi-order sky coverage and set operations|https://github.com/cds-astro/mocpy|Existing Rust core is relevant to reuse.
Regions|Python|Pixel and sky regions with astronomy file interfaces|https://github.com/astropy/regions|A shared selection primitive.
Reproject|Python|Reprojection and mosaicking of astronomical images|https://github.com/astropy/reproject|Flux conservation and surface brightness semantics matter.
ndcube|Python|Multidimensional data with coordinates, masks and uncertainties|https://github.com/sunpy/ndcube|Useful shared cube object model.
specutils|Python|Spectral data containers, manipulation and analysis|https://github.com/astropy/specutils|Shared spectral abstraction.
synphot|Python|Synthetic photometry and passband integration|https://github.com/spacetelescope/synphot_refactor|Calibration data are separate dependencies.
dust_extinction|Python|Extinction laws for astronomical sources|https://github.com/karllark/dust_extinction|Physical assumptions and valid wavelength domains differ by law.
dustmaps|Python|Interfaces to Galactic dust maps|https://github.com/gregreen/dustmaps|Code does not include every map dataset.
Skyfield|Python|Ephemerides, positions and observation geometry|https://github.com/skyfielders/python-skyfield|Depends on ephemeris and time data.
SPICE Toolkit|Fortran; C; IDL and MATLAB interfaces|Spacecraft and planetary geometry, time and reference frames|https://naif.jpl.nasa.gov/naif/toolkit.html|Kernels are essential data products; SpiceyPy is a Python binding.
IVOA standards|Specifications; multiple implementations|VOTable, TAP, ADQL, SAMP, UWS, SIAv2, SSA, ObsCore, MOC and HiPS|https://www.ivoa.net/documents/|Standards inventory; not counted as an executable package.
''')

group('02 Integrated suites and legacy environments', '''
AMUSE|Python with C, C++ and Fortran community codes|Couple gravitational dynamics, stellar evolution, hydrodynamics and radiative transfer|https://github.com/amusecode/amuse|Closest architectural analogue for multiphysics integration; its code list warns it can be stale.
Starlink Software Collection|Fortran; C; C++; Java; Perl; Tcl/Tk|Image and spectral reduction, coordinates, data models and pipeline tools|https://github.com/Starlink/starlink|Longstanding suite with current community and observatory stewardship.
IRAF Community Distribution|SPP; C; Fortran; CL|General astronomical image and spectral reduction|https://iraf-community.github.io/|Institutional development ended; community maintenance continues. Do not label simply dead.
PyRAF|Python|Python access to IRAF tasks and command workflows|https://github.com/iraf-community/pyraf|Community continuation of a legacy interface.
STSDAS and TABLES|SPP; CL; C; Python|Historical HST calibration and analysis tasks under IRAF|https://www.stsci.edu/institute/software_hardware/stsdas|Legacy lineage; modern HST and Python successors must be considered per task.
ESO MIDAS|Fortran; C; MIDAS command language|Munich Image Data Analysis System for images and spectra|https://www.eso.org/sci/software/esomidas/|Legacy suite worth mining; current support status requires source-specific review.
IDL|Proprietary IDL runtime|Array-oriented numerical analysis and visualization environment|https://www.nv5geospatialsoftware.com/products/IDL|Likely the language recalled by the user; Sullivan attribution not established.
IDL Astronomy Library|IDL|FITS I/O, astrometry, photometry and numerical utilities|https://github.com/wlandsman/IDLAstro|NASA website frozen in 2022 points to this repository; runtime and library licenses differ.
GNU Data Language|C++; IDL-compatible language|Open implementation of much of the IDL language and numerical environment|https://github.com/gnudatalanguage/gdl|Compatibility reference; do not assume complete IDL equivalence.
SolarSoft|Primarily IDL; other languages|Solar mission analysis libraries, calibration and common utilities|https://www.lmsal.com/solarsoft/|Large distributed ecosystem; not one homogeneous package.
PDL|Perl; C|N-dimensional numerical arrays and scientific computing in Perl|https://pdl.perl.org/|Relevant numerical foundation behind Perl scientific workflows.
ORAC-DR and PICARD|Perl; external Starlink tasks|Recipe-based instrument reduction and post-processing|https://github.com/Starlink/ORAC-DR|Current and historical pipeline knowledge; include calibration selection rules.
Gnuastro|C; command line|Astronomical arithmetic, detection, measurement and image processing|https://www.gnu.org/software/gnuastro/|Integrated command-line suite with reusable libraries.
''')

group('03 Image reduction and measurement', '''
Photutils|Python|Source detection, segmentation, aperture and PSF photometry, backgrounds|https://github.com/astropy/photutils|Major Astropy ecosystem capability.
ccdproc|Python|CCD calibration, combination, uncertainty and mask processing|https://github.com/astropy/ccdproc|Detector correction primitives.
Astro-SCRAPPY|Cython; Python|Cosmic-ray rejection based on L.A.Cosmic|https://github.com/astropy/astroscrappy|Preserve noise and saturation assumptions.
Source Extractor|C|Detect sources, deblend and measure image catalogues|https://github.com/astromatic/sextractor|Also known as SExtractor.
SEP|C; Python|Source extraction library derived from Source Extractor|https://github.com/sep-developers/sep|Related implementation lineage; do not count algorithm twice.
PSFEx|C|Model spatially varying point-spread functions|https://github.com/astromatic/psfex|Pairs with Source Extractor catalogues.
SCAMP|C|Astrometric and photometric calibration of imaging surveys|https://github.com/astromatic/scamp|Global catalogue and exposure calibration.
SWarp|C|Resample and coadd astronomical images|https://github.com/astromatic/swarp|Separate coaddition covariance from per-pixel variance.
Astrometry.net|C; Python|Blind astrometric image solving|https://github.com/dstndstn/astrometry.net|Index catalogues are separate resources.
DAOPHOT and ALLSTAR|Fortran|Crowded-field stellar PSF photometry|https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/community/STETSON/daophot/|Historical and specialist capability; distribution and rights need checking.
DOLPHOT|C|Stellar photometry with instrument-specific PSFs and calibrations|https://americano.dolphinsim.com/dolphot/|Important HST and other instrument lineage.
DoPHOT|Fortran|Automated PSF stellar photometry|https://github.com/AbhijitSaha/DoPhot|Legacy algorithm family with later implementations.
GALFIT|C|Parametric two-dimensional galaxy image decomposition|https://users.obs.carnegiescience.edu/peng/work/galfit/galfit.html|Executable availability is not the same as reusable source availability.
Imfit|C++|Galaxy image fitting and structural models|https://github.com/perwin/imfit|Alternative to GALFIT with public implementation.
The Tractor|Python; native extensions|Forward modelling and forced photometry across images|https://github.com/dstndstn/tractor|Models sources in the individual observation space.
scarlet|Python|Constrained multiband source separation and deblending|https://github.com/pmelchior/scarlet|Research deblending model family.
statmorph|Python|Nonparametric galaxy morphology measurements|https://github.com/vrodgom/statmorph|Gini, asymmetry and other morphology statistics.
HOTPANTS|C|PSF-matched image subtraction|https://github.com/acbecker/hotpants|Alard-Lupton difference imaging lineage.
ISIS image subtraction|C; scripts|Image registration and difference photometry|https://www.iap.fr/useriap/alard/package.html|Distinct from the X-ray package named ISIS.
Montage|C; Python interfaces|Astronomical image mosaics and background matching|https://montage.ipac.caltech.edu/|Mosaic workflow reference.
''')

group('04 Spectroscopy and observatory pipelines', '''
PypeIt|Python|Optical and near-infrared spectroscopic reduction|https://github.com/pypeit/PypeIt|Broad multi-instrument pipeline.
Gemini DRAGONS|Python|Recipe-driven Gemini data reduction|https://github.com/GeminiDRSoftware/DRAGONS|Includes instrument data models and calibration workflows.
AstroData|Python|Common interface for astronomical instrument data products|https://github.com/GeminiDRSoftware/astrodata|Separate reusable component of Gemini ecosystem.
JWST pipeline|Python; compiled dependencies|Detector correction, calibration, imaging and spectroscopy|https://github.com/spacetelescope/jwst|Requires instrument reference files and CRDS context.
HSTCAL|C|Hubble instrument calibration tasks|https://github.com/spacetelescope/hstcal|Includes ACS, WFC3 and STIS calibration lineage.
DrizzlePac|Python; C|HST image alignment, distortion correction and drizzle combination|https://github.com/spacetelescope/drizzlepac|Sampling and correlated-noise semantics matter.
romancal|Python|Roman WFI imaging calibration and processing|https://github.com/spacetelescope/romancal|Roman pipeline development; distinguish from completed survey use.
stpipe and stcal|Python|Shared pipeline framework and detector calibration algorithms|https://github.com/spacetelescope/stcal|Component relationship to JWST and Roman.
CRDS|Python|Select and manage calibration reference data|https://github.com/spacetelescope/crds|Scientific reproducibility requires pinning reference context.
Rubin Science Pipelines|Python; C++|Instrument signature removal, calibration, coadds, detection and measurements|https://pipelines.lsst.io/|Includes task graph and Butler data management architecture.
Rubin Butler|Python|Dataset registry, storage and provenance for processing pipelines|https://github.com/lsst/daf_butler|Reusable data-management design reference.
ESO CPL and EsoRex|C|Common pipeline library and recipe execution|https://www.eso.org/sci/software/cpl/|Instrument recipes are additional packages.
ESO Reflex|Java; Python; external recipes|Graphical reduction workflows for ESO instruments|https://www.eso.org/sci/software/esoreflex/|Workflow and recipe orchestration layer.
ESO EDPS and PyCPL|Python; C bindings|Data organisation and pipeline scheduling with Python recipe access|https://www.eso.org/sci/software/pipe_aem_main.html|Newer ESO orchestration ecosystem.
ESO instrument recipes|Mostly C via CPL; varies|MUSE, X-shooter, UVES, ESPRESSO, CRIRES, FORS, KMOS, SPHERE and VLTI reduction|https://www.eso.org/sci/software/pipelines/|Family entry; each instrument needs its own future inventory.
DESI spectroscopic pipeline|Python; native dependencies|Spectral extraction, calibration, redshifts and survey processing|https://github.com/desihub/desispec|Includes relationships to specter and redrock.
SDSS idlspec2d|IDL; supporting scripts|SDSS and BOSS optical spectroscopy reduction|https://github.com/sdss/idlspec2d|Substantial IDL survey pipeline lineage.
MaNGA DAP|Python|Integral-field spectral analysis and derived galaxy maps|https://github.com/sdss/mangadap|Separate data analysis from upstream reduction.
MPDAF|Python|MUSE integral-field data analysis|https://github.com/musevlt/mpdaf|Cube operations, spectra and images.
pPXF|Python; historical IDL|Stellar and gas kinematics and population fitting from spectra|https://pypi.org/project/ppxf/|Retain template resolution, regularization and noise conventions.
GANDALF|IDL; Python ports|Joint stellar continuum and emission-line analysis|https://ascl.net/1708.012|Legacy lineage; verify chosen port separately.
linetools|Python|Spectroscopic analysis of absorption lines and systems|https://github.com/linetools/linetools|Atomic line and redshift utilities.
SpecPro|IDL|Interactive spectrum inspection and redshift identification|https://specpro.caltech.edu/|Proprietary runtime ecosystem.
L.A.Cosmic original|IDL|Laplacian cosmic-ray identification|https://www.astro.yale.edu/dokkum/lacosmic/|Legacy algorithm reference alongside Astro-SCRAPPY.
''')

group('05 X ray and gamma ray astronomy', '''
HEASoft and FTOOLS|C; C++; Fortran; Perl; Tcl; Python|High-energy mission reduction and general FITS analysis|https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/|Large suite; component entries below are not independent suites.
XSPEC and PyXspec|C++; Fortran models; Tcl; Python|Forward-folded spectral fitting and physical emission models|https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/|Responses, backgrounds and model tables are essential.
XSELECT|Fortran; command interface|Filter event lists and extract images, spectra and light curves|https://heasarc.gsfc.nasa.gov/docs/software/ftools/xselect/|HEASoft component.
XRONOS|Fortran; C|High-energy time-series analysis|https://heasarc.gsfc.nasa.gov/docs/xanadu/xronos/|Legacy and maintained-suite component lineage.
XIMAGE|Fortran; C; Tcl|High-energy image analysis and source detection|https://heasarc.gsfc.nasa.gov/docs/xanadu/ximage/|HEASoft component.
HEASoft mission tasks|Mixed compiled languages; Perl; Python|NICERDAS, NuSTARDAS, Swift, Suzaku, IXPE and XRISM instrument processing|https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/versions.html|Family entry with explicit component/version evidence.
CIAO|C; C++; Python; other compiled components|Chandra calibration, event processing, detection and analysis|https://cxc.cfa.harvard.edu/ciao/|Includes mission-specific CALDB dependencies.
Sherpa|Python; C++; Fortran|Statistical modelling and fitting of spectra and images|https://github.com/sherpa/sherpa|Mission-independent library also distributed with CIAO.
XMM Newton SAS|C++; Fortran; Perl; Python interfaces|XMM event calibration and science products|https://www.cosmos.esa.int/web/xmm-newton/sas|Instrument calibrations and observation metadata are separate.
SPEX|Fortran; Python interface|High-resolution plasma spectral modelling and fitting|https://www.sron.nl/astrophysics-spex|Distribution and model database terms need separate review.
ISIS spectral analysis|C; S-Lang|Interactive high-resolution X-ray spectroscopy|https://space.mit.edu/cxc/isis/|Different ISIS from optical difference imaging.
Stingray|Python|Spectral timing, variability, cross spectra and light curves|https://github.com/StingraySoftware/stingray|Event timing and noise-statistics capability.
HENDRICS|Python|High-energy timing workflows built on Stingray|https://github.com/StingraySoftware/HENDRICS|Pipeline layer; avoid double-counting its numerical base.
Gammapy|Python|Gamma-ray event, map, spectrum and likelihood analysis|https://github.com/gammapy/gammapy|Observation and instrument-response data model reference.
ctools and GammaLib|C++; Python|Cherenkov gamma-ray analysis and general high-energy likelihood library|https://cta.irap.omp.eu/ctools/|Related library and tool suite.
Fermitools|C++; Python|Fermi LAT event selection, response and likelihood analysis|https://github.com/fermi-lat/Fermitools-conda|Mission suite distribution repository.
Fermipy|Python|Higher-level Fermi LAT analysis workflows|https://github.com/fermiPy/fermipy|Uses Fermitools; not an independent response implementation.
threeML|Python|Multi-instrument likelihood modelling|https://github.com/threeML/threeML|Plugin architecture for multi-messenger analysis.
Naima|Python|Nonthermal radiation modelling and inference|https://github.com/zblz/naima|Particle populations to observed spectra.
ctapipe|Python|Cherenkov telescope low-level event reconstruction|https://github.com/cta-observatory/ctapipe|Distinct from high-level Gammapy analysis.
eSASS|Mixed; consult distribution|eROSITA event and survey analysis|https://erosita.mpe.mpg.de/edr/DataAnalysis/|Release-specific accessibility and calibration dependencies.
ixpeobssim|Python|X-ray polarimetry simulation and analysis|https://github.com/lucabaldini/ixpeobssim|Include Stokes statistics and modulation response.
SOXS and pyXSIM|Python|Synthetic X-ray observations and photons from simulations|https://github.com/lynx-x-ray-observatory/soxs|Couples plasma emission, responses and simulation outputs.
MARX|C|Chandra observation and ray-trace simulation|https://space.mit.edu/cxc/marx/|Instrument simulator and validation reference.
''')

group('06 Radio interferometry and spectral cubes', '''
CASA|C++; Python|Radio calibration, imaging, deconvolution and single-dish analysis|https://casadocs.readthedocs.io/en/stable/|Full suite; preserve Measurement Set and measurement-equation semantics.
casacore|C++|Tables, measures, coordinates and radio astronomy data primitives|https://github.com/casacore/casacore|Foundational library distinct from CASA tasks.
AIPS|Fortran; C; POPS|Radio interferometry calibration and imaging|https://www.aips.nrao.edu/|Longstanding NRAO system; age alone does not imply abandoned.
ParselTongue|Python|Python scripting of AIPS|https://www.jive.eu/jivewiki/doku.php?id=parseltongue:parseltongue|Wrapper lineage and workflow compatibility.
MIRIAD|Fortran; C|Radio interferometric reduction and analysis|https://www.atnf.csiro.au/computing/software/miriad/|Legacy ecosystem with observatory-specific branches.
GILDAS|Fortran; C; SIC|CLASS spectroscopy, CLIC calibration and MAPPING imaging|https://www.iram.fr/IRAMFR/GILDAS/|CLASS here is unrelated to the cosmology Boltzmann code.
WSClean|C++|Wide-field radio synthesis imaging and deconvolution|https://wsclean.readthedocs.io/|Modern imager; dependent on shared radio data libraries.
AOFlagger|C++; Lua|Radio-frequency interference detection and flagging|https://aoflagger.readthedocs.io/|Selection and flag provenance matter.
DP3|C++; Python plugins|Visibility preprocessing, averaging and calibration|https://dp3.readthedocs.io/|LOFAR and SKA lineage.
DDFacet|Python; C/C++|Direction-dependent wide-field radio imaging|https://github.com/saopicc/DDFacet|Facet and beam-aware imaging.
killMS|Python; native kernels|Direction-dependent radio calibration|https://github.com/saopicc/killMS|Companion to DDFacet.
CubiCal|Python; accelerated kernels|Radio interferometric calibration|https://github.com/ratt-ru/CubiCal|Gain solutions and calibration equations.
QuartiCal|Python; Numba; Dask|Scalable radio interferometric calibration|https://github.com/ratt-ru/QuartiCal|Related modern calibration ecosystem.
MeqTrees|C++; Python|Numerical solution of measurement equations|https://github.com/ratt-ru/meqtrees-timba|Flexible equation-tree architecture.
RASCIL and SKA SDP components|Python; C/C++ components|Simulation, calibration, imaging and radio data models|https://developer.skao.int/projects/rascil/en/latest/|Parts migrated to ska-sdp-datamodels and ska-sdp-func-python.
OSKAR|C++; CUDA; Python|Radio interferometer and beam simulation|https://github.com/OxfordSKA/OSKAR|Instrument forward modelling.
Stimela and CARACal|Python; containerized tools|Compose radio calibration and imaging workflows|https://github.com/caracal-pipeline/caracal|Orchestration, not a distinct physical model.
PyBDSF|Python; Fortran|Radio image source detection and characterization|https://github.com/lofar-astron/PyBDSF|Former PyBDSM naming lineage.
Aegean|Python|Radio source finding and catalogue measurement|https://github.com/PaulHancock/Aegean|Includes background estimation and related tools.
SoFiA 2|C|Neutral-hydrogen spectral-cube source finding|https://github.com/SoFiA-Admin/SoFiA-2|Three-dimensional detection and reliability.
BBarolo|C++|Three-dimensional tilted-ring galaxy kinematic modelling|https://github.com/editeodoro/Bbarolo|HI and other resolved spectral cubes.
TiRiFiC|C|Tilted-ring fitting directly to spectral cubes|https://github.com/gigjozsa/tirific|Alternative forward model and fitting strategy.
spectral-cube and radio-beam|Python|WCS-aware cube analysis and radio beam operations|https://github.com/radio-astro-tools/spectral-cube|Beam and brightness-temperature equivalencies are essential.
GBTIDL|IDL|Green Bank Telescope single-dish spectral analysis|https://gbtidl.nrao.edu/|IDL observatory workflow.
dysh|Python|Single-dish radio spectral analysis|https://github.com/GreenBankObservatory/dysh|Python ecosystem for GBT and single-dish workflows.
eht-imaging|Python|VLBI imaging and modelling for black-hole-scale sources|https://github.com/achael/eht-imaging|Visibility-domain likelihoods and closure quantities.
SMILI|Python; C; Fortran|Sparse-model VLBI image reconstruction|https://github.com/astrosmili/smili|Alternative EHT imaging approach.
Difmap|C|VLBI visibility editing, imaging and model fitting|https://sites.astro.caltech.edu/~mcs/difmap/|Longstanding specialist tool.
DiFX|C++; Python; other components|Software VLBI correlation|https://github.com/difx/difx|Correlation is a distinct upstream capability.
pyuvdata and pyuvsim|Python|Radio visibility, beam and calibration formats; simulations|https://github.com/RadioAstronomySoftwareGroup/pyuvdata|Shared objects and interoperability for low-frequency arrays.
''')

group('07 Pulsars fast transients and timing', '''
PRESTO|C; Python|Pulsar searches, acceleration searches and analysis|https://github.com/scottransom/presto|Filterbank and time-domain signal-processing chain.
PSRCHIVE|C++; Python|Pulsar archives, calibration and profile analysis|https://psrchive.sourceforge.net/|Polarization and profile data model.
DSPSR|C++|Pulsar coherent dedispersion and signal processing|https://dspsr.sourceforge.net/|Streaming baseband and folded data.
TEMPO|Fortran|Pulsar timing model fitting|https://tempo.sourceforge.net/|Legacy reference with scientific importance.
TEMPO2|C++; C; plugins|Precision pulsar timing|https://bitbucket.org/psrsoft/tempo2/|Different numerical and parameter conventions require validation.
PINT pulsar timing|Python|Pulsar timing models and fitting|https://github.com/nanograv/PINT|Not the unrelated Python units package named Pint.
enterprise|Python|Pulsar timing array noise and gravitational-wave inference|https://github.com/nanograv/enterprise|Correlated timing-residual likelihoods.
PSRSIGSIM|Python|Pulsar signal and observation simulation|https://github.com/PsrSigSim/PsrSigSim|Synthetic-data validation.
SIGPROC|C; Fortran|Filterbank processing and pulsar search utilities|https://sigproc.sourceforge.net/|Legacy file formats and algorithms remain useful.
Heimdall|C++; CUDA|GPU single-pulse transient searches|https://sourceforge.net/projects/heimdall-astro/|Fast-radio-burst search lineage.
dsptools and baseband|Python|Baseband data access and processing|https://github.com/mhvk/baseband|Precise timestamps and sample formats.
''')

group('08 Strong weak and microlensing', '''
lenstronomy|Python|Strong-lens imaging, mass models, time delays and cosmography|https://github.com/lenstronomy/lenstronomy|Broad reference for reusable lens primitives.
PyAutoLens|Python|Automated strong-lens modelling and source reconstruction|https://github.com/Jammy2211/PyAutoLens|Related PyAutoGalaxy and PyAutoFit ecosystem.
Herculens|Python; JAX|Differentiable strong-lens modelling|https://github.com/Herculens/herculens|Gradient-based pixel and parametric models.
caustics|Python; PyTorch|Differentiable gravitational-lensing simulations|https://github.com/Ciela-Institute/caustics|Accelerator-oriented lens forward models.
jaxtronomy|Python; JAX|JAX implementations of lenstronomy-style lens calculations|https://github.com/lenstronomy/JAXtronomy|Related implementation lineage.
Lenstool|C; supporting interfaces|Parametric strong and weak cluster lens modelling|https://projets.lam.fr/projects/lenstool/wiki|Cluster-scale mass reconstruction.
glafic|C|Strong-lensing calculations and mass modelling|https://github.com/oguri/glafic2|Point and extended sources; cluster and galaxy applications.
gravlens and lensmodel|Compiled code; language requires inspection|Lens calculations, image solving and model optimization|https://www.physics.rutgers.edu/~keeton/gravlens/|Legacy specialist code; distribution rights and source access need review.
GLEE|Mixed; inspect distribution|Strong-lens mass and source reconstruction|https://wwwmpa.mpa-garching.mpg.de/~suyu/glee/|Access and license may be restricted; do not assume downloadable source.
GLASS|Python; native components|Free-form gravitational-lens mass modelling|https://github.com/jpcoles/glass|Distinct from the cosmological simulation package also called GLASS.
WSLAP plus|Fortran; IDL analysis lineage|Free-form cluster lens reconstruction|https://www.ifca.unican.es/users/jdiego/LensExplorer/|Specialist legacy family; source availability requires follow-up.
PyCS3|Python|Time-delay estimation from lensed quasar light curves|https://github.com/COSMOGRAIL/PyCS|Time series and microlensing variability treatment.
hierArc|Python|Hierarchical strong-lens cosmography inference|https://github.com/lenstronomy/hierArc|Population and nuisance assumptions sit above image modelling.
SNTD|Python|Analysis of strongly lensed supernova light curves|https://github.com/sntd/sntd|Connects transient and lensing modules.
MulensModel|Python; C/C++ kernels|Microlensing light curves and event modelling|https://github.com/rpoleski/MulensModel|Finite-source and binary-lens methods.
pyLIMA|Python|Microlensing event modelling and fitting|https://github.com/ebachelet/pyLIMA|Independent fitting and workflow reference.
VBMicrolensing|C++; Python bindings|Fast binary and multiple-lens magnification calculations|https://github.com/valboz/VBMicrolensing|Successor family to VBBinaryLensing.
GalSim|C++; Python|Realistic astronomical image and weak-lensing simulations|https://github.com/GalSim-developers/GalSim|PSFs, shear, noise and detector effects; not an N-body code.
ngmix|Python; compiled acceleration|Galaxy and PSF model fitting and shear estimation|https://github.com/esheldon/ngmix|Weak-lensing measurement algorithms.
metadetect|Python|Detection-aware metacalibration for weak lensing|https://github.com/esheldon/metadetect|Selection effects are part of the estimator.
TreeCorr|C++; Python|Two-point and three-point correlation functions|https://github.com/rmjarvis/TreeCorr|Shared large-scale structure and shear statistics.
LensTools|Python; native components|Weak-lensing maps, simulations and analysis|https://github.com/apetri/LensTools|Map-level statistics and ray-tracing workflows.
''')

group('09 Cosmology CMB and large scale structure', '''
CAMB|Fortran; Python|Linear cosmological perturbations and CMB/matter predictions|https://github.com/cmbant/CAMB|Boltzmann solver; distinct from inference packages.
CLASS cosmology|C; Python classy wrapper|Background and perturbation evolution and spectra|https://github.com/lesgourg/class_public|Distinct from IRAM CLASS spectroscopy.
Cobaya|Python|Theory/likelihood composition and Bayesian sampling|https://github.com/CobayaSampler/cobaya|Framework can share CAMB or CLASS engines.
CosmoMC|Fortran; Python tooling|Cosmological Monte Carlo parameter inference|https://github.com/cmbant/CosmoMC|Historical and specialist inference lineage.
MontePython|Python|Cosmological likelihood and parameter inference|https://github.com/brinckmann/montepython_public|Typically coupled to CLASS.
CosmoSIS|Python; C; C++; Fortran modules|Modular cosmological parameter-estimation pipelines|https://github.com/joezuntz/cosmosis|Useful multi-language component architecture.
CCL and pyccl|C; Python|Cosmological predictions for large-scale structure|https://github.com/LSSTDESC/CCL|Core Cosmology Library.
Colossus|Python|Cosmology, halo profiles, masses and large-scale structure|https://bitbucket.org/bdiemer/colossus/|Useful halo-physics primitives.
jax-cosmo|Python; JAX|Differentiable background and projected large-scale-structure spectra|https://github.com/DifferentiableUniverseInitiative/jax_cosmo|Approximations and probe scope differ from full Boltzmann solvers.
CosmoPower and CosmoPower-JAX|Python; TensorFlow or JAX|Neural emulators of cosmological predictions|https://github.com/alessiospuriomancini/cosmopower|Learned model domain and training artefacts are part of the scientific product.
CLASS-PT|C; Python|Perturbation-theory galaxy clustering predictions|https://github.com/Michalychforever/CLASS-PT|Extends CLASS; shared lineage.
PyBird|Python|Effective-field-theory large-scale-structure predictions|https://github.com/pierrexyz/pybird|Nonlinear galaxy clustering analysis.
velocileptors|Python|Perturbative large-scale-structure and redshift-space models|https://github.com/sfschen/velocileptors|Alternative perturbative implementation.
Corrfunc|C; Python|Fast pair counting for correlation functions|https://github.com/manodeep/Corrfunc|Reusable performance-critical statistics kernel.
nbodykit|Python; MPI; native dependencies|Large-scale structure catalogues, meshes and statistics|https://github.com/bccp/nbodykit|Historical dependency stack; inspect maintenance before adoption.
pycorr and pypower|Python; native dependencies|Survey correlation functions and power-spectrum estimation|https://github.com/cosmodesi/pycorr|DESI ecosystem; pypower is a separate companion repository.
NaMaster|C; Python pymaster|Masked-sky angular power-spectrum estimation|https://github.com/LSSTDESC/NaMaster|Mode-coupling and spin-field treatment.
pixell|Python; native libraries|CMB maps, curved-sky transforms and analysis|https://github.com/simonsobs/pixell|Map geometry and Fourier conventions.
TOAST|C++; Python|Time-ordered CMB data simulation and mapmaking|https://github.com/hpc4cmb/toast|Distributed observation-level workflows.
sotodlib|Python|Simons Observatory time-stream analysis|https://github.com/simonsobs/sotodlib|Instrument-specific analysis atop common numerical infrastructure.
Commander|Fortran; supporting languages|Bayesian CMB component separation|https://github.com/Cosmoglobe/Commander|Foreground and instrument modelling.
SMICA and Planck likelihood code|Mixed compiled code; Python interfaces|CMB component separation and likelihood evaluation|https://pla.esac.esa.int/|Family and archive entry; individual source accessibility varies.
lenspyx|Python; native transforms|Curved-sky lensing deflections and map remapping|https://github.com/carronj/lenspyx|CMB lensing simulation.
21cmFAST|C; Python|Semi-numerical reionization and 21 cm simulations|https://github.com/21cmfast/21cmFAST|Distinct early-universe modelling assumptions.
21CMMC|Python|Bayesian inference around 21cmFAST|https://github.com/21cmfast/21CMMC|Shared simulation engine.
Cobaya GetDist companion|Python|Posterior-chain summaries and visualization|https://github.com/cmbant/getdist|Canonical project name is GetDist; works beyond Cobaya.
''')

group('10 Supernovae and transient ecosystems', '''
SNANA|C; Fortran; Python; historical Perl|Supernova simulation, fitting and cosmology workflows|https://github.com/RickKessler/SNANA|Survey selection, calibration and bias correction are major capabilities.
SNCosmo|Python|Transient light-curve models, synthetic photometry and fitting|https://github.com/sncosmo/sncosmo|Shared passband and model data dependencies.
BayeSN|Python; JAX; NumPyro|Hierarchical Type Ia optical/NIR emission and distance modelling|https://github.com/bayesn/bayesn|Different population and dust assumptions from SALT workflows.
SALTShaker|Python|Train SALT-family supernova spectral models|https://github.com/djones1040/SALTShaker|Training and fitting are distinct stages.
SNooPy|Python; native components|Supernova light-curve fitting and distance methods|https://github.com/obscode/snpy|Carnegie Supernova Project lineage.
SNID|Fortran|Supernova spectral classification and redshift matching|https://people.lam.fr/blondin.stephane/software/snid/|Template library is essential scientific input.
Superfit|IDL|Supernova spectral-template fitting|https://github.com/dahowell/superfit|Legacy proprietary-language workflow.
SuperNNova|Python; PyTorch|Neural supernova light-curve classification|https://github.com/supernnova/SuperNNova|Training populations and calibration determine applicability.
MOSFiT|Python|Physical transient light-curve modelling and inference|https://github.com/guillochon/MOSFiT|Modular model architecture.
Redback|Python|Bayesian transient and multi-messenger modelling|https://github.com/nikhil-sarin/redback|Links physical models and inference engines.
TARDIS|Python; Numba; native components|Monte Carlo supernova spectral synthesis|https://github.com/tardis-sn/tardis|Atomic data and radiative-transfer assumptions are separate dependencies.
SEDONA|C++|Time-dependent radiation transport for transients|https://github.com/dnkasen/pubsed|Specialist simulation reference.
SNEC|Fortran|Supernova radiation-hydrodynamics and light curves|https://stellarcollapse.org/index.php/SNEC.html|Explosion modelling and microphysics inputs.
SkyPortal|Python; JavaScript|Collaborative transient data and follow-up management|https://github.com/skyportal/skyportal|Application/workflow capability, not a physical solver.
TOM Toolkit|Python; JavaScript|Target observation management and telescope integration|https://github.com/TOMToolkit/tom_base|Operational scheduling and provenance.
Fink|Python; Spark|Transient alert processing and scientific enrichment|https://github.com/astrolabsoftware/fink-broker|Distributed broker; external services and trained models.
Lasair|Python; service stack|Transient alert ingestion, filtering and discovery|https://github.com/lsst-uk/lasair-lsst|Service architecture rather than standalone numerical library.
ALeRCE|Python; service stack|Alert classification and time-domain discovery|https://github.com/alercebroker|Organization-level family entry.
''')

group('11 Stars populations and spectral synthesis', '''
MESA|Fortran|Detailed stellar structure and evolution with microphysics|https://github.com/MESAHub/mesa|Equation-of-state, opacity and nuclear-network layers matter independently.
GYRE|Fortran|Stellar oscillation and asteroseismic calculations|https://github.com/rhdtownsend/gyre|Can consume MESA stellar models.
ADIPLS|Fortran|Adiabatic stellar pulsation calculations|https://users-phys.au.dk/jcd/adipack.n/|Longstanding asteroseismic reference.
SSE and BSE|Fortran|Rapid single and binary stellar evolution|https://astronomy.swin.edu.au/~jhurley/|Legacy fitting-formula lineage used by later population codes.
COSMIC|Python; Fortran|Compact-binary population synthesis|https://github.com/COSMIC-PopSynth/COSMIC|Builds on SSE/BSE lineage.
COMPAS|C++; Python analysis|Rapid binary population synthesis|https://github.com/TeamCOMPAS/COMPAS|Alternative binary-evolution assumptions.
POSYDON|Python; MESA grids|Population synthesis with detailed stellar/binary evolution grids|https://github.com/POSYDON-code/POSYDON|Grid generation and interpolation are distinct from inference.
binary_c|C|Binary population synthesis and nucleosynthesis|https://binary_c.gitlab.io/|Specialist evolution framework.
BPASS|Mixed; stellar model and data distribution|Binary population synthesis and spectral predictions|https://bpass.auckland.ac.nz/|Separate downloadable grids from code access.
FSPS and python-fsps|Fortran; Python|Flexible stellar population synthesis|https://github.com/cconroy20/fsps|Stellar libraries and isochrones affect predictions.
Starburst99|Fortran; supporting tools|Population synthesis for young stellar systems|https://www.stsci.edu/science/starburst99/docs/default.htm|Legacy and model-grid reference.
Prospector|Python|Bayesian stellar-population and galaxy SED inference|https://github.com/bd-j/prospector|Commonly uses FSPS; shared models imply dependence.
Bagpipes|Python|Galaxy spectral-energy-distribution and star-formation-history fitting|https://github.com/ACCarnall/bagpipes|Model assumptions and priors explicit.
CIGALE|Python|Multiwavelength galaxy SED modelling|https://cigale.lam.fr/|Energy-balance model family.
MAGPHYS|Fortran; model libraries|Galaxy SED fitting with energy balance|https://www.iap.fr/magphys/|Distribution and model-library access need separate review.
LePhare|C++; Python; historical Fortran|Photometric redshifts and SED fitting|https://github.com/lephare-photoz/lephare|Track current rewrite and historical algorithms separately.
EAZY|Python; historical C|Photometric redshift fitting|https://github.com/gbrammer/eazy-py|Template and prior dependence.
kcorrect|Python; C; historical IDL interfaces|K corrections and galaxy SED reconstruction|https://github.com/blanton144/kcorrect|Long-lived multi-language lineage.
MOOG|Fortran|LTE stellar line analysis and spectral synthesis|https://www.as.utexas.edu/~chris/moog.html|Classic stellar abundance software.
Turbospectrum|Fortran|Stellar spectral synthesis|https://github.com/bertrandplez/Turbospectrum2019|Atmosphere and line-list resources required.
SME and PySME|C/C++; IDL or Python interfaces|Stellar spectral fitting and parameter estimation|https://github.com/AWehrhahn/SME|Spectroscopy Made Easy lineage.
iSpec|Python; external synthesis engines|Stellar spectral analysis with multiple synthesis backends|https://github.com/marblestation/iSpec|Integration example and algorithm comparison environment.
ATLAS and SYNTHE|Fortran|Stellar atmospheres and line-blanketed synthetic spectra|https://wwwuser.oats.inaf.it/castelli/sources.html|Legacy numerical and line-list lineage.
TLUSTY and SYNSPEC|Fortran|Non-LTE stellar atmospheres and spectra|https://tlusty.oca.eu/|Distinct physical approximation regime.
PHOENIX|Fortran; mixed infrastructure|Stellar and planetary atmosphere radiative transfer|https://www.hs.uni-hamburg.de/EN/For/ThA/phoenix/|Public grids do not imply full code is freely available.
''')

group('12 Radiation chemistry and atomic physics', '''
Cloudy|C++|Photoionization, plasma microphysics and emitted spectra|https://gitlab.nublado.org/cloudy/cloudy|Code and atomic/molecular data both required.
XSTAR|Fortran|Photoionized gas and X-ray spectral models|https://heasarc.gsfc.nasa.gov/lheasoft/xstar/xstar.html|Also distributed with HEASoft.
CHIANTI and ChiantiPy|IDL; Python; atomic database|Optically thin plasma emission and diagnostics|https://www.chiantidatabase.org/|Model database and two language ecosystems.
AtomDB and PyAtomDB|Atomic database; Python; compiled emissivity tools|X-ray plasma emissivities and spectral modelling|https://www.atomdb.org/|Data versions are part of the scientific model.
PyNeb|Python|Emission-line diagnostics and nebular abundances|https://github.com/Morisset/PyNeb_devel|Atomic coefficients and uncertainty propagation.
RADMC-3D|Fortran; Python tools|Dust and line radiative transfer|https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/|Monte Carlo and ray tracing.
MCFOST|Fortran|Three-dimensional dust and line radiative transfer|https://github.com/cpinte/mcfost|Disk and circumstellar applications.
SKIRT|C++|Three-dimensional Monte Carlo radiative transfer|https://github.com/SKIRT/SKIRT9|Dust, gas and synthetic observations.
Hyperion|Fortran; Python|Dust continuum radiative transfer|https://github.com/hyperion-rt/hyperion|Radiation/density grid abstractions.
LIME|C|Non-LTE molecular line transfer|https://github.com/lime-rt/lime|Irregular-grid radiative transfer.
RADEX|Fortran|Non-LTE molecular excitation and line intensities|https://personal.sron.nl/~vdtak/radex/index.shtml|Escape-probability approximation; not full spatial transport.
MOCASSIN|Fortran|Three-dimensional photoionization and dust transfer|https://mocassin.nebulousresearch.org/|Specialist Monte Carlo lineage.
TORUS|Fortran|Radiation transport and radiation hydrodynamics|https://www.astro.ex.ac.uk/people/th2/torus_html/|Model and distribution scope require inspection.
KROME|Python generator; Fortran output|Chemical reaction networks and thermal microphysics|https://bitbucket.org/tgrassi/krome/|Generated solver architecture.
Grackle|C; Fortran; Python|Cooling, chemistry and radiative heating for simulations|https://github.com/grackle-project/grackle|Reusable microphysics library.
DESPOTIC|Python|Interstellar cloud thermal and chemical modelling|https://bitbucket.org/krumholz/despotic/|Reduced physical cloud models.
pynucastro|Python; generated C++ and Fortran|Nuclear reaction networks and rate handling|https://github.com/pynucastro/pynucastro|Rates and network generation are separate capabilities.
''')

group('13 Dynamics hydrodynamics and cosmological simulations', '''
GADGET-4|C++|Cosmological N-body and hydrodynamical simulations|https://wwwmpa.mpa-garching.mpg.de/gadget4/|Separate historical GADGET-2/3 branches and descendants.
GADGET-2|C|Tree/SPH cosmological simulations|https://wwwmpa.mpa-garching.mpg.de/gadget/|Legacy public reference with many descendants.
AREPO|C|Moving-mesh hydrodynamics and gravity|https://arepo-code.org/|Public version and collaboration extensions may differ.
GIZMO|C|Gravity and multiple hydrodynamic methods|https://www.tapir.caltech.edu/~phopkins/Site/GIZMO.html|Public branch and optional physics vary.
SWIFT|C|Task-based cosmological gravity and hydrodynamics|https://github.com/SWIFTSIM/SWIFT|Parallel execution and subgrid physics framework.
RAMSES|Fortran|Adaptive-mesh cosmological hydrodynamics and gravity|https://github.com/ramses-organisation/ramses|AMR and particle-mesh methods.
Enzo|C++; Fortran; Python|Adaptive-mesh astrophysical fluid simulations|https://github.com/enzo-project/enzo-dev|Longstanding AMR code.
Enzo-E|C++|Scalable adaptive-mesh astrophysics|https://github.com/enzo-project/enzo-e|Distinct implementation using Cello infrastructure.
FLASH|Fortran; C|Adaptive-mesh multiphysics and astrophysical explosions|https://flash.rochester.edu/site/flashcode/|Access and redistribution terms require review.
Athena and Athena++|C; C++|Astrophysical hydrodynamics and MHD|https://github.com/PrincetonUniversity/athena|Athena++ is distinct from the original C implementation.
AthenaK|C++; Kokkos|Portable accelerated astrophysical fluid and relativity simulations|https://github.com/IAS-Astrophysics/athenak|GPU-capable successor family.
PLUTO|C|Astrophysical fluid, MHD and relativistic flows|https://plutocode.ph.unito.it/|Multiple physical/numerical modules.
MPI-AMRVAC|Fortran|Adaptive-mesh hydrodynamics, MHD and related physics|https://github.com/amrvac/amrvac|Extensible PDE framework.
Pencil Code|Fortran; Python|High-order compressible MHD and multiphysics|https://github.com/pencil-code/pencil-code|Turbulence and dynamo applications.
Phantom|Fortran|Smoothed-particle hydrodynamics and astrophysical dynamics|https://github.com/danieljprice/phantom|Disks and other Lagrangian-fluid applications.
ChaNGa|C++; Charm++|Parallel cosmological N-body and SPH simulations|https://github.com/N-BodyShop/changa|Task/runtime design reference.
PKDGRAV3|C; C++|Large cosmological N-body gravity simulations|https://github.com/CLS-project/PKDGRAV3|Performance-critical gravity implementation.
REBOUND|C; Python|N-body orbital integration|https://github.com/hannorein/rebound|Multiple integrators and collision treatments.
REBOUNDx|C; Python|Additional forces and physical effects for REBOUND|https://github.com/dtamayo/reboundx|Extension mechanism reference.
NBODY6 and NBODY6++GPU|Fortran; C/C++; CUDA|Direct collisional star-cluster dynamics|https://github.com/nbodyx/Nbody6ppGPU|Legacy and accelerated family; regularization and close encounters.
PeTar|C++|Star-cluster dynamics with hybrid gravitational methods|https://github.com/lwang-astro/PeTar|Coupling dynamics and stellar evolution.
CMC|C; Python interfaces|Monte Carlo star-cluster evolution|https://github.com/ClusterMonteCarlo/CMC-COSMIC|Related COSMIC evolution dependencies.
galpy|Python; C|Galactic orbits, potentials and distribution functions|https://github.com/jobovy/galpy|Useful reusable dynamics core.
gala|Python; C/Cython|Galactic dynamics and orbit analysis|https://github.com/adrn/gala|Astropy-compatible phase-space objects.
AGAMA|C++; Python|Galaxy dynamics, potentials, actions and distribution functions|https://github.com/GalacticDynamics-Oxford/Agama|Alternative dynamics implementation.
NEMO|C; C++; Fortran tools|Stellar dynamics toolbox and simulation analysis|https://github.com/teuben/nemo|Legacy integrated toolbox with ongoing stewardship.
yt|Python; Cython|Analysis and visualization of volumetric simulation data|https://github.com/yt-project/yt|Reader/analysis layer rather than a simulation engine.
pynbody|Python; C/Cython|N-body and hydrodynamic simulation analysis|https://github.com/pynbody/pynbody|Particle data, units, halos and derived quantities.
''')

group('14 Relativity spacetime and black holes', '''
Einstein Toolkit|C; C++; Fortran; Python|Numerical relativity, initial data, spacetime evolution and diagnostics|https://einsteintoolkit.org/|A suite of separately maintained thorns; Cactus and Carpet/CarpetX infrastructure.
Cactus and CarpetX|C++; C; Fortran components|Modular scientific simulation framework and AMR driver|https://github.com/eschnett/CarpetX|Infrastructure related to Einstein Toolkit; not an independent gravity theory.
SpECTRE|C++; Python|Parallel relativistic astrophysics using high-order numerical methods|https://github.com/sxs-collaboration/spectre|Distinct from the collaboration's SpEC code.
SpEC|C++; other components|Spectral numerical relativity and compact-binary simulations|https://www.black-holes.org/code/SpEC.html|Code accessibility differs from public waveform availability.
GRChombo|C++|Adaptive-mesh numerical relativity|https://github.com/GRChombo/GRChombo|Einstein-equation evolution and matter systems.
NRPy|Python; generated C/CUDA|Symbolic generation of numerical-relativity solvers|https://github.com/nrpy/nrpy|Useful symbolic-to-numerical architecture reference.
EinsteinPy|Python|Symbolic relativity, geodesics and spacetime calculations|https://github.com/einsteinpy/einsteinpy|Not a general replacement for numerical-relativity evolution codes.
xAct|Wolfram Language|Tensor algebra, perturbations and differential geometry|https://www.xact.es/|Open package ecosystem with proprietary Mathematica runtime.
Cadabra|C++; Python|Symbolic tensor algebra and field-theory calculations|https://cadabra.science/|Open symbolic environment.
SageManifolds|Python; SageMath|Differential geometry, tensor fields and relativity|https://sagemanifolds.obspm.fr/|Open symbolic alternative with numerical interfaces.
GRTensorIII|Maple|Tensor calculations in general relativity|https://github.com/grtensor/grtensor|Legacy GRTensorII lineage and proprietary Maple runtime.
Black Hole Perturbation Toolkit|Wolfram Language; Python; C/C++; others|Perturbations, Kerr orbits, waveforms and special functions|https://bhptoolkit.org/|Umbrella; inspect individual packages and data dependencies.
Teukolsky and KerrGeodesics|Wolfram Language|Black-hole perturbation solutions and Kerr orbital quantities|https://bhptoolkit.org/Teukolsky/|Components of BHPToolkit; not full nonlinear spacetime evolution.
FastEMRIWaveforms|Python; C++; CUDA|Extreme-mass-ratio-inspiral waveforms|https://github.com/BlackHolePerturbationToolkit/FastEMRIWaveforms|Specialized waveform approximations and model domains.
qnm|Python|Black-hole quasinormal modes|https://github.com/duetosymmetry/qnm|Perturbative ringdown spectrum.
GYOTO|C++; Python; Yorick|General-relativistic ray tracing|https://github.com/gyoto/Gyoto|Metric plugins and radiative models.
ipole|C|Polarized general-relativistic radiative transfer|https://github.com/moscibrodzka/ipole|Post-processing of relativistic fluid simulations.
grtrans|Fortran; Python|General-relativistic polarized radiation transport|https://github.com/jadexter/grtrans|Metric/radiation conventions need explicit comparison.
BHAC|Fortran|General-relativistic magnetohydrodynamics|https://bhac.science/|Compact-object plasma simulations.
HARM and iharm3D|C|General-relativistic magnetohydrodynamics|https://github.com/AFD-Illinois/iharm3d|HARM lineage; scope varies by implementation.
KHARMA|C++; Kokkos; Parthenon|Portable general-relativistic MHD|https://github.com/AFD-Illinois/kharma|Modern implementation lineage.
gevolution|C++|Relativistic cosmological N-body simulations|https://github.com/gevolution-code/gevolution-1.2|Weak-field relativistic cosmology; not interchangeable with full NR.
CosmoLattice|C++|Expanding-universe classical field simulations|https://cosmolattice.net/|Early-universe nonlinear fields.
''')

group('15 Gravitational waves and multi messenger inference', '''
LALSuite|C; Python bindings|Gravitational-wave waveforms and analysis algorithms|https://git.ligo.org/lscsoft/lalsuite|Umbrella for LAL, LALSimulation, inference and other components.
PyCBC|Python; C/CUDA components|Compact-binary searches and inference|https://github.com/gwastro/pycbc|Often uses LAL waveform models.
Bilby|Python|Bayesian inference for gravitational waves and other data|https://github.com/bilby-dev/bilby|Likelihood and sampler interface layer.
GWpy|Python|Gravitational-wave detector time-series analysis|https://github.com/gwpy/gwpy|Data quality and time-frequency analysis.
GWOSC client|Python|Access public gravitational-wave strain and event data|https://github.com/gwosc/client|Service/data interface.
GstLAL|C; Python; GStreamer|Streaming gravitational-wave analysis|https://git.ligo.org/lscsoft/gstlal|Low-latency pipeline architecture.
cWB|C++; ROOT|Coherent burst gravitational-wave searches|https://gwburst.gitlab.io/|Waveform-independent search family.
BayesWave|C; supporting scripts|Bayesian signal and glitch reconstruction|https://git.ligo.org/lscsoft/bayeswave|Transient noise and signal separation.
RIFT|Python; native dependencies|Rapid gravitational-wave parameter inference|https://git.ligo.org/rapidpe-rift/rift|Alternative inference architecture.
PESummary|Python|Posterior summaries and comparison reports|https://github.com/pesummary/pesummary|Provenance and result presentation layer.
ligo.skymap|Python; C|Sky localization maps and multi-messenger utilities|https://git.ligo.org/lscsoft/ligo.skymap|Probability on the sphere and distance information.
gwpopulation|Python; NumPy/CuPy/JAX options|Hierarchical compact-object population inference|https://github.com/ColmTalbot/gwpopulation|Selection and detection efficiency are central.
gwsurrogate|Python; native components|Numerical-relativity waveform surrogate evaluation|https://github.com/sxs-collaboration/gwsurrogate|Training model files have separate provenance.
PyRing|Python|Black-hole ringdown inference|https://git.ligo.org/lscsoft/pyring|Tests of remnant and perturbation models.
''')

group('16 Exoplanets planetary and solar system science', '''
Lightkurve|Python|Kepler and TESS light-curve processing|https://github.com/lightkurve/lightkurve|Time-series extraction and systematics correction.
batman|C; Python|Exoplanet transit light curves|https://github.com/lkreidberg/batman|Finite-source occultation and limb-darkening models.
starry|Python; C++|Analytic occultation and surface-map light curves|https://github.com/rodluger/starry|Differentiable map/geometry modelling.
exoplanet|Python; probabilistic-programming dependencies|Exoplanet transit and orbital inference|https://github.com/exoplanet-dev/exoplanet|Inspect current backend versions before integration.
juliet|Python|Joint transit and radial-velocity inference|https://github.com/nespinoza/juliet|Combines existing models and samplers.
RadVel|Python|Radial-velocity orbit fitting|https://github.com/California-Planet-Search/radvel|Instrument offsets and stellar variability.
allesfitter|Python|Joint exoplanet and stellar-system light-curve fitting|https://github.com/MNGuenther/allesfitter|Integrates several numerical/model engines.
PHOEBE|Python; C/C++|Eclipsing-binary light curves and stellar geometry|https://github.com/phoebe-project/phoebe2|Wilson-Devinney lineage with modern model architecture.
Wilson Devinney code|Fortran|Eclipsing binary light and radial-velocity modelling|https://www.astro.ufl.edu/~wilson/|Legacy physical-model reference.
petitRADTRANS|Python; Fortran components|Exoplanet atmospheric spectra and retrieval|https://github.com/wbalmer/petitRADTRANS|Opacity tables can dominate scientific data requirements.
TauREx|Python|Exoplanet atmospheric retrieval|https://github.com/ucl-exoplanets/TauREx3_public|Modular atmospheric forward models.
Exo-Transmit|C|Exoplanet transmission spectra|https://github.com/elizakempton/Exo_Transmit|Opacity and chemistry tables.
VPLanet|C; Python tooling|Coupled planetary evolution modules|https://github.com/VirtualPlanetaryLaboratory/vplanet|Multiphysics coupling reference.
Tudat|C++; Python|Astrodynamics, propagation and orbit determination|https://github.com/tudat-team/tudat|Spaceflight and planetary dynamics.
Orekit|Java|Orbit propagation, frames, time and estimation|https://www.orekit.org/|Important non-Python astrodynamics ecosystem.
GMAT|C++; scripting language|Mission analysis and trajectory design|https://gmat.atlassian.net/wiki/|NASA mission-analysis application.
OpenOrb|Fortran; Python bindings|Small-body orbit computation and uncertainty propagation|https://github.com/oorb/oorb|Asteroid orbit workflows.
sbpy|Python|Small-body ephemerides, photometry and activity analysis|https://github.com/NASA-Planetary-Science/sbpy|Astropy ecosystem for planetary astronomy.
USGS ISIS|C++; Python tooling|Planetary image calibration, mapping and geometry|https://github.com/DOI-USGS/ISIS3|Third distinct ISIS name in this inventory.
''')

group('17 Solar plasma and space physics', '''
SunPy|Python|Solar data discovery, maps, coordinates and time-series analysis|https://github.com/sunpy/sunpy|Major domain ecosystem analogous to Astropy.
aiapy|Python|SDO AIA calibration and analysis|https://github.com/LM-SAL/aiapy|Instrument-specific layer.
sunkit-image|Python|Solar image processing methods|https://github.com/sunpy/sunkit-image|Shared solar image analysis.
sunkit-spex|Python|Solar X-ray spectral analysis|https://github.com/sunpy/sunkit-spex|Relationship to legacy SolarSoft spectroscopy.
PlasmaPy|Python|Plasma parameters, particles and analysis|https://github.com/PlasmaPy/PlasmaPy|General plasma foundation beyond astronomy.
SpacePy|Python; C; Fortran|Space physics coordinates, radiation belts and analysis|https://github.com/spacepy/spacepy|External empirical models and datasets.
PySPEDAS|Python|Space-physics mission data and analysis|https://github.com/spedas/pyspedas|Python counterpart to IDL SPEDAS.
SPEDAS|IDL|Space Physics Environment Data Analysis System|https://spedas.org/|Proprietary-runtime mission ecosystem.
FISM and solar model archives|Mixed; models and data|Solar irradiance and space-weather model products|https://lasp.colorado.edu/lisird/|Data/model boundary; not one reusable codebase.
BOUT++|C++|Plasma-fluid PDE simulations|https://github.com/boutproject/BOUT-dev|General physics and fusion overlap.
Smilei|C++; Python|Particle-in-cell plasma simulation|https://github.com/SmileiPIC/Smilei|Kinetic plasma rather than fluid MHD.
WarpX|C++; Python; AMReX|Accelerated particle-in-cell simulations|https://github.com/BLAST-WarpX/warpx|Electromagnetic particles and fields.
EPOCH|Fortran|Particle-in-cell plasma simulations|https://github.com/epochpic/epoch|Alternative kinetic implementation.
Vlasiator|C++|Hybrid-Vlasov space-plasma simulations|https://github.com/fmihpc/vlasiator|Distribution-function physics and large-scale simulation.
''')

group('18 Julia scientific and astronomy ecosystems', '''
JuliaAstro ecosystem|Julia|Astronomical types, FITS, coordinates, cosmology, time and images|https://juliaastro.org/|Umbrella; individual packages have independent status.
AstroImages.jl|Julia|Astronomical image arrays, FITS and WCS-aware display|https://github.com/JuliaAstro/AstroImages.jl|Image/data model reference.
FITSIO.jl|Julia; CFITSIO|FITS reading and writing|https://github.com/JuliaAstro/FITSIO.jl|Wrapper, not an independent FITS engine.
SkyCoords.jl|Julia|Astronomical coordinate systems and transformations|https://github.com/JuliaAstro/SkyCoords.jl|Validate frames and precision scope.
AstroTime.jl|Julia|Astronomical time scales and epochs|https://github.com/JuliaAstro/AstroTime.jl|Precision time representations.
Cosmology.jl|Julia|Cosmological background quantities and distances|https://github.com/JuliaAstro/Cosmology.jl|Do not confuse similarly named repositories or assume full Boltzmann support.
UnitfulAstro.jl|Julia|Astronomical units for Unitful quantities|https://github.com/JuliaAstro/UnitfulAstro.jl|Shared units extension.
PSFModels.jl|Julia|Point-spread-function models and fitting|https://github.com/JuliaAstro/PSFModels.jl|Imaging primitive.
Photometry.jl|Julia|Astronomical photometry methods|https://github.com/JuliaAstro/Photometry.jl|Scope and maintenance require package-level review.
Transits.jl|Julia|Exoplanet transit and occultation calculations|https://github.com/JuliaAstro/Transits.jl|Geometry and differentiability reference.
Octofitter.jl|Julia|Joint orbital fitting and inference|https://github.com/sefffal/Octofitter.jl|Direct imaging, astrometry and related orbital observations.
VLBISkyModels.jl|Julia|Sky intensity models and interferometric observables|https://github.com/EHTJulia/VLBISkyModels.jl|EHT and VLBI forward modelling.
Comrade.jl|Julia|Bayesian VLBI imaging and modelling|https://github.com/ptiede/Comrade.jl|Probabilistic radio-imaging ecosystem.
Gradus.jl|Julia|General-relativistic ray tracing and spectral models|https://github.com/astro-group-bristol/Gradus.jl|Spacetime geodesics and accretion-disk observables.
Krang.jl|Julia|Kerr ray tracing and black-hole image models|https://github.com/dchang10/Krang.jl|Specialist relativistic modelling.
Skylight.jl|Julia|Relativistic ray tracing and radiation transfer|https://github.com/joaquinpelle/Skylight.jl|Repository identity must be verified; retain as candidate if unresolved.
Bolt.jl|Julia|Differentiable cosmological Boltzmann integration|https://github.com/xzackli/Bolt.jl|Research implementation; scope differs from CAMB/CLASS.
SymBoltz.jl|Julia|Symbolic-numeric Einstein-Boltzmann solver|https://github.com/hersle/SymBoltz.jl|Direct model equations and automatic differentiation.
''')

group('19 General numerical inference and physics foundations', '''
NumPy and SciPy|Python; C; C++; Fortran|Arrays, linear algebra, integration, optimization, statistics and special functions|https://scipy.org/|Foundational capability families, not astronomy-specific.
SymPy|Python|Symbolic mathematics and code generation|https://github.com/sympy/sympy|Algebraic foundation for physics tooling.
JAX|Python; compiled XLA backend|Automatic differentiation, array transformations and accelerators|https://github.com/jax-ml/jax|Execution/model layer rather than a science model.
Diffrax|Python; JAX|Differentiable ordinary and stochastic differential equations|https://github.com/patrick-kidger/diffrax|Useful solver architecture.
SciML DifferentialEquations.jl|Julia|ODE, SDE, DAE, delay and related numerical solvers|https://github.com/SciML/DifferentialEquations.jl|Large solver ecosystem with individually versioned components.
ModelingToolkit.jl|Julia|Symbolic system modelling and numerical code generation|https://github.com/SciML/ModelingToolkit.jl|Equation-to-solver architecture reference.
Dedalus|Python; compiled numerical libraries|Spectral PDE solvers and physical simulations|https://github.com/DedalusProject/dedalus|Equations and boundary conditions specified at high level.
FEniCSx|C++; Python|Finite-element PDE discretization and solution|https://fenicsproject.org/|Code generation and variational formulations.
Firedrake|Python; generated/native kernels|Finite-element PDE solution|https://www.firedrakeproject.org/|Alternative variational solver architecture.
deal.II|C++|Finite-element numerical methods|https://www.dealii.org/|General multiphysics infrastructure.
MFEM|C++|Modular finite-element methods|https://mfem.org/|High-order methods and accelerators.
PETSc|C; Fortran; Python interfaces|Parallel linear and nonlinear systems and time stepping|https://petsc.org/|Numerical engine used by many domain applications.
SUNDIALS|C; C++; Fortran interfaces|ODE, DAE and nonlinear solvers with sensitivities|https://sundials.readthedocs.io/|Stiff systems and sensitivity analysis.
AMReX|C++; Fortran; Python interfaces|Block-structured adaptive mesh refinement|https://github.com/AMReX-Codes/amrex|Infrastructure for astrophysics and plasma codes.
Kokkos|C++|Portable parallel numerical kernels|https://github.com/kokkos/kokkos|Execution portability reference.
emcee|Python|Ensemble Markov chain Monte Carlo|https://github.com/dfm/emcee|Sampler, not a likelihood or physical model.
dynesty|Python|Nested sampling and evidence estimation|https://github.com/joshspeagle/dynesty|Posterior and evidence workflow.
UltraNest|Python; Cython|Nested sampling and diagnostics|https://github.com/JohannesBuchner/UltraNest|Alternative sampler and validation tooling.
MultiNest and PyMultiNest|Fortran; Python|Nested sampling and multimodal evidence calculation|https://github.com/farhanferoz/MultiNest|Legacy native sampler with wrappers.
PolyChord|Fortran; C++; Python|High-dimensional nested sampling|https://github.com/PolyChord/PolyChordLite|Sampler implementation family.
NumPyro|Python; JAX|Probabilistic models and gradient-based Bayesian inference|https://github.com/pyro-ppl/numpyro|Autodifferentiation and model compilation dependency.
PyMC|Python; PyTensor|Probabilistic programming and Bayesian inference|https://github.com/pymc-devs/pymc|Alternative modelling ecosystem.
Stan|Stan language; C++|Compiled probabilistic models and inference|https://mc-stan.org/|Independent scientific modelling language.
iminuit and Minuit2|Python; C++|Numerical minimization and uncertainty estimation|https://github.com/scikit-hep/iminuit|Optimizer diagnostics and likelihood assumptions must remain visible.
celerite and tinygp|Python; C++; JAX variants|Gaussian-process time-series and covariance modelling|https://github.com/dfm/celerite2|Family entry; separate algorithms and supported kernels.
''')

group('20 General physics quantum and particle software', '''
Geant4|C++|Particle transport through matter and detector simulation|https://geant4.web.cern.ch/|Relevant to instrument response and particle astrophysics.
ROOT and RooFit|C++; Python|Particle-physics data processing, statistics and likelihood models|https://root.cern/|Large integrated analysis ecosystem.
CRPropa|C++; Python|High-energy cosmic-ray propagation|https://github.com/CRPropa/CRPropa3|Astroparticle and magnetic-field modelling.
GALPROP|C++; other components|Galactic cosmic-ray transport and diffuse emission|https://galprop.stanford.edu/|Data, source populations and propagation assumptions matter.
CORSIKA|Fortran; C++ successor|Extensive air-shower simulation|https://www.iap.kit.edu/corsika/|Distinguish CORSIKA 7 and CORSIKA 8 implementations.
QuTiP|Python; Cython|Open quantum-system dynamics|https://github.com/qutip/qutip|General quantum physics capability.
QuantumOptics.jl|Julia|Quantum optics and open-system simulation|https://github.com/qojulia/QuantumOptics.jl|Julia counterpart in a related capability area.
Kwant|Python; C/C++|Quantum transport in tight-binding systems|https://kwant-project.org/|General condensed-matter physics.
Meep|C++; Python; Scheme|Finite-difference time-domain electromagnetics|https://github.com/NanoComp/meep|Useful for optics and instrument modelling.
MPB|C; Scheme; Python interfaces|Photonic band-structure calculation|https://github.com/NanoComp/mpb|Frequency-domain electromagnetic eigenproblems.
LAMMPS|C++|Classical molecular dynamics and particle simulation|https://github.com/lammps/lammps|General physics scope extension.
GROMACS|C++; C|Molecular dynamics|https://gitlab.com/gromacs/gromacs|Primarily molecular/biophysical domain.
OpenMM|C++; Python|Molecular simulation and programmable particle forces|https://github.com/openmm/openmm|Reusable GPU simulation architecture.
Quantum ESPRESSO|Fortran; C|Electronic structure and density-functional theory|https://www.quantum-espresso.org/|Materials physics, beyond core astronomy.
GPAW and ASE|Python; C; numerical libraries|Electronic structure and atomistic simulation workflows|https://gpaw.readthedocs.io/|Related but distinct calculator and environment.
FeynCalc|Wolfram Language|Symbolic quantum-field-theory calculations|https://github.com/FeynCalc/feyncalc|Open package with proprietary runtime.
FeynRules|Wolfram Language|Lagrangian models and particle-physics interaction rules|https://feynrules.irmp.ucl.ac.be/|Model-generation capability.
''')

group('21 Visualization archives and operations', '''
SAOImage DS9|C++; Tcl/Tk|FITS image display, regions and interactive analysis|https://ds9.si.edu/site/Home.html|Widely integrated viewer; do not reduce scope to plotting.
TOPCAT and STILTS|Java|Catalogue exploration, crossmatching and table processing|https://www.star.bris.ac.uk/~mbt/topcat/|GUI and command-line pair on STIL library.
Aladin Desktop and Aladin Lite|Java; JavaScript|Sky atlas, catalogue overlays and survey visualization|https://aladin.cds.unistra.fr/|Desktop and browser implementations differ.
Glue|Python|Linked multidimensional scientific data exploration|https://github.com/glue-viz/glue|Cross-dataset selection and visualization.
Ginga|Python|Astronomical image viewer and plugin framework|https://github.com/ejeschke/ginga|Embedded or standalone user interface.
Jdaviz|Python; web widgets|Interactive imaging, spectral and cube analysis|https://github.com/spacetelescope/jdaviz|STScI analysis application ecosystem.
CARTA|C++; TypeScript|Interactive large radio-image and cube visualization|https://cartavis.org/|Client/server architecture for large data.
VisIt|C++; Python|Large scientific-simulation visualization|https://visit-dav.github.io/visit-website/|General simulation visualization.
ParaView and VTK|C++; Python|Scientific visualization, meshes and volume rendering|https://www.paraview.org/|General infrastructure and application pair.
yt OpenSpace and WorldWide Telescope ecosystem|Mixed|Interactive astronomical and simulation visualization|https://www.openspaceproject.com/|OpenSpace is the source for this row; other named systems require their own review.
Astroplan|Python|Observation planning, observability and scheduling|https://github.com/astropy/astroplan|Constraints depend on sites, ephemerides and time conventions.
INDI|C++; XML protocol|Astronomical instrument control|https://github.com/indilib/indi|Hardware/operations boundary.
ASCOM and Alpaca|Windows/.NET; network protocol|Telescope and instrument interoperability|https://ascom-standards.org/|Protocol compatibility matters more than implementation language.
RTS2|C++; Python|Robotic telescope control and scheduling|https://github.com/RTS2/rts2|Operations rather than science inference.
''')

group('22 Existing Rust scientific building blocks', '''
ndarray|Rust|N-dimensional arrays and numerical operations|https://github.com/rust-ndarray/ndarray|Candidate foundation; benchmark needed, not chosen automatically.
nalgebra|Rust|Linear algebra and geometric transformations|https://github.com/dimforge/nalgebra|Candidate foundation.
faer|Rust|Dense and sparse linear algebra|https://github.com/sarah-ek/faer-rs|Candidate numerical engine.
uom|Rust|Type-safe units of measurement|https://github.com/iliekturtles/uom|Astronomical equivalencies, logarithmic units and runtime units need separate design.
hifitime|Rust|High-precision epochs, durations and time scales|https://github.com/nyx-space/hifitime|Potential astronomy/astrodynamics foundation.
ANISE|Rust; Python bindings|Spacecraft geometry, ephemerides and reference frames|https://github.com/nyx-space/anise|Existing Rust astrodynamics capability.
Nyx|Rust|Astrodynamics and orbit determination|https://github.com/nyx-space/nyx|Existing domain engine.
rust-fitsio|Rust; CFITSIO binding|FITS data access|https://github.com/simonrw/rust-fitsio|FFI binding; not pure Rust FITS implementation.
CDS HEALPix Rust|Rust|HEALPix geometry and spatial indexing|https://github.com/cds-astro/cds-healpix-rust|Existing astronomy kernel.
CDS MOC Rust|Rust|Multi-order spatial and space-time coverage|https://github.com/cds-astro/cds-moc-rust|Core related to MOCPy.
RustFFT|Rust|Fast Fourier transforms|https://github.com/ejmahler/RustFFT|Shared signal-processing kernel.
argmin|Rust|Numerical optimization|https://github.com/argmin-rs/argmin|Candidate optimizer framework.
diffsol|Rust|Differential-equation solvers|https://github.com/martinjrobins/diffsol|Assess stiff solvers and sensitivity support against required physics.
Burn|Rust|Tensor computation and machine learning|https://github.com/tracel-ai/burn|Potential accelerator/model backend, not a physics package.
''')

group('23 Optical instruments adaptive optics and interferometry', '''
HCIPy|Python|Wave optics, atmospheric turbulence, adaptive optics and coronagraphs|https://github.com/ehpor/hcipy|End-to-end high-contrast instrument models.
POPPY|Python|Physical-optics propagation and point-spread functions|https://github.com/spacetelescope/poppy|Reusable optical propagation engine.
STPSF and WebbPSF|Python|Space-telescope instrument PSF simulation|https://github.com/spacetelescope/stpsf|STPSF is the successor name; optical calibration data are required.
dLux|Python; JAX|Differentiable optical-system modelling|https://github.com/LouisDesdoigts/dLux|Wavefront and telescope inference.
COMPASS|C++; CUDA; Python|Adaptive-optics system simulation|https://github.com/ANR-COMPASS/shesha|GPU wavefront and control simulations.
Soapy|Python|Adaptive-optics simulation|https://github.com/AOtools/soapy|Alternative modular AO implementation.
AOtools|Python|Adaptive-optics numerical utilities|https://github.com/AOtools/aotools|Atmospheric screens, wavefronts and related primitives.
VIP|Python|High-contrast imaging and faint-companion analysis|https://github.com/vortex-exoplanet/VIP|Post-processing, detection and characterization.
pyKLIP|Python|PSF subtraction and forward modelling for direct imaging|https://bitbucket.org/pyKLIP/pyklip/|Karhunen-Loeve image-projection algorithms.
PMOIRED|Python|Optical-interferometric data modelling|https://github.com/amerand/PMOIRED|OIFITS observables and parametric models.
MiRA|Yorick; C components|Image reconstruction from optical interferometry|https://github.com/emmt/MiRA|Optimization under incomplete Fourier information.
SQUEEZE|C++|Optical-interferometric image reconstruction|https://github.com/fabienbaron/squeeze|Alternative reconstruction and regularization strategies.
LITpro|Java client; numerical services|Parametric optical-interferometry model fitting|https://www.jmmc.fr/english/tools/proposal-preparation/litpro/|Distributed service access and model implementation differ.
OIFITS and OIFITSlib|Standard; C implementations|Exchange optical-interferometric observations and uncertainties|https://www.jmmc.fr/oifits/|Data-standard entry; implementations need separate audit.
''')

group('24 Additional legacy R and domain specialist software', '''
ProFound|R; C++|Source detection, segmentation and photometry|https://github.com/asgr/ProFound|Important astronomy ecosystem outside Python.
ProFit|R; C++|Bayesian galaxy image modelling|https://github.com/ICRAR/ProFit|Galaxy fitting and structural inference.
ProSpect|R; C++|Galaxy spectral-energy-distribution modelling|https://github.com/asgr/ProSpect|Stellar populations, dust and fitting.
celestial|R|Astronomical coordinate and cosmology utilities|https://github.com/asgr/celestial|Alternative-language foundation.
KAPPA|Fortran; C|General image and NDF manipulation and analysis|https://starlink.eao.hawaii.edu/docs/sun95.htx/sun95.html|Starlink component; useful mature numerical task catalogue.
SMURF|C; Fortran|Submillimetre time-stream mapmaking and reduction|https://starlink.eao.hawaii.edu/docs/sun258.htx/sun258.html|JCMT detector and mapmaking knowledge.
FIGARO and DIPSO|Fortran; C components|Spectral reduction, fitting and interactive analysis|https://starlink.eao.hawaii.edu/docs/sun86.htx/sun86.html|Legacy Starlink spectral-analysis lineage; FIGARO documentation cited.
CCDPACK|Fortran; C|CCD calibration, registration and combination|https://starlink.eao.hawaii.edu/docs/sun139.htx/sun139.html|Legacy pipeline algorithms and parameter semantics.
CUPID|C; Fortran|Identify and characterize emission clumps in images and cubes|https://starlink.eao.hawaii.edu/docs/sun255.htx/sun255.html|Multiple detection algorithms with different assumptions.
PGPLOT|Fortran; C bindings|Scientific plotting used by legacy astronomy software|https://sites.astro.caltech.edu/~tjp/pgplot/|Preserve diagnostic plots when recovering old codes.
SM SuperMongo|C; command language|Interactive scientific plotting and scripting|https://www.astro.princeton.edu/~rhl/sm/|Historical astronomy analysis environment.
MPFIT|IDL; C ports|Nonlinear least-squares fitting based on MINPACK|https://pages.physics.wisc.edu/~craigm/idl/idl.html|Shared optimizer lineage behind many older analyses.
IDLSPEC2D and idlutils utilities|IDL; C support|SDSS analysis utilities, spectra, masks and fitting|https://github.com/sdss/idlutils|Utilities distinct from the idlspec2d pipeline entry.
VARTOOLS|C|Astronomical light-curve analysis|https://www.astro.princeton.edu/~jhartman/vartools.html|Broad time-series task collection.
VaST|C; shell; other components|Variability search from imaging time series|https://github.com/kirxkirx/vast|Legacy-friendly photometric workflow.
''')

group('25 Additional simulation and emerging modelling capabilities', '''
GRTeclyn|C++; AMReX|Accelerated numerical relativity|https://github.com/GRTLCollaboration/GRTeclyn|Successor to GRChombo; preserve historical code as a comparison target.
Castro|C++; Fortran; AMReX|Compressible astrophysical radiation hydrodynamics|https://github.com/AMReX-Astro/Castro|Stellar explosions and microphysics coupling.
MAESTROeX|C++; Fortran; AMReX|Low-Mach-number stellar hydrodynamics|https://github.com/AMReX-Astro/MAESTROeX|Filters acoustic dynamics; different validity regime from compressible codes.
MUSIC|C++|Multiscale cosmological initial conditions|https://github.com/ohahn/MUSIC|Initial-condition generation is separate from evolution.
ROCKSTAR|C|Phase-space halo identification and merger analysis|https://bitbucket.org/gfcstanford/rockstar/|Halo definition changes downstream catalogues.
AHF|C|Adaptive Meshigauss halo finding|https://popia.ft.uam.es/AHF/|Alternative halo definitions and substructure.
VELOCIraptor|C++|Phase-space structure and halo identification|https://github.com/pelahi/VELOCIraptor-STF|Particle catalogues and merger workflows.
NIFTy|Python; JAX components|Information field theory and Bayesian field reconstruction|https://github.com/NIFTy-PPL/NIFTy|Astronomy and general inverse problems.
JaxPM|Python; JAX|Differentiable particle-mesh cosmological simulations|https://github.com/DifferentiableUniverseInitiative/JaxPM|Gradient-based initial-condition and cosmology inference.
JAX-GalSim|Python; JAX|Differentiable astronomical image simulation|https://github.com/GalSim-developers/JAX-GalSim|Project documentation explicitly labels early development and directs scientific use to reference GalSim.
AstroML|Python|Astronomical statistics, machine learning and examples|https://github.com/astroML/astroML|Useful algorithm and example collection.
sbi|Python; PyTorch|Simulation-based Bayesian inference|https://github.com/sbi-dev/sbi|Validation must include simulator support and coverage.
nuSQuIDS|C++; Python bindings|Neutrino evolution and oscillations in matter|https://github.com/arguelles/nuSQuIDS|Astroparticle and general physics capability.
SNEWPY|Python|Supernova-neutrino models, flavor transformations and detector predictions|https://github.com/SNEWS2/snewpy|Links source simulations to SNOwGLoBES detector calculations.
SNOwGLoBES|C; data and scripts|Supernova-neutrino detector event-rate calculations|https://github.com/SNOwGLoBES/snowglobes|Cross sections and detector tables are essential inputs.
IceTray|C++; Python|IceCube event processing framework|https://github.com/icecube/icetray-public|Public framework does not imply every collaboration module is public.
Bifrost|Fortran; IDL/Python analysis|Radiation-MHD solar atmosphere simulations|https://github.com/ITA-Solar/Bifrost|Source and build availability depend on public distribution.
MURaM|C++; Fortran lineage|Solar and stellar radiative magnetohydrodynamics|https://www2.mps.mpg.de/projects/solar-mhd/muram_site/|Distribution conditions and branches require checking.
RH|C; Python tools|Non-LTE radiative transfer and solar spectral modelling|https://github.com/ITA-Solar/rh|Atomic models and formal-solver assumptions.
STiC|C++; Python|Non-LTE spectropolarimetric inversion|https://github.com/jaimedelacruz/stic|Infers atmospheric properties from polarized spectra.
NICOLE|Fortran; Python tools|Non-LTE inversion and polarized spectral synthesis|https://github.com/hsocasnavarro/NICOLE|Alternative solar-atmosphere inversion lineage.
AstroLib.jl|Julia|Astronomical and astrophysical utility algorithms|https://github.com/JuliaAstro/AstroLib.jl|Explicit alternative-language analogue of historical utilities.
Trixi.jl|Julia|High-order hyperbolic PDE solvers with adaptive meshes|https://github.com/trixi-framework/Trixi.jl|General conservation-law physics engine.
GLoW|Python; C; Cython|Wave-optics gravitational lensing and amplification factors|https://github.com/glow-astro/GLoW|Diffraction and interference extend geometric-optics lensing; connects directly to gravitational-wave inference.
''')

group('26 Shared data storage execution and numerical infrastructure', '''
HDF5|C; Fortran; C++ interfaces|Hierarchical scientific arrays, metadata, chunking and compression|https://github.com/HDFGroup/hdf5|Shared storage infrastructure beneath many simulation and observational workflows.
h5py|Python; Cython; HDF5 binding|NumPy-oriented interface to HDF5 datasets|https://github.com/h5py/h5py|Wrapper around HDF5; not an independent storage engine.
netCDF and netCDF4-Python|C; Fortran; Python|Self-describing multidimensional scientific datasets|https://github.com/Unidata/netcdf4-python|Some netCDF variants use HDF5; preserve schema and dimension semantics.
Zarr|Specification; Python and other implementations|Chunked compressed N-dimensional array storage|https://github.com/zarr-developers/zarr-python|Relevant to object storage and large cubes; not a scientific observation model by itself.
xarray|Python|Named dimensions, coordinates and labelled multidimensional arrays|https://github.com/pydata/xarray|Useful shared data abstraction; can use Dask and several storage backends.
Dask|Python|Task graphs, parallel arrays and distributed data processing|https://github.com/dask/dask|Execution and chunk scheduling overlap with bespoke scientific pipeline frameworks.
Apache Arrow and Parquet|C++; Rust; Java; Python; format specifications|Columnar memory interchange and analytical table storage|https://arrow.apache.org/|Data interchange is distinct from astronomy metadata semantics; Arrow and Parquet serve different roles.
Polars|Rust; Python|Columnar table processing, lazy queries and joins|https://github.com/pola-rs/polars|Existing Rust engine relevant to catalogue processing; astronomy types still need integration.
pandas|Python; Cython; native dependencies|In-memory table cleaning, joins, grouping and statistics|https://github.com/pandas-dev/pandas|Overlaps with Astropy tables at generic operations, with different units and metadata behaviour.
Vaex|Python; C++|Large tabular data exploration and aggregation|https://github.com/vaexio/vaex|Catalogue-scale processing; evaluate data model and execution differences.
DuckDB|C++; Python and other interfaces|Embedded analytical SQL over files and tables|https://github.com/duckdb/duckdb|Candidate common catalogue-query service, not an astronomy-specific package.
fsspec|Python|Common interfaces to local and remote filesystems|https://github.com/fsspec/filesystem_spec|Centralize access, caching and range requests; filesystem access is separate from scientific format decoding.
MPI and mpi4py|Standard; C/Fortran implementations; Python bindings|Distributed communication and collective numerical execution|https://github.com/mpi4py/mpi4py|Bindings and MPI implementations are separate layers, not independent physics engines.
FFTW|C; Fortran interfaces|Fast Fourier transforms and reusable execution plans|https://www.fftw.org/|Shared numerical kernel across signal, map and simulation workflows.
ducc|C++; Python bindings|FFT, spherical-harmonic transforms, gridding and related numerical primitives|https://github.com/mreineck/ducc|Includes development lineage related to libsharp; useful radio and CMB overlap.
libsharp|C|Spherical-harmonic transforms for scientific maps|https://github.com/Libsharp/libsharp|Shared transform engine; compare normalization and spin conventions with alternatives.
SHTns|C; Python interfaces|Spherical-harmonic transforms|https://bitbucket.org/nschaeff/shtns/|Alternative implementation and execution architecture.
BLAS and LAPACK|Fortran; C interfaces; multiple implementations|Dense linear algebra, matrix factorizations and eigenproblems|https://www.netlib.org/lapack/|Specifications, reference implementation and optimized implementations must be distinguished.
OpenBLAS|C; Fortran; assembly|Optimized BLAS and LAPACK-related numerical kernels|https://github.com/OpenMathLib/OpenBLAS|Multiple frontend packages may call the same underlying library.
SuiteSparse|C; C++; MATLAB interfaces|Sparse matrix factorization and graph-based numerical kernels|https://github.com/DrTimothyAldenDavis/SuiteSparse|Shared core for inverse problems and sparse solvers.
GNU Scientific Library|C|Special functions, integration, interpolation, random numbers, optimization and ODE solvers|https://www.gnu.org/software/gsl/|Broad overlapping numerical foundation used by many compiled scientific packages.
MINPACK|Fortran; multiple ports and wrappers|Nonlinear least squares and nonlinear equations|https://www.netlib.org/minpack/|Legacy algorithm lineage behind several fitting interfaces; ports are not independent methods.
ODEPACK|Fortran|Ordinary differential equations including stiff and nonstiff systems|https://www.netlib.org/odepack/|Legacy LSODE/LSODA family; modern wrapper availability does not change ancestry.
QUADPACK|Fortran|Adaptive one-dimensional numerical integration|https://www.netlib.org/quadpack/|Legacy quadrature algorithms behind multiple scientific interfaces.
FITPACK and DIERCKX|Fortran|Spline approximation, fitting and interpolation|https://www.netlib.org/dierckx/|Important interpolation lineage; basis, smoothing and extrapolation conventions matter.
FFTPACK|Fortran; later translations|Fourier, sine and cosine transforms|https://www.netlib.org/fftpack/|Legacy transform reference, not identical to modern FFTW or pocketfft implementations.
Cuba|C; Fortran and other interfaces|Multidimensional numerical integration|https://feynarts.de/cuba/|Algorithms have different stochastic and deterministic assumptions.
FFTLog|Fortran; later ports|Fast logarithmic-grid Hankel and Fourier-Bessel transforms|https://jila.colorado.edu/~ajsh/FFTLog/|Cosmological power-spectrum and correlation-function transformations.
mcfit|Python|Cosmological integral transforms based on FFTLog|https://github.com/eelregit/mcfit|Related algorithm family with reusable transform kernels.
FAST-PT|Python|Fast perturbation-theory convolution integrals|https://github.com/JoeMcEwen/FAST-PT|Overlaps with cosmological perturbation calculations at kernel level, not every physical approximation.
fitsio Python|Python; CFITSIO binding|Read and write FITS data through CFITSIO|https://github.com/esheldon/fitsio|Shares an underlying engine with Julia FITSIO and Rust fitsio bindings.
ArviZ|Python|Posterior diagnostics, summaries and inference-data interchange|https://github.com/arviz-devs/arviz|Unify sampler outputs without hiding sampler-specific weights or diagnostics.
ChainConsumer|Python|Posterior-chain summaries and comparative plotting|https://github.com/Samreay/ChainConsumer|Overlaps with GetDist, ArviZ and domain-specific posterior summaries.
Snakemake|Python|Dependency-aware scientific workflow execution|https://github.com/snakemake/snakemake|General orchestration reference; scientific calibration rules remain domain-specific.
Nextflow|Groovy; Java|Portable workflow and task execution|https://github.com/nextflow-io/nextflow|General infrastructure comparison; inclusion does not assert widespread cosmology adoption.
''')

# Corrections made after primary-page identity and lineage review.
corrections = {
    'dsptools and baseband': {'name': 'baseband'},
    'Cobaya GetDist companion': {'name': 'GetDist'},
    'yt OpenSpace and WorldWide Telescope ecosystem': {'name': 'OpenSpace', 'capability': 'Interactive visualization of astronomical datasets and space missions', 'note': 'Scientific visualization application; distinct from yt and WorldWide Telescope.'},
    'Aladin Desktop and Aladin Lite': {'language': 'Java; JavaScript; Rust/WebAssembly'},
    'SPEX': {'source_url': 'https://www.sron.nl/en/pillars/science/astrophysics/', 'note': 'SRON describes SPEX as GPLv3 open source with source distributed through Zenodo.'},
    'hierArc': {'source_url': 'https://github.com/TDCOSMO/hierArc'},
    'SNTD': {'source_url': 'https://github.com/jpierel14/sntd'},
    'GWOSC client': {'source_url': 'https://github.com/gwpy/gwosc'},
    'TauREx': {'source_url': 'https://taurex3.readthedocs.io/en/latest/user/'},
    'petitRADTRANS': {'source_url': 'https://gitlab.com/mauricemolli/petitRADTRANS'},
    'GRChombo': {'note': 'Primary repository identifies GRChombo as predecessor to GRTeclyn; retain both implementations.'},
    'Gradus.jl': {'note': 'GitHub repository states that development moved to Codeberg; follow its link for current source.'},
    'Skylight.jl': {'note': 'Primary repository identity confirmed for general-relativistic ray tracing and radiative transfer.'},
    'ESO MIDAS': {'note': 'ESO lists it as legacy software while still advertising a 25MAY release; legacy does not mean abandoned.'},
    'Teukolsky and KerrGeodesics': {'source_url': 'https://github.com/BlackHolePerturbationToolkit/Teukolsky'},
    'FISM and solar model archives': {'name': 'LISIRD solar model and data archive', 'note': 'Data/model service entry; not an independent numerical software package.'},
    'IDLSPEC2D and idlutils utilities': {'name': 'idlutils'},
    'pPXF': {'source_url': 'https://users.physics.ox.ac.uk/~cappellari/software/'},
    'AtomDB and PyAtomDB': {'source_url': 'https://github.com/AtomDB/pyatomdb'},
    'MUSIC': {'source_url': 'https://www.oca.eu/en/olivier-hahn/1914-music'},
    'AHF': {'capability': 'Adaptive Mesh Investigations of Galaxy Assembly halo finder'},
    'Bifrost': {'source_url': 'https://ita-solar.github.io/Bifrost/Start_here/', 'note': 'Official documentation is public; source-repository access and distribution require separate confirmation.'},
    'GLEE': {'source_url': 'https://arxiv.org/abs/2209.03094', 'language': 'Implementation language not confirmed in this survey', 'note': 'Author paper documents strong-lens modelling and comparisons; public source-distribution status not established.'},
    'WSLAP plus': {'source_url': 'https://arxiv.org/abs/1304.2393', 'language': 'Implementation language not confirmed in this survey', 'note': 'Primary author paper documents WSLAP+ reconstruction; surviving source distribution requires follow-up.'},
    'COMPASS': {'note': 'This GitHub repository is marked obsolete and points to an Observatoire de Paris GitLab successor; preserve lineage when acquiring source.'},
    'SoFiA 2': {'note': 'The GitHub landing page points to a GitLab successor repository; acquisition must follow the migration.'},
    'Tudat': {'note': 'The original C++ repository is archived and points to tudatpy, which now contains the C++ core; follow the current distribution.'},
    'faer': {'note': 'GitHub landing page states that development migrated to Codeberg and GitHub remains a mirror.'},
    'ADIPLS': {'note': 'The cited landing page explicitly describes an obsolete adipack version and links to its replacement; retain the lineage.'},
}
for record in records:
    record.update(corrections.get(record['name'], {}))

# Stable registry fallbacks preserve software identity when an original site is
# unavailable. They do not establish source availability or current maintenance.
registry_fallbacks = {
    'DAOPHOT and ALLSTAR': '1104.011', 'DOLPHOT': '1608.013',
    'Montage': '1010.036', 'SpecPro': '1404.014',
    'ParselTongue': '1208.020', 'Lenstool': '1102.004',
    'binary_c': '2307.035', 'PHOENIX': '1010.056',
    'RADEX': '1010.075', 'TORUS': '1404.006', 'GIZMO': '1410.003',
    'GALPROP': '1010.028', 'MPFIT': '1208.019', 'AHF': '1102.009',
    'STSDAS and TABLES': '1206.003',
    'L.A.Cosmic original': '1207.005',
}
for record in records:
    if record['name'] in registry_fallbacks:
        record['original_source_url'] = record['source_url']
        record['source_url'] = 'https://ascl.net/' + registry_fallbacks[record['name']]
        record['source_kind'] = 'ASCL registry record; original site also retained'
    else:
        record['source_kind'] = 'Author paper' if 'arxiv.org' in record['source_url'] else 'Project, institution or standards source'

assert len({(r['domain'], r['name']) for r in records}) == len(records)
for i, r in enumerate(records, 1):
    r['id'] = f'ASTRO-{i:04d}'
    r['evidence_status'] = 'Primary source nominated; retrieval and content review tracked separately'

(ROOT / 'catalogue.json').write_text(json.dumps(records, indent=2, ensure_ascii=False) + '\n')
print(json.dumps({'entries': len(records), 'domains': dict(Counter(r['domain'] for r in records))}, indent=2))
