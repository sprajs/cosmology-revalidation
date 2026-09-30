# Where cosmology software repeats work

This companion groups the survey by jobs rather than package names. It covers the full path from raw files and archives through calibration, measurement, physical calculations, simulation, inference and inspection. The 36 families below identify **capability overlap and consolidation opportunities**. They are not measurements of redundant runtime or proof that all listed programs use the same algorithms.

The useful distinction is between four kinds of overlap:

1. **One engine behind several interfaces.** CFITSIO is exposed through Python fitsio, Julia FITSIO and Rust fitsio. Calling these separate packages does not mean three independently developed FITS engines.
2. **Related code or algorithm ancestry.** SEP and Source Extractor, Astro-SCRAPPY and L.A.Cosmic, SOFA and ERFA, and MINPACK-related fitters belong to identifiable lineages.
3. **Independent implementations of similar calculations.** CAMB and CLASS, several lens modelling packages, radio imagers and stellar-population fitters overlap in purpose. Their independence can be valuable for checking results.
4. **Repeated workflow work.** Loading the same observations, converting arrays, computing the same background distances, convolving the same response or rebuilding a covariance factorization can happen at application boundaries. Whether a particular workflow actually repeats it requires profiling and dependency inspection.

The individual projects and primary sources are linked in each row. The consolidation column is a survey synthesis, not a claim that an upstream package endorses this design.

## End-to-end scope

| Stage | Inputs and repeated operations to inventory | Relevant families |
|---|---|---|
| Acquire and decode | Archive queries, remote reads, FITS/HDF5/Zarr/ASDF, catalogue files, event streams, visibilities, simulation snapshots | Astroquery, PyVO, CFITSIO, h5py, casacore, yt, fsspec |
| Represent and calibrate | Units, epochs, coordinates, masks, uncertainty, detector responses, reference selection | Astropy, ERFA/SOFA, Starlink, CRDS, observatory pipelines |
| Measure | Sources, spectra, time series, redshifts, shapes, map statistics and selection functions | Photutils, PypeIt, HEASoft, CASA, GalSim, TreeCorr |
| Predict | Background cosmology, transfer functions, matter power, lensing, populations, radiation and instruments | CAMB, CLASS, CCL, lens codes, FSPS, radiative-transfer codes |
| Simulate | Initial conditions, particles, meshes, gravity, hydro/MHD, spacetime, mock observations | MUSIC, GADGET, AREPO, RAMSES, Einstein Toolkit, AMUSE |
| Infer and compare | Theory calls, likelihoods, covariance, priors, optimization, posterior sampling and diagnostics | Cobaya, CosmoSIS, SNANA, Bilby, samplers, GetDist, ArviZ |
| Inspect and reproduce | Linked views, fit diagnostics, intermediate artifacts, workflow execution and lineage | DS9, TOPCAT, Glue, CARTA, Rubin Butler, ORAC-DR, workflow engines |

## Capability overlap map

### 01. FITS reading and writing

[Astropy](https://github.com/astropy/astropy), [CFITSIO](https://heasarc.gsfc.nasa.gov/fitsio/), [FITSIO.jl](https://github.com/JuliaAstro/FITSIO.jl), [rust-fitsio](https://github.com/simonrw/rust-fitsio), [fitsio Python](https://github.com/esheldon/fitsio)

**Relationship:** Both independent interfaces and wrappers around CFITSIO.

**Potential shared work:** One shared reader/writer layer, header handling, memory mapping and compression service.

**Keep explicit:** Astropy is not simply another CFITSIO wrapper; FITS scaling, nulls, masks, extension types and ownership differ.

### 02. Scientific arrays and file containers

[ASDF](https://github.com/asdf-format/asdf), [Starlink HDS and NDF](https://starlink.eao.hawaii.edu/starlink), [ndcube](https://github.com/sunpy/ndcube), [HDF5](https://github.com/HDFGroup/hdf5), [h5py](https://github.com/h5py/h5py), [netCDF and netCDF4-Python](https://github.com/Unidata/netcdf4-python), [Zarr](https://github.com/zarr-developers/zarr-python), [xarray](https://github.com/pydata/xarray)

**Relationship:** Related storage and data-model roles at different layers.

**Potential shared work:** Common in-memory arrays, named dimensions and adapters; avoid repeated materialization and conversions.

**Keep explicit:** HDF5/Zarr/ASDF/NDF are not interchangeable schemas. Preserve units, WCS, quality and variance rather than only array values.

### 03. Tables, joins, catalogue queries and columnar data

[Astropy](https://github.com/astropy/astropy), [TOPCAT and STILTS](https://www.star.bris.ac.uk/~mbt/topcat/), [Apache Arrow and Parquet](https://arrow.apache.org/), [Polars](https://github.com/pola-rs/polars), [pandas](https://github.com/pandas-dev/pandas), [Vaex](https://github.com/vaexio/vaex), [DuckDB](https://github.com/duckdb/duckdb)

**Relationship:** Overlapping generic table operations with different execution models.

**Potential shared work:** Shared typed catalogue objects, joins, filters and query execution.

**Keep explicit:** Sky-coordinate joins differ from ordinary equijoins; nulls, masks, units and integer precision matter.

### 04. Archive access and local/remote caching

[Astroquery](https://github.com/astropy/astroquery), [PyVO](https://github.com/astropy/pyvo), [IVOA standards](https://www.ivoa.net/documents/), [GWOSC client](https://github.com/gwpy/gwosc), [fsspec](https://github.com/fsspec/filesystem_spec)

**Relationship:** Complementary access layers with overlapping connection/cache concerns.

**Potential shared work:** One resource manager, credentials boundary, download cache and dataset identity scheme.

**Keep explicit:** VO discovery, mission archive APIs and filesystem range requests require separate adapters.

### 05. Units, constants, epochs and positional transformations

[Astropy](https://github.com/astropy/astropy), [ERFA and PyERFA](https://github.com/liberfa/erfa), [IAU SOFA](https://www.iausofa.org/), [Starlink AST](https://starlink.eao.hawaii.edu/starlink/AST), [Skyfield](https://github.com/skyfielders/python-skyfield), [SPICE Toolkit](https://naif.jpl.nasa.gov/naif/toolkit.html), [IDL Astronomy Library](https://github.com/wlandsman/IDLAstro), [SkyCoords.jl](https://github.com/JuliaAstro/SkyCoords.jl), [AstroTime.jl](https://github.com/JuliaAstro/AstroTime.jl), [UnitfulAstro.jl](https://github.com/JuliaAstro/UnitfulAstro.jl), [uom](https://github.com/iliekturtles/uom), [hifitime](https://github.com/nyx-space/hifitime), [ANISE](https://github.com/nyx-space/anise)

**Relationship:** Shared lineage and partially overlapping independent implementations.

**Potential shared work:** Common physical-quantity, epoch, frame and ephemeris services.

**Keep explicit:** SOFA/ERFA ancestry; Earth-orientation tables; leap seconds; observer location; precision and frame conventions.

### 06. WCS, reprojection and mosaicking

[Astropy](https://github.com/astropy/astropy), [GWCS](https://github.com/spacetelescope/gwcs), [WCSLIB](https://www.atnf.csiro.au/people/mcalabre/WCS/), [Starlink AST](https://starlink.eao.hawaii.edu/starlink/AST), [Reproject](https://github.com/astropy/reproject), [SWarp](https://github.com/astromatic/swarp), [Montage](https://ascl.net/1010.036), [DrizzlePac](https://github.com/spacetelescope/drizzlepac)

**Relationship:** Wrappers plus independent reprojection and coaddition implementations.

**Potential shared work:** Reuse coordinate transforms, resampling plans and exposure geometry.

**Keep explicit:** Distortion models, flux versus surface brightness, correlated noise and pixel-area conventions.

### 07. Sky indexing, coverage and spherical transforms

[HEALPix and healpy](https://healpix.sourceforge.io/), [astropy-healpix](https://github.com/astropy/astropy-healpix), [MOCPy and CDS MOC](https://github.com/cds-astro/mocpy), [NaMaster](https://github.com/LSSTDESC/NaMaster), [pixell](https://github.com/simonsobs/pixell), [CDS HEALPix Rust](https://github.com/cds-astro/cds-healpix-rust), [CDS MOC Rust](https://github.com/cds-astro/cds-moc-rust), [ducc](https://github.com/mreineck/ducc), [libsharp](https://github.com/Libsharp/libsharp), [SHTns](https://bitbucket.org/nschaeff/shtns/)

**Relationship:** Related indexing algorithms, some wrappers and different transform engines.

**Potential shared work:** Shared sky geometry and transform backends, with reusable mask/window operators.

**Keep explicit:** HEALPix indexing is not a spherical-harmonic transform; MOC coverage is not a science map; spin, ordering and normalization differ.

### 08. Pipeline orchestration and calibration selection

[ORAC-DR and PICARD](https://github.com/Starlink/ORAC-DR), [Gemini DRAGONS](https://github.com/GeminiDRSoftware/DRAGONS), [stpipe and stcal](https://github.com/spacetelescope/stcal), [CRDS](https://github.com/spacetelescope/crds), [Rubin Science Pipelines](https://pipelines.lsst.io/), [Rubin Butler](https://github.com/lsst/daf_butler), [ESO CPL and EsoRex](https://www.eso.org/sci/software/cpl/), [ESO Reflex](https://www.eso.org/sci/software/esoreflex/), [ESO EDPS and PyCPL](https://www.eso.org/sci/software/pipe_aem_main.html), [Stimela and CARACal](https://github.com/caracal-pipeline/caracal), [Dask](https://github.com/dask/dask), [Snakemake](https://github.com/snakemake/snakemake), [Nextflow](https://github.com/nextflow-io/nextflow)

**Relationship:** Different task systems repeatedly manage dependencies, files and execution.

**Potential shared work:** One task graph, artifact catalogue and cache mechanism; retain domain recipes.

**Keep explicit:** A generic scheduler does not replace instrument calibration rules or scientific data-quality gates.

### 09. Detector correction, masking and cosmic-ray removal

[IRAF Community Distribution](https://iraf-community.github.io/), [STSDAS and TABLES](https://ascl.net/1206.003), [ESO MIDAS](https://www.eso.org/sci/software/esomidas/), [ccdproc](https://github.com/astropy/ccdproc), [Astro-SCRAPPY](https://github.com/astropy/astroscrappy), [JWST pipeline](https://github.com/spacetelescope/jwst), [HSTCAL](https://github.com/spacetelescope/hstcal), [romancal](https://github.com/spacetelescope/romancal), [L.A.Cosmic original](https://ascl.net/1207.005), [CCDPACK](https://starlink.eao.hawaii.edu/docs/sun139.htx/sun139.html)

**Relationship:** Independent tasks plus shared algorithm lineage.

**Potential shared work:** Reusable bias/dark/flat, noise, masking and rejection operations.

**Keep explicit:** Detector models, gain, saturation, ramp sampling and propagated uncertainty are instrument-dependent.

### 10. Detection, segmentation and deblending

[Photutils](https://github.com/astropy/photutils), [Source Extractor](https://github.com/astromatic/sextractor), [SEP](https://github.com/sep-developers/sep), [scarlet](https://github.com/pmelchior/scarlet), [PyBDSF](https://github.com/lofar-astron/PyBDSF), [Aegean](https://github.com/PaulHancock/Aegean), [SoFiA 2](https://github.com/SoFiA-Admin/SoFiA-2), [ProFound](https://github.com/asgr/ProFound), [CUPID](https://starlink.eao.hawaii.edu/docs/sun255.htx/sun255.html)

**Relationship:** Shared ancestry in some pairs and distinct methods elsewhere.

**Potential shared work:** Common image/cube detection primitives, segmentation objects and catalogue outputs.

**Keep explicit:** SEP derives from Source Extractor; scarlet-style source models, radio detection and cube clumps are not identical methods.

### 11. PSF, aperture and crowded-field photometry

[Photutils](https://github.com/astropy/photutils), [PSFEx](https://github.com/astromatic/psfex), [DAOPHOT and ALLSTAR](https://ascl.net/1104.011), [DOLPHOT](https://ascl.net/1608.013), [DoPHOT](https://github.com/AbhijitSaha/DoPhot), [The Tractor](https://github.com/dstndstn/tractor), [PSFModels.jl](https://github.com/JuliaAstro/PSFModels.jl), [Photometry.jl](https://github.com/JuliaAstro/Photometry.jl), [VaST](https://github.com/kirxkirx/vast)

**Relationship:** Partially overlapping measurement tasks.

**Potential shared work:** Shared PSF objects, aperture integrals, fitting primitives and uncertainty outputs.

**Keep explicit:** Spatial PSFs, crowding, multiband constraints, forced positions and detector calibration change the inference.

### 12. Galaxy image modelling and morphology

[GALFIT](https://users.obs.carnegiescience.edu/peng/work/galfit/galfit.html), [Imfit](https://github.com/perwin/imfit), [The Tractor](https://github.com/dstndstn/tractor), [scarlet](https://github.com/pmelchior/scarlet), [statmorph](https://github.com/vrodgom/statmorph), [PyAutoLens](https://github.com/Jammy2211/PyAutoLens), [ProFit](https://github.com/ICRAR/ProFit)

**Relationship:** Overlapping image forward models; different reconstruction and fitting methods.

**Potential shared work:** Common image formation, PSF convolution, source profiles and fit-data contracts.

**Keep explicit:** Parametric structure, nonparametric morphology, deblending and lens-source reconstruction remain distinct scientific operations.

### 13. Image subtraction and time-series detection

[HOTPANTS](https://github.com/acbecker/hotpants), [ISIS image subtraction](https://www.iap.fr/useriap/alard/package.html), [Rubin Science Pipelines](https://pipelines.lsst.io/), [VaST](https://github.com/kirxkirx/vast)

**Relationship:** Similar workflow goals and related difference-imaging algorithms.

**Potential shared work:** Shared registration, convolution, difference-image uncertainty and detection products.

**Keep explicit:** Kernel basis, noise model, reference-image construction and detection selection must remain explicit.

### 14. Spectral reduction, extraction and line fitting

[specutils](https://github.com/astropy/specutils), [PypeIt](https://github.com/pypeit/PypeIt), [Gemini DRAGONS](https://github.com/GeminiDRSoftware/DRAGONS), [DESI spectroscopic pipeline](https://github.com/desihub/desispec), [SDSS idlspec2d](https://github.com/sdss/idlspec2d), [MaNGA DAP](https://github.com/sdss/mangadap), [MPDAF](https://github.com/musevlt/mpdaf), [pPXF](https://users.physics.ox.ac.uk/~cappellari/software/), [GANDALF](https://ascl.net/1708.012), [linetools](https://github.com/linetools/linetools), [SpecPro](https://ascl.net/1404.014), [FIGARO and DIPSO](https://starlink.eao.hawaii.edu/docs/sun86.htx/sun86.html)

**Relationship:** Repeated generic primitives inside instrument and survey workflows.

**Potential shared work:** Common spectrum objects, wavelength grids, masks, response functions and extraction/fitting interfaces.

**Keep explicit:** Pixel covariance, instrumental line spread, sky subtraction and redshift/template assumptions.

### 15. Synthetic photometry, extinction and SED transformations

[synphot](https://github.com/spacetelescope/synphot_refactor), [dust_extinction](https://github.com/karllark/dust_extinction), [dustmaps](https://github.com/gregreen/dustmaps), [SNCosmo](https://github.com/sncosmo/sncosmo), [FSPS and python-fsps](https://github.com/cconroy20/fsps), [Prospector](https://github.com/bd-j/prospector), [Bagpipes](https://github.com/ACCarnall/bagpipes), [CIGALE](https://cigale.lam.fr/), [MAGPHYS](https://www.iap.fr/magphys/), [LePhare](https://github.com/lephare-photoz/lephare), [EAZY](https://github.com/gbrammer/eazy-py), [kcorrect](https://github.com/blanton144/kcorrect), [ProSpect](https://github.com/asgr/ProSpect)

**Relationship:** Repeated passband integration and attenuation operations across different physical models.

**Potential shared work:** Share passbands, extinction laws, interpolation and magnitude/flux conversions.

**Keep explicit:** Rest/observer frame, photon/energy weighting, zero-points, wavelength ranges and model validity.

### 16. High-energy event filtering, response folding and spectra

[HEASoft and FTOOLS](https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/), [XSPEC and PyXspec](https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/), [XSELECT](https://heasarc.gsfc.nasa.gov/docs/software/ftools/xselect/), [CIAO](https://cxc.cfa.harvard.edu/ciao/), [Sherpa](https://github.com/sherpa/sherpa), [XMM Newton SAS](https://www.cosmos.esa.int/web/xmm-newton/sas), [SPEX](https://www.sron.nl/en/pillars/science/astrophysics/), [ISIS spectral analysis](https://space.mit.edu/cxc/isis/), [Gammapy](https://github.com/gammapy/gammapy), [ctools and GammaLib](https://cta.irap.omp.eu/ctools/), [Fermitools](https://github.com/fermi-lat/Fermitools-conda), [threeML](https://github.com/threeML/threeML)

**Relationship:** Overlapping analysis operations in mission and modelling ecosystems.

**Potential shared work:** Shared event lists, good-time intervals, responses, likelihood interfaces and model evaluation.

**Keep explicit:** Response definitions, Poisson backgrounds, plasma models and mission calibration are not automatically interchangeable.

### 17. Timing, period searches and stochastic signals

[Astropy](https://github.com/astropy/astropy), [XRONOS](https://heasarc.gsfc.nasa.gov/docs/xanadu/xronos/), [Stingray](https://github.com/StingraySoftware/stingray), [HENDRICS](https://github.com/StingraySoftware/HENDRICS), [PRESTO](https://github.com/scottransom/presto), [TEMPO](https://tempo.sourceforge.net/), [TEMPO2](https://bitbucket.org/psrsoft/tempo2/), [PINT pulsar timing](https://github.com/nanograv/PINT), [Lightkurve](https://github.com/lightkurve/lightkurve), [celerite and tinygp](https://github.com/dfm/celerite2), [VARTOOLS](https://www.astro.princeton.edu/~jhartman/vartools.html)

**Relationship:** Related signal-analysis primitives with domain-specific timing models.

**Potential shared work:** Common time-series, Fourier, windowing, period-search and covariance operations.

**Keep explicit:** Pulsar barycentric timing is more than a periodogram; irregular sampling, exposure and noise models differ.

### 18. Radio calibration and visibility imaging

[CASA](https://casadocs.readthedocs.io/en/stable/), [casacore](https://github.com/casacore/casacore), [AIPS](https://www.aips.nrao.edu/), [MIRIAD](https://www.atnf.csiro.au/computing/software/miriad/), [GILDAS](https://www.iram.fr/IRAMFR/GILDAS/), [WSClean](https://wsclean.readthedocs.io/), [DP3](https://dp3.readthedocs.io/), [DDFacet](https://github.com/saopicc/DDFacet), [killMS](https://github.com/saopicc/killMS), [CubiCal](https://github.com/ratt-ru/CubiCal), [QuartiCal](https://github.com/ratt-ru/QuartiCal), [MeqTrees](https://github.com/ratt-ru/meqtrees-timba), [RASCIL and SKA SDP components](https://developer.skao.int/projects/rascil/en/latest/), [Stimela and CARACal](https://github.com/caracal-pipeline/caracal)

**Relationship:** Alternative pipelines and solvers around related measurement equations.

**Potential shared work:** One visibility data model, Jones operators, flagging interfaces, gridding and calibration products.

**Keep explicit:** Direction-dependent gains, beam models, weighting, wide-field approximations and solver degeneracies.

### 19. VLBI and optical-interferometric reconstruction

[eht-imaging](https://github.com/achael/eht-imaging), [SMILI](https://github.com/astrosmili/smili), [Difmap](https://sites.astro.caltech.edu/~mcs/difmap/), [VLBISkyModels.jl](https://github.com/EHTJulia/VLBISkyModels.jl), [Comrade.jl](https://github.com/ptiede/Comrade.jl), [PMOIRED](https://github.com/amerand/PMOIRED), [MiRA](https://github.com/emmt/MiRA), [SQUEEZE](https://github.com/fabienbaron/squeeze), [LITpro](https://www.jmmc.fr/english/tools/proposal-preparation/litpro/), [OIFITS and OIFITSlib](https://www.jmmc.fr/oifits/)

**Relationship:** Related inverse problems with different observables and regularizers.

**Potential shared work:** Shared Fourier-domain operators, closure quantities and reconstruction interfaces.

**Keep explicit:** Radio visibilities, optical squared visibilities and closure phases carry different calibration and likelihood information.

### 20. Strong-lens potentials, ray tracing and source inference

[lenstronomy](https://github.com/lenstronomy/lenstronomy), [PyAutoLens](https://github.com/Jammy2211/PyAutoLens), [Herculens](https://github.com/Herculens/herculens), [caustics](https://github.com/Ciela-Institute/caustics), [jaxtronomy](https://github.com/lenstronomy/JAXtronomy), [Lenstool](https://ascl.net/1102.004), [glafic](https://github.com/oguri/glafic2), [gravlens and lensmodel](https://www.physics.rutgers.edu/~keeton/gravlens/), [GLEE](https://arxiv.org/abs/2209.03094), [GLASS](https://github.com/jpcoles/glass), [WSLAP plus](https://arxiv.org/abs/1304.2393)

**Relationship:** Independent implementations, ports and competing reconstruction methods.

**Potential shared work:** Common lens potentials, deflections, Jacobians, multiple-plane geometry and image formation.

**Keep explicit:** Source bases, priors, mass-sheet degeneracy, boundary conditions and supported mass models.

### 21. Lensing time delays, microlensing and coherent waves

[PyCS3](https://github.com/COSMOGRAIL/PyCS), [hierArc](https://github.com/TDCOSMO/hierArc), [SNTD](https://github.com/jpierel14/sntd), [MulensModel](https://github.com/rpoleski/MulensModel), [pyLIMA](https://github.com/ebachelet/pyLIMA), [VBMicrolensing](https://github.com/valboz/VBMicrolensing), [GLoW](https://github.com/glow-astro/GLoW)

**Relationship:** Connected problems; not one interchangeable solver family.

**Potential shared work:** Reuse distances, trajectories, source profiles, lens geometry and likelihood plumbing.

**Keep explicit:** Time-delay cosmography, finite-source microlensing and wave-optics diffraction require different numerical engines.

### 22. Weak-lensing shapes, simulation and calibration

[GalSim](https://github.com/GalSim-developers/GalSim), [ngmix](https://github.com/esheldon/ngmix), [metadetect](https://github.com/esheldon/metadetect), [TreeCorr](https://github.com/rmjarvis/TreeCorr), [LensTools](https://github.com/apetri/LensTools), [JAX-GalSim](https://github.com/GalSim-developers/JAX-GalSim)

**Relationship:** Shared lineage in some tools and distinct measurement/calibration stages.

**Potential shared work:** Common PSF/image operators, shear conventions, catalogues, masks and correlation primitives.

**Keep explicit:** Metacalibration, detection effects and selection response are part of the estimator; JAX-GalSim is labelled early development.

### 23. Background cosmology, distances, ages and growth

[Astropy](https://github.com/astropy/astropy), [CAMB](https://github.com/cmbant/CAMB), [CLASS cosmology](https://github.com/lesgourg/class_public), [CCL and pyccl](https://github.com/LSSTDESC/CCL), [Colossus](https://bitbucket.org/bdiemer/colossus/), [jax-cosmo](https://github.com/DifferentiableUniverseInitiative/jax_cosmo), [Cosmology.jl](https://github.com/JuliaAstro/Cosmology.jl), [Bolt.jl](https://github.com/xzackli/Bolt.jl), [SymBoltz.jl](https://github.com/hersle/SymBoltz.jl)

**Relationship:** Many packages expose overlapping background calculations.

**Potential shared work:** One versioned background solution and interpolation service per cosmological model and parameter state.

**Keep explicit:** Neutrinos, radiation, curvature, dark energy and numerical tolerances must be in the cache key and interface.

### 24. Einstein–Boltzmann, power-spectrum and perturbation calculations

[CAMB](https://github.com/cmbant/CAMB), [CLASS cosmology](https://github.com/lesgourg/class_public), [CCL and pyccl](https://github.com/LSSTDESC/CCL), [jax-cosmo](https://github.com/DifferentiableUniverseInitiative/jax_cosmo), [CosmoPower and CosmoPower-JAX](https://github.com/alessiospuriomancini/cosmopower), [CLASS-PT](https://github.com/Michalychforever/CLASS-PT), [PyBird](https://github.com/pierrexyz/pybird), [velocileptors](https://github.com/sfschen/velocileptors), [Bolt.jl](https://github.com/xzackli/Bolt.jl), [SymBoltz.jl](https://github.com/hersle/SymBoltz.jl), [FFTLog](https://jila.colorado.edu/~ajsh/FFTLog/), [mcfit](https://github.com/eelregit/mcfit), [FAST-PT](https://github.com/JoeMcEwen/FAST-PT)

**Relationship:** Independent solvers, ports, approximations and trained surrogates.

**Potential shared work:** Reuse background histories, transfer functions, transform kernels and requested observable grids.

**Keep explicit:** Linear/nonlinear regimes, EFT terms, gauge, neutrinos, normalization and emulator training support. A surrogate is not independent validation of its training solver.

### 25. Correlations, power estimators and survey windows

[TreeCorr](https://github.com/rmjarvis/TreeCorr), [Corrfunc](https://github.com/manodeep/Corrfunc), [nbodykit](https://github.com/bccp/nbodykit), [pycorr and pypower](https://github.com/cosmodesi/pycorr), [NaMaster](https://github.com/LSSTDESC/NaMaster), [pixell](https://github.com/simonsobs/pixell)

**Relationship:** Different estimators with overlapping pair-count and Fourier primitives.

**Potential shared work:** Shared catalogue selections, pair counts, transforms, binning and window operators.

**Keep explicit:** Estimator normalization, random catalogues, masks, shot noise, weights and covariance conventions.

### 26. Cosmological likelihood and parameter pipelines

[Cobaya](https://github.com/CobayaSampler/cobaya), [CosmoMC](https://github.com/cmbant/CosmoMC), [MontePython](https://github.com/brinckmann/montepython_public), [CosmoSIS](https://github.com/joezuntz/cosmosis), [SMICA and Planck likelihood code](https://pla.esac.esa.int/), [21CMMC](https://github.com/21cmfast/21CMMC), [GetDist](https://github.com/cmbant/getdist)

**Relationship:** Different orchestration systems, sometimes calling the same theory engines or datasets.

**Potential shared work:** Common theory-product cache, parameter graph, likelihood interfaces and chain interchange.

**Keep explicit:** Shared data and calibration create covariance; identical theory engines are not independent scientific checks.

### 27. Supernova and transient modelling

[SNANA](https://github.com/RickKessler/SNANA), [SNCosmo](https://github.com/sncosmo/sncosmo), [BayeSN](https://github.com/bayesn/bayesn), [SALTShaker](https://github.com/djones1040/SALTShaker), [SNooPy](https://github.com/obscode/snpy), [SNID](https://people.lam.fr/blondin.stephane/software/snid/), [Superfit](https://github.com/dahowell/superfit), [SuperNNova](https://github.com/supernnova/SuperNNova), [MOSFiT](https://github.com/guillochon/MOSFiT), [Redback](https://github.com/nikhil-sarin/redback)

**Relationship:** Partly overlapping light-curve, spectral, training and inference work.

**Potential shared work:** Shared observations, filters, redshifts, time grids, model adapters and selection simulations.

**Keep explicit:** SALT-family fitting, BayeSN, template classification and physical transients are different models; preserve training and calibration lineage.

### 28. Stellar populations and galaxy SED inference

[MESA](https://github.com/MESAHub/mesa), [SSE and BSE](https://astronomy.swin.edu.au/~jhurley/), [COSMIC](https://github.com/COSMIC-PopSynth/COSMIC), [COMPAS](https://github.com/TeamCOMPAS/COMPAS), [POSYDON](https://github.com/POSYDON-code/POSYDON), [binary_c](https://ascl.net/2307.035), [BPASS](https://bpass.auckland.ac.nz/), [FSPS and python-fsps](https://github.com/cconroy20/fsps), [Starburst99](https://www.stsci.edu/science/starburst99/docs/default.htm), [Prospector](https://github.com/bd-j/prospector), [Bagpipes](https://github.com/ACCarnall/bagpipes), [CIGALE](https://cigale.lam.fr/), [MAGPHYS](https://www.iap.fr/magphys/), [ProSpect](https://github.com/asgr/ProSpect)

**Relationship:** Nested physical layers with overlapping fitting and population interfaces.

**Potential shared work:** Reusable stellar libraries, population integration, dust/emission components and observation projection.

**Keep explicit:** Stellar evolution, binary populations, synthesis and SED fitting are distinct levels; IMF and evolutionary tracks remain explicit.

### 29. Radiative transfer, ionization and spectral microphysics

[TARDIS](https://github.com/tardis-sn/tardis), [SEDONA](https://github.com/dnkasen/pubsed), [MOOG](https://www.as.utexas.edu/~chris/moog.html), [Turbospectrum](https://github.com/bertrandplez/Turbospectrum2019), [SME and PySME](https://github.com/AWehrhahn/SME), [ATLAS and SYNTHE](https://wwwuser.oats.inaf.it/castelli/sources.html), [TLUSTY and SYNSPEC](https://tlusty.oca.eu/), [PHOENIX](https://ascl.net/1010.056), [Cloudy](https://gitlab.nublado.org/cloudy/cloudy), [XSTAR](https://heasarc.gsfc.nasa.gov/lheasoft/xstar/xstar.html), [CHIANTI and ChiantiPy](https://www.chiantidatabase.org/), [AtomDB and PyAtomDB](https://github.com/AtomDB/pyatomdb), [PyNeb](https://github.com/Morisset/PyNeb_devel), [RADMC-3D](https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/), [MCFOST](https://github.com/cpinte/mcfost), [SKIRT](https://github.com/SKIRT/SKIRT9), [Hyperion](https://github.com/hyperion-rt/hyperion), [LIME](https://github.com/lime-rt/lime), [RADEX](https://ascl.net/1010.075), [MOCASSIN](https://mocassin.nebulousresearch.org/), [TORUS](https://ascl.net/1404.006), [RH](https://github.com/ITA-Solar/rh)

**Relationship:** Related physics with different regimes and numerical formulations.

**Potential shared work:** Share atomic/molecular data access, opacity/emissivity objects and transport interfaces.

**Keep explicit:** LTE/non-LTE, scattering, polarization, velocity fields, geometry and equilibrium assumptions define separate solvers.

### 30. N-body gravity, hydro and magnetohydrodynamics

[AMUSE](https://github.com/amusecode/amuse), [GADGET-4](https://wwwmpa.mpa-garching.mpg.de/gadget4/), [GADGET-2](https://wwwmpa.mpa-garching.mpg.de/gadget/), [AREPO](https://arepo-code.org/), [GIZMO](https://ascl.net/1410.003), [SWIFT](https://github.com/SWIFTSIM/SWIFT), [RAMSES](https://github.com/ramses-organisation/ramses), [Enzo](https://github.com/enzo-project/enzo-dev), [Enzo-E](https://github.com/enzo-project/enzo-e), [FLASH](https://flash.rochester.edu/site/flashcode/), [Athena and Athena++](https://github.com/PrincetonUniversity/athena), [AthenaK](https://github.com/IAS-Astrophysics/athenak), [PLUTO](https://plutocode.ph.unito.it/), [MPI-AMRVAC](https://github.com/amrvac/amrvac), [Pencil Code](https://github.com/pencil-code/pencil-code), [Phantom](https://github.com/danieljprice/phantom), [ChaNGa](https://github.com/N-BodyShop/changa), [PKDGRAV3](https://github.com/CLS-project/PKDGRAV3), [REBOUND](https://github.com/hannorein/rebound), [NBODY6 and NBODY6++GPU](https://github.com/nbodyx/Nbody6ppGPU), [PeTar](https://github.com/lwang-astro/PeTar), [AMReX](https://github.com/AMReX-Codes/amrex), [Castro](https://github.com/AMReX-Astro/Castro), [MAESTROeX](https://github.com/AMReX-Astro/MAESTROeX), [Trixi.jl](https://github.com/trixi-framework/Trixi.jl)

**Relationship:** Independent numerical methods with shared mesh/particle, time-step and communication needs.

**Potential shared work:** Common particle/mesh objects, initial conditions, gravity and microphysics interfaces, checkpointing and diagnostics.

**Keep explicit:** SPH, meshless, moving mesh, AMR, direct N-body and low-Mach methods must not be collapsed into one undocumented approximation.

### 31. Simulation readers, halo finding and dynamical analysis

[galpy](https://github.com/jobovy/galpy), [gala](https://github.com/adrn/gala), [AGAMA](https://github.com/GalacticDynamics-Oxford/Agama), [NEMO](https://github.com/teuben/nemo), [yt](https://github.com/yt-project/yt), [pynbody](https://github.com/pynbody/pynbody), [MUSIC](https://www.oca.eu/en/olivier-hahn/1914-music), [ROCKSTAR](https://bitbucket.org/gfcstanford/rockstar/), [AHF](https://ascl.net/1102.009), [VELOCIraptor](https://github.com/pelahi/VELOCIraptor-STF)

**Relationship:** Analysis layers overlap on readers, units, selections and particle data.

**Potential shared work:** One snapshot reader layer, spatial indexing, units and output catalogue format.

**Keep explicit:** Halo definitions, unbinding, merger-tree construction and dynamics models can change scientific results.

### 32. Symbolic geometry, geodesics and numerical relativity

[Einstein Toolkit](https://einsteintoolkit.org/), [Cactus and CarpetX](https://github.com/eschnett/CarpetX), [SpECTRE](https://github.com/sxs-collaboration/spectre), [SpEC](https://www.black-holes.org/code/SpEC.html), [GRChombo](https://github.com/GRChombo/GRChombo), [NRPy](https://github.com/nrpy/nrpy), [EinsteinPy](https://github.com/einsteinpy/einsteinpy), [xAct](https://www.xact.es/), [Cadabra](https://cadabra.science/), [SageManifolds](https://sagemanifolds.obspm.fr/), [GRTensorIII](https://github.com/grtensor/grtensor), [Black Hole Perturbation Toolkit](https://bhptoolkit.org/), [GYOTO](https://github.com/gyoto/Gyoto), [ipole](https://github.com/moscibrodzka/ipole), [grtrans](https://github.com/jadexter/grtrans), [BHAC](https://bhac.science/), [HARM and iharm3D](https://github.com/AFD-Illinois/iharm3d), [KHARMA](https://github.com/AFD-Illinois/kharma), [Gradus.jl](https://github.com/astro-group-bristol/Gradus.jl), [Krang.jl](https://github.com/dchang10/Krang.jl), [Skylight.jl](https://github.com/joaquinpelle/Skylight.jl), [GRTeclyn](https://github.com/GRTLCollaboration/GRTeclyn)

**Relationship:** Several physically related but numerically distinct levels.

**Potential shared work:** Share tensor conventions, metric definitions, geometric primitives and diagnostics where valid.

**Keep explicit:** Symbolic equations, fixed-spacetime ray tracing and evolving Einstein equations are different jobs; gauge and formulation must remain visible.

### 33. Gravitational-wave signals, likelihoods and populations

[LALSuite](https://git.ligo.org/lscsoft/lalsuite), [PyCBC](https://github.com/gwastro/pycbc), [Bilby](https://github.com/bilby-dev/bilby), [GWpy](https://github.com/gwpy/gwpy), [GstLAL](https://git.ligo.org/lscsoft/gstlal), [cWB](https://gwburst.gitlab.io/), [BayesWave](https://git.ligo.org/lscsoft/bayeswave), [RIFT](https://git.ligo.org/rapidpe-rift/rift), [PESummary](https://github.com/pesummary/pesummary), [ligo.skymap](https://git.ligo.org/lscsoft/ligo.skymap), [gwpopulation](https://github.com/ColmTalbot/gwpopulation), [gwsurrogate](https://github.com/sxs-collaboration/gwsurrogate), [PyRing](https://git.ligo.org/lscsoft/pyring)

**Relationship:** Shared waveform libraries plus alternative search and inference pipelines.

**Potential shared work:** Common strain/PSD data, waveform calls, detector responses, likelihood products and posterior formats.

**Keep explicit:** Search detection, event inference and population selection are different levels; calibration and waveform approximations matter.

### 34. Generic optimization, sampling and posterior summaries

[Sherpa](https://github.com/sherpa/sherpa), [Cobaya](https://github.com/CobayaSampler/cobaya), [Bilby](https://github.com/bilby-dev/bilby), [emcee](https://github.com/dfm/emcee), [dynesty](https://github.com/joshspeagle/dynesty), [UltraNest](https://github.com/JohannesBuchner/UltraNest), [MultiNest and PyMultiNest](https://github.com/farhanferoz/MultiNest), [PolyChord](https://github.com/PolyChord/PolyChordLite), [NumPyro](https://github.com/pyro-ppl/numpyro), [PyMC](https://github.com/pymc-devs/pymc), [Stan](https://mc-stan.org/), [iminuit and Minuit2](https://github.com/scikit-hep/iminuit), [MPFIT](https://ascl.net/1208.019), [GNU Scientific Library](https://www.gnu.org/software/gsl/), [MINPACK](https://www.netlib.org/minpack/), [ArviZ](https://github.com/arviz-devs/arviz), [ChainConsumer](https://github.com/Samreay/ChainConsumer)

**Relationship:** Some shared optimizer ancestry, many independent algorithms and repeated adapters.

**Potential shared work:** Common parameter transforms, likelihood/gradient interfaces, diagnostics and sample schemas.

**Keep explicit:** Least squares, MCMC, nested sampling, variational inference and optimization are not equivalent statistical procedures.

### 35. Dense/sparse linear algebra, FFTs, quadrature and ODEs

[NumPy and SciPy](https://scipy.org/), [Diffrax](https://github.com/patrick-kidger/diffrax), [SciML DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl), [Dedalus](https://github.com/DedalusProject/dedalus), [FEniCSx](https://fenicsproject.org/), [Firedrake](https://www.firedrakeproject.org/), [deal.II](https://www.dealii.org/), [MFEM](https://mfem.org/), [PETSc](https://petsc.org/), [SUNDIALS](https://sundials.readthedocs.io/), [ndarray](https://github.com/rust-ndarray/ndarray), [nalgebra](https://github.com/dimforge/nalgebra), [faer](https://github.com/sarah-ek/faer-rs), [RustFFT](https://github.com/ejmahler/RustFFT), [argmin](https://github.com/argmin-rs/argmin), [diffsol](https://github.com/martinjrobins/diffsol), [FFTW](https://www.fftw.org/), [ducc](https://github.com/mreineck/ducc), [BLAS and LAPACK](https://www.netlib.org/lapack/), [OpenBLAS](https://github.com/OpenMathLib/OpenBLAS), [SuiteSparse](https://github.com/DrTimothyAldenDavis/SuiteSparse), [GNU Scientific Library](https://www.gnu.org/software/gsl/), [MINPACK](https://www.netlib.org/minpack/), [ODEPACK](https://www.netlib.org/odepack/), [QUADPACK](https://www.netlib.org/quadpack/), [FITPACK and DIERCKX](https://www.netlib.org/dierckx/), [FFTPACK](https://www.netlib.org/fftpack/), [Cuba](https://feynarts.de/cuba/)

**Relationship:** Frequently shared underlying engines or historical algorithm families beneath many frontends.

**Potential shared work:** Central numerical services, memory layout, thread pools, transform/factorization plans and tolerance conventions.

**Keep explicit:** Backend choice, conditioning, precision, stiffness, convergence criteria and deterministic execution; common ancestry reduces independence.

### 36. Visualization and analysis applications

[PDL](https://pdl.perl.org/), [SAOImage DS9](https://ds9.si.edu/site/Home.html), [TOPCAT and STILTS](https://www.star.bris.ac.uk/~mbt/topcat/), [Aladin Desktop and Aladin Lite](https://aladin.cds.unistra.fr/), [Glue](https://github.com/glue-viz/glue), [Ginga](https://github.com/ejeschke/ginga), [Jdaviz](https://github.com/spacetelescope/jdaviz), [CARTA](https://cartavis.org/), [VisIt](https://visit-dav.github.io/visit-website/), [ParaView and VTK](https://www.paraview.org/), [OpenSpace](https://www.openspaceproject.com/), [PGPLOT](https://sites.astro.caltech.edu/~tjp/pgplot/), [SM SuperMongo](https://www.astro.princeton.edu/~rhl/sm/)

**Relationship:** Repeated readers, selections, plotting and interactive inspection across applications.

**Potential shared work:** One linked observation/catalogue viewer backed by shared objects, selections and provenance.

**Keep explicit:** Preserve image, spectrum, cube, sky-map, event and simulation-specific interactions rather than forcing a single generic plot.

## Repetition to inspect in actual workflows

The following are concrete places to look for repeated execution. They are hypotheses for a future workflow audit, not measured findings from this survey.

| Repeated work | What could be reused | What makes reuse valid |
|---|---|---|
| Every application reads and converts the same file | Decoded, typed observation objects and lazy views | Source bytes, schema, selected HDUs/columns, scaling, masks and units |
| Multiple packages recalculate distances or expansion histories | A cached cosmological background solution | Full cosmological parameters, physical assumptions, constants and tolerance |
| Several likelihoods separately call the same Boltzmann solver | A theory-product dependency graph | Model state, requested observables, grids, nonlinear prescription and solver settings |
| Each fitter independently interpolates the same templates and filters | Reusable template grids and passband operators | Model version, calibration, wavelength frame, integration convention and interpolation |
| PSF convolutions or FFT plans are repeatedly rebuilt | Prepared operators and transform plans | Kernel/image geometry, boundary conditions, normalization, precision and device |
| The same covariance is factorized for every likelihood call | Cached factorization or linear operator | Covariance truly parameter-independent; same ordering, masks, normalization and numerical method |
| Data are copied between C/Fortran/Python/Julia/Rust containers | Shared memory buffers and explicit ownership | Compatible strides, types, alignment, lifetime, masks, units and device placement |
| Survey masks, pair counts and window operators are rebuilt | Versioned geometry and selection products | Identical objects, weights, cuts, estimator definitions and random catalogues |
| Each pipeline downloads or parses calibration assets separately | A shared versioned calibration/data cache | Exact release, validity interval, instrument configuration and checksum |
| Every application converts and plots posterior chains differently | Common weighted-sample and diagnostic schema | Correct weights, transforms, log-density conventions, priors and sampler metadata |

## Consolidation conclusion

The greatest broadly shared surface is below the individual scientific models: data access, scientific objects, calibration assets, geometry, numerical primitives, task execution and inference interfaces. Above it, a unified software unit can expose multiple named model and solver implementations. This combines the repeated plumbing while preserving alternative algorithms and the ability to compare them.

No speedup estimate is claimed. A library dependency graph, representative end-to-end workflows and profiling would be needed to establish how much computation, memory movement or I/O is actually duplicated.

See the [full catalogue](CATALOGUE.md), [overall analysis](REPORT.md), [source evidence](EVIDENCE.md), and [machine-readable overlap map](overlap-map.json).
