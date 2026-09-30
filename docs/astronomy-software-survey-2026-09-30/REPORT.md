# Cosmology, astronomy and computational-physics software survey

The follow-on [light design for an agent-operated engine](LIGHT-DESIGN.md) applies the later scope decisions: a thin CLI, shared compiled physics, optional data adapters, no visualization layer, and C++/CUDA as the proposed implementation choice. It supersedes the broad architecture sketch below for the intended product; the survey itself remains a reference inventory.

This survey follows the full cosmology workflow: obtaining and loading data, representing observations, calibrating instruments, measuring signals, calculating physical predictions, running simulations, fitting models and inspecting results. It identifies capabilities spread across applications and packages, including repeated infrastructure and overlapping scientific calculations, to inform a possible unified Rust suite. It includes modern packages, observatory systems, historical programs and specialist research codes. Legacy software remains in scope even when it is difficult to install, no longer maintained, or available only as documentation or binaries. Project size is not used as an exclusion criterion.

The result has two levels: a curated catalogue of 528 entries across 26 groups, and a separate title-and-link index of 4,260 records from the Astrophysics Source Code Library (4,105 registered and 155 submitted). The curated catalogue includes individual packages, related software families, umbrella suites, and a few essential standards and data services. Its count is therefore **not a count of independent engines**. The ASCL index is a breadth and discovery resource; its thousands of records have not each received the same review as the curated list. The collection date is 30 September 2026, Europe/London.

Start with [the overlap map](OVERLAP-MAP.md) for 36 families of repeated capabilities across 376 curated identities. It distinguishes shared engines, related implementations, alternative methods and potentially repeated pipeline computation. Open [the searchable catalogue](catalogue.html) to filter both collections. [CATALOGUE.md](CATALOGUE.md) contains the full curated list in reading order. [catalogue.json](catalogue.json) and [ascl-index.json](ascl-index.json) are machine-readable. [EVIDENCE.md](EVIDENCE.md) records source retrieval limits and validation counts.

The main architectural conclusion is a proposal drawn from the survey: one coherent product would benefit from a common system for physical quantities, observations, models, solvers, uncertainty, calibration, and provenance, with specialist engines above it. Astropy, AMUSE, CASA, HEASoft, Starlink, and SciML each provide important pieces of this design. None is a complete substitute for the others. This is a capability inventory and design synthesis, not an implementation or a measured performance comparison.

## The principal ecosystems to study

| Ecosystem | What it contributes to the prospective suite | Source |
|---|---|---|
| Astropy and its coordinated and affiliated packages | Common scientific objects, astronomy conventions, file access and interoperable analysis libraries | [Astropy](https://github.com/astropy/astropy), [ecosystem registry](https://www.astropy.org/affiliated/) |
| Starlink | A broad task system with mature image, spectral, coordinate, quality, variance and pipeline conventions | [Starlink collection](https://github.com/Starlink/starlink) |
| IRAF and ESO MIDAS | Historical reduction algorithms, task semantics, interactive workflows and instrument expertise | [IRAF community](https://iraf-community.github.io/), [MIDAS](https://www.eso.org/sci/software/esomidas/) |
| CASA and casacore | Radio measurement equations, calibration, visibilities, imaging and radio data structures | [CASA documentation](https://casadocs.readthedocs.io/en/stable/), [casacore](https://github.com/casacore/casacore) |
| HEASoft, CIAO and SAS | Event-based observations, response matrices, spectral/timing analysis and mission calibration | [HEASoft](https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/), [CIAO](https://cxc.cfa.harvard.edu/ciao/), [SAS](https://www.cosmos.esa.int/web/xmm-newton/sas) |
| SunPy, SolarSoft and SPEDAS | Solar and space-physics coordinates, mission readers, time series and instrument analysis | [SunPy](https://github.com/sunpy/sunpy), [SolarSoft](https://www.lmsal.com/solarsoft/), [SPEDAS](https://spedas.org/) |
| AMUSE | Coupling independently implemented physical solvers through common quantities, particles, grids and interfaces | [AMUSE](https://github.com/amusecode/amuse), [architecture paper](https://arxiv.org/abs/1307.3016) |
| Einstein Toolkit and related numerical-relativity codes | Initial data, spacetime evolution, matter coupling, horizons and waveform diagnostics | [Einstein Toolkit](https://einsteintoolkit.org/), [SpECTRE](https://github.com/sxs-collaboration/spectre), [GRTeclyn](https://github.com/GRTLCollaboration/GRTeclyn) |
| SciML, Dedalus and finite-element ecosystems | Equation specification, discretization, solver composition, sensitivities and execution | [SciML](https://github.com/SciML/DifferentialEquations.jl), [Dedalus](https://github.com/DedalusProject/dedalus), [FEniCS](https://fenicsproject.org/) |
| Rubin and STScI pipelines | Task graphs, instrument data models, calibration reference management and reproducible processing | [Rubin pipelines](https://pipelines.lsst.io/), [Butler](https://github.com/lsst/daf_butler), [CRDS](https://github.com/spacetelescope/crds) |

Astropy's official registry distinguishes coordinated packages from affiliated packages: they do not all have the same governance or release cycle. The registry is also undergoing a review-process transition, so its displayed affiliated-package list should not be mistaken for the complete astronomy Python ecosystem. The catalogue consequently includes relevant packages regardless of formal affiliation. [Astropy registry](https://www.astropy.org/affiliated/)

## Astropy capability analysis

Astropy is a necessary reference for the proposed common layer, but its name hides many separate contracts. A Rust counterpart would need explicit behaviour for each of these areas. This matrix is a proposed decomposition of documented Astropy capabilities, not a claim that every item should be a one-to-one port. [Astropy source and documentation](https://github.com/astropy/astropy)

| Capability area | Semantics to preserve or explicitly redesign | Related implementations |
|---|---|---|
| Units and physical quantities | Dimensional algebra, conversions, equivalencies, logarithmic quantities, magnitude systems and array behaviour | Astropy units; Unitful/UnitfulAstro; Rust uom |
| Physical constants | Named reference sets, uncertainties and version identity rather than silently changing constants | Astropy constants; domain-specific physical tables |
| Time and epochs | UTC, TAI, TT, TDB and UT1; leap seconds; Earth orientation; precision; location and barycentric corrections | Astropy time; SOFA/ERFA; SPICE; AstroTime; hifitime |
| Coordinates and reference frames | Celestial, terrestrial, observer and instrument frames; velocities and proper motions; transformation metadata | Astropy coordinates; SOFA/ERFA; Starlink AST; SPICE/ANISE |
| WCS | Pixel/world transformations, projection conventions, distortions, invertibility and axis identity | WCSLIB; Astropy WCS; GWCS; Starlink AST |
| Tables and catalogues | Typed columns, units, masks, joins, groups, metadata, large-table access and spatial matching | Astropy tables; TOPCAT/STIL; Arrow-based infrastructure |
| Images and cubes | Units, masks, quality flags, variance, WCS, provenance and lazy/chunked storage | NDData; ndcube; NDF; spectral-cube; observatory data models |
| File and protocol support | FITS and VOTable plus ASDF, HDF5, radio Measurement Sets, event formats and VO services | CFITSIO; ASDF; casacore; PyVO; Astroquery |
| Models and fitting | Parameters, bounds, compound models, evaluation, derivatives and fitting interfaces | Astropy modeling; Sherpa; iminuit; SciML |
| Statistical and signal utilities | Robust estimates, clipping, period searches, convolution, uncertainty representations and numerical conventions | Astropy stats/timeseries/convolution; SciPy; Stingray; VARTOOLS |
| Cosmological backgrounds | Distances, expansion rates, ages, volume elements and model parameter conventions | Astropy cosmology; Cosmology.jl; CCL; Colossus |
| Interoperability | Data-object interfaces, serialization, reader/writer registries, modelling conventions and provenance | Astropy ecosystem interfaces; IVOA standards |

The explicit dependencies are important. Astropy's FITS WCS interface and its broader high-level WCS interface serve different scopes; a new implementation should avoid assuming that every instrument transformation can be encoded in FITS WCS keywords. [Astropy WCS documentation](https://docs.astropy.org/en/latest/wcs/index.html)

Photutils, specutils, ccdproc, Astroquery, PyVO, Reproject, Regions, GWCS and ndcube belong in the same analysis because much of the end-user capability lives outside the core repository. Their source links and descriptions appear individually in the catalogue.

## Lensing deserves several connected engines

The survey identifies at least seven distinct capability families within gravitational lensing. Combining them is scientifically useful, but treating them as interchangeable would erase important assumptions.

| Capability | Representative software | What a common implementation would need |
|---|---|---|
| Strong-lens images and mass models | lenstronomy, PyAutoLens, Lenstool, glafic, GLEE, gravlens/lensmodel | Lens potentials and deflections, multiple planes, image solvers, PSFs, source reconstruction, likelihoods and parameter transformations |
| Differentiable lens modelling | Herculens, caustics, jaxtronomy | Derivatives through ray tracing, image formation and reconstruction; accelerator execution |
| Time-delay cosmography | PyCS3, hierArc, lenstronomy, SNTD | Time-series delays, lens distances, stellar kinematics, line-of-sight structure, microlensing and population assumptions |
| Microlensing | MulensModel, pyLIMA, VBMicrolensing | Binary/multiple lenses, caustics, finite sources, limb darkening, orbital motion and parallax |
| Weak-lensing measurement | GalSim, ngmix, metadetect, TreeCorr, LensTools | PSF correction, shape estimators, detection/selection effects, shear calibration, maps and correlations |
| CMB lensing and spherical maps | lenspyx, HEALPix, pixell, NaMaster | Spin fields, spherical transforms, remapping, masks and mode-coupling corrections |
| Wave-optics gravitational lensing | GLoW | Frequency-dependent amplification, diffraction and interference for coherent signals, including gravitational waves |

These groupings are design proposals grounded in the individual project scopes. See [lenstronomy](https://github.com/lenstronomy/lenstronomy), [hierArc](https://github.com/TDCOSMO/hierArc), [VBMicrolensing](https://github.com/valboz/VBMicrolensing), [metadetect](https://github.com/esheldon/metadetect) and [lenspyx](https://github.com/carronj/lenspyx).

GLoW supplies a useful bridge between lensing and gravitational-wave analysis: its Python interface wraps a C numerical core for wave-optics amplification. [GLoW project description](https://glow-astro.org/glow-code.html)

Lensing also connects directly to optical simulation. GalSim, POPPY, HCIPy, STPSF, dLux, adaptive-optics models and detector calibration supply the image-formation machinery that a lens inference engine needs. Radio/VLBI lens observations additionally require visibility-domain and closure-quantity likelihoods rather than only image-plane fitting.

## Einstein equations and relativistic physics

Four levels should be represented explicitly in a unified suite:

1. **Symbolic geometry:** metrics, connections, tensors, curvature, field equations and perturbative expansions. Reference projects include xAct, Cadabra, SageManifolds, GRTensor and EinsteinPy.
2. **Propagation in a prescribed spacetime:** timelike/null geodesics, ray tracing, polarization and radiative transfer. Reference projects include GYOTO, ipole, grtrans, Gradus, Krang and Skylight.
3. **Perturbative dynamics and waveforms:** black-hole perturbations, quasinormal modes, EMRI waveforms and cosmological perturbations. BHPToolkit, Teukolsky, qnm, FastEMRIWaveforms, CAMB, CLASS, Bolt and SymBoltz occupy different parts of this level.
4. **Evolving spacetime coupled to matter:** numerical relativity, relativistic fluids/MHD, initial-data constraints, horizons and gravitational-wave extraction. Einstein Toolkit, SpECTRE, SpEC, GRChombo/GRTeclyn, BHAC, HARM/iharm3D, KHARMA and AthenaK are relevant references.

This is a conceptual classification, not a ranking. A geodesic integrator does not implement the same problem as nonlinear Einstein-equation evolution. Likewise, linear Einstein–Boltzmann calculations and weak-field relativistic cosmological N-body simulations have different validity regimes. Their common interfaces should expose those regimes rather than conceal them. [Einstein Toolkit capabilities](https://einsteintoolkit.org/capabilities.html), [BHPToolkit](https://bhptoolkit.org/), [SymBoltz](https://github.com/hersle/SymBoltz.jl), [gevolution](https://github.com/gevolution-code/gevolution-1.2)

A concrete lineage finding is that the public GRChombo repository now identifies it as the predecessor to GRTeclyn; the latter uses AMReX. Both are useful in a recovery or reimplementation programme because an older implementation can provide an independent comparison case. [GRChombo](https://github.com/GRTLCollaboration/GRChombo), [GRTeclyn](https://github.com/GRTLCollaboration/GRTeclyn)

## Languages and historical knowledge

**IDL is the strongest match for the proprietary language described in the request.** Its official name is Interactive Data Language. The IDL Astronomy Library, SolarSoft, GBTIDL, SPEDAS, Superfit, SpecPro and SDSS utilities provide concrete examples of astronomy software in that environment. The search did not establish a sufficiently direct primary source for the particular claim about Mark Sullivan's personal use, so that attribution remains unresolved. [IDL product](https://www.nv5geospatialsoftware.com/products/IDL), [IDL Astronomy Library](https://github.com/wlandsman/IDLAstro), [Superfit](https://github.com/dahowell/superfit)

Other language ecosystems must remain visible:

- **Fortran:** stellar evolution, radiation transport, stellar spectra, pulsation, older dynamics, cosmology, pulsar timing and observatory libraries. MESA, GYRE, CAMB, MOOG, Turbospectrum, RADMC-3D, MCFOST and TEMPO are distinct examples.
- **Perl:** scientific arrays through PDL and substantial recipe orchestration through ORAC-DR/PICARD. Some observatory suites also contain Perl utilities. This is often where workflow and calibration-selection knowledge lives.
- **SPP and CL:** IRAF's implementation and command-language ecosystem. Understanding its task interfaces and data conventions matters more than assuming it is ordinary C or Fortran.
- **Java:** TOPCAT/STILTS, Aladin Desktop and Orekit demonstrate important catalogue, visualization and astrodynamics capabilities.
- **R:** ProFound, ProFit, ProSpect and celestial cover source extraction, galaxy structure, spectral modelling and astronomy utilities.
- **Julia:** JuliaAstro, Comrade, Octofitter, Gradus, Skylight, Bolt, SymBoltz and SciML provide both specialist applications and alternative numerical architectures.
- **Wolfram Language and Maple:** xAct, parts of BHPToolkit, FeynCalc/FeynRules and GRTensor preserve symbolic-physics capabilities. Public package source and proprietary runtime availability are different questions.
- **S-Lang, Yorick, Tcl and Scheme:** ISIS X-ray spectroscopy, MiRA/GYOTO-related workflows, legacy interactive systems and electromagnetic solvers show why the search cannot be limited to the largest languages.

Each example is linked in the curated catalogue. Language fields identify the principal implementation/interface family, not audited line-count percentages; generated code, wrappers and native dependencies are stated where material.

Legacy status requires evidence. IRAF lost continuous institutional development but has a community distribution. ESO lists MIDAS among legacy software while its own page advertises a 25MAY release. The NASA IDL-library site was frozen and directs users to GitHub. AIPS has current maintenance notes. Conversely, a web page or package index that still loads does not establish current maintenance. [IRAF](https://iraf-community.github.io/), [ESO legacy list](https://www.eso.org/sci/software/legacy.html), [MIDAS](https://www.eso.org/sci/software/esomidas/), [IDL archive notice](https://asd.gsfc.nasa.gov/archive/idlastro/idlfaq.html), [AIPS](https://www.aips.nrao.edu/)

## Proposed structure of the Rust software unit

The following module boundaries are my synthesis of the surveyed capabilities. They are not upstream project recommendations. One installable product can expose a single workspace, command line, library API and user interface while preserving explicit internal modules and scientific model identities.

| Proposed module | Responsibilities | Primary reference families |
|---|---|---|
| Quantities and frames | Units, constants, epochs, coordinate transformations and ephemerides | Astropy, SOFA/ERFA, AST, SPICE, hifitime, ANISE |
| Scientific data | Arrays, masks, variance/covariance, tables, images, cubes, events, visibilities, particles and meshes | Astropy, NDF, casacore, ndcube, yt |
| Formats and archives | Readers/writers, schema evolution, VO protocols, service adapters and catalogue matching | CFITSIO, ASDF, PyVO, Astroquery, TOPCAT |
| Calibration and pipelines | Detector response, reference selection, task execution, caching and provenance | CRDS, stpipe/stcal, Rubin Butler, ESO CPL/EDPS, ORAC-DR |
| Measurement | Detection, extraction, photometry, spectroscopy, timing, shapes and image reconstruction | Photutils, PypeIt, DAOPHOT, XSPEC, Stingray, radio imaging packages |
| Instrument simulation | Optics, atmosphere, detector noise, response matrices, beams and correlators | GalSim, HCIPy, POPPY, OSKAR, Geant4, MARX |
| Geometry and gravity | Potentials, orbits, lens planes, relativistic metrics and geodesics | galpy, REBOUND, lensing packages, GYOTO, Gradus |
| Physical evolution | Stellar evolution, hydrodynamics, MHD, kinetic plasmas and numerical relativity | MESA, AMUSE, fluid codes, PIC codes, Einstein Toolkit |
| Radiation and microphysics | Opacities, atomic/molecular processes, dust, chemistry, nuclear reactions and transfer | Cloudy, CHIANTI, AtomDB, RADMC-3D, SKIRT, Grackle, pynucastro |
| Cosmological predictions | Backgrounds, perturbations, structure growth, reionization, foregrounds and observables | CAMB, CLASS, CCL, 21cmFAST, CMB packages |
| Numerical methods | Linear algebra, FFTs, sparse systems, integration, PDE discretization, automatic differentiation and optimizers | SciPy, SciML, PETSc, SUNDIALS, FEniCS, Dedalus |
| Statistical inference | Likelihoods, priors, hierarchical models, selection functions, samplers and validation | Sherpa, Cobaya, CosmoSIS, Bilby, Stan, NumPyro, nested samplers |
| Exploration and operations | Visualization, diagnostics, scheduling, instrument control and collaboration | DS9, CARTA, Glue, Jdaviz, Astroplan, TOM, INDI/ASCOM |
| General physics extensions | Quantum systems, particle transport, electronic structure, electromagnetics and molecular dynamics | QuTiP, QuantumOptics, Geant4, Meep, Quantum ESPRESSO, LAMMPS |

The shared object model should preserve the following distinctions:

- **Observed data versus inferred quantities.** Detector counts, calibrated fluxes, reconstructed sources and posterior samples should not become indistinguishable arrays.
- **Variance versus covariance.** Reprojection, coaddition, extraction, calibration and shared model training can introduce correlations that a per-element error bar cannot describe.
- **Model versus solver versus fit.** The same physical equations can have several numerical solvers; the same solver can support several likelihoods and priors.
- **Code versus scientific assets.** Response matrices, ephemerides, opacity tables, atomic rates, passbands, stellar libraries, trained emulators and calibration files need versioned identities.
- **Physical convention versus numerical representation.** Coordinate frames, metric signature, velocity definition, Stokes convention, magnitude system, distance convention and normalization must be machine-readable.
- **Selection versus observation.** Detection thresholds, masks, trigger criteria and follow-up selection belong in the model when the analysis depends on them.

These are proposed interface requirements inferred from the combined scope, not features claimed to be implemented by one existing package.

## How the existing Rust projects fit

The inventory includes ndarray, nalgebra, faer, uom, RustFFT, argmin, diffsol and Burn as candidate numerical building blocks. It also includes domain-specific Rust systems: hifitime, ANISE, Nyx, CDS HEALPix and CDS MOC. Their inclusion is discovery, not a recommendation based on benchmarks or numerical validation.

Rust is already inside some astronomy interfaces: CDS's HEALPix and MOC implementations are used from other languages, and Aladin Lite describes Rust core libraries. The Rust fitsio crate, by contrast, wraps CFITSIO; using it would not make FITS support entirely Rust-native. This distinction should be explicit when deciding whether the final product means a Rust application with external engines or an implementation whose computational core is fully Rust. The inventory supports either end goal. [CDS HEALPix](https://github.com/cds-astro/cds-healpix-rust), [CDS MOC](https://github.com/cds-astro/cds-moc-rust), [Aladin Lite](https://github.com/cds-astro/aladin-lite), [rust-fitsio](https://github.com/simonrw/rust-fitsio), [ANISE](https://github.com/nyx-space/anise)

For a fully Rust implementation, existing engines remain valuable as specifications, test oracles and sources of published algorithms. Agreement between ports of a shared implementation is weaker evidence than agreement between independent algorithms. A future validation programme should record ancestry as well as numerical residuals.

## Recovery records for old software

For each selected legacy component, the follow-on record should identify the algorithm, scientific conventions, source revision or surviving manual, source-access conditions, representative input/output cases, calibration/model data, and successor implementations. This is a proposed research record, not work claimed complete for all entries.

Particularly valuable recovery targets include IRAF task packages, MIDAS applications, Starlink spectral tasks, DAOPHOT/DoPHOT, IDL survey and spectral-analysis utilities, older Fortran stellar/transport codes, and non-Python lens-model tools. Where source is unavailable, a manual or paper can still define capabilities worth reimplementing, but it cannot establish permission to copy an implementation or supply missing test data.

## Evidence and coverage limits

The curated descriptions are concise capability summaries based on official repositories, project/observatory documentation, author pages, and selected software papers or ASCL records. Source retrieval records retain timestamps, redirects and checksums of extracted text. HTTP success is reported separately from substantive page content. Some successful responses contain only metadata, JavaScript shells, directory listings or redirects; some older sites are inaccessible. These remain visible in the evidence report.

The ASCL breadth index captures all records in the registry's browse listing at collection, including submitted records. It does not assert that every submission is accepted, every code remains downloadable, or every entry is currently used. Titles and record links are retained; abstracts and source code are not bulk reproduced in the delivered index. Counts from the registry and curated catalogue must not be added as if disjoint.

This is a broad, inspectable sweep rather than a literal census of every astronomy program ever written. Remaining gaps include private collaboration code, unpublished scripts, complete per-instrument recipe inventories, every Julia/R/IDL extension, proprietary binary internals, and exhaustive general-physics coverage. The language fields are orientation labels, not full source audits. Licenses, buildability, adoption, numerical accuracy and maintenance have not been individually audited for every entry. No software was installed, benchmarked or scientifically validated as part of this survey.

The search was guided by the project context in this workspace, but software facts were checked against external sources during this sweep. The catalogue is the starting inventory; the proposed module structure is the analysis that makes it useful for the Rust-suite objective.
