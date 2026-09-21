#!/usr/bin/env python3
"""Build a version/hash-pinned literature map with explicit reading/reproduction limits."""
from pathlib import Path
import json,hashlib,datetime,re
ROOT=Path(__file__).resolve().parents[3];OUT=ROOT/'phase2/literature'

def source(path,url=None):
 p=ROOT/path
 return {'path':path,'url':url,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size}

catalog=json.loads((ROOT/'catalog/papers.json').read_text());entries={}
for p in catalog:
 id=p['arxiv_id'];pdf=p['pdf_paths'][0] if p['pdf_paths'] else 'papers/pdf/popovic2024-dust-published.pdf';stem=Path(pdf).stem;txt=pdf.replace('/pdf/','/text/').replace('.pdf','.txt')
 entries[id]={'id':id,'title':p['metadata']['title'][0],'source_url':'https://arxiv.org/abs/'+id,'exact_version':stem,'sources':[source(pdf,'https://arxiv.org/pdf/'+stem)]+([source(txt)] if (ROOT/txt).exists() else []),'reading_status':'phase1-focused-summary-carried','methods_read':[],'claim':'See pinned paper; no additional phase2 numerical claim.','data_level':'Published source and linked products inventoried in phase1.','likelihood':'Not re-audited in phase2 for this supporting source.','assumptions':{'population':'Not re-extracted in phase2.','dust':'Not re-extracted in phase2.','age':'Not re-extracted in phase2.','colour_stretch':'Not re-extracted in phase2.','selection':'Not re-extracted in phase2.'},'priors':'Not re-extracted; do not infer an uninformative prior from this field.','sample_reuse':'No independence claim. Crossmatch needed before joint use.','known_critiques':[],'counterevidence':[],'reproduction':{'status':'source-map-only-in-phase2','data_code':['docs/literature-and-data-status.md'],'exact_limits':['This entry preserves coverage of the earlier source review; phase2 method-level re-reading is not claimed.']},'preregistered_tests':[],'roles':['supporting']}
for a in json.loads((OUT/'sources/acquisition.json').read_text()):
 if 'error' in a:continue
 id=a['arxiv_id'];entries[id]={'id':id,'title':a['title'],'source_url':'https://arxiv.org/abs/'+id,'exact_version':a['version'],'sources':a['files'],'reading_status':'focused-methods','methods_read':[],'claim':'','data_level':'','likelihood':'','assumptions':{},'priors':'','sample_reuse':'','known_critiques':[],'counterevidence':[],'reproduction':{'status':'source-and-method-review','data_code':[],'exact_limits':['Author cosmological result not independently reproduced in this phase2 literature task.']},'preregistered_tests':[],'roles':[]}

def setentry(id,claim,data,likelihood,sections,population,dust,age,cs,selection,priors,reuse,critique,counter,code,limits,tests,roles):
 e=entries[id];e.update(claim=claim,data_level=data,likelihood=likelihood,reading_status='focused-methods',methods_read=sections,assumptions={'population':population,'dust':dust,'age':age,'colour_stretch':cs,'selection':selection},priors=priors,sample_reuse=reuse,known_critiques=critique,counterevidence=counter,roles=roles,preregistered_tests=tests)
 e['reproduction']={'status':'source-and-method-review','data_code':code,'exact_limits':limits}

DES='sources/repos/des-science__DES-SN5YR@1.3'; PPLUS='sources/repos/PantheonPlusSH0ES__DataRelease'
setentry('2401.02945','DES uses 1,635 photometric DES events plus 194 low-z SNe after all analysis cuts; 1,499 DES objects have nominal P(Ia)>0.5, but the likelihood retains the broader mixture.',
 'SMP epoch fluxes, SALT3 model, released fitted parameters/distances, classification probabilities, 25 selected Ia/non-Ia mocks, covariance and PIPPIN inputs.',
 'SALT3 flux fitting; BEAMS Ia/contaminant mixture with BBC4D bias corrections; binned nuisance estimation followed by unbinned distance Gaussian cosmological likelihood, with BEAMS-adjusted uncertainties.',
 ['Sections 3.1-3.6; Eqs 1-12','Sections 4.1-4.3; 5.1 and Table 4; 6.4.3; 6.5'],
 'SNANA survey simulations; intrinsic-colour Gaussian plus dust reddening; stretch distribution depends on host mass; volumetric rate sets generated redshifts.',
 'BS21/P21 forward model with separate intrinsic colour, intrinsic beta, RV and exponential E(B-V); host-dependent dust distributions; effective beta in standardization is distinct.',
 'No explicit progenitor age latent variable in nominal analysis; residual host mass/colour step may absorb other astrophysics.',
 'alpha,beta,gamma fitted by BBC; intrinsic-colour and stretch population conditioning in simulation; raw covariance required for new joint inference.',
 'Pipeline detection + host-spec-z + SALT cuts + outlier removal + valid bias correction + intersection of systematic variants. Selected mocks are not all generated events.',
 'Survey/population calibration and systematic scales supply assumptions; nuisance fitting precedes cosmology; exact cosmological priors must follow release config per model.',
 'Shares low-z anchors and calibration/training with Pantheon+; DES objects overlapping spectroscopic samples have different selection conditioning.',
 ['2408.07175','2411.05299','2510.13121'],['2501.06664','2601.13785'],[DES,'phase2/official/inputs/SNDATA_ROOT'],
 ['No completed end-to-end raw hierarchical cosmology in this literature task. Exact selected-mock p0 file equality and final fit selection must be verified.','Official distance covariance cannot simply be attached to a different raw-SALT hierarchy.'],
 ['SEL-01','SEL-02','SEL-03','SEL-04','POP-01','CONT-01','CAL-01','OVL-01'],['DES-original','methods','selection'])
setentry('2401.02929','Original DES5YR cosmology uses the calibrated and bias-corrected release; SN-only and external-probe combinations answer different questions.',
 '1829-event distance vector and covariance; raw-to-distance methods delegated to Vincenzi et al.',
 'Gaussian distance likelihood in specified FLRW models; external CMB/BAO combinations explicitly separate.',
 ['Cosmological model and results sections; source companion 2401.02945 Sections 3-6'],
 'Inherits DES-SN5YR simulation population and nuisance treatment.','Inherits P21/BS21 dust corrections.','No fitted free progenitor-evolution law.','SALT3 standardization; nuisance fit inherited.','All Table 4 selection stages in companion paper apply.',
 'Model-specific flatness/dark-energy parameterization and external likelihood priors; use release configs, not abstract intervals.',
 'Original and Dovekie are overlapping reductions, not independent surveys.','2408.07175 2510.13121'.split(),['2501.06664','2601.13785'],[DES+'/5_COSMOLOGY'],
 ['This source map does not certify a new raw-photometry reproduction.'],['CAL-01','COS-01','OVL-01'],['DES-original','cosmology'])
setentry('2406.05046','DES public light curves support raw-flux reanalysis but the cosmology sample is already selected.',
 'Scene-model photometry and difference-imaging products with metadata; public release companion.',
 'Photometry extraction and validation rather than new cosmological likelihood.',
 ['Data-release/photometry product descriptions; release READMEs checked against actual HEAD/PHOT tables'],
 'Sample ascertainment and astrophysical types must be retained in modeling.','Calibration/Milky-Way treatment lives in released metadata and models.','Host age is not supplied as a direct progenitor-age observation.','Flux-to-SALT estimates require an independently pinned model and passbands.','Released sample boundaries differ between full transients and cosmology SMP.',
 'No standalone cosmological prior.','Reuses the DES5YR observations underlying 2401.02929.',[],['2401.02945'],[DES+'/0_DATA'],
 ['Epoch-level covariance/model training/calibration choices must be audited; raw flux availability does not by itself release an all-generated selection denominator.'],['SEL-01','CAL-01'],['DES-original','data-release'])
setentry('1811.02379','Image-fake-calibrated simulations model DES observing conditions and selection-related distance bias.',
 'SNANA catalogue simulations based on 10,000 overlaid fake light curves, cadence/noise/PSF/zero-points and detection efficiency.',
 'Forward simulation and bias-correction validation; cosmology applied after SALT fitting and sample cuts.',
 ['Simulation construction and detection-efficiency sections; release efficiency documentation'],
 'Generated SN Ia distributions and rates supplied to simulator.','Alternative intrinsic-scatter models affect distance biases.','No direct progenitor-age inference.','Multidimensional bias dependence includes colour/stretch.','Image recovery calibrates per-epoch detection; later quality and follow-up selection still required.',
 'Input simulation cosmology/scatter/populations are assumptions to vary.','DES3YR detection calibration reused in DES5YR; fake overlays are not a new cosmological SN sample.',
 ['Selected-distribution agreement is insufficient to identify a misspecified population/selection combination.'],['DES5YR uses host-based redshift selection rather than identical DES3YR spectroscopy.'],
 ['phase2/official/inputs/SNDATA_ROOT/models/searcheff/SEARCHEFF_PIPELINE_DES.DAT'],
 ['Per-epoch detection table lacks injection-level counts/uncertainties in this release; joint fit efficiency needs simulation.'],['SEL-01','SEL-03'],['DES-simulation','selection'])
setentry('2012.07180','Host-redshift efficiency is measured relative to galaxies satisfying OzDES targeting criteria, conditional on host brightness, colour, field and discovery epoch.',
 'Targeted host-galaxy sample with/without successful redshifts; DES difference imaging and CC/Ia simulations.',
 'Empirical binwise host-z efficiency plus forward SN/host simulations; validation of observed distributions and contaminants.',
 ['Sections 2.3-2.4; 3.1-3.2.3; Figs 1-2'],
 'SN rates depend on galaxy stellar mass/SFR; Ia and CC have different host populations.','Type-dependent extinction and luminosity/template assumptions enter forward simulation.','No direct age likelihood; host demographics are not progenitor delays.','SALT2 loose/tight cuts; host colour is observed g-r, not SN SALT c.','Success denominator is targeted-eligible hosts, not every galaxy or generated SN. Redshift-from-SN spectra excluded from this efficiency branch.',
 'Rates/templates imported from prior studies; host-z kernel is empirically estimated.','Same DES observing sample; distinct reduction or host subsample does not create independence.',
 ['Target eligibility/host association and quality-selection conditioning must not be omitted.'],['2402.18690 demonstrates the host-z requirement removes a systematically different SN population.'],
 ['phase2/official/inputs/SNDATA_ROOT/models/searcheff/SEARCHEFF_zHOST_DES-SN5YR.DAT'],
 ['Published efficiency grid does not expose per-cell success/target counts; finite calibration uncertainty cannot be estimated from grid alone.'],['SEL-01','SEL-02','CONT-01'],['DES-simulation','selection','contamination'])
setentry('2307.13696','Host mismatches can generate wrong redshifts and must be propagated through the selection and classification pipeline.',
 'DES5YR host catalogue and SNANA host-matching simulations.','Directional-light-radius association, simulated mismatches, and cosmology comparisons.',
 ['Host-matching and simulation methodology; data availability'],
 'Host library and SN-host rate weighting determine match population.','Inherits simulated SN dust/scatter.','Host age association inherits host-identification uncertainty.','Wrong host changes z and host covariates jointly.','Match quality and redshift availability correlate with selection.',
 'Host catalogue and matching model supply prior structure.','Same DES host and observing sample.',[],['2401.02945 includes mismatch-related systematic checks.'],[DES+'/1_SIMULATIONS'],
 ['A marginal mismatch fraction is not an independent redshift-error PDF for each event.'],['SEL-02','CONT-01'],['DES-hosts','redshift'])
setentry('2111.10382','Small simulated contamination biases are conditional on the tested non-Ia template families and cosmological priors.',
 'DES-like Ia/non-Ia light curves, SuperNNova and other classifier variants.',
 'BEAMS/BBC mixture, redshift-binned distances and w/CPL fits on simulations.',
 ['Classifier/BEAMS methodology and simulation-validation design'],
 'Ia plus several non-Ia template families and population rates.','Inherited light-curve/scatter assumptions.','No explicit age drift.','SALT quality cuts reshape contamination.','Classifier accuracy is evaluated after survey and quality selection.',
 'Flat universe; published w validation uses Gaussian Ωm=0.311±0.010; this is not SN-only evidence.',
 'Training/test simulations and DES5YR analysis machinery overlap.',
 ['Classifier probabilities can be miscalibrated outside training population support.'],['Tests multiple non-Ia models and classifiers rather than assuming perfect purity.'],[DES+'/3_CLASSIFICATION',DES+'/7_PIPPIN_FILES'],
 ['Reported bounded w biases do not certify arbitrary population/selection alternatives.'],['CONT-01','POP-01'],['DES-classification','contamination'])
for id,claim,sections in [
 ('2201.11142','SuperNNova constructs a DES photometric Ia sample using DES-like simulations and redshift-aware classification.',['Preprocessing, loose cuts, training/validation sections']),
 ('2402.18690','Classification without host redshift recovers a broader detected population; missing host-z events differ in colour/stretch/host brightness.',['Photometric-only classification; photometric redshift and sample-comparison sections'])]:
 setentry(id,claim,'DES difference-imaging light curves and classifier probabilities.','Probabilistic photometric classifier; this catalogue is not itself a selection-corrected cosmological likelihood.',sections,
 'Training simulations define class populations.','Dust affects colours used by classifier.','No directly measured progenitor-age variable.','Classification and photo-z use the same flux histories as SALT fitting.','Detection, transient quality and redshift-availability masks must remain explicit.',
 'Training class balance/calibration rather than a cosmological prior.','Same DES observations as original cosmology, with broader/different masks.',
 ['P(Ia|photometry,z) cannot be treated as independent additional data beside the same light-curve likelihood.'],['Redshift-free sample provides a withheld population diagnostic.'],[DES+'/3_CLASSIFICATION'],
 ['No independent survey; no all-generated injection denominator supplied by class-probability catalogue.'],['CONT-01','SEL-02'],['DES-classification'])
for id,claim,sections,role in [
 ('2104.07795','SALT3 learns an empirical spectral time-series model; training data and calibration enter subsequent light-curve inference.',['Training and model-validation sections'],'light-curve-model'),
 ('2301.10644','SALT2/SALT3 are compared after identical training; this isolates framework differences conditional on training/calibration.',['Training, DES3YR application, simulation comparisons'],'light-curve-model'),
 ('2112.03864','Fragilistic cross-calibration propagates correlated passband/zero-point uncertainty and model retraining.',['Calibration solution and uncertainty propagation'],'calibration')]:
 setentry(id,claim,'Training light curves/spectra, passbands and calibration products.','Empirical SED training/calibration fit; subsequent cosmology depends on fitted-parameter covariance.',sections,
 'Training sample coverage controls interpolation and extrapolation.','Empirical colour law does not uniquely separate dust and intrinsic colour.','No identification of physical progenitor age from SALT training.','SALT x1/c are empirical coordinates; coordinate changes across trained models must be checked.','Training sample and validation selection differ from cosmology selection.',
 'Training regularization/calibration constraints are part of model, not absent assumptions.','Training/calibration is shared by multiple SN compilations.',
 ['Good standardization residuals do not establish population invariance or physical dust identification.'],['Independent training/calibration variants provide sensitivity checks.'],[DES+'/2_LCFIT_MODEL',PPLUS+'/Pantheon+_Data/2_CALIBRATION'],
 ['Full retraining and covariance reconstruction not executed by this literature task.'],['CAL-01','OVL-01'],[role])
setentry('1710.00845','The PS1 spectroscopic sample joins other surveys in Pantheon and requires survey-specific selection corrections.',
 '365 PS1 spectroscopically confirmed light curves; 279 useful-distance PS1 objects in the 1048-SN Pantheon compilation.',
 'SALT2 standardization, simulation bias correction and covariance-distance cosmology.',
 ['Survey/selection and calibration methods; sample construction'],
 'Survey-specific parent distributions and intrinsic scatter.','Colour/brightness scatter models tested.','No direct free progenitor-age evolution fit.','x1/c standardization and population/selection corrections.','Spectroscopic follow-up selects a different PS1 population from photometric PS1.',
 'External CMB/BAO combinations distinguished from SN-only; exact priors depend on cosmological model.',
 'PS1 objects/calibration recur in Pantheon+ and training; not an independent check if reused.',
 ['Heterogeneous low-z selection/calibration dominate some uncertainty.'],['1710.00846 uses broader photometric sample and BEAMS sensitivity tests.'],['https://github.com/dscolnic/Pantheon'],
 ['No new PS1 raw photometry fit or object crossmatch here.'],['OVL-01','SEL-01'],['Pan-STARRS','Pantheon'])
setentry('1710.00846','Photometric PS1 cosmology marginalizes a contaminant distribution and tests classification priors; apparent robustness is conditional on those variants.',
 '1169 PS1 events and 195 low-z SNe.','Three-Gaussian BEAMS: Ia high/low host mass plus CC, with redshift control points and distance covariance passed to cosmology.',
 ['Sections 2.3, 3.1; Table 4'],
 'SN Ia and CC distributions with redshift dependence.','G10/C11-type scatter and colour-standardization variants.','No explicit age latent.','alpha,beta,host step fitted; stretch/colour selection included via simulation.','Host-z and live-SN spectroscopic selection separately modeled.',
 'Table 4: alpha .155±.05, beta 2.947±.50, mass step .07±.07; classifier renormalization A=1±.2 and shift S=0±.2; broad contaminant priors.',
 'Reuses PS1 survey/calibration and nearby anchors; not independent of Pantheon PS1 subset.',
 ['A flexible contaminant mean can absorb cosmological effects; class priors from same flux require conditional treatment.'],['Four classification priors and alternative CC models tested.'],['https://github.com/djones1040/BEAMS'],
 ['Published model not executed here; PS1 parent-selection assets not yet acquired into phase2.'],['CONT-01','OVL-01'],['Pan-STARRS','contamination'])
setentry('2112.03863','Pantheon+ publishes light curves and cross-survey metadata; duplicate observations and common calibration affect independence.',
 'Raw light curves and metadata for 1701 light curves/1550 distinct SNe in published Pantheon+ framework.','Data release and SALT-based quality cuts; cosmology likelihood in companion Brout paper.',
 ['Dataset, selection and release-product sections'],
 'Heterogeneous discovery/follow-up populations.','Colour law and scatter inferred in companion analyses.','Host mass is not progenitor age.','SALT2 summary measurements shared with cosmology products.','Survey-specific quality/redshift cuts and duplicated light curves.',
 'Not a standalone cosmological posterior.','Repeated light curves and shared objects must be counted at event rather than table-row grain.',[],['2202.04077'],[PPLUS],
 ['No raw hierarchical Pantheon+ refit in phase2 literature work.'],['OVL-01'],['Pantheon+','data-release'])
setentry('2202.04077','Pantheon+ already encodes host/dust/scatter effects in bias corrections and error estimates; a small residual mass step does not mean absent host correction.',
 'SALT2-B22 distances, bias-correction products, full covariance and optional Cepheid-calibrator data.',
 'BBC-derived standardization and covariance Gaussian distance likelihood; optional SH0ES calibrator likelihood separately.',
 ['Sections 2.1-2.3, 3.1-3.3; Eqs 1-3; Table 2'],
 'Heterogeneous survey populations modeled with SNANA and dust/scatter prescriptions.','BS21 nominal plus P21/G10/C11 variants; Milky-Way correction and model retraining.','No direct free age law; host-correlated scatter corrected by simulation.','Effective beta and alpha; full SALT covariance plus simulation-derived uncertainty scaling.','Survey detection and spectroscopic follow-up, quality cuts, valid bias corrections.',
 'Cosmological model priors and external CMB/BAO separate; Cepheid calibration must not silently enter an uncalibrated SN-only fit.',
 'DES, PS1, low-z and Union training/calibrations overlap.','2411.05299 2510.13121 2601.19424'.split(),['2601.13785','2501.06664'],[PPLUS],
 ['Modern corrected residuals cannot be interpreted as raw luminosities. Full raw-fit covariance recalibration not reproduced here.'],['POP-01','CAL-01','OVL-01'],['Pantheon+','methods'])
setentry('2112.04456','Dust2Dust forward-models colour/scatter populations rather than estimating a redshift-only brightness correction.',
 'SNANA forward simulations and observed colour/HR distributions.','Simulation-driven fit of intrinsic colour, beta, dust RV and reddening distributions.',
 ['Model/distribution definitions and inference design (supported by DES 4.2 and released GENPDF)'],
 'Separate intrinsic-colour Gaussian and exponential reddening; host-dependent RV/tau.','Dust and intrinsic-colour latent variables; inference depends on their chosen functional forms.','No direct progenitor-delay likelihood.','Colour-luminosity relation combines beta_int and dust extinction.','Selection is applied to generated simulated population before observed comparisons.',
 'Parameterized distribution families and fitted population hyperparameters; tabulated simulation PDF is not an assumption-free prior.',
 'Inputs recur in Pantheon+/DES corrections; inference is not independent of those samples.',
 ['Dust/age/intrinsic-colour alternatives can be observationally degenerate.'],['Reproduces multiple distribution summaries rather than only mean residual.'],['phase2/official/inputs/SNDATA_ROOT/models/population_pdf/DES-SN5YR'],
 ['Exact original DES mock GENPDF filename differs from distributed table; equality not established.'],['POP-01','SEL-04'],['dust','population'])
setentry('1610.04677','BBC is a practical approximate two-stage method; its own motivation distinguishes it from a full parent-population likelihood.',
 'Simulated Ia+CC light curves and SALT fits.','BEAMS mixture with simulation bias corrections and binned distances, then cosmology.',
 ['Introduction; Section 5 likelihood/bias correction; Section 8 validation discussion'],
 'Simulated parent populations enter corrections.','Intrinsic scatter prescriptions assumed.','No direct age inference.','Simultaneous nuisance fitting; corrected mean distances.','Monte Carlo calibration includes selection and fit bias.',
 'Flat nuisance ranges in first stage; cosmological priors in second stage, including fixed Ωm in some validation tests.',
 'Foundational method reused in DES/Pantheon, not independent data.',
 ['Approximate likelihood depends on validation domain and simulation adequacy.'],['Simulation validation measures bias under multiple tested scenarios.'],['https://github.com/RickKessler/SNANA'],
 ['No proof for arbitrary out-of-family populations from nominal-mock success.'],['POP-01','SEL-03'],['methods','BBC'])
for id,claim,sections in [
 ('1507.01602','UNITY jointly models latent light-curve populations, nonlinear standardization, outliers and selection.',['Model and selection-likelihood sections']),
 ('2311.12098','UNITY1.5 extends joint population/selection/systematic inference to Union3; survey-depth priors matter especially for sparse nearby samples.',['Sections 2.2 and 4.5; prior/selection discussion; Appendix B is identified but not fully transcribed'])]:
 setentry(id,claim,'SALT summary vectors/covariance and heterogeneous SN sample metadata.','Bayesian hierarchical likelihood with latent populations, normalization for detection, outlier mixture and systematic nuisance parameters.',sections,
 'Population parameters modeled rather than fixed to one universal observed Gaussian.','Colour population/standardization and residual covariance; not unique physical dust identification.','Host dependence considered; a direct progenitor-delay law is a separate extension.','Nonlinear relations and latent measurement error.','Efficiency model enters normalization; survey thresholds and widths are uncertain.',
 'For Union3, survey-depth priors approximately ±0.5 mag, tested at ±0.25 mag; prior table must be checked for complete reproduction.',
 'Large overlap with earlier compilations; changed likelihood is not independent new photons.',
 ['Prior sensitivity and simplified selection function can matter; verify with injections.'],['Framework permits population drift and forward selection normalization.'],['https://github.com/rubind/Unity'],
 ['No independent UNITY execution; complete prior transcription not claimed.'],['POP-01','SEL-03','PRIOR-01'],['methods','hierarchical'])
setentry('1506.01354','Marginal acceleration significance is obtained under global independent Gaussian latent M,x1,c populations and the specified cosmological model.',
 'JLA SALT2 summaries with statistical/systematic covariance.','Analytic Gaussian latent integral including determinant; profile maximum likelihood over cosmology and population parameters (Eqs 4-10).',
 ['Methods Eqs 4-10 and profile-likelihood section'],
 'One global independent Gaussian each for M,x1,c; population means/variances fitted.','No explicit nonnegative dust component.','No explicit redshift-dependent progenitor age.','x1/c population independent of redshift and survey.','Uses JLA corrected input but no full explicit generated selection likelihood.',
 'Profile likelihood rather than declared Bayesian posterior; latent distribution assumptions still constrain inference; curvature/flatness choices change acceleration test.',
 'Same JLA observations as mainstream JLA/Rubin-Hayden comparison.',
 ['1610.08972','1912.02191'],['1912.04257'],['https://arxiv.org/abs/1506.01354'],
 ['No new JLA execution; global observed-population assumption is not a valid default for selected DES5YR.'],['POP-01','COS-01'],['acceleration-critical','population'])
setentry('1610.08972','Allowing survey/redshift-dependent x1/c distributions increases acceleration significance relative to a fixed-population hierarchy.',
 'Same JLA SALT2 summaries.','Hierarchical latent Gaussian model sampled with HMC; cosmological/kinematic variants and q0 inference.',
 ['Statistical model; redshift-independent/dependent distributions; other cosmological models'],
 'Mean x1 and c vary linearly with redshift by discovery sample (HST mean constant); more flexible variant checked.','No unique physical dust model.','Host-age correlations motivate population flexibility but are not measured progenitor ages.','Selected distributions vary with z; fixed global shrinkage can remove real standardization.','Observed population drift is flexible; not a full physical generated selection model.',
 'Flat cosmological/kinematic priors in reported variants; curvature assumptions materially alter significance.',
 'Reanalysis of JLA, not independent events.',
 ['1912.04257 disputes aspects of frame/correction-based evidence in related analyses.'],['Demonstrates model-induced latent shrinkage using observed x1/c trends.'],['https://arxiv.org/abs/1610.08972'],
 ['Population flexibility alone does not identify cosmology if arbitrary luminosity drift is allowed.'],['POP-01','COS-01','PRIOR-01'],['acceleration-counterargument','population'])
for id,claim,crit,counter in [
 ('1912.02191','Velocity-frame and population treatment can spuriously generate a deceleration dipole.',['1912.04257'],['1610.08972']),
 ('1912.04257','The response argues peculiar-velocity corrections and population choices change isotropic-acceleration significance.',['1912.02191'],['1506.01354'])]:
 setentry(id,claim,'JLA redshifts, sky positions and light-curve parameters.','Kinematic monopole/dipole reanalysis with alternate frames/population corrections.',['Frame/velocity treatment and model-comparison sections'],
 'Population assumptions compared with related hierarchical analyses.','No independent physical dust identification.','No new direct progenitor-age likelihood.','Global versus sample/redshift-dependent summary populations.','Survey anisotropy/selection must be propagated, not treated as isotropic sampling.',
 'Frame, velocity model, curvature/kinematic truncation and nuisance restrictions are conditional assumptions.',
 'Same JLA events; correlated argumentative exchange.',crit,counter,['https://arxiv.org/abs/'+id],
 ['No independent phase2 directional refit; see phase1 directional audit for newer Pantheon+ debate.'],['COS-01','POP-01'],['acceleration-critical','frames'])
setentry('2408.07175','A small differential low/high-z offset can materially change evolving-dark-energy inference; source v3 comments on the DES response.',
 'Common DES5YR/Pantheon+ distance objects and external cosmological datasets.','Covariance-distance cosmological fits and inter-compilation magnitude-offset comparison.',
 ['Common-object comparison; bias-correction discussion; v3 response and conclusion'],
 'Questions adequacy of host/scatter parameterization.','Dust/scatter assumptions feed simulation bias corrections.','Local age/metallicity/SFR are possible omitted variables, not measured corrections here.','Distances already standardized under different models.','Different spectroscopic/photometric selection changes expected bias corrections for identical events.',
 'Cosmological model and CMB/BAO assumptions; arbitrary offset is sensitivity intervention, not a measured new calibration likelihood.',
 'Shared events make difference measurements sensitive but not independent cosmology samples.',
 ['2501.06664'],['A residual sensitivity concern remains even after accounting for selection differences.'],[DES,PPLUS],
 ['Offset subtraction does not by itself identify a causal calibration defect.'],['CAL-01','OVL-01','SEL-01'],['DES-critical','calibration'])
setentry('2501.06664','DES explains much of the matched-object offset through scatter/host revisions and differing selection corrections.',
 'Common DES5YR/Pantheon+ events; source methods use 118 Foundation and 145 DES common objects.','Factorial distance/bias-correction comparisons within the DES framework.',
 ['Sections 2-5; Table 1; methodological appendix'],
 'Reverts DES scatter/population and host estimates to Pantheon+ alternatives.','Dust/scatter changes contribute to bias corrections.','No direct test of arbitrary progenitor luminosity evolution.','Different trained models and nuisance corrections considered.','Same SN has different expected bias correction when conditioning on different selection mechanisms.',
 'Comparison conditional on BBC simulations; no universal proof against all systematics.',
 'Explicit common-event comparison; neither catalogue supplies independent photons for those objects.',
 ['2408.07175'],['Measured factorial changes explain why per-object corrected distances need not coincide.'],[DES,PPLUS],
 ['Full historical intermediate pipeline products not all reproduced in this literature task.'],['CAL-01','SEL-01','OVL-01'],['DES-counterargument','selection'])
setentry('2503.14738','DESI DR2 BAO measures relative distances to the sound horizon; CMB/BBN calibration and SN compilation choices must be reported separately.',
 'BAO DM/rd, DH/rd or DV/rd measurements/covariance plus public cosmological likelihood configs/chains.','BAO distance likelihood; separate full CMB or early-universe compressed-prior combinations.',
 ['Section IV A-B and supernova discussion; Appendix A purpose checked'],
 'Galaxy/quasar tracer selection and reconstruction differ from SN populations.','SN dust inherited only when combining SN likelihood.','No direct SN progenitor-age measurement.','Not applicable to BAO observable.','Tracer/BAO reconstruction selection; SN sample normalization remains a separate likelihood.',
 'BBN Ωb h²=.02218±.00055 for standard neutrino case; θ*=1.04110±.00053 in 100θ units for acoustic-only variant; full CMB differs and requires perturbations.',
 'DR1/DR2 overlap; do not multiply DR1 and DR2 as independent. Pantheon+/Union3/DESY5 compared separately because overlap.',
 ['SN likelihood choice changes dynamical-dark-energy preference.'],['BAO provides geometry independent of supernova luminosity, conditional on rd treatment.'],['data/bao/desi-dr2-reference','sources/repos/CobayaSampler__bao_data'],
 ['Compressed CMB or Ωm priors cannot be labeled a full CMB reproduction.'],['COS-01','PRIOR-01','OVL-01'],['DESI','external-probes'])
setentry('1807.06209','Planck cosmological constraints are inferred through a physical early-universe/perturbation model and foreground likelihood.',
 'CMB temperature/polarization spectra and lensing likelihoods.','Hybrid low-ell and high-ell likelihoods plus lensing; Boltzmann prediction and nuisance marginalization.',
 ['Sections 2.1-2.3; baseline theoretical model and power-spectrum likelihood'],
 'Primordial perturbation and matter model rather than SN population.','Galactic/extragalactic foreground dust differs from SN host extinction.','No SN age information.','Not applicable.','Sky masks, frequency selection and foreground cleaning rather than SN efficiency.',
 'Baseline flat six-parameter ΛCDM, adiabatic nearly power-law scalar spectrum; standard neutrino setup with one .06eV massive species; extensions alter conclusions.',
 'Planck spectra reused in many combined analyses; CMB posterior priors are not independent repeated observations.',
 ['Late-time model changes can invalidate simplistic fixed-Ωm compression.'],['Full likelihood or model-valid early-universe compression makes assumptions auditable.'],['https://pla.esac.esa.int/','sources/repos/CobayaSampler__cobaya'],
 ['This task does not execute a Planck likelihood; reference-only external-prior audit.'],['COS-01','PRIOR-01'],['CMB','external-probes'])

# Central age and updated-standardization exchange: carry exact phase1 method audits,
# explicitly distinguishing them from a new full-text reanalysis.
carry={
 '2411.05299':('Age-luminosity correlation claim from revised host ages.','Published paper plus correction; G11/R19 overlap; corrected and original residual estimands differ.','docs/investigations/age-signal-audit.md',['2601.13785']),
 '2510.13121':('Age-evolution correction changes late-time cosmological inference.','Published methods/priors and CSFH/DTD read in phase1; external CMB and self-consistent correction remain limited.','docs/literature-and-data-status.md',['2601.13785']),
 '2601.13785':('Modern host/dust corrections weaken an additional age residual and reduce inferred progenitor evolution.','Published Appendix A; matched sample and W22 simulation examined; cosmic short-DTD number not fully recovered.','docs/literature-and-data-status.md',['2605.21586']),
 '2605.12596':('Age correction is proposed as origin of environmental steps.','Published HTML and distinct 175-object cut/70-young-host correction verified in phase1; mask/posterior inputs missing.','docs/literature-and-data-status.md',['2601.13785']),
 '2605.21586':('Reply argues broad redshift mixing dilutes age relation and transport preserves correction.','Published reply contains mock absent from arXiv v1; exact author mock/85-percent transport still unreproduced.','docs/literature-and-data-status.md',['2601.13785']),
 '2604.16597':('TITAN host SFH modeling changes progenitor-delay estimates.','Host SFH/deconvolution methodology differs from measured HR-age cosmology; public final posterior/selection products absent at phase1 check.','docs/literature-and-data-status.md',['2605.21586']),
 '2606.09650':('Directional Pantheon+ age-corrected analysis reports deceleration.','Full latent C2 and corrected/reversed-bias C1 differ; phase1 directional methods and current responses audited.','docs/investigations/directional-audit.md',['1912.02191']),
 '2511.07517':('Recalibration/dust revisions shift DES cosmology.','Dovekie v3 uses changed calibration and Milky-Way treatment; matched factorial comparison performed in phase1.','docs/investigations/standardization-audit.md',['2408.07175']),
 '2506.05471':('Dovekie reassesses cross-survey calibration uncertainty.','Calibration and standardization source audit performed in phase1.','docs/literature-and-data-status.md',['2408.07175']),
 '2004.10206':('Host-dependent dust can explain correlated luminosity-scatter and mass-step behavior.','BS21 model is a physical latent explanation tested on distributions, not unique causal identification.','docs/literature-and-data-status.md',['2411.05299']),
 '2207.05583':('Galaxy-driven forward model connects host demographics, progenitor delays and SN properties.','Host population/SFH/DTD and selection define transport; selected simulator unavailable in complete original form.','docs/literature-and-data-status.md',['2605.21586']),
 '2002.12382':('Observed host-age correlations need not overturn supernova acceleration.','Counterargument concerns latent errors, covariance and correction transport, not denial of all environment dependence.','docs/literature-and-data-status.md',['2411.05299'])}
published={'2411.05299':['papers/pdf/chung2025-published.pdf','papers/pdf/chung2025-correction-published.pdf'],'2510.13121':['papers/pdf/son2025-published.pdf'],'2601.13785':['papers/pdf/wiseman2026-published.pdf'],'2605.21586':['papers/pdf/chung2026-published.pdf']}
for id,(claim,limit,doc,counter) in carry.items():
 e=entries[id];e['claim']=claim;e['reading_status']='phase1-method-audit-carried';e['methods_read']=[{'scope':'phase1 audit carried; not new phase2 full re-reading','document':'docs/literature-and-data-status.md','details':limit}]
 e['data_level']='Host/corrected-residual or calibration/reanalysis level; see phase1 linked exact source/input audit.'
 e['likelihood']='Published source-specific likelihood already audited in phase1; this entry does not replace its exact methods.'
 e['assumptions']={'population':'Host/SN population transport must be propagated jointly with selection.','dust':'Modern residuals may already include dust/host correction; avoid adding a second correction without redefining estimand.','age':'Host mean stellar age and SN progenitor delay are distinct; conditional host-SFH/DTD and selection required.','colour_stretch':'Existing standardization changes residual-age slope; new raw hierarchy must fit jointly.','selection':'Residual sample/host-age ascertainment and low/high-z composition matter.'}
 e['priors']='Source-specific priors are in the prior audit; no universal age-slope prior endorsed here.'
 e['sample_reuse']='G11/R19 overlap and modern compilation reuse recorded in phase1; no independent-sample label inferred.'
 e['counterevidence']=counter;e['known_critiques']=[limit]
 e['reproduction']={'status':'phase1-partial-reproductions-carried','data_code':['docs/literature-and-data-status.md'],'exact_limits':[limit,'No new raw-data cosmology result from this literature-map task.']}
 e['preregistered_tests']=['POP-01','PRIOR-01','OVL-01'];e['roles']=['age-exchange' if id[0:2] in ['24','25','26'] else 'dust-age-methods']
 for path in published.get(id,[]):
  if (ROOT/path).exists():e['sources'].append(source(path));e['exact_version']+=' + published journal version/correction pinned by hash'

# Structured hypothesis tests, fixed before any phase2 cosmological comparison.
tests=[
 {'id':'SEL-01','hypothesis':'The simulation sample is an adequate generated denominator.','design':'Count README generated/written and DUMP/HEAD rows; inspect applied masks and all selection stages.','pass_condition':'Generated failures available or exact p0, selection-identical importance-ratio proof with support validated.','status':'Selected-only diagnosis complete; full denominator likelihood not yet validated.'},
 {'id':'SEL-02','hypothesis':'Host-z selection is ignorable after conditioning on available host properties.','design':'Evaluate empirical r, g-r, field, year grid; test missing-host/host-z and low/high host-colour residual distributions on broader detected catalogue.','pass_condition':'Predeclared held-out distribution agreement under propagated grid/targeting uncertainty.','status':'Grid reader/audit complete; held-out adequacy test pending.'},
 {'id':'SEL-03','hypothesis':'Candidate inference recovers injected cosmology under identical quality cuts.','design':'Generate broad-support parent injections including failures; run same flux fitter, convergence, x1/c uncertainty, host and probability cuts; split train/validation by mock.','pass_condition':'Coverage and bias tested on held-out non-ΛCDM/population variants, not only nominal model.','status':'Preregistered; not executed by literature task.'},
 {'id':'SEL-04','hypothesis':'Importance ratios over released selected mocks identify candidate normalization.','design':'Check exact generated density, true host conditioning, effective vs intrinsic SALT mapping, absolute continuity, stratified ESS and mock-to-mock variance.','pass_condition':'Reject non-overlap; require ESS>=1000 overall and >=100 per reported conditional cell, max normalized weight<.01; tighten until result stable to Monte Carlo error.','status':'Luminosity-support singularity measured; colour-only stress diagnostics computed; no cosmology certification.'},
 {'id':'POP-01','hypothesis':'One redshift-independent Gaussian x1/c population is adequate.','design':'Compare joint mB,x1,c distributions and host/z conditioning for Gaussian, dust-tail and redshift-drift models using held-out posterior predictive checks.','pass_condition':'Report failed features and cosmology shifts; never select model only by acceleration preference.','status':'Preregistered.'},
 {'id':'CONT-01','hypothesis':'Nominal photometric classifier probabilities remain calibrated under candidate populations.','design':'Use Ia+CC mixtures conditional on the same photometry/host/z; alternative CC templates/classifiers and classification-odds calibration; retain same-flux dependence.','pass_condition':'Held-out purity/calibration and inferred parameter recovery; no double counting of BEAMS posterior probability.','status':'Preregistered.'},
 {'id':'CAL-01','hypothesis':'Calibration/model changes leave inference stable within propagated uncertainty.','design':'Reproduce official release first; vary calibration/model retraining and Milky-Way law with correlated nuisance/covariance; compare common-event raw fits and corrected distances separately.','pass_condition':'All released systematic variants or explicit limits; no diagonal-only reinterpretation of correlated calibration.','status':'Official baseline agent owns execution.'},
 {'id':'PRIOR-01','hypothesis':'Data identify luminosity evolution independently of cosmology and nuisance priors.','design':'Report likelihood profiles/degeneracies, proper priors and width sensitivity; separate physical host-age model from arbitrary redshift luminosity shift.','pass_condition':'If broad luminosity drift is degenerate with μ(z), report non-identifiability instead of a sigma claim.','status':'Preregistered.'},
 {'id':'COS-01','hypothesis':'Evidence about q0 is equivalent to evidence for evolving dark energy.','design':'Report q0<0, ΛCDM departure, and redshift-dependent acceleration separately; SN-only, BAO+free rd, early-universe priors, full CMB distinct.','pass_condition':'Correct estimand/model/priors for every result; external evidence never labeled SN-only.','status':'Preregistered.'},
 {'id':'OVL-01','hypothesis':'Comparison catalogues or model variants are independent datasets.','design':'Object-level crossmatch plus shared calibration/training/host-estimator dependency graph; use covariance or disjoint validation subsets.','pass_condition':'No double-counting of shared events or external datasets; unresolved aliases identified.','status':'Phase1 overlaps carried; phase2 raw sample audit independent.'}]
additional=[
 ('2606.19429','BayeSN with Dovekie calibration','v1','sources/updates/2026-09-20-ztf/originals/2606.19429v1.pdf','sources/updates/2026-09-20-ztf/extracted/2606.19429v1.txt','G26 training and DES validation conventions differ from a completed selection-corrected cosmology; training and test roles must be kept distinct.'),
 ('2605.06799','ZTF BayeSN environmental study','v3, 17 September 2026','sources/updates/2026-09-20-ztf/originals/2605.06799v3.pdf','sources/updates/2026-09-20-ztf/extracted/2605.06799v3.txt','New intrinsic-colour flexibility; 932 text objects, 929 residual rows, 944 linked master rows and 933 lite light curves need exact masks; calibration not yet cosmology ready.'),
 ('2601.19424','Union3.1 self-consistent host properties','2601.19424v1','sources/updates/2026-09-20-standardization/2601.19424.pdf','sources/updates/2026-09-20-standardization/2601.19424.txt','Host-mass revisions matter; claimed correction sign conflicts with separately released reconstruction and cannot be resolved without executed author vector.'),
 ('2607.24443','Pantheon+ host-mass correction reconstruction','v2, 15 September 2026','sources/updates/2026-09-20-standardization/2607.24443.pdf','sources/updates/2026-09-20-standardization/2607.24443.txt','114-row reconstruction reproduced in phase1; public new-minus-old sign conflicts with Hoyt textual sign; covariance held fixed.'),
 ('2608.02484','Sah response to acceleration/deceleration reanalysis','v1','runs/directional/sources/2608.02484.pdf','runs/directional/sources/2608.02484.txt','Response to simplified directional reanalysis; compare exact likelihood and coordinate convention rather than abstract conclusions.'),
 ('2607.20570','Ray revised acceleration reanalysis','v2, 3 August 2026','runs/directional/sources/2607.20570.pdf','runs/directional/sources/2607.20570.txt','v2 acknowledges coordinate error and withdraws original strong deceleration; author also cannot recover original full-sample baseline.'),
 ('2609.12083','Kim September host-age indicator study','v1','sources/updates/2026-09-20-ztf/originals/searches/2609.12083v1.pdf','sources/updates/2026-09-20-ztf/extracted/2609.12083v1.txt','Host-age indicator methods do not provide an independent released residual-age cosmological likelihood.'),
 ('2609.16972','Kelsey September host-SED study','v1','sources/updates/2026-09-20-ztf/originals/searches/2609.16972v1.pdf','sources/updates/2026-09-20-ztf/extracted/2609.16972v1.txt','Host SED/local versus global property inference does not directly observe a SN progenitor delay or independently establish an HR-age cosmology.')]
import copy
for id,title,version,pdf,txt,limit in additional:
 e=copy.deepcopy(entries['2604.16597']);e.update(id=id,title=title,source_url='https://arxiv.org/abs/'+id,exact_version=version,sources=[source(pdf,'https://arxiv.org/abs/'+id),source(txt)],claim=limit,methods_read=[{'scope':'phase1 method audit carried; no new phase2 full rereading','document':'docs/literature-and-data-status.md','details':limit}],reading_status='phase1-method-audit-carried',known_critiques=[limit],counterevidence=[],roles=['updated-source','phase1-continuity'])
 e['reproduction']={'status':'phase1-audit-carried','data_code':['docs/literature-and-data-status.md'],'exact_limits':[limit]};entries[id]=e

ledger=[]
for e in entries.values():
 for s in e['sources']:
  if s['path'] not in {r['path'] for r in ledger}:ledger.append(s)
for p in sorted((OUT/'sources').glob('*2fe0f56.c')):ledger.append(source(str(p.relative_to(ROOT)),'https://raw.githubusercontent.com/RickKessler/SNANA/2fe0f56/src/'+p.name.replace('-2fe0f56','')))
result={'schema_version':'1.0','built_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'DES raw-data likelihood dependencies, original and updated calibration, Pan-STARRS/Pantheon+, acceleration-critical arguments and counterarguments, DESI/CMB; broad source map with explicit reading depth, not exhaustive publication census.','selection_contract':'phase2/literature/selection-contract.json','entries':list(entries.values()),'preregistered_tests':tests,'reading_status_counts':{s:sum(e['reading_status']==s for e in entries.values()) for s in sorted({e['reading_status'] for e in entries.values()})}}
(OUT/'literature-map.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
(OUT/'source-hash-ledger.json').write_text(json.dumps(ledger,indent=2)+'\n')
print(len(entries),result['reading_status_counts'])
