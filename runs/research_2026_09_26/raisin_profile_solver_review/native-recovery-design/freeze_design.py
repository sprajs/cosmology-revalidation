"""Freeze metadata-only timing pilot; no native calls or distance scoring."""
from pathlib import Path
import csv, gzip, hashlib, json

ROOT = Path(__file__).resolve().parents[4]
OUT = Path(__file__).resolve().parent
BASE = ROOT / 'runs/research_2026_09_26/raisin_simulation_execution_review'
TAB = BASE / 'nominal_fitres'

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def table(p, marker):
    opener = gzip.open if p.suffix == '.gz' else open
    cols = None
    with opener(p, 'rt') as f:
        for line in f:
            if line.startswith('VARNAMES:'):
                cols = line.split()[1:]
            elif line.startswith(marker + ':'):
                vals = line.split()[1:]
                assert len(vals) == len(cols)
                yield dict(zip(cols, vals))

def dump(p, x):
    assert not p.exists(), p
    p.write_text(json.dumps(x, indent=2) + '\n')

lcpaths = [p for p in (TAB / 'lcplot-feasibility').rglob('*.gz')]
assert len(lcpaths) == 1, lcpaths
lcpath = lcpaths[0]
nir = {r['CID']: r for r in table(TAB / 'nir.FITRES.gz', 'SN')}
joint_ids = {r['CID'] for r in table(TAB / 'optnir.FITRES.gz', 'SN')}
epochs = {}
for r in table(lcpath, 'OBS'):
    if r['DATAFLAG'] == '1':
        epochs.setdefault(r['CID'], []).append(r)
ids = [str(i) for i in range(1, 501)]
assert set(ids) == set(epochs)
ledger = []
for cid in ids:
    r = nir[cid]
    ledger.append({'CID': cid, 'in_optnir': cid in joint_ids,
                   'n_J': sum(e['BAND'] == 'J' for e in epochs[cid]),
                   'n_H': sum(e['BAND'] == 'H' for e in epochs[cid]),
                   'zHEL': r['zHEL'], 'zCMB': r['zCMB'], 'zHD': r['zHD'],
                   'MWEBV': r['MWEBV'], 'SIM_LIBID': r['SIM_LIBID'],
                   'pilot': int(cid) <= 8})
with (OUT / 'all500-membership.csv').open('x', newline='') as f:
    w = csv.DictWriter(f, fieldnames=ledger[0]); w.writeheader(); w.writerows(ledger)
sources = [Path(__file__), TAB / 'nir.FITRES.gz', TAB / 'optnir.FITRES.gz', lcpath,
    TAB / 'acquisition.json', TAB / 'timing-audit-protocol.json',
    TAB / 'lcplot-feasibility/acquisition.json', TAB / 'lcplot-feasibility/inventory.json',
    ROOT / 'runs/research_2026_09_26/raisin_simulation_assets/timing/REFAC_DES_RAISIN_NIR.nml',
    ROOT / 'runs/research_2026_09_26/raisin_simulation_assets/timing/acquisition.json',
    BASE / 'author/fit/sim/REFAC_DES_RAISIN_optnir.nml',
    BASE / 'author/output/fit_nir/DES_RAISIN_NIR_SIM/SUBMIT.INFO',
    BASE / 'author/output/fit_all/DES_RAISIN_OPTNIR_SIM/SUBMIT.INFO',
    BASE / 'author/SIMLOGS_DES_SPEC/SIMGEN_MASTER_DESSPEC.INPUT',
    ROOT / 'runs/research_2026_09_26/raisin_sign_source/sim/simlibs/DES_RAISIN.simlib',
    ROOT / 'sources/repos/RickKessler__SNANA@v11_04k/src/snlc_fit.car',
    ROOT / 'sources/repos/RickKessler__SNANA@v11_04k/src/snana.car',
    ROOT / 'sources/repos/RickKessler__SNANA@v11_04k/src/sntools_output_text.c',
    ROOT / 'phase2/official/build/SNANA-v11_04k/bin/snlc_fit.exe',
    ROOT / 'sources/repos/djones1040__RAISIN_DataRelease@a383c4b/kcor/kcor_DES_NIR.fits',
    OUT / 'all500-membership.csv']
model = ROOT / 'sources/repos/djones1040__RAISIN_DataRelease@a383c4b/model/snoopy.B18'
sources += sorted(p for p in model.rglob('*') if p.is_file())
p = {
  'state': 'Frozen before reconstructed baseline or paired timing distance outcomes; design author runs no native fits.',
  'purpose': 'Bounded approximate archived-NIR input replay and conditional timing intervention; no population or cosmology correction.',
  'pilot_CIDs': ids[:8], 'all_available_CIDs': ids,
  'missing_optnir_in_all500': [cid for cid in ids if cid not in joint_ids],
  'missing_optnir_global': sorted(set(nir) - joint_ids, key=int),
  'selection': 'All first500 remain ledger. First8 numeric CIDs only are pilot, chosen by processing order, no replacement on failure. No paired response for other492 without later expansion protocol.',
  'native': {'binary': 'phase2/official/build/SNANA-v11_04k/bin/snlc_fit.exe',
    'nominal_nml': 'runs/research_2026_09_26/raisin_simulation_assets/timing/REFAC_DES_RAISIN_NIR.nml',
    'model': str(model.relative_to(ROOT)), 'MWlaw': 'historical default94; commented99 stays inactive',
    'parameters': 'stretch1, AV0, RV1.518, fixed peak; distance only free; NFIT_ITERATION3, native covariance feedback retained',
    'paths': 'Work-local NML and short relative TEXTFILE_PREFIX/private/version/override paths; hash every copy; no pinned source edits.',
    'threads': 1},
  'reconstruction': [
    'Only DATAFLAG1 J/H data rows. Preserve signs, error tokens, order and multiplicity. DATAFLAG0 and SIM_FLUXCAL never become measured photometry.',
    'Separate full original cadence unknown from exported accepted rows. No invented rejected epochs; any PHOTFLAG/default header value synthesized must be labelled and source-gated.',
    'Map observer phase/redshift from zHEL; retain zCMB,zHD,VPEC separately for distance-reference bookkeeping. Verify actual native REDSHIFT_FIT and output values rather than assigning zHD as the spectral redshift.',
    'Use archived printed MWEBV; preserve RV_MW3.1 and nominal USE_MWCOR false/OPT_MWEBV1. SIM_LIBID→SIMLIB coordinates may supply RA/DEC only with exact unique block join; do not derive dust from sky maps or substitute observed SN metadata.',
    'Baseline PEAKMJD is archived PKMJDINI. Record header parsing float32 coordinates and printed PKMJD; do not reconstruct MJD from the coarser Tobs column.',
    'SNFIT LCPLOT error is JEP_DATAFLUX_ERR=preprocessed SNLC_FLUXCAL_ERRTOT×FUDGE_DATAERR_SCALE, not full model/MW covariance. Default scale1 and exact NML. Trace preprocessing/default maps; feed exported error only if no duplicated transform is applied; verify outgoing data-error intervals. Do not subtract model error or set fullC=diag(exportederror²).',
    'Before any baseline score, freeze execution script, generated NML/header files, exact model/KCOR/binary hashes and all allowed path-only differences. If metadata/error convention cannot be source-closed, stop.'
  ],
  'printed_input_intervals': {
    'MJD_days': 0.0005, 'peak_days': 0.00005,
    'flux_error_and_header_tokens': 'For each printed decimal/scientific token, half one unit in its last displayed digit, using Decimal; never a relative five-digit blanket approximation.',
    'source_storage': 'Additionally record REAL*4 epoch-offset and peak parsing; use source float32 round-trip. LCPLOT MJD is offset REAL*4 plus MJDOFF, so decimal bounds alone do not assert original full precision.',
    'uncertain_variables': 'Each epoch MJD, FLUXCAL, FLUXCAL_ERR; header peak, zHEL,zCMB,MWEBV. VPEC/zHD source map kept distinct; no physical measurement uncertainty added.',
    'screen': 'After baseline, for first8 only test each nonzero rounding coordinate at ±half-last-digit, plus joint directions maximizing/minimizing local D and dataQ changes. Use paired identical source transformations and output quantization intervals. Report the L1 first-order envelope and measured nonlinear remainder as a sensitivity screen, not a rigorous all-corners interval proof. Stop if this exceeds fixed engineering caps or if any masks change.'
  },
  'baseline_gate': {
    'copies': 'Two isolated baseline directories with byte-identical reconstruction: require bitwise identical FITRES SN rows and LCPLOT output after excluding explicit path/comment metadata only (every exclusion frozen first).',
    'all_eight': 'Require finite output, ERRFLAG0, exact archived NDOF, accepted J/H (band,printedMJD) multiset equality, source float32 fixed parameter closure and zero fixed-parameter errors. Failures stay in ledger; no replacement.',
    'data_error': 'Outgoing data flux/error must agree with archived token intervals plus explicitly calculated source float32 round-trip interval; no extra fitted tolerance.',
    'DLMAG_absolute_cap_mag': 0.001,
    'FITCHI2_absolute_cap': 'max(0.01,0.001*abs(archived FITCHI2)); exported FITCHI2 is dataQ, not full minimand.',
    'meaning': 'These are predeclared engineering closure caps, NOT a derived rounding bound or proof of historical execution. Both baseline discrepancies and independent rounding screen must pass. Never widen caps after outcomes.',
    'native_state': 'Record iteration log, epoch masks, fixed peak and covariance settings. No ranking across state-dependent objectives by dataQ alone.'
  },
  'paired_after_gate': {
    'arms': ['INI: archived neartruth PKMJDINI', 'JOINT: sameCID archived optical+NIR fitted PKMJD'],
    'change': 'Only PEAKMJD header. Identical reconstructed NIR flux/error/epoch vector, redshift, MW, model and settings. No new noise/scatter; no optical shape/AV import.',
    'primary_estimand': 'Perobject deltaD=D_JOINT−D_INI, then unweighted first8 mean and empirical dispersion; all native failures/mask changes preserved. This is conditional replay, not a population bias estimator.',
    'mask_gate': 'Report native cut-aware result only on available exported rows. If mask differs between arms, retain it but do not mix into common-mask summary or claim full-cadence counterfactual; separate new fixed-mask amendment would be required.',
    'rounding_stability': 'Replay the jointD-oriented signed endpoint rounding directions in both arms to screen stability of deltaD; archive peak token uncertainty is included in each arm separately.',
    'shared_photometry_assumption': 'Same VERSION_PHOTOMETRY, CID and printed truth support pairing but do not prove crossbranch epoch/noise identity. JOINT timing may share JH noise with INI; retain this assumed coupling, never treat JOINT time as exogenous.',
    'convergence': 'Do not equate ERRFLAG0 with uniqueness. With fixed shape/color/time, record true objective/prior state and verify baseline copies. Distinct amplitude starts are meaningful only if logs confirm distinct post-initialization starts; otherwise mark untested instead of claiming multistart. A larger run requires an audited one-dimensional distance profile if instability appears.'
  },
  'resource_cap': {'single_worker_wall_seconds': 600,
    'stage_rule': 'Batch deterministic rounding clones to avoid process initialization per scalar perturbation. Checkpoint source gate, baseline, rounding, paired stage separately. No 500-object expansion. If cap reached, preserve completed states and report incomplete gates, not pass.',
    'forecast': '16 baseline object fits plus at most~400 rounding/control object fits and8 JOINT; actual pilot row counts fix exact count before execution. Native process overhead must be amortized. Stop after first process if forecast exceeds cap.'},
  'interpretation': [
    'Archived joint peak itself starts near truth, defaults to initial −4:2:+4d grid search, and uses5d prior recentered after each prior fit. It is a source-supported simulated alternative, not reconstruction of the unknown observed header creator.',
    'Active generator OIR.J19 scatter is already in the same simulated light curves. Paired timing response may correlate with existing distance residuals/scatter; do not add its RMS independently or call it an omitted physical scatter term.',
    'No SIM_DLMAG residual or MSE scoring before a separate truth/scatter definition gate. Future variance difference=Var(deltaD)+2Cov(D_INI−truth,deltaD); MSE difference=2E[(D_INI−truth)deltaD]+E[deltaD²].',
    'First500 is an output-order cap; first8 is an engineering pilot. Missing thirteen global joint fits and CID91 in available500 are retained. No population/cosmology extrapolation.'
  ],
  'native_sign_recovery': 'Defer. Preserve old1536-nativefit proposal unchanged. A distinct finite-library experiment may use only deterministic521×65 grid core, not observed-adaptive nodes; source-chosen truths need fresh native verification. Details in companion report.',
  'inputs_sha256': {str(x.relative_to(ROOT)): sha(x) for x in sources}
}
dump(OUT / 'timing-pilot-protocol.json', p)
print(json.dumps({'protocol_sha256': sha(OUT / 'timing-pilot-protocol.json'),
                  'pilot': ids[:8], 'missing_in500': p['missing_optnir_in_all500'],
                  'input_files': len(sources)}))
