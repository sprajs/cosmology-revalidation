"""Metadata-only CSP Jdw relabeling design. No fits or distance outcomes."""
import csv, gzip, hashlib, json, re
from decimal import Decimal
from pathlib import Path

ROOT=Path(__file__).resolve().parents[5]
OUT=Path(__file__).resolve().parent
R=ROOT/'runs/research_2026_09_26'
REL=ROOT/'sources/repos/djones1040__RAISIN_DataRelease@a383c4b'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,x):
    assert not p.exists(),p
    p.write_text(json.dumps(x,indent=2)+'\n')

membership=R/'raisin_differential/frozen-membership.csv'
lineage=R/'csp_magnitude_lineage/released-physical-filter-ledger.csv'
members=[x for x in csv.DictReader(membership.open()) if x['survey']=='CSP']
assert len(members)==42 and len({x['CID'] for x in members})==42
lines=list(csv.DictReader(lineage.open()))
cols=None;raw={}
with gzip.open(OUT/'CSP-nominal.FITRES.gz','rt') as f:
    for line in f:
        if line.startswith('VARNAMES:'): cols=line.split()[1:]
        if line.startswith('SN:'):
            row=dict(zip(cols,line.split()[1:]));assert row['CID'] not in raw;raw[row['CID']]=row
assert all(m['CID'] in raw for m in members)
roster=[];changes=[];ambiguous=[];phot_paths=[]
for m in members:
    path=ROOT/m['photometry_path'];phot_paths.append(path)
    text=path.read_text();peaks=re.findall(r'^PEAKMJD:\s+(\S+)',text,re.M)
    assert peaks and len(set(peaks))==1
    n=inside=0
    for x in lines:
        if x['CID']!=m['CID']: continue
        if x['raw_physical_filter']=='Jdw' and x['physical_unique']=='True' and x['converter_identity']=='True':
            assert x['band']=='J'
            count=int(x['multiplicity']);t=Decimal(x['time'])+Decimal(53000)
            phase=(float(t)-float(peaks[0]))/(1+float(m['zHEL']))
            n+=count;inside+=count*(-15<=phase<=45)
            changes.append({'CID':m['CID'],'absolute_MJD':str(t),'old_band':'J','new_band':'j','multiplicity':count,'header_phase':phase,'inside_header_window':-15<=phase<=45,'photometry_path':m['photometry_path']})
        elif x['physical_unique']!='True' and 'Jdw' in x['raw_physical_filter'].split('|'):
            ambiguous.append(x)
    roster.append({'CID':m['CID'],'n_known_Jdw_rows':n,'n_known_Jdw_header_window':inside,
                   'zHEL':m['zHEL'],'PEAKMJD':peaks[0], 'photometry_path':m['photometry_path'],
                   'archived_raw_present':m['CID'] in raw,'affected':n>0})
for name,rows in [('cohort42.csv',roster),('known-relabel-rows.csv',changes),('ambiguous-rows.csv',ambiguous)]:
    with (OUT/name).open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
affected=sorted(x['CID'] for x in roster if x['affected'])
unaffected=sorted(x['CID'] for x in roster if not x['affected'])
paths=[Path(__file__),membership,lineage,OUT/'CSP-nominal.FITRES.gz',OUT/'acquisition.json',
       ROOT/'phase2/official/build/SNANA-v11_04k/bin/snlc_fit.exe',
       REL/'lcfitting/REFAC_CSP_RAISIN_nir_sys.nml',REL/'kcor/kcor_CSPDR3_BD17.fits',
       REL/'kcor/kcor_CSPDR3_BD17.input',REL/'vpec/vpec_baseline_raisin.list',
       R/'raisin_differential/mass-threshold/code/output/fit_nir_sys/CSP_RAISIN/SUBMIT.INFO']
paths+=phot_paths+sorted(p for p in (REL/'model/snoopy.B18').rglob('*') if p.is_file())
p={
 'state':'Frozen metadata-only design before native baseline or changed-filter distance outcomes; no execution by design author.',
 'question':'Conditional raw NIR distance response to interpreting uniquely identified WIRC J rows in the paper RC2/WIRC natural-system approximation rather than released RC1 assignment.',
 'membership':'Exactly all42 final-selected CSP objects from previously frozen RAISIN differential membership; no new quality or phase selection.',
 'counts':{'cohort':42,'known_affected_objects':len(affected),'known_relabel_rows':sum(x['n_known_Jdw_rows'] for x in roster),'header_phase_affected_objects':sum(x['n_known_Jdw_header_window']>0 for x in roster),'header_phase_rows':sum(x['n_known_Jdw_header_window'] for x in roster),'ambiguous_Jdw_rows':sum(int(x['multiplicity']) for x in ambiguous)},
 'affected_CIDs':affected,'unaffected_CIDs':unaffected,
 'pilot_CIDs':[affected[0],unaffected[0]],
 'input_relabel':'Only listed physical-unique, converter-exact Jdw rows: J→j. Match exact Decimal MJD and oldband with full multiplicity against row ledger; fail ambiguity. Keep all other row tokens, header, epochs, MAG/ERR and FLUXCAL/ERR byte-identical. Relabel all188 known rows, including outside final phasewindow, preserving native initialization pathway.',
 'primary_ambiguity_rule':'Keep four2006kf J|Jdw rows asJ; do not assign hidden physical identity. Optional separately labelled allfour→j stress arm for2006kf, after primary. It bounds this discrete assignment choice only, not all calibration uncertainty.',
 'observation_operator':'j uses released RC2 throughput and BD17 natural magnitude8.4192, matching paper WIRC/RC2 grouping. This is a counterfactual operator choice, not a pure flux-unit transformation; no extra scalar magnitude or flux offset and no second S-correction.',
 'native':{'source':'Existingv11_04kbuild, matching acquired rawFITRES header; exact released NIR NML with path/output-format substitutions only; fixedstretch1 AV0 RV1.518 peak, empirical SNooPy.B18, sourceMWdefault94, native covariance iteration retained. Do not borrow v11_03c observedDES baseline.',
           'paths':'Short work-local relative output/version/override paths and byte-identical local VPEC copy; hash allgeneratedinputs/runner beforebaseline.',
           'mask':'Same available full photometry; require baseline and relabel accepted physical-row multisets equal after undoing onlyJ→j in identity. No automatic row deletion. Header phase counts in roster are metadata, not certified native masks.'},
 'baseline_gates':{'pilot':'Firstaffected2004ef and firstunaffected2005hc, two isolated identical nominal copies before relabel. Require bitwise FITRES science rows/LCPLOT, ERRFLAG0, fixedparameterclosure.',
                   'archive':'Compare native rawD to acquired uncorrected FITOPT000, never cosmology bias/mass-corrected table. Require absoluteDdifference≤0.001mag, exact archivedNDOF and dataQdifference≤max(.01,.001*abs(archivedQ)); caps frozen before comparisons, no widening. Record allredshift/MW/peak/header agreement.',
                   'support':'Verify actualnative evaluated rest phases/filter wavelengths and mean clipping, including iteration1; preserve any extrapolationflag. FixedC amplitude proof only where native means scale multiplicatively and stay within valid flux/magnitude support.',
                   'full_cohort':'Only after pilotpass and runtimeforecast: all42nominal+nominalcopy. Any objectfailure remains in ledger; no silent reducedsample. Quantitative42anchor summary waits for all42baselinegates.'},
 'optimizer_gate':{'meaning':'ERRFLAG0 and positivecovariance alone insufficient. One freeD is simpler than failed multishape profiles, but still has native state-dependentC.',
                   'test':'For each pilot arm require actual distinct post-initialization amplitude starts separated by±.2mag around archivedrawD, if source supports this. Logs must show distinct objective entry coordinates. If native initialization erases differences, do not call it multistart.',
                   'fixed_state_check':'At final native mask/C, independently verify positive-amplitude quadratic optimum and nativeD±.15 mean multiplicativity with actualfloat32coordinates; export actualnativeobjective and covariance without SALT x0 conversion. FittedD versus fixed-state optimum≤.001mag and independentQ closure≤1e−6. This validates final-state objective only; require lasttwo iterationDchange≤.001mag and compare covariance states.',
                   'failure':'If any needed exporter/gate is unavailable, freeze a separate sourceinstrumentation/profile amendment before newresponse outcomes. Do not run a broad automatic multishape solver or conceal untestedstability.'},
 'arms':['nominal','nominal_copy','known_Jdw_to_j','optional2006kf_all_four_ambiguous_to_j'],
 'estimands':['Perobject deltaD=rawD_relabel−rawD_nominal withmask/support/stateflags.','Unweighted mean among32knownaffected; separately 42-object mean retaining10zero-change controls.','Conditional fixed-sample high-minus-low distance contrast response is minus42-objectmean deltaD if allhigh-z objects held unchanged. This is pre-BBC processing response, not a correction or cosmologyfit.','Phase and Jdw-row-count descriptors from frozenmetadata only; do not select based on response.'],
 'control':'Ten unaffected objects must be byte-identical input/output between nominal and relabel; any nonzeroD signals implementation/state nondeterminism and fails gate.',
 'resource':'Singleworker,600active seconds including pilot; preserve perprocess timing/checkpoints. Pilot2objects before all42; forecast fullnominalcopies+primary and stabilitygate. If exceeded report estimate/incomplete, no reducedoutcome-selected cohort.',
 'scope_limits':['RC2 is paper-supported approximation toWIRC, not exactWIRCthroughput. A physicalWIRC operator requires separately calibratedprimary and regeneratedKCOR.','Fixed empirical template and its training/calibration inheritance remain; no retraining or populationbias inferred.','Final42 membership is held fixed; detection/followup/lightcurve selection and nativecovariance generation changes are outside this counterfactual.','No cosmology or correctedsystematicmatrix before coherent bias simulations, selection and covariance regeneration. No addition of syntheticphaseRMS as independenterror.'],
 'inputs_sha256':{str(x.relative_to(ROOT)):sha(x) for x in paths}}
dump(OUT/'protocol.json',p)
print(json.dumps({'protocol_sha256':sha(OUT/'protocol.json'),'counts':p['counts'],'pilot':p['pilot_CIDs'],'n_raw_fitres':len(raw)}))
