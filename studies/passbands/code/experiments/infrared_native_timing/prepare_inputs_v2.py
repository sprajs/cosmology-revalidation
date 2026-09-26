from pathlib import Path
import json,hashlib,re,csv,shutil,numpy as np
R=Path('/home/szymon/Documents/ChatGPT/supernova');P=R/'runs/research_2026_09_26/raisin_profile_solver_review/csp-filter-interpretation/native-design/timing-followon-design/prospective-engineering';I=P/'inputs-v2';I.mkdir(exist_ok=False)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
src=R/'runs/research_2026_09_26/raisin_sign_source/sim/simlibs/DES_RAISIN.simlib';s=src.read_text();b=re.search(r'^LIBID:\s+11\s.*?^END_LIBID:\s*11.*?$',s,re.M|re.S).group();p=float(np.float32(57707.8));z=float(np.float32(.453));rv=float(np.float32(1.518))
selected=[];ledger=[];ordinal=0
for ln,line in enumerate(s.splitlines(),1):
 if line not in b.splitlines() or not line.startswith('S:'):continue
 # This source uses unique exposure lines; assert exact occurrence below.
 t=line.split();phase=(float(t[1])-p)/(1+z);keep=-15<=phase<=45;ordinal+=1
 ledger.append(dict(source_line=ln,block_ordinal=ordinal,exposure_id=t[2],band=t[3],mjd=t[1],truth_rest_phase=phase,keep=keep,reason='supported metadata cadence' if keep else 'outside prospective truth phase window'))
 if keep:selected.append(line)
assert ordinal==735 and len(selected)==117
head=s.split('BEGIN LIBGEN')[0];head=re.sub(r'NLIBID:\s+\d+','NLIBID: 1',head)
bhead='\n'.join(list(__import__('itertools').takewhile(lambda line:not line.startswith('S:'),b.splitlines())))+'\n';bhead=re.sub(r'NOBS:\s+735','NOBS: 117',bhead);bhead=re.sub(r'REDSHIFT:\s+\S+\s+PEAKMJD:\s+\S+',f'REDSHIFT: {z:.17g}   PEAKMJD: {p:.17g}',bhead)
(I/'cadence.simlib').write_text(head+'BEGIN LIBGEN\n\n'+bhead+'\n'.join(selected)+'\nEND_LIBID: 11\n\nEND_OF_SIMLIB: 1 ENTRIES\n')
with (I/'cadence-ledger.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(ledger[0]));w.writeheader();w.writerows(ledger)
model=R/'sources/repos/djones1040__RAISIN_DataRelease@a383c4b/model/snoopy.B18';kcor=R/'sources/repos/djones1040__RAISIN_DataRelease@a383c4b/kcor/kcor_DES_NIR.fits';emap=R/'runs/research_2026_09_26/raisin_simulation_assets/author/sim/inputs/DES/DES3YR_SIM_ERRORFUDGES.DAT'
assert all(x.exists() for x in [model,kcor,emap])
(I/'snoopy.B18').symlink_to(model,target_is_directory=True);(I/'kcor.fits').symlink_to(kcor);(I/'error-map.dat').symlink_to(emap)
base=f'''# Controlled prospective experiment, not historical simulation reproduction.
SIMLIB_FILE: ../../inputs-v2/cadence.simlib
SIMLIB_NREPEAT: 8
SIMLIB_IDLOCK: 1
USE_SIMLIB_PEAKMJD: 1
USE_SIMLIB_REDSHIFT: 1
NGEN_LC: 0
NGENTOT_LC: {{attempts}}
CIDOFF: 0
GENVERSION: PTE
GENSOURCE: RANDOM
GENMODEL: ../../inputs-v2/snoopy.B18
GENFILTERS: grizJH
KCOR_FILE: ../../inputs-v2/kcor.fits
GENMAG_SMEAR_MODELNAME: NONE
GENMAG_SMEAR: 0
GENMODEL_ERRSCALE: 0
WEAKLENS_DSIGMADZ: 0
GENRANGE_PEAKMJD: 57700 57720
GENSIGMA_SEARCH_PEAKMJD: 0
GENRANGE_TREST: -15 45
GENRANGE_REDSHIFT: 0.4 0.5
GENSIGMA_REDSHIFT: 0
GENSIGMA_VPEC: 0
VPEC_ERR: 0
VEL_CMBAPEX: 0
OPT_MWEBV: 1
OPT_MWCOLORLAW: 94
RV_MWCOLORLAW: 3.1
GENSIGMA_MWEBV_RATIO: 0
GENSIGMA_MWEBV: 0
GENSHIFT_MWEBV: 0
GENSCALE_MWEBV: 1
RANSEED: 26092671
SMEARFLAG_FLUX: {{smear}}
SMEARFLAG_ZEROPT: 0
FLUXERRMODEL_FILE: ../../inputs-v2/error-map.dat
APPLY_SEARCHEFF_OPT: 0
SEARCHEFF_SPEC_FILE: NONE
SEARCHEFF_PIPELINE_FILE: NONE
SEARCHEFF_PIPELINE_LOGIC_FILE: NONE
APPLY_CUTWIN_OPT: 0
CUTWIN_NEPOCH: 1 -999
NEWMJD_DIF: 0.4
OMEGA_MATTER: 0.3
OMEGA_LAMBDA: 0.7
W0_LAMBDA: -1
H0: 70
FORMAT_MASK: 32
DNDZ: POWERLAW 1.7E-5 2.11
SIMGEN_DUMP: 7 CID ZCMB STRETCH AV DLMAG PEAKMJD LIBID
GENPEAK_STRETCH: 1
GENRANGE_STRETCH: 1 1
GENSIGMA_STRETCH: 0 0
GENPEAK_RV: {rv:.17g}
GENRANGE_RV: {rv:.17g} {rv:.17g}
GENSIGMA_RV: 0 0
GENRANGE_AV: 0 0
GENTAU_AV: 0.2
PATH_SNDATA_SIM: {{output}}
'''
for branch in ['original','ledger','noiseless']:
 w=P/'generation-v2'/branch;w.mkdir(parents=True,exist_ok=False);(w/'output').mkdir()
 (w/'sim.input').write_text(base.format(attempts=1 if branch=='noiseless' else 8,smear=0 if branch=='noiseless' else 1,output=str(R/'phase2/pte/generation-v2'/branch/'output')))
# Source-derived fitter template; actual generated native FITS files are read directly.
fit='''&SNLCINP
 PRIVATE_DATA_PATH = '{data_path}'
 VERSION_PHOTOMETRY = 'PTE'
 KCOR_FILE = '../../inputs-v2/kcor.fits'
 NFIT_ITERATION = {iterations}
 INTERP_OPT = 1
 MXLC_PLOT = 8
 USE_MINOS = F
 SNTABLE_LIST = 'FITRES(text:key) LCPLOT(text:key)'
 TEXTFILE_PREFIX = 'fit'
 LDMP_SNFAIL = T
 USE_MWCOR = F
 H0_REF = 70
 OLAM_REF = 0.7
 OMAT_REF = 0.3
 W0_REF = -1
 SNCID_LIST = 0
 CUTWIN_CID = 0, 20000000
 cutwin_redshift = 0.01, 2
 cutwin_Nepoch = 1
 RV_MWCOLORLAW = 3.1
 OPT_MWCOLORLAW = 94
 OPT_MWEBV = 1
 MWEBV_SCALE = 1
 MWEBV_SHIFT = 0
 EPCUT_SNRMIN = ''
 ABORT_ON_NOEPOCHS = F
&END
&FITINP
 FITMODEL_NAME = '../../inputs-v2/snoopy.B18'
 OPT_PRIOR_AV = 0
 PRIOR_MJDSIG = 5
 PRIOR_LUMIPAR_RANGE = -5, 5
 INIVAL_SHAPE = 1
 INISTP_SHAPE = 0
 INIVAL_AV = 0
 INISTP_AV = 0
 INIVAL_RV = 1.518
 INISTP_RV = 0
 INISTP_DLMAG = 0.5
 INISTP_PEAKMJD = {peak_step}
 INIVAL_PEAKMJD_SHIFT = {peak_start_shift}
 OPT_XTMW_ERR = 1
 OPT_COVAR_FLUX = 0
 TREST_REJECT = -20, 70
 NGRID_PDF = 0
 FUDGEALL_ITER1_MAXFRAC = 0.02
 FILTLIST_FIT = '{bands}'
&END
'''
(I/'fit-template.nml').write_text(fit)
# schema terminology amendment, binary unchanged
ip=P/'instrumentation';sc=json.loads((ip/'schema.json').read_text());sc['PROSP_NOISE'][sc['PROSP_NOISE'].index('gauss_S')]='gauss_SZ';(ip/'schema-v2.json').write_text(json.dumps(sc,indent=2)+'\n')
d={'source_simlib':str(src.relative_to(R)),'source_simlib_sha256':sha(src),'chosen_LIBID':11,'source_block_rows':735,'selected_rows':117,'selected_band_counts':{b:sum(x['keep'] and x['band']==b for x in ledger) for b in 'grizJH'},'source_peak':57707.8,'generation_peak':p,'peak_shift_days':p-57707.8,'source_redshift':.453,'generation_redshift':z,'VEL_CMBAPEX':0,'redshift_change_reason':'Representable native-R4 controlled redshift and removal of frame/velocity nuisance; not historical reproduction','RV':rv,'selection_phase_rest_days':[-15,45],'fit_guard_phase_days':[-20,70],'model_empirical_phase_support':[-20,70],'filter_or_epoch_value_edits':'Exposure lines unchanged. Only one block selected; NOBS/NLIBID/peak/redshift/global end count updated. No photometry outcomes exist.','noise_map_scope':'Unchanged author map retained; actual full FIELD DES16E1dcx matches no map, and J/H absent. Must verify zero F variance. No field remapping.','files':{str(f.relative_to(R)):sha(f) for f in [I/'cadence.simlib',I/'cadence-ledger.csv',I/'fit-template.nml',ip/'schema-v2.json',kcor,emap]+list(model.glob('*'))+list((P/'generation-v2').glob('*/sim.input')) if f.is_file()}}
(I/'manifest.json').write_text(json.dumps(d,indent=2)+'\n');print(json.dumps({'inputs_manifest_sha256':sha(I/'manifest.json'),'files':len(d['files'])}))
