#!/usr/bin/env python3
"""Prepare predeclared DES calibrated-flux injections with all-generated DUMPs.
The P21 run uses released 2024 assets, not an asserted byte-exact 2022 replay.
"""
from pathlib import Path
import argparse, hashlib,json,re,datetime,gzip
ROOT=Path(__file__).resolve().parents[3]
SNDATA=ROOT/'phase2/official/inputs/SNDATA_ROOT'
OUT=ROOT/'phase2/literature/simulations'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write_immutable(p,text):
 p.parent.mkdir(parents=True,exist_ok=True)
 if p.exists() and p.read_text()!=text:raise RuntimeError(f'Refusing changed existing input {p}')
 if not p.exists():p.write_text(text)

parser=argparse.ArgumentParser();parser.add_argument('--n',type=int,default=26518);parser.add_argument('--seed',type=int,default=213849);parser.add_argument('--tag',default='pilot01');parser.add_argument('--compat-alpha',action='store_true');args=parser.parse_args()
base=SNDATA/'sample_input_files/DES-SN5YR/base_files/sim'
survey=(base/'sim_des5yr_survey.input').read_text().split('SIMGEN_DUMP:')[0]
survey=re.sub(r'^RANSEED:.*$',f'RANSEED: {args.seed}',survey,flags=re.M)
survey=re.sub(r'^OMEGA_MATTER:.*$','OMEGA_MATTER: 0.315',survey,flags=re.M)
survey=re.sub(r'^OMEGA_LAMBDA:.*$','OMEGA_LAMBDA: 0.685',survey,flags=re.M)
ia=(base/'sim_ia_salt_des5yr.input').read_text()
ia=re.sub(r'^GENPDF_FILE:.*$','',ia,flags=re.M)
ia=ia.replace('$PLASTICC_MODELS/SNIa_Extrap_LateTime_2expon.TEXT',str(ROOT/'phase2/literature/sources/SNIa_Extrap_LateTime_2expon.TEXT'))
cols='CID GENZ GALZTRUE LIBID RA DEC MWEBV MU PEAKMJD PEAKMAG_g PEAKMAG_r PEAKMAG_i PEAKMAG_z SNRMAX_g SNRMAX_r SNRMAX_i SNRMAX_z SNRMAX SNRMAX2 SNRMAX3 NOBS TRESTMIN TRESTMAX TGAPMAX CUTMASK SIM_EFFMASK GALID LOGMASS_TRUE LOGSFR_TRUE r_obs_auto obs_gr_auto SALT2mB SALT2x1 SALT2c SALT2alpha SALT2beta AV RV MAGSMEAR_COH'.split()
models={
 'P21':{'map':'P21','extra':'','limits':['Released 2024 renamed survey/population assets; source equivalence to original 2022 mock inputs not established.']},
 'P21_dmplus010':{'map':'P21','extra':'GENMAG_OFF_GLOBAL: 0.10\n','limits':['Physical flux injection +0.10 mag in addition to population MAG_OFFSET; not importance-reweighted source magnitudes.']},
 'P21_dmminus010':{'map':'P21','extra':'GENMAG_OFF_GLOBAL: -0.10\n','limits':['Physical flux injection -0.10 mag in addition to population MAG_OFFSET; not importance-reweighted source magnitudes.']},
 'P21_dmz020':{'map':'P21','extra':'ZVARIATION_POLY: GENMAG_OFF_GLOBAL 0,0.20\n','limits':['Physical flux injection +0.20*z mag; fixed cosmology retained only for generating survey flux.']},
 'BS21':{'map':'BS21','extra':'GENPEAK_SALT2ALPHA: 0.15\nGENPEAK_SALT2BETA: 1.98\nGENSIGMA_SALT2BETA: 0.35 0.35\nGENRANGE_SALT2BETA: 0.4 3.0\nGENPEAK_SALT2c: -0.084\nGENSIGMA_SALT2c: 0.042 0.042\nGENRANGE_SALT2c: -0.5 0.5\nGENMAG_OFF_GLOBAL: -0.12\n','limits':['Released approximate BS21 RV/EBV/x1 grids plus published intrinsic colour/beta moments (arXiv2004.10206v2, section4).','Beta and colour truncation ranges and common -0.12 offset are explicit controlled choices; exact historical BS21 job configuration not recovered.','Use free magnitude offsets in validation; do not interpret global magnitude difference as cosmology.']},
 'G10':{'map':'G10','extra':'GENMAG_SMEAR_MODELNAME: G10\nGENPEAK_SALT2ALPHA: 0.15\nGENPEAK_SALT2BETA: 3.1\n','limits':['Released G10 mass-dependent x1/c grids, G10 wavelength scatter and beta3.1 from arXiv2004.10206v2 section3.1; common alpha0.15.','Controlled public-assets variant, not certified exact original G10 systematic job.']},
}
runs=[]
compat_map=None
if args.compat_alpha:
 src=SNDATA/'models/population_pdf/DES-SN5YR/DES-SN5YR_DES_S3_P21.DAT.gz'
 compat_map=OUT/'inputs'/'DES-SN5YR_DES_S3_P21_constant_alpha_relocated.DAT'
 original=gzip.open(src,'rt').read()
 derived=re.sub(r'^GENPEAK_SALT2ALPHA:.*$', '# Constant alpha0.15 relocated to input to avoid historical zero-sigma map bug.', original,flags=re.M)
 write_immutable(compat_map,derived)
 write_immutable(OUT/'inputs'/'constant-alpha-compatibility.json',json.dumps({'source':str(src.relative_to(ROOT)),'source_sha256':sha(src),'derived':str(compat_map.relative_to(ROOT)),'derived_sha256':sha(compat_map),'change':'Remove GENPEAK_SALT2ALPHA0.15 from PDF and set same exact constant in simulation input. Historical2023 generator fails constructing zero-sigmaGaussian grid; no added scatter.','diagnostic':'phase2/official/diagnostics/simulation-segfault-gdb.log'},indent=2)+'\n')
for name,model in models.items():
 version=f'PH2_{args.tag}_{name}'
 inp=OUT/'inputs'/f'{version}.input'
 text=f'# Prepared before cosmology comparison; all generated attempts requested in DUMP.\nGENVERSION: {version}\nGENPREFIX: {version}\nPATH_SNDATA_SIM: {OUT/"outputs"}\nFORMAT_MASK: 32\nNGENTOT_LC: {args.n}\n'+survey+'\n'+ia+'\n'+f'GENPDF_FILE: $SNDATA_ROOT/models/population_pdf/DES-SN5YR/DES-SN5YR_DES_S3_{model["map"]}.DAT\n'+model['extra']+'HOSTLIB_STOREVAR: LOGMASS_TRUE,LOGSFR_TRUE,r_obs_auto,obs_gr_auto\n'+f'SIMGEN_DUMPALL: {len(cols)}\n  '+' '.join(cols)+'\n'
 if compat_map and model['map']=='P21':
  text=text.replace('$SNDATA_ROOT/models/population_pdf/DES-SN5YR/DES-SN5YR_DES_S3_P21.DAT',str(compat_map))+'GENPEAK_SALT2ALPHA: 0.15\n'
 write_immutable(inp,text)
 runs.append({'name':name,'version':version,'input':str(inp.relative_to(ROOT)),'input_sha256':sha(inp),'generated_requested':args.n,'seed':args.seed,'population_map':f'DES-SN5YR_DES_S3_{model["map"]}.DAT.gz','limits':model['limits']})
prereg={'purpose':'Forward model/selection validation at calibrated epoch flux level; no cosmological parameter estimation.','sample_design':'One initial common-seed realization per model; fixed number of generated attempts, not fixed accepted count.','runs':runs,'W22':{'status':'blocked-exact-input-missing','reason':'W22 x1 PDF needs SN_age; nominal released DES host library lacks SN_age. No substitute age is invented.'},'pairing':'Same seed and survey source does not guarantee identical host/cadence draws when model RNG consumption differs. Audit CID/redshift/LIBID/host equality before claiming paired uncertainty reduction.','selection_stages':['Generated DUMP attempt','Pipeline+hostz+minimal generated cuts, FITS HEAD membership','Identical calibrated-flux SALT fitter','Predeclared measured quality subset','Any additional final analysis membership separately declared'],'heldout_tests':['Compare joint fitted c,x1 and host mass within redshift/field bins on held-out CID split.','Compare within-z magnitude-colour-host patterns after estimating one free magnitude intercept per z bin on training split.','Use same measured quality cuts and covariance handling on observations and simulations.','Report simulated counts and multinomial/Monte Carlo uncertainties; one realization is a pilot, not decisive model rejection.','Measure selection response to global plus/minus0.10 and redshift0.20*z luminosity injections; no singular-proposal importance approximation.'],'nonclaims':['No independent SN_age measurement from mean host age.','No validation from preferred LCDM fit; injection cosmology is a flux-generation input only.','No final-cosmology or complete contamination validation from Ia-only mocks.']}
write_immutable(OUT/'inputs'/f'preregistration-{args.tag}.json',json.dumps(prereg,indent=2)+'\n')
(OUT/'outputs').mkdir(exist_ok=True);(OUT/'logs').mkdir(exist_ok=True)
print(json.dumps({'runs':len(runs),'outputs':str(OUT/'outputs'),'preregistered':str(OUT/'inputs'/f'preregistration-{args.tag}.json')}))
