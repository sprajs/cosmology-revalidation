from astropy.io import fits
from pathlib import Path
import json,hashlib
r=Path(__file__).resolve().parents[3];out=Path(__file__).resolve().parent/'noise_variant_design';base=r/'phase2/literature/simulations/outputs';arms=['P21','P21_rho000','P21_rho090','P21_noisetrue120'];h={}
for arm in arms:
 v='PH2_pilot02_'+arm;a=fits.getdata(base/v/(v+'_HEAD.FITS'),1);h[arm]={int(x['SNID']):x for x in a}
cols=['RA','DEC','MWEBV','MWEBV_ERR','REDSHIFT_HELIO','REDSHIFT_HELIO_ERR','REDSHIFT_FINAL','REDSHIFT_FINAL_ERR','VPEC','VPEC_ERR','HOSTGAL_LOGMASS','HOSTGAL_LOGMASS_ERR','SIM_SALT2x0','SIM_SALT2x1','SIM_SALT2c','SIM_AV','SIM_RV','SIM_LIBID'];checks=[]
for arm in arms[1:]:
 common=sorted(set(h[arm])&set(h['P21']));diff={c:int(sum(h[arm][i][c]!=h['P21'][i][c] for i in common)) for c in cols};assert not any(diff.values());checks.append({'arm':arm,'common_written':len(common),'differences':diff})
f=out/'header-identity.json';f.write_text(json.dumps({'checks':checks,'source_rng_fact':'snlc_sim.c254 fills fixed random lists before each generated attempt; sntools.c6307-6345 fills exactly3*MXSTORE_RAN values independent of later rejection; GENFLUX_DRIVER draws epochnoise before any rho-dependent mixing. Raw variate arrays absent.','rng_source_sha256':hashlib.sha256((r/'phase2/official/build/SNANA-2fe0f56/src/sntools.c').read_bytes()).hexdigest(),'check_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n');print('Header metadata identical for',sum(q['common_written'] for q in checks),'paired object cases')
