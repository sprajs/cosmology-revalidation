"""Post-run verification of saved calibration artifacts; does not rewrite them."""
from datetime import datetime,timezone
from pathlib import Path
import hashlib,json,sys
import numpy as np,pandas as pd
from scipy.linalg import solve_triangular
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'runs/research_2026_09_26/deps'))
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy_healpix import HEALPix

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'runs/research_2026_09_26'
DEST=BASE/'calibration_verification'
CASES={
 'foreground_response/nonlinear_check':'foreground_refit_check.py',
 'calibration_stars':'calibration_star_extinction.py',
 'calibration_stars/sn_transfer':'stellar_correction_transfer.py',
}
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
checks={};snapshots={}
for name,source in CASES.items():
 path=BASE/name;manifest=json.loads((path/'manifest.json').read_text())
 source_path=ROOT/'scripts/research_2026_09_26'/source
 assert sha(source_path)==manifest['source_sha256']==sha(path/'executed_source.py')
 inputs=manifest.get('inputs_sha256',manifest.get('input_sha256',{}))
 for rel,digest in inputs.items():assert sha(ROOT/rel)==digest,(name,rel)
 snapshots[name]={'original_manifest_sha256':sha(path/'manifest.json'),
                  'source_sha256':sha(source_path),
                  'saved_files_sha256':{str(p.relative_to(ROOT)):sha(p) for p in sorted(path.iterdir()) if p.is_file()},
                  'input_hashes_verified':len(inputs)}

# Independently check local map metadata, bit meanings and the 64 saved samples.
readme=(ROOT/'data/dust/csfd-v2/readme.txt').read_text()
assert "Bit  0" not in readme # readme describes bits in a simple table below
assert "0    'LSS_corr'" in readme and "2    'cosmology'" in readme
mapdir=ROOT/'data/dust/csfd-v2'
for name in ('sfd_ebv','csfd_ebv','mask'):
 with fits.open(mapdir/(name+'.fits'),memmap=True) as hd:
  assert int(hd[1].header['NSIDE'])==2048 and hd[1].header['ORDERING'].strip()=='RING'
maprows=pd.read_csv(BASE/'foreground_response/map_samples.csv',dtype={'CID':str}).set_index('CID')
sky=SkyCoord(maprows.RA.to_numpy()*u.deg,maprows.DEC.to_numpy()*u.deg,frame='icrs').galactic
hp=HEALPix(nside=2048,order='ring',frame='galactic');pixels,weights=hp.bilinear_interpolation_weights(sky.l,sky.b);nearest=hp.lonlat_to_healpix(sky.l,sky.b)
with fits.open(mapdir/'mask.fits',memmap=True) as hd:bits=hd[1].data.field(0).reshape(-1).astype(int)
reliable=np.all((bits[pixels]&4)!=0,axis=0);footprint=np.all((bits[pixels]&1)!=0,axis=0)
assert np.array_equal(reliable,maprows.reliable.to_numpy())
assert np.array_equal(footprint,maprows.correction_footprint.to_numpy())
assert np.array_equal(bits[nearest],maprows.mask_nearest.to_numpy())
headpath=ROOT/'sources/repos/des-science__DES-SN5YR@1.3/0_DATA/DES-SN5YR_DES/DES-SN5YR_DES_HEAD.FITS.gz'
with fits.open(headpath) as hd:
 heads=hd[1].data;lookup={str(s).strip():i for i,s in enumerate(heads['SNID'])}
 raw=np.array([heads['MWEBV'][lookup[c]] for c in maprows.index]);head_dtype=str(heads['MWEBV'].dtype)
head_diff=float(np.max(abs(raw-maprows.HEAD_MWEBV.to_numpy())))
assert head_diff<1e-8
checks['map_and_head']={'map_nside':2048,'map_order':'RING','reliable_objects':int(reliable.sum()),'correction_footprint_objects':int(footprint.sum()),'HEAD_MWEBV_dtype':head_dtype,'HEAD_to_saved_csv_max_abs_difference':head_diff,'direct_replacement_identity_max_abs_error':float(np.max(abs(maprows.delta_ebv_direct_replacement-(maprows.delta_ebv_cleaning+maprows.delta_ebv_resampling))))}

# Check the star source join, distinct sky positions, fixed selection and gray mode.
base=ROOT/'sources/repos/bap37__Dovekie';obs=pd.read_csv(base/'output_observed_apermags+AV/DES_observed.csv')
sd=pd.read_csv(base/'output_synthetic_magsaper/synth_DES_shift_0.000.txt',sep=r'\s+')
sp=pd.read_csv(base/'output_synthetic_magsaper/synth_PS1_shift_0.000.txt',sep=r'\s+')
joined=sd.merge(sp,on=['standard','standard_catagory'],validate='one_to_one')
assert len(joined)==len(sd)==len(sp)==451 and not obs[['RA','DEC']].duplicated().any()
selection=pd.read_csv(BASE/'calibration_stars/star_selection.csv')
assert len(selection)==len(obs)==1255 and int(selection.selected_both.sum())==1178
stars=pd.read_csv(BASE/'calibration_stars/extinction_branch_changes.csv');stars=stars[stars.arm=='fixed_intersection']
delta={lib:part.set_index('band').loc[list('griz'),'delta_mean_correction_stored_minus_zero'].to_numpy() for lib,part in stars.groupby('library')}
checks['star_join_and_cuts']={'observed_rows':len(obs),'unique_observed_coordinates':True,'synthetic_DES_rows':len(sd),'synthetic_PS1_rows':len(sp),'synthetic_joined_rows':len(joined),'fixed_intersection_rows':int(selection.selected_both.sum()),'heldout_sky_cells':int(selection.loc[selection.selected_both & selection.test,'cell'].nunique()),'fixed_vectors_gray_mmag':{lib:float(values.mean()*1000) for lib,values in delta.items()}}

# Recompute all stored foreground and stellar projected shifts/gains with an
# SVD orthogonal projector independent of the generation code's QR projector.
a=np.load(ROOT/'runs/salt_dust_audit/flux_response/matrices.npz')
fr=pd.read_csv(BASE/'foreground_response/object_responses.csv',dtype={'CID':str})
sr=pd.read_csv(BASE/'calibration_stars/sn_transfer/object_responses.csv',dtype={'CID':str})
D=np.array([1,.16087,-3.1178,0.]);maxerr={k:0. for k in ('foreground_shift','foreground_gain','stellar_shift','stellar_gain','gray_response')}
for cid,row in maprows.iterrows():
 pre=cid+'__';J=a[pre+'jacobian_flux'];G=a[pre+'nuisance_flux'];res=a[pre+'flux_observed']-a[pre+'flux_model']
 for weighting in ('measurement_only','measurement_plus_model'):
  L=np.linalg.cholesky(a[pre+weighting+'_covariance']);jw=solve_triangular(L,J,lower=True);Q,_,_=np.linalg.svd(jw,full_matrices=False)
  rp=solve_triangular(L,res,lower=True);rp=rp-Q@(Q.T@rp)
  response=a[pre+weighting+'_response']
  maxerr['gray_response']=max(maxerr['gray_response'],float(np.max(abs(response[:,5:].sum(axis=1)-[1,0,0,0]))))
  for mode in ('cleaning','resampling','direct_replacement'):
   e=row['delta_ebv_'+mode];v=solve_triangular(L,G[:,4]*e,lower=True);vp=v-Q@(Q.T@v)
   saved=fr[(fr.CID==cid)&(fr.weighting==weighting)&(fr['mode']==mode)].iloc[0]
   shift=-float(D@response[:,4])*e;gain=float(rp@vp-.5*(vp@vp))
   maxerr['foreground_shift']=max(maxerr['foreground_shift'],abs(saved.pre_BBC_delta_standardized_mag-shift))
   maxerr['foreground_gain']=max(maxerr['foreground_gain'],abs(saved.fixed_map_loglike_gain-gain))
  for lib,vector in delta.items():
   v=solve_triangular(L,G[:,5:]@vector,lower=True);vp=v-Q@(Q.T@v)
   saved=sr[(sr.CID==cid)&(sr.weighting==weighting)&(sr.library==lib)].iloc[0]
   shift=float(D@response[:,5:]@vector);gain=float(-rp@vp-.5*(vp@vp))
   maxerr['stellar_shift']=max(maxerr['stellar_shift'],abs(saved.delta_standardized_preBBC-shift))
   maxerr['stellar_gain']=max(maxerr['stellar_gain'],abs(saved.residual_loglike_gain_if_applied-gain))
assert len(fr)==384 and len(sr)==256 and max(maxerr.values())<1e-7
checks['independent_svd_reprojection']={'foreground_rows':len(fr),'stellar_rows':len(sr),'max_abs_errors':maxerr,'model_foreground_sign':'fixed observed flux, model derivative G: fitted parameter shift -R delta; gain r_perp dot G_perp delta - norm^2/2','data_calibration_sign':'positive observed magnitude delta, data derivative G: fitted parameter shift +R delta; gain -r_perp dot G_perp delta - norm^2/2'}

nonlinear=pd.read_csv(BASE/'foreground_response/nonlinear_check/paired_refits.csv',dtype={'CID':str})
summary=json.loads((BASE/'foreground_response/nonlinear_check/summary.json').read_text())
assert len(nonlinear)==12 and nonlinear.CID.nunique()==6 and set(nonlinear.target)=={'noiseless_reference','actual_observed'}
assert nonlinear.baseline_success.all() and nonlinear.alternative_success.all() and summary['all_success']
checks['nonlinear_refit_gate']={'rows':len(nonlinear),'unique_objects':nonlinear.CID.nunique(),'all_optimizer_success':True,'max_absolute_linear_error_by_target':nonlinear.groupby('target').difference_mag.apply(lambda x:float(x.abs().max())).to_dict(),'max_baseline_optimality':float(nonlinear.baseline_optimality.max()),'max_alternative_optimality':float(nonlinear.alternative_optimality.max())}

DEST.mkdir(exist_ok=False)
manifest={'status':'post-run verification of current saved checkpoint; output hashes were not captured by the original three run manifests','verified_utc':datetime.now(timezone.utc).isoformat(),'verifier_source_sha256':sha(Path(__file__)),'runs':snapshots,'checks':checks}
(DEST/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps({'runs':list(snapshots),'checks':list(checks),'max_reprojection_error':max(maxerr.values())},indent=2))
