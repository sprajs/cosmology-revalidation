"""Targeted holdout mean-edge tangent sensitivity, preserving primary residual outputs."""
from pathlib import Path
import os,sys,re,json,hashlib,subprocess
from concurrent.futures import ThreadPoolExecutor
import numpy as np,pandas as pd
from scipy.linalg import solve_triangular
ROOT=Path(__file__).resolve().parents[3];P=Path(__file__).resolve().parent;factor=.5 if '--half' in sys.argv else 1.;O=P/('hold_jh' if factor==.5 else 'hold_jm');E=ROOT/'runs/research_2026_09_26/simulation_holdout_control'
EXE=ROOT/'phase2/official/build/SNANA-audit-v3/bin/snlc_fit.exe';PY=ROOT/'phase2/env-official/bin/python';EXPORT=ROOT/'scripts/research_2026_09_26/export_flux_objectives.py'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
# Fixed source-identified mean-edge objects, all4 across both laws; no residual score selection.
selection={'P21':[19675,14533],'G10':[980,6801]}
assert json.loads((E/'native-edge-selection.json').read_text())
selection_hash=sha(E/'native-edge-selection.json')
steps={1:.001*factor,2:.0001*factor,3:.01*factor}
if not O.exists():
 O.mkdir();protocol={'selection':selection,'rule':'Maximum mean mismatch then maximum nonpositive-count remaining object perarm, fixed by saved implementation ledgers; scorer computation already underway but no aggregate inspected','steps_x1_c_t0':steps,'laws':[-99],'native_coordinates':'Nominal free-fit final x0,x1,c,t0 plus/minus only one declared increment; LFIXPAR_ALL=T for derivative measurement','fixed_frozen_inputs':'Same version,KCOR,model,ordinary cuts; require exact accepted band/time/data/error identity with free baseline','covariance':'Nominal free-fit C retained for tangent and projection comparison; varied C not used','comparison':'Python modern F99 tangent as actual Sol scoring approximation, law-matched Python diagnostic, exact native finite derivative; native amplitude exact','scope':'4 implementation-selected holdout objects; other residual scores unchanged; not a global derivative-error bound or selection null','selection_sha256':selection_hash,'mask_amendment':'Freeze original accepted rawPHOT rows via rejection bit29; no reclipping or phasecuts in derivative-only exports. Original freefit/C/scores unchanged. Prior two attempts preserved.','boundary_derivative':'CID14533 x1 inward difference because native hardbound−5; repeat with halfstep. Othercolumns central.','source_sha256':sha(__file__),'binary_sha256':sha(EXE)};(O/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
else:raise RuntimeError('preserve existing derivatives')
from astropy.io import fits
from collections import defaultdict,deque
maskled=[]
for arm,ids in selection.items():
 version='HE_'+arm;datadir=O/'inputs'/version;datadir.mkdir(parents=True);orig=ROOT/'phase2/literature/simulations/outputs'/('PH2_pilot02_'+arm);stem='PH2_pilot02_'+arm
 with fits.open(orig/(stem+'_HEAD.FITS')) as h, fits.open(orig/(stem+'_PHOT.FITS')) as ph:
  photo=ph[1].data;flags=photo['PHOTFLAG'].copy();photo['PHOTFLAG'][:]=np.bitwise_or(flags,1<<29)
  for cid in ids:
   head=[x for x in h[1].data if int(x['SNID'])==cid][0];lo=int(head['PTROBS_MIN'])-1;hi=int(head['PTROBS_MAX']);index=defaultdict(deque)
   for rawidx in range(lo,hi):
    v=photo[rawidx];key=(float(v['MJD']),str(v['BAND']).strip(),float(v['FLUXCAL']),float(v['FLUXCALERR']));index[key].append(rawidx)
   a=np.load(E/arm/f'objectives/objective_{cid}.npz');selected=[]
   for t,b,f,e in zip(a['MJD'],a['band'],a['data_flux'],a['data_fluxerr']):
    key=(float(t),str(b),float(f),float(e));assert len(index[key]);rawidx=index[key].popleft();photo['PHOTFLAG'][rawidx]=int(flags[rawidx]) & ~(1<<29);selected.append(rawidx)
   maskled.append(dict(arm=arm,CID=cid,accepted_raw_zero_based=selected))
  h.writeto(datadir/(version+'_HEAD.FITS'));ph.writeto(datadir/(version+'_PHOT.FITS'))
 (datadir/(version+'.LIST')).write_text(version+'_HEAD.FITS\n');(datadir/(version+'.README')).write_text((orig/(stem+'.README')).read_text())
(O/'mask-ledger.json').write_text(json.dumps(maskled,indent=2)+'\n')
def freeze_nml(nml,arm):
 nml=re.sub(r"(?m)^\s*PRIVATE_DATA_PATH\s*=.*$",f"    PRIVATE_DATA_PATH = '{O/'inputs'}'",nml)
 nml=re.sub(r"(?m)^\s*VERSION_PHOTOMETRY\s*=.*$",f"    VERSION_PHOTOMETRY = 'HE_{arm}'",nml)
 nml=re.sub(r'(?m)^\s*OPT_SNCID_LIST\s*=.*$', '    OPT_SNCID_LIST = 3',nml)
 for key,value in [('CUTWIN_TREST','-999.,999.'),('FITWIN_TREST','-999.,999.'),('DELCHI2_REJECT','1.0E9')]:nml=re.sub(r'(?m)^\s*'+key+r'\s*=.*$',f'    {key} = {value}',nml)
 return nml.replace('&SNLCINP','&SNLCINP\n    PHOTFLAG_MSKREJ = 536870912',1)
jobs=[]
for arm,ids in selection.items():
 for law in ['approx_minus99']:
  baseline=E/arm
  for j,h in steps.items():
   for sign in [-1,1]:
    d=O/arm/law/f'p{j}_{"plus" if sign>0 else "minus"}';d.mkdir(parents=True);seed=d/'seed.FITRES';lines=['VARNAMES: CID PKMJD x0 x1 c']
    for cid in ids:
     pars=np.load(baseline/f'objectives/objective_{cid}.npz')['parameters_x0_x1_c_t0'].copy();pars[j]+=0. if (cid==14533 and j==1 and sign==-1) else sign*h;x0,x1,c,t0=pars;lines.append(f'SN: {cid} {t0:.17g} {x0:.17g} {x1:.17g} {c:.17g}')
    seed.write_text('\n'.join(lines)+'\n');nml=freeze_nml((baseline/'fit.nml').read_text(),arm);nml=re.sub(r"(?m)^\s*SNCID_LIST_FILE\s*=.*$",f"    SNCID_LIST_FILE = '{seed}'",nml);nml=re.sub(r"(?m)^\s*TEXTFILE_PREFIX\s*=.*$",f"    TEXTFILE_PREFIX = '{d/'fit'}'",nml);nml=nml.replace('&FITINP','&FITINP\n    LFIXPAR_ALL = T',1);(d/'fit.nml').write_text(nml);(d/'cids.txt').write_text('\n'.join(map(str,ids))+'\n');jobs.append(d)
# Independent zero-increment fixed-coordinate closure before derivatives.
for arm,ids in selection.items():
 d=O/arm/'baseline';d.mkdir(parents=True);seed=d/'seed.FITRES';lines=['VARNAMES: CID PKMJD x0 x1 c']
 for cid in ids:
  pars=np.load(E/arm/f'objectives/objective_{cid}.npz')['parameters_x0_x1_c_t0'];x0,x1,c,t0=pars;lines.append(f'SN: {cid} {t0:.17g} {x0:.17g} {x1:.17g} {c:.17g}')
 seed.write_text('\n'.join(lines)+'\n');nml=freeze_nml((E/arm/'fit.nml').read_text(),arm);nml=re.sub(r"(?m)^\s*SNCID_LIST_FILE\s*=.*$",f"    SNCID_LIST_FILE = '{seed}'",nml);nml=re.sub(r"(?m)^\s*TEXTFILE_PREFIX\s*=.*$",f"    TEXTFILE_PREFIX = '{d/'fit'}'",nml);nml=nml.replace('&FITINP','&FITINP\n    LFIXPAR_ALL = T',1);(d/'fit.nml').write_text(nml);(d/'cids.txt').write_text('\n'.join(map(str,ids))+'\n');jobs.append(d)
env=os.environ.copy();env.update(SNANA_DIR=str(ROOT/'phase2/official/build/SNANA-audit-v3'),SNDATA_ROOT=str(ROOT/'phase2/official/inputs/SNDATA_ROOT'),LD_LIBRARY_PATH=str(ROOT/'phase2/official/build/sysroot/usr/lib'),OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
def run(d):
 with (d/'fit.log').open('x') as f:r=subprocess.run([str(EXE),str(d/'fit.nml')],cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT)
 assert r.returncode==0,d
 r=subprocess.run([str(PY),str(EXPORT),'--log',str(d/'fit.log'),'--output',str(d/'objectives'),'--expected-cids',str(d/'cids.txt')],cwd=ROOT,env=env,capture_output=True,text=True);(d/'export.log').write_text(r.stdout+r.stderr);assert r.returncode==0,(d,r.stderr)
 return str(d.relative_to(O))
with ThreadPoolExecutor(max_workers=2) as pool:
 for d in pool.map(run,jobs):print('finished',d,flush=True)
# Comparison is separate so a failed mask gate leaves source artifacts intact.
sys.path.insert(0,str(ROOT));from scripts.salt_dust_audit import flux_response as fr
from scripts.salt_dust_audit.snana_extinction import SnanaExtinction
ext=SnanaExtinction(ROOT/'runs/salt_dust_audit/snana_extinction/libsnana_extinction.so')
class OldLaw(fr.VariableF99):
 def propagate(self,wave,flux,phase=None):
  ebv,rv=self._parameters;return flux*10**(-.4*ext(wave,rv,ebv,option=-99))
model,bands,paths,zp=fr.build_model();old=fr.sncosmo.Model(source=model.source,effects=[OldLaw(),fr.VariableF99()],effect_names=['mw','host'],effect_frames=['obs','rest']);offsets={str(r['Filter Name'])[-1]:float(r['Primary Mag']) for r in zp}
basis=np.load(P/'validation1020/frozen-discovery-coefficients.npz');coeff=basis['gauge_griz']@basis['basis_mean'];rows=[];arrays={}
for arm,ids in selection.items():
 for law in ['approx_minus99']:
  for cid in ids:
   a=np.load(E/arm/f'objectives/objective_{cid}.npz');baseline0=np.load(O/arm/'baseline'/f'objectives/objective_{cid}.npz');
   for key in ['MJD','band','data_flux','data_fluxerr','model_flux','parameters_x0_x1_c_t0']:assert np.array_equal(a[key],baseline0[key]),(arm,cid,key)
   J=np.zeros((len(a['MJD']),4));J[:,0]=-fr.K*a['model_flux'];maxpar=0.
   for j,h in steps.items():
    pp=[]
    for sign in [-1,1]:
     b=np.load(O/arm/law/f'p{j}_{"plus" if sign>0 else "minus"}'/f'objectives/objective_{cid}.npz')
     for key in ['MJD','band','data_flux','data_fluxerr']:assert np.array_equal(a[key],b[key]),(arm,law,cid,j,sign,key)
     expect=a['parameters_x0_x1_c_t0'].copy();expect[j]+=0. if (cid==14533 and j==1 and sign==-1) else sign*h;err=float(np.max(abs(expect-b['parameters_x0_x1_c_t0'])));maxpar=max(maxpar,err);assert err<1e-12
     pp.append(b['model_flux'])
    J[:,j]=(pp[1]-pp[0])/(h if (cid==14533 and j==1) else 2*h)
   C=a['frozen_flux_covariance'];L=np.linalg.cholesky(C);jw=solve_triangular(L,J,lower=True);U,s,V=np.linalg.svd(jw,full_matrices=False);assert (s>s[0]*1e-10).sum()==4
   raw=solve_triangular(L,-fr.K*a['model_flux']*np.array([coeff['griz'.index(b)] for b in a['band']]),lower=True);v=raw-U@(U.T@raw);rw=solve_triangular(L,a['data_flux']-a['model_flux'],lower=True)
   native_a=float(v@rw);native_I=float(v@v);x0,x1,c,t0=a['parameters_x0_x1_c_t0'];bp=np.array([bands[b] for b in a['band']],dtype=object);conv=np.array([10**(-.4*(.27+offsets[b])) for b in a['band']]);z=float(a['zHEL'][0]);ebv=float(a['MWEBV'][0])
   for label,m in [('actual_Python_exact99',model),('law_matched_Python',old if law=='approx_minus99' else model)]:
    def f(th):
     m.set(z=z,t0=t0+th[3],x0=x0*np.exp(-fr.K*th[0]),x1=x1+th[1],c=c+th[2],mwebv=ebv,mwrv=3.1,hostebv=0,hostrv=3.1);return conv*m.bandflux(bp,a['MJD'],zp=27.5,zpsys='ab')
    jj=J.copy()
    for j,h in {1:.001,2:.0001,3:.01}.items():
     e=np.eye(4)[j]*h;jj[:,j]=(f(e)-f(-e))/(2*h)
    w=solve_triangular(L,jj,lower=True);uu,ss,_=np.linalg.svd(w,full_matrices=False);vp=raw-uu@(uu.T@raw)
    row=dict(arm=arm,law=law,CID=cid,zHEL=z,epochs=len(v),comparison=label,parameter_max_error=maxpar,relative_whitened_derivative_norm=float(np.linalg.norm(w-jw)/np.linalg.norm(jw)),projector_spectral_difference=float(np.linalg.norm(U@U.T-uu@uu.T,2)),fixed_vector_difference_norm=float(np.linalg.norm(vp-v)),fixed_vector_relative_difference=float(np.linalg.norm(vp-v)/np.linalg.norm(v)),native_information=native_I,native_matched=native_a,Python_minus_native_information=float(vp@vp-native_I),Python_minus_native_matched=float(vp@rw-native_a));row['Python_minus_native_gain']=row['Python_minus_native_matched']-.5*row['Python_minus_native_information'];rows.append(row)
   arrays[f'{arm}_{law}_{cid}_J_native']=J
pd.DataFrame(rows).to_csv(O/'comparison.csv',index=False);np.savez_compressed(O/'native-jacobians.npz',**arrays)
result={'status':'PASS: exact same accepted rows and coordinates for all native finite differences','checks':rows,'source_sha256':sha(__file__),'protocol_sha256':sha(O/'protocol.json'),'scope':'Finite selected-coordinate tangent control, not full sample error bound'};(O/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(pd.DataFrame(rows).to_string(index=False),flush=True)
