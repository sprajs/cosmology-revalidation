"""One-object native mean/C export gate, not a replacement scientific fitter."""
from pathlib import Path
import json,re,os,subprocess,time,argparse
import numpy as np
import paired_refit as p
from check_results import table
O=p.OUT;G=O/'native-profile-gate-v2';CID='DES16C1cim'
BUILD=p.ROOT/'phase2/official/build/SNANA-audit-v3';EXE=BUILD/'bin/snlc_fit.exe'

def export(d):
    flux={};inv={};ct=None;setup=None;covset={}
    s=(d/'fit.log').read_text();assert 'ENDING PROGRAM GRACEFULLY.' in s
    for line in s.splitlines():
      a=line.split()
      if not a:continue
      if a[0]=='PHASE2_FLUX:':
        assert a[1]==CID;flux[int(a[2])]=[a[3]]+[float(x) for x in a[4:]]
      elif a[0]=='PHASE2_COVINV:':inv[int(a[2])]=[float(x) for x in a[3:]]
      elif a[0]=='PHASE2_OBJECTIVE:':ct=(int(a[2]),np.array([float(x) for x in a[3:]]))
      elif a[0]=='PHASE2_SETUP:':setup=np.array([float(x) for x in a[2:]])
      elif a[0]=='PHASE2_COVSET:':covset[int(a[2])]=[float(x) for x in a[3:]]
    assert ct is not None
    n,v=ct;ff=[flux[i] for i in range(1,n+1)];a=np.array([r[1:] for r in ff]);W=np.array([inv[i] for i in range(1,n+1)])
    assert W.shape==(n,n)
    C=np.linalg.inv(W);res=a[:,4]-a[:,2];err=float(res@W@res+v[1]-v[0]);assert abs(err)<1e-7
    assert np.linalg.eigvalsh(C).min()>0
    # Generic native slots are DLMAG,STRETCH,AV,PKMJD for SNooPy. No SALT conversion.
    q=dict(parameters=np.array(v[2:6]),total_chi2=v[0],prior_chi2=v[1],data_chi2=v[0]-v[1],band=np.array([r[0] for r in ff]),
      MJD=a[:,0],rest_phase=a[:,1],model_flux=a[:,2],model_magerr=a[:,3],data_flux=a[:,4],data_fluxerr=a[:,5],zHEL=a[:,6],MWEBV=a[:,7],C=C,W=W,setup=setup,
      covset=np.array([covset[i] for i in range(1,n+1)]),objective_reconstruction_error=err)
    np.savez_compressed(d/'native.npz',**q)
    return q

def config(d,version='RSR_B',fixed=None):
    d.mkdir(parents=True,exist_ok=False)
    s=(O/'fits/full/fixed/B/fit.nml').read_text()
    s=re.sub(r'(?m)^\s*SNCCID_LIST\s*=.*$',f" SNCCID_LIST = '{CID}'",s)
    s=re.sub(r'(?m)^\s*TEXTFILE_PREFIX\s*=.*$',f" TEXTFILE_PREFIX = '{d}/fit'",s)
    s=re.sub(r'(?m)^\s*VERSION_PHOTOMETRY\s*=.*$',f" VERSION_PHOTOMETRY = '{version}'",s)
    if fixed is not None:
      s=s.replace('  &FITINP','  &FITINP\n LFIXPAR_ALL = T\n INISTP_SHAPE = 0\n INISTP_AV = 0\n INISTP_DLMAG = 0')
      for key,val in zip(['INIVAL_DLMAG','INIVAL_SHAPE','INIVAL_AV','INIVAL_PEAKMJD'],fixed):
        s=re.sub(r'(?m)^\s*'+key+r'\s*=.*$', '',s)
        s=s.replace('  &FITINP',f'  &FITINP\n {key} = {val:.17g}')
    (d/'fit.nml').write_text(s)
    return d/'fit.nml'

def run(d):
    env=os.environ.copy();env.update(SNANA_DIR=str(BUILD),SNDATA_ROOT=str(p.ROOT/'phase2/official/inputs/SNDATA_ROOT'),LD_LIBRARY_PATH=str(p.ROOT/'phase2/official/build/sysroot/usr/lib'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    assert not(d/'fit.log').exists();t=time.monotonic()
    with(d/'fit.log').open('x') as f:
      rr=subprocess.run([str(EXE),str(d/'fit.nml')],cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=180)
    p.save(d/'execution.json',dict(returncode=rr.returncode,seconds=time.monotonic()-t,log_sha256=p.sha(d/'fit.log'),nml_sha256=p.sha(d/'fit.nml')))
    assert rr.returncode==0
    print(json.dumps(dict(run=d.name,seconds=time.monotonic()-t)),flush=True)
    return export(d)

def main():
    G.mkdir(exist_ok=False)
    protocol=dict(status='Frozen after original engineering list-reader abort; parent authorized one object',CID=CID,arm='B',peak='fixed to released header',law='modern explicit99',
      reader_amendment='SNCID_LIST_FILE stores SALT names only; original abort retained. Native NML INIVAL keys are used instead, with explicit float32 rounding before seed checks. No model or observation changes.',
      reference='Rerun original start1 fixed-peak B with output-only audit binary; compare current binary parameters/mask.',
      subsequent='Hold exact reference coordinates except DLMAG offsets ±.15; exact fixed-coordinate nominal and time-sorted/reversed identical row variants.',
      objective='Reference C and accepted mask will be retained for later conditional geometry. Other runs C are diagnostic only, never substituted pointwise.',
      gates=dict(reference_parameter_absolute_tolerance=[1e-5,1e-5,1e-5,1e-4],objective_relative_tolerance=1e-6,multiplicative_model_relative_tolerance=2e-5,
      model_permutation_relative_tolerance=2e-5,raw_flux_error_identity='exact native exported values; match by band and exact MJD, preserve duplicates',mask='same multiset across all mean probes'),
      amplitude='Only after native multiplicativity gate: a=b/q, D=Dref−2.5log10(a), a>0; fixed C. No claim of a proper flat-distance marginal likelihood.',
      input_hashes={str(f.relative_to(p.ROOT)):p.sha(f) for f in [Path(__file__),EXE,BUILD/'src/snlc_fit.F90',O/'fits/full/fixed/B/fit.nml',O/'fits/full/fixed/B/fit.FITRES.TEXT',O/'data/RSR_B'/f'{CID}.snana.dat',O/'data/RSR_Bsort'/f'{CID}.snana.dat',O/'data/RSR_Breverse'/f'{CID}.snana.dat']})
    p.save(G/'protocol.json',protocol)
    config(G/'reference');ref=run(G/'reference')
    original=next(r for r in table(O/'fits/full/fixed/B/fit.FITRES.TEXT') if r['CID']==CID)
    rerun=table(G/'reference/fit.FITRES.TEXT')[0]
    changes={k:float(rerun[k])-float(original[k]) for k in ['DLMAG','STRETCH','AV','PKMJD','FITCHI2','NDOF']}
    assert np.allclose([changes[k] for k in ['DLMAG','STRETCH','AV','PKMJD']],0,rtol=0,atol=[1e-5,1e-5,1e-5,1e-4]),changes
    assert abs(changes['FITCHI2'])<1e-6*max(1,float(original['FITCHI2'])) and changes['NDOF']==0
    results={};nom=None
    for name,delta,version in [('fixed',0.,'RSR_B'),('distance_plus',.15,'RSR_B'),('distance_minus',-.15,'RSR_B'),('sorted',0.,'RSR_Bsort'),('reversed',0.,'RSR_Breverse')]:
      pars=ref['parameters'].copy();pars[0]+=delta;pars=pars.astype(np.float32).astype(float)
      config(G/name,version,pars);v=run(G/name)
      # Exact seed/parameter closure is measured separately from displayed FITRES precision.
      assert np.max(abs(v['parameters']-pars))<1e-9,(name,v['parameters'],pars)
      if name=='fixed':nom=v
      key=lambda q:np.lexsort((q['data_fluxerr'],q['data_flux'],q['MJD'],q['band']))
      i=key(nom);j=key(v)
      for k in ['band','MJD','data_flux','data_fluxerr']:
        assert np.array_equal(nom[k][i],v[k][j]),(name,k)
      actual_delta=float(v['parameters'][0]-nom['parameters'][0])
      predicted=nom['model_flux'][i]*10**(-.4*actual_delta)
      norm=max(np.max(abs(predicted)),1e-30);relative=float(np.max(abs(v['model_flux'][j]-predicted))/norm)
      assert relative<2e-5,(name,relative)
      results[name]=dict(requested_offset_DLMAG=delta,actual_offset_DLMAG=actual_delta,mean_scaling_max_relative=relative,max_parameter_seed_error=float(np.max(abs(v['parameters']-pars))),
        C_relative_to_reference=float(np.linalg.norm(v['C'][np.ix_(j,j)]-ref['C'][np.ix_(key(ref),key(ref))])/np.linalg.norm(ref['C'])),
        objective_reconstruction_error=float(v['objective_reconstruction_error']))
    y=ref['data_flux'];h=nom['model_flux'];W=ref['W'];b=float(h@W@y);q=float(h@W@h);a=b/q
    assert a>0
    D=float(nom['parameters'][0]-2.5*np.log10(a));r=y-a*h
    out=dict(status='One-object native mean/export and amplitude proof gate passed; no global minimization yet',CID=CID,epochs=len(y),reference_current_minus_audit=changes,
      reference_parameters_DLMAG_STRETCH_AV_PKMJD=ref['parameters'].tolist(),probes=results,
      analytic_positive_amplitude=dict(a=a,DLMAG=D,data_chi2=float(r@W@r),normal_equation_residual=float(h@W@r)),
      limitation='Uses one explicitly frozen native C and fixed peak/mask. All later C changes are ignored by design. No historical likelihood or sign-selection bias inference.')
    p.save(G/'result.json',out)
    p.save(G/'manifest.json',dict(protocol_sha256=p.sha(G/'protocol.json'),all_files_sha256={str(f.relative_to(G)):p.sha(f) for f in G.rglob('*') if f.is_file() and f.name!='manifest.json'}))
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
