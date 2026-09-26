"""Six-object M20 conditional forward/information gate; no observed outcome fitting.

Run with owned .venv, CPU affinity <=2. Source mutation occurs only in an
in-memory subclass to refine integration resolution; pinned files stay unchanged.
"""
import os
os.sched_setaffinity(0, sorted(os.sched_getaffinity(0))[:2])
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('MKL_NUM_THREADS','1')
os.environ.setdefault('XLA_FLAGS','--xla_cpu_multi_thread_eigen=false')
from pathlib import Path
import sys,csv,json,time,inspect,textwrap,hashlib
import numpy as np
from astropy.cosmology import FlatLambdaCDM
OUT=Path(__file__).resolve().parent; ROOT=OUT.parents[2]
sys.path.insert(0,str(OUT/'official-code'))
from calibration_bridge import MAP,raw_phot
import bayesn.bayesn_model as module
from bayesn import SEDmodel
import jax
import jax.numpy as jnp


def save(name,r): (OUT/name).write_text(json.dumps(r,indent=2)+'\n')


def main():
    start=time.monotonic()
    protocol={'parent_protocol_sha256':hashlib.sha256((OUT/'protocol.json').read_bytes()).hexdigest(),
              'scope':'Conditional forward operator and finite-derivative Fisher only; raw flux values never enter objective or derivatives.',
              'resolution_sequence':[300,600,1200,2400],
              'fiducial':{'AV':.3,'RV':3.1,'theta':0.,'epsilon_white':'42zero values',
                          'D':'flatLCDM(73.24,.28) at frozen zHD, used only for signal scale in design'},
              'derivative_steps':{'D':1e-4,'AV':1e-4,'RV':1e-3,'theta':1e-4,'eps_white':1e-4},
              'noise':'Quoted FLUXCALERR squared; measurement covariance diagonal; no SALT modelC or SN outcome residual.',
              'masks':'fixed header PEAKMJD, strict -10<phase<40; excludeu; full1%-peak throughput interval must lie in3000..18500A rest range; no S/N cut',
              'gates':{'successive_max_error_over_quoted_sigma':.1,'relative_derivative_halfstep_error':1e-4,
                       'gray_gauge_relative_flux_error':1e-12,'AV0_RV_derivative_max':1e-12},
              'Fisher':'dustAV,RV with nuisanceD,theta,42eps; compare originalDprior vs noDprior, same unit normal theta/eps prior; separately profile all nuisances without intrinsic priors',
              'physical_support':'OriginalRV>=.5 kept in parent model; law negative-extinction support independently recorded in protocol-addendum.json; no clamping.',
              'resource':'2CPU affinity; no posterior sampling; stop after finite sixobject gate.'}
    if (OUT/'forward-protocol.json').exists():
        assert json.loads((OUT/'forward-protocol.json').read_text())==protocol
    else:save('forward-protocol.json',protocol)
    pilots=list(csv.DictReader((OUT/'pilot-membership.csv').open()))
    bridge=list(csv.DictReader((OUT/'calibration-bridge.csv').open()))
    bridge={(r['survey'],r['raw_band']):r for r in bridge}
    objects=[]; rowledger=[]
    for p in pilots:
        head,rows=raw_phot(ROOT/p['photometry_path']); z=float(p['zHEL']); peak=float(head['PEAKMJD'][0])
        keep=[]
        for i,r in enumerate(rows):
            b=bridge[(p['survey'],r['band'])]; curve=np.loadtxt(OUT/'release-filters'/f"{b['custom_filter']}.dat")
            active=curve[:,1]>.01*curve[:,1].max(); lo,hi=curve[active,0][[0,-1]]/(1+z)
            phase=(r['MJD']-peak)/(1+z)
            reasons=[]
            if r['band']=='u':reasons.append('M20 u not validated')
            if not(-10<phase<40):reasons.append('phase outside open(-10,40)')
            if lo<3000 or hi>18500:reasons.append('1percent throughput support outside model')
            rowledger.append(dict(CID=p['CID'],row=i,band=r['band'],phase=phase,error=r['error'],
                                  rest_support_low=lo,rest_support_high=hi,keep=not reasons,reason=';'.join(reasons)))
            if not reasons:keep.append(dict(phase=phase,error=r['error'],band=b['custom_filter']))
        assert keep and all(r['error']>0 for r in keep)
        zerr=float(head.get('REDSHIFT_FINAL',head['REDSHIFT_HELIO'])[2])
        objects.append(dict(CID=p['CID'],z=z,zHD=float(p['zHD']),zerr=zerr,ebv=float(head['MWEBV'][0]),rows=keep))
    with (OUT/'forward-row-mask.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rowledger[0]));w.writeheader();w.writerows(rowledger)
    n=len(objects);nobs=max(len(o['rows']) for o in objects)
    t=np.zeros((nobs,n)); errors=np.ones_like(t); mask=np.zeros_like(t); bands=np.full((nobs,n),'NULL_BAND',object)
    for j,o in enumerate(objects):
        for i,r in enumerate(o['rows']):t[i,j]=r['phase'];errors[i,j]=r['error'];bands[i,j]=r['band'];mask[i,j]=1
    D=FlatLambdaCDM(H0=73.24,Om0=.28).distmod([o['zHD'] for o in objects]).value
    pars=np.zeros((n,46));pars[:,0]=D;pars[:,1]=.3;pars[:,2]=3.1
    # Sole implementation change: controlled quadrature refinement.
    source=textwrap.dedent(inspect.getsource(SEDmodel._setup_band_weights))
    assert source.count('self.spectrum_bins = 300')==1
    source=source.replace('self.spectrum_bins = 300','self.spectrum_bins = self.audit_spectrum_bins')
    scope={};exec(compile(source,'<audit_resolution_override>','exec'),module.__dict__,scope)
    class Refined(SEDmodel):_setup_band_weights=scope['_setup_band_weights']
    (OUT/'resolution-method.py').write_text(source)
    rows=[];prev=None;last=None
    for bins in protocol['resolution_sequence']:
        Refined.audit_spectrum_bins=bins
        tic=time.monotonic();model=Refined(load_model='M20_model',num_devices=1,filter_yaml=str(OUT/'release-filter-config.yaml'))
        weights=model._calculate_band_weights(jnp.array([o['z'] for o in objects]),jnp.array([o['ebv'] for o in objects]))
        bi=jnp.array([[model.band_dict[b] for b in row] for row in bands])
        hs=jnp.array([19+np.floor(t),19+np.ceil(t),np.remainder(t,1)])
        jt=model.J_t_map(t.flatten(order='F'),model.tau_knots,model.KD_t).reshape((nobs,n,6),order='F').transpose(1,2,0)
        @jax.jit
        def forward(p):
            eps=(model.L_Sigma@p[:,4:].T).T.reshape((n,7,6),order='F')
            eps=jnp.zeros((n,9,6)).at[:,1:-1,:].set(eps)
            return model.get_flux_batch(model.M0,p[:,3],p[:,1],model.W0,model.W1,eps,p[:,0],p[:,2],bi,jnp.array(mask),jt,hs,weights)
        pred=np.asarray(forward(jnp.array(pars))); assert np.isfinite(pred).all()
        diff=np.zeros(n) if prev is None else np.max(np.abs(pred-prev)/errors*mask,axis=0)
        rows.append(dict(bins=bins,seconds=time.monotonic()-tic,max_delta_over_error=float(diff.max()),per_object_delta_over_error=list(diff)))
        np.savez_compressed(OUT/f'forward-{bins}.npz',prediction=pred,errors=errors,mask=mask,parameters=pars,
                            phase=t,CID=np.array([o['CID'] for o in objects]))
        print(json.dumps(rows[-1]),flush=True)
        prev=pred;last=(model,forward,pred)
    save('forward-resolution.json',rows)
    model,forward,pred=last
    # Derivatives are simultaneous same-coordinate shifts for the6independent objects.
    steps=np.r_[1e-4,1e-4,1e-3,np.repeat(1e-4,43)]
    def jac(steps):
        cols=[]
        for k,h in enumerate(steps):
            shift=np.zeros_like(pars);shift[:,k]=h
            cols.append((np.asarray(forward(jnp.array(pars+shift)))-np.asarray(forward(jnp.array(pars-shift))))/(2*h))
        return np.stack(cols,axis=-1)
    tic=time.monotonic();J=jac(steps);Jh=jac(steps/2)
    scaled_error=float(np.linalg.norm((J-Jh)/errors[:,:,None])/max(1.,np.linalg.norm(Jh/errors[:,:,None])))
    zpars=pars.copy();zpars[:,1]=0
    zhi=zpars.copy();zlo=zpars.copy();zhi[:,2]+=.1;zlo[:,2]-=.1
    av0deriv=(np.asarray(forward(jnp.array(zhi)))-np.asarray(forward(jnp.array(zlo))))/.2
    gray=pars.copy();gray[:,0]+=.15
    graypred=np.asarray(forward(jnp.array(gray)))
    expected=pred*10**(-.4*.15)
    gauge=float(np.max(np.abs(graypred-expected)/np.maximum(abs(expected),1e-200)*mask))
    results=[]
    for i,o in enumerate(objects):
        use=mask[:,i]>0; A=Jh[use,i,:]/errors[use,i,None]
        signal=pred[use,i]/errors[use,i]
        dust=A[:,1:3]; nuisance=A[:,np.r_[0,3:46]]
        U,s,V=np.linalg.svd(nuisance,full_matrices=False);tol=max(nuisance.shape)*np.finfo(float).eps*s[0]
        U=U[:,s>tol];P=dust-U@(U.T@dust)
        free=P.T@P
        r=dict(CID=o['CID'],N_epochs=int(use.sum()),D_fiducial=float(D[i]),
               unregularized_nuisance_rank=U.shape[1],unregularized_dust_information=free.tolist(),
               unregularized_dust_eigenvalues=np.linalg.eigvalsh(free).tolist())
        sigmaext=5/(o['zHD']*np.log(10))*np.sqrt(o['zerr']**2+(150/299792.458)**2)
        for label,prec in [('no_distance',0.),('external_distance',1/(sigmaext**2+.088**2))]:
            B=nuisance.T@nuisance+np.diag(np.r_[prec,np.ones(43)])
            eff=dust.T@dust-dust.T@nuisance@np.linalg.solve(B,nuisance.T@dust)
            cov=np.linalg.inv(eff)
            r[label]=dict(information=eff.tolist(),eigenvalues=np.linalg.eigvalsh(eff).tolist(),
                          sigma_AV_mag=float(np.sqrt(cov[0,0])),sigma_RV=float(np.sqrt(cov[1,1])),
                          correlation=float(cov[0,1]/np.sqrt(cov[0,0]*cov[1,1])))
        r['sigma_ext_mag']=float(sigmaext)
        # Noiseless recovery of amplitude only; f is arbitrary reference, y is generated.
        amp=signal@(signal*10**(-.4*.15))/(signal@signal)
        r['noiseless_amplitude_recovery_error']=float(amp-10**(-.4*.15))
        results.append(r)
    np.savez_compressed(OUT/'forward-geometry-arrays.npz',prediction=pred,J=J,Jhalf=Jh,parameters=pars,
                        errors=errors,mask=mask,AV0_RV_derivative=av0deriv,CID=np.array([o['CID'] for o in objects]))
    save('forward-geometry-result.json',dict(status='No measured-flux likelihood or dust fit; all Fisher values at prescribed latent coordinates.',
        objects=results,halfstep_relative_derivative_error=scaled_error,AV0_RV_derivative_max=float(abs(av0deriv).max()),
        gray_scaling_relative_error=gauge,derivative_seconds=time.monotonic()-tic,total_seconds=time.monotonic()-start,
        final_resolution_gate_pass=rows[-1]['max_delta_over_error']<.1,
        gates_pass=(rows[-1]['max_delta_over_error']<.1 and scaled_error<1e-4 and abs(av0deriv).max()<1e-12 and gauge<1e-12)))
    print(json.dumps({'geometry_seconds':time.monotonic()-tic,'total_seconds':time.monotonic()-start,
                      'halfstep_error':scaled_error,'AV0_null':float(abs(av0deriv).max()),'gray_error':gauge}),flush=True)


if __name__=='__main__':main()
