"""Native spectral operators and conditional formed-mass age feasibility.

The noise ellipsoid and finite stellar/dust family are assumptions, not a
posterior or a measured supernova luminosity-evolution correction.
"""
import json
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.linalg import solve_triangular
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import lsq_linear
from scipy.stats import chi2
from astropy.cosmology import FlatLambdaCDM
import clarabel
import extinction

ROOT = Path(__file__).resolve().parents[4]
CODE = Path(__file__).parent
WORK = ROOT / '.work/unified-cosmology/calibrated-host-physics'
INPUT = ROOT / '.work/unified-cosmology/calibrated-hosts'
BANDS = [(3850.,3950.),(4000.,4100.),(4041.6,4079.75),(4083.5,4122.25),(4128.5,4161.)]
CAMERAS = [('B',-np.inf,5780.),('R',5780.,7570.),('Z',7570.,np.inf)]
C_KM_S = 299792.458
COSMO = FlatLambdaCDM(H0=70., Om0=.3)


def resolution(data):
    """Author FITS DIA convention, descending offsets; no normalization."""
    ndiag, n = data.shape
    assert ndiag % 2 == 1
    return sparse.dia_matrix((data, np.arange(ndiag//2, -(ndiag//2)-1, -1)), shape=(n,n)).tocsr()


def interpolation(wave, target):
    """Piecewise-linear density samples, with strict wavelength support."""
    assert target.min() > wave[0] and target.max() < wave[-1]
    right = np.searchsorted(wave,target)
    frac = (target-wave[right-1])/(wave[right]-wave[right-1])
    rows = np.repeat(np.arange(len(target)),2)
    return sparse.csr_matrix((np.column_stack([1-frac,frac]).ravel(),
        (rows,np.column_stack([right-1,right]).ravel())),shape=(len(target),len(wave)))


def velocity_smooth(wave, spectra, sigma):
    """Doppler convolution conserving integral Llambda dlambda locally."""
    if sigma == 0:
        return spectra
    lw = np.arange(np.log(3300.),np.log(4700.)+1e-4,1e-4)
    fine = np.exp(lw)
    local = (wave>3400)&(wave<4600)
    flat = spectra.reshape(-1,len(wave))
    answer = flat.copy()
    for j,row in enumerate(flat):
        luminosity_per_log_wave = np.interp(fine,wave,row)*fine
        convolved = gaussian_filter1d(luminosity_per_log_wave,sigma/C_KM_S/1e-4,mode='nearest')/fine
        answer[j,local] = np.interp(wave[local],fine,convolved)
    return answer.reshape(spectra.shape)


class Library:
    def __init__(self, refined=False):
        self.path = WORK/('ssp-library-refined.npz' if refined else 'ssp-library.npz')
        self.data = np.load(self.path)
        self.wave = self.data['wavelength_A']
        self.ages = self.data['ages_Gyr']
        self.metals = self.data['metallicity_logZ_solar']
        self.base = self.data['spectra_Lsun_per_A_per_formed_Msun']
        lines = self.data['emission_line_wavelength_A']
        self.lines = np.unique([x for x in lines if any(lo<=x<=hi for lo,hi in BANDS)])
        self.cache = {}

    def grid(self, z, scenario):
        velocity = {'velocity150':150.,'velocity300':300.}.get(scenario,0.)
        if velocity not in self.cache:
            self.cache[velocity] = velocity_smooth(self.wave,self.base,velocity)
        spectra = self.cache[velocity]
        ages = self.ages.copy()
        if scenario == 'clock_LCDM':
            limit = min(COSMO.age(z).value,ages[-1])
            keep = ages<=limit
            selected = spectra[:,keep]
            selected_ages = ages[keep]
            if selected_ages[-1] < limit:
                right = np.searchsorted(ages,limit)
                w = np.log(limit/ages[right-1])/np.log(ages[right]/ages[right-1])
                selected = np.concatenate([selected,((1-w)*spectra[:,right-1]+w*spectra[:,right])[:,None]],axis=1)
                selected_ages = np.r_[selected_ages,limit]
            ages,spectra = selected_ages,selected
        tau = [0.,.5,1.,2.] + ([4.] if scenario == 'dust_tau4' else [])
        rows,parameters = [],[]
        for iz,metal in enumerate(self.metals):
            for t in tau:
                for power in ([-.7] if t==0 else [-.7,-1.]):
                    rows.extend(spectra[iz]*np.exp(-t*(self.wave/5500.)**power))
                    parameters.extend([dict(age_Gyr=float(age),logZ=float(metal),tauV=t,power=power) for age in ages])
        return np.array(rows),parameters


def observed_operator(arrays,z,lines,mask_velocity):
    """Same signed means/covariance as data; emission masks alter both sides."""
    parts=[]
    for camera,lo,hi in CAMERAS:
        wave=arrays[camera+'_WAVELENGTH'].astype(float)
        flux=arrays[camera+'_FLUX'].astype(float)
        ivar=arrays[camera+'_IVAR'].astype(float)
        rest=wave/(1+z)
        valid=(wave>=lo)&(wave<hi)&(arrays[camera+'_MASK']==0)&np.isfinite(flux)&np.isfinite(ivar)&(ivar>0)
        if mask_velocity:
            for line in lines:
                valid &= np.abs(np.log(rest/line))*C_KM_S > mask_velocity
        widths=np.array([np.maximum(0.,np.minimum(rest+.4/(1+z),b)-np.maximum(rest-.4/(1+z),a))*valid for a,b in BANDS])
        parts.append(dict(camera=camera,wave=wave,flux=np.where(valid,flux,0.),
            variance=np.divide(1.,ivar,out=np.zeros_like(ivar),where=valid),widths=widths,valid=valid))
    retained=sum(p['widths'].sum(axis=1) for p in parts)
    assert np.all(retained>0), retained
    y=np.zeros(5);cov=np.zeros((5,5))
    for p in parts:
        weight=p['widths']/retained[:,None]
        weight[:2]*=p['wave'][None,:]**2/2.99792458e18*1e12
        p['weight']=weight
        y+=weight@p['flux'];cov+=(weight*p['variance'])@weight.T
    return y,cov,parts,retained


def problem(targetid,lib,scenario='primary'):
    meta=json.loads((INPUT/'spectra'/f'{targetid}.json').read_text())
    arrays=np.load(ROOT/meta['array_file'])
    z=meta['desi_z']
    mask_velocity=0 if scenario=='original_free_emission' else (800 if scenario=='mask800' else 400)
    y,cov,parts,widths=observed_operator(arrays,z,lib.lines,mask_velocity)
    contributing=[p['camera'] for p in parts if p['widths'].sum()>0]
    kernel=np.zeros((5,len(lib.wave)))
    gas=np.zeros((5,len(lib.lines)))
    resolution_rows=[]
    for p in parts:
        if not p['widths'].sum():continue
        cam,wave,weight=p['camera'],p['wave'],p['weight']
        native=resolution(arrays[cam+'_RESOLUTION'].astype(float))
        op=sparse.eye(len(wave),format='csr') if scenario=='identity_resolution' else native
        back=np.asarray(op.T@weight.T).T
        support=np.any(back!=0,axis=0)
        assert support.any()
        obs=wave[support];rest=obs/(1+z)
        foreground_scale=1. if scenario=='foreground1' else .86
        fore=10**(-.4*extinction.fitzpatrick99(obs,3.1*meta['fibermap']['EBV']*foreground_scale,3.1))
        factor=fore/(1+z)
        if scenario.startswith('tilt_'):
            factor*=np.exp((1 if scenario=='tilt_plus' else -1)*.02*(rest-4000)/100.)
        if scenario.startswith('camera_') and len(contributing)>1 and cam==contributing[-1]:
            factor*=1+(0.03 if scenario=='camera_plus' else -.03)
        kernel+=(back[:,support]*factor)@interpolation(lib.wave,rest)
        # A delta-line density has unit integral on the native observed grid;
        # its arbitrary nuisance amplitude is not counted as stellar mass.
        for j,line in enumerate(lib.lines):
            center=line*(1+z)
            if wave[0]<center<wave[-1]:
                right=np.searchsorted(wave,center)
                frac=(center-wave[right-1])/(wave[right]-wave[right-1])
                delta=np.zeros(len(wave));delta[right-1:right+1]=[(1-frac)/.8,frac/.8]
                gas[:,j]+=weight@(op@delta)
        relevant=np.any(weight!=0,axis=0)
        rowsum=np.asarray(native.sum(axis=1)).ravel()[relevant]
        resolution_rows.extend(rowsum.tolist())
    spectra,parameters=lib.grid(z,scenario)
    matrix=kernel@spectra.T
    # ONE common scale for every stellar column preserves formed-mass ratios.
    scale=float(np.median(matrix[2:]))
    assert scale>0
    matrix/=scale
    ages=np.array([p['age_Gyr'] for p in parameters]);mass=np.ones(len(ages))
    if scenario=='original_free_emission':
        keep=np.linalg.norm(gas,axis=0)>0
        matrix=np.column_stack([matrix,gas[:,keep]])
        ages=np.r_[ages,np.zeros(keep.sum())];mass=np.r_[mass,np.zeros(keep.sum())]
        parameters.extend([dict(emission_line_A=float(x),stellar_mass=0.) for x in lib.lines[keep]])
    if scenario=='discrepancy3percent':cov+=np.diag((.03*y)**2)
    info=dict(targetid=str(targetid),z=z,scenario=scenario,
        band_flux=y.tolist(),band_covariance=cov.tolist(),retained_rest_width_A=widths.tolist(),
        cameras=contributing,emission_mask_km_s=mask_velocity,stellar_columns=int(mass.sum()),
        emission_columns=int((mass==0).sum()),age_ceiling_Gyr=float(ages.max()),
        native_resolution_row_sum_range=[min(resolution_rows),max(resolution_rows)],
        common_stellar_scale=scale)
    return matrix,y,cov,ages,mass,parameters,info


def whiten(matrix,y,cov):
    chol=np.linalg.cholesky(cov)
    return solve_triangular(chol,matrix,lower=True),solve_triangular(chol,y,lower=True)


def nonnegative(matrix,y):
    scale=np.linalg.norm(matrix,axis=0)
    assert np.all(scale>0)
    fit=lsq_linear(matrix/scale,y,bounds=(0,np.inf),method='bvls',tol=1e-11,max_iter=2000)
    coefficients=fit.x/scale
    objective=float(np.sum((matrix@coefficients-y)**2))
    if not fit.success or fit.optimality>1e-6:
        raise RuntimeError(f'NNLS failed status={fit.status}, KKT={fit.optimality}')
    return coefficients,objective,float(fit.optimality)


def bound(aw,yw,ages,mass,threshold,maximize=False):
    """Generalized Charnes-Cooper: m.u=1; gas entries have m=0."""
    n,b=aw.shape[1],len(yw)
    stellar=mass>0
    limits=[float(ages[stellar].min()),float(ages[stellar].max())]
    j=np.flatnonzero(stellar)[np.argmax(ages[stellar]) if maximize else np.argmin(ages[stellar])]
    gas=np.flatnonzero(~stellar)
    base=np.zeros(n)
    if len(gas):base[gas],_,_=nonnegative(aw[:,gas],yw)
    residual=yw-aw@base
    zero_chi2=float(residual@residual)
    if zero_chi2<threshold-1e-8:
        aa=float(aw[:,j]@aw[:,j]);ay=float(aw[:,j]@residual)
        amplitude=(ay+np.sqrt(ay**2+aa*(threshold-zero_chi2)))/aa/2
        coeff=base.copy();coeff[j]=amplitude
        return dict(age_Gyr=float(ages[j]),primal_age_Gyr=float(ages[j]),duality_gap_Gyr=0.,
            chi2=float(np.sum((aw@coeff-yw)**2)),coefficients=coeff.tolist(),
            status='AnalyticGasOrZeroInterior',retries=[],mass_normalization_error=0.)
    scales=np.linalg.norm(aw,axis=0);scales/=np.median(scales[stellar])
    aa=aw/scales
    eq=sparse.csc_matrix(np.r_[mass/scales,0.][None,:])
    nonneg=-sparse.eye(n+1,format='csc')
    soc=sparse.csc_matrix(np.vstack([np.r_[np.zeros(n),-np.sqrt(threshold)],np.column_stack([-aa,yw])]))
    constraint=sparse.vstack([eq,nonneg,soc],format='csc')
    rhs=np.r_[1.,np.zeros(n+1+b+1)]
    cones=[clarabel.ZeroConeT(1),clarabel.NonnegativeConeT(n+1),clarabel.SecondOrderConeT(b+1)]
    sign=-1 if maximize else 1
    q=np.r_[sign*ages/scales,0.]
    settings=clarabel.DefaultSettings();settings.verbose=False
    settings.tol_gap_abs=settings.tol_gap_rel=settings.tol_feas=1e-9
    rejected=[]
    for attempt,reg in enumerate([1e-8,1e-10,1e-6]):
        settings.static_regularization_constant=reg;settings.max_iter=200 if attempt==0 else 500
        solution=clarabel.DefaultSolver(sparse.csc_matrix((n+1,n+1)),q,constraint,rhs,cones,settings).solve()
        status=str(solution.status)
        if status not in ('Solved','AlmostSolved'):
            rejected.append(dict(status=status,reason='solver_status',regularization=reg));continue
        x=np.array(solution.x);u=x[:-1]/scales;t=x[-1]
        mass_error=abs(float(mass@u)-1.)
        gap=abs(solution.obj_val-solution.obj_val_dual)
        if t<=0 or min(u)<-1e-7 or mass_error>1e-6:
            rejected.append(dict(status=status,reason='simplex',regularization=reg,t=float(t),min_u=float(min(u)),mass_error=mass_error));continue
        u=np.maximum(u,0.);u/=mass@u
        coeff=u/t;objective=float(np.sum((aw@coeff-yw)**2))
        if objective>threshold+2e-5 or gap>1e-5:
            rejected.append(dict(status=status,reason='primal_or_gap',regularization=reg,chi2=objective,gap=float(gap)));continue
        return dict(age_Gyr=float(np.clip(sign*solution.obj_val_dual,*limits)),
            primal_age_Gyr=float(ages@u),duality_gap_Gyr=float(gap),chi2=objective,
            coefficients=coeff.tolist(),status=status,retries=rejected,mass_normalization_error=mass_error)
    return dict(status='NumericalFailure',retries=rejected,age_Gyr=None)


def solve(matrix,y,cov,ages,mass):
    aw,yw=whiten(matrix,y,cov)
    coeff,objective,kkt=nonnegative(aw,yw)
    threshold=float(chi2.ppf(.95,len(y)))
    report=dict(best_chi2=objective,threshold=threshold,observed_rank=len(y),nnls_kkt=kkt,
        compatible=objective<=threshold,status='incompatible' if objective>threshold else 'compatible',
        best_coefficients=coeff.tolist(),lower=None,upper=None)
    if report['compatible']:
        report['lower']=bound(aw,yw,ages,mass,threshold)
        report['upper']=bound(aw,yw,ages,mass,threshold,maximize=True)
        if any(report[s]['age_Gyr'] is None for s in ['lower','upper']):report['status']='numerical_failure'
    return report
