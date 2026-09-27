"""Photon-counting photometry and convex, SFH-prior-free age bounds.

These are conditional inverse-problem bounds, not cosmological corrections.
"""
from pathlib import Path
import json
import numpy as np
from scipy import sparse
from scipy.optimize import lsq_linear
from astropy.cosmology import FlatLambdaCDM
import clarabel
import extinction

ROOT = Path(__file__).resolve().parents[4]
WORK = ROOT / '.work/physical-ages'
COSMO = FlatLambdaCDM(H0=70, Om0=.3)
C_A_S = 2.99792458e18


class Library:
    def __init__(self, path=WORK / 'ssp-library.npz', response_path=None):
        self.data = np.load(path)
        self.wave = self.data['wavelength_A']
        self.ages = self.data['ages_Gyr']
        self.spectra = self.data['spectra_Lsun_per_A_per_formed_Msun']
        self.mass = self.data['surviving_mass_fraction']
        self.bands=list('ugriz') if response_path is None else ['u','g','r','i','z','j','h','ks']
        response=np.load(response_path) if response_path is not None else self.data
        self.filters={b:(response[b+'_wavelength_A'],response[b+'_throughput']) if response_path else
            (response[f'sdss_{b}_wavelength'],response[f'sdss_{b}_throughput']) for b in self.bands}

    def at_redshift(self, z):
        """Keep SSPs younger than the clock; interpolate one SSP at its ceiling."""
        limit = min(float(COSMO.age(z).value), float(self.ages[-1]))
        keep = self.ages <= limit
        ages, spectra, mass = self.ages[keep], self.spectra[:, keep], self.mass[:, keep]
        if ages[-1] < limit:
            j = np.searchsorted(self.ages, limit)
            weight = np.log(limit / self.ages[j-1]) / np.log(self.ages[j] / self.ages[j-1])
            spectra = np.concatenate([spectra, ((1-weight)*self.spectra[:, j-1] + weight*self.spectra[:, j])[:, None]], axis=1)
            mass = np.concatenate([mass, ((1-weight)*self.mass[:, j-1] + weight*self.mass[:, j])[:, None]], axis=1)
            ages = np.append(ages, limit)
        return ages, spectra, mass

    def photometry(self, spectra, z, ebv, tau=0., power=-.7, integration='gauss'):
        """Return band luminosity densities up to ONE common distance factor.

        Interpolate rest Lnu, not Llambda, to observed wavelengths.
        fnu_obs=(1+z)Lnu_rest/(4*pi*DL^2).
        The omitted distance factor is common to all bands and all SSPs.
        Gauss integration resolves BOTH passband and redshifted spectral knots.
        Native/dense rules independently test the quadrature/discretization.
        """
        rest_nu = spectra * self.wave**2 / C_A_S
        host = np.exp(-tau * (self.wave / 5500.)**power)
        rest_nu = rest_nu * host
        answer = []
        for band in self.bands:
            fw,ft=self.filters[band]
            wave,quadrature=observer_quadrature(self.wave,fw,z,integration)
            throughput = np.interp(wave, fw, ft, left=0., right=0.)
            foreground = 10**(-.4*extinction.fitzpatrick99(wave.astype(float), 3.1*ebv, 3.1))
            kernel = throughput * foreground / wave
            observed = np.array([np.interp(wave/(1+z), self.wave, row) for row in rest_nu.reshape(-1, len(self.wave))])
            integral = (observed@(kernel*quadrature)) / (quadrature@(throughput/wave))
            answer.append(integral.reshape(spectra.shape[:-1])*(1+z))
        return np.stack(answer, axis=-2)

    def dn_bands(self, spectra, tau=0., power=-.7):
        """Rest-frame narrow D4000 mean Fnu, before Galactic extinction."""
        fnu = spectra * self.wave**2 / C_A_S * np.exp(-tau*(self.wave/5500.)**power)
        result = []
        for lo, hi in [(3850., 3950.), (4000., 4100.)]:
            wave = np.unique(np.r_[lo, self.wave[(self.wave>lo)&(self.wave<hi)], hi])
            vals = np.array([np.interp(wave, self.wave, row) for row in fnu.reshape(-1, len(self.wave))])
            result.append((np.trapz(vals, wave, axis=1)/(hi-lo)).reshape(spectra.shape[:-1]))
        return np.stack(result, axis=-2)


def trapezoid_weights(x):
    result = np.empty(len(x))
    result[0], result[-1] = (x[1]-x[0])/2, (x[-1]-x[-2])/2
    result[1:-1] = (x[2:]-x[:-2])/2
    return result


def observer_quadrature(rest_wave,filter_wave,z,method='gauss'):
    lo,hi=filter_wave.min(),filter_wave.max()
    if method=='native':
        wave=rest_wave[(rest_wave>=lo)&(rest_wave<=hi)]
        return wave,trapezoid_weights(wave)
    if method=='dense':
        wave=np.unique(np.r_[np.arange(lo,hi,.1),hi])
        return wave,trapezoid_weights(wave)
    shifted=rest_wave*(1+z)
    knots=np.unique(np.r_[filter_wave,shifted[(shifted>lo)&(shifted<hi)]])
    nodes,weights=np.polynomial.legendre.leggauss(4)
    half=np.diff(knots)/2;middle=(knots[1:]+knots[:-1])/2
    return (middle[:,None]+half[:,None]*nodes).ravel(),(half[:,None]*weights).ravel()


class Grid:
    """One BLAS matrix multiplication replaces thousands of interpolations."""
    def __init__(self, library=None, design=None):
        self.lib = library or Library()
        self.design = design or json.loads((Path(__file__).parent/'design.json').read_text())
        self.parameters, spectra, mass = [], [], []
        for iz, metal in enumerate(self.lib.data['metallicity']):
            for tau in self.design['host_dust_tauV']:
                for power in self.design['host_dust_power_indices'][:1 if tau == 0 else 2]:
                    self.parameters.append(dict(metallicity=float(metal), tau=float(tau), power=float(power)))
                    spectra.append(self.lib.spectra[iz]*np.exp(-tau*(self.lib.wave/5500.)**power))
                    mass.append(self.lib.mass[iz])
        self.spectra = np.asarray(spectra)
        self.mass = np.asarray(mass)
        self.nu_flat = (self.spectra*self.lib.wave**2/C_A_S).reshape(-1, len(self.lib.wave))
        self.dn = self.lib.dn_bands(self.spectra)

    def kernels(self, z, ebv):
        kernels = []
        for band in self.lib.bands:
            fw,ft=self.lib.filters[band]
            wave,quadrature=observer_quadrature(self.lib.wave,fw,z)
            response = np.interp(wave, fw, ft)/wave
            extinction_factor = 10**(-.4*extinction.fitzpatrick99(wave.astype(float), 3.1*ebv, 3.1))
            weights = quadrature*response*extinction_factor/(quadrature@response)*(1+z)
            rest = wave/(1+z)
            assert rest.min() > self.lib.wave.min() and rest.max() < self.lib.wave.max()
            right = np.searchsorted(self.lib.wave, rest); left = right-1
            frac = (rest-self.lib.wave[left])/(self.lib.wave[right]-self.lib.wave[left])
            kernel = np.zeros(len(self.lib.wave))
            np.add.at(kernel, left, weights*(1-frac)); np.add.at(kernel, right, weights*frac)
            kernels.append(kernel)
        return np.asarray(kernels)

    def at_redshift(self, z, ebv):
        matrices = (self.nu_flat@self.kernels(z, ebv).T).reshape(len(self.parameters), len(self.lib.ages), len(self.lib.bands)).transpose(0,2,1)
        limit = min(float(COSMO.age(z).value), float(self.lib.ages[-1]))
        keep = self.lib.ages <= limit
        ages, flux, dn, mass = self.lib.ages[keep], matrices[:,:,keep], self.dn[:,:,keep], self.mass[:,keep]
        if ages[-1] < limit:
            j = np.searchsorted(self.lib.ages, limit)
            weight = np.log(limit/self.lib.ages[j-1])/np.log(self.lib.ages[j]/self.lib.ages[j-1])
            def extend(array, original):
                return np.concatenate([array, ((1-weight)*original[...,j-1]+weight*original[...,j])[...,None]], axis=-1)
            flux, dn, mass = extend(flux, matrices), extend(dn, self.dn), extend(mass, self.mass)
            ages = np.append(ages, limit)
        scale = np.median(flux[:,:,np.argmin(abs(ages-1))])
        return ages, flux/scale, dn/scale, mass

    def observed_dn(self, z, ebv):
        """Fnu means with foreground extinction at the observed wavelength.

        Ratio is independent of distance and of the common (1+z) factor.
        Use complete observed band coverage; masked partial-band indices are
        not silently compared with an unmasked model.
        """
        kernels=[]
        for lo,hi in [(3850.,3950.),(4000.,4100.)]:
            wave=np.unique(np.r_[lo,self.lib.wave[(self.lib.wave>lo)&(self.lib.wave<hi)],hi])
            fore=10**(-.4*extinction.fitzpatrick99((wave*(1+z)).astype(float),3.1*ebv,3.1))
            weights=trapezoid_weights(wave)*fore/(hi-lo)
            right=np.searchsorted(self.lib.wave,wave);left=right-1
            frac=(wave-self.lib.wave[left])/(self.lib.wave[right]-self.lib.wave[left])
            kernel=np.zeros(len(self.lib.wave));np.add.at(kernel,left,weights*(1-frac));np.add.at(kernel,right,weights*frac)
            kernels.append(kernel)
        values=(self.nu_flat@np.asarray(kernels).T).reshape(len(self.parameters),len(self.lib.ages),2).transpose(0,2,1)
        limit=min(float(COSMO.age(z).value),float(self.lib.ages[-1]));keep=self.lib.ages<=limit
        result=values[:,:,keep]
        if self.lib.ages[keep][-1]<limit:
            j=np.searchsorted(self.lib.ages,limit);weight=np.log(limit/self.lib.ages[j-1])/np.log(self.lib.ages[j]/self.lib.ages[j-1])
            result=np.concatenate([result,((1-weight)*values[:,:,j-1]+weight*values[:,:,j])[:,:,None]],axis=-1)
        return result/np.median(result[:,0,:])


def feasible_fit(matrix, flux, error):
    """Nonnegative least squares, with column conditioning mapped back to mass."""
    aw, yw = matrix/error[:, None], flux/error
    scales = np.sqrt(np.sum(aw**2, axis=0))
    result=lsq_linear(aw/scales,yw,bounds=(0,np.inf),method='bvls',tol=1e-11,max_iter=1000)
    if not result.success or result.optimality>1e-6:
        raise RuntimeError(f'Bounded least-squares did not converge: {result.message}; KKT {result.optimality}')
    coeff=result.x/scales
    # Recompute in original coordinates; never trust only a solver's residual.
    return coeff, float(np.sum(((matrix@coeff-flux)/error)**2))


def age_bound(matrix, flux, error, ages, threshold, maximize=False, inequalities=None, diagnostics=None):
    """Charnes-Cooper SOCP: ||W(Au-y*t)|| <= sqrt(T)*t, 1'u=1.

    u=m/(1'm), t=1/(1'm). Extra homogeneous inequalities G*m<=0
    become G*u<=0. Returned solutions are checked in original units.
    """
    n, b = len(ages), len(flux)
    zero_chi2=float(np.sum((flux/error)**2))
    if inequalities is None and zero_chi2<threshold-1e-8:
        # A nondetection containing zero in its interior admits arbitrarily
        # small positive mass at every included age. Do not select it away.
        j=int(np.argmax(ages) if maximize else np.argmin(ages))
        column=matrix[:,j]/error;yw=flux/error
        aa=float(column@column);ay=float(column@yw)
        upper=(ay+np.sqrt(ay**2+aa*(threshold-zero_chi2)))/aa
        amplitude=upper/2
        weights=np.zeros(n);weights[j]=1.
        measured=float(np.sum(((matrix[:,j]*amplitude-flux)/error)**2))
        assert amplitude>0 and measured<=threshold+1e-8
        return dict(age=float(ages[j]),primal_value=float(ages[j]),duality_gap=0.,chi2=measured,weights=weights.tolist(),amplitude=amplitude,status='AnalyticZeroInterior',inequality_max=None,retries=0)
    aw, yw = matrix/error[:, None], flux/error
    scales=np.linalg.norm(aw,axis=0);scales/=np.median(scales)
    aw=aw/scales
    eq = np.r_[1/scales, 0.][None, :]
    nonneg = -np.eye(n+1)
    soc = np.vstack([np.r_[np.zeros(n), -np.sqrt(threshold)], np.column_stack([-aw, yw])])
    blocks, rhs = [eq, nonneg], [np.ones(1), np.zeros(n+1)]
    cones = [clarabel.ZeroConeT(1), clarabel.NonnegativeConeT(n+1)]
    if inequalities is not None:
        g = np.column_stack([inequalities/scales, np.zeros(len(inequalities))])
        # Normalize inequalities to avoid arbitrary luminosity units.
        g /= np.maximum(np.linalg.norm(g, axis=1, keepdims=True), 1e-300)
        blocks.append(g); rhs.append(np.zeros(len(g)))
        cones.append(clarabel.NonnegativeConeT(len(g)))
    blocks.append(soc); rhs.append(np.zeros(b+1)); cones.append(clarabel.SecondOrderConeT(b+1))
    settings = clarabel.DefaultSettings(); settings.verbose = False
    settings.tol_gap_abs = settings.tol_gap_rel = settings.tol_feas = 1e-9
    settings.max_iter = 200
    q = np.r_[(-1 if maximize else 1)*ages/scales, 0.]
    # A solver status is not the acceptance gate. AlmostSolved can violate
    # original-unit simplex/photometry constraints after column unscaling.
    # Try the same numerical regularizations after ANY rejected candidate;
    # never loosen either solver tolerances or physical acceptance checks.
    regularizations = [settings.static_regularization_constant, 1e-10, 1e-6]
    for retries, regularization in enumerate(regularizations):
        settings.static_regularization_constant = regularization
        settings.max_iter = 200 if retries == 0 else 500
        solver = clarabel.DefaultSolver(sparse.csc_matrix((n+1,n+1)), q,
            sparse.csc_matrix(np.vstack(blocks)), np.concatenate(rhs), cones, settings)
        sol = solver.solve()
        reason = None
        if str(sol.status) not in ('Solved', 'AlmostSolved'):
            reason = {'reason': 'status', 'status': str(sol.status)}
        else:
            x = np.asarray(sol.x); u, t = x[:-1]/scales, x[-1]
            if t <= 0 or np.min(u) < -1e-7 or abs(u.sum()-1) > 1e-6:
                reason = {'reason': 'simplex', 't': float(t),
                    'min_u': float(u.min()), 'sum_u': float(u.sum())}
            else:
                u = np.maximum(u, 0.); u /= u.sum()
                coeff = u/t
                chi2 = float(np.sum(((matrix@coeff-flux)/error)**2))
                slack = None if inequalities is None else float(np.max(
                    inequalities@u / np.maximum(np.linalg.norm(inequalities, axis=1), 1e-300)))
                gap = abs(sol.obj_val-sol.obj_val_dual)
                if chi2 > threshold+2e-5 or (slack is not None and slack > 1e-7):
                    reason = {'reason': 'primal', 'chi2': chi2,
                        'threshold': float(threshold), 'slack': slack}
                elif gap > 1e-5:
                    reason = {'reason': 'duality_gap', 'gap': float(gap)}
                else:
                    # Dual objective gives a conservative outward endpoint;
                    # the original-unit primal mixture witnesses attainability.
                    outer = (-1 if maximize else 1)*sol.obj_val_dual
                    outer = np.clip(outer, float(np.min(ages)), float(np.max(ages)))
                    return dict(age=float(outer), primal_value=float(ages@u),
                        duality_gap=float(gap), chi2=chi2, weights=u.tolist(),
                        amplitude=float(1/t), status=str(sol.status),
                        inequality_max=slack, retries=retries)
        if diagnostics is not None:
            diagnostics.append(dict(reason, attempt=retries,
                static_regularization=float(regularization)))
        if str(sol.status) == 'PrimalInfeasible':
            # A certified infeasible cell is distinct from a numerical failure.
            return None
    return None


def fieller(red, blue, variance_red, variance_blue, covariance=0., zscore=1.959963984540054):
    """Connected Fieller interval only; return None for unbounded/disjoint set."""
    a = blue**2-zscore**2*variance_blue
    b = -2*red*blue+2*zscore**2*covariance
    c = red**2-zscore**2*variance_red
    disc = b*b-4*a*c
    if a <= 0 or disc < 0:
        return None
    return ((-b-np.sqrt(disc))/(2*a), (-b+np.sqrt(disc))/(2*a))
