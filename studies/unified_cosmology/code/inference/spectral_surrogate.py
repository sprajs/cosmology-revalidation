"""Experimental local CMB-spectrum surrogate with exact CAMB background.

This is an acceleration proposal, not a replacement scientific likelihood.
Posterior intervals require separately recorded held-out likelihood errors and
exact CAMB importance correction. Outside the training envelope, use CAMB.
No training arrays or third-party code are committed to the repository.
"""
from pathlib import Path
from itertools import combinations_with_replacement
import json
import os
import numpy as np
import camb
from cobaya.theory import Theory

COORDINATES = ['100thetaMC','ombh2','omch2','logA','ns','tau','w','wa']
SPECTRA = ['tt','ee','bb','te','pp']


def polynomial_exponents(dimension=8,degree=3):
    rows = [np.zeros(dimension,dtype=int)]
    for order in range(1,degree+1):
        for indices in combinations_with_replacement(range(dimension),order):
            rows.append(np.bincount(indices,minlength=dimension))
    return np.array(rows)


def polynomial(x,exponents):
    return np.prod(np.asarray(x)[...,None,:]**exponents,axis=-1)


def coordinates(parameters,background):
    return np.array([100*background.cosmomc_theta(),parameters.ombh2,parameters.omch2,
                     np.log(1e10*parameters.InitPower.As),parameters.InitPower.ns,
                     parameters.Reion.optical_depth,parameters.DarkEnergy.w,parameters.DarkEnergy.wa])


class SpectralSurrogate(Theory):
    model_file: str = ''
    input_params = ['H0','ombh2','omch2','As','ns','tau','w','wa']
    extra_args: dict = {}
    envelope: float = 4.
    exact_mode: bool = False

    def initialize(self):
        with np.load(self.model_file,allow_pickle=False) as f:
            self.centre = f['centre'];self.chol = f['coordinate_cholesky']
            self.exponents = f['exponents'];self.coefficients = f['coefficients']
            self.output_scale = f['output_scale'];self.length = int(f['length'])
        self.required_cls = {};self.exact_calls = 0;self.surrogate_calls = 0

    def get_can_provide_params(self):
        return ['rdrag','omegam']

    def must_provide(self,**requirements):
        for key,value in requirements.items():
            if key=='Cl':
                for name,lmax in value.items():
                    name = name.lower()
                    if name not in SPECTRA+['et'] or lmax>=self.length:
                        raise ValueError(f'Unsupported spectrum request {name}:{lmax}')
                    self.required_cls[name] = max(lmax,self.required_cls.get(name,0))
            elif key not in ['Hubble','angular_diameter_distance','comoving_radial_distance','rdrag','omegam']:
                raise ValueError(f'Unsupported surrogate product {key}')

    def get_can_provide(self):
        return ['Cl','Hubble','angular_diameter_distance','comoving_radial_distance']

    def calculate(self,state,want_derived=True,**parameters):
        extra = dict(self.extra_args)
        extra['lmax'] = max(extra.get('lmax',0),max(self.required_cls.values(),default=0))
        try:
            pars = camb.set_params(**parameters,**extra)
            background = camb.get_background(pars)
        except camb.CAMBError as error:
            self.record_failure(parameters,error)
            return False
        point = coordinates(pars,background)
        x = np.linalg.solve(self.chol,point-self.centre)
        fallback = self.exact_mode or np.max(abs(x))>self.envelope
        if not fallback:
            flattened = (polynomial(x,self.exponents)@self.coefficients)*self.output_scale
            cls = dict(zip(SPECTRA,flattened.reshape(len(SPECTRA),self.length)))
            fallback = any(np.any(~np.isfinite(cls[key])) for key in SPECTRA)
            fallback |= any(np.any(cls[key][2:]<=0) for key in ['tt','ee','pp'])
        if fallback:
            try:
                full = camb.get_results(pars)
            except camb.CAMBError as error:
                self.record_failure(parameters,error)
                return False
            cmb = full.get_cmb_power_spectra(lmax=self.length-1,CMB_unit='muK',raw_cl=False)
            cls = {name:cmb['total'][:,index].copy() for name,index in [('tt',0),('ee',1),('bb',2),('te',3)]}
            cls['pp'] = cmb['lens_potential'][:,0].copy()
            background = full
            self.exact_calls += 1
        else:
            self.surrogate_calls += 1
        state['background'] = background
        state['Cl'] = cls
        state['derived'] = {'rdrag':background.get_derived_params()['rdrag'],
                            'omegam':pars.omegam}
        return True

    def record_failure(self,parameters,error):
        path = Path(self.model_file).parent/f'surrogate-camb-failures-{os.getpid()}.jsonl'
        with path.open('a') as file:
            file.write(json.dumps({'parameters':parameters,'error':str(error)},default=float)+'\n')

    def get_Cl(self,ell_factor=False,units='FIRASmuK2'):
        cls = {name:value.copy() for name,value in self.current_state['Cl'].items()}
        if units in ['FIRASmuK2','muK2']:
            factor = 1.
        elif units=='K2':
            factor = 1e-12
        elif units=='1':
            factor = 1/(2.7255e6)**2
        else:
            raise ValueError(f'Unsupported CMB unit {units}')
        ell = np.arange(self.length)
        ell2 = ell*(ell+1)
        for name in ['tt','ee','bb','te']:
            cls[name] *= factor
            if not ell_factor:
                cls[name][1:] *= 2*np.pi/ell2[1:]
        if not ell_factor:
            cls['pp'][1:] *= 2*np.pi/ell2[1:]**2
        cls['et'] = cls['te'].copy();cls['ell'] = ell
        return cls

    def get_Hubble(self,z,units='km/s/Mpc'):
        # Cobaya returns a length-one array even for scalar redshift requests.
        # Its BAO vector assembly relies on this shape contract.
        h = np.atleast_1d(self.current_state['background'].hubble_parameter(z))
        if units=='1/Mpc':return h/299792.458
        if units=='km/s/Mpc':return h
        raise ValueError(units)

    def get_angular_diameter_distance(self,z):
        return np.atleast_1d(self.current_state['background'].angular_diameter_distance(z))

    def get_comoving_radial_distance(self,z):
        return np.atleast_1d(self.current_state['background'].comoving_radial_distance(z))


def replace_theory(info,model_file,exact=False):
    old = info['theory'].pop('camb')
    info['theory']['spectral_surrogate'] = {'external':SpectralSurrogate,
        'model_file':str(Path(model_file).resolve()),'extra_args':old['extra_args'],
        'exact_mode':exact}
    return info
