"""Independent calibrated-flux integration. No SNANA library/fitter calls.

The direct path uses scipy tensor splines, polynomial colour law, explicit
photon integrals and the calibration primary SED. The sncosmo path separately
uses its native SALT3 interpolator/integrator. Shared assets are intentional.
"""
from __future__ import annotations
import gzip
import hashlib
import importlib.metadata
import json
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.io import fits
import extinction
import sncosmo
from scipy.interpolate import RectBivariateSpline, RegularGridInterpolator

ROOT = Path(__file__).resolve().parents[3]
BUNDLE = ROOT / 'phase2/official/portable_pilot'
ASSETS = BUNDLE / 'assets'
OUT = ROOT / 'phase2/independent_flux'
MAG_OFFSET = .27
LB, LV = 4302.57, 5428.55


def grid(path):
    a = np.loadtxt(path)
    t, w = np.unique(a[:, 0]), np.unique(a[:, 1])
    assert len(a) == len(t) * len(w)
    return t, w, a[:, 2].reshape(len(t), len(w))


class Engine:
    def __init__(self):
        self.direct = {}
        for k in (0, 1):
            t, w, v = grid(ASSETS / f'salt3_template_{k}.dat.gz')
            self.direct[k] = RectBivariateSpline(t, w, v, kx=3, ky=3)
        self.tmin, self.tmax, self.wmin, self.wmax = t[0], t[-1], w[0], w[-1]
        self.err = {}
        for key, filename in [('00','variance_0'), ('11','variance_1'), ('01','covariance_01')]:
            t, w, v = grid(ASSETS / f'salt3_lc_model_{filename}.dat.gz')
            self.err[key] = RegularGridInterpolator((t,w),v,bounds_error=True)
        self.cd_w, self.cd = np.loadtxt(ASSETS / 'salt3_color_dispersion.dat.gz',unpack=True)
        with gzip.open(ASSETS / 'salt3_color_correction.dat.gz','rt') as f:
            lines = f.read().split()
        n = int(lines[0]); coefs = np.array(lines[1:n+1],dtype=float)
        self.clpoly = np.r_[0,1-coefs.sum(),coefs]
        with fits.open(ASSETS / 'calib_DES-SN5YR_DES.fits.gz') as f:
            self.wave = np.array(f['FilterTrans'].data.field(0),float)
            self.trans = {b:np.array(f['FilterTrans'].data.field(j+1),float) for j,b in enumerate('griz')}
            self.primary = np.array(f['PrimarySED'].data['AB'],float)
            self.primary_mag = {r['Filter Name'][-1]: float(r['Primary Mag']) for r in f['ZPoff'].data}
        self.bands = {}
        for b,tr in self.trans.items():
            nz = np.where(tr > 0)[0]
            sl = slice(max(0,nz[0]-1), min(len(tr),nz[-1]+2))
            self.bands[b] = sncosmo.Bandpass(self.wave[sl],tr[sl],name='des5yr_exact_'+b)
        kwargs = {'m0file':'salt3_template_0.dat.gz','m1file':'salt3_template_1.dat.gz',
                  'clfile':'salt3_color_correction.dat.gz','cdfile':'salt3_color_dispersion.dat.gz',
                  'lcrv00file':'salt3_lc_model_variance_0.dat.gz',
                  'lcrv11file':'salt3_lc_model_variance_1.dat.gz',
                  'lcrv01file':'salt3_lc_model_covariance_01.dat.gz'}
        handles = {k:gzip.open(ASSETS/v,'rt') for k,v in kwargs.items()}
        try:
            source = sncosmo.SALT3Source(**handles,name='DES5YR-local')
        finally:
            for f in handles.values(): f.close()
        self.snmodel = sncosmo.Model(source=source,effects=[sncosmo.F99Dust(3.1)],
                                    effect_names=['mw'],effect_frames=['obs'])

    def colorlaw(self,w):
        x = (np.asarray(w)-LB)/(LV-LB)
        lo,hi = (np.array([2800.,8000.])-LB)/(LV-LB)
        xc = np.clip(x,lo,hi)
        p = np.polynomial.polynomial.polyval(xc,self.clpoly)
        dp = np.polynomial.polynomial.polyval(xc,np.polynomial.polynomial.polyder(self.clpoly))
        return -(p+(x-xc)*dp)

    def prepare(self,z,ebv,step=5.):
        result = {}
        for band,tr in self.trans.items():
            nz = np.where(tr>0)[0]
            lo,hi = self.wave[max(0,nz[0]-1)], self.wave[min(len(tr)-1,nz[-1]+1)]
            wave = np.linspace(lo,hi,int(round((hi-lo)/step))+1)
            trans = np.interp(wave,self.wave,tr)
            # Trapezoid is equivalent to rectangular for zero end transmission.
            dw = wave[1]-wave[0]
            rest = wave/(1+z)
            good = (rest>self.wmin)&(rest<self.wmax)
            primary = np.interp(wave,self.wave,self.primary)
            denom = np.sum(primary*trans*wave)*dw
            assert denom>0
            mw = 10**(-.4*extinction.fitzpatrick99(wave,ebv*3.1,3.1))
            weights = 1e-12/(1+z)*wave*trans*mw*dw/denom*10**(.4*(27.5-self.primary_mag[band]-MAG_OFFSET))
            cl = self.colorlaw(rest)
            f99 = extinction.fitzpatrick99(rest,3.1,3.1) - extinction.fitzpatrick99(np.array([LB]),3.1,3.1)[0]
            result[band] = {'wave':wave[good], 'rest':rest[good], 'weights':weights[good],
                            'cl':cl[good], 'f99':f99[good],
                            'norm_weights':(wave*trans)[good]/np.sum(wave*trans),
                            'meanrest':np.sum(wave*trans)/np.sum(trans)/(1+z)}
        return result

    def flux(self,p,bands,times,z,prepared,family='salt',interpolation='scipy'):
        """p=(ln x0,x1,c,t0[,d]), times observer MJD, flux at zeropoint27.5."""
        amp,x1,c,t0 = p[:4]
        times,bands = np.asarray(times),np.asarray(bands)
        phase = (times-t0)/(1+z)
        if np.min(phase)<self.tmin or np.max(phase)>self.tmax:
            raise ValueError('Phase extrapolation forbidden')
        ans = np.empty(len(times))
        for b in np.unique(bands):
            ix = np.where(bands==b)[0]; g=prepared[b]
            if interpolation=='scipy':
                shape = self.direct[0](phase[ix],g['rest'],grid=True) + x1*self.direct[1](phase[ix],g['rest'],grid=True)
            else:
                shape = (self.snmodel.source._model['M0'](phase[ix],g['rest']) + x1*self.snmodel.source._model['M1'](phase[ix],g['rest']))/1e-12
            color = c if family!='phase_colour' else c+p[4]*np.tanh(phase[ix,None]/20)
            law = g['f99'] if family=='f99' else g['cl']
            ans[ix] = np.exp(amp)*np.sum(shape * 10**(-.4*color*law) * g['weights'],axis=1)
        return ans

    def sncosmo_flux(self,p,bands,times,z,ebv,covariance=False):
        self.snmodel.set(z=z,t0=p[3],x0=np.exp(p[0]),x1=p[1],c=p[2],mwebv=ebv)
        obs_bands = np.array([self.bands[b] for b in bands],dtype=object)
        scales = np.array([10**(-.4*(MAG_OFFSET+self.primary_mag[b])) for b in bands])
        if covariance:
            fl,cv = self.snmodel.bandfluxcov(obs_bands,times,zp=27.5,zpsys='ab')
            return fl*scales, cv*np.outer(scales,scales)
        return self.snmodel.bandflux(obs_bands,times,zp=27.5,zpsys='ab')*scales

    def model_covariance(self,p,bands,times,z,prepared,mode='snana_formula'):
        """Independent analytic construction, not claimed exact frozen SNANA C.

        The snana_formula branch retains actual source's colour/z-dependent
        flux_train denominator despite its comment describing c=0. It applies
        linear error-map interpolation, same-band colour-dispersion covariance,
        and the maximum-with-0.005mag floor. MW uncertainty added separately.
        """
        phase=(np.asarray(times)-p[3])/(1+z); bands=np.asarray(bands)
        f=self.flux(p,bands,times,z,prepared,interpolation='sncosmo')
        snake=np.zeros(len(f)); cd=np.zeros(len(f))
        for b in np.unique(bands):
            ix=np.where(bands==b)[0]; g=prepared[b]
            shape=(self.snmodel.source._model['M0'](phase[ix],g['rest'])+p[1]*self.snmodel.source._model['M1'](phase[ix],g['rest']))/1e-12
            if mode=='snana_formula':
                shape=shape*10**(-.4*p[2]*g['cl'])/(1+z)
            train=shape@g['norm_weights']
            coords=np.column_stack([phase[ix],np.full(len(ix),g['meanrest'])])
            v=self.err['00'](coords)+2*p[1]*self.err['01'](coords)+p[1]**2*self.err['11'](coords)
            if np.any(v<0):
                raise ValueError('Negative model variance; must audit before choosing fallback')
            snake[ix]=v/train**2
            cd[ix]=np.interp(g['meanrest'],self.cd_w,self.cd)
        total=np.maximum(snake+cd**2,(.005*np.log(10)/2.5)**2)
        cv=np.diag(f**2*(total-cd**2))
        cv+=np.outer(f*cd,f*cd)*(bands[:,None]==bands[None,:])
        return f,cv,snake,cd


def read_inputs():
    return (pd.read_csv(BUNDLE/'calibrated_observations.csv',dtype={'CID':str}),
            pd.read_csv(BUNDLE/'snana_parameters.csv',dtype={'CID':str}).set_index('CID'),
            json.loads((BUNDLE/'double_parameters_and_hessian.json').read_text()))


def mask_observations(observations,which='published'):
    ref=pd.read_csv(BUNDLE/f'{which}_predictions_and_mask.csv',dtype={'CID':str})
    ref=ref[ref.DATA_MODEL==1]
    matched=[]
    for r in ref.itertuples():
        q=observations[(observations.CID==r.CID)&(observations.BAND==r.BAND)]
        d=abs(q.MJD-r.MJD)
        if len(d)==0 or d.min()>.003: raise ValueError(f'Mask unmatched {r}')
        if np.sum(d<.003)!=1: raise ValueError(f'Mask ambiguous {r}')
        matched.append(d.idxmin())
    assert len(set(matched))==len(matched)
    q=observations.loc[matched].copy().sort_values(['CID','MJD','BAND'])
    q['epoch_fold']=[int(hashlib.sha256(f'independent-flux-v1|{cid}|{int(np.floor(t+.5))}'.encode()).hexdigest(),16)%5 for cid,t in zip(q.CID,q.MJD)]
    return q


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def provenance(name,outputs,extra_inputs=()):
    inputs=[p for p in BUNDLE.rglob('*') if p.is_file()]
    inputs += list(extra_inputs)
    code=list(Path(__file__).parent.glob('*.py'))
    record={'inputs_sha256':{str(p.relative_to(ROOT)):sha(p) for p in inputs},
            'code_sha256':{str(p.relative_to(ROOT)):sha(p) for p in code},
            'outputs_sha256':{str(p.relative_to(ROOT)):sha(p) for p in outputs},
            'packages':{p:importlib.metadata.version(p) for p in ['numpy','scipy','sncosmo','astropy','extinction','iminuit','pandas']}}
    (OUT/f'{name}_manifest.json').write_text(json.dumps(record,indent=2)+'\n')
