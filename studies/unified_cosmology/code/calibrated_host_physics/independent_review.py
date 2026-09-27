"""Independent DESI DIA-scatter, masked-band and forward-spectrum checks."""
from pathlib import Path
import hashlib
import json
import numpy as np
import extinction
import model

ROOT=Path(__file__).resolve().parents[4]
WORK=ROOT/'.work/unified-cosmology/calibrated-hosts'
OUTPUT=ROOT/'studies/unified_cosmology/results/calibrated_host_physics/independent-review.json'
RNG=np.random.default_rng(928167)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def scatter(diagonals, vector):
    """Construct y_i=sum_j R_ij x_j directly from FITS DIA column storage."""
    nd,n=diagonals.shape
    output=np.zeros(n)
    for row in range(nd):
        offset=nd//2-row
        columns=np.arange(max(0,offset),min(n,n+offset))
        output[columns-offset]+=diagonals[row,columns]*vector[columns]
    return output


def direct(arrays,z,lines,velocity,synthetic=None,ebv=0.):
    # Accumulate unnormalized native-pixel contributions over all arms before
    # normalizing each window. No sparse matrix or producer weights are used.
    widths=np.zeros(5); raw=np.zeros(5); covariance=np.zeros((5,5))
    predicted=None if synthetic is None else np.zeros((5,len(synthetic[1])))
    for camera,cutlo,cuthi in [('B',-np.inf,5780.),('R',5780.,7570.),('Z',7570.,np.inf)]:
        wave=arrays[camera+'_WAVELENGTH'].astype(float)
        flux=arrays[camera+'_FLUX'].astype(float)
        ivar=arrays[camera+'_IVAR'].astype(float)
        valid=(arrays[camera+'_MASK']==0)&np.isfinite(flux)&np.isfinite(ivar)&(ivar>0)&(wave>=cutlo)&(wave<cuthi)
        rest=wave/(1+z)
        if velocity:
            valid &= np.all(np.abs(np.log(rest[:,None]/lines[None,:]))*299792.458>velocity,axis=1)
        rawweights=np.zeros((5,len(wave)))
        for b,(lo,hi) in enumerate(model.BANDS):
            lower=np.maximum((wave-.4)/(1+z),lo)
            upper=np.minimum((wave+.4)/(1+z),hi)
            overlap=np.maximum(upper-lower,0.)*valid
            widths[b]+=overlap.sum()
            rawweights[b]=overlap*(wave**2/2.99792458e18*1e12 if b<2 else 1.)
        raw+=np.sum(rawweights*np.where(valid,flux,0.)[None,:],axis=1)
        for b in range(5):
            for c in range(5):
                covariance[b,c]+=np.sum(rawweights[b,valid]*rawweights[c,valid]/ivar[valid])
        if synthetic is not None:
            fore=10**(-.4*extinction.fitzpatrick99(wave,3.1*ebv*.86,3.1))/(1+z)
            for j,spectrum in enumerate(synthetic[1]):
                density=np.interp(rest,synthetic[0],spectrum)*fore
                pixels=scatter(arrays[camera+'_RESOLUTION'].astype(float),density)
                predicted[:,j]+=np.sum(rawweights*pixels[None,:],axis=1)
    out=(raw/widths,covariance/np.outer(widths,widths),widths)
    return out if synthetic is None else out+(predicted/widths[:,None],)


class SyntheticLibrary:
    def __init__(self,lines):
        self.wave=np.arange(3300.,4700.01,.2)
        self.lines=lines
        self.spectra=np.array([np.ones(len(self.wave)),1+.0002*(self.wave-4000),1-.45*np.exp(-.5*((self.wave-4102)/1.2)**2)])
    def grid(self,z,scenario):
        return self.spectra,[{'age_Gyr':x} for x in [.1,5.,10.]]


def relative_error(x,y):
    return float(np.max(abs(x-y)/np.maximum(1.,abs(y))))


def main():
    source=Path(model.__file__)
    initial_hash=sha(source)
    with np.load(model.WORK/'ssp-library.npz')as library:
        all_lines=library['emission_line_wavelength_A']
        lines=np.unique([x for x in all_lines if any(lo<=x<=hi for lo,hi in model.BANDS)])
    synthetic=SyntheticLibrary(lines)
    aggregate={'DIA_scatter_max_relative_error':0.,'masked_flux_max_relative_error':0.,'masked_covariance_max_relative_error':0.,'retained_width_max_absolute_error':0.,'synthetic_forward_max_relative_error':0.}
    inputs={};records=[]
    for path in sorted((WORK/'spectra').glob('*.json')):
        meta=json.loads(path.read_text());data_path=ROOT/meta['array_file']
        inputs[str(path.relative_to(ROOT))]=sha(path);inputs[str(data_path.relative_to(ROOT))]=sha(data_path)
        arrays=np.load(data_path);row={'targetid':meta['targetid'],'cases':[]}
        for camera in ['B','R','Z']:
            diagonal=arrays[camera+'_RESOLUTION'].astype(float)
            x=RNG.normal(size=diagonal.shape[1])
            err=relative_error(model.resolution(diagonal)@x,scatter(diagonal,x))
            aggregate['DIA_scatter_max_relative_error']=max(aggregate['DIA_scatter_max_relative_error'],err)
        for scenario,velocity in [('original_free_emission',0),('primary',400),('mask800',800)]:
            measured,cov,_,width=model.observed_operator(arrays,meta['desi_z'],lines,velocity)
            yy,cc,ww,predicted=direct(arrays,meta['desi_z'],lines,velocity,(synthetic.wave,synthetic.spectra),meta['fibermap']['EBV'])
            matrix,y,c,ages,mass,parameters,info=model.problem(meta['targetid'],synthetic,scenario)
            errors={'masked_flux_max_relative_error':relative_error(measured,yy),'masked_covariance_max_relative_error':relative_error(cov,cc),'retained_width_max_absolute_error':float(np.max(abs(width-ww))),'synthetic_forward_max_relative_error':relative_error(matrix[:,:3]*info['common_stellar_scale'],predicted)}
            for k,v in errors.items():aggregate[k]=max(aggregate[k],v)
            assert np.allclose(y,measured,rtol=1e-11,atol=1e-11)
            assert np.allclose(c,cov,rtol=1e-11,atol=1e-11)
            assert np.all(mass[:3]==1) and np.array_equal(ages[:3],[.1,5.,10.])
            if scenario=='original_free_emission':
                assert np.all(mass[3:]==0) and np.all(ages[3:]==0)
            row['cases'].append({'scenario':scenario,**errors})
        records.append(row)
    # A gas-only interior admits arbitrarily small positive stellar mass of
    # either extreme age. Its age ratio is otherwise undefined at zero mass.
    matrix=np.array([[1.,.4,1.],[.3,1.2,1.]])
    observed=np.ones(2);ages=np.array([.1,10.,0.]);mass=np.array([1.,1.,0.])
    lower=model.bound(matrix,observed,ages,mass,1.)
    upper=model.bound(matrix,observed,ages,mass,1.,maximize=True)
    assert lower['age_Gyr']==.1 and upper['age_Gyr']==10.
    assert lower['status']==upper['status']=='AnalyticGasOrZeroInterior'
    assert lower['chi2']<=1 and upper['chi2']<=1
    assert sha(source)==initial_hash
    assert max(aggregate.values())<1e-8,aggregate
    result={'schema':'calibrated-host-physics-independent-review-v1','status':'passed','spectra':len(records),'masked_scenarios_per_spectrum':3,'seed':928167,'metrics':aggregate,'gas_only_semantics':{'age_interval_Gyr':[lower['age_Gyr'],upper['age_Gyr']],'meaning':'Gas coefficients carry zero stellar mass; a gas-only interior leaves the full positive-stellar-mass age range unidentified.'},'scope':'Independent DIA column scatter, same native data/model camera and emission masks, overlapping-band formal diagonal-IVAR covariance, three analytic spectral shapes including a narrow absorption line, and gas-mass ratio semantics. No stellar-library or astrophysical completeness validation.','records':records,'source_sha256':{str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),source,model.WORK/'desispec-resolution.py']},'input_sha256':inputs}
    OUTPUT.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items()if k not in ['records','source_sha256','input_sha256']},indent=2))


if __name__=='__main__':main()
