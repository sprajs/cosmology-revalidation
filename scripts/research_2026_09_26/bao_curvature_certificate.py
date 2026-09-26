"""Conservative Gaussian marginal-interval certificate, arbitrary FLRW curvature."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.stats import norm
from scipy.integrate import quad
from bao_shape import ROOT, INPUT, OUT
from bao_recent_counterexample import geometry

PROTOCOL = ROOT/'docs/research-2026-09-26/bao-curvature-protocol.md'


def certificate(b, B, g, D):
    a=b/(1+b);t=np.log1p(b)
    if min(B,g,D)<=0:
        return {'valid':False,'reason':'nonpositive confidence limit'}
    upper_dm=np.sqrt(B**3/(a*g));x=t*g/D
    if D<=upper_dm or x>=np.pi/2:
        return {'valid':False,'reason':'branch/domain gate','upper_DM_BGS':upper_dm,'phase_lower':x}
    lower_B=(a*g*(D*np.sin(x))**2)**(1/3)
    return {'valid':True,'upper_DM_BGS':upper_dm,'phase_lower':x,
            'minimum_BGS':lower_B,'deficit':lower_B-B,
            'closed_to_flat_bound_ratio':np.sinc(x/np.pi)**(2/3)}


def validate():
    rng=np.random.default_rng(2026092641)
    edges=np.array([0,.295,.51,.706,.934,1.321,1.484,2.33])
    z=edges[1:];b=z[0]
    counts={'flat':0,'open':0,'closed':0,'closed_high_after_turnover':0}
    worst=0.;quadrature=0.
    for repeat in range(5000):
        q=rng.uniform(0,4,len(z));scale=np.exp(rng.uniform(np.log(20),np.log(200)))
        dc,dh=geometry(z,np.log(scale),q,edges)
        radius=dc[-1]/rng.uniform(.05,np.pi-.02)
        for label in ['flat','open','closed']:
            dm=dc if label=='flat' else radius*(np.sinh(dc/radius) if label=='open' else np.sin(dc/radius))
            B=(b*dm[0]**2*dh[0])**(1/3)
            for j in range(1,len(z)):
                g=(1+z[j])*dh[j]
                test=certificate(b,B,g,dm[-1])
                if not test['valid']:continue
                counts[label]+=1
                if label=='closed' and dc[-1]/radius>np.pi/2:
                    counts['closed_high_after_turnover']+=1
                worst=max(worst,test['deficit'])
                assert test['deficit']<1e-8,(repeat,label,test)
        if repeat<100:
            te=np.log1p(edges)
            def invE(x):
                t=np.log1p(x)
                return np.exp(-sum(v*max(0,min(t,hi)-lo) for v,lo,hi in zip(q,te[:-1],te[1:])))/(1+x)
            integral=scale*quad(invE,0,b,epsabs=1e-11,epsrel=1e-11)[0]
            quadrature=max(quadrature,abs(integral-dc[0]))
    assert quadrature<1e-8
    x=1e-6
    assert abs(np.sinc(x/np.pi)**(2/3)-1)<1e-12
    return dict(random_histories=5000,seed=2026092641,valid_certificates=counts,
                maximum_false_deficit=worst,quadrature_max_error=quadrature,
                flat_limit_error=float(abs(np.sinc(x/np.pi)**(2/3)-1)))


def main():
    out=OUT/'curvature_certificate';out.mkdir(exist_ok=False)
    (out/'executed_source.py').write_bytes(Path(__file__).read_bytes())
    (out/'protocol.md').write_bytes(PROTOCOL.read_bytes())
    inp=[INPUT/'desi_gaussian_bao_ALL_GCcomb_mean.txt',INPUT/'desi_gaussian_bao_ALL_GCcomb_cov.txt']
    data=pd.read_csv(inp[0],sep=r'\s+',comment='#',names=['z','value','kind'])
    cov=np.loadtxt(inp[1]);Bidx=int(data.index[data.kind=='DV_over_rs'][0]);b=float(data.loc[Bidx,'z'])
    high=data[data.kind=='DM_over_rs'].sort_values('z').iloc[-1];Didx=int(high.name)
    B=float(data.loc[Bidx,'value']);sigmaB=np.sqrt(cov[Bidx,Bidx]);D=float(high.value);sigmaD=np.sqrt(cov[Didx,Didx])
    assert np.isclose(b,.295) and np.isclose(high.z,2.33)
    verification=validate();rows=[]
    for idx,row in data[data.kind=='DH_over_rs'].iterrows():
        g=(1+row.z)*row.value;sigmag=(1+row.z)*np.sqrt(cov[idx,idx])
        central=certificate(b,B,g,D)
        def at(k):return certificate(b,B+k*sigmaB,g-k*sigmag,D-k*sigmaD)
        def deficit(k):
            test=at(k)
            return test.get('deficit',-1.) if test['valid'] else -1.
        kmax=0.
        if deficit(0)>0:
            upper=.25
            while deficit(upper)>0 and upper<20:upper+=.25
            kmax=brentq(deficit,0,upper,xtol=1e-12)
            assert at(kmax)['valid']
        rows.append(dict(z=float(row.z),g=float(g),sigma_g=float(sigmag),central=central,
                         certificate_k=kmax,one_radial_union_bound=min(1.,3*norm.sf(kmax)),
                         six_radial_union_bound=min(1.,8*norm.sf(kmax)),
                         confidence_limits_at_boundary=dict(B_upper=B+kmax*sigmaB,g_lower=g-kmax*sigmag,D_lower=D-kmax*sigmaD),
                         boundary=at(kmax)))
    best=max(rows,key=lambda row:row['certificate_k'])
    result=dict(scope='Post-result exploratory confidence certificate; no curvature prior, no cosmology posterior or present-q localization.',
                assumptions=['Homogeneous/isotropic metric FLRW geometry','positive expansion','one constant comoving ruler','closed geometry on positive first-antipode branch through z2.33','quoted Gaussian marginal errors correct; no extra mean biases or variance'],
                BGS=dict(z=b,value=B,sigma=float(sigmaB)),highest_transverse=dict(z=float(high.z),value=D,sigma=float(sigmaD)),
                radial_results=rows,best=best,verification=verification)
    (out/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    (out/'manifest.json').write_text(json.dumps(dict(inputs_sha256={str(p.relative_to(ROOT)):sha(p) for p in inp+[PROTOCOL,Path(__file__),Path(__file__).with_name('bao_recent_counterexample.py'),Path(__file__).with_name('bao_shape.py')]},outputs_sha256={p.name:sha(p) for p in out.iterdir() if p.is_file()}),indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
