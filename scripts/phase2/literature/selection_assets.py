#!/usr/bin/env python3
"""Read released DES selection/population assets; never infer efficiency from selected mocks.

GENPDF densities are table-defined and normalized numerically on their generated
coordinate. Host-z uses the release table. Neither function supplies SALT-fit cuts.
"""
from pathlib import Path
import gzip, numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.stats import truncnorm
ROOT=Path(__file__).resolve().parents[3]
AUX=ROOT/'phase2/official/inputs/SNDATA_ROOT'

def text(path):
    return gzip.open(path,'rt').read() if str(path).endswith('.gz') else Path(path).read_text()

def grids(path, row='PDF:'):
    out=[]; names=None; rows=[]
    for line in text(path).splitlines():
        w=line.split()
        if not w:continue
        if w[0]=='VARNAMES:':
            if names is not None:out.append((names,np.array(rows)))
            names=w[1:]; rows=[]
        elif w[0]==row:rows.append([float(x) for x in w[1:len(names)+1]])
    if names is not None:out.append((names,np.array(rows)))
    return out

def regular(names, data):
    axes=[np.unique(data[:,i]) for i in range(len(names)-1)]
    shape=tuple(map(len,axes)); vals=np.full(shape,np.nan)
    for row in data:
        vals[tuple(np.searchsorted(a,v) for a,v in zip(axes,row[:-1]))]=row[-1]
    if np.isnan(vals).any():raise ValueError('Incomplete grid')
    return axes,vals

class PopulationPDF:
    """Conditional generated-variable densities; outside generated support => zero.

    Context (host mass/redshift) is clipped at table edges as GENPDF_OPTMASK=1
    permits in the historical simulator. Reproduction still requires historical
    reference-file equality and validation against the original generator.
    """
    def __init__(self,path=AUX/'models/population_pdf/DES-SN5YR/DES-SN5YR_DES_S3_P21.DAT.gz'):
        self.path=Path(path);self.maps={}
        for names,data in grids(path):
            axes,vals=regular(names,data)
            norm=np.trapezoid(vals,axes[0],axis=0)
            self.maps[names[0]]=(names,axes,RegularGridInterpolator(axes,vals,bounds_error=False,fill_value=0),norm)
    def density(self,name,value,**context):
        names,axes,interp,norm=self.maps[name]
        vals=np.broadcast_arrays(np.asarray(value),*[np.asarray(context[n]) for n in names[1:-1]])
        pts=np.stack(vals,axis=-1).reshape(-1,len(axes)).astype(float)
        for i,axis in enumerate(axes[1:],1):pts[:,i]=np.clip(pts[:,i],axis[0],axis[-1])
        if len(axes)==1:denom=float(norm)
        else:denom=RegularGridInterpolator(axes[1:],norm,bounds_error=True)(pts[:,1:])
        return (interp(pts)/denom).reshape(vals[0].shape)
    def log_joint(self,c_intrinsic,x1,rv,ebv,beta,logmass,z_host):
        ps=[self.density('SALT2c',c_intrinsic),self.density('SALT2x1',x1,LOGMASS=logmass),self.density('RV',rv,LOGMASS=logmass),self.density('EBV',ebv,ZTRUE=z_host,LOGMASS=logmass),truncnorm.pdf(beta,(.4-2.07646589)/.21728215,(3-2.07646589)/.21728215,loc=2.07646589,scale=.21728215)]
        with np.errstate(divide='ignore'):return np.sum([np.log(p) for p in ps],axis=0)

class HostEfficiency:
    def __init__(self,path=AUX/'models/searcheff/SEARCHEFF_zHOST_DES-SN5YR.DAT'):
        self.path=Path(path);self.maps=[];rows=[];epoch=None;fields=None
        def finish():
            if rows:
                axes,vals=regular(['r_obs_auto','obs_gr_auto','HOSTEFF'],np.array(rows))
                self.maps.append((epoch,fields,axes,RegularGridInterpolator(axes,vals,bounds_error=True)))
        for line in text(path).splitlines():
            w=line.split()
            if not w:continue
            if w[0]=='PEAKMJD_RANGE:':finish();rows=[];epoch=tuple(map(float,w[1:3]))
            elif w[0]=='FIELDLIST:':fields=w[1].split('+')
            elif w[0]=='HOSTEFF:':rows.append(list(map(float,w[1:4])))
        finish()
    def probability(self,r,gr,field,peakmjd,shift_mag=0):
        for (lo,hi),fields,axes,interp in self.maps:
            if field in fields and lo<=peakmjd<hi:
                # Positive shift samples the table at a fainter host magnitude.
                point=[np.clip(r+shift_mag,axes[0][0],axes[0][-1]),np.clip(gr,axes[1][0],axes[1][-1])]
                return float(np.clip(interp([point])[0],0,1))
        raise ValueError(f'No host-efficiency map for {field=} {peakmjd=}')

class DetectionEfficiency:
    def __init__(self,path=AUX/'models/searcheff/SEARCHEFF_PIPELINE_DES.DAT'):
        self.path=Path(path);self.maps={};filt=None
        for line in text(path).splitlines():
            w=line.split()
            if not w:continue
            if w[0]=='FILTER:':filt=w[1];self.maps[filt]=[]
            elif w[0]=='SNR:':self.maps[filt].append(tuple(map(float,w[1:3])))
        self.maps={k:np.array(v) for k,v in self.maps.items()}
    def probability(self,band,snr):
        tab=self.maps[band]
        # Historical GETEFF_PIPELINE_DETECT sets negative SNR to zero efficiency.
        return np.where(np.asarray(snr)<=0,0,np.interp(snr,tab[:,0],tab[:,1]))
    def trigger_probability(self,mjd,band,snr,separation=.4):
        """DES two independent detection-period successes, grouped per SNANA 2fe0f56.
        Call with valid, nonsaturated observations and the simulator's SNR;
        undefined-model epochs/extra PHOTPROB cuts must be removed upstream.
        This is conditional on cadence and SNR, not a population selection rate.
        """
        order=np.argsort(mjd); t=np.asarray(mjd)[order]
        p=np.array([self.probability(str(b).strip(),s) for b,s in zip(np.asarray(band)[order],np.asarray(snr)[order])],float)
        if len(t)==0:return 0.
        groups=[];start=t[0];q=1.
        for ti,pi in zip(t,p):
            if ti-start>separation:groups.append(1-q);start=ti;q=1.
            q*=1-pi
        groups.append(1-q)
        p0,p1=1.,0.
        for pg in groups:p0,p1=p0*(1-pg),p1*(1-pg)+p0*pg
        return float(np.clip(1-p0-p1,0,1))

def importance_diagnostics(logweights):
    """Weights only; does not certify correct p0 or matching sample conditioning."""
    x=np.asarray(logweights,float)
    if len(x)==0 or not np.isfinite(x).any():raise ValueError('No support')
    a=np.max(x);w=np.exp(x-a);total=w.sum()
    return {'n':len(w),'ess':float(total**2/(w@w)),'ess_fraction':float(total**2/(w@w)/len(w)),'max_normalized_weight':float(w.max()/total),'log_mean_weight':float(a+np.log(total/len(w)))}
