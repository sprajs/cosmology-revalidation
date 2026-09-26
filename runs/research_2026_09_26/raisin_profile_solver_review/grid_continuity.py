"""Execute exact source-extracted native SNooPy grid interpolation on fixed coordinates.

No SN photometry is read. A minimal adapter loads the released FITS arrays
with the native reader's one-based/float32 conventions. This is a rest-grid
component check, not an observer-flux/KCOR or complete objective closure.
"""
from pathlib import Path
import ctypes as ct
import hashlib
import json
import subprocess
import time
import numpy as np
from astropy.io import fits

ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).resolve().parent
SRC=ROOT/'phase2/official/build/SNANA-current/src'
MODEL=ROOT/'sources/repos/djones1040__RAISIN_DataRelease@a383c4b/model/snoopy.B18/SNooPy_B18.fits'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
save=lambda p,v:Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')

def main():
    start=time.monotonic()
    native=(SRC/'genmag_snoopy.c').read_text()
    interp=native[native.index('void gridinterp_snoopy('):]
    reader=(SRC/'sntools_modelgrid_read.c').read_text()
    index=reader[reader.index('int INDEX_GRIDGEN('):reader.index('// end of INDEX_GRIDGEN')]
    header=(SRC/'sntools_modelgrid.h').read_text()
    prefix='''#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MXPATHLEN 1024
#define MXFILTINDX 100
#define SEV_FATAL 1
typedef void fitsfile;
char c1err[2048],c2err[2048];
void errmsg(int a,int b,char*c,char*d,char*e){fprintf(stderr,"%s: %s %s\\n",c,d,e);abort();}
'''
    wrapper='''
SNGRID_DEF SNGRID_SNOOPY;
void load_arrays(int ns,float*s,int nt,float*t,int*ptr,short*mag,short*err,float pack){
 memset(&SNGRID_SNOOPY,0,sizeof(SNGRID_SNOOPY));
 int pars[2]={IPAR_GRIDGEN_SHAPEPAR,IPAR_GRIDGEN_TREST};
 int sizes[2]={ns,nt};float*vals[2]={s,t};
 for(int q=0;q<2;q++){int p=pars[q],n=sizes[q];
  SNGRID_SNOOPY.NBIN[p]=n;
  for(int j=0;j<n;j++)SNGRID_SNOOPY.VALUE[p][j+1]=vals[q][j];
  SNGRID_SNOOPY.VALMIN[p]=vals[q][0];SNGRID_SNOOPY.VALMAX[p]=vals[q][n-1];
  float dif=vals[q][n-1]-vals[q][0];SNGRID_SNOOPY.BINSIZE[p]=dif/(float)(n-1);
 }
 SNGRID_SNOOPY.PTR_GRIDGEN_LC=ptr;
 SNGRID_SNOOPY.I2GRIDGEN_LCMAG=mag;SNGRID_SNOOPY.I2GRIDGEN_LCERR=err;
 GRIDGEN_I2LCPACK=pack;
}
'''
    cfile=OUT/'source-extracted-grid.c'
    cfile.write_text(prefix+header+wrapper+index+interp)
    with fits.open(MODEL) as h:
        shape=h['LUMI-GRID'].data.field(0).astype('f4')
        phase=h['TREST-GRID'].data.field(0).astype('f4')
        ptr=np.r_[np.int32(0),h['PTR_I2LCMAG'].data.field(0).astype('i4')]
        mag=np.r_[np.int16(0),h['I2LCMAG'].data.field(0).astype('i2')]
        err=np.r_[np.int16(0),h['I2LCMAG'].data.field(1).astype('i2')]
        pack=float(h[0].header['I2LCPACK'])
    phases=np.array([-7.,0.,10.,30.,45.])
    protocol={'scope':'Pure rest-frame released grid; no measured SN flux/objective outcomes',
              'sources_sha256':{str(p.relative_to(ROOT)):sha(p) for p in [SRC/'genmag_snoopy.c',SRC/'sntools_modelgrid_read.c',SRC/'sntools_modelgrid.h',MODEL,Path(__file__),cfile]},
              'phases':phases.tolist(),'filters':'BVugriYJH','shape':'all129 internal algorithm cell boundaries',
              'one_sided_steps_cell_fraction':[1e-3,1e-4],
              'controls':'Unchanged C function bodies; independent Python bilinear implementation; full two-sided grid boundary scan; unsupported rest magnitudes reported, not clipped'}
    assert not (OUT/'grid-protocol.json').exists()
    save(OUT/'grid-protocol.json',protocol)
    command=['cc','-O2','-fPIC','-shared',str(cfile),'-o',str(OUT/'grid.so')]
    subprocess.run(command,check=True,capture_output=True)
    lib=ct.CDLL(str(OUT/'grid.so'))
    fp=np.ctypeslib.ndpointer(dtype='f4',flags='C_CONTIGUOUS')
    dp=np.ctypeslib.ndpointer(dtype='f8',flags='C_CONTIGUOUS')
    ip=np.ctypeslib.ndpointer(dtype='i4',flags='C_CONTIGUOUS')
    sp=np.ctypeslib.ndpointer(dtype='i2',flags='C_CONTIGUOUS')
    lib.load_arrays.argtypes=[ct.c_int,fp,ct.c_int,fp,ip,sp,sp,ct.c_float]
    lib.gridinterp_snoopy.argtypes=[ct.c_int,ct.c_double,ct.c_int,dp,dp,dp]
    lib.load_arrays(len(shape),shape,len(phase),phase,ptr,mag,err,pack)
    sb=float(np.float32(np.float32(shape[-1]-shape[0])/np.float32(len(shape)-1)))
    tb=float(np.float32(np.float32(phase[-1]-phase[0])/np.float32(len(phase)-1)))
    bounds=float(shape[0])+np.arange(1,len(shape)-1)*sb
    def native_eval(s):
        m=np.empty((9,len(phases)));e=np.empty_like(m)
        for band in range(9):lib.gridinterp_snoopy(band,s,len(phases),phases,m[band],e[band])
        return np.stack([m,e])
    def python_eval(s):
        # Independent direct four-corner expression, with exact native indexes.
        k=min(max(int((s-float(shape[0]))/sb),0),len(shape)-2)
        rs=(s-float(shape[k]))/sb
        ans=np.empty((2,9,len(phases)))
        for b in range(9):
            for j,t in enumerate(phases):
                l=min(max(int((t-float(phase[0]))/tb),0),len(phase)-2)
                rt=(t-float(phase[l]))/tb
                for v,table in enumerate([mag,err]):
                    a=table[ptr[k+1]+b*len(phase)+2+l]/pack
                    b0=table[ptr[k+1]+b*len(phase)+3+l]/pack
                    c=table[ptr[k+2]+b*len(phase)+2+l]/pack
                    d=table[ptr[k+2]+b*len(phase)+3+l]/pack
                    ans[v,b,j]=(1-rs)*((1-rt)*a+rt*b0)+rs*((1-rt)*c+rt*d)
        return ans
    mismatches=[]; scan=[]
    for k,bound in enumerate(bounds):
        for frac in [1e-3,1e-4]:
            step=sb*frac
            ml2,ml,mr,mr2=[native_eval(bound+x*step) for x in [-2,-1,1,2]]
            lm=2*ml-ml2;rm=2*mr-mr2
            left=(ml-ml2)/step;right=(mr2-mr)/step
            if frac==1e-3:
                mismatches.append(float(np.max(np.abs(ml-python_eval(bound-step)))))
            scan.append(dict(k=k+1,bound=bound,frac=frac,limit_left=lm,limit_right=rm,
                             jump=rm-lm,derivative_left=left,derivative_right=right))
    stack=lambda key:np.stack([x[key] for x in scan])
    jumps=stack('jump');dl=stack('derivative_left');dr=stack('derivative_right')
    values=stack('limit_left')
    worst=np.unravel_index(np.argmax(abs(jumps[:,0])),jumps[:,0].shape)
    worstder=np.unravel_index(np.argmax(abs(dr[:,0]-dl[:,0])),jumps[:,0].shape)
    results={'protocol_sha256':sha(OUT/'grid-protocol.json'),'compiled_binary_sha256':sha(OUT/'grid.so'),
             'command':command,'compiler':subprocess.check_output(['cc','--version'],text=True).splitlines()[0],
             'python_bilinear_max_difference':max(mismatches),
             'shape_float32_step':sb,'shape_knots':len(shape),'phase_knots':len(phase),
             'stored_vs_algorithm_boundary_max_abs':float(np.max(abs(bounds-shape[1:-1]))),
             'max_magnitude_limit_jump':float(np.max(abs(jumps[:,0]))),
             'max_magerr_limit_jump':float(np.max(abs(jumps[:,1]))),
             'max_magnitude_derivative_jump_mag_per_stretch':float(np.max(abs(dr[:,0]-dl[:,0]))),
             'worst_jump':{'boundary':scan[worst[0]]['bound'],'band':'BVugriYJH'[worst[1]],'phase':phases[worst[2]],'jump_mag':jumps[worst[0],0,worst[1],worst[2]]},
             'worst_derivative_jump':{'boundary':scan[worstder[0]]['bound'],'band':'BVugriYJH'[worstder[1]],'phase':phases[worstder[2]],'left':dl[worstder[0],0,worstder[1],worstder[2]],'right':dr[worstder[0],0,worstder[1],worstder[2]]},
             'rest_magnitudes_outside_native_valid_minus30_minus10':int(np.sum((values[:,0]<=-30)|(values[:,0]>=-10))),
             'half_decade_step_limit_jump_change_max':float(np.max(abs(jumps[::2]-jumps[1::2]))),
             'seconds':time.monotonic()-start,
             'interpretation':'Piecewise bilinear magnitude/error component only. Does not establish observer-flux continuity, covariance-state closure or any likelihood minimum.'}
    np.savez_compressed(OUT/'grid-scan.npz',phases=phases,shape=shape,boundaries=bounds,
                        limits_left=stack('limit_left'),limits_right=stack('limit_right'),jumps=jumps,derivatives_left=dl,derivatives_right=dr)
    save(OUT/'grid-result.json',results)
    print(json.dumps(results,indent=2))

if __name__=='__main__':main()
