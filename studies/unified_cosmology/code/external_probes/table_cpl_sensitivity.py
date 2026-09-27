"""Conditional actual-CMB data scan across CAMB's analytic CPL boundary.

Uses CAMB's documented tabulated-w interface, no global/source modification.
Fixed physical densities/H0/primordial spectrum/calibration/foregrounds: this is
not a marginalized posterior or a profile-likelihood exclusion of a domain.
"""
import os
os.environ.setdefault('CLIPY_NOJAX','1')
import json
import time
import sys
import ctypes
from contextlib import contextmanager
import numpy as np
import pandas as pd
import camb
from cobaya.model import get_model
from scipy.linalg import solve_triangular
from adapter import external_info,reference_point,PACKAGES,WORK
from acquire import sha,HERE,RESULTS
from validate import vector


@contextmanager
def capture_native_output(path):
    sys.stdout.flush();sys.stderr.flush();ctypes.CDLL(None).fflush(None)
    saved=[os.dup(i) for i in (1,2)]
    with path.open('w') as f:
        try:
            os.dup2(f.fileno(),1);os.dup2(f.fileno(),2)
            yield
        finally:
            sys.stdout.flush();sys.stderr.flush();ctypes.CDLL(None).fflush(None)
            for i,s in zip((1,2),saved):os.dup2(s,i);os.close(s)


def main():
    m=get_model(external_info('full','cpl'))
    m.add_requirements({'CAMBdata':None})
    point=reference_point(m)
    baseline=m.logposterior(point)
    nuisance_inputs=m.parameterization.to_input(point)
    pars=m.provider.get_CAMBdata().Params.copy()
    full=m.likelihood['planck_2018_highl_plik.TTTEEE']
    import clipy
    clik=clipy.clik(str(PACKAGES/'data/planck_2018/baseline/plc_3.0/hi_l/plik/plik_rd12_HM_v22b_TTTEEE.clik'))
    base=PACKAGES/'data/bao_data/desi_bao_dr2'
    data=pd.read_csv(base/'desi_gaussian_bao_ALL_GCcomb_mean.txt',sep=r'\s+',comment='#',names=['z','value','kind'])
    L=np.linalg.cholesky(np.loadtxt(base/'desi_gaussian_bao_ALL_GCcomb_cov.txt'))

    def evaluate(r):
        ps=r.get_cmb_power_spectra(CMB_unit='muK')['total']
        ell=np.arange(len(ps)); factor=ell*(ell+1)/(2*np.pi);factor[:2]=1
        dls={k:ps[:,j] for k,j in [('tt',0),('ee',1),('bb',2),('te',3)]}
        dls['pp']=r.get_lens_potential_cls(lmax=len(ps)-1)[:,0]
        cls={k:v/factor for k,v in dls.items() if k!='pp'}
        cls['pp']=dls['pp']/(factor**2*2*np.pi)
        likes={'highl':float(clik(vector(clik,cls,nuisance_inputs)))}
        for k in ['TT','EE']:
            likes[k]=float(m.likelihood['planck_2018_lowl.'+k].log_likelihood(dls[k.lower()],point['A_planck']))
        likes['lensing']=float(m.likelihood['planck_2018_lensing.native'].log_likelihood(dls,A_planck=point['A_planck']))
        rd=r.get_derived_params()['rdrag']; pred=[]
        for row in data.itertuples():
            dm=(1+row.z)*r.angular_diameter_distance(row.z);dh=299792.458/r.hubble_parameter(row.z)
            pred.append({'DM_over_rs':dm/rd,'DH_over_rs':dh/rd,'DV_over_rs':(row.z*dm**2*dh)**(1/3)/rd}[row.kind])
        q=solve_triangular(L,data.value.to_numpy()-pred,lower=True)
        likes['BAO']=float(-.5*q@q)
        return likes

    original=evaluate(m.provider.get_CAMBdata())
    expected=dict(zip(['highl','TT','EE','lensing','BAO'],map(float,baseline.loglikes)))
    assert max(abs(original[k]-expected[k]) for k in original)<1e-6, (original,expected)
    rows=[]
    cases=[(-1.,0.,4000)]
    cases += [(w,early-w,n) for w in [-.8,-1.,-2.]
              for early in [-.02,0.,.02] for n in ([2000,4000] if early==.02 else [4000])]
    for w,wa,n in cases:
        start=time.monotonic(); p=pars.copy();a=np.geomspace(1e-9,1,n)
        row={'w':w,'wa':wa,'w_plus_wa':w+wa,'table_nodes':n}
        try:
            p.DarkEnergy.set_w_a_table(a,w+wa*(1-a))
            native_log=WORK/f'table-cpl-native-{w}-{wa}-{n}.log'
            with capture_native_output(native_log):
                r=camb.get_results(p); likes=evaluate(r)
            warning_lines=[s for s in native_log.read_text().splitlines() if 'WARNING' in s.upper()]
            row.update(status='finite_theory',rdrag_Mpc=float(r.get_derived_params()['rdrag']),
                native_warning_count=len(warning_lines),
                native_warning_messages=sorted(set(warning_lines)),
                numerical_warning_free=not warning_lines,
                OmegaDE_z1100=float(r.get_Omega('de',1100.)),
                loglikes={k:v if np.isfinite(v) else None for k,v in likes.items()},
                nonfinite_likelihoods=[k for k,v in likes.items() if not np.isfinite(v)],
                delta_chi2_from_reference_by_component={k:float(-2*(v-original[k])) if np.isfinite(v) else None for k,v in likes.items()})
        except Exception as e:
            row.update(status='theory_or_likelihood_error',exception=type(e).__name__,message=str(e))
        row['seconds']=time.monotonic()-start;rows.append(row)
        print(json.dumps(row),flush=True)
    out={'status':'completed_conditional_scan','reference_point':{k:float(v) for k,v in point.items()},
         'reference_loglikes':original,'direct_reconstruction_max_difference':max(abs(original[k]-expected[k]) for k in original),
         'rows':rows,'active_baseline_unchanged':True,
         'scope':'Finite fixed-parameter actual spectra/lensing/BAO likelihood evaluations across the analytic w+wa=0 implementation boundary. This does not optimize other cosmological or nuisance parameters and cannot establish integrated posterior probability outside the boundary. Tabulated PPF is a phenomenological extrapolation, not a new physical dark-energy model.',
         'code_sha256':sha(__file__),'adapter_sha256':sha(HERE/'adapter.py')}
    (RESULTS/'table-cpl-sensitivity.json').write_text(json.dumps(out,indent=2)+'\n')

if __name__=='__main__':main()
