#!/usr/bin/env python3
"""Inspect published SNANA GENPDF distributions and extinction responses.

This integrates probability density maps, using the multilinear map interpolation
of SNANA, before survey selection. It does not regenerate BBC or fit light curves.
"""
from __future__ import annotations
import csv
import gzip
import hashlib
import json
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import brentq
from snana_extinction import ROOT, SnanaExtinction

INPUT = ROOT/'phase2/official/inputs/SNDATA_ROOT/models/population_pdf/DES-SN5YR'
OUT = ROOT/'runs/salt_dust_audit/snana_population'
WAVES=np.array([3500.,4400.,5500.,6500.,8000.])


def parse(path):
    txt=gzip.open(path,'rt').read() if path.suffix=='.gz' else path.read_text()
    maps={}; cols=None; config={}
    for line in txt.splitlines():
        if line.startswith('VARNAMES:'):
            cols=line.split()[1:];maps[cols[0]]={'columns':cols,'rows':[]}
        elif line.startswith('PDF:'):
            maps[cols[0]]['rows'].append([float(x) for x in line.split()[1:]])
        elif line.startswith(('GENPEAK_','GENSIGMA_','GENRANGE_','MAG_OFFSET:')):
            key, val=line.split(':',1);config[key]=val.strip()
    for name, m in maps.items():
        a=np.array(m.pop('rows'));axes=[np.unique(a[:,j]) for j in range(a.shape[1]-1)]
        vals=np.full(tuple(len(x) for x in axes),np.nan)
        indexes=tuple(np.searchsorted(axes[j],a[:,j]) for j in range(len(axes)))
        vals[indexes]=a[:,-1]
        assert np.all(np.isfinite(vals)), name
        m.update(axes=axes, values=vals, interp=RegularGridInterpolator(axes,vals,bounds_error=True))
    return maps,config


def quadrature(m, params, n=4, breaks=()):
    # Interval Gauss-Legendre exactly integrates linearly interpolated pdf moments
    # through degree 6 at n=4. Split intervals at thresholds before tail integrals.
    x=m['axes'][0]
    edges=np.unique(np.r_[x,[b for b in breaks if x[0]<b<x[-1]]])
    t,w=leggauss(n);mid=(edges[1:]+edges[:-1])/2;half=np.diff(edges)/2
    q=(mid[:,None]+half[:,None]*t).ravel(); qw=(half[:,None]*w).ravel()
    arr=np.column_stack([q]+[np.full(len(q),params[k]) for k in m['columns'][1:-1]])
    density=m['interp'](arr);pw=density*qw;norm=pw.sum();assert norm>0
    keep=pw>0
    return q[keep],pw[keep]/norm


def write_csv(path, rows):
    with path.open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    ext=SnanaExtinction(ROOT/'runs/salt_dust_audit/snana_extinction/libsnana_extinction.so')
    paths=sorted(INPUT.glob('DES-SN5YR_*_S3_P21*.DAT.gz'))+sorted(INPUT.glob('DES-SN5YR_*_S3_BS21.DAT.gz'))
    rows=[];responses=[];maps_meta=[];configrows=[]
    zero_roots=[]
    for wave in WAVES:
        for option in [99,-99]:
            f=lambda rv:float(ext(wave,rv,option=option))
            if f(.01)*f(8)<0:
                zero_roots.append({'wave_A':float(wave),'option':option,'RV_zero':float(brentq(f,.01,8))})
    breaks=[1,1.4,2,6]+[r['RV_zero'] for r in zero_roots]
    # Exact monochromatic response evaluated at quadrature nodes and normalized.
    # Some nodes are RV<2; retain that extrapolation so the released grid is audited.
    cache={}
    for path in paths:
        maps,config=parse(path)
        parts=path.name.split('_');survey=parts[1];variant=parts[3].split('.')[0]
        maps_meta.append({'file':str(path.relative_to(ROOT)),'axes':{k:{col:axis.tolist() for col,axis in zip(m['columns'][:-1],m['axes'])} for k,m in maps.items()}, 'config':config})
        configrows.append(dict(survey=survey,variant=variant,**config))
        if any('COLOR' in m['columns'] for m in maps.values()):
            maps_meta[-1]['integration_status']='Not integrated: conditions on host COLOR rather than mass; no fabricated color-to-mass mapping.'
            continue
        for mass in [9.,9.5,10.,11.]:
            for z in [.05,.1,.5,1.]:
                pars={'LOGMASS':mass,'ZTRUE':z}
                rv,p=quadrature(maps['RV'],pars,breaks=breaks)
                ebv,pe=quadrature(maps['EBV'],pars,breaks=[.3])
                mean_e=float(ebv@pe);mean_rv=float(rv@p)
                row={'survey':survey,'variant':variant,'logmass':mass,'z':z,'mean_RV':mean_rv,'sd_RV':float(np.sqrt((rv**2)@p-mean_rv**2)),'P_RV_lt1':float(p[rv<1].sum()),'P_RV_lt1p4':float(p[rv<1.4].sum()),'P_RV_lt2':float(p[rv<2].sum()),'P_RV_gt6':float(p[rv>6].sum()),'mean_EBV':mean_e,'P_EBV_gt0p3':float(pe[ebv>.3].sum()),'mean_AV_independent_conditional':mean_e*mean_rv}
                if 'SALT2c' in maps:
                    c,pc=quadrature(maps['SALT2c'],pars); cm=float(c@pc)
                    row['mean_c_int']=cm;row['sd_c_int']=float(np.sqrt(c**2@pc-cm**2))
                else:row['mean_c_int']=row['sd_c_int']=None
                key=tuple(rv)
                if key not in cache:
                    cache[key]=(np.stack([ext(WAVES,float(r)) for r in rv]),np.stack([ext(WAVES,float(r),option=-99) for r in rv]))
                exact,approx=cache[key]
                neg_exact=(exact<0).any(axis=1);neg_approx=(approx<0).any(axis=1)
                row['P_negative_A_any_5_waves_exact']=float(p[neg_exact].sum())
                row['P_negative_A_any_5_waves_approx']=float(p[neg_approx].sum())
                rows.append(row)
                for j,wave in enumerate(WAVES):
                    responses.append({'survey':survey,'variant':variant,'logmass':mass,'z':z,'wavelength_A':wave,'mean_A_exact_mag':float(mean_e*(p@exact[:,j])),'mean_A_approx_mag':float(mean_e*(p@approx[:,j])),'mean_exact_minus_approx_mag':float(mean_e*(p@(exact[:,j]-approx[:,j]))),'P_A_exact_lt0':float(p[exact[:,j]<0].sum()),'P_A_approx_lt0':float(p[approx[:,j]<0].sum())})
    write_csv(OUT/'conditional_moments.csv',rows);write_csv(OUT/'conditional_extinction_response.csv',responses)
    (OUT/'grid_axes_and_config.json').write_text(json.dumps(maps_meta,indent=2)+'\n')
    # Scenario response matrix: each row is a nominal DES (mass,ZTRUE,wavelength)
    # point, columns change one documented assumption. Not a posterior covariance.
    lookup={(r['variant'],r['logmass'],r['z'],r['wavelength_A']):r for r in responses if r['survey']=='DES'}
    scenario_rows=[]
    for row in responses:
        if row['survey']!='DES' or row['variant']!='P21':continue
        m,z,w=row['logmass'],row['z'],row['wavelength_A'];nom=row['mean_A_exact_mag']
        out={'logmass':m,'ZTRUE':z,'wavelength_A':w,'F99_exact_minus_approx':row['mean_exact_minus_approx_mag']}
        for v in ['P21sys1','P21sys2','P21sys3','BS21']:
            out[v+'_minus_P21']=lookup[(v,m,z,w)]['mean_A_exact_mag']-nom
        scenario_rows.append(out)
    write_csv(OUT/'scenario_response_matrix.csv',scenario_rows)
    columns=list(scenario_rows[0])[3:];R=np.array([[r[c]for c in columns]for r in scenario_rows])
    norms=np.linalg.norm(R,axis=0);unit=R/norms
    gram=unit.T@unit
    write_csv(OUT/'scenario_cosine_matrix.csv',[dict(column=c,**dict(zip(columns,row)))for c,row in zip(columns,gram)])
    U,sv,Vt=np.linalg.svd(unit,full_matrices=False)
    (OUT/'scenario_geometry.json').write_text(json.dumps({'columns':columns,'interpretation':'Euclidean cosine and SVD of unweighted parent-grid extinction scenario vectors. These are collinearities of assumption responses, not posterior parameter correlations or an uncertainty covariance. Systematic scenarios are not independent Gaussian draws.', 'unit_column_singular_values':sv.tolist(),'unit_column_condition_number':float(sv[0]/sv[-1]),'right_singular_vectors':Vt.tolist(),'norms_mag':norms.tolist()},indent=2)+'\n')
    # Numerical convergence check at nominal high-mass branch, where low RV matters.
    maps,_=parse(INPUT/'DES-SN5YR_DES_S3_P21.DAT.gz');pars={'LOGMASS':11.,'ZTRUE':.5}
    expectation=[]
    for n in [4,8,16]:
        rv,p=quadrature(maps['RV'],pars,n=n,breaks=breaks);e,pe=quadrature(maps['EBV'],pars,n=n)
        f=np.stack([ext(WAVES,float(r)) for r in rv]);expectation.append(float(e@pe)*(p@f))
    convergence=float(np.max(np.abs(expectation[0]-expectation[-1])))
    assert convergence<1e-8
    result={'classification':'Deterministic integrations of published parent simulation probability grids. No detection/cosmology selection and no BBC correction or empirical confidence interval.', 'quadrature':'4-point Gauss-Legendre within every PDF segment; multilinear host/redshift interpolation; normalize conditional PDF after interpolation. RV threshold tails integrated after splitting at threshold.', 'conditional_independence':'Separate RV and EBV grids imply independent draws at fixed mass/redshift for these inputs, so mean A = mean EBV * mean[A(RV,EBV=1)].', 'selection_warning':'These parent means/tails cannot be substituted for selected-sample means or cosmological corrections.', 'extinction_zero_roots':zero_roots, 'negative_extinction_warning':'P_negative_A tests only the 5 recorded monochromatic wavelengths. This diagnoses F99/approx-law extrapolation into low-RV support; it is not a claim about negative broadband extinction after fitting.', 'convergence_max_abs_mag_4_vs16':convergence,'nominal_DES_examples':[r for r in rows if r['survey']=='DES' and r['variant']=='P21' and r['z'] in [.05,.5]], 'files_hashed':[{'path':str(p.relative_to(ROOT)), 'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}for p in paths]}
    (OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    files=[Path(__file__).resolve(),ROOT/'scripts/salt_dust_audit/snana_extinction.py',ROOT/'sources/repos/RickKessler__SNANA/src/sntools_genPDF.c',ROOT/'sources/repos/RickKessler__SNANA/src/snlc_sim.c']+paths+sorted(OUT.glob('*.csv'))+[OUT/'results.json',OUT/'grid_axes_and_config.json',OUT/'scenario_geometry.json']
    (OUT/'manifest.json').write_text(json.dumps({'files':[{'path':str(p.relative_to(ROOT)), 'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size}for p in files]},indent=2)+'\n')
    print(json.dumps({'convergence_max_abs_mag':convergence,'nominal_DES_examples':result['nominal_DES_examples']},indent=2))

if __name__=='__main__':main()
