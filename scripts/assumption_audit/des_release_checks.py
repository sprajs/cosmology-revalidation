#!/usr/bin/env python3
"""Focused release consistency checks; never changes upstream artifacts."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.linalg import eigh
from scipy.optimize import minimize_scalar
from astropy.table import Table
from des_covariance import BASE, ROOT, OUT, precision, inverse, projected, mu


def main():
    used = [BASE/'STATONLY.npz', BASE/'STAT+SYS.npz']
    cs = inverse(precision(used[0]))
    ct = inverse(precision(used[1]))
    w = 1/np.sqrt(np.diag(cs))
    increments = {}
    for f in sorted((BASE/'SingleSYS_CovMatrix').glob('*.npz')):
        used.append(f)
        increments[f.stem] = inverse(precision(f))-cs
    summed = sum(increments.values())
    residual = ct-cs-summed+increments['CALSPEC']
    purecal = increments['CAL_SALT3']-increments['CALSPEC']
    evcal = eigh(w[:,None]*purecal*w[None,:],eigvals_only=True)
    fullsys = ct-cs
    closure = {'identity':'Ctotal-Cstat = sum(single_i-Cstat) - (C_CALSPEC-Cstat)',
        'max_abs_remainder_mag2':float(np.max(abs(residual))),
        'relative_stat_whitened_frobenius_remainder':float(np.linalg.norm(w[:,None]*residual*w[None,:])/np.linalg.norm(w[:,None]*fullsys*w[None,:])),
        'CAL_SALT3_minus_CALSPEC_positive_rank_at_2e_6':int(np.sum(evcal>2e-6)),
        'CAL_SALT3_minus_CALSPEC_top_eigenvalues':evcal[-11:].tolist()}
    # Quantisation is diagnosed from the released covariance itself; we do not
    # attribute its provenance to a particular unreleased generation step.
    v = increments['VPEC']
    qv = np.round(v, 5)
    ev, vec = eigh(v)
    evw = eigh(w[:,None]*v*w[None,:],eigvals_only=True)
    diag = np.maximum(0,np.diag(qv))
    cauchy_excess = np.abs(qv)-np.sqrt(diag[:,None]*diag[None,:])
    np.fill_diagonal(cauchy_excess,0)
    i,j = np.unravel_index(np.argmax(cauchy_excess),v.shape)
    hdpath = BASE/'DES-Dovekie_HD.csv'; used.append(hdpath)
    d = pd.read_csv(hdpath,sep=r'\s+',comment='#',dtype={'CID':str})
    selected = [i,j]
    vpec = {'min_eigenvalue_mag2':float(ev[0]),'min_stat_whitened_eigenvalue':float(evw[0]),
        'negative_stat_whitened_eigenvalues_below_minus_2e_6':int(np.sum(evw< -2e-6)),
        'max_abs_deviation_from_5_decimal_quantisation_mag2':float(np.max(abs(v-qv))),
        'violating_pair':[{'CID':str(d.CID.iloc[k]),'IDSURVEY':int(d.IDSURVEY.iloc[k]),'zHD':float(d.zHD.iloc[k])} for k in selected],
        'violating_pair_matrix_mag2':qv[np.ix_(selected,selected)].tolist(),
        'violating_pair_Cauchy_excess_mag2':float(cauchy_excess[i,j]),
        'violating_pair_min_eigenvalue_using_unrounded_release':float(eigh(v[np.ix_(selected,selected)],eigvals_only=True)[0]),
        'total_min_stat_whitened_eigenvalue':float(eigh(w[:,None]*ct*w[None,:],eigvals_only=True,subset_by_index=[0,0])[0])}
    # This arbitrary PSD completion quantifies scale; it is not a replacement
    # uncertainty model or a recommended correction to the published data.
    repair = (vec[:,ev<0]*(-ev[ev<0]))@vec[:,ev<0].T
    fitrows=[]
    z,zh,y = d.zHD.to_numpy(),d.zHEL.to_numpy(),d.MU.to_numpy()
    for label,cov in [('published_total',ct),('illustrative_VPEC_negative_eigenvalues_clipped',ct+repair)]:
        pr=projected(inverse(cov))
        fit=minimize_scalar(lambda om:(y-mu(om,z,zh))@pr@(y-mu(om,z,zh)),bounds=(.1,.6),method='bounded',options={'xatol':1e-11})
        deriv=(mu(fit.x+1e-5,z,zh)-mu(fit.x-1e-5,z,zh))/2e-5
        fitrows.append({'label':label,'Omega_m':float(fit.x),'local_Fisher_sigma_Omega_m':float(1/np.sqrt(deriv@pr@deriv))})
    vpec['illustrative_PSD_sensitivity_not_a_physical_refit']=fitrows
    # Build table of published prior draws; no light curves are re-fit.
    rows=[]
    sysbase = BASE.parent/'2_LCFIT_MODEL/SALT3.DOVEKIE-SYS'
    for f in sorted(sysbase.glob('*/SALT3.INFO')):
        used.append(f)
        for lineno,line in enumerate(f.read_text().splitlines(),start=1):
            parts=line.split()
            if parts and parts[0] in ['WAVESHIFT:','MAGSHIFT:']:
                rows.append({'realisation':f.parent.name,'type':parts[0][:-1],
                    'survey':parts[1],'filter':parts[2],'value':float(parts[3]),
                    'source':str(f.relative_to(ROOT)),'line':lineno})
    draw=pd.DataFrame(rows)
    draw.to_csv(OUT/'calibration_prior_draws.csv',index=False)
    pairs=[]
    for filt in ['PS1-g','PS1-r','PS1-i','PS1-z']:
        s=draw[(draw.type=='WAVESHIFT')&(draw['filter']==filt)].pivot(index='realisation',columns='survey',values='value')
        a,b=s.PS1MD.to_numpy(),s.FOUNDATION.to_numpy()
        pairs.append({'filter':filt,'n':len(s),'number_different':int(np.sum(a!=b)),
            'empirical_correlation_across_released_draws':float(np.corrcoef(a,b)[0,1]),
            'rms_difference_Angstrom':float(np.sqrt(np.mean((a-b)**2))),
            'first_PS1MD_Angstrom':float(a[0]),'first_Foundation_Angstrom':float(b[0])})
    # Reproduce exactly the table reader used in the shipped likelihood.
    parsed = Table.read(hdpath,format='ascii.csv')
    reader = {'shipped_format':'ascii.csv','parsed_columns':parsed.colnames,
        'zHD_available': 'zHD' in parsed.colnames,
        'default_data_path_exists':(BASE/'DES-SN5YR_HD.csv').exists(),
        'default_covmat_path_exists':(BASE/'STAT+SYS.txt.gz').exists()}
    result={'calibration_overlap':closure,'vpec_increment':vpec,'PS1_Foundation_wavelength_draws':pairs,
        'shipped_likelihood_reader':reader,
        'finite_draw_projected_variance_fractional_SD_iid_Gaussian':{'9':float(np.sqrt(2/9)),'10':float(np.sqrt(2/10)),'3':float(np.sqrt(2/3))}}
    (OUT/'release_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
    for f in ['sources/repos/bap37__Dovekie/surfaces-dovekie.py','sources/repos/bap37__Dovekie/DOVEKIE_DEFS.yml',
        'sources/repos/bap37__Dovekie/templates/new_kcor_templates/Foundation.input','sources/repos/bap37__Dovekie/templates/new_kcor_templates/PS1SN.input',
        'sources/repos/des-science__DES-SN5YR/4_DISTANCES_COVMAT/DES-Dovekie-SN_Likelihood.py']:
        used.append(ROOT/f)
    used.append(ROOT/'scripts/assumption_audit/des_covariance.py')
    codeout=OUT/'code'
    codeout.mkdir(exist_ok=True)
    (codeout/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    (OUT/'release_checks_manifest.json').write_text(json.dumps({'script_sha256':sha(Path(__file__)),
        'code_snapshot':str((codeout/Path(__file__).name).relative_to(ROOT)),
        'inputs':{str(f.relative_to(ROOT)):sha(f) for f in used},
        'outputs':{str((OUT/f).relative_to(ROOT)):sha(OUT/f) for f in ['release_checks.json','calibration_prior_draws.csv']}},indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__': main()
