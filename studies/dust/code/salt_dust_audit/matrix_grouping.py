#!/usr/bin/env python3
"""Resolve whether released systematic groups overlap; do not sum labels blindly."""
import json
import numpy as np
from scipy.linalg import eigh
import matrix_audit as audit


def spectrum(k,sig):
    kw=k/sig[:,None]/sig[None,:]
    ev=eigh(kw,eigvals_only=True)
    tol=max(1e-12,max(abs(ev))*1e-5)
    return {'stat_whitened_min_eigenvalue':float(ev[0]),'stat_whitened_max_eigenvalue':float(ev[-1]),
            'positive_rank_at_relative_1e_minus5':int(sum(ev>tol)),
            'negative_rank_at_relative_1e_minus5':int(sum(ev< -tol)),
            'stat_whitened_frobenius':float(np.linalg.norm(ev)),
            'stat_whitened_eigenvalues':ev.tolist()}


def main():
    out=audit.ROOT/'runs/salt_dust_audit/matrix_grouping'
    out.mkdir(parents=True,exist_ok=True)
    result={}
    for label,cfg in audit.RELEASES.items():
        base=audit.ROOT/cfg['root'];hd=audit.read_table(base/cfg['hd']);n=len(hd)
        audit.INPUTS.add(base.parent/'7_PIPPIN_FILES/D5yr_analysis.yml')
        suffix='.txt.gz' if label=='original' else '.npz'
        stat=audit.load_matrix(base/('STATONLY'+suffix),n)
        total=audit.load_matrix(base/('STAT+SYS'+suffix),n)
        if label=='original':
            stat+=np.diag(hd.MUERR_FINAL.to_numpy()**2)
        else:
            total-=stat
        sig=np.sqrt(stat.diagonal())
        summ=np.zeros((n,n));cal=None;calspec=None
        for p in sorted((base/'SingleSYS_CovMatrix').glob('*'+suffix)):
            k=audit.load_matrix(p,n)
            if label=='Dovekie':k-=stat
            summ+=k
            if p.name in ['CALIBplusSALT3.txt.gz','CAL_SALT3.npz']:cal=k
            if p.name.startswith('CALSPEC.'):calspec=k
        r=total-(summ-calspec)
        rr=spectrum(r,sig)
        totalnorm=np.linalg.norm(total/sig[:,None]/sig[None,:])
        result[label]={
            'calibration_group_minus_CALSPEC':spectrum(cal-calspec,sig),
            'total_sys_minus_sum_except_separate_CALSPEC':rr,
            'residual_relative_total_sys_stat_whitened_frobenius':rr['stat_whitened_frobenius']/totalnorm,
            'residual_trace_mag2':float(np.trace(r)),
            'interpretation':'CALSPEC subtraction is a covariance grouping diagnostic, not evidence of independent physical priors or a correction to data.'}
        print(label,'relative corrected sum residual',result[label]['residual_relative_total_sys_stat_whitened_frobenius'],flush=True)
    from pathlib import Path
    audit.INPUTS.add(Path(__file__).resolve())
    audit.INPUTS.add(Path(audit.__file__).resolve())
    audit.INPUTS.add(audit.ROOT/'sources/repos/RickKessler__SNANA/util/create_covariance.py')
    result['input_sha256']={str(p.relative_to(audit.ROOT)):audit.sha(p) for p in sorted(audit.INPUTS)}
    (out/'results.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
