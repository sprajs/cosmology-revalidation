"""Check Pantheon+ grouped systematics against the final statistical baseline.

For the two one-variation MW products, recover the rank-one systematic from
an isolated high-z reference row. Its sign is arbitrary; only covariance is
used. Agreement of separately recovered baselines is the independent check.
"""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
import pandas as pd
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize_scalar

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts/cosmology'))
from core import mu


def digest(p):
    with open(p,'rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def load(p):
    a=np.loadtxt(p);n=int(a[0]);assert a.size==1+n*n
    c=a[1:].reshape(n,n);return (c+c.T)/2


def main():
    base=ROOT/'sources/repos/PantheonPlusSH0ES__DataRelease/Pantheon+_Data/4_DISTANCES_AND_COVAR'
    p=base/'Pantheon+SH0ES.dat';stat_path=base/'Pantheon+SH0ES_STATONLY.cov'
    f=pd.read_csv(p,sep=r'\s+');s=load(stat_path)
    keep=f.zHD.to_numpy()>.01
    unique=~f.CID.duplicated(keep=False).to_numpy()
    inputs=[p,stat_path,ROOT/'scripts/cosmology/core.py',Path(__file__)]
    result={'classification':'Confirmed product baseline mismatch; main analysis uses the final total covariance and is unaffected.',
            'groups':{}, 'fits':{}}
    reconstructed={};variants={'final_statonly':s}
    for label in ['MWEBV','MWCOLORLAW']:
        path=base/f'sytematic_groupings/Pantheon+SH0ES_122221_{label}.cov';inputs.append(path)
        w=load(path);d=w-s
        valid=(f.zHD.to_numpy()>.1)&unique
        k=int(np.argmax(np.where(valid,np.diag(d),-np.inf)))
        v=d[k,:]/np.sqrt(d[k,k]);b=w-np.outer(v,v)
        reconstructed[label]=b
        delta=b-s
        bad=abs(delta)>1e-6
        same=f.CID.to_numpy()[:,None]==f.CID.to_numpy()[None,:]
        result['groups'][label]={
            'reference_row':k,'reference_CID':str(f.CID.iloc[k]),
            'maximum_statistical_baseline_difference_mag2':float(abs(delta).max()),
            'diagonal_entries_differing_above_1e_6':int(np.count_nonzero(np.diag(bad))),
            'offdiagonal_entries_differing_above_1e_6':int(np.count_nonzero(bad)-np.count_nonzero(np.diag(bad))),
            'entries_between_different_literal_CIDs_above_1e_6':int(np.count_nonzero(bad&~same)),
            'selected_rows_touched':int(np.count_nonzero(np.any(bad[np.ix_(keep,keep)],axis=1))),
            'examples':[]}
        for i,j in [(47,48),(48,49)]:
            result['groups'][label]['examples'].append({'rows':[i,j],'CIDs':[str(f.CID.iloc[i]),str(f.CID.iloc[j])],
                'final_stat':float(s[i,j]),'group_total':float(w[i,j]),'recovered_group_stat':float(b[i,j])})
        variants[label+'_as_shipped']=w
        variants[label+'_on_final_statistical_baseline']=s+np.outer(v,v)
    residual=reconstructed['MWEBV']-reconstructed['MWCOLORLAW']
    result['baseline_crosscheck_max_mag2']=float(abs(residual).max())
    assert abs(residual).max()<2e-7, 'The one-mode recovery must leave the same baseline in both products'
    # Differing entries beyond rounding provide an inspectable row-pair record.
    delta=reconstructed['MWEBV']-s
    ii,jj=np.where(np.triu(abs(delta)>1e-6))
    table=pd.DataFrame({'i':ii,'j':jj,'CID_i':f.CID.to_numpy()[ii],'CID_j':f.CID.to_numpy()[jj],
                        'final_stat':s[ii,jj],'group_stat':reconstructed['MWEBV'][ii,jj],
                        'difference':delta[ii,jj]})
    for label,c in variants.items():
        cc=c[np.ix_(keep,keep)]
        inv=cho_solve(cho_factor(cc,lower=True),np.eye(keep.sum()))
        u=inv.sum(axis=1);a=inv-np.outer(u,u)/u.sum()
        def fun(om):
            r=f.loc[keep,'m_b_corr'].to_numpy()-mu(f.loc[keep,'zHD'].to_numpy(),[om],'lcdm',f.loc[keep,'zHEL'].to_numpy())[0]
            return float(r@a@r)
        opt=minimize_scalar(fun,bounds=(.05,.7),method='bounded',options={'xatol':1e-11})
        h=1e-4;curv=(fun(opt.x+h)-2*fun(opt.x)+fun(opt.x-h))/h**2
        assert opt.success and curv>0
        result['fits'][label]={'Om_profile':float(opt.x),'local_sigma':float(np.sqrt(2/curv)),
                              'chisq':float(opt.fun),'rows':int(keep.sum())}
    result['interpretation']=[
        'Systematic group files cannot be treated as final STATONLY plus the named systematic.',
        'Subtracting final STATONLY produces a non-PSD difference because the statistical baselines differ; it is not evidence of negative physical systematic variance.',
        'Two independent rank-one reconstructions leave a common background to printing precision.',
        'Rank-one sign is not determined; these are covariance-only comparisons, not signed dust corrections.',
        'The grouped covariance has older/mismatched statistical terms; the tests do not identify which published author results used those files.',
        'The profile fits isolate this release-file trap; they are not full STAT+SYS cosmology estimates.',
    ]
    out=ROOT/'runs/assumption_audit';out.mkdir(parents=True,exist_ok=True)
    dest=out/'covariance-group-differences.csv';table.to_csv(dest,index=False)
    target=out/'covariance-groups.json';target.write_text(json.dumps(result,indent=2)+'\n')
    (out/'covariance-groups-manifest.json').write_text(json.dumps({
        'inputs_sha256':{str(p.relative_to(ROOT)):digest(p) for p in inputs},
        'outputs_sha256':{str(p.relative_to(ROOT)):digest(p) for p in [target,dest]},
    },indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
