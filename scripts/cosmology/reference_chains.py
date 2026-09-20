"""Derived quantities from ORIGINAL DESI chains; not a likelihood reproduction."""
import json
import numpy as np
import pandas as pd
from core import ROOT, manifest, qvalue


def stats(x,w):
    idx=np.argsort(x); xx=x[idx]; ww=w[idx]; c=(np.cumsum(ww)-.5*ww)/ww.sum()
    m=np.average(x,weights=w)
    return dict(mean=float(m),sd=float(np.sqrt(np.average((x-m)**2,weights=w))),
                **{str(p):float(np.interp(p,c,xx)) for p in [.025,.16,.5,.84,.975]})


def main():
    out=ROOT/'runs/cosmology/desi-reference-chains';out.mkdir(parents=True,exist_ok=True)
    inputs=[];result={}; histories=[]
    for d in sorted((ROOT/'data/bao/desi-dr2-reference/cobaya/base_w_wa').iterdir()):
        label='BAO+CMB'+('+Pantheon+' if 'pantheonplus' in d.name else '+DESY5' if 'desy5sn' in d.name else '')
        result[label]={}
        files=sorted(d.glob('chain.[1-4].txt'));inputs+=files
        if not files: continue
        for burn in [0.,.1,.3,.5]:
            parts=[];perchain=[]
            for p in files:
                columns=p.open().readline().strip().lstrip('#').split()
                a=np.loadtxt(p);a=a[int(len(a)*burn):]
                df=pd.DataFrame(a,columns=columns)
                parts.append(df)
                q=.5+1.5*df.w.to_numpy()*(1-df.omegam.to_numpy())
                perchain.append({'chain':p.name,'q0_mean':float(np.average(q,weights=df.weight)),
                                 'P_q0_negative':float(np.average(q<0,weights=df.weight))})
            df=pd.concat(parts,ignore_index=True);w=df.weight.to_numpy()
            q=.5+1.5*df.w.to_numpy()*(1-df.omegam.to_numpy())
            result[label][f'row_burn_{burn}']=dict(rows=len(df),sum_weights=float(w.sum()),
                q0=stats(q,w),P_q0_negative=float(np.average(q<0,weights=w)),
                w0=stats(df.w.to_numpy(),w),wa=stats(df.wa.to_numpy(),w),
                Om=stats(df.omegam.to_numpy(),w),per_chain=perchain)
            if burn==.3:
                for z in np.linspace(0,2,101):
                    qz=qvalue([z],df[['omegam','w','wa']].to_numpy())[:,0]
                    s=stats(qz,w)
                    histories.append({'sample':label,'z':z,'q16':s['0.16'],'median':s['0.5'],
                                      'q84':s['0.84'],'P_q_negative':float(np.average(qz<0,weights=w))})
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    pd.DataFrame(histories).to_csv(out/'q_histories.csv',index=False)
    manifest(out,'Derived q(z) from published DESI samples, not independent refit',inputs,
             {'burn_in_row_fractions':[0,.1,.3,.5],'display_burn':.3,'q_approximation':'flat low-z matter+DE; radiation neglected in q calculation'},
             [out/'summary.json',out/'q_histories.csv'])
    for label,r in result.items():
        x=r['row_burn_0.3'];print(label,x['q0'],x['P_q0_negative'])

if __name__=='__main__':main()
