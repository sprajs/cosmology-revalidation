"""Initial 13D modern-CMB proposal only, not a scientific prior."""
import json
import numpy as np
from adapter import WORK,ROOT
from acquire import sha,RESULTS

def main():
    rows={}
    for model in ['lcdm','cpl']:
        source=WORK/f'proposal-{model}-dovekie.covmat'
        names=source.read_text().splitlines()[0].lstrip('# ').split()
        c=np.loadtxt(source)
        use=['H0','ombh2','omch2','logA','ns','tau']+(['w','wa'] if model=='cpl' else [])+['A_planck']
        ind=[names.index(k) for k in use]
        extras={'P_act':.003,'Tcal':.0036,'Ecal':.0095,'A_fg':.3}
        out=np.zeros((len(use)+len(extras),)*2)
        out[:len(use),:len(use)]=c[np.ix_(ind,ind)]
        for i,sd in enumerate(extras.values(),len(use)):out[i,i]=sd*sd
        np.linalg.cholesky(out)
        dest=WORK/f'proposal-modern-{model}-dovekie.covmat'
        np.savetxt(dest,out,header=' '.join(use+list(extras)))
        rows[model]={'path':str(dest.relative_to(ROOT)),'sha256':sha(dest),
                     'source_sha256':sha(source),'dimension':len(out)}
    record={'models':rows,'status':'positive_definite',
         'scope':'Subset of transformed Planck+SN proposal retaining cosmology/A_planck correlations; independent initial scales for new calibration/foreground parameters. Does not use modern CMB posterior and does not change scientific prior.',
         'code_sha256':sha(__file__)}
    (RESULTS/'modern-proposal.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2))

if __name__=='__main__':main()
