"""Independent reading of output-only native states; no fitter imports.

This checks each recorded fixed-state quadratic. Multiplicative-mean probes
and post-initialization starts are separate gates and are not presumed passed.
"""
from pathlib import Path
import argparse,hashlib,json
import numpy as np

def parse(path):
    rows=[];state=None;records=[];entry=[]
    for line in path.read_text().splitlines():
        a=line.split()
        if not a:continue
        if a[0]=='CSP_ENTRY:':entry.append(a[1:])
        elif a[0]=='CSP_ROW:':
            cid,it,idx,epoch,band=a[1:6];v=list(map(float,a[6:]));assert len(v)==13
            rows.append(dict(CID=cid,iteration=int(it),index=int(idx),epoch=int(epoch),band=band,values=v))
        elif a[0]=='CSP_OBJECTIVE:':
            cid,it,n,cov=a[1:5];v=list(map(float,a[5:]));assert len(v)==10
            assert len(rows)==int(n),(path,cid,it,len(rows),n)
            assert [r['index'] for r in rows]==list(range(1,int(n)+1))
            assert all(r['CID']==cid and r['iteration']==int(it) for r in rows)
            state=dict(CID=cid,iteration=int(it),n=int(n),cov=cov=='T',objective=v,rows=rows,weights=[]);rows=[]
        elif a[0] in ['CSP_WROW:','CSP_WDIAG:']:
            assert state is not None
            cid,it,ix=a[1:4];assert cid==state['CID'] and int(it)==state['iteration']
            assert int(ix)==len(state['weights'])+1
            v=list(map(float,a[4:]));assert len(v)==(state['n'] if state['cov'] else 1)
            state['weights'].append(v)
            if len(state['weights'])==state['n']:records.append(state);state=None
    assert not rows and state is None
    return records,entry

def main():
    ap=argparse.ArgumentParser();ap.add_argument('log',type=Path);ap.add_argument('out',type=Path);args=ap.parse_args()
    args.out.mkdir(exist_ok=False,parents=True)
    records,entry=parse(args.log);assert records
    result=[]
    for record in records:
        arr=np.array([r['values'] for r in record['rows']]);m=arr[:,2];y=arr[:,4]
        W=np.array(record['weights']);W=W if record['cov'] else np.diag(W[:,0])
        symmetric=float(np.max(abs(W-W.T)));scale=float(np.max(abs(W)))
        assert symmetric<=max(1e-12,1e-10*scale)
        chol=np.linalg.cholesky(W)
        residual=y-m;ourQ=float(np.dot(chol.T@residual,chol.T@residual))
        Q,prior,sigma,D,shape,AV,peak,priorcenter,searchpeak,mjdoff=record['objective']
        nativeQ=Q-prior-sigma
        den=float(m@W@m);a=float(m@W@y/den);Dopt=D-2.5*np.log10(a) if a>0 else None
        result.append(dict(CID=record['CID'],iteration=record['iteration'],n=record['n'],full_covariance=record['cov'],D=D,fixed_shape=shape,fixed_AV=AV,absolute_peak=peak,prior_center=priorcenter,search_peak=searchpeak,native_total_Q=Q,native_prior_Q=prior,native_sigma_Q=sigma,independent_data_Q=ourQ,native_data_Q=nativeQ,absolute_Q_gap=abs(ourQ-nativeQ),fixed_C_amplitude=a,fixed_C_optimum_D=Dopt,delta_D_optimum=None if Dopt is None else Dopt-D,positive_definite=True,symmetry_gap=symmetric,epochs=[r['epoch'] for r in record['rows']],bands=[r['band'] for r in record['rows']],phase_range=[float(min(arr[:,1])),float(max(arr[:,1]))],rest_wavelength_range=[float(min(arr[:,8])),float(max(arr[:,8]))],inside_rest_wavelength_limits=bool(np.all((arr[:,8]>=arr[:,11])&(arr[:,8]<=arr[:,12]))),weight_matrix=W.tolist()))
    compact=[{k:v for k,v in r.items() if k!='weight_matrix'} for r in result]
    last={}
    for r in result:last.setdefault(r['CID'],{})[r['iteration']]=r
    convergence={}
    for cid,byiteration in last.items():
        ordered=[byiteration[i] for i in sorted(byiteration)];a,b=ordered[-2:]
        convergence[cid]=dict(last_two_iterations=[a['iteration'],b['iteration']],D_change=b['D']-a['D'],same_epoch_mask=a['epochs']==b['epochs'] and a['bands']==b['bands'],last_fixed_C_D_gap=b['delta_D_optimum'],last_Q_gap=b['absolute_Q_gap'])
        if convergence[cid]['same_epoch_mask']:
            A=np.array(a['weight_matrix']);B=np.array(b['weight_matrix']);convergence[cid]['relative_weight_Frobenius_change']=float(np.linalg.norm(B-A)/np.linalg.norm(A))
    answer=dict(status='Independent state check; does not certify unexecuted multiplicativity/real-start/support gates.',log=str(args.log.resolve()),log_sha256=hashlib.sha256(args.log.read_bytes()).hexdigest(),recorded_states=len(records),entry_records=entry,states=compact,last_two_iteration_checks=convergence)
    (args.out/'result.json').write_text(json.dumps(answer,indent=2,allow_nan=False)+'\n')
    (args.out/'executed-source.py').write_bytes(Path(__file__).read_bytes())
    (args.out/'manifest.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in args.out.iterdir() if p.is_file() and p.name!='manifest.json'},indent=2)+'\n')
    print(json.dumps({'states':len(records),'last_two':convergence},indent=2))
if __name__=='__main__':main()
