"""Paired generated-attempt selection response with explicit failed attempts."""
from pathlib import Path
import argparse,hashlib,json
import numpy as np,pandas as pd

ROOT=Path(__file__).resolve().parents[3];OUT=ROOT/'phase2/forward_discrimination'
MODELS=['P21','BS21','G10','P21_dmplus010','P21_dmminus010','P21_dmz020','P21_rho000','P21_rho090','P21_noisetrue120']
EDGES=np.array([0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1,1.2,2])

def contrast(a,b,cluster):
    a=np.asarray(a,dtype=float);b=np.asarray(b,dtype=float);d=a-b;n=len(d);delta=d.mean();se=d.std(ddof=1)/np.sqrt(n)
    ag=pd.DataFrame({'d':d-delta,'cluster':cluster}).groupby('cluster').d.sum();g=len(ag)
    sec=np.sqrt(g/(g-1)*np.sum(ag.to_numpy()**2))/n if g>1 else np.nan
    p=a.mean();z=1.959963984540054;den=1+z*z/n;ctr=(p+z*z/(2*n))/den;rad=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return {'n_generated':n,'n_selected':int(a.sum()),'fraction':float(p),'wilson95':[float(ctr-rad),float(ctr+rad)],'n_selected_P21':int(b.sum()),'delta_vs_P21':float(delta),'paired_standard_error':float(se),'paired_ci95':[float(delta-z*se),float(delta+z*se)],'libid_cluster_standard_error':float(sec),'libid_cluster_ci95':[float(delta-z*sec),float(delta+z*sec)],'n_libid_clusters':g,'newly_selected_vs_P21':int(((a==1)&(b==0)).sum()),'lost_vs_P21':int(((a==0)&(b==1)).sum())}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--require-fits',action='store_true');args=ap.parse_args()
    root=ROOT/'phase2/hierarchy/forward';summ=json.loads((root/'summary.json').read_text());ready={m['model']:m['fit_complete'] for m in summ['models']}
    if args.require_fits:assert all(ready.values())
    frames={m:pd.read_csv(root/(m+'-attempts.csv.gz')) for m in MODELS};base=frames['P21'];out={};inputs=[Path(__file__),root/'summary.json']
    for m,f in frames.items():
        inputs.append(root/(m+'-attempts.csv.gz'))
        for c in ['generated_attempt_index','GENZ','LIBID','GALID','PEAKMJD']:np.testing.assert_array_equal(base[c],f[c])
        assert len(f)==26518 and f.generated_attempt_index.is_unique
        stages=['written']+(['fit_pass','basic_quality_pass'] if ready[m] and ready['P21'] else [])
        out[m]={}
        for stage in stages:
            byz=[]
            for lo,hi in zip(EDGES[:-1],EDGES[1:]):
                pick=(f.GENZ>=lo)&(f.GENZ<hi)
                if pick.sum()>1:byz.append({'z_min':lo,'z_max':hi,**contrast(f.loc[pick,stage],base.loc[pick,stage],f.loc[pick,'LIBID'])})
            out[m][stage]={'overall':contrast(f[stage],base[stage],f.LIBID),'by_redshift':byz}
    (OUT/'selection-response.json').write_text(json.dumps({'models':out,'fit_complete':ready,'scope':'Only fixed forward selection experiments. Common design stream paired by generated attempt, not reused CID. No normalization for historical BBC/classifier/systematic-intersection cohort. One seed; intervals quantify conditional Monte Carlo sampling, not physical calibration uncertainty.'},indent=2)+'\n')
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();(OUT/'selection-manifest.json').write_text(json.dumps({'inputs_sha256':{str(p.relative_to(ROOT)):sha(p) for p in inputs},'output_sha256':sha(OUT/'selection-response.json')},indent=2)+'\n')
    for m,x in out.items():print(m,{s:(r['overall']['n_selected'],r['overall']['delta_vs_P21'],r['overall']['paired_ci95']) for s,r in x.items()})

if __name__=='__main__':main()
