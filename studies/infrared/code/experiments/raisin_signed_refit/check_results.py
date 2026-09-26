"""Post-fit native acceptance, pair and archived-baseline ledger; no input changes."""
from pathlib import Path
import csv,gzip,json,re,hashlib
import numpy as np
from scipy.optimize import linear_sum_assignment
import paired_refit as p
O=p.OUT

def table(f,key='SN:'):
    if not f.exists() and f.name=='fit.FITRES.TEXT' and (f.parent/'fit').exists():
        # Native chooses literal prefix if the absolute path contains a dot.
        f=f.parent/'fit'
        assert '# TABLE_NAME:    FITRES' in f.read_text()
    s=gzip.open(f,'rt').read() if str(f).endswith('.gz') else f.read_text();rows=[]
    for line in s.splitlines():
        a=line.split()
        if not a:continue
        if a[0]=='VARNAMES:':cols=a[1:]
        elif a[0]==key:rows.append(dict(zip(cols,a[1:])))
    return rows

def gate(run,arm):
    fr=table(run/'fit.FITRES.TEXT');lc=table(run/'fit.LCPLOT.TEXT','OBS:')
    log=(run/'fit.log').read_text();cov={}
    for line in log.splitlines():
        a=line.split()
        if a and a[0]=='PHASE2_HESSIAN:':
            v=np.array(list(map(float,a[2:])));cov[a[1]]=v[-16:].reshape(4,4)
    allrows=[]
    for r in fr:
        cid=r['CID'];h,raw=p.phot(O/'data'/f'RSR_{arm}'/f'{cid}.snana.dat')
        data=[x for x in lc if x['CID']==cid and x['DATAFLAG']=='1']
        assert len(data)==int(float(r['NDOF']))+(3 if '/fixed/' in str(run) else 4),(run,cid,'ndof')
        maxdt=maxdf=maxde=0.;ambiguous=0
        for b in 'griz':
            lr=[x for x in data if x['BAND']==b];rr=[x for x in raw if x['FLT']==b]
            if not lr:continue
            dist=np.abs(np.array([float(x['MJD']) for x in lr])[:,None]-np.array([float(x['MJD']) for x in rr])[None,:])
            i,j=linear_sum_assignment(dist);assert len(i)==len(lr) and np.max(dist[i,j])<.0021,(run,cid,b,'MJD closure')
            ambiguous+=int(np.sum(np.sum(dist<.0021,axis=1)>1))
            for u,v in zip(i,j):
                x,y=lr[u],rr[v];df=abs(float(x['FLUXCAL'])-float(y['FLUXCAL']));de=abs(float(x['FLUXCAL_ERR'])-float(y['FLUXCALERR']))
                # Native LCPLOT uses five significant digits, unlike raw quoted inputs.
                tolflux=max(1e-8,abs(float(y['FLUXCAL']))*5.1e-5);tolerr=max(1e-8,abs(float(y['FLUXCALERR']))*5.1e-5)
                assert df<tolflux and de<tolerr,(run,cid,b,'flux/error printed closure',df,de)
                maxdt=max(maxdt,dist[u,v]);maxdf=max(maxdf,df);maxde=max(maxde,de)
        C=cov[cid];use=[0,1,2] if '/fixed/' in str(run) else [0,1,2,3]
        eig=np.linalg.eigvalsh((C+C.T)[np.ix_(use,use)]/2)
        info={k:float(r[k]) for k in ['PKMJDINI','PKMJD','PKMJDERR','STRETCH','STRETCHERR','AV','AVERR','RV','RVERR','DLMAG','DLMAGERR','NDOF','FITCHI2','FITPROB','ERRFLAG_FIT']}
        info.update(CID=cid,arm=arm,timing='fixed' if '/fixed/' in str(run) else 'free',accepted=len(data),accepted_negative=sum(float(x['FLUXCAL'])<0 for x in data),accepted_ambiguous_date_rows=ambiguous,max_native_date_error_day=float(maxdt),max_printed_flux_error=maxdf,max_printed_quoted_error_difference=maxde,covariance_min_eigenvalue=float(eig.min()),within_stretch_grid=.7<=info['STRETCH']<=1.3,relative_covariance_asymmetry=float(np.max(abs(C-C.T))/max(abs(C).max(),1e-30)))
        allrows.append(info)
    return allrows

def main():
    runs=[]
    for timing in ['free','fixed']:
        for arm in ['R','A','B','H']:runs.extend(gate(O/'fits/full'/timing/arm,arm))
    assert len(runs)==80
    p.csvsave(O/'native-fit-ledger.csv',runs)
    lookup={(r['timing'],r['arm'],r['CID']):r for r in runs};pairs=[]
    for timing in ['free','fixed']:
        for cid in sorted({r['CID'] for r in runs}):
            for before,after,label in [('A','B','author_signed_minus_positive'),('R','H','hybrid_minus_release'),('R','A','authorpositive_minus_release')]:
                a=lookup[timing,before,cid];b=lookup[timing,after,cid]
                pairs.append(dict(CID=cid,timing=timing,comparison=label,**{'delta_'+k:b[k]-a[k] for k in ['PKMJD','STRETCH','AV','DLMAG','FITCHI2','accepted','accepted_negative']},both_inside_stretch_grid=a['within_stretch_grid'] and b['within_stretch_grid']))
    p.csvsave(O/'paired-shifts.csv',pairs)
    archive={r['CID']:r for r in table(O/'author-optical-FITOPT000.FITRES.gz')};baseline=[]
    for cid in sorted({r['CID'] for r in runs}):
        a=lookup['free','R',cid]
        if cid in archive:baseline.append(dict(CID=cid,**{'delta_'+k:a[k]-float(archive[cid][k]) for k in ['PKMJD','STRETCH','AV','DLMAG','NDOF','FITCHI2']}))
    p.csvsave(O/'modern-minus-author-baseline.csv',baseline)
    order=[]
    for arm in ['Bsort','Breverse']:
        d=O/'fits/order/free'/arm
        if not (d/'fit.FITRES.TEXT').exists() and not (d/'fit').exists():continue
        r=gate(d,arm)[0];a=lookup['free','B',r['CID']]
        order.append(dict(arm=arm,CID=r['CID'],**{'delta_'+k:r[k]-a[k] for k in ['PKMJD','STRETCH','AV','DLMAG','FITCHI2','accepted','accepted_negative']}))
    summary={}
    for timing in ['free','fixed']:
        for label in sorted({r['comparison'] for r in pairs}):
            rr=[r for r in pairs if r['timing']==timing and r['comparison']==label]
            summary[timing+'/'+label]={k:{'median':float(np.median([r['delta_'+k] for r in rr])),'minimum':min(r['delta_'+k] for r in rr),'maximum':max(r['delta_'+k] for r in rr)} for k in ['DLMAG','AV','STRETCH','PKMJD']}
    p.save(O/'result.json',dict(status='Native acceptance and quoted-input gates pass; conditional fit sensitivity, not bias-corrected cosmology.',N_fits=len(runs),all_native_ERRFLAG_zero=all(r['ERRFLAG_FIT']==0 for r in runs),covariance_nonpositive=[r for r in runs if r['covariance_min_eigenvalue']<=0],outside_stretch_grid=[r for r in runs if not r['within_stretch_grid']],order_checks=order,summary=summary,all_input_manifest_sha256=p.sha(O/'input-manifest.json'),source_sha256=p.sha(Path(__file__))))
    print(json.dumps({'Nfits':len(runs),'ERRFLAGS':sorted({r['ERRFLAG_FIT'] for r in runs}),'outside_grid':[(r['CID'],r['arm'],r['timing']) for r in runs if not r['within_stretch_grid']],'order':order,'summary':summary},indent=2))

if __name__=='__main__':main()
