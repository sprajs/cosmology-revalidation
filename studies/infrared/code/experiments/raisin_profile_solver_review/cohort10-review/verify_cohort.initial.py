"""Independent cohort saved-output arithmetic. No native process or owner imports."""
from pathlib import Path
from collections import Counter,defaultdict,deque
import csv,hashlib,importlib.util,json,sys,time
import numpy as np
from scipy.linalg import cho_factor,cho_solve

ROOT=Path(__file__).resolve().parents[4];OUT=Path(__file__).resolve().parent
BASE=ROOT/'runs/research_2026_09_26/astra_design/raisin_signed_refit';C=BASE/'fixed-c-profile/cohort10'
LABELS=['Banchor_A','Banchor_B','Aanchor_A','Aanchor_B'];LM={s:i for i,s in enumerate(LABELS)}
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
save=lambda p,x:Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def csvrows(p):
    with p.open() as f:return [{k:(float(v) if k not in ['stage','metric'] and v else v) for k,v in r.items()} for r in csv.DictReader(f)]
def compare(a,b,tol=1e-8):
    if isinstance(a,dict):
        for key,value in a.items():compare(value,b[key],tol)
    elif isinstance(a,(list,tuple)):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y,tol)
    elif isinstance(a,(bool,str)) or a is None:assert a==b,(a,b)
    else:assert abs(a-b)<=tol,(a,b)

def parse_prefix(path,nom,cid):
    flux={};W={};p=None
    with path.open() as f:
        for line in f:
            s=line.split()
            if not s:continue
            if s[0]=='RAISIN_READY:':n=int(s[1]);break
            if s[0].startswith('PHASE2_') and s[0] in ['PHASE2_FLUX:','PHASE2_COVINV:','PHASE2_OBJECTIVE:']:
                assert s[1]==cid
                if s[0]=='PHASE2_FLUX:':flux[int(s[2])]=(s[3],np.array(s[4:],float))
                elif s[0]=='PHASE2_COVINV:':W[int(s[2])]=np.array(s[3:],float)
                else:p=np.array(s[5:9],float)
        else:raise AssertionError('no READY prefix')
    assert n==len(nom['MJD']) and np.array_equal(p,nom['parameters'])
    a=np.array([flux[j][1] for j in range(1,n+1)])
    fields=dict(band=np.array([flux[j][0] for j in range(1,n+1)]),MJD=a[:,0],model_flux=a[:,2],data_flux=a[:,4],data_fluxerr=a[:,5],W=np.array([W[j] for j in range(1,n+1)]))
    for key,value in fields.items():assert np.array_equal(value,nom[key]),key
    return dict(epochs=n,parameters=p.tolist(),exact_fields=list(fields))

def verify(cid,g,cohort):
    start=time.monotonic();r=json.loads((g/'result.json').read_text());p=json.loads((g/'protocol.json').read_text())
    manifest=json.loads((g/'manifest.json').read_text())['files_sha256'];hash_notes=[]
    for f,h in manifest.items():
        if sha(g/f)==h:continue
        assert f=='runner.log',(cid,f)
        raw=(g/f).read_bytes();tail=(json.dumps(r,indent=2)+'\n').encode()
        assert raw.endswith(tail) and hashlib.sha256(raw[:-len(tail)]).hexdigest()==h
        hash_notes.append(dict(file=f,reason='Manifest hashed before final result print; exact earlier prefix verifies.',appended_bytes=len(tail),closed_sha256=sha(g/f),prefix_sha256=h))
    for f,h in p['hashes'].items():assert sha(ROOT/f)==h,(cid,f)
    engines=list(BASE.glob('cohort_profile_engine*.py'));matched=[f for f in engines if sha(f)==p['engine_sha256']];assert len(matched)==1
    assert r['protocol_sha256']==sha(g/'protocol.json') and r['cohort_protocol_sha256']==sha(C/'protocol.json')
    z=np.load(g/'native-profiles.npz');ref=np.load(g/'reference/native.npz');nom=np.load(g/'fixed_reference/native.npz');final=np.load(g/'stream/native.npz');alt=np.load(g/'A_covariance_on_B/native.npz');accepted=np.load(g/'A_acceptance/native.npz')
    for key in ['band','MJD','data_flux','data_fluxerr']:
        for other in [nom,final,alt,z]:assert np.array_equal(ref[key],other[key]),(cid,key)
    for key in ['parameters','W','C','model_flux']:assert np.array_equal(nom[key],final[key]),(cid,key)
    assert np.array_equal(nom['parameters'],z['reference_parameters']) and np.array_equal(alt['parameters'],p['alternate_anchor_parameters'])
    prefix=parse_prefix(g/'stream/fit.log',nom,cid)
    keys=lambda a:list(zip(a['band'],a['MJD'],a['data_flux'],a['data_fluxerr']))
    bk=keys(ref);ak=keys(accepted);queues=defaultdict(deque);groups=defaultdict(list)
    for j,key in enumerate(bk):queues[key].append(j);groups[key].append(j)
    ai=np.array([queues[key].popleft() for key in ak]);assert np.array_equal(ai,z['A_indices'])
    assert len(ai)==cohort['expected_masks'][cid]['A'] and len(bk)==cohort['expected_masks'][cid]['B']
    assert Counter(ak)==Counter(key for key in bk if key[2]>0)
    duperrors=[]
    for indices in groups.values():
        for index in indices[1:]:
            perm=np.arange(len(bk));perm[indices[0]],perm[index]=perm[index],perm[indices[0]]
            for full in [ref['C'],alt['C']]:
                error=float(np.max(abs(full[np.ix_(perm,perm)]-full))/np.max(abs(full)));assert error<1e-12;duperrors.append(error)
    spec=importlib.util.spec_from_file_location('ownreader',OUT.parent/'native-proof-review/verify_native_proof.py');reader=importlib.util.module_from_spec(spec);spec.loader.exec_module(reader);reader.CID=cid
    covariance_errors=[]
    for folder,export in [('reference',ref),('fixed_reference',nom),('A_covariance_on_B',alt)]:
        raw=reader.parse_log(g/folder/'fit.log');assert np.array_equal(raw['W'],export['W'])
        inv=cho_solve(cho_factor(raw['W'],lower=True),np.eye(raw['n']));error=float(np.max(abs(inv-export['C'])));assert error<1e-8;covariance_errors.append(error)
        residual=export['data_flux']-export['model_flux'];objective=float(residual@raw['W']@residual+export['prior_chi2']);assert abs(objective-export['total_chi2'])<1e-7
    assert np.array_equal(z['covariance_B'],ref['C']) and np.array_equal(z['covariance_A_anchor'],alt['C'])
    xy=z['coordinates'];H=z['model_means'];saved=z['metric_results'];y=z['data_flux'];Dref=z['reference_parameters'][0];n=len(xy)
    lookup={tuple(v):i for i,v in enumerate(xy)};assert len(lookup)==n
    assert np.isfinite(H).all();re=np.empty_like(saved);information=np.empty((n,4));metricerrors={};support=[]
    for k,label in enumerate(LABELS):
        ix=ai if k%2==0 else np.arange(len(y));full=ref['C'] if k<2 else alt['C'];sub=full[np.ix_(ix,ix)]
        assert np.max(abs(sub-sub.T))/np.max(abs(sub))<1e-12
        W=cho_solve(cho_factor(sub,lower=True),np.eye(len(ix)));normal=0.
        for st in range(0,n,4096):
            en=min(n,st+4096);h=H[st:en][:,ix];wh=h@W;q=np.einsum('ij,ij->i',wh,h);a=(wh@y[ix])/q
            assert np.all(q>0) and np.all(a>0)
            delta=y[ix]-a[:,None]*h;Q=np.einsum('ij,ij->i',delta@W,delta);D=Dref-2.5*np.log10(a)
            re[st:en,k]=np.column_stack([Q,D,a]);information[st:en,k]=q
            normal=max(normal,float(np.max(abs(np.einsum('ij,ij->i',wh,delta))/np.maximum(1,q))))
        error=np.max(abs(re[:,k]-saved[:,k]),axis=0);assert np.max(error)<1e-8,(cid,label,error)
        scale=re[:,k,2,None]*H;valid=(H>1e-5).all() and (H<1e9).all() and (scale>1e-5).all() and (scale<1e9).all() and (re[:,k,1]>10).all() and (re[:,k,1]<60).all()
        support.append(bool(valid));metricerrors[label]=dict(Q_D_amplitude_max_error=error.tolist(),normal_equation_scaled_max=normal,C_min_eigenvalue=float(np.linalg.eigvalsh(sub).min()))
        assert normal<1e-10
    profiles=csvrows(g/'shape-profiles.csv');branches=csvrows(g/'all-AV-branches.csv');boundaries=csvrows(g/'algorithm-boundary-probes.csv');ampchecks=csvrows(g/'competitive-amplitude-checks.csv')
    for rows in [profiles,branches,boundaries]:
        for row in rows:
            value=re[lookup[(row['shape'],row['AV'])],LM[row['metric']]]
            assert np.max(abs(value-np.array([row['Q'],row['DLMAG'],row['amplitude']])))<1e-8
    knots=np.array(p['domain']['shape_knots']);assert np.array_equal(knots,cohort['same_method']['domain']['shape_knots'])
    for stage,nsub,nav in [('coarse',2,33),('fine',4,65)]:
        shapes=np.unique(np.concatenate([np.linspace(a,b,nsub+1) for a,b in zip(knots[:-1],knots[1:])]))
        avs=np.linspace(-1,2,nav)
        for s in shapes:assert all((s,a) in lookup for a in avs)
        for label in LABELS:
            rows=[x for x in profiles if x['stage']==stage and x['metric']==label];assert np.array_equal([x['shape'] for x in rows],shapes)
            for x in rows:
                Q=np.array([re[lookup[(x['shape'],a)],LM[label],0] for a in avs]);assert x['Q']<=Q.min()+1e-8
                assert abs(x['AV_low_Q']-Q[0])<1e-8 and abs(x['AV_high_Q']-Q[-1])<1e-8
    step=float(np.float32(np.float32(knots[-1]-knots[0])/np.float32(len(knots)-1)))
    expected=np.array([knots[0]+i*step+off for i in range(1,len(knots)-1) for off in [-1e-8,0.,1e-8]])
    for label in LABELS:assert np.array_equal([x['shape'] for x in boundaries if x['metric']==label],expected)
    summaries={};geometry={}
    for k,label in enumerate(LABELS):
        rows=[x for x in profiles if x['metric']==label];fine=[x for x in rows if x['stage']=='fine'];coarse=[x for x in rows if x['stage']=='coarse']
        best=min(rows,key=lambda x:x['Q']);co=min(coarse,key=lambda x:x['Q']);compare(r['metrics'][label]['finest'],best)
        assert best['Q']-re[:,k,0].min()<1e-3
        avedge=min(min(x['AV_low_Q'],x['AV_high_Q']) for x in fine)-best['Q'];sedge=min(fine[0]['Q'],fine[-1]['Q'])-best['Q'];dq=co['Q']-best['Q'];dd=co['DLMAG']-best['DLMAG']
        summaries[label]=dict(coarse_minus_finest_Q=dq,coarse_minus_finest_DLMAG=dd,coarse_fine_tolerance_pass=abs(dq)<1e-3 and abs(dd)<1e-3,AV_edge_within9=avedge<=9,shape_edge_within1=sedge<=1)
        compare(summaries[label],r['metrics'][label]);summaries[label].update(AV_edge_delta_Q=avedge,shape_edge_delta_Q=sedge,minimum=best)
        modes=[]
        for j,x in enumerate(fine):
            if (j==0 or x['Q']<=fine[j-1]['Q']) and (j==len(fine)-1 or x['Q']<=fine[j+1]['Q']):
                modes.append(dict(x,delta_Q=x['Q']-best['Q']))
                if 0<j<len(fine)-1:assert any(fine[j-1]['shape']<=z1['shape']<=fine[j+1]['shape'] for z1 in rows if z1['stage']=='shape_polish_final')
        levels={}
        for threshold in [1.,4.,9.]:
            selected=re[:,k,0]<=best['Q']+threshold;a=re[selected,k,2];half=np.sqrt(np.maximum(0,best['Q']+threshold-re[selected,k,0])/information[selected,k]);low=a-half;high=a+half
            minD=10**(-.4*(60-Dref));maxD=10**(-.4*(10-Dref));nativeMin=np.max(1e-5/H[selected],axis=1);nativeMax=np.min(1e9/H[selected],axis=1)
            edgeD=bool(np.any(low<=minD)|np.any(high>=maxD));edgeN=bool(np.any(low<=nativeMin)|np.any(high>=nativeMax))
            low=np.maximum(np.maximum(low,minD),nativeMin);high=np.minimum(np.minimum(high,maxD),nativeMax);assert np.all(low<=high)
            dlo=Dref-2.5*np.log10(high);dhi=Dref-2.5*np.log10(low);sel=xy[selected];ridge=[x['DLMAG'] for x in fine if x['Q']<=best['Q']+threshold]
            levels[str(threshold)]=dict(amplitude_inclusive_distance_range=[float(dlo.min()),float(dhi.max())],shape_profile_ridge_distance_range=[min(ridge),max(ridge)],finite_native_grid_coordinates=int(selected.sum()),shape_box_touched=bool(np.any(abs(sel[:,0]-knots[0])<5e-7)|np.any(abs(sel[:,0]-knots[-1])<5e-7)),AV_box_touched=bool(np.any(abs(sel[:,1]+1)<1e-9)|np.any(abs(sel[:,1]-2)<1e-9)),distance_box_touched=edgeD,native_magnitude_box_touched=edgeN)
        geometry[label]=dict(minimum_Q=best['Q'],cached_minimum_Q=float(re[:,k,0].min()),local_modes=modes,level_sets=levels)
    paired={}
    for anchor in ['B','A']:
        a=anchor+'anchor_A';b=anchor+'anchor_B';amin=summaries[a]['minimum'];bmin=summaries[b]['minimum'];match=min(geometry[b]['local_modes'],key=lambda x:abs(x['shape']-amin['shape']))
        paired[anchor]=dict(global_minimum_delta_DLMAG=bmin['DLMAG']-amin['DLMAG'],matched_B_mode_nearest_A_shape=match,matched_mode_delta_DLMAG=match['DLMAG']-amin['DLMAG'])
        assert abs(paired[anchor]['global_minimum_delta_DLMAG']-r['delta_DLMAG_B_minus_A'][anchor])<1e-10
    geom=json.loads((g/'branch-level-geometry.json').read_text())
    for f,h in geom['input_hashes'].items():assert sha(g/f)==h
    compare(geometry,geom['metrics']);compare(paired,geom['paired'])
    assert geom['source_sha256']==sha(BASE/'cohort_profile_geometry_v2.py')
    numerical=bool(r['computational_completion'] and all(support) and all(x['coarse_fine_tolerance_pass'] and not x['AV_edge_within9'] and not x['shape_edge_within1'] for x in summaries.values()))
    assert numerical==r['numerical_gate_pass'] and len(ampchecks)==r['competitive_amplitude_checks']
    assert max(x['relative_error'] for x in ampchecks)<2e-8 and r['total_oracle_calls']<=100000
    zero=None
    if len(ai)==len(y):
        zero={anchor:float(np.max(abs(re[:,2*k]-re[:,2*k+1]))) for k,anchor in enumerate(['B','A'])};assert max(zero.values())<1e-9
        assert max(abs(x['global_minimum_delta_DLMAG']) for x in paired.values())<1e-9
    output=dict(CID=cid,status='PASS: independent arithmetic and declared gates verified',numerical_gate_pass=numerical,active_path=str(g.relative_to(C)),manifest_entries=len(manifest),manifest_log_append_ledger=hash_notes,engine=matched[0].name,epochs=dict(A=len(ai),B=len(y)),negative_rows=int(np.sum(y<0)),reference_nominal_final_identity=True,prefix=prefix,duplicate_covariance_swap_checks=len(duperrors),max_duplicate_swap_error=max(duperrors,default=0.),raw_C_inverse_max_error=max(covariance_errors),vectors=n,metric_errors=metricerrors,metrics=summaries,geometry=geometry,paired=paired,zero_change_control=zero,seconds=time.monotonic()-start,input_hashes={f:sha(g/f) for f in ['manifest.json','protocol.json','result.json','native-profiles.npz','shape-profiles.csv','branch-level-geometry.json','stream/fit.log']})
    save(OUT/f'{cid}.json',output);return output

def main():
    cohort=json.loads((C/'protocol.json').read_text());index=json.loads((C/'active-attempts.json').read_text())['attempts'];out={};pending=[]
    for f,h in cohort['hashes'].items():assert sha(ROOT/f)==h,f
    requested=sys.argv[1:] or cohort['cohort']
    for cid in requested:
        g=C/index[cid]
        if not all((g/f).exists() for f in ['manifest.json','result.json','branch-level-geometry.json','native-profiles.npz']):pending.append(cid);continue
        prior=OUT/f'{cid}.json'
        if prior.exists():
            cached=json.loads(prior.read_text())
            if cached.get('checker_sha256')==sha(Path(__file__)) and all(sha(g/f)==h for f,h in cached['input_hashes'].items()):out[cid]=cached;continue
        r=verify(cid,g,cohort);r['checker_sha256']=sha(Path(__file__));save(prior,r);out[cid]=r
        print(json.dumps(dict(CID=cid,numerical_gate_pass=r['numerical_gate_pass'],paired={k:v['global_minimum_delta_DLMAG'] for k,v in r['paired'].items()},seconds=r['seconds'])),flush=True)
    complete=set(out)==set(cohort['cohort']) and not pending
    overall=dict(status='Full cohort independently checked' if complete else 'Partial review; no whole-cohort certification',completed=list(out),pending=pending,all_ten_complete=complete,numerical_gate_pass_all=complete and all(r['numerical_gate_pass'] for r in out.values()),zero_controls={cid:r['zero_change_control'] for cid,r in out.items() if r['zero_change_control'] is not None},pairs={cid:r['paired'] for cid,r in out.items()},cohort_protocol_sha256=sha(C/'protocol.json'),active_index_sha256=sha(C/'active-attempts.json'),checker_sha256=sha(Path(__file__)))
    save(OUT/'result.json',overall)
if __name__=='__main__':main()
