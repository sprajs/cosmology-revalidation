"""Independent saved-array/raw-log review. No native execution or fit imports."""
from pathlib import Path
import csv, hashlib, importlib.util, json, time
import numpy as np
from scipy.linalg import cho_factor, cho_solve
from astropy.io import fits

ROOT=Path(__file__).resolve().parents[4]
OUT=Path(__file__).resolve().parent
BASE=ROOT/'runs/research_2026_09_26/astra_design/raisin_signed_refit'
G=BASE/'fixed-c-profile/global-profile'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(name,value): (OUT/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
def rows(name):
    with (G/name).open() as f:
        return [{k:(float(v) if k not in ('stage','metric') and v else v) for k,v in r.items()} for r in csv.DictReader(f)]
def main():
    start=time.monotonic()
    manifest=json.loads((G/'manifest.json').read_text())['files_sha256']
    assert all(sha(G/p)==h for p,h in manifest.items())
    protocol=json.loads((G/'protocol.json').read_text()); amendment=json.loads((G/'amendment-mean-domain.json').read_text())
    assert sha(BASE/'global_fixed_profile.py')==amendment['runner_sha256']
    native=np.load(G/'native-profiles.npz'); reported=json.loads((G/'result.json').read_text())
    coords=native['coordinates']; H=native['model_means']; saved=native['metric_results']; ai=native['A_indices']
    y=native['data_flux']; e=native['data_fluxerr']; par=native['reference_parameters']; n=len(coords)
    assert H.shape==(n,74) and saved.shape==(n,4,3)
    assert np.isfinite(H).all() and np.min(H)>1e-5 and np.max(H)<1e9
    lookup={tuple(c):i for i,c in enumerate(coords)}; assert len(lookup)==n
    labels=['Banchor_A','Banchor_B','Aanchor_A','Aanchor_B']; lm={s:i for i,s in enumerate(labels)}
    spec=importlib.util.spec_from_file_location('own_raw_reader',OUT.parent/'native-proof-review/verify_native_proof.py')
    reader=importlib.util.module_from_spec(spec);spec.loader.exec_module(reader)
    reference=np.load(BASE/'native-profile-gate-v2/reference/native.npz')
    nominal=np.load(BASE/'native-profile-gate-v2/fixed/native.npz')
    alternate=np.load(G/'A_covariance_on_B/native.npz')
    accepted=np.load(G/'A_acceptance/native.npz')
    assert np.array_equal(par,nominal['parameters'])
    assert np.array_equal(alternate['parameters'],protocol['alternate_anchor_parameters'])
    for key in ['band','MJD','data_flux','data_fluxerr']:
        assert np.array_equal(native[key],reference[key]) and np.array_equal(native[key],alternate[key])
        assert np.array_equal(native[key][ai],accepted[key])
    assert np.array_equal(np.sort(ai),np.flatnonzero(y>0)) and len(ai)==65
    assert np.array_equal(native['covariance_B'],reference['C'])
    assert np.array_equal(native['covariance_A_anchor'],alternate['C'])
    raw_cov=[]
    for directory,z in [(BASE/'native-profile-gate-v2/reference',reference),(G/'A_covariance_on_B',alternate),(G/'A_acceptance',accepted)]:
        raw=reader.parse_log(directory/'fit.log')
        assert np.array_equal(raw['W'],z['W'])
        C=cho_solve(cho_factor(raw['W'],lower=True),np.eye(raw['n']))
        err=float(np.max(abs(C-z['C'])));assert err<1e-9
        raw_cov.append(dict(directory=str(directory.relative_to(ROOT)),inverse_raw_precision_max_error=err))
    recomputed=np.empty_like(saved); infos=np.empty((n,4)); matrices={}; covchecks={}
    for k,label in enumerate(labels):
        full=native['covariance_B'] if k<2 else native['covariance_A_anchor']
        ix=ai if k%2==0 else np.arange(74)
        C=full[np.ix_(ix,ix)]
        W=cho_solve(cho_factor(C,lower=True),np.eye(len(ix)))
        matrices[label]=(ix,W)
        errors=[]; normals=[]
        for st in range(0,n,4096):
            en=min(n,st+4096); h=H[st:en][:,ix]; v=y[ix]
            wh=h@W; q=np.einsum('ij,ij->i',wh,h); b=wh@v
            a=b/q; D=par[0]-2.5*np.log10(a);r=v-a[:,None]*h
            Q=np.einsum('ij,ij->i',r@W,r)
            recomputed[st:en,k]=np.column_stack([Q,D,a]);infos[st:en,k]=q
            normals.extend(abs(np.einsum('ij,ij->i',wh,r))/np.maximum(1,q))
        assert np.all(infos[:,k]>0)
        assert np.all(recomputed[:,k,2]>0) and np.all((recomputed[:,k,1]>10)&(recomputed[:,k,1]<60))
        scaled=H*recomputed[:,k,2,None];assert scaled.min()>1e-5 and scaled.max()<1e9
        covchecks[label]=dict(min_eigenvalue=float(np.linalg.eigvalsh(C).min()),asymmetry=float(np.max(abs(C-C.T))/np.max(abs(C))),normal_equation_max_scaled=float(max(normals)),
            max_Q_D_amplitude_error=np.max(abs(recomputed[:,k]-saved[:,k]),axis=0).tolist(),scaled_mean_range=[float(scaled.min()),float(scaled.max())])
    assert np.max(abs(recomputed-saved))<1e-9
    profiles=rows('shape-profiles.csv'); branches=rows('all-AV-branches.csv'); bounds=rows('algorithm-boundary-probes.csv'); amps=rows('competitive-amplitude-checks.csv')
    csv_errors={}
    for filename,rr in [('shape-profiles',profiles),('all-AV-branches',branches),('algorithm-boundary-probes',bounds)]:
        dif=[]
        for r in rr:
            j=lookup[(r['shape'],r['AV'])];k=lm[r['metric']]
            dif.append(abs(np.array([r['Q'],r['DLMAG'],r['amplitude']])-recomputed[j,k]))
        csv_errors[filename]=dict(rows=len(rr),max_Q_D_amplitude_error=np.max(dif,axis=0).tolist())
        assert np.max(dif)<1e-9
    # Check raw native stream, independent of NPZ and CSV materialization.
    seq=0;baseline=None; pending=None; raw_amp_errors=[];point_errors=[]
    with (G/'stream/fit.log').open() as f:
        for line in f:
            if not line.startswith('RAISIN_MEAN:'):continue
            q=line.split();seq+=1;assert int(q[1])==seq and int(q[2])==74
            v=np.array(list(map(float,q[3:])));p,h=v[:4],v[4:]
            if seq==1:baseline=h.copy();assert np.array_equal(p,par)
            elif seq<=n+1:
                j=seq-2;expected=par.copy();expected[1:3]=coords[j]
                assert np.array_equal(p,expected) and np.array_equal(h,H[j])
            elif seq<=n+1+2*len(amps):
                t=seq-(n+2);r=amps[t//2];expected=par.copy();expected[1:3]=r['shape'],r['AV']
                if t%2==0:
                    assert np.array_equal(p,expected);pending=h.copy()
                else:
                    expected[0]=r['DLMAG'];assert np.array_equal(p,expected)
                    a=10**(-.4*(r['DLMAG']-par[0]));pred=a*pending
                    err=float(np.max(abs(h-pred))/np.max(abs(pred)))
                    # Recompute with the amplitude saved for this exact coordinate.
                    a_saved=saved[lookup[(r['shape'],r['AV'])],lm[r['metric']],2]
                    err_saved=float(np.max(abs(h-a_saved*pending))/np.max(abs(a_saved*pending)))
                    assert abs(err_saved-r['relative_error'])<1e-14
                    assert err<2e-8
                    raw_amp_errors.append(err);point_errors.append(float(np.max(abs((h-pred)/pred))))
            else:
                assert seq==n+2+2*len(amps) and np.array_equal(p,par) and np.array_equal(h,baseline)
    assert seq==reported['total_oracle_calls']==n+2+2*len(amps)
    final=np.load(G/'stream/native.npz')
    for key in ['band','MJD','data_flux','data_fluxerr','C','W']:assert np.array_equal(final[key],nominal[key])
    # Native knots, full two/four subdivision grids and AV scans.
    model=ROOT/'sources/repos/djones1040__RAISIN_DataRelease@a383c4b/model/snoopy.B18/SNooPy_B18.fits'
    knots=np.asarray(fits.getdata(model,'LUMI-GRID').field(0),float)
    assert np.array_equal(knots,protocol['domain']['shape_knots'])
    gridchecks={}; modes={}; minima={}; intervals={}; allmodes=[]
    for stage,subdiv,nav in [('coarse',2,33),('fine',4,65)]:
        sg=np.unique(np.concatenate([np.linspace(a,b,subdiv+1) for a,b in zip(knots[:-1],knots[1:])]))
        ag=np.linspace(-1,2,nav)
        for s in sg:
            assert all((s,a) in lookup for a in ag)
        for label in labels:
            rr=[r for r in profiles if r['stage']==stage and r['metric']==label]
            assert np.array_equal([r['shape'] for r in rr],sg)
            for r in rr:
                values=np.array([recomputed[lookup[(r['shape'],a)],lm[label],0] for a in ag])
                assert r['Q']<=values.min()+1e-9
                assert abs(r['AV_low_Q']-values[0])<1e-9 and abs(r['AV_high_Q']-values[-1])<1e-9
        gridchecks[stage]=dict(shapes=len(sg),AV_points=nav,all_grid_coordinates_saved=True)
    step=float(np.float32(np.float32(knots[-1]-knots[0])/np.float32(len(knots)-1)))
    expected=np.array([knots[0]+i*step+o for i in range(1,len(knots)-1) for o in [-1e-8,0.,1e-8]])
    for label in labels:
        rr=[r for r in bounds if r['metric']==label]
        assert np.array_equal([r['shape'] for r in rr],expected)
        fine=[r for r in profiles if r['metric']==label and r['stage']=='fine']
        for r in rr:
            near=min(fine,key=lambda f:abs(f['shape']-r['shape']))
            assert near['AV']==r['AV']
    # Retain every sampled shape-profile local mode and descriptive level-set envelopes.
    for k,label in enumerate(labels):
        rr=[r for r in profiles if r['metric']==label]; fi=[r for r in rr if r['stage']=='fine'];co=[r for r in rr if r['stage']=='coarse']
        best=min(rr,key=lambda r:r['Q']);actual=reported['metrics'][label]['finest'];assert all(best[key]==value for key,value in actual.items())
        j=int(np.argmin(recomputed[:,k,0]));allbest=dict(shape=float(coords[j,0]),AV=float(coords[j,1]),Q=float(recomputed[j,k,0]),DLMAG=float(recomputed[j,k,1]))
        assert best['Q']-allbest['Q']<1e-3
        coarse=min(co,key=lambda r:r['Q']);assert abs(coarse['Q']-best['Q'])<1e-3 and abs(coarse['DLMAG']-best['DLMAG'])<1e-3
        avedge=min(min(r['AV_low_Q'],r['AV_high_Q']) for r in fi)-best['Q'];sedge=min(fi[0]['Q'],fi[-1]['Q'])-best['Q']
        assert avedge>9 and sedge>1
        local=[]
        for t in range(1,len(fi)-1):
            if fi[t]['Q']<=fi[t-1]['Q'] and fi[t]['Q']<=fi[t+1]['Q']:
                lo,hi=fi[t-1]['shape'],fi[t+1]['shape'];select=(coords[:,0]>=lo)&(coords[:,0]<=hi)
                jj=np.flatnonzero(select)[np.argmin(recomputed[select,k,0])]
                row=dict(metric=label,fine_shape=fi[t]['shape'],shape=float(coords[jj,0]),AV=float(coords[jj,1]),Q=float(recomputed[jj,k,0]),delta_Q=float(recomputed[jj,k,0]-best['Q']),DLMAG=float(recomputed[jj,k,1]))
                assert any(lo<=r['shape']<=hi for r in rr if r['stage']=='shape_polish_final')
                local.append(row);allmodes.append(row)
        local.sort(key=lambda r:r['Q']);modes[label]=local
        minima[label]=dict(reported_minimum=best,minimum_of_all_saved_probes=allbest,reported_minus_all_saved_Q=best['Q']-allbest['Q'],coarse_delta_Q=coarse['Q']-best['Q'],coarse_delta_D=coarse['DLMAG']-best['DLMAG'],AV_edge_delta_Q=avedge,shape_edge_delta_Q=sedge,number_sampled_shape_local_minima=len(local))
        levels={}
        for delta in [1.,4.,9.]:
            select=recomputed[:,k,0]<=best['Q']+delta
            a=recomputed[select,k,2];dq=np.maximum(0,best['Q']+delta-recomputed[select,k,0]);half=np.sqrt(dq/infos[select,k])
            assert np.all(a>half)
            dlo=par[0]-2.5*np.log10(a+half);dhi=par[0]-2.5*np.log10(a-half)
            assert np.all((dlo>10)&(dhi<60))
            assert np.min((a-half)[:,None]*H[select])>1e-5 and np.max((a+half)[:,None]*H[select])<1e9
            levels[str(delta)]=dict(saved_coordinates=int(select.sum()),profiled_amplitude_D_range=[float(recomputed[select,k,1].min()),float(recomputed[select,k,1].max())],including_analytic_amplitude_level_set_D_range=[float(dlo.min()),float(dhi.max())],shape_domain_edge_within_level=sedge<=delta,AV_domain_edge_within_level=avedge<=delta,interpretation='Descriptive fixed-quadratic level set over saved shape/AV coordinates; not a calibrated confidence interval or posterior. Envelopes may be incomplete between saved coordinates or where a domain boundary is reached.')
        intervals[label]=levels
    # Marginal/conditional block decomposition at both best B shape branches.
    ni=np.setdiff1d(np.arange(74),ai);C=native['covariance_B'];CA=C[np.ix_(ai,ai)];WA=cho_solve(cho_factor(CA,lower=True),np.eye(65))
    R=C[np.ix_(ni,ai)]@WA;V=C[np.ix_(ni,ni)]-R@C[np.ix_(ai,ni)];WV=cho_solve(cho_factor(V,lower=True),np.eye(9));decomposition=[]
    for mode in modes['Banchor_B']:
        j=lookup[(mode['shape'],mode['AV'])];r=y-recomputed[j,1,2]*H[j];cn=r[ni]-R@r[ai]
        qp=float(r[ai]@WA@r[ai]);qn=float(cn@WV@cn);qt=float(recomputed[j,1,0]);assert abs(qp+qn-qt)<1e-9
        decomposition.append(dict(**mode,positive_marginal_Q=qp,negative_conditional_Q=qn,sum_minus_full_Q=qp+qn-qt))
    # Independently reproduce the parent's four-point diagnostic, without its code.
    parentpath=ROOT/'runs/research_2026_09_26/raisin_profile_decomposition/result.json'
    parent=json.loads(parentpath.read_text());parent_diffs={}
    for name,r0 in parent['points'].items():
        j=lookup[(r0['shape'],r0['AV'])];k=0 if name.endswith('A_amplitude') else 1
        a=recomputed[j,k,2];r=y-a*H[j];cn=r[ni]-R@r[ai]
        own=dict(amplitude=a,DLMAG=recomputed[j,k,1],Q_positive_marginal=float(r[ai]@WA@r[ai]),Q_negative_conditional=float(cn@WV@cn),Q_full=float(r@matrices['Banchor_B'][1]@r))
        parent_diffs[name]={key:float(own[key]-r0[key]) for key in own}
    assert max(abs(v) for r in parent_diffs.values() for v in r.values())<1e-9
    save('parent-decomposition-check.json',dict(source_sha256=sha(parentpath),differences=parent_diffs,max_absolute_error=max(abs(v) for r in parent_diffs.values() for v in r.values())))
    pairs={anchor: minima[anchor+'anchor_B']['reported_minimum']['DLMAG']-minima[anchor+'anchor_A']['reported_minimum']['DLMAG'] for anchor in ['A','B']}
    for anchor,value in pairs.items():assert abs(value-reported['delta_DLMAG_B_minus_A'][anchor])<1e-12
    save('all-local-modes.json',modes);save('level-sets.json',intervals);save('block-decomposition.json',decomposition)
    result=dict(status='Independent saved global profile checks PASS; no native fits run',hashes_checked=len(manifest),native_vectors=n,raw_native_calls_verified=seq,accepted_rows=dict(A=65,B=74,negative_B=9),all_C_data_masks_match_raw_exports=True,raw_covariance_checks=raw_cov,metric_reconstruction=covchecks,CSV_checks=csv_errors,
      raw_competitive_amplitude_checks=len(amps),max_competitive_amplitude_error=float(max(raw_amp_errors)),max_competitive_pointwise_amplitude_error=float(max(point_errors)),native_replay_exact=True,
      grid_coverage=gridchecks,native_float32_shape_step=step,algorithm_boundary_probes_per_metric=len(expected),all_boundary_coordinates_exact=True,minima=minima,delta_DLMAG_B_minus_A=pairs,
      covariance_anchor_pair_difference=pairs['A']-pairs['B'],level_sets=intervals,local_modes=modes,
      native_mean_range=[float(H.min()),float(H.max())],numerical_gates_recomputed=True,
      limitation='Finite restricted-domain coverage and explicit boundary probes, not a mathematical continuum-global proof. Fitted minima are not uniquely distance-identifying. Different A/B row objectives are not a likelihood ratio; C is state frozen, negative AV is empirical, no selected likelihood or physical correction.',seconds=time.monotonic()-start,
      input_sha256={str(p.relative_to(ROOT)):sha(p) for p in [G/'manifest.json',G/'protocol.json',G/'amendment-mean-domain.json',G/'native-profiles.npz',G/'stream/fit.log',model,Path(__file__)]})
    save('result.json',result);print(json.dumps(dict(status=result['status'],seconds=result['seconds'],pairs=pairs,metric_errors={k:v['max_Q_D_amplitude_error'] for k,v in covchecks.items()},modes=modes),indent=2))
if __name__=='__main__':main()
