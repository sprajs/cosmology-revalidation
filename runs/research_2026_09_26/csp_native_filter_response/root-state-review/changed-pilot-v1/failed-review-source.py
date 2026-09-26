"""Root independent arithmetic/scaling/start checks of all native pilot probes."""
from pathlib import Path
import argparse,json,hashlib,csv
from collections import Counter
import numpy as np
from review_csp_native_states import parse
ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'runs/research_2026_09_26/csp_native_filter_response'
I=ROOT/'runs/research_2026_09_26/raisin_profile_solver_review/csp-filter-interpretation/native-design/instrumentation'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--probe',type=Path,default=I/'probes-v3');ap.add_argument('--name',default='nominal-probes');ap.add_argument('--changed',action='store_true');a=ap.parse_args()
    out=BASE/'root-state-review'/a.name;out.mkdir(exist_ok=False)
    fits=a.probe/'fits';states={};entries={};qgaps=[];hashes={}
    for d in sorted(fits.iterdir()):
        if not d.is_dir():continue
        path=d/'fit.log';records,ent=parse(path);states[d.name]=records;entries[d.name]=ent;hashes[str(path.relative_to(ROOT))]=sha(path)
        for r in records:
            arr=np.array([row['values'] for row in r['rows']]);W=np.array(r['weights']);W=W if r['cov'] else np.diag(W[:,0]);L=np.linalg.cholesky(W)
            residual=arr[:,4]-arr[:,2];Q=float(np.sum((L.T@residual)**2));gap=abs(Q-(r['objective'][0]-r['objective'][1]-r['objective'][2]));qgaps.append(gap);assert gap<1e-6
            assert np.all(arr[:,2]>0) and np.all(np.isfinite(arr))
            assert arr[:,1].min()>=-20 and arr[:,1].max()<=70
            assert np.all((arr[:,8]>=arr[:,11])&(arr[:,8]<=arr[:,12]))
    def final(arm,cid):return [r for r in states[arm] if r['CID']==cid][-1]
    def identity(r):return [(x['epoch'],x['band'],x['values'][0]) for x in r['rows']]
    def model(r):return np.array([x['values'][2] for x in r['rows']])
    ids=['2004ef','2005hc'];starts=[];scaling=[];fixed=[]
    for cid in ids:
        base=final('nominal',cid);D=base['objective'][3]
        row=np.array([x['values'] for x in base['rows']]);W=np.array(base['weights']);m=row[:,2];y=row[:,4]
        amplitude=(m@W@y)/(m@W@m);delta=-2.5*np.log10(amplitude);assert abs(delta)<=.001
        rr=[r for r in states['nominal'] if r['CID']==cid];last=rr[-1];previous=[r for r in rr if r['iteration']<last['iteration']][-1]
        iteration_delta=last['objective'][3]-previous['objective'][3];assert abs(iteration_delta)<=.001 and identity(last)==identity(previous)
        fixed.append(dict(CID=cid,fixed_C_D_gap=float(delta),last_iteration_D_change=iteration_delta))
        en=[r for r in entries['nominal'] if r[0]==cid and r[1]=='1'];assert len(en)==1;baseentry=float(en[0][5])
        for arm,offset in [('start_minus',-.2),('start_plus',.2)]:
            r=final(arm,cid);e=[r for r in entries[arm] if r[0]==cid and r[1]=='1'];assert len(e)==1
            shift=float(e[0][5])-baseentry;dd=r['objective'][3]-D
            assert abs(shift-offset)<1e-10 and abs(dd)<=.001 and identity(r)==identity(base)
            starts.append(dict(CID=cid,arm=arm,entry_shift=shift,final_D_change=dd))
        center=final(cid+'_center',cid)
        for suffix in ['minus','plus']:
            r=final(cid+'_'+suffix,cid);assert identity(r)==identity(center)
            actual_delta=r['objective'][3]-center['objective'][3];scale=10**(-.4*actual_delta)
            gap=float(np.max(abs(model(r)/(model(center)*scale)-1)))
            assert gap<1e-6
            scaling.append(dict(CID=cid,arm=suffix,actual_delta=actual_delta,max_relative_gap=gap))
    original=BASE/'fits/nominal';current=fits/'nominal'
    def science(path):return '\n'.join(l for l in path.read_text().splitlines() if l.startswith('VARNAMES:') or l.startswith('SN:') and (not a.changed or l.split()[1]=='2005hc'))
    defaultscience=science(original/'fit.FITRES.TEXT')==science(current/'fit.FITRES.TEXT')
    def plot(path):return '\n'.join(l for l in path.read_text().splitlines() if not a.changed or l.startswith('VARNAMES:') or l.startswith('OBS:') and l.split()[1]=='2005hc')
    defaultplot=plot(original/'fit.LCPLOT.TEXT')==plot(current/'fit.LCPLOT.TEXT')
    assert defaultscience and defaultplot
    physical_checks=[];responses=[]
    if a.changed:
        old,_=parse(I/'probes-v3/fits/nominal/fit.log')
        ledger=BASE/'cohort-preparation/changed-row-ledger.csv'
        changed={(r['CID'],float(np.float32(r['MJD'])),r['old_band']):r['new_band'] for r in csv.DictReader(ledger.open()) if r['scope']=='pilot'}
        assert len(changed)==6
        def physical(r,relabel):
            result=[]
            for row in r['rows']:
                v=row['values'];band=row['band']
                if relabel:band=changed.get((r['CID'],v[0],band),band)
                result.append((band,v[0],v[4],v[5]))
            return Counter(result)
        assert len(old)==len(states['nominal'])==6
        for before,after in zip(old,states['nominal']):
            assert (before['CID'],before['iteration'])==(after['CID'],after['iteration'])
            equal=physical(before,True)==physical(after,False);assert equal
            physical_checks.append(dict(CID=before['CID'],iteration=before['iteration'],n=before['n'],physical_multiset_equal_after_declared_relabel=equal))
        for cid in ids:
            before=[r for r in old if r['CID']==cid][-1];after=final('nominal',cid)
            responses.append(dict(CID=cid,delta_D=after['objective'][3]-before['objective'][3]))
    result=dict(pass_nominal_numerical_gate=True,scope='Nominal pilot numerical gates only. Changed-pilot and full42 archive/state/mask gates still required; no physical goodness-of-fit validation.',native_binary_sha256=sha(I/'SNANA-v11_04k-output/bin/snlc_fit.exe'),default_science_equal=defaultscience,default_LCPLOT_equal=defaultplot,callbacks=sum(len(r) for r in states.values()),max_Q_gap=max(qgaps),fixed_state=fixed,distinct_starts=starts,mean_scaling=scaling,reviewed_log_hashes=hashes,source_entry_proof='Agent source excerpts: snana.car FITPAR_PREP→MNFIT_DRIVER unchangedINIVAL→MNPARM→MinuitU. Exported entries test actual post-initialization values.',still_pending='Same gates for changedpilot, then all42 unchanged nativebaseline/copy andchangedacceptedmask/support/state verification.')
    if a.changed:
        result.pop('pass_nominal_numerical_gate')
        result.update(pass_changed_pilot_numerical_gate=True,scope='Changed two-object pilot numerical gates only. Default equality refers to unchanged control 2005hc; no full42 release or physical goodness-of-fit validation.',physical_row_checks=physical_checks,pilot_processing_response=responses,still_pending='Full42 baseline has already failed separate frozen archived-Q and convergence gates. No full42 changed-response release follows from this pilot.')
    (out/'result.json').write_text(json.dumps(result,indent=2)+'\n');(out/'executed-source.py').write_bytes(Path(__file__).read_bytes());(out/'parser-source.py').write_bytes((ROOT/'scripts/research_2026_09_26/review_csp_native_states.py').read_bytes())
    (out/'manifest.json').write_text(json.dumps({p.name:sha(p) for p in out.iterdir() if p.is_file() and p.name!='manifest.json'},indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='reviewed_log_hashes'},indent=2))
if __name__=='__main__':main()
