"""Independent raw-table/photometry verification; no native fits or imports."""
from pathlib import Path
from collections import Counter
import csv
import hashlib
import json
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / 'runs/research_2026_09_26/raisin_nir_timing_sensitivity'
OUT = RUN / 'root-independent-review'
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()


def table(path, prefix):
    cols = None; rows = []
    for line in path.read_text().splitlines():
        bits = line.split()
        if bits and bits[0] == 'VARNAMES:': cols = bits[1:]
        if bits and bits[0] == prefix:
            assert cols is not None and len(bits)-1 == len(cols)
            rows.append(dict(zip(cols, bits[1:])))
    return rows


def epochs(path):
    grouped = {}
    for row in table(path, 'OBS:'):
        if row['DATAFLAG'] != '1': continue
        assert row['BAND'] in 'JH'
        grouped.setdefault(row['CID'], Counter()).update([
            (row['BAND'], row['MJD'], row['FLUXCAL'], row['FLUXCAL_ERR'])])
    return grouped


def main():
    OUT.mkdir(exist_ok=True)
    protocol=json.loads((RUN/'protocol-shortvpec.json').read_text())
    reported=json.loads((RUN/'result.json').read_text())
    hashes={**protocol['inputs_sha256'],**reported['inputs_sha256']}
    for p,h in hashes.items(): assert sha(ROOT/p)==h,p
    baseline={r['CID']:r for r in table(RUN/'fits-shortvpec/baseline/fit.FITRES.TEXT','SN:')}
    bm=epochs(RUN/'fits-shortvpec/baseline/fit.LCPLOT.TEXT')
    cids=protocol['cohort']; assert set(baseline)==set(cids)
    results={}; row_ledger=[]; input_checks=0; parsed=0; mask_checks=0
    for job in protocol['jobs']:
        label=job['condition']; folder=(ROOT/job['nml']).parent
        data=ROOT/job['private_data_path']/'DES_RAISIN'
        fits={r['CID']:r for r in table(folder/'fit.FITRES.TEXT','SN:')}
        masks=epochs(folder/'fit.LCPLOT.TEXT')
        assert set(fits)==set(cids)
        deltas=[]
        for cid in cids:
            p=data/(cid+'.snana.dat'); nominal=RUN/'data/baseline/DES_RAISIN'/(cid+'.snana.dat')
            raw=p.read_text().splitlines(); base=nominal.read_text().splitlines()
            ispeak=lambda s:s.strip().startswith('PEAKMJD:')
            assert [s for s in raw if not ispeak(s)] == [s for s in base if not ispeak(s)]
            peaks=[float(s.split()[1]) for s in raw if ispeak(s)]
            assert len(peaks)==2 and peaks[0]==peaks[1]
            expected=float(f'{np.float32(peaks[0]):.4f}')
            r=fits[cid]; parsed+=1;input_checks+=1
            for key,value in [('PKMJD',expected),('PKMJDINI',expected),('PKMJDERR',0.),('STRETCH',1.),('STRETCHERR',0.),('AV',0.),('AVERR',0.),('RV',1.518),('RVERR',0.),('ERRFLAG_FIT',0.)]:
                assert float(r[key])==value,(cid,label,key,r[key],value)
            assert masks[cid]==bm[cid];mask_checks+=1
            assert float(r['NDOF'])==sum(masks[cid].values())-1
            delta=float(r['DLMAG'])-float(baseline[cid]['DLMAG']);deltas.append(delta)
            row_ledger.append({'condition':label,'CID':cid,'DLMAG':float(r['DLMAG']),'delta_DLMAG':delta,'accepted_epochs':sum(masks[cid].values())})
        results[label]={'mean':float(np.mean(deltas)),'median':float(np.median(deltas)),
                        'minimum':float(min(deltas)),'maximum':float(max(deltas)),
                        'max_abs':float(max(abs(x) for x in deltas))}
        if label in reported['summary']:
            rr=reported['summary'][label]
            for ours,theirs in [('mean','mean_delta_DLMAG'),('median','median_delta_DLMAG'),('minimum','min_delta_DLMAG'),('maximum','max_delta_DLMAG'),('max_abs','max_abs_delta_DLMAG')]:
                assert abs(results[label][ours]-rr[theirs])<1e-13
    assert results['baseline']['max_abs']==results['baseline_copy']['max_abs']==0
    derivatives=[]
    for cid in cids:
        values={r['condition']:r['delta_DLMAG'] for r in row_ledger if r['CID']==cid}
        derivatives.append({'CID':cid,'central_slope_halfday_mag_per_day':values['plus0p5']-values['minus0p5'],
            'central_slope_oneday_mag_per_day':(values['plus1']-values['minus1'])/2,
            'second_difference_halfday_mag_per_day2':(values['plus0p5']+values['minus0p5'])/.25,
            'second_difference_oneday_mag_per_day2':values['plus1']+values['minus1']})
    for name,rows in [('fit-ledger.csv',row_ledger),('timing-response-differences.csv',derivatives)]:
        with (OUT/name).open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    result={'passed':True,'verified_hash_entries':len(hashes),'parsed_native_records':parsed,
        'input_only_peak_change_checks':input_checks,'unchanged_mask_checks':mask_checks,'summaries':results,
        'finite_difference_scope':'Descriptive conditional response. No timing-error distribution or new variance/bias is assigned.',
        'protocol_sha256':sha(RUN/'protocol-shortvpec.json'),'reviewer_source_sha256':sha(Path(__file__))}
    (OUT/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    (OUT/'executed_source.py').write_bytes(Path(__file__).read_bytes())
    (OUT/'manifest.json').write_text(json.dumps({'files_sha256':{str(p.relative_to(OUT)):sha(p) for p in OUT.iterdir() if p.is_file() and p.name!='manifest.json'}},indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__': main()
