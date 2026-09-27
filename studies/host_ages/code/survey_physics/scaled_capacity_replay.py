#!/usr/bin/env python3
"""Preserve capacity failures and replay with an independently reviewed repair."""
import argparse
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import re
import numpy as np
from common import RESULTS, HERE, sha
from scaled_bbc import SCALE, BBC_EXE, execute_bbc, fitres

REPLAY_ROOT = SCALE / "capacity-replay"


def task(spec):
    label,old,inputs,success=spec
    dest=REPLAY_ROOT/label
    record=execute_bbc(inputs[0],inputs[1],dest,True)
    evidence={'label':label, 'original_log':{'path':str(old/'bbc.log'),'sha256':sha(old/'bbc.log')},
              'original_input':{'path':str(old/'bbc.input'),'sha256':sha(old/'bbc.input')},
              'original_success':success,'replay':record}
    if success:
        assert record['graceful']
        f,g=fitres(old/'bbc.FITRES'),fitres(dest/'bbc.FITRES')
        assert f.equals(g),label
        evidence['science_table_exact_numeric_identity']=True
        evidence['science_rows_byte_identity']=[x for x in (old/'bbc.FITRES').read_text().splitlines() if x.startswith('SN:')]==[x for x in (dest/'bbc.FITRES').read_text().splitlines() if x.startswith('SN:')]
        evidence['covariance_byte_identity']=(old/'bbc.COV').read_bytes()==(dest/'bbc.COV').read_bytes()
        assert evidence['science_rows_byte_identity'] and evidence['covariance_byte_identity']
    else:
        original=(old/'bbc.log').read_text()
        capacity = 'NBIN_SIGINT=100 exceeds bound MXSTORE_PULL=100' in original
        evidence['original_failure_type']='scan_capacity' if capacity else 'zero_MAD_pull'
        assert capacity or 'Invalid stdPull = 0.000000' in original
        evidence['original_failure_tail']=original[-2500:]
        assert record['graceful'] if capacity else not record['graceful'],label
        if not capacity: assert 'Invalid stdPull = 0.000000' in record['failure_tail']
    scans=[]
    for line in record['scatter_capacity_scan_diagnostics']:
        fields={k:float(v) for k,v in re.findall(r'(N|start|mean|STD|steps|result|lower_floor)=([^ ]+)',line)}
        assert all(np.isfinite(v) for v in fields.values())
        fields['requested_capacity']=max(100,int(np.ceil((fields['start']+.3)/.01))+2)
        fields['cell_context']=line.split(' N=')[0]
        scans.append(fields)
    evidence['expanded_scan_cells']=scans
    print(label,record['graceful'],len(scans),flush=True)
    return evidence


def main():
    global REPLAY_ROOT
    parser=argparse.ArgumentParser()
    parser.add_argument('--restore-fixtures',action='store_true',
                        help='Reconstruct four baseline and twenty nominal k012 pilot cases using the original executable and original seeds, in a separate fixture directory.')
    args=parser.parse_args()
    if args.restore_fixtures:
        from scaled_restore_capacity import restore
        specs=restore()
        REPLAY_ROOT=SCALE/'capacity-fixtures/repaired'
    else:
        specs=[]
        prior=json.loads((RESULTS/'scaled-bbc-k012-common-fixed.json').read_text())
        for name,r in prior['cases'].items():
            old=Path(r['folder'])
            txt=(old/'bbc.input').read_text()
            get=lambda k:re.search(r'^'+k+r'\s*=\s*(\S+)',txt,re.M).group(1)
            specs.append(('baseline-'+name,old,(get('datafile'),get('simfile_biascor')),True))
        for old in sorted((SCALE/'bootstrap/k012-common-fixed/joint').glob('r*/nominal')):
            if 'FATAL ERROR ABORT' not in (old/'bbc.log').read_text():continue
            txt=(old/'bbc.input').read_text()
            get=lambda k:re.search(r'^'+k+r'\s*=\s*(\S+)',txt,re.M).group(1)
            specs.append(('pilot-'+old.parent.name,old,(get('datafile'),get('simfile_biascor')),False))
    assert len(specs)==10
    with ThreadPoolExecutor(max_workers=3) as pool: records=list(pool.map(task,specs))
    out={'status':'passed','code_sha256':sha(__file__),'patch':{'path':str(HERE/'scaled_scatter_capacity.patch'),'sha256':sha(HERE/'scaled_scatter_capacity.patch')},
         'build_record_sha256':sha(RESULTS/'scaled-native-capacity-build.json'), 'executable':{'path':str(BBC_EXE),'sha256':sha(BBC_EXE)},
         'scope':'Only storage capacity and diagnostics changed. Original pilot has one capacity abort and five zero-MAD failures. Capacity abort repairs; zero-MAD failures persist. All four successful science tables/covariances exactly identical.',
         'unchanged_science':'same .01 step, −.3 lower limit, MAD, covariance floor, weights, interpolation and support criteria', 'records':records, 'restored_fixtures':args.restore_fixtures}
    if args.restore_fixtures:
        out['restoration_manifest']={'path':str(SCALE/'capacity-fixtures/restoration.json'),'sha256':sha(SCALE/'capacity-fixtures/restoration.json')}
        out['restoration_code_sha256']=sha(HERE/'scaled_restore_capacity.py')
    (RESULTS/('scaled-native-capacity-restoration.json' if args.restore_fixtures else 'scaled-native-capacity-replay.json')).write_text(json.dumps(out,indent=2)+'\n')

if __name__=='__main__':main()
