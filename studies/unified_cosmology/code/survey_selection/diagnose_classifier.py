#!/usr/bin/env python3
"""Quantify and certify the timing-input explanation after the full comparison."""
import json,sys
import numpy as np
import pandas as pd
from common import ROOT,WORK,RESULTS,sha

def main():
    source=RESULTS/'classifier-full.json';rec=json.loads(source.read_text());tables={}
    for name,r in rec['results'].items():
        p=ROOT/r['path'];assert sha(p)==r['rows_sha256'];d=pd.read_csv(p,dtype={'CID':str}).set_index('CID');assert d.index.is_unique and len(d)==17733;tables[name]=d
    baseline=tables['public_HEAD_PEAKMJD'];native=tables['native_OPT_SETPKMJD16_PKMJDINI'].loc[baseline.index]
    assert np.array_equal(baseline.pIa_released,native.pIa_released)
    assert baseline.status.eq('predicted').all() and native.status.eq('predicted').all()
    b=(baseline.pIa_reconstructed-baseline.pIa_released).abs();n=(native.pIa_reconstructed-native.pIa_released).abs()
    tol=5.1e-5;assert np.all(n<=tol)
    top=b.idxmax();delta=native.used_peak-baseline.used_peak
    sample=pd.read_csv(WORK/'classifier-closure.csv',dtype={'CID':str}).set_index('CID')
    mismatch=sample.index[(sample.pIa_reconstructed-sample.pIa_released).abs()>tol]
    result={'status':'observed_classifier_inference_reproduced_to_release_precision','code_sha256':sha(__file__),'source_result_sha256':sha(source),'all_released_objects':len(n),'nominal_rounding_tolerance':.00005,'additional_float32_allowance':.000001,'combined_tolerance':tol,'native_max_absolute_probability_difference':float(n.max()),'native_rows_outside_tolerance':int((n>tol).sum()),'baseline_rows_outside_tolerance':int((b>tol).sum()),'original_256_mismatches':len(mismatch),'original_256_mismatches_repaired':int(n.loc[mismatch].le(tol).sum()),'baseline_gate_disagreements':int((baseline.pIa_reconstructed.gt(.5)!=baseline.pIa_released.gt(.5)).sum()),'native_gate_disagreements':int((native.pIa_reconstructed.gt(.5)!=native.pIa_released.gt(.5)).sum()),'peak_change_gt_1day':int(delta.abs().gt(1).sum()),'largest_baseline_mismatch':{'CID':top,'released_probability':float(baseline.loc[top,'pIa_released']),'header_peak_probability':float(baseline.loc[top,'pIa_reconstructed']),'native_peak_probability':float(native.loc[top,'pIa_reconstructed']),'header_peak_MJD':float(baseline.loc[top,'used_peak']),'native_peak_MJD':float(native.loc[top,'used_peak']),'shift_days':float(delta.loc[top])},'interpretation':['This verifies deterministic observed-data preprocessing and inference with released weights at published probability precision.','It does not establish classifier truth accuracy, CC-prior population, retraining under new physics, or the selection-conditioned cosmological likelihood.','Peak estimates were generated from published OPT_SETPKMJD16 and PHOTFLAG1016 settings before comparison; no peak was optimized against a released probability.'],'inputs':{str((WORK/'classifier-full-design.json').relative_to(ROOT)):sha(WORK/'classifier-full-design.json'),str((WORK/'classifier-closure.csv').relative_to(ROOT)):sha(WORK/'classifier-closure.csv')}}
    (RESULTS/'classifier-diagnosis.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
