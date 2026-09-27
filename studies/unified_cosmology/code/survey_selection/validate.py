#!/usr/bin/env python3
"""Verify frozen survey inputs, identities and the declared scientific boundaries."""
import functools,json
from pathlib import Path
import numpy as np
import pandas as pd
from common import ROOT,HERE,WORK,RESULTS,sha

@functools.lru_cache(None)
def digest(path):return sha(path)

def main():
    producers={'acquisition':'acquire','followup-acquisition':'acquire_followup','processing-source':'acquire_processing','dovekie-interface':'dovekie','observed-survey-checks':'observations','ozdes-crossmatch':'ozdes','native-assets':'assets','calibration-registry':'calibration','classifier-closure':'classifier_closure','classifier-full':'classifier_full','classifier-diagnosis':'diagnose_classifier','dataprep':'dataprep','des3yr-acquisition':'des3yr','des3yr-support':'des3yr_checks','des3yr-interface':'des3yr_likelihood','geometry-review':'review_geometry','prior-classifier-audit':'prior_classifier_audit'}
    count=0;records={};failures=[]
    def check(p,expected):
        nonlocal count
        p=Path(p);p=p if p.is_absolute() else ROOT/p
        count+=1
        if not p.exists():failures.append({'path':str(p),'error':'missing'})
        elif digest(str(p))!=expected:failures.append({'path':str(p),'error':'hash mismatch'})
    def walk(node):
        if isinstance(node,list):
            for v in node:walk(v)
        if isinstance(node,dict):
            if 'path' in node and 'sha256' in node:check(node['path'],node['sha256'])
            for k,v in node.items():
                if k in ['inputs','outputs','output','dependencies','source_sha256','reviewed_source_sha256','native_object_sha256'] and isinstance(v,dict):
                    for p,h in v.items():
                        if isinstance(h,str) and len(h)==64:check(p,h)
                walk(v)
    for label,producer in producers.items():
        p=RESULTS/(label+'.json');d=json.loads(p.read_text());records[label]=d;check(HERE/(producer+'.py'),d['code_sha256']);walk(d)
    for r in records['acquisition']['new_bundle']['extracted_assets']+records['native-assets']['new_assets']:
        check(WORK/'SNDATA_ROOT_2026-04-10'/r['member'],r['sha256'])
    full=records['classifier-full'];computed={}
    for name,r in full['results'].items():
        p=ROOT/r['path'];check(p,r['rows_sha256']);f=pd.read_csv(p,dtype={'CID':str});assert len(f)==17733 and f.CID.is_unique
        assert f.status.eq('predicted').all();error=(f.pIa_reconstructed-f.pIa_released).abs()
        countclose=int(error.le(5.1e-5).sum());gate=int((f.pIa_reconstructed.gt(.5)!=f.pIa_released.gt(.5)).sum())
        assert countclose==r['all']['within_4decimal_rounding'] and gate==r['all']['gate_0_5_disagreement']
        computed[name]={'rows':len(f),'within_rounding_and_float32_allowance':countclose,'classification_gate_disagreements':gate}
    assert computed['native_OPT_SETPKMJD16_PKMJDINI']['within_rounding_and_float32_allowance']==17733
    assert computed['public_HEAD_PEAKMJD']['classification_gate_disagreements']==35
    plan=json.loads((WORK/'classifier-full-design.json').read_text());assert plan['native_source_sha256']==sha(WORK/'native-dataprep/dataprep.SNANA.TEXT')
    cross=pd.read_csv(WORK/'followup/dovekie-ozdes-clean.csv');assert len(cross)==1089 and cross.oz_OzDES_ID.nunique()==1088
    assert cross.nearest_arcsec.le(1).all() and cross.second_arcsec.gt(1).all() and cross.redshift_delta.abs().le(.003).all()
    prior=records['prior-classifier-audit'];assert prior['campaigns']==41 and prior['classified_occurrences']==55431
    assert all(c['stored_peak_vs_native_PKMJDINI_max_days']<1e-9 and not c['uses_simulation_truth_peak_bit2048'] for c in prior['checks'])
    assert records['dovekie-interface']['rows']==1820 and records['observed-survey-checks']['physical_ID_join']['detected_without_SMP']==11930
    evidence=records['des3yr-support']['spectroscopic_evidence']
    assert evidence['DES_with_at_least_one_SNIa_label']==185 and evidence['DES_with_only_provisional_SNIa_question_mark_labels']==22
    for name,sample in records['des3yr-interface']['samples'].items():
        assert sample['likelihood_rows']==20 and sample['physical_SNe_represented']==(329 if name=='combined' else 207)
        matrices={}
        for typ,entry in sample['assembly'].items():
            check(entry['NPZ'],entry['sha256']);a=np.load(ROOT/entry['NPZ']);matrices[typ]=a['covariance']
            assert len(a['CID'])==20 and len(np.unique(a['CID']))==20
            assert np.max(abs(a['covariance']@a['precision']-np.eye(20)))<1e-9
            np.linalg.cholesky(a['covariance'])
        b=WORK/'des3yr/05-COSMOLOGY/COSMOLOGY_INPUTS';label='DES+LOWz' if name=='combined' else 'DESonly'
        raw=np.loadtxt(b/f'lcparam_{label}.txt');systematic=np.loadtxt(b/f'sys_{label}_ALLSYS.txt')[1:].reshape(20,20)
        assert np.array_equal(matrices['stat'],np.diag(raw[:,5]**2))
        assert np.allclose(matrices['total'],matrices['stat']+systematic,atol=1e-12,rtol=0)
    docs=[HERE/'README.md',ROOT/'studies/unified_cosmology/notes/survey-selection.md']
    result={'status':'passed' if not failures else 'failed','code_sha256':sha(__file__),'file_hash_checks':count,'unique_files_hashed':digest.cache_info().currsize,'failures':failures,'classifier_recomputed':computed,'scientific_boundaries':{'released_distance_joint_likelihood':'available; conditional on released upstream calibration, bias and contamination treatment','observed_classifier_reproduction':'all 17733 at published precision plus stated float32 allowance','full_physical_joint_selection':'not closed by these checks','DES3YR_spectroscopic_alternative':'20 released likelihood bins represent 329 observed SNe; 22 of 207 DES objects have only provisional SNIa? labels; historical model conditions retained'},'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in sorted(HERE.glob('*.py'))},'documentation_sha256':{str(p.relative_to(ROOT)):sha(p) for p in docs},'result_sha256':{str(p.relative_to(ROOT)):sha(p) for p in sorted(RESULTS.glob('*.json')) if p.name!='validation.json'}}
    (RESULTS/'validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':result['status'],'hash_checks':count,'failures':failures},indent=2));assert not failures
if __name__=='__main__':main()
