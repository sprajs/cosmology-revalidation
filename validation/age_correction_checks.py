#!/usr/bin/env python3
"""Verify the executed age studies without promoting conditional tests to physics."""
from pathlib import Path
import argparse,ast,hashlib,json,re
from datetime import datetime,timezone
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'studies/host_ages'
TOPICS=['age_reconciliation','age_recovery','environment_validation','population_transport']

def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--archive',type=Path,default=Path.home()/'.local/share/cosmology-revalidation/archive/17487bf659fcbdeeea072221492bac14b04a0a85');args=ap.parse_args()
    verified={}
    def check(p,h):
        p=Path(p);actual=sha(p);assert actual==h,(str(p),h,actual);verified[str(p)]=h
    def mapping(d,base=ROOT):
        for p,h in d.items():check(base/p,h)
    a=read(BASE/'results/age_reconciliation/manifest.json')
    mapping(a['inputs'],args.archive);mapping(a['outputs'],ROOT/'.work/age-reconciliation')
    check(BASE/'code/age_reconciliation/run.py',a['code_sha256']);check(BASE/'code/age_reconciliation/design.json',a['design_sha256']);check(BASE/'results/age_reconciliation/summary.json',a['summary_sha256'])
    source=read(BASE/'results/age_reconciliation/public-input-check.json');check(BASE/'code/age_reconciliation/check_public_inputs.py',source['code_sha256'])
    aa=read(BASE/'results/age_reconciliation/summary.json')
    assert aa['original_quality_reconstruction']['quality_rows']==175 and aa['original_quality_reconstruction']['quality_young_lt4']==70
    assert aa['counts']['duplicate_G11_R19_SNe']==33
    assert aa['validation']['eiv_block_vs_conditional_nll_difference']<1e-10
    for f in ['run.json','titan-run.json']:
        d=read(BASE/'results/population_transport'/f)
        for k in ['input_sha256','code_sha256','output_sha256','additional_input_sha256']:
            if k in d:mapping(d[k])
    dp=read(BASE/'results/population_transport/validation.json');assert dp['passed'];check(BASE/'code/population_transport/validate.py',dp['validator_sha256'])
    groups=read(BASE/'results/population_transport/group-feasibility.json')
    mapping(groups['input_sha256']);mapping(groups['code_sha256'])
    assert not groups['new_age_slope_fit_performed']
    check(ROOT/'.work/population-transport/ztf-titan-tully-crosswalk.csv',groups['crosswalk_sha256'])
    assert groups['validation']['independent_haversine_separation_max_error_arcsec']<1e-7
    assert not groups['group_width_colour_age_design']['age_identifiable_after_group_intercepts_width_colour']
    for item in read(BASE/'results/population_transport/group-inputs.json')['files']:check(ROOT/item['path'],item['sha256'])
    for f in ['cosmology-recovery.json','flux-recovery.json']:
        d=read(BASE/'results/age_recovery'/f);r=d['provenance'];mapping(r['code_sha256']);mapping(r['input_sha256']);check(BASE/'code/age_recovery'/('cosmology-protocol.json' if f.startswith('cosmology') else 'protocol.json'),r['protocol_sha256'])
        check(ROOT/'.work/age-recovery'/('cosmology-draws.csv' if f.startswith('cosmology') else 'flux-event-ledger.csv'),r.get('draw_sha256',r.get('event_ledger_sha256')))
    flux=read(BASE/'results/age_recovery/flux-recovery.json');assert flux['eligible_event_draws']==1728
    assert all(x['success']==x['selected'] and x['boundaries']==0 and x['max_two_start_chi2_gap']<1e-7 for x in flux['arms'])
    response=read(BASE/'results/age_recovery/flux-redshift-response.json');check(BASE/'results/age_recovery/flux-recovery.json',response['input_record_sha256']);check(BASE/'code/age_recovery/flux_response.py',response['code_sha256'])
    for v in response['arms'].values():
        C=np.array(v['paired_noise_covariance_mag2']);assert np.max(abs(C-C.T))<1e-14 and np.linalg.eigvalsh(C).min()>-1e-14
        assert np.allclose(C/12,v['covariance_of_estimated_mean_mag2'],rtol=0,atol=1e-15)
    avail=read(BASE/'results/age_recovery/survey-availability.json');check(BASE/'code/age_recovery/survey_availability.py',avail['code_sha256']);assert avail['objects_checked']==250
    independent=read(BASE/'results/age_recovery/independent-checks.json');assert independent['passed']
    mapping(independent['input_sha256']);mapping(independent['code_sha256'])
    # Environmental records include their own source paths because archive and
    # freshly downloaded inputs reside in distinct acquisition roots.
    for p in sorted((BASE/'results/environment_validation').glob('*.json')):
        d=read(p)
        for item in d.get('inputs',[]):
            path=Path(item['path'])
            if not path.is_absolute():path=ROOT/path
            check(path,item['sha256'])
        for k in ['dependencies_sha256']:
            for name,h in d.get(k,{}).items():check((ROOT if len(Path(name).parts)>1 else BASE/'code/environment_validation')/name,h)
        if 'raisin_selected_header_sha256' in d:mapping(d['raisin_selected_header_sha256'])
        code={'dustpedia-summary.json':'analyse.py','ztf-summary.json':'ztf.py','titan-ztf-summary.json':'titan.py','titan-calibration.json':'calibrate.py','robustness.json':'robustness.py','validation.json':'validate.py','coverage.json':'coverage.py','dovekie-summary.json':'dovekie.py','coverage-validation.json':'validate_coverage.py'}.get(p.name)
        if code and 'code_sha256' in d:check(BASE/'code/environment_validation'/code,d['code_sha256'])
        if 'helper_sha256' in d:check(BASE/'code/environment_validation/analyse.py',d['helper_sha256'])
        if 'reader_sha256' in d:check(ROOT/'lib/records.py',d['reader_sha256'])
    ev=read(BASE/'results/environment_validation/validation.json')
    assert ev['status'].startswith('pass for the declared')
    mapping(ev['result_sha256'],BASE/'results/environment_validation')
    coverage=read(BASE/'results/environment_validation/coverage-validation.json')
    assert coverage['status'].startswith('pass for source identity')
    mapping(coverage['result_sha256'],BASE/'results/environment_validation')
    astro=[]
    for topic in TOPICS:
        for p in (BASE/'code'/topic).glob('*.py'):ast.parse(p.read_text());astro.append(p)
    docs=[ROOT/'README.md',ROOT/'docs/age-correction-results.md',ROOT/'docs/age-correction-evidence.md',ROOT/'docs/age-correction-plan.md',ROOT/'provenance/README.md',ROOT/'validation/README.md',BASE/'README.md']+list((BASE/'notes').glob('*-results.md'))+list((BASE/'code/environment_validation').glob('README.md'))
    linkcount=0
    for p in docs:
        for target in re.findall(r'\]\(([^)]+)\)',p.read_text()):
            if '://' in target or target.startswith('#'):continue
            target=target.split('#')[0].split()[0].strip('<>')
            assert (p.parent/target).exists(),(str(p),target);linkcount+=1
    artifact_paths=astro+docs+[p for topic in TOPICS for p in (BASE/'code'/topic).glob('*.json')]+[p for topic in TOPICS for p in (BASE/'results'/topic).glob('*.json')]
    inventory={str(p.relative_to(ROOT)):sha(p) for p in sorted(set(artifact_paths))}
    report={'schema':'age-correction-execution-validation-v1','checked_utc':datetime.now(timezone.utc).isoformat(),'status':'passed','scope':'Current code/input/output identities, independent numerical records, selected-fit gates, finite covariance structure and documentation links. This does not validate unavailable host likelihoods, survey selection, blinding provenance, physical age interpretation or high-redshift transport.','file_identities_checked':len(verified),'python_files_parsed':len(astro),'local_links_verified':linkcount,'artifact_sha256':inventory,'validator_sha256':sha(Path(__file__)),'results':{'recovered_original_quality_sample':175,'recovered_young_quality_sample':70,'paired_cosmology_realizations_per_amplitude':500,'flux_event_intervention_draws':1728,'unavailable_Dovekie_files_at_checked_endpoints':avail['objects_checked']-avail['verified'],'new_physical_cosmology_claim':False}}
    (ROOT/'validation/reports/age-correction-execution.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='artifact_sha256'},indent=2))
if __name__=='__main__':main()
