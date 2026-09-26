"""Verify and package the completed cross-model audit, not its astrophysics.

This is a record-integrity/completeness check. Scientific claims depend on the
named numerical tests and primary sources, not on passing this inventory.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import re

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'runs/assumption_audit'
DOC=ROOT/'docs/assumption-audit'
SCRIPT=ROOT/'scripts/assumption_audit'


def strict(path):
    def reject(x):raise ValueError(f'Nonfinite JSON constant {x} in {path}')
    return json.loads(path.read_text(),parse_constant=reject)


def digest(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def main():
    required={
        'host_age_original_source_and_all_archived_draws':[
            'docs/assumption-audit/pantheon.md',
            'runs/assumption_audit/pantheon/age-semantics.json',
            'runs/assumption_audit/pantheon/host-dust-contract.json',
            'runs/assumption_audit/pantheon/r19-global-validated-age-summary.json',
            'runs/assumption_audit/pantheon/r19-global-prior-violations.json',
            'runs/assumption_audit/pantheon/r19-global-validated-age-audit.csv',
            'runs/assumption_audit/pantheon/r19-sample-crosswalk.json',
            'runs/assumption_audit/pantheon/r19-age-audit.png',
            'runs/assumption_audit/des/r19-prior-peer-review.json',
            'runs/assumption_audit/age-independent-crosscheck.json'],
        'pantheon_covariance_foreground_and_shared_uncertainty':[
            'docs/assumption-audit/foreground-and-covariance.md',
            'runs/assumption_audit/covariance-groups.json',
            'runs/assumption_audit/bao/group-peer-review.json',
            'runs/assumption_audit/foreground-map-summary.json',
            'runs/assumption_audit/foreground-peer-review.json',
            'runs/assumption_audit/shared-uncertainty.json'],
        'des_dust_calibration_age_and_error_components':[
            'docs/assumption-audit/des.md',
            'runs/assumption_audit/des/covariance_results.json',
            'runs/assumption_audit/des/release_checks.json'],
        'bao_joint_probe_and_cmb_assumptions':[
            'docs/assumption-audit/bao.md',
            'runs/assumption_audit/bao/release_check.json',
            'runs/assumption_audit/bao/verification.json',
            'runs/assumption_audit/bao/wa-10-seed210921/summary.json',
            'runs/assumption_audit/bao/wa-10-seed210922/summary.json',
            'runs/assumption_audit/bao/joint-wa-10-seed210926/summary.json',
            'runs/assumption_audit/bao/joint-wa-3-seed210927/summary.json'],
        'cross_model_synthesis_and_scientific_limits':['docs/assumption-audit/README.md','docs/assumption-audit/HANDOFF.md'],
    }
    for paths in required.values():
        for name in paths:
            p=ROOT/name
            assert p.is_file() and p.stat().st_size>0,name
            if p.suffix=='.json':strict(p)
    age=strict(OUT/'pantheon/r19-global-validated-age-summary.json')
    detail=strict(OUT/'pantheon/r19-global-prior-violations.json')
    assert age['hosts']==len(detail)==103
    assert len({v['CID'] for v in detail})==103
    assert age['draws_total']==sum(v['draws_total'] for v in detail)
    assert age['draws_valid']+age['draws_excluded']==age['draws_total']
    assert age['draws_valid']==sum(v['draws_used'] for v in detail)
    assert age['draws_excluded']==sum(v['draws_excluded'] for v in detail)
    # These checks ensure references point to actual evidence, not just a plan.
    assert strict(OUT/'age-independent-crosscheck.json')['max_abs_difference_Gyr']<1e-8
    assert strict(OUT/'covariance-groups.json')['baseline_crosscheck_max_mag2']<2e-7
    assert strict(OUT/'shared-uncertainty.json')['checks']['explicit_amplitude_minimum_vs_marginal_quadratic_max']<1e-7
    assert strict(OUT/'des/release_checks.json')['calibration_overlap']['relative_stat_whitened_frobenius_remainder']<2e-6
    assert strict(OUT/'bao/verification.json')['cpl_relevant_ast_identical_after_removing_qbins_branches']
    links=[]
    for p in DOC.glob('*.md'):
        for link in re.findall(r'\]\(([^)]+)\)',p.read_text()):
            if '://' in link or link.startswith('#'):continue
            link=link.split('#')[0]
            target=(p.parent/link).resolve()
            assert target.exists(),(p,link)
            links.append(str(target.relative_to(ROOT)))
    # All producing-code/input/output hash maps retained by the investigators.
    # Older source versions may be resolved ONLY by an exact preserved copy.
    cache={}
    def cached(p):
        key=(str(p),p.stat().st_size,p.stat().st_mtime_ns)
        if key not in cache:cache[key]=digest(p)
        return cache[key]
    maps={'inputs','input_sha256','inputs_sha256','outputs','outputs_sha256',
          'sources_sha256','executed_script_hashes_captured_at_start'}
    verified=[];historical=[]
    manifests=sorted(p for p in OUT.rglob('*.json') if
        'manifest' in p.name and p.name!='final-manifest.json' and
        not any(x in p.parts for x in ['preliminary','sources','primary']))
    # Include pointwise independent reviews that carry their own input digests.
    manifests += [OUT/'age-independent-crosscheck.json',OUT/'bao/group-peer-review.json']
    for p in manifests:
        content=strict(p)
        hashmaps=[content[key] for key in maps & content.keys() if isinstance(content[key],dict)]
        # Some small manifests are themselves a flat path-to-digest map.
        if content and all(isinstance(v,str) and re.fullmatch('[0-9a-f]{64}',v) for v in content.values()):
            hashmaps.append(content)
        for hashmap in hashmaps:
            for path,want in hashmap.items():
                if not isinstance(want,str) or not re.fullmatch('[0-9a-f]{64}',want):continue
                candidates=[ROOT/path,p.parent/path]
                current=next((q for q in candidates if q.is_file()),None)
                if current is not None and cached(current)==want:
                    verified.append({'manifest':str(p.relative_to(ROOT)),'path':path,'sha256':want})
                    continue
                # Exact hash-addressed imported source snapshots are retained.
                snapshots=list((OUT/'bao/provenance').glob(want+'.py'))+list((ROOT/'runs/cosmology/provenance').glob(want+'.py'))
                snapshots += [q for q in (OUT/'pantheon/provenance').glob('*.py') if cached(q)==want]
                found=next((q for q in snapshots if cached(q)==want),None)
                assert found is not None,('No current or preserved exact bytes',str(p),path,want)
                historical.append({'manifest':str(p.relative_to(ROOT)),'recorded_path':path,
                                   'exact_snapshot':str(found.relative_to(ROOT)),'sha256':want})
    # Inventory small deliverables; bulk sources/chains stay local and are
    # bound by the verified upstream manifests rather than copied into Git.
    artifacts=[*DOC.glob('*.md'),*SCRIPT.glob('*.py')]
    artifacts += [p for p in OUT.rglob('*') if p.is_file() and p.suffix in ['.json','.csv','.md','.png']
                  and not any(x in p.parts for x in ['preliminary','sources','primary','provenance','code'])
                  and p.name not in ['final-manifest.json','final-verification.json']]
    for p in artifacts:
        if p.suffix=='.json':strict(p)
    inventory={str(p.relative_to(ROOT)):cached(p) for p in sorted(set(artifacts))}
    final={'created_utc':datetime.now(timezone.utc).isoformat(),
           'scope':'Audit of DES/Pantheon/BAO/CMB, extinction/dust/host ages, correlation and uncertainty assumptions outside separate SALT/SNM correction work',
           'requirements_and_evidence':required,'artifacts_sha256':inventory,
           'verification_role':'Integrity and coverage of audit record; not a proof of the physical truth of any model.',
           'bulk_data_policy':'Original data, downloaded maps/archives and posterior chains remain local. Preserve them and the referenced manifests for reproduction.',
           'concurrent_source_edit':'A separate task changed qbins-only core code during this audit. Original/current snapshots and CPL equivalence are recorded in bao/verification.json; no unrelated edits were reverted.',
           'not_claimed':['An official corrected C25 age catalogue','A measured replacement age-luminosity law','A full foreground/SN/host/BAO pipeline refit','A corrected CMB likelihood','A unique acceleration or deceleration conclusion']}
    (OUT/'final-manifest.json').write_text(json.dumps(final,indent=2,allow_nan=False)+'\n')
    check={'status':'PASS','authoritative_host_count':age['hosts'],'authoritative_archive_draw_count':age['draws_total'],
           'verified_hash_references':len(verified),'historical_exact_snapshot_references':historical,
           'checked_local_links':len(links),'inventoried_artifacts':len(inventory),
           'final_manifest_sha256':digest(OUT/'final-manifest.json'),
           'qualification':'Scientific checks have their own scopes and caveats; this check only binds their record.'}
    (OUT/'final-verification.json').write_text(json.dumps(check,indent=2,allow_nan=False)+'\n')
    print(json.dumps(check,indent=2))


if __name__=='__main__':main()
