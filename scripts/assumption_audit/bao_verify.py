"""Verify BAO audit records and scope a concurrent core edit precisely.

This verifies bytes and AST equivalence, not the physical truth of a model.
Historical manifests are preserved, including the documented capture-timing limit.
"""
from pathlib import Path
from datetime import datetime, timezone
import ast
import copy
import hashlib
import json

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'runs/assumption_audit/bao'
OLD = '18cc923c532e96791492ae2f72c4c8a8f8fc9cdfface9dfdff8f29cd4178eecc'
NEW = '83f852a3d51e851ed9030ee8f9e5b6c74fae5893bb5c12891a3fe1a8eaa1bc8c'


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def normalize(tree):
    tree = copy.deepcopy(tree)
    tree.body = [n for n in tree.body if not
                 (isinstance(n, ast.FunctionDef) and n.name == '_check_qbin_domain')]
    for n in tree.body:
        if isinstance(n, ast.FunctionDef) and n.name in ('efunc', 'integral', 'qvalue'):
            # Remove only branches whose exact guard is model == 'qbins'.
            n.body = [b for b in n.body if not
                      (isinstance(b, ast.If) and ast.dump(b.test) ==
                       ast.dump(ast.parse("model == 'qbins'", mode='eval').body))]
    return ast.dump(tree, include_attributes=False)


def main():
    snapshots = OUT / 'provenance'
    for p in snapshots.glob('*.py'):
        assert sha(p) == p.stem, p
    trees = [ast.parse((snapshots / (h + '.py')).read_text()) for h in (OLD, NEW)]
    exact = {}
    for name in ('Pantheon', 'BAO', 'mu', 'interpolation_basis'):
        nodes = [next(n for n in t.body if getattr(n, 'name', None) == name) for t in trees]
        exact[name] = ast.dump(nodes[0]) == ast.dump(nodes[1])
        assert exact[name], name
    unchanged = normalize(trees[0]) == normalize(trees[1])
    assert unchanged
    verified = []; archived = []
    for manifest in sorted(OUT.glob('**/manifest.json')):
        record = json.loads(manifest.read_text())
        sh = record['script_sha256']
        assert sha(snapshots / (sh + '.py')) == sh, manifest
        for key in ('inputs_sha256', 'outputs_sha256'):
            for name, want in record[key].items():
                candidates = [ROOT / name, manifest.parent / name]
                actual = next((p for p in candidates if p.is_file() and sha(p) == want), None)
                if actual is None:
                    snapshot = snapshots / (want + '.py')
                    assert snapshot.is_file() and sha(snapshot) == want, (manifest, name, want)
                    archived.append(dict(manifest=str(manifest.relative_to(ROOT)),
                                         recorded_path=name, exact_snapshot=str(snapshot.relative_to(ROOT)),
                                         sha256=want))
        verified.append(str(manifest.relative_to(ROOT)))
    # Byte identity is deliberately not used to infer when a running process imported code.
    record = dict(
        created_utc=datetime.now(timezone.utc).isoformat(),
        verifier_sha256=sha(Path(__file__)),
        verified_prior_run_manifests=[p for p in verified if p != str((OUT/'manifest.json').relative_to(ROOT))],
        release_audit_manifest_verified=True,
        historical_exact_snapshot_references=archived,
        old_core_sha256=OLD, current_core_sha256=sha(ROOT/'scripts/cosmology/core.py'),
        compared_core_sha256=[OLD, NEW],
        whole_module_ast_identical_except_qbins_branches_and_helper=unchanged,
        cpl_relevant_ast_identical_after_removing_qbins_branches=unchanged,
        exact_ast_matches=exact,
        cpl_scope='The two whole modules are identical after removing only exact model==qbins '
                  'branches in efunc, integral, qvalue, and the added _check_qbin_domain helper.',
        provenance_timing_limit='Initial audit scripts hashed input code at completion. '
            'Two long joint runs overlap a concurrent qbins-only edit and their manifest '
            'records the later whole-file hash. Old and new snapshots are retained; '
            'the CPL-relevant AST is identical. Later script revisions capture imported code bytes.',
        runs_with_completion_time_core_hash=['joint-wa-10-seed210926', 'joint-wa-3-seed210927'],
        report_sha256=sha(ROOT/'docs/assumption-audit/bao.md'))
    (OUT/'verification.json').write_text(json.dumps(record, indent=2)+'\n')
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
