"""Restore the exact input tree for the bounded BayeSN held-out experiment.

Downloads remain outside Git. This preserves historical fixtures and pinned
author assets; it neither parses measured flux nor creates a Python environment.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import urllib.parse
import urllib.request

ROOT = Path(__file__).resolve().parents[3]
HISTORY = '17487bf659fcbdeeea072221492bac14b04a0a85'
BAYESN = '08eef9e54188f601506ef9f7f47fc65f9f52a946'
RAISIN = 'a383c4bd03c9fbfd64bf5bfda525aec38d32b39c'
GATE = ROOT/'studies/infrared/results/bayesn-nir-forward-gate.json'
GATE_SHA256 = '21d4bcfdf3ef63be82e363466313dc1f038ff2e7d69609a63fc71044d3e830f2'


def digest(data):
    return hashlib.sha256(data).hexdigest()


def source_url(name):
    bayesn = 'runs/research_2026_09_26/bayesn_distance_identification/official-code/'
    raisin = 'sources/repos/djones1040__RAISIN_DataRelease@a383c4b/'
    if name.startswith(bayesn):
        repository, commit, path = 'bayesn/bayesn', BAYESN, name[len(bayesn):]
        kind = 'pinned_author_asset'
    elif name.startswith(raisin):
        repository, commit, path = 'djones1040/RAISIN_DataRelease', RAISIN, name[len(raisin):]
        kind = 'pinned_author_asset'
    else:
        repository, commit, path = 'sprajs/cosmology-revalidation', HISTORY, name
        kind = 'historical_execution_fixture'
    return f'https://raw.githubusercontent.com/{repository}/{commit}/{urllib.parse.quote(path,safe="/")}', kind


def restore(destination, prefer_local_git):
    assert digest(GATE.read_bytes())==GATE_SHA256
    inputs = json.loads(GATE.read_text())['execution_design']['input_sha256']
    assert not destination.exists(), 'Use a fresh tree; never overwrite a prior restoration.'
    destination.mkdir(parents=True)
    records = []
    plan = dict(source_sha256=digest(Path(__file__).read_bytes()), parent_gate_sha256=GATE_SHA256,
                destination=str(destination), prefer_local_git=prefer_local_git,
                scope='Exact byte restoration only. No measured-flux parsing, model evaluation, environment installation or inference.',
                expected_files=len(inputs), files=records)
    try:
        for name, expected in inputs.items():
            relative = Path(name)
            assert not relative.is_absolute() and '..' not in relative.parts
            url, kind = source_url(name)
            data = None; method = 'https'
            if prefer_local_git:
                result = subprocess.run(['git','cat-file','blob',f'{HISTORY}:{name}'],cwd=ROOT,capture_output=True)
                if result.returncode==0:
                    data = result.stdout; method = 'local_historical_git_blob'
            if data is None:
                request = urllib.request.Request(url,headers={'User-Agent':'cosmology-revalidation-input-restoration'})
                with urllib.request.urlopen(request,timeout=60) as response:
                    data = response.read()
            actual = digest(data)
            assert actual==expected, f'Identity mismatch; no substitute accepted: {name}'
            output = destination/relative
            output.parent.mkdir(parents=True,exist_ok=True)
            with output.open('xb') as stream:
                stream.write(data)
            records.append(dict(path=name,sha256=actual,bytes=len(data),url=url,source_kind=kind,method=method))
        assert len(records)==187
        plan['status']='passed_exact_input_restoration'
    except BaseException as error:
        plan['status']='failed_partial_restoration_preserved'; plan['error']=repr(error)
        raise
    finally:
        (destination/'restoration.json').write_text(json.dumps(plan,indent=2,allow_nan=False)+'\n')
    return plan


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--destination',type=Path,required=True)
    parser.add_argument('--prefer-local-git',action='store_true',help='Reuse historical Git blobs when present; fetch absent assets from pinned URLs.')
    args=parser.parse_args()
    result=restore(args.destination.resolve(),args.prefer_local_git)
    print(json.dumps({key:value for key,value in result.items() if key!='files'},indent=2))
