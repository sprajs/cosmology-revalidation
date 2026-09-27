#!/usr/bin/env python3
"""Build a separately identified BBC executable with scan-storage-only repair.

The configured native tree is a dependency from recover.py/native_checks.py;
its simulation binaries and all previous executable identities remain intact.
"""
import json
from pathlib import Path
import shutil
import subprocess
from common import HERE, WORK, RESULTS, sha, native_env

BASE = WORK / 'SNANA-legacy-eff'
DEST = WORK / 'scaled-bbc/native-capacity'


def main():
    assert not DEST.exists(), 'Refuse to overwrite an existing build identity'
    shutil.copytree(BASE, DEST, ignore=shutil.ignore_patterns('.git'))
    for p in DEST.rglob('Makefile'):
        p.write_text(p.read_text().replace(str(BASE), str(DEST)))
    patch = HERE / 'scaled_scatter_capacity.patch'
    subprocess.run(['patch', '-p1', '-i', str(patch)], cwd=DEST, check=True)
    log = DEST / 'capacity-build.log'
    with log.open('w') as stream:
        subprocess.run(['make', '-C', 'src', 'SALT2mu', '-j2'], cwd=DEST,
                       env=native_env(), stdout=stream, stderr=subprocess.STDOUT, check=True)
    result = {
        'status':'built',
        'scope':'Storage capacity and diagnostics only; unchanged scan step, lower floor, objective, interpolation and native scientific support gates.',
        'patch':{'path':str(patch),'sha256':sha(patch)},
        'previous_high_snr_capacity_patch':{'path':str(HERE/'bbc_bounds.patch'),'sha256':sha(HERE/'bbc_bounds.patch')},
        'code_sha256':sha(__file__),
        'build_log':{'path':str(log),'sha256':sha(log)},
        'original_executable':{'path':str(BASE/'bin/SALT2mu.exe'),'sha256':sha(BASE/'bin/SALT2mu.exe')},
        'executable':{'path':str(DEST/'bin/SALT2mu.exe'),'sha256':sha(DEST/'bin/SALT2mu.exe')},
        'source':{str(p):sha(p) for p in [BASE/'src/sntools.c', DEST/'src/sntools.c', BASE/'src/SALT2mu.c',DEST/'src/SALT2mu.c']},
        'unchanged_SALT2mu_source':sha(BASE/'src/SALT2mu.c')==sha(DEST/'src/SALT2mu.c'),
    }
    (RESULTS/'scaled-native-capacity-build.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__': main()
