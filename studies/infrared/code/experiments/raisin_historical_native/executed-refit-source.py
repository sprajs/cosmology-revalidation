"""Source-version-matched RAISIN native controls; preserve all fit outcomes."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import re
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'runs/research_2026_09_26/raisin_historical_native'
MODERN = ROOT / 'runs/research_2026_09_26/astra_design/raisin_signed_refit'
BUILD = ROOT / 'phase2/official/build/SNANA-v11_04k'
SOURCE = ROOT / 'sources/repos/RickKessler__SNANA@v11_04k'
REL = ROOT / 'sources/repos/djones1040__RAISIN_DataRelease@a383c4b'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def prepare():
    assert not (OUT / 'fit-protocol.json').exists()
    previous = json.loads((MODERN / 'input-manifest.json').read_text())
    for name, expected in previous['inputs_sha256'].items():
        assert sha(ROOT / name) == expected, name
    acquisition = json.loads((OUT / 'acquisition.json').read_text())
    for name, expected in acquisition['source_files_sha256'].items():
        assert sha(SOURCE / name) == expected, name
    paths = [Path(__file__), ROOT / 'scripts/research_2026_09_26/build_raisin_historical.sh',
             OUT / 'acquisition.json', OUT / 'build.log', BUILD / 'bin/snlc_fit.exe',
             REL / 'vpec/vpec_baseline_raisin.list', MODERN / 'input-manifest.json',
             MODERN / 'author-optical-FITOPT000.FITRES.gz']
    jobs = []
    for stage in ('engineering', 'engineering_copy', 'full'):
        timings = ['free'] if stage != 'full' else ['free', 'fixed']
        arms = ['R'] if stage == 'engineering' else ['Rcopy'] if stage == 'engineering_copy' else ['R', 'A', 'B', 'H']
        starts = [1.0] if stage != 'full' else [1.0, .85, 1.15]
        for timing in timings:
            for arm in arms:
                for start in starts:
                    source_stage = stage if stage != 'engineering_copy' else 'engineering_copy'
                    base = MODERN / 'fits' / source_stage / timing / arm / 'fit.nml'
                    laws = [94, 99] if stage == 'full' and timing == 'free' and arm == 'R' else [94]
                    for law in laws:
                        target = OUT / 'fits' / stage / str(law) / timing / arm / str(start)
                        target.mkdir(parents=True, exist_ok=True)
                        text = base.read_text()
                        # Inherited published NML leaves the MW law unspecified. The historical default is 94.
                        text = re.sub(r'(?m)^ OPT_MWCOLORLAW = 99\n', '' if law == 94 else ' OPT_MWCOLORLAW = 99\n', text)
                        text = re.sub(r"(?m)^\s*TEXTFILE_PREFIX\s*=.*$", " TEXTFILE_PREFIX = '" + str(target / 'fit') + "'", text)
                        text = re.sub(r"(?m)^\s*HEADER_OVERRIDE_FILE\s*=.*$", " HEADER_OVERRIDE_FILE = '" + str(REL / 'vpec/vpec_baseline_raisin.list') + "'", text)
                        text = re.sub(r'(?m)^\s*INIVAL_SHAPE\s*=.*$', ' INIVAL_SHAPE = ' + str(start), text)
                        (target / 'fit.nml').write_text(text)
                        paths.append(target / 'fit.nml')
                        jobs.append({'stage': stage, 'law': law, 'timing': timing, 'arm': arm, 'start': start,
                                     'nml': str((target / 'fit.nml').relative_to(ROOT))})
    env = os.environ.copy(); env['LD_LIBRARY_PATH'] = str(ROOT / 'phase2/official/build/sysroot/usr/lib')
    tools = {name: subprocess.run(command, env=env, capture_output=True, text=True).stdout
             for name, command in [('ldd', ['ldd', str(BUILD / 'bin/snlc_fit.exe')]), ('gfortran', [str(ROOT / 'phase2/official/build/sysroot/usr/bin/gfortran'), '--version']), ('gcc', ['gcc', '--version'])]}
    dump(OUT / 'build-environment.json', tools); paths.append(OUT / 'build-environment.json')
    protocol = {
        'status': 'Frozen before any source-version-matched fit outcome',
        'purpose': 'Reproduce author baseline with the tagged source version; determine whether modern start instability also occurs there. No global optimum or physical correction inferred from native success flags.',
        'source_tag': acquisition['tag'], 'source_commit': acquisition['commit'],
        'known_results': 'Sign omission and modern240-fit numerical instability are known; this is a registered source-version follow-up.',
        'primary_MW_law': 'Unspecified in original NML, historical default94. Old-version free-R law99 controls retain the modern explicit law for a version-vs-law comparison.',
        'data_arms': 'Existing immutable R/A/B/H files and source row ledger, with Rcopy engineering reproducibility control.',
        'fit_sequence': 'First engineering R and byte-identical Rcopy, then 240 primary fits (10SN x4arms x2timings x3starts) plus30 old-code/law99 free-R controls.',
        'engineering_gates': 'Actual released grid/KCOR initialization, one finite converged object, metadata and RV closure, identical copied-input fitted outputs to native printed precision; retain failure rather than silently retune.',
        'comparison': 'Author raw FITOPT000, all paired-distance/AV/stretch/peak changes, data FITCHI2, accepted masks and start spreads. Print rounding sets minimum resolution. No FITCHI2 ranking across evolving covariances or prior centers.',
        'environment_limit': 'Tagged source built with current local compiler/libraries; historical executable/toolchain and exact execution not reproduced.',
        'jobs': jobs, 'cohort': previous['cohort'],
        'inputs_sha256': {str(p.relative_to(ROOT)): sha(p) for p in paths},
        'parent_inputs_sha256': previous['inputs_sha256'],
    }
    dump(OUT / 'fit-protocol.json', protocol)
    (OUT / 'executed-refit-source.py').write_bytes(Path(__file__).read_bytes())
    print(json.dumps({'jobs': len(jobs), 'protocol_sha256': sha(OUT / 'fit-protocol.json')}))


def run(stage):
    protocol = json.loads((OUT / 'fit-protocol.json').read_text())
    for hashes in (protocol['inputs_sha256'], protocol['parent_inputs_sha256']):
        for name, expected in hashes.items():
            assert sha(ROOT / name) == expected, name
    env = os.environ.copy()
    env.update(SNANA_DIR=str(BUILD), SNDATA_ROOT=str(ROOT / 'phase2/official/inputs/SNDATA_ROOT'),
               LD_LIBRARY_PATH=str(ROOT / 'phase2/official/build/sysroot/usr/lib'),
               OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    for job in protocol['jobs']:
        if job['stage'] != stage:
            continue
        nml = ROOT / job['nml']; directory = nml.parent
        assert not (directory / 'fit.log').exists(), 'Never overwrite executed fit'
        start = time.monotonic()
        with (directory / 'fit.log').open('x') as handle:
            process = subprocess.run([str(BUILD / 'bin/snlc_fit.exe'), str(nml)], cwd=directory,
                                     env=env, stdout=handle, stderr=subprocess.STDOUT, timeout=900)
        result = {'returncode': process.returncode, 'seconds': time.monotonic() - start,
                  'protocol_sha256': sha(OUT / 'fit-protocol.json'), 'log_sha256': sha(directory / 'fit.log'), **job}
        dump(directory / 'execution.json', result)
        print(json.dumps(result), flush=True)
        assert process.returncode == 0, 'Preserve failure and diagnose before retry'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('action', choices=['prepare', 'run']); parser.add_argument('--stage', default='engineering')
    args = parser.parse_args()
    prepare() if args.action == 'prepare' else run(args.stage)
