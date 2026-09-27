"""Build isolated high-resolution FSPS bindings and pin existing model inputs."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import urllib.request

CODE=Path(__file__).resolve().parent
ROOT=CODE.parents[3]
WORK=ROOT/'.work/unified-cosmology/calibrated-host-physics'
RESULT=ROOT/'studies/unified_cosmology/results/calibrated_host_physics'
FSPS=ROOT/'.work/physical-ages/fsps'


def sha(path):
    with Path(path).open('rb') as file:
        return hashlib.file_digest(file,'sha256').hexdigest()


def compiler():
    """Use an available compiler, or extract the pinned Arch frontend locally."""
    existing=os.environ.get('FC') or shutil.which('gfortran')
    if existing:return existing
    folder=WORK/'toolchain'; folder.mkdir(exist_ok=True)
    name='gcc-fortran-16.2.1+r23+gd564253eb6c8-1-x86_64.pkg.tar.zst'
    url='https://archive.archlinux.org/packages/g/gcc-fortran/'+name
    expected='008443629461abc92ff3ea6d363a5d951e2b314998a4165b4d6f0f7dc5a97007'
    for suffix in ['', '.sig']:
        path=folder/(name+suffix)
        if not path.exists():
            with urllib.request.urlopen(url+suffix) as response:path.write_bytes(response.read())
    package=folder/name; assert sha(package)==expected
    signature=subprocess.run(['gpgv','--keyring','/etc/pacman.d/gnupg/pubring.gpg',
        str(package)+'.sig',str(package)],capture_output=True,text=True,check=True)
    version=subprocess.check_output(['gcc','-dumpfullversion'],text=True).strip()
    assert version=='16.2.1','Provide FC for a different host toolchain.'
    executable=folder/'usr/bin/gfortran'
    if not executable.exists():subprocess.run(['bsdtar','-xf',str(package),'-C',str(folder)],check=True)
    # Relocated GCC also needs its matching host linker plugin/startup objects.
    # Only local symlinks are created; no system compiler files are changed.
    support=folder/'usr/lib/gcc/x86_64-pc-linux-gnu/16'
    for path in Path('/usr/lib/gcc/x86_64-pc-linux-gnu/16').iterdir():
        destination=support/path.name
        if not destination.exists():destination.symlink_to(path,target_is_directory=path.is_dir())
    record={'url':url,'sha256':sha(package),'signature_sha256':sha(str(package)+'.sig'),
        'signature_verification':signature.stderr,'system_gcc':version,
        'scope':'Locally extracted frontend; existing matching GCC/runtime used; no host package changes.'}
    record['host_linker_support_sha256']={str(p):sha(p) for p in
        Path('/usr/lib/gcc/x86_64-pc-linux-gnu/16').iterdir() if p.is_file()}
    (folder/'toolchain.json').write_text(json.dumps(record,indent=2)+'\n')
    return str(executable)


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--build-environment',action='store_true')
    args=parser.parse_args(); WORK.mkdir(parents=True,exist_ok=True); RESULT.mkdir(parents=True,exist_ok=True)
    design=json.loads((CODE/'design.json').read_text())
    commit=subprocess.check_output(['git','-C',str(FSPS),'rev-parse','HEAD'],text=True).strip()
    assert commit==design['library']['source_commit']
    for path,digest in design['input_sha256'].items(): assert sha(ROOT/path)==digest,path
    # Validate the original complete tracked stellar-data inventory without changing it.
    old=ROOT/'.work/physical-ages/stellar-file-hashes.json'
    original=json.loads((ROOT/'studies/host_ages/results/physical_ages/acquisition.json').read_text())
    assert sha(old)==original['stellar_inventory_sha256']
    inventory=json.loads(old.read_text())
    for relative,expected in inventory.items():assert sha(FSPS/relative)==expected,relative
    with urllib.request.urlopen('https://pypi.org/pypi/fsps/0.5.0/json') as response:
        metadata=json.load(response)
    source=next(row for row in metadata['urls'] if row['packagetype']=='sdist')
    package=WORK/source['filename']
    if not package.exists():
        with urllib.request.urlopen(source['url']) as response:
            package.write_bytes(response.read())
    assert sha(package)==source['digests']['sha256']
    record={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'source_package':{'url':source['url'],'path':str(package.relative_to(ROOT)),
                          'sha256':sha(package)},
        'FSPS_data_commit':commit,'prior_stellar_inventory_sha256':sha(old),
        'verified_prior_stellar_files':len(inventory),
        'design_sha256':sha(CODE/'design.json'),'code_sha256':sha(__file__),
        'compiler_flags':design['library']['compile_flags'],
        'primary_build_documentation':'https://python-fsps.readthedocs.io/en/stable/installation/',
        'primary_resolution_documentation':f'https://github.com/cconroy20/fsps/blob/{commit}/SPECTRA/C3K/readme.md',
        'data_files_sha256':{str(path.relative_to(ROOT)):sha(path) for path in
            sorted((FSPS/'SPECTRA/C3K/c3k_hr').glob('*')) if path.is_file()}}
    if args.build_environment:
        python=WORK/'.venv/bin/python'
        if not python.exists():
            subprocess.run(['uv','venv','--python',str(ROOT/'.venv/bin/python'),str(WORK/'.venv')],check=True)
        fc=compiler()
        env=dict(os.environ,FFLAGS=design['library']['compile_flags'],SPS_HOME=str(FSPS),FC=fc,
                 UV_CACHE_DIR=str(WORK/'uv-cache'))
        dependencies=['numpy==1.26.4','scipy==1.15.3','astropy==5.3.4','pandas==2.2.3',
                      'clarabel==0.11.1','extinction==0.4.9','requests==2.32.5']
        subprocess.run(['uv','pip','install','--python',str(python),*dependencies],check=True,env=env)
        subprocess.run(['uv','pip','install','--python',str(python),'--no-cache','--no-binary=fsps',
                        '--reinstall-package=fsps',str(package)],check=True,env=env)
        check="import fsps,json; p=fsps.StellarPopulation(zcontinuous=1); print(json.dumps({'fsps':fsps.__version__,'libraries':[x.decode() for x in p.libraries]}))"
        measured=json.loads(subprocess.check_output([str(python),'-c',check],text=True,env=env).splitlines()[-1])
        assert 'c3k_hr' in [name.lower() for name in measured['libraries']],measured
        record['built_environment']=measured
        record['compiler_path']=fc
        record['compiler_version']=subprocess.check_output([fc,'--version'],text=True).splitlines()[0]
        if (WORK/'toolchain/toolchain.json').exists():
            record['local_toolchain']=json.loads((WORK/'toolchain/toolchain.json').read_text())
        freeze=subprocess.check_output(['uv','pip','freeze','--python',str(python)],text=True)
        freeze='\n'.join('fsps==0.5.0' if line.startswith('fsps @') else line for line in freeze.splitlines())+'\n'
        (CODE/'requirements-lock.txt').write_text(freeze)
        record['requirements_sha256']=sha(CODE/'requirements-lock.txt')
        extension=list((WORK/'.venv').glob('lib/python*/site-packages/fsps/_fsps*.so'))
        assert len(extension)==1
        record['compiled_extension_sha256']=sha(extension[0])
    failures=[WORK/'acquire-initial.py',WORK/'acquire-no-fortran.log',WORK/'acquire-missing-lto-plugin.log']
    record['initial_missing_compiler_attempt_sha256']={str(p.relative_to(ROOT)):sha(p) for p in failures if p.exists()}
    (RESULT/'acquisition.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({k:v for k,v in record.items() if k!='data_files_sha256'},indent=2))


if __name__=='__main__':
    main()
