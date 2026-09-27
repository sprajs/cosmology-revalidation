"""Build isolated matched CAMB1.6.6 spline controls; no cosmological calls."""
import argparse
import base64
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile
import urllib.request
import zipfile

ROOT = Path(__file__).resolve().parents[4]
WORK = ROOT / '.work/unified-cosmology/camb-spline-audit'
REFERENCE = ROOT / '.work/unified-cosmology/external-probes/.modern-venv'
ASSETS = {
    'camb-1.6.6-py3-none-manylinux_2_27_x86_64.manylinux_2_28_x86_64.whl':
        ('pypi', 'bb8e0c85b3e33d78ca0676c84b0b60b54dfa72d4e520a805166cab0458e81094'),
    'camb-1.6.6.tar.gz': ('pypi', '9856202a5c05570256e52377b20431891c7b08b2e9c334e141fd08d2a085516f'),
    'camb-1.6.6-tag.tar.gz': ('https://codeload.github.com/cmbant/CAMB/tar.gz/3ef0272d6f7ba1231128872e56e6d4c12af8267b',
        '7c7c8c3e1161d42f193eaf9690e818b54e01004f65342449042267055dfd2de0'),
    'forutils-841f06.tar.gz': ('https://codeload.github.com/cmbant/forutils/tar.gz/841f06d5356877e90437f94b2d976dd98f7f923c',
        '701f98b02c8f691638ce4404e945d9c324190034b65ac5ada3389311fa9fb73f'),
    'upstream.patch': ('https://github.com/cmbant/CAMB/commit/56a95f78fdd72709c5af6b668141ec0c617777c3.patch',
        '25c3b93cc241ac0ccb611644a74f61f57cea0507526d448938c1df8efd310c72'),
}
OLD = 'call spline_def(State%Transfer_Times, scaling, State%num_transfer_redshifts,ddScaling)'
NEW = 'call spline_def(State%Transfer_Times(1:State%num_transfer_redshifts), scaling, State%num_transfer_redshifts,ddScaling)'


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def rel(path):
    path = Path(path).resolve()
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


def fetch(work):
    work.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request('https://pypi.org/pypi/camb/1.6.6/json', headers={'User-Agent':'cosmology-revalidation'})
    path = work/'pypi-1.6.6.json'
    if not path.exists(): path.write_bytes(urllib.request.urlopen(request, timeout=60).read())
    registry = json.loads(path.read_text())
    records = []
    for name,(url,expected) in ASSETS.items():
        if url == 'pypi':
            row = next(x for x in registry['urls'] if x['filename'] == name)
            assert row['digests']['sha256'] == expected
            url = row['url']
        path = work/name
        if not path.exists():
            data = urllib.request.urlopen(url, timeout=60).read()
            assert hashlib.sha256(data).hexdigest() == expected
            path.write_bytes(data)
        assert sha(path) == expected
        records.append({'path':rel(path), 'sha256':expected, 'url':url, 'bytes':path.stat().st_size})
    return records


def members(path):
    with tarfile.open(path) as stream:
        return {x.name.split('/',1)[1]:stream.extractfile(x).read() for x in stream.getmembers()
                if x.isfile() and '/' in x.name}


def installed_provenance(work, reference):
    site = next((reference/'lib').glob('python*/site-packages'))
    wheel = next(work.glob('*.whl')); archive = zipfile.ZipFile(wheel)
    files = {}
    for name in archive.namelist():
        if name.endswith('/') or name.endswith('/RECORD'): continue
        path = site/name
        assert path.is_file() and path.read_bytes() == archive.read(name), name
        files[rel(path)] = sha(path)
    dist = site/'camb-1.6.6.dist-info'
    for row in csv.reader(io.StringIO((dist/'RECORD').read_text())):
        if row[1]:
            algorithm, expected = row[1].split('=',1)
            assert algorithm == 'sha256'
            assert base64.urlsafe_b64encode(bytes.fromhex(sha(site/row[0]))).decode().rstrip('=') == expected
    return {'wheel_path':rel(wheel), 'wheel_sha256':sha(wheel), 'matching_wheel_payload_files':len(files),
        'installed_file_sha256':files, 'RECORD_sha256':sha(dist/'RECORD'),
        'direct_url_present':(dist/'direct_url.json').exists(), 'installer':(dist/'INSTALLER').read_text().strip(),
        'binary_comment':subprocess.check_output(['readelf','-p','.comment',str(site/'camb/camblib.so')],text=True),
        'interpretation':'Exact public wheel bytes and matching release sources; the wheel does not embed a complete reproducible build attestation.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work',type=Path,default=WORK)
    parser.add_argument('--reference',type=Path,default=REFERENCE)
    parser.add_argument('--compiler',type=Path,default=ROOT/'.work/unified-cosmology/calibrated-host-physics/toolchain/usr/bin/gfortran')
    args = parser.parse_args(); work = args.work.resolve(); reference=args.reference.resolve()
    assert work.is_relative_to(ROOT/'.work') and work != reference
    assert not (work/'build-record.json').exists(), 'A completed build record is immutable; use a fresh work directory.'
    assets=fetch(work); original=installed_provenance(work,reference)
    source=members(work/'camb-1.6.6.tar.gz'); tag=members(work/'camb-1.6.6-tag.tar.gz'); utils=members(work/'forutils-841f06.tar.gz')
    common=[p for p in source if p in tag and p.startswith(('camb/','fortran/'))]
    assert common and all(source[p]==tag[p] for p in common)
    utility_common=[p for p in utils if 'forutils/'+p in source]
    assert utility_common and all(utils[p]==source['forutils/'+p] for p in utility_common)
    compiler=args.compiler.resolve(); assert compiler.is_file()
    env=os.environ.copy(); env.update({k:'1' for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']})
    env['PATH']=str(args.compiler.parent.resolve())+os.pathsep+env['PATH'];env['MAKEFLAGS']='-j1'
    toolchain={'path':rel(compiler),'sha256':sha(compiler),'version':subprocess.check_output([str(compiler),'--version'],text=True),
               'make':subprocess.check_output(['make','--version'],text=True).splitlines()[0]}
    for name in ['f951','libgfortran.so','libquadmath.so','libgomp.so','libgcc_s.so']:
        path=Path(subprocess.check_output([str(compiler),'-print-prog-name='+name if name=='f951' else '-print-file-name='+name],text=True).strip())
        if path.is_file():toolchain[name]={'path':rel(path),'sha256':sha(path)}
    reference_site=next((reference/'lib').glob('python*/site-packages'))
    builds={}
    for arm in ['control','fixed']:
        tree=work/arm; assert not tree.exists(), 'Partial build preserved; choose a fresh work path.'
        tree.mkdir()
        for name,data in source.items():
            path=tree/name; assert path.resolve().is_relative_to(tree)
            path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(data)
        path=tree/'fortran/cmbmain.f90';text=path.read_text();assert text.count(OLD)==1
        if arm=='fixed':path.write_text(text.replace(OLD,NEW))
        source_hash={p:sha(tree/p) for p in source}
        log=work/(arm+'-build.log')
        command=['make','-j1','python','PYCAMB_OUTPUT_DIR='+str(tree/'camb')+'/', 'CLUSTER_SAFE=1']
        with log.open('x') as out:subprocess.run(command,cwd=tree/'fortran',env=env,stdout=out,stderr=subprocess.STDOUT,check=True,timeout=900)
        assert (tree/'camb/camblib.so').is_file()
        venv=work/('.'+arm+'-venv')
        subprocess.run(['uv','venv','--python',str(reference/'bin/python'),str(venv)],check=True,env=env)
        site=next((venv/'lib').glob('python*/site-packages'))
        (site/'00-local-camb.pth').write_text(str(tree)+'\n')
        (site/'99-reference-dependencies.pth').write_text(str(reference_site)+'\n')
        probe='import camb,json,importlib.metadata as m;print(json.dumps({"camb_file":camb.__file__,"version":camb.__version__,"versions":{d.metadata["Name"]:d.version for d in m.distributions() if d.metadata["Name"]}}))'
        info=json.loads(subprocess.check_output([str(venv/'bin/python'),'-c',probe],env=env,text=True))
        assert Path(info['camb_file']).resolve()==tree/'camb/__init__.py' and info['version']=='1.6.6'
        builds[arm]={'source_sha256':source_hash,'binary_path':rel(tree/'camb/camblib.so'),'binary_sha256':sha(tree/'camb/camblib.so'),
            'environment_path':rel(venv),'dependency_import_scope':'Isolated CAMB overlay; remaining imports read-only from pinned original environment.',
            'environment_pth_sha256':{rel(p):sha(p) for p in site.glob('*.pth')},'module_probe':info,
            'command':command,'log_path':rel(log),'log_sha256':sha(log),
            'ldd':subprocess.check_output(['ldd',str(tree/'camb/camblib.so')],text=True)}
    differences=[p for p in source if builds['control']['source_sha256'][p]!=builds['fixed']['source_sha256'][p]]
    assert differences==['fortran/cmbmain.f90'] or set(differences)=={'fortran/cmbmain.f90'}
    assert installed_provenance(work,reference)==original
    record={'schema':'isolated-camb-spline-build-v1','status':'built_matched_controls_no_physical_calls','physical_calls':0,
        'source_sha256':{rel(__file__):sha(__file__)},'assets':assets,'installed_provenance':original,
        'sdist_tag_identical_source_files':len(common),'sdist_submodule_identical_files':len(utility_common),
        'compiler':toolchain,'builds':builds,'changed_source_files':differences,
        'patch':{'old':OLD,'new':NEW,'scope':'Only explicit spline abscissa slice; adapted upstream56a95f fix to1.6.6 spline_def name.'}}
    (work/'build-record.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({'status':record['status'],'build_record':rel(work/'build-record.json')}))


if __name__=='__main__':main()
