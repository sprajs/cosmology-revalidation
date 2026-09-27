"""Complete upstream setup.py's data-copy step after isolated make builds."""
import argparse
import json
from pathlib import Path
import shutil
import zipfile

import camb_spline_build as build


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work',type=Path,default=build.WORK)
    args=parser.parse_args();work=args.work.resolve()
    output=work/'package-record.json';assert not output.exists()
    record_path=work/'build-record.json';record=json.loads(record_path.read_text())
    for path,expected in record['source_sha256'].items():assert build.sha(build.ROOT/path)==expected
    archive=zipfile.ZipFile(next(work.glob('*.whl')))
    data={name:archive.read(name) for name in archive.namelist()
          if name.startswith('camb/') and name.endswith(('.py','.dat'))}
    rows={}
    for arm,row in record['builds'].items():
        tree=(build.ROOT/row['binary_path']).parent.parent
        template='HighLExtrapTemplate_lenspotentialCls.dat'
        source=tree/'fortran'/template;destination=tree/'camb'/template
        assert source.read_bytes()==data['camb/'+template]
        assert not destination.exists()
        shutil.copyfile(source,destination)
        for name,expected in data.items():assert (tree/name).read_bytes()==expected,name
        rows[arm]={'python_and_data_files_equal_public_wheel':len(data),
                   'file_sha256':{build.rel(tree/name):build.sha(tree/name) for name in data},
                   'binary_sha256':build.sha(build.ROOT/row['binary_path'])}
        assert rows[arm]['binary_sha256']==row['binary_sha256']
    value={'schema':'isolated-camb-spline-package-v1','status':'matched_python_data_and_different_declared_binaries',
           'physical_calls':0,'source_sha256':{build.rel(__file__):build.sha(__file__)},
           'original_build_record':build.rel(record_path),'original_build_sha256':build.sha(record_path),
           'scope':'Replays setup.py make_library template-copy step, omitted by direct make; original build receipt is preserved.',
           'packages':rows}
    output.write_text(json.dumps(value,indent=2)+'\n')
    print(json.dumps({'status':value['status'],'report':build.rel(output)}))


if __name__=='__main__':main()
