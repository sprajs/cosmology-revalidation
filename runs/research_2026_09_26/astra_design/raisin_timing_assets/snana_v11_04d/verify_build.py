from pathlib import Path
import json,hashlib,difflib,os,subprocess
P=Path(__file__).parent;R=Path.cwd();source=P/'source';build=P/'build';a=json.loads((P/'acquisition.json').read_text());changes=[];patch=[]
for r in a['files']:
 p=source/r['path'];assert hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'];b=build/r['path']
 if p.read_bytes()!=b.read_bytes():
  changes.append(r['path']);patch+=list(difflib.unified_diff(p.read_text().splitlines(True),b.read_text().splitlines(True),fromfile='source/'+r['path'],tofile='build/'+r['path']))
(P/'build-generated-source.patch').write_text(''.join(patch));assert set(changes)<= {'src/sntools_output.h','src/genmag_PySEDMODEL.h'},changes
binary=build/'bin/snlc_fit.exe';assert binary.is_file();env=os.environ.copy();gsl=R/'phase2/official/build/sysroot/usr';env['LD_LIBRARY_PATH']=str(gsl/'lib')+(':'+env['LD_LIBRARY_PATH'] if 'LD_LIBRARY_PATH'in env else '');env['PATH']=str(gsl/'bin')+':'+env['PATH']
ldd=subprocess.check_output(['ldd',str(binary)],text=True,env=env);assert 'not found' not in ldd
res={'pass':True,'untouched_source_files_verified':len(a['files']),'build_modified_source_files':changes,'scientific_source_changes':False,'binary':str(binary.resolve()),'binary_bytes':binary.stat().st_size,'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),'ldd':ldd,'gfortran':subprocess.check_output(['gfortran','--version'],text=True,env=env).splitlines()[0],'gcc':subprocess.check_output(['gcc','--version'],text=True,env=env).splitlines()[0],'build_log_sha256':hashlib.sha256((P/'build.log').read_bytes()).hexdigest(),'scope':'Build/link gate only. No native fits and no original2021 binary/environment claim.'};(P/'build-verification.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps({k:v for k,v in res.items() if k!='ldd'},indent=2));print('PATCH:',''.join(patch))
