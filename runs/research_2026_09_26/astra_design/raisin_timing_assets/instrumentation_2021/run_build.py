from pathlib import Path
import subprocess,time,json,hashlib
P=Path(__file__).parent;start=time.monotonic()
with (P/'build.log').open('w') as f:r=subprocess.run(['timeout','--signal=TERM','--kill-after=5','180s','bash',str(P/'build.sh')],stdout=f,stderr=subprocess.STDOUT)
result={'exit_code':r.returncode,'elapsed_seconds':time.monotonic()-start,'command':['timeout','--signal=TERM','--kill-after=5','180s','bash',str(P/'build.sh')],'protocol_sha256':hashlib.sha256((P/'build-protocol.json').read_bytes()).hexdigest(),'build_script_sha256':hashlib.sha256((P/'build.sh').read_bytes()).hexdigest()};(P/'build-result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
