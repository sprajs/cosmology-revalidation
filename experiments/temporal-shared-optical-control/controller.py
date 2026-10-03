#!/usr/bin/env python3
"""One frozen SDK experiment: bounded transport, original169 comparisons, immutable receipts."""
import argparse
import datetime
from decimal import Decimal as D, getcontext
import io
import hashlib
import json
import os
from pathlib import Path
import re
import resource
import selectors
import signal
import stat
import subprocess
import sys
import tarfile
import time
import verify_sdk as V

getcontext().prec = 80
PACKET = Path(__file__).absolute().parent
ROOT = PACKET.parents[1]
RESULTS = ROOT / 'results' / 'temporal-shared-optical-control'
SLUG = re.compile(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}\Z')
GIT = '/usr/bin/git'
CONSUMED_AMENDMENT_SHA = '919ce4560762b2617f9c17ae63e322ad5fdfa627b644e6d26f0b20fb9ea91b3e'
SCHEDULE = ('repro-initial-head','repro-initial-branch','repro-initial-status','repro-initial-tree','repro-source-archive',
            'sdk-initial-head','sdk-initial-status','compiler-version','discovery','runtime-before','sdk-precompile-head','sdk-precompile-status',
            'compile','loader-before','sdk-prenative-head','sdk-prenative-status','native','sdk-prereference-head','sdk-prereference-status','reference',
            'loader-after-finally','runtime-after-finally','sdk-final-head','sdk-final-status','repro-final-head','repro-final-branch','repro-final-status','repro-final-tree')


class ChildFailure(RuntimeError):
    pass


class StoreLimit(ValueError):
    pass


class Capture:
    """Retain the exact consumed buffer independently of mutable producer files."""
    def __init__(self, attempt, artifacts, buffers):
        self.attempt = attempt; self.artifacts = artifacts; self.buffers = buffers
        self.raw = {}; self.parsed = {}; self.failures = []
        self.quota = None; self.consumed=set()

    def retain(self, label, blob, live, parsed=False, complete=True, observed_bytes=None):
        key = label
        if key in self.buffers: raise ValueError('duplicate capture label')
        self.buffers[key] = blob
        identity = {'path': str(Path(live).relative_to(self.attempt)), 'bytes': len(blob), 'sha256': V.digest(blob), 'complete': complete,
                    'observed_bytes':len(blob) if observed_bytes is None else observed_bytes}
        self.raw[key] = identity
        if parsed:
            self.consumed.add(key)
            if self.quota is not None: self.quota(len(blob))
            self.artifacts[key] = V.write_new(self.attempt / (label + '.consumed.json'), blob)
            try:
                value = V.parse(blob)
                self.parsed[key] = value
                return value
            except Exception as exc:
                self.failures.append({'stage': label + '-parse', 'kind': type(exc).__name__, 'message': str(exc)[:2048]})
        return None

    def terminal(self):
        """Every captured file is checked, even after failure; original objects stay intact."""
        for key, identity in list(self.raw.items()):
            for role, item in [('producer', identity), ('consumed', self.artifacts.get(key) if key in self.consumed else None)]:
                if item is None: continue
                try:
                    blob = V.read(self.attempt / item['path'])
                    if len(blob) != item['bytes'] or V.digest(blob) != item['sha256']:
                        raise ValueError('raw bytes differ from consumed buffer')
                except Exception as exc:
                    self.failures.append({'stage': 'raw-output-drift', 'label': key, 'role': role, 'path': item['path'], 'expected':dict(item), 'kind': type(exc).__name__, 'message': str(exc)[:2048]})
                    if role == 'consumed':
                        try:
                            recovered = key + '.recovered.consumed.json'
                            if self.quota is not None: self.quota(len(self.buffers[key]))
                            self.artifacts[key + '-recovered'] = V.write_new(self.attempt / recovered, self.buffers[key])
                            self.artifacts[key] = self.artifacts[key + '-recovered']
                        except Exception as recovery:
                            self.failures.append({'stage': 'consumed-buffer-recovery', 'label': key, 'kind': type(recovery).__name__, 'message': str(recovery)[:2048]})
                    elif key not in self.consumed:
                        try:
                            if self.quota is not None: self.quota(len(self.buffers[key]))
                            recovery=V.write_new(self.attempt/(key+'.recovered.log'),self.buffers[key])
                            self.artifacts[key+'-recovered']=recovery
                            self.failures[-1]['recovery_artifact']=recovery
                        except Exception as recovery:
                            self.failures.append({'stage':'consumed-buffer-recovery','label':key,'kind':type(recovery).__name__,'message':str(recovery)[:2048]})
        return not self.failures


class Runner:
    def __init__(self, attempt, request, artifacts, capture):
        self.attempt = attempt; self.request = request; self.artifacts = artifacts; self.capture = capture
        self.children = []; self.used = set(); self.skipped = {}; self.stop_cause = 'required preceding evidence was not earned'
        self.env = {'PATH': '/usr/bin:/bin', 'LANG': 'C', 'LC_ALL': 'C', 'TZ': 'UTC',
                    'TMPDIR': str(attempt / 'tmp'), **request['resources']['environment'],
                    **request['resources']['git_fixed_environment']}
        self.tool_baseline = V.file_identity(GIT, 67108864)[1]
        self.capture.quota=lambda extra: self.room(extra)

    def store_bytes(self):
        total = 0
        for path in self.attempt.rglob('*'):
            info = path.lstat()
            if stat.S_ISLNK(info.st_mode): raise ValueError('attempt contains symlink')
            if stat.S_ISREG(info.st_mode): total += info.st_size
        return total

    def total_room(self, extra=0):
        if self.store_bytes() + extra > self.request['resources']['attempt_store_bytes']:
            raise StoreLimit('attempt store quota')

    def room(self, extra=0):
        # Every ordinary writer shares the same intrinsic terminal reservation.
        self.total_room(extra+65536)

    def skip(self,label,cause):
        if label not in self.request['resources']['schedule'] or label in self.used:
            raise ValueError('cannot skip undeclared/already attempted label')
        self.skipped[label]=cause

    def consumption_map(self,label,args,git):
        """Cheap consumed-file maps, independently retaining errors; no311 scan."""
        files={}; errors=[]; wanted={}
        if label in ('compiler-version','discovery','compile','native','reference','loader-before','loader-after-finally'):
            sdk=self.request['sdk']; root=Path(sdk['artifacts_root'])
            for name in ('manifest','archive','cli','compiler','standard_library'):
                item=sdk[name]; path=Path(item['path'])
                if not path.is_absolute(): path=root/path
                if name in ('compiler','standard_library'): path=path.resolve()
                wanted['sdk/'+name]=(path,268435456,item['sha256'])
            item=self.request['native_dependencies_proposed']
            wanted['sdk/loader']=(Path(item['loader_path']).resolve(),67108864,item['loader_sha256'])
        wanted['tool/executable']=(Path(args[0]).resolve(),67108864,self.tool_baseline['sha256'] if git else None)
        if label in ('runtime-before','reference','runtime-after-finally'):
            item=self.request['reference_runtime']
            wanted['reference/python']=(Path(item['requested_executable']).resolve(),67108864,item['executable_sha256'])
            script=Path(args[1])
            for name,key in [('reference.py','reference'),('controller.py','controller'),('verify_sdk.py','sdk_verifier'),('consumer.cpp','consumer')]:
                wanted['source/'+name]=(script.parent/name,V.MIB,self.request['source_ports'][key]['sha256'])
        if label=='reference':
            wanted['input/request']=(self.attempt/'request.json',V.MIB,None)
            wanted['input/contract']=(self.attempt/'contract.snapshot.json',V.MIB,V.CONTRACT_SHA)
        if label=='compile':
            for arg in args:
                if arg.endswith('.cpp'): wanted['source/consumer']=(Path(arg),V.MIB,V.CONSUMER_SHA)
            for item in self.request['sdk']['headers']:
                wanted['header/'+item['path']]=(Path(self.request['sdk']['source_root'])/item['path'],V.MIB,item['sha256'])
        if label in ('loader-before','loader-after-finally'):
            wanted['compiled/consumer']=(Path(args[-1]),67108864,None)
        for name,(path,limit,pin) in wanted.items():
            try: files[name]=V.file_identity(path,limit,pin)[1]
            except Exception as exc:
                files[name]=None; errors.append({'path':str(path),'role':name,'kind':type(exc).__name__,'message':str(exc)[:2048]})
        return {'files':files,'errors':errors}

    def launch(self, label, args, git=False, allow=(0,), json_output=False):
        if label not in self.request['resources']['schedule'] or label in self.used:
            raise ValueError('undeclared/repeated direct launch: ' + label)
        if len(self.used) >= 28: raise ValueError('direct child count bound')
        if not args or any(type(arg) is not str for arg in args): raise ValueError('literal argument list required')
        self.used.add(label)
        if git: args = [GIT if arg == 'git' and i == 0 else arg for i,arg in enumerate(args)]
        r = self.request['resources']; cpu = 30 if git else 180; wall = 30 if git else 180
        file_limit = 4194304 if label == 'repro-source-archive' else (67108864 if label == 'compile' else V.MIB)
        def limits():
            resource.setrlimit(resource.RLIMIT_AS,(1073741824,1073741824))
            resource.setrlimit(resource.RLIMIT_CPU,(cpu,cpu))
            resource.setrlimit(resource.RLIMIT_FSIZE,(file_limit,file_limit))
            os.umask(0o077)
        receipt={'label':label,'argv':args,'cwd':str(self.attempt),'environment':self.env,'limits':{'cpu_seconds_per_process':cpu,'wall_seconds_process_group':wall,'address_space_bytes':1073741824,'regular_file_bytes':file_limit,'stdout_bytes':V.MIB,'stderr_bytes':V.MIB},
                 'launched':False,'status':'launch-failed','returncode':None,'wall_seconds':None,'wait4':None,'streams':{},'artifacts_before':None,'artifacts_after':None}
        self.children.append(receipt)
        receipt['started_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
        start=time.monotonic(); streams={'stdout':bytearray(),'stderr':bytearray()}; observed={'stdout':0,'stderr':0}; overflow=set()
        cause=None; proc=None; selector=None; usage=None; exit_status=None; streams_complete=False
        def read_ready(wait):
            for key,_ in selector.select(wait):
                part=os.read(key.fileobj.fileno(),65536)
                if not part: selector.unregister(key.fileobj); key.fileobj.close(); continue
                name=key.data; observed[name]+=len(part); remaining=V.MIB-len(streams[name]); streams[name].extend(part[:remaining])
                if len(part)>remaining: overflow.add(name)
        def kill_group():
            if proc is not None:
                try: os.killpg(proc.pid,signal.SIGKILL)
                except ProcessLookupError: pass
        try:
            # Refuse before a launch whose bounded known outputs would consume
            # the terminal reserve; compiler scratch still has the store watchdog.
            child_output=67108864 if label=='compile' else 6291456 if label=='repro-source-archive' else V.MIB if label=='reference' else 0
            self.room(child_output+4*V.MIB+(V.MIB if json_output else 0)+(V.MIB if label=='reference' else 0))
            before=self.consumption_map(label,args,git); receipt['artifacts_before']=before
            self.artifact(label+'-artifacts-before',before)
            if before['errors']: raise V.IdentityError('before-consumption identity refused',before)
            proc=subprocess.Popen(args,cwd=self.attempt,env=self.env,stdin=subprocess.DEVNULL,stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                                  preexec_fn=limits,start_new_session=True)
            receipt['launched']=True
            selector=selectors.DefaultSelector()
            for name, pipe in [('stdout',proc.stdout),('stderr',proc.stderr)]:
                os.set_blocking(pipe.fileno(),False); selector.register(pipe,selectors.EVENT_READ,name)
            while selector.get_map() or exit_status is None:
                if time.monotonic()-start>wall and cause is None: cause='timeout'
                try: self.room()
                except Exception: cause=cause or 'attempt-store-limit'
                if label=='reference' and (self.attempt/'reference.json').exists():
                    if (self.attempt/'reference.json').lstat().st_size>V.MIB: cause=cause or 'independent-output-limit'
                if cause is not None: kill_group()
                read_ready(0.05)
                if overflow: cause=cause or 'output-limit'
                if exit_status is None:
                    pid,code,rusage=os.wait4(proc.pid,os.WNOHANG)
                    if pid:
                        exit_status=code; usage=rusage; proc.returncode=os.waitstatus_to_exitcode(code)
            streams_complete=True
            receipt['returncode']=proc.returncode
            receipt['wait4']={'cpu_user_seconds':usage.ru_utime,'cpu_system_seconds':usage.ru_stime,'maximum_rss_kib':usage.ru_maxrss,
                              'scope':'OS wait4 attribution for reaped direct child; not an enforced aggregate descendant CPU budget'}
            receipt['status']=cause or ('completed' if proc.returncode in allow else 'nonzero-exit')
        except BaseException as exc:
            receipt['exception']={'kind':type(exc).__name__,'message':str(exc)[:2048]}
            receipt['status']='interrupted' if proc is not None else 'launch-failed'
            # A reaped parent can still have descendants holding our pipes.
            kill_group()
            if proc is not None:
                if exit_status is None:
                    try:
                        _,code,usage=os.wait4(proc.pid,0); exit_status=code; proc.returncode=os.waitstatus_to_exitcode(code)
                    except ChildProcessError: pass
                receipt['returncode']=proc.returncode
                if selector is not None:
                    deadline=time.monotonic()+1
                    try:
                        while selector.get_map() and time.monotonic()<deadline: read_ready(0.05)
                        streams_complete=not selector.get_map()
                    except BaseException as drain:
                        receipt['drain_failure']={'kind':type(drain).__name__,'message':str(drain)[:2048]}
        finally:
            if receipt['status']!='completed': kill_group()
            if selector is not None: selector.close()
            if proc is not None:
                for pipe in (proc.stdout,proc.stderr):
                    if pipe is not None and not pipe.closed: pipe.close()
            receipt['wall_seconds']=time.monotonic()-start
            receipt['finished_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
            receipt['process_status']=receipt['status']
            # Pair the same consumed identities immediately, on every outcome.
            try:
                after=self.consumption_map(label,args,git); receipt['artifacts_after']=after
                self.artifact(label+'-artifacts-after',after)
                if after['errors']: raise V.IdentityError('after-consumption identity refused',after)
                if receipt['artifacts_before'] is not None: V.exact(after,receipt['artifacts_before'],'consumed identities before/after')
            except Exception as exc:
                receipt['identity_failure']={'kind':type(exc).__name__,'message':str(exc)[:2048],'details':getattr(exc,'details',None)}
                if receipt['status']=='completed': receipt['status']='identity-failed'
            for name in ('stdout','stderr'):
                if not receipt['launched']: continue  # Unearned streams remain null.
                blob=bytes(streams[name]); path=self.attempt/(label+'.'+name+'.log')
                identity={'path':path.name,'bytes':len(blob),'sha256':V.digest(blob),'complete':streams_complete and name not in overflow,'observed_bytes':observed[name]}
                receipt['streams'][name]=identity
                try:
                    self.room(len(blob)); artifact=V.write_new(path,blob)
                    artifact.update(complete=identity['complete'],observed_bytes=observed[name]); self.artifacts[label+'-'+name]=artifact
                except Exception as exc:
                    receipt['capture_failure']={'kind':type(exc).__name__,'message':str(exc)[:2048]}; receipt['status']='capture-failed'
                try: self.capture.retain(label+'-'+name,blob,path,parsed=json_output and name=='stdout',complete=identity['complete'],observed_bytes=observed[name])
                except Exception as exc:
                    receipt['capture_failure']={'kind':type(exc).__name__,'message':str(exc)[:2048]}; receipt['status']='capture-failed'
            if usage is not None and receipt['wait4'] is None:
                receipt['wait4']={'cpu_user_seconds':usage.ru_utime,'cpu_system_seconds':usage.ru_stime,'maximum_rss_kib':usage.ru_maxrss,'scope':'OS wait4 attribution'}
        if receipt['status']!='completed': raise ChildFailure(label+': '+receipt['status'])
        return bytes(streams['stdout'])

    def artifact(self, label, value, limit=V.MIB):
        blob=V.encode(value); self.room(len(blob)); identity=V.write_new(self.attempt/(label+'.json'),blob,limit)
        self.artifacts[label]=identity; return identity


def source_map(root, tree):
    """Check every committed source leaf against Git's recorded object identity."""
    entries={}; total=0; errors=[]
    for raw in tree.split(b'\0'):
        if not raw: continue
        metadata,name=raw.split(b'\t',1); mode,kind,oid=metadata.decode('ascii').split(' ')
        name=name.decode('utf-8','strict'); path=Path(name)
        if path.is_absolute() or '..' in path.parts or mode not in ('100644','100755') or kind!='blob': raise ValueError('unsupported committed source entry')
        try:
            if total>2097152: raise ValueError('remaining source leaves unavailable after2MiB bound')
            V.walk_directory((root/path).parent); blob=V.read(root/path)
            total+=len(blob)
            if total>2097152: raise ValueError('committed source tree exceeds2MiB')
            object_bytes=b'blob '+str(len(blob)).encode()+b'\0'+blob
            actual_oid=hashlib.sha1(object_bytes).hexdigest()
            entries[name]={'mode':mode,'git_blob':oid,'actual_git_blob':actual_oid,'bytes':len(blob),'sha256':V.digest(blob)}
            if actual_oid!=oid: errors.append({'path':name,'kind':'ValueError','message':'working source differs from committed blob'})
        except Exception as exc:
            entries[name]=None; errors.append({'path':name,'kind':type(exc).__name__,'message':str(exc)[:2048]})
    if not entries: raise ValueError('empty committed source')
    if errors: raise V.IdentityError('committed source map refused',{'sources':entries,'errors':errors})
    return entries


def repro_map(root, phase, run, initial_head=None):
    base=['git',*V.GIT_OPTIONS,'-C',str(root)]
    errors=[]
    def collect(name,args,allow=(0,)):
        try: return run('repro-'+phase+'-'+name,base+args,git=True,allow=allow)
        except BaseException as exc:
            errors.append({'stage':name,'kind':type(exc).__name__,'message':str(exc)[:2048]}); return None
    head_raw=collect('head',['rev-parse','--verify','HEAD'])
    branch_raw=collect('branch',['symbolic-ref','--quiet','--short','HEAD'],(0,1))
    status=collect('status',['status','--porcelain=v1','-z','--untracked-files=all'])
    def decode(raw,name):
        if raw is None: return None
        try:
            text=raw.decode('utf-8','strict').strip()
            if name=='head' and re.fullmatch('[0-9a-f]{40}',text) is None: raise ValueError('committed HEAD token differs')
            return text
        except Exception as exc:
            errors.append({'stage':name+'-parse','kind':type(exc).__name__,'message':str(exc)[:2048]}); return None
    head=decode(head_raw,'head'); branch=decode(branch_raw,'branch')
    tree_revision=head if head is not None else initial_head
    tree=None; sources=None
    if tree_revision is not None:
        tree=collect('tree',['ls-tree','-r','-z','--full-tree',tree_revision])
    else:
        errors.append({'stage':'tree','kind':'NotStarted','message':'current and initial committed HEAD unavailable'})
        owner=getattr(run,'__self__',None)
        if isinstance(owner,Runner): owner.skip('repro-'+phase+'-tree','current and initial committed HEAD unavailable after HEAD metadata refusal')
    if tree is not None:
        try: sources=source_map(root,tree)
        except Exception as exc:
            details=getattr(exc,'details',None) or {}; sources=details.get('sources')
            errors.append({'stage':'committed-source','kind':type(exc).__name__,'message':str(exc)[:2048],'details':details.get('errors')})
    if status is not None and status!=b'': errors.append({'stage':'clean-source','kind':'ValueError','message':'Repro tracked/untracked status is dirty'})
    if initial_head is not None and head is not None and head!=initial_head: errors.append({'stage':'HEAD','kind':'ValueError','message':'Repro committed HEAD changed'})
    return {'head':head,'branch':branch,'tracked_and_untracked_status':status.decode('utf-8','backslashreplace') if status is not None else None,
            'tree_revision':tree_revision,'tree_sha256':V.digest(tree) if tree is not None else None,'sources':sources,'errors':errors}


def admit_source(value):
    if value['errors']: raise V.IdentityError('committed source admission refused',{'errors':value['errors']})
    return value


def unpack_snapshot(attempt, source):
    archive=attempt/'source.tar'; blob=V.read(archive,4194304); total=0; files={}
    with tarfile.open(fileobj=io.BytesIO(blob),mode='r:') as tar:
        members=tar.getmembers()
        if len(members)>4096: raise ValueError('archive member limit')
        for entry in members:
            relative=Path(entry.name)
            if relative.is_absolute() or '..' in relative.parts or entry.issym() or entry.islnk(): raise ValueError('escaped/link archive member')
            target=attempt/'source'/relative
            if entry.isdir(): target.mkdir(parents=True,exist_ok=True); V.walk_directory(target); continue
            if not entry.isfile() or entry.size>V.MIB or entry.name not in source['sources']: raise ValueError('unsupported/extra archive member')
            total+=entry.size
            if total>2097152: raise ValueError('archive source tree bound')
            target.parent.mkdir(parents=True,exist_ok=True); content=tar.extractfile(entry).read(V.MIB+1)
            wanted=source['sources'][entry.name]
            V.exact(len(content),wanted['bytes']); V.exact(V.digest(content),wanted['sha256'])
            files[entry.name]=V.write_new(target,content,mode=0o555 if wanted['mode']=='100755' else 0o444)
    V.exact(set(files),set(source['sources']),'archive complete inventory')
    os.chmod(archive,0o444)
    return {'path':'source.tar','bytes':len(blob),'sha256':V.digest(blob)}


def loader_map(blob, request):
    result=[]
    for line in blob.decode('utf-8','strict').splitlines():
        text=line.strip()
        if not text: continue
        if text.startswith('linux-vdso.so.'): continue  # virtual mapping: explicitly outside resolved-file scope
        if 'not found' in text: raise ValueError('loader dependency unavailable')
        if ' => ' in text: path=text.split(' => ',1)[1].split(' (',1)[0]
        else: path=text.split(' (',1)[0]
        p=Path(path)
        if not p.is_absolute(): raise ValueError('unrecognized loader resolution line')
        _,identity=V.file_identity(p.resolve(),67108864); result.append(identity)
    unique={item['resolved_path']:item for item in result}
    if len(unique)!=len(result) or not unique or len(unique)>256: raise ValueError('loader dependency inventory')
    hashes={item['sha256'] for item in unique.values()}
    if request['sdk']['standard_library']['sha256'] not in hashes: raise ValueError('compiled consumer resolved standard library differs from SDK')
    if request['native_dependencies_proposed']['loader_sha256'] not in hashes: raise ValueError('resolved ELF loader differs from pin')
    return {'scope':request['native_dependencies_proposed']['scope'],'libraries':[unique[key] for key in sorted(unique)]}


def source_port_hashes(packet, contract):
    result={name:V.file_identity(packet/name)[1]['sha256'] for name in ('consumer.cpp','reference.py','controller.py','verify_sdk.py')}
    result['contract.snapshot.json']=V.file_identity(contract)[1]['sha256']; return result


def controller_runtime_collect(module_items,maps_text,status_text):
    """Bound discovery and consumption; retain modules and lexical mappings apart."""
    paths=set(); modules=[]; mapped=set()
    def add(path):
        path=Path(path)
        if not path.is_absolute(): raise ValueError('runtime origin must be absolute')
        resolved=path.resolve(strict=True)
        if not stat.S_ISREG(resolved.stat().st_mode): raise ValueError('runtime origin is not regular')
        paths.add(resolved)
        if len(paths)>4096: raise ValueError('controller runtime file count')
        return str(resolved)
    for path in (Path(sys.executable),Path(__file__).absolute(),Path(V.__file__).absolute()): add(path)
    for name,module in module_items:
        if len(modules)>=4096: raise ValueError('controller module count')
        spec=getattr(module,'__spec__',None); origin=getattr(spec,'origin',None); cached=getattr(module,'__cached__',None)
        if origin is None: origin=getattr(module,'__file__',None)
        namespace=[]; locations=getattr(spec,'submodule_search_locations',None)
        if origin in ('built-in','frozen'): kind=origin; recorded=None
        elif origin is not None: kind='source-or-extension'; recorded=add(origin)
        elif locations is not None:
            kind='namespace'; recorded=None
            for location in locations:
                if len(namespace)>=32: raise ValueError('controller namespace count')
                path=Path(location).resolve(strict=True); V.walk_directory(path); namespace.append(str(path))
        else:
            # A module without a file must have an explicit builtin/frozen identity.
            raise ValueError('controller module origin unavailable: '+name)
        modules.append({'name':name,'kind':kind,'origin':recorded,'namespace_paths':sorted(set(namespace))})
        if cached is not None and Path(cached).exists(): add(cached)
    for line in maps_text.splitlines():
        parts=line.split(maxsplit=5)
        if len(parts)==6 and '.so' in parts[5]:
            lexical=re.sub(r'\\([0-7]{3})',lambda match:chr(int(match.group(1),8)),parts[5])
            if not lexical.startswith('/') or lexical.endswith(' (deleted)'): raise ValueError('mapped runtime library unavailable')
            resolved=add(lexical); mapped.add((lexical,resolved))
            if len(mapped)>256: raise ValueError('controller mapped library count')
    inventory=[]; total=0
    for path in sorted(paths):
        # Check actual stat before reading; the bounded read also checks drift.
        size=path.stat().st_size
        if size>67108864-total: raise ValueError('controller runtime byte quota before read')
        _,item=V.file_identity(path,67108864-total)
        total+=item['bytes']; inventory.append(item)
        if total>67108864: raise ValueError('controller runtime actual byte quota')
    threads=[line.split() for line in status_text.splitlines() if line.startswith('Threads:')]
    if threads!=[['Threads:','1']]: raise ValueError('controller actual thread count differs')
    if getcontext().prec!=80: raise ValueError('controller Decimal precision differs')
    return {'version':sys.version,'executable':str(Path(sys.executable).resolve()),'dont_write_bytecode':sys.dont_write_bytecode,'actual_thread_count':1,
            'decimal_precision':80,'modules':sorted(modules,key=lambda item:item['name']),
            'mapped_libraries':[{'path':path,'resolved_path':resolved} for path,resolved in sorted(mapped)],'inventory':inventory}


def controller_runtime():
    """Actual controller code/runtime map; stdlib only, no reference imports."""
    if len(sys.modules)>4096: raise ValueError('controller module count')
    return controller_runtime_collect(sorted(sys.modules.items()),Path('/proc/self/maps').read_text(),Path('/proc/self/status').read_text())
def compare(native, ref):
    checks = []
    def check(name, passed, **details):
        checks.append(dict(name=name, passed=bool(passed), **details))
    primary = native["primary"]
    check("one-immutable-owner", native["prepared_status"]=="ok" and native["grid_count"]==1 and native["band_count"]==4)
    check("all-primary-rows", primary["status"]=="ok" and primary["required_rows_admitted"] and len(primary["rows"])==8)
    max_photon = D(0)
    for i, row in enumerate(primary["rows"]):
        n = row["photons"]["value"]
        a = D(ref["photons"][i]); ae = D(ref["photon_errors"][i])
        budget = D("1e-300")+D("2e-12")*abs(a)
        used = abs(n-a)+ae
        max_photon = max(max_photon, used/budget)
        check("photons-"+str(i), row["admission_status"]=="ok" and row["photons"]["availability"]==1 and row["photons"]["status"]=="ok" and n>0 and used<=budget,
              native=n, reference=a, earned_reference_error=ae, total_used=used, allowance=budget)
        check("mask-"+str(i), row["mean_flux"]["availability"]==0 and row["energy"]["availability"]==0)
        e = (i%4)//2
        check("identity-coverage-"+str(i), row["grid_index"]==0 and row["band_index"]==2*(i//4)+(i%2) and row["source_epoch_second"]==10 and
              row["observer_duration_second"]==(2 if e==0 else 3) and row["covered_observer_second"]==2 and row["coverage"]==("full" if e==0 else "partial") and
              abs(row["covered_fraction"]-(D(1) if e==0 else D(2)/3))<D("1e-18"))
        f = ref["frequency_rows"][i]
        f32, f64 = D(f["frequency32"]), D(f["frequency64"])
        e32, e64 = D(f["frequency32_error"]), D(f["frequency64_error"])
        fb = D("1e-300")+D("2e-13")*abs(a)
        check("frequency-polynomial-"+str(i), abs(f64-a)+e64+ae<=fb,
              used=abs(f64-a)+e64+ae, allowance=fb)
        check("frequency-refinement-"+str(i), abs(f64-f32)+e64+e32<=fb,
              used=abs(f64-f32)+e64+e32, allowance=fb)
    max_log = D(0);max_ref = D(0);witness = D(0);correlation_receipts=[]
    for f in ref["records"]:
        s, r = f["sigma_index"], f["record_index"]
        actual = native["detector_controls"][s]["records"][r]
        check(f"record-status-{s}-{r}", actual["status"]=="ok")
        for quantity, error in (("log_joint","joint_log_error"),("log_event","event_log_error"),("log_value","value_log_error")):
            v = D(f[quantity]);earned = D(f[error]);n = actual[quantity]
            budget = D("2e-10")*(1+abs(v));used = abs(n-v)+earned
            max_log=max(max_log,used/budget);max_ref=max(max_ref,earned/budget)
            check(f"{quantity}-{s}-{r}", used<=budget and earned<=D("0.05")*budget,
                  native=n, reference=v, earned_reference_log_error=earned, total_used=used, allowance=budget)
            estimate=actual[error]
            check(f"empirical-{quantity}-{s}-{r}", estimate<=D("5e-12")+D("1e-8")*(1+abs(n)), estimate=estimate)
        if r==1:
            check(f"same-selection-denominator-{s}", abs(actual["log_value"]-(actual["log_joint"]-actual["log_event"]))<=D("1e-17"))
        if r==0:witness=max(witness,min(abs(actual["record_relative_difference"]),abs(D(f["record_relative_difference"]))-D(f["record_relative_difference_error"])))
        correlation_receipts.append({"sigma_index":s,"record_index":r,
                                     "native_record_relative_difference":actual["record_relative_difference"],
                                     "reference_record_relative_difference":f["record_relative_difference"],
                                     "record_earned_reference_relative_error":f["record_relative_difference_error"],
                                     "native_event_relative_difference":actual["event_relative_difference"],
                                     "reference_event_relative_difference":f["event_relative_difference"],
                                     "event_earned_reference_relative_error":f["event_relative_difference_error"],
                                     "reference_record_signed_difference":f["record_signed_difference"],
                                     "reference_record_signed_difference_error":f["record_signed_difference_error"],
                                     "reference_event_signed_difference":f["event_signed_difference"],
                                     "reference_event_signed_difference_error":f["event_signed_difference_error"]})
    check("meaningful-common-state-numerator-witness",witness>D("1e-3"), max_absolute_relative_difference=witness)
    for noise, control in enumerate(native["detector_controls"]):
        check(f"state-channel-count-{noise}",len(control["attempts"])==8)
        for i, attempt in enumerate(control["attempts"]):
            check(f"state-channel-map-{noise}-{i}",attempt["state_index"]==i//4 and attempt["channel_index"]==i%4)
            for a, batch in enumerate(attempt["batches"]):
                check(f"required-attempt-{noise}-{i}-{a}",batch["status"]=="ok" and len(batch["rows"])==6 and all(row["status"]=="ok" for row in batch["rows"]))
    for control in native["actual_temporal_controls"]:
        name, b = control["id"], control["native"]
        if name=="no-time-overlap":
            row=b["rows"][0]
            check(name, row["admission_status"]=="ok" and row["coverage"]=="no_overlap" and row["photons"]["value"]==0 and row["observer_duration_second"]==2)
        elif name in ("zero-duration","invalid-redshift"):
            check(name,len(b["rows"])==8 and b["rows"][0]["admission_status"]=="outside_domain" and not b["required_rows_admitted"] and control["joint_withheld"])
        elif name=="required-state-row-refusal":
            check(name,len(b["rows"])==8 and b["rows"][7]["admission_status"]!="ok" and not b["required_rows_admitted"] and control["joint_withheld"])
        elif name=="temporal-work-refusal":
            check(name,b["status"]=="ok" and len(b["rows"])==8 and all(row["admission_status"]=="work_limit" and row["photons"]["availability"]!=1 for row in b["rows"]) and control["joint_withheld"])
        else:
            check(name,b["status"]=="work_limit" and control["joint_withheld"])
    for control in native["actual_detector_controls"]:
        name, b = control["id"], control["native"]
        if name=="detector-QE-domain": check(name,b["status"]=="outside_domain")
        elif name=="detector-work-refusal": check(name,b["status"]=="ok" and b["rows"][0]["status"]=="conditioning_budget_exceeded" and b["rows"][0]["log_value"] is None)
        elif name=="selected-nondetection": check(name,b["rows"][0]["status"]=="invalid_input")
        else:check(name,b["status"]=="ok" and b["rows"][0]["status"]=="ok" and b["rows"][0]["zero_probability"] and b["rows"][0]["log_value"] is None)
    for control in native["actual_external_reducer_controls"]:
        r=control["reducer"]
        expected={"required-QE-refusal":"outside_domain","required-work-refusal":"conditioning_budget_exceeded","selected-censored-refusal":"invalid_input"}[control["id"]]
        check(control["id"],r["status"]==expected and all(r[q] is None for q in ("log_joint","log_event","log_value","product_per_exposure_log_record","product_per_exposure_log_event","record_relative_difference","event_relative_difference")) and control["retained_state_masses"]==[D("0.25"),D("0.75")])
    counted=sum(b["poisson_terms"] for control in native["detector_controls"] for attempt in control["attempts"] for b in attempt["batches"] if b is not None)
    check("primary-poisson-work-accounting",native["primary_poisson_terms"]==counted)
    counted+=sum(c["native"]["poisson_terms"] for c in native["actual_detector_controls"])
    counted+=sum(c["replacement_batch"]["poisson_terms"] for c in native["actual_external_reducer_controls"] if c["replacement_batch"] is not None)
    check("all-poisson-work-accounting",native["all_poisson_terms"]==counted)
    check("all-temporal-work-accounting",native["all_temporal_segment_work"]==primary["segment_work"]+sum(c["native"]["segment_work"] for c in native["actual_temporal_controls"]))
    return {"accepted": all(c["passed"] for c in checks), "checks": checks,
            "check_count":len(checks), "failed":[c for c in checks if not c["passed"]],
            "maximum_photon_budget_fraction":max_photon,"maximum_log_budget_fraction":max_log,
            "maximum_reference_log_budget_fraction":max_ref,"correlation_witness":witness,"correlation_receipts":correlation_receipts}



def terminal_record(attempt, record, artifacts, quota, total_quota=None):
    """Encode before exclusive create; a bounded fallback survives receipt overflow."""
    try:
        blob=V.encode(record)
        if len(blob)>V.MIB: raise ValueError('terminal receipt exceeds1MiB')
        quota(len(blob)); V.write_new(attempt/'record.json',blob)
        return record
    except Exception as exc:
        if (attempt/'record.json').exists(): raise
        selected=[]; omitted=0
        for key,item in sorted(artifacts.items()):
            path=item.get('path') if type(item) is dict else None
            valid=type(path) is str and len(path.encode())<=256 and not Path(path).is_absolute() and '..' not in Path(path).parts
            valid=valid and type(item.get('bytes')) is int and item['bytes']>=0 and type(item.get('sha256')) is str and V.SHA_PATTERN.fullmatch(item['sha256']) is not None
            if len(selected)<128 and valid:
                selected.append({'label':key[:128],'path':path,'bytes':item['bytes'],'sha256':item['sha256']})
            else: omitted+=1
        pinned=record.get('request_sha256')
        if pinned is None and type(record.get('request')) is dict: pinned=record['request'].get('sha256')
        fallback={'schema_version':1,'interface_id':'temporal-shared-optical-terminal-refusal/v1','status':'failed','numerical_accepted':False,
                  'request_sha256':pinned,'reason':'terminal-serialization-or-quota','message':str(exc).encode()[:2048].decode('utf-8','ignore'),
                  'artifacts':selected,'unavailable_artifact_count':omitted,'seal':'regular files0444; executable0555; exclusive terminal'}
        blob=V.encode(fallback)
        while len(blob)>65536 and selected:
            selected.pop(); omitted+=1; fallback['unavailable_artifact_count']=omitted; blob=V.encode(fallback)
        if len(blob)>65536: raise ValueError('bounded terminal fallback exceeded')
        # Reserve64KiB before work starts; fallback uses that reserved store space.
        (total_quota if total_quota is not None else quota)(len(blob))
        V.write_new(attempt/'record.json',blob,65536)
        return fallback


def apply_capture_gate(record,capture):
    """A previously passing numerical object remains evidence after raw drift."""
    capture.terminal(); record['failures'].extend(capture.failures)
    if record['failures']: record['numerical_accepted']=False
    record['status']='completed' if record['numerical_accepted'] else 'failed'


def seal(attempt,terminal=True):
    for path in attempt.rglob('*'):
        info=path.lstat()
        if stat.S_ISLNK(info.st_mode): raise ValueError('cannot seal symlink')
        if stat.S_ISREG(info.st_mode): os.chmod(path,0o555 if path.name=='consumer' else 0o444)
    if terminal:
        for path in sorted((p for p in attempt.rglob('*') if p.is_dir()),key=lambda p:len(p.parts),reverse=True): os.chmod(path,0o555)
        os.chmod(attempt,0o555)


def parsed_receipt(label,capture,admitted=None,error=None):
    value=capture.parsed.get(label); identity=capture.artifacts.get(label)
    if label not in capture.buffers: return None
    parse_errors=[failure for failure in capture.failures if failure.get('stage')==label+'-parse']
    schema=(value.get('interface_id') or value.get('schema')) if type(value) is dict else None
    return {'consumed_bytes_sha256':V.digest(capture.buffers[label]),'schema':schema,
            'admitted':bool(admitted) if admitted is not None else False,'admission_error':error or (parse_errors[0] if parse_errors else None),
            'value_artifact':identity}


def captured_identity(attempt,item,error=None):
    if item is None: return None
    complete=item.get('complete',True)
    return {'path':item['path'],'resolved_path':str(attempt/item['path']),'bytes':item['bytes'],'sha256':item['sha256'],
            'status':'captured' if complete and error is None else 'refused','error':error,'observed_bytes':item.get('observed_bytes',item['bytes']),'truncated':not complete}


def not_started_record(label,cause,environment=None):
    return {'name':label,'command':None,'environment':environment,'limits':None,'started_utc':None,'finished_utc':None,
            'launch_status':'not_started','exit_code':None,'signal':None,'termination_cause':'not_started','wall_seconds':None,
            'user_cpu_seconds':None,'system_cpu_seconds':None,'maximum_rss_kib':None,'stdout':None,'stderr':None,'independent_output':None,'parsed_output':None,
            'parse_error':{'kind':'NotStarted','message':cause}}


def child_records(runner,capture,native_admitted,reference_admitted):
    result=[]
    V.exact(runner.request['resources']['schedule'],list(SCHEDULE),'fixed28 schedule')
    causes={'timeout':'timeout','output-limit':'stdout_limit','independent-output-limit':'independent_output_limit','attempt-store-limit':'store_limit','launch-failed':'launch_failure','interrupted':'interrupted','capture-failed':'store_limit'}
    by_label={child['label']:child for child in runner.children}
    for label in SCHEDULE:
        if label not in by_label:
            result.append(not_started_record(label,runner.skipped.get(label,runner.stop_cause),runner.env))
            continue
        child=by_label[label]
        code=child['returncode']; usage=child['wait4']; label=child['label']
        signal_number=-code if code is not None and code<0 else None
        termination=child.get('process_status',child['status'])
        if termination=='completed' and child['status']=='capture-failed': termination='capture-failed'
        cause=causes.get(termination,'signal' if signal_number else 'exit')
        if termination=='capture-failed' and child.get('capture_failure',{}).get('kind')!='StoreLimit': cause='interrupted'
        if termination=='output-limit' and child['streams'].get('stderr',{}).get('complete') is False: cause='stderr_limit'
        admission=native_admitted if label=='native' else reference_admitted if label=='reference' else child['status']=='completed'
        result.append({'name':label,'command':child['argv'],'environment':child['environment'],'limits':child['limits'],
                       'started_utc':child.get('started_utc'),'finished_utc':child.get('finished_utc'),'launch_status':'started' if child['launched'] else 'launch_failed',
                       'exit_code':code,'signal':signal_number,'termination_cause':cause,'wall_seconds':child['wall_seconds'],
                       'user_cpu_seconds':usage['cpu_user_seconds'] if usage else None,'system_cpu_seconds':usage['cpu_system_seconds'] if usage else None,
                       'maximum_rss_kib':usage['maximum_rss_kib'] if usage else None,
                       'stdout':captured_identity(runner.attempt,child['streams'].get('stdout'),child.get('capture_failure')),
                       'stderr':captured_identity(runner.attempt,child['streams'].get('stderr'),child.get('capture_failure')),
                       'independent_output':captured_identity(runner.attempt,capture.raw.get('reference-output')) if label=='reference' else None,
                       'parsed_output':parsed_receipt(label+'-stdout',capture,admission,child.get('identity_failure') or child.get('capture_failure')),
                       'parse_error':next((e for e in capture.failures if e.get('stage')==label+'-stdout-parse'),child.get('identity_failure') or child.get('capture_failure') or child.get('exception'))})
    return result


def output_inventory(attempt):
    result={}; total=0
    for path in sorted(attempt.rglob('*')):
        info=path.lstat()
        if stat.S_ISLNK(info.st_mode): raise ValueError('output inventory symlink')
        if not stat.S_ISREG(info.st_mode): continue
        relative=str(path.relative_to(attempt))
        if relative=='record.json': continue
        if len(result)>=4096: raise ValueError('output inventory count')
        total+=info.st_size
        if total>268435456: raise ValueError('output inventory store bound')
        blob=V.read(path,67108864 if path.name=='consumer' else 4194304 if path.name=='source.tar' else V.MIB)
        result[relative]={'bytes':len(blob),'sha256':V.digest(blob),'mode':stat.S_IMODE(info.st_mode)}
    return result


def final_receipt(internal,slug,artifacts,capture,inventory,native_admitted,reference_admitted):
    errors=internal['failures']; gates={}
    def gate(name,earned,stages=()):
        relevant=[error for error in errors if any(error['stage'].startswith(stage) for stage in stages)]
        gates[name]={'status':'failed' if relevant else 'passed' if earned else 'not_earned',
                     'cause':relevant[0] if relevant else None if earned else 'required evidence was not earned'}
    gate('execution',internal.get('execution_complete',False),('execution',))
    gate('source_admission',internal.get('source_admitted',False),('source-admission',))
    gate('source_unchanged','repro-final' in artifacts,('source-terminal','input-terminal','artifact-terminal-source-archive'))
    gate('sdk_admission',internal.get('sdk_admitted',False),('sdk-admission',))
    gate('sdk_unchanged','sdk-final' in artifacts,('sdk-terminal','sdk-artifact-terminal'))
    gate('controller_runtime_unchanged','controller-runtime-after' in artifacts,('controller-runtime-terminal',))
    gate('reference_runtime_admission','runtime-before-admitted' in artifacts)
    gate('reference_runtime_unchanged','runtime-after-admitted' in artifacts,('runtime-terminal',))
    gate('compiled_consumer_unchanged','compiled-consumer' in artifacts,('compiled-consumer-terminal',))
    gate('native_dependencies_unchanged','loader-after-resolution' in artifacts,('native-terminal-identity',))
    gate('raw_output_unchanged',bool(capture.raw),('raw-output-drift','consumed-buffer-recovery','artifact-terminal'))
    gate('native_output_admission',native_admitted,('native-terminal-admission','native-output-admission'))
    gate('reference_output_admission',reference_admitted,('reference-terminal-admission','reference-output-admission'))
    gate('original169_comparisons',internal['comparison_accepted_before_terminal'],('original169-comparisons',))
    required=list(gates.values())
    numerical=internal['numerical_accepted'] and all(item['status']=='passed' for item in required)
    gates['numerical']={'status':'passed' if numerical else 'failed','cause':None if numerical else 'one or more execution/admission/identity/original169 gates unavailable or failed'}
    gates['inference']={'status':'blocked','cause':internal['inference_status']}
    gates['interpretation']={'status':'blocked','cause':'synthetic control does not qualify a physical/event-calibrated likelihood'}
    inventory_sha=V.digest(json.dumps(inventory,sort_keys=True,separators=(',',':'),allow_nan=False).encode()) if inventory is not None else None
    qualification={'scope':internal['qualification'],'physical_qualification':False,'observational_qualification':False,'original_input_certificate':False,
                   'consumed_buffer_amendment_sha256':CONSUMED_AMENDMENT_SHA,'direct_launches':internal.get('direct_launches',0),'maximum_direct_launches':28,
                   'cpu_scope':'RLIMIT_CPU per-process; actual wait4 OS attribution, compiler descendants under same serial group wall/store watchdog',
                   'native_dependencies_scope':'loader resolution only; no claim of native process self-maps','original169_comparison_semantics':'unchanged'}
    return {'schema_version':1,'interface_id':'temporal-shared-optical-attempt-record/v1','created_utc':internal['started_utc'],'finished_utc':internal['ended_utc'],'attempt':slug,
            'request':{'sha256':internal['request_sha256'],'value_artifact':artifacts.get('request'),'contract_artifact':artifacts.get('contract')},
            'source_before':artifacts.get('repro-initial'),'source_after':artifacts.get('repro-final'),
            'sdk_before':artifacts.get('sdk-initial'),'sdk_precompile':artifacts.get('sdk-precompile'),'sdk_prenative':artifacts.get('sdk-prenative'),'sdk_prereference':artifacts.get('sdk-prereference'),'sdk_after':artifacts.get('sdk-final'),
            'controller_runtime_before':artifacts.get('controller-runtime-before'),'controller_runtime_after':artifacts.get('controller-runtime-after'),
            'reference_runtime_before':artifacts.get('runtime-before-admitted'),'reference_runtime_after':artifacts.get('runtime-after-admitted'),
            'subprocesses':artifacts.get('children'),'native':parsed_receipt('native-stdout',capture,native_admitted),'reference':parsed_receipt('reference-output',capture,reference_admitted),
            'comparisons':internal['comparison_artifact'],'output_inventory':inventory,'gates':gates,'errors':errors,'qualification':qualification,
            'seal':{'status':'failed' if inventory is None else 'sealed-regular-output-files','file_permissions':{'regular':'0444','compiled_consumer':'0555'},'inventory_sha256':inventory_sha,'errors':[error for error in errors if error['stage'].startswith('seal')]}}


def execute(slug):
    if type(slug) is not str or SLUG.fullmatch(slug) is None: raise ValueError('single bounded attempt slug required')
    V.walk_directory(ROOT); V.walk_directory(PACKET)
    results=RESULTS
    if not results.exists(): results.mkdir(parents=True)
    V.walk_directory(results)
    attempt=results/slug; attempt.mkdir(mode=0o700,exist_ok=False); (attempt/'tmp').mkdir(mode=0o700)
    artifacts={}; buffers={}; capture=Capture(attempt,artifacts,buffers); runner=None; request=None; baseline=None; snapshot_packet=None
    runtime=None; own_before=None; sdk=None; loader=None; executable=None; hashes=None; comparisons=None; native=None; reference=None
    native_admitted=False; reference_admitted=False; phase='execution'
    record={'schema_version':1,'interface_id':'temporal-shared-optical-attempt/v1','status':'failed','numerical_accepted':False,
            'inference_status':'blocked-original-event-frame-reduction-calibration-selection-covariance-gaps',
            'qualification':'empirical-synthetic-fixed-temporal-shared-optical-control/v1','consumed_buffer_amendment_sha256':CONSUMED_AMENDMENT_SHA,
            'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'request_sha256':None,'failures':[],'artifacts':artifacts,'children':[],
            'source_roots':{'reproducible':str(ROOT),'sdk_source':None,'sdk_artifacts':None},'comparison_artifact':None,'comparison_accepted_before_terminal':False}
    def fail(stage,exc): record['failures'].append({'stage':stage,'kind':type(exc).__name__,'message':str(exc)[:2048],'details':getattr(exc,'details',None)})
    def guard(stage,call):
        try: return call()
        except BaseException as exc: fail(stage,exc); return None
    def hash_inputs():
        wanted=source_port_hashes(snapshot_packet,attempt/'contract.snapshot.json')
        V.exact(wanted,hashes,'snapshot source ports')
        V.exact(V.digest(V.read(attempt/'request.json')),record['request_sha256'],'attempt request')
        V.exact(source_port_hashes(PACKET,PACKET/'contract.snapshot.json'),hashes,'current source ports')
        V.exact(V.digest(V.read(PACKET/'request.json')),record['request_sha256'],'current request')
        return wanted
    try:
        if not sys.dont_write_bytecode: raise ValueError('controller requires Python -B or PYTHONDONTWRITEBYTECODE=1')
        raw=V.read(PACKET/'request.json'); record['request_sha256']=V.digest(raw)
        artifacts['request']=V.write_new(attempt/'request.json',raw)
        request=V.parse(raw); V.request_admission(request,PACKET)
        # Public experiment metadata binds the final executable request bytes.
        metadata=V.parse(V.read(PACKET/'experiment.json'))
        if metadata.get('request_sha256')!=V.digest(raw): raise ValueError('experiment metadata does not pin exact request bytes')
        record['request_sha256']=V.digest(raw); record['source_roots']['sdk_source']=request['sdk']['source_root']; record['source_roots']['sdk_artifacts']=request['sdk']['artifacts_root']
        artifacts['contract']=V.write_new(attempt/'contract.snapshot.json',V.read(PACKET/'contract.snapshot.json'))
        runner=Runner(attempt,request,artifacts,capture); record['children']=runner.children
        phase='source-admission'
        baseline=repro_map(ROOT,'initial',runner.launch); runner.artifact('repro-initial',baseline); admit_source(baseline); record['source_admitted']=True
        phase='execution'
        runner.launch('repro-source-archive',[GIT,*V.GIT_OPTIONS,'-C',str(ROOT),'archive','--format=tar','--output='+str(attempt/'source.tar'),baseline['head']],git=True)
        artifacts['source-archive']=unpack_snapshot(attempt,baseline)
        snapshot_packet=attempt/'source'/'experiments'/'temporal-shared-optical-control'
        hashes=source_port_hashes(snapshot_packet,attempt/'contract.snapshot.json'); hash_inputs()
        V.request_admission(request,snapshot_packet)
        phase='sdk-admission'
        sdk=V.fingerprint(request,'initial',runner.launch); runner.artifact('sdk-initial',sdk); V.admit_fingerprint(sdk); record['sdk_admitted']=True
        phase='execution'
        compiler=Path(request['sdk']['compiler']['path']).resolve()
        version=runner.launch('compiler-version',[str(compiler),'--version'])
        runner.artifact('compiler-version-identity',{'executable':V.file_identity(compiler,67108864)[1],'version_sha256':V.digest(version)})
        cli=Path(request['sdk']['artifacts_root'])/request['sdk']['cli']['path']
        discovery=V.parse(runner.launch('discovery',[str(cli),'describe','--json'],json_output=True))
        V.discovery_admission(discovery,request); runner.artifact('discovery-parsed',discovery)
        runtime=V.parse(runner.launch('runtime-before',[request['reference_runtime']['requested_executable'],str(snapshot_packet/'reference.py'),'--runtime-fingerprint'],json_output=True))
        V.runtime_admission(runtime,request); runner.artifact('runtime-before-admitted',runtime)
        own_before=controller_runtime(); runner.artifact('controller-runtime-before',own_before)
        hash_inputs()
        phase='sdk-admission'
        current=V.fingerprint(request,'precompile',runner.launch); runner.artifact('sdk-precompile',current); V.admit_fingerprint(current); V.exact(current,sdk,'SDK precompile')
        phase='execution'
        native_path=attempt/'consumer'
        args=[str(compiler),*request['sdk']['compiler_flags'],'-I'+str(Path(request['sdk']['source_root'])/'cpp/include'),str(snapshot_packet/'consumer.cpp'),str(Path(request['sdk']['artifacts_root'])/request['sdk']['archive']['path']),'-o',str(native_path)]
        runner.launch('compile',args); executable=V.file_identity(native_path,67108864)[1]; os.chmod(native_path,0o555); runner.artifact('compiled-consumer',executable)
        loader=loader_map(runner.launch('loader-before',[request['native_dependencies_proposed']['loader_path'],'--list',str(native_path)]),request); runner.artifact('loader-before-resolution',loader)
        phase='sdk-admission'
        current=V.fingerprint(request,'prenative',runner.launch); runner.artifact('sdk-prenative',current); V.admit_fingerprint(current); V.exact(current,sdk,'SDK prenative'); hash_inputs()
        phase='execution'
        V.exact(V.file_identity(native_path,67108864)[1],executable,'compiled native consumption')
        native=V.parse(runner.launch('native',[str(native_path)],json_output=True))
        phase='native-output-admission'
        native_complete=V.native_admission(native,request); native_admitted=native_complete
        if not native_complete: raise ValueError('native required owner refused; scientific comparison unavailable')
        phase='sdk-admission'
        current=V.fingerprint(request,'prereference',runner.launch); runner.artifact('sdk-prereference',current); V.admit_fingerprint(current); V.exact(current,sdk,'SDK prereference'); hash_inputs()
        phase='execution'
        reference_stdout=runner.launch('reference',[request['reference_runtime']['requested_executable'],str(snapshot_packet/'reference.py'),'--request',str(attempt/'request.json'),'--output',str(attempt/'reference.json')],json_output=True)
        summary=V.parse(reference_stdout); V.closed(summary,['interface_id','status','request_sha256','reference_json_sha256','reference_json_bytes'])
        V.exact(summary['interface_id'],'temporal-shared-optical-reference-summary/v1'); V.exact(summary['request_sha256'],record['request_sha256'])
        ref_raw=V.read(attempt/'reference.json'); reference=capture.retain('reference-output',ref_raw,attempt/'reference.json',parsed=True)
        V.exact(summary['reference_json_sha256'],V.digest(ref_raw)); V.exact(summary['reference_json_bytes'],len(ref_raw)); V.exact(summary['status'],reference['status'])
        phase='reference-output-admission'
        reference_complete=V.reference_admission(reference,request,record['request_sha256'],hashes,runtime); reference_admitted=reference_complete
        if not reference_complete: raise ValueError('reference refused; scientific comparison unavailable')
        phase='original169-comparisons'
        comparisons=compare(native,reference['payload']); V.exact(comparisons['check_count'],169,'original comparison count')
        record['comparison_artifact']=runner.artifact('comparisons',comparisons)
        record['comparison_accepted_before_terminal']=comparisons['accepted']
        if not comparisons['accepted']: raise ValueError('original numerical comparison gates refused')
        record['numerical_accepted']=True
    except BaseException as exc:
        fail(phase,exc)
        if runner is not None: runner.stop_cause=phase+': '+type(exc).__name__+': '+str(exc)[:2048]
    finally:
        # Recover independently written refusal/partial JSON only when not already consumed.
        # Never replace accepted in-memory values with a reread of mutable producer files.
        if 'reference-output' not in capture.buffers and (attempt/'reference.json').exists():
            def ingest_reference():
                blob=V.read(attempt/'reference.json'); return capture.retain('reference-output',blob,attempt/'reference.json',parsed=True)
            partial_ref=guard('reference-partial-capture',ingest_reference)
            if reference is None: reference=partial_ref
        if runner is not None:
            if native is None: native=capture.parsed.get('native-stdout')
            if native is not None:
                earned_native=guard('native-terminal-admission',lambda: V.native_admission(native,request))
                native_child=next((item for item in runner.children if item['label']=='native'),None)
                if earned_native is True and native_child is not None and native_child['status']=='completed' and capture.raw.get('native-stdout',{}).get('complete'):
                    native_admitted=True
            if executable is not None:
                def final_loader():
                    current=loader_map(runner.launch('loader-after-finally',[request['native_dependencies_proposed']['loader_path'],'--list',str(attempt/'consumer')]),request)
                    runner.artifact('loader-after-resolution',current); V.exact(current,loader,'loader dependency resolution before/after')
                guard('native-terminal-identity',final_loader)
                guard('compiled-consumer-terminal',lambda: V.exact(V.file_identity(attempt/'consumer',67108864)[1],executable,'compiled native terminal identity'))
            if runtime is not None and snapshot_packet is not None:
                def final_runtime():
                    current=V.parse(runner.launch('runtime-after-finally',[request['reference_runtime']['requested_executable'],str(snapshot_packet/'reference.py'),'--runtime-fingerprint'],json_output=True))
                    V.runtime_admission(current,request); runner.artifact('runtime-after-admitted',current); V.exact(current,runtime,'reference runtime preflight/terminal')
                guard('runtime-terminal-identity',final_runtime)
                def final_controller_runtime():
                    current=controller_runtime(); runner.artifact('controller-runtime-after',current); V.exact(current,own_before,'controller runtime before/after')
                guard('controller-runtime-terminal',final_controller_runtime)
            def final_sdk():
                current=V.fingerprint(request,'final',runner.launch); runner.artifact('sdk-final',current)
                V.admit_fingerprint(current)
                if sdk is not None: V.exact(current,sdk,'SDK final complete map')
            guard('sdk-terminal-identity',final_sdk)
            guard('sdk-artifact-terminal',lambda: V.artifact_checks(request))
            def final_source():
                current=repro_map(ROOT,'final',runner.launch,baseline['head'] if baseline is not None else None); runner.artifact('repro-final',current); admit_source(current)
                if baseline is not None: V.exact(current,baseline,'committed source before/after')
                if hashes is not None: hash_inputs()
                if snapshot_packet is not None and baseline is not None:
                    for name,item in baseline['sources'].items():
                        blob=V.read(attempt/'source'/name); V.exact(V.digest(blob),item['sha256'],'archived source terminal')
                V.exact(V.file_identity(GIT,67108864)[1],runner.tool_baseline,'Git executable before/after')
            guard('source-terminal-identity',final_source)
            if hashes is not None: guard('input-terminal-identity',hash_inputs)
            if reference is not None and runtime is not None and hashes is not None:
                earned_reference=guard('reference-terminal-admission',lambda: V.reference_admission(reference,request,record['request_sha256'],hashes,runtime))
                reference_child=next((item for item in runner.children if item['label']=='reference'),None)
                if earned_reference is True and reference_child is not None and reference_child['status']=='completed' and capture.raw.get('reference-output',{}).get('complete'):
                    reference_admitted=True
            guard('store-terminal',lambda: runner.room())
        apply_capture_gate(record,capture)
        record['ended_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
        record['raw_captures']=capture.raw
        record['retained_native_object']=artifacts.get('native-stdout')
        record['retained_reference_object']=artifacts.get('reference-output')
        record['seal']='regular files0444; executable0555; source directories555; exclusive new terminal record'
        # Full child records are sealed independently; compact terminal contains identities.
        if runner is not None:
            guard('children-receipt',lambda: runner.artifact('children',child_records(runner,capture,native_admitted,reference_admitted)))
            record['children']={'artifact':artifacts.get('children'),'direct_launches':sum(child['launched'] for child in runner.children),'reserved_schedule_labels':len(runner.used),'scope':request['resources']['child_count_scope']}
            record['direct_launches']=sum(child['launched'] for child in runner.children)
            record['execution_complete']=len(runner.children)==28 and all(child['status']=='completed' for child in runner.children)
            for name,item in list(artifacts.items()):
                def check_artifact(item=item):
                    path=attempt/item['path']; blob=V.read(path,4194304 if path.name=='source.tar' else V.MIB)
                    V.exact(len(blob),item['bytes'],'artifact bytes'); V.exact(V.digest(blob),item['sha256'],'artifact hash')
                guard('artifact-terminal-'+name,check_artifact)
        else:
            def unstarted_children():
                cause=phase+': '+record['failures'][0]['message'] if record['failures'] else 'request/source admission unavailable'
                blob=V.encode([not_started_record(label,cause) for label in SCHEDULE])
                # No admitted runner exists; the frozen literal store cap still applies.
                total=sum(path.lstat().st_size for path in attempt.rglob('*') if stat.S_ISREG(path.lstat().st_mode))
                if total+len(blob)+65536>268435456: raise StoreLimit('unstarted child receipt store quota')
                artifacts['children']=V.write_new(attempt/'children.json',blob)
            guard('children-receipt',unstarted_children)
        if record['failures']: record['numerical_accepted']=False; record['status']='failed'
        guard('seal-regular-files',lambda: seal(attempt,False))
        inventory=guard('seal-output-inventory',lambda: output_inventory(attempt))
        if record['failures']: record['numerical_accepted']=False; record['status']='failed'
        final=final_receipt(record,slug,artifacts,capture,inventory,native_admitted,reference_admitted)
        result=terminal_record(attempt,final,artifacts,runner.room if runner else lambda extra:None,runner.total_room if runner else None)
        seal(attempt)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument('--attempt',required=True)
    args=parser.parse_args(); result=execute(args.attempt)
    accepted=result.get('gates',{}).get('numerical',{}).get('status')=='passed'
    print(json.dumps({'status':'completed' if accepted else 'failed','numerical_accepted':accepted,'record':str(RESULTS/args.attempt/'record.json')},sort_keys=True))
    return 0 if accepted else 1


if __name__=='__main__':
    raise SystemExit(main())
