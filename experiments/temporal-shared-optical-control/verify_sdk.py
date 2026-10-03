#!/usr/bin/env python3
"""Admission helpers for the one frozen temporal/shared-optical SDK control.

No numerical kernels, process launcher, installation or CLI inference lives here.
The controller supplies its bounded, recorded launcher for the fixed Git calls.
"""
from decimal import Decimal
import hashlib
import json
import os
from pathlib import Path
import re
import stat

MIB = 1048576
CONTRACT_SHA = 'e3a1ceb57d96f599cf82d938462d03feef5d424ad9d2c70f82316d63c34ba6d7'
CONSUMER_SHA = '46e17d76d0b03390048fdc66ba97ca7012d8eda215beb9ab2aeea96a92905ab8'
DECIMAL_PATTERN = re.compile(r'-?(0|[1-9][0-9]*)(\.[0-9]+)?([eE][+-]?[0-9]+)?\Z')
SHA_PATTERN = re.compile(r'[0-9a-f]{64}\Z')
STATUSES = {'ok', 'invalid_input', 'nonfinite_input', 'overflow', 'work_limit',
            'outside_domain', 'singular', 'not_positive_definite', 'conditioning_budget_exceeded'}
GIT_OPTIONS = ['-c', 'core.fsmonitor=false', '-c', 'core.preloadIndex=false',
               '-c', 'index.threads=1', '--no-pager']


class IdentityError(ValueError):
    def __init__(self,message,details):
        super().__init__(message); self.details=details


def exact(actual, expected, label='identity'):
    if type(actual) is not type(expected) or actual != expected:
        raise ValueError(label + ' differs in type/value')
    if type(expected) is dict:
        for key in expected: exact(actual[key], expected[key], label + '/' + key)
    elif type(expected) is list:
        for index, (a, b) in enumerate(zip(actual, expected)):
            exact(a, b, label + '/' + str(index))


def closed(value, fields, label='object'):
    if type(value) is not dict or set(value) != set(fields):
        raise ValueError(label + ' closed fields differ')
    return value


def uint(value, maximum=4294967295):
    if type(value) is not int or not 0 <= value <= maximum:
        raise ValueError('bounded integer required')
    return value


def number(value, lower=None, upper=None):
    if type(value) not in (int, Decimal):
        raise ValueError('exact JSON numeric token required')
    result = Decimal(value)
    if not result.is_finite() or (lower is not None and result < lower) or (upper is not None and result > upper):
        raise ValueError('finite numeric domain differs')
    return result


def decimal_string(value, nonnegative=False):
    if type(value) is not str or len(value) > 256 or DECIMAL_PATTERN.fullmatch(value) is None:
        raise ValueError('strict finite decimal string required')
    result = Decimal(value)
    if not result.is_finite() or abs(result.adjusted()) > 1000 or (nonnegative and result < 0):
        raise ValueError('reference decimal domain differs')
    return result


def sha(value):
    if type(value) is not str or SHA_PATTERN.fullmatch(value) is None:
        raise ValueError('SHA256 required')
    return value


def array(value, size=None, maximum=None):
    if type(value) is not list or (size is not None and len(value) != size) or (maximum is not None and len(value) > maximum):
        raise ValueError('ordered array size differs')
    return value


def status(value):
    if value not in STATUSES or type(value) is not str:
        raise ValueError('unknown native status')


def boolean(value):
    if type(value) is not bool:
        raise ValueError('boolean required')


def parse(blob):
    if type(blob) is not bytes or len(blob) > MIB:
        raise ValueError('bounded JSON bytes required')
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError('duplicate JSON key: ' + key)
            result[key] = value
        return result
    def token(text):
        if len(text) > 256:
            raise ValueError('numeric token limit')
        value = Decimal(text)
        if not value.is_finite() or abs(value.adjusted()) > 1000:
            raise ValueError('nonfinite/unbounded numeric token')
        return value
    def integer(text):
        if len(text) > 256:
            raise ValueError('integer token limit')
        return int(text)
    def bad_constant(text):
        raise ValueError('nonfinite JSON constant: ' + text)
    value = json.loads(blob.decode('utf-8', errors='strict'), object_pairs_hook=pairs,
                       parse_float=token, parse_int=integer, parse_constant=bad_constant)
    count = 0
    def bound(item, depth):
        nonlocal count
        count += 1
        if count > 65536 or depth > 32:
            raise ValueError('JSON structure quota')
        if type(item) is str and len(item.encode('utf-8')) > 8192:
            raise ValueError('JSON string quota')
        if type(item) is dict:
            for key, child in item.items():
                bound(key, depth + 1); bound(child, depth + 1)
        elif type(item) is list:
            for child in item:
                bound(child, depth + 1)
    bound(value, 0)
    return value


def encode(value):
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False,
                       default=lambda item: str(item) if type(item) is Decimal else _unsupported(item)) + '\n').encode('utf-8')


def _unsupported(item):
    raise TypeError('unsupported receipt type ' + type(item).__name__)


def digest(blob):
    return hashlib.sha256(blob).hexdigest()


def read(path, limit=MIB):
    path = Path(path)
    walk_directory(path.parent)
    descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    try:
        info = os.fstat(descriptor)
        if not stat.S_ISREG(info.st_mode) or info.st_size > limit:
            raise ValueError('regular bounded file required: ' + str(path))
        chunks, total = [], 0
        while True:
            chunk = os.read(descriptor, min(65536, limit + 1 - total))
            if not chunk:
                break
            chunks.append(chunk); total += len(chunk)
            if total > limit:
                raise ValueError('file grew beyond bound: ' + str(path))
        after = os.fstat(descriptor)
        current = path.stat(follow_symlinks=False)
        key = lambda item: (item.st_dev, item.st_ino, item.st_size, item.st_mtime_ns)
        if key(info) != key(after) or key(after) != key(current) or total != after.st_size:
            raise ValueError('file changed during read: ' + str(path))
        return b''.join(chunks)
    finally:
        os.close(descriptor)


def file_identity(path, limit=MIB, expected=None):
    blob = read(path, limit)
    identity = {'path': str(Path(path).absolute()), 'resolved_path': str(Path(path).resolve()),
                'bytes': len(blob), 'sha256': digest(blob)}
    if expected is not None and identity['sha256'] != expected:
        raise IdentityError('file hash differs: ' + str(path),{'expected_sha256':expected,'actual':identity})
    return blob, identity


def walk_directory(path):
    path = Path(path).absolute()
    for part in [path, *path.parents]:
        info = part.lstat()
        if not stat.S_ISDIR(info.st_mode) or stat.S_ISLNK(info.st_mode):
            raise ValueError('nonsymlink directory required: ' + str(part))
    return path


def write_new(path, blob, limit=MIB, mode=0o444):
    if type(blob) is not bytes or len(blob) > limit:
        raise ValueError('bounded output bytes required')
    path = Path(path); walk_directory(path.parent)
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
    try:
        view = memoryview(blob)
        while view:
            size = os.write(descriptor, view)
            if size <= 0:
                raise OSError('short output write')
            view = view[size:]
        os.fsync(descriptor); os.fchmod(descriptor, mode)
    finally:
        os.close(descriptor)
    written=read(path,limit)
    if len(written)!=len(blob) or digest(written)!=digest(blob): raise ValueError('new output bytes differ after write')
    directory=os.open(path.parent,os.O_RDONLY|os.O_DIRECTORY|os.O_NOFOLLOW)
    try: os.fsync(directory)
    finally: os.close(directory)
    return {'path': path.name, 'bytes': len(blob), 'sha256': digest(blob)}


def request_admission(value, packet):
    expected = parse(REQUEST_PROTOTYPE_TEXT.encode('utf-8'))
    expected.update(request_id='temporal-shared-optical-reproducible-request/v1',
                    interface_id='temporal-shared-optical-bounded-controller/v1',
                    status='reviewed-native-synthetic-control', execution=None)
    closed(value, expected, 'request')
    for name, filename in [('reference', 'reference.py'), ('controller', 'controller.py'), ('sdk_verifier', 'verify_sdk.py')]:
        candidate = sha(value['source_ports'][name]['sha256'])
        file_identity(Path(packet) / filename, MIB, candidate)
        expected['source_ports'][name]['sha256'] = candidate
    if value != expected:
        raise ValueError('frozen request contract differs')
    # Recursive exact types exclude bool masquerading as integer even under equality.
    def typed(actual, wanted):
        exact(type(actual).__name__, type(wanted).__name__, 'request scalar/container type')
        if type(actual) is dict:
            closed(actual, wanted)
            for key in wanted: typed(actual[key], wanted[key])
        elif type(actual) is list:
            array(actual, len(wanted))
            for a, b in zip(actual, wanted): typed(a, b)
        else: exact(actual, wanted)
    typed(value, expected)
    file_identity(Path(packet) / 'contract.snapshot.json', MIB, CONTRACT_SHA)
    file_identity(Path(packet) / 'consumer.cpp', MIB, CONSUMER_SHA)
    return value


def fingerprint(request, phase, run):
    sdk=request['sdk']; root=Path(sdk['source_root']); errors=[]
    def collect(label,call):
        try: return call()
        except Exception as exc:
            errors.append({'stage':label,'kind':type(exc).__name__,'message':str(exc)[:2048],'details':getattr(exc,'details',None)})
            return None
    collect('source-root',lambda: walk_directory(root))
    head=collect('HEAD',lambda:run('sdk-'+phase+'-head',['git',*GIT_OPTIONS,'-C',str(root),'rev-parse','--verify','HEAD'],git=True).decode().strip())
    raw_status=collect('tracked-status',lambda:run('sdk-'+phase+'-status',['git',*GIT_OPTIONS,'-C',str(root),'status','--porcelain=v1','-z','--untracked-files=no'],git=True))
    collect('HEAD-pin',lambda:exact(head,sdk['revision']))
    if raw_status is not None: collect('tracked-status-pin',lambda:exact(raw_status,b''))
    manifest_file=collect('manifest-file',lambda:file_identity(root/sdk['manifest']['path'],MIB))
    manifest=collect('manifest-JSON',lambda:parse(manifest_file[0])) if manifest_file else None
    manifest_identity=manifest_file[1] if manifest_file else None
    if manifest_identity is not None: collect('manifest-pin',lambda:exact(manifest_identity['sha256'],sdk['manifest']['sha256']))
    if manifest is not None and type(manifest) is not dict:
        errors.append({'stage':'manifest-object','kind':'ValueError','message':'manifest object required','details':None}); manifest=None
    build_id=None; expected_sources={}
    if manifest is not None:
        collect('manifest-HEAD',lambda:exact(manifest['git_head'],head)); collect('manifest-status',lambda:exact(manifest['git_status'],''))
        build={k:v for k,v in manifest.items() if k not in ('build_id','git_head','git_status')}
        build_id=collect('build-encoding',lambda:digest(json.dumps(build,sort_keys=True,separators=(',',':'),allow_nan=False).encode()))
        collect('build-pin',lambda:exact(build_id,sdk['build_id'])); collect('manifest-build',lambda:exact(manifest['build_id'],build_id))
        expected_sources=manifest.get('sources',{})
        if type(expected_sources) is not dict: errors.append({'stage':'manifest-sources','kind':'ValueError','message':'source map required','details':None}); expected_sources={}
        if 'compiler_executable_digest' in manifest: collect('compiler-manifest',lambda:exact(manifest['compiler_executable_digest'],sdk['compiler']['sha256']))
        if 'standard_library_digest' in manifest: collect('stdlib-manifest',lambda:exact(manifest['standard_library_digest'],sdk['standard_library']['sha256']))
    paths=set(sdk['source_explicit'])
    def enumerate_paths():
        for pattern in sdk['source_patterns']:
            for path in root.glob(pattern):
                paths.add(str(path.relative_to(root)))
                if len(paths)>512: raise ValueError('SDK source enumeration bound')
    collect('source-enumeration',enumerate_paths)
    collect('source-closed-set',lambda:exact(paths,set(expected_sources))); collect('source-count',lambda:exact(len(paths),311))
    sources={}
    for name in sorted(paths):
        relative=Path(name)
        def source_file(relative=relative):
            if relative.is_absolute() or '..' in relative.parts: raise ValueError('escaped SDK source')
            return file_identity(root/relative,MIB)[1]
        sources[name]=collect('source-file/'+name,source_file)
        if sources[name] is not None: collect('source-pin/'+name,lambda name=name:exact(sources[name]['sha256'],expected_sources.get(name)))
    headers={}
    def enumerate_headers():
        for path in (root/'cpp/include').rglob('*'):
            if path.is_dir() and not path.is_symlink(): continue
            name=str(path.relative_to(root))
            if len(headers)>=512: raise ValueError('header inventory bound')
            headers[name]=collect('header/'+name,lambda path=path:file_identity(path,MIB)[1])
    collect('header-enumeration',enumerate_headers)
    documentary={}
    for entry in sdk['headers']+sdk['guides']:
        name=entry['path']; documentary[name]=collect('documentary/'+name,lambda name=name:file_identity(root/name,MIB)[1])
        if documentary[name] is not None: collect('documentary-pin/'+name,lambda name=name,entry=entry:exact(documentary[name]['sha256'],entry['sha256']))
    artifacts={}
    for name in ('manifest','archive','cli','compiler','standard_library'):
        item=sdk[name]; path=Path(item['path'])
        if not path.is_absolute(): path=Path(sdk['artifacts_root'])/path
        if name in ('compiler','standard_library'): path=path.resolve()
        artifacts[name]=collect('artifact/'+name,lambda path=path:file_identity(path,268435456)[1])
        if artifacts[name] is not None: collect('artifact-pin/'+name,lambda name=name,item=item:exact(artifacts[name]['sha256'],item['sha256']))
    loader=request['native_dependencies_proposed']; artifacts['loader']=collect('loader',lambda:file_identity(Path(loader['loader_path']).resolve(),67108864)[1])
    if artifacts['loader'] is not None: collect('loader-pin',lambda:exact(artifacts['loader']['sha256'],loader['loader_sha256']))
    return {'head':head,'tracked_status':raw_status.decode('utf-8','backslashreplace') if raw_status is not None else None,'build_id':build_id,
            'manifest':manifest_identity,'sources':sources,'headers':headers,'guides_and_selected_headers':documentary,'artifacts':artifacts,'errors':errors}


def admit_fingerprint(value):
    if value['errors']: raise IdentityError('SDK complete-map admission refused',{'errors':value['errors']})
    return value


def artifact_checks(request):
    sdk = request['sdk']; root = Path(sdk['artifacts_root'])
    identities = {}
    for name in ('manifest', 'archive', 'cli', 'compiler', 'standard_library'):
        item = sdk[name]; path = Path(item['path'])
        if not path.is_absolute(): path = root / path
        # Explicit requested compiler/stdlib aliases are pinned resolved identities.
        if name in ('compiler', 'standard_library'): path = path.resolve()
        identities[name] = file_identity(path, 268435456, item['sha256'])[1]
    exact(identities['compiler']['sha256'], request['sdk']['compiler']['sha256'])
    loader = request['native_dependencies_proposed']
    identities['loader'] = file_identity(Path(loader['loader_path']).resolve(), 67108864, loader['loader_sha256'])[1]
    return identities


# Exact source-approved request template. Only the four approved metadata values
# and the three final reviewed script SHAs differ in the executable request.
REQUEST_PROTOTYPE_TEXT = r'''
{
  "schema_version": 1,
  "request_id": "temporal-shared-optical-reproducible-request/v2-proposal",
  "status": "source-only-proposal-unexecutable-until-parent-freeze",
  "execution": null,
  "packet": "temporal-shared-optical-control",
  "interface_id": "temporal-shared-optical-bounded-controller/v2-proposal",
  "origin": {
    "kind": "development_smoke",
    "component_repository": "Irreducible",
    "component_revision": "bd2dff2d7d84fec45e384326d8e0e9d2aa1e9d09",
    "component_root": "/home/szymon/.codex/worktrees/all19-optical/irreducible",
    "source_directory": "evidence/next-wave-temporal-optical-sdk",
    "source_tracking": "ignored component files; each exact byte identity is independent of tracked HEAD",
    "admission_path": ".work/temporal-optical-admission-20261003/admission.json",
    "admission_sha256": "1ee77abd937be8342f3be946c17453a162f75e9cd315a9a53b05f8fe2db401cc",
    "preserved_directory": ".work/temporal-optical-admission-20261003/preserved",
    "original_sources": [
      {
        "path": "contract-v2.json",
        "bytes": 26350,
        "sha256": "e3a1ceb57d96f599cf82d938462d03feef5d424ad9d2c70f82316d63c34ba6d7",
        "role": "exact current original design copied as contract.snapshot.json"
      },
      {
        "path": "consumer.cpp",
        "bytes": 17548,
        "sha256": "46e17d76d0b03390048fdc66ba97ca7012d8eda215beb9ab2aeea96a92905ab8",
        "role": "unchanged native SDK consumer and external reducer"
      },
      {
        "path": "reference.py",
        "bytes": 12857,
        "sha256": "f912a6f3b0b351b91ece9d26af959200618776908d1842bb64caa63ed87ebe0e",
        "role": "original scientific arithmetic for I/O-only fork review"
      },
      {
        "path": "controller.py",
        "bytes": 16090,
        "sha256": "4f76eca4f9ec4f01ead57978a7031f6dbf28e07eaf14c7e7e286424a1412e535",
        "role": "original169 comparison logic and budgets for semantic preservation"
      },
      {
        "path": "verify_sdk.py",
        "bytes": 4566,
        "sha256": "3354c6fbb7b8b2ec201351b859a3b3e3774cead6db17bb27a015e7bebcd70a68",
        "role": "original full311 SDK fingerprint review ancestry"
      }
    ],
    "historical_result": {
      "attempt": "attempt03",
      "result_sha256": "cc2d51012ee7167205719adcd99d6eea040e00fec696590b4b39c7dcbc26a740",
      "comparison_sha256": "35ddcde6d6cc6b94901db69c5f4d77a283d411a34e218b02e69594fdef3bcc2d",
      "reported_checks": 169,
      "new_execution_evidence": false
    }
  },
  "original_contract": {
    "path": "contract.snapshot.json",
    "bytes": 26350,
    "sha256": "e3a1ceb57d96f599cf82d938462d03feef5d424ad9d2c70f82316d63c34ba6d7",
    "all_numerical_model_endpoint_budgets": "unchanged exact original contract bytes"
  },
  "source_ports": {
    "consumer": {
      "path": "consumer.cpp",
      "bytes": 17548,
      "sha256": "46e17d76d0b03390048fdc66ba97ca7012d8eda215beb9ab2aeea96a92905ab8",
      "input": "no arguments or input JSON; exact source declares synthetic design",
      "output_schema": "external-temporal-optical-sdk-control/v1"
    },
    "reference": {
      "path": "reference.py",
      "sha256": null,
      "interface_id": "temporal-shared-optical-reference-transport/v1",
      "scientific_payload_schema": "original-high-precision-temporal-optical-reference/v1",
      "change_scope": "transport/runtime/failure accounting only; original arithmetic unchanged"
    },
    "controller": {
      "path": "controller.py",
      "sha256": null,
      "arguments": [
        "--attempt",
        "SLUG"
      ],
      "comparison_source_ancestry_sha256": "4f76eca4f9ec4f01ead57978a7031f6dbf28e07eaf14c7e7e286424a1412e535",
      "comparison_semantics": "original169 checks and all budgets unchanged"
    },
    "sdk_verifier": {
      "path": "verify_sdk.py",
      "sha256": null,
      "change_scope": "bounded full inventory/raw CLI/runtime/source admission only"
    },
    "closed_output_contract": {
      "source_owner": "/root/reference_owner",
      "documentary_proposal_sha256": "791ea3981b8eef42eec3a5622e0193cf8ecfe06a5f28472dda032c50f1f1c266",
      "embedded_schema_key": "shape_contract",
      "ignored_proposal_file_required_at_runtime": false
    },
    "pin_rule": "final request pins all new scripts; experiment metadata pins final request bytes; source admission requires clean committed blobs, avoiding cyclic script/request hash literals"
  },
  "sdk": {
    "source_root": "/home/szymon/.codex/worktrees/box-target-integrated/irreducible",
    "artifacts_root": "/home/szymon/.codex/worktrees/box-target-integrated/irreducible",
    "revision": "6f869532c1951ed1afd9f2506b5d05c6cfd03c82",
    "tracked_status": "clean",
    "build_id": "f22c25423cfb9cbac3c2b91a4e514b13ce604e92e7010f55a9aa0bdd42f40f59",
    "manifest": {
      "path": "build/build-manifest-release.json",
      "sha256": "f2e4d6a22257c13f8ab46de51cdfc92654d9857a82a7bc3ea5b8c0edf234f062"
    },
    "archive": {
      "path": "build/native-release/libirred_core.a",
      "sha256": "1cb2b85ad292334f3b5041d669187b04a6fca0916bca8c41070b0ff878499dff"
    },
    "cli": {
      "path": "target/release/irred",
      "sha256": "09b5bfd05ce6c057423b1d41f6f81f8e4d0db44dfd35424fdccd0b34fd9b22d9",
      "discovery": "fixed describe --json; no temporal operation inferred"
    },
    "compiler": {
      "path": "/usr/bin/c++",
      "sha256": "f04191f6a7b2cd7d9a62e1745872b8a6088791e5af6955488c69c9b2c4668bc9"
    },
    "standard_library": {
      "path": "/usr/lib/gcc/x86_64-pc-linux-gnu/16/libstdc++.so",
      "sha256": "f5fc7380f2ae46fa4053a64be04e7b98109f1066a4bbfff3c37042488aa0be0e"
    },
    "source_inventory_count": 311,
    "source_patterns": [
      "src/**/*.rs",
      "tests/**/*.rs",
      "cpp/**/*.cpp",
      "cpp/**/*.h",
      "cpp/**/*.hpp",
      "cpp/**/*.inc",
      "cpp/**/*.cmake",
      "schema/*.json",
      "tools/*.py"
    ],
    "source_explicit": [
      "Cargo.toml",
      "Cargo.lock",
      "build.rs",
      "cpp/CMakeLists.txt"
    ],
    "headers": [
      {
        "path": "cpp/include/irred/temporal_photometry.hpp",
        "sha256": "78df6251066b9272482f79a43acd3f4a48a7128e2594e8d5e4004a7d198ab151"
      },
      {
        "path": "cpp/include/irred/detector_selection.hpp",
        "sha256": "2973496bcf85f51ef42a1f5ab2be319019aa34076d7603b8e2d6eb19ea5a1dbd"
      },
      {
        "path": "cpp/include/irred/photometry_calibration.hpp",
        "sha256": "0a825b6fd7bb659491801dec667e08e3c52b739ea2c36dff1c4a55113a1980fd"
      },
      {
        "path": "cpp/include/irred/photometry.hpp",
        "sha256": "0361c9e43696bbc49791e4b9b1331da6fcae430273cdb2b71dad3a74da64dec4"
      }
    ],
    "guides": [
      {
        "path": "docs/temporal-photometry.md",
        "sha256": "68b93d0278bd148970a9656fb49f84aa040dfadd3fc656c04540c416c8a6ef67"
      },
      {
        "path": "docs/detector-selection.md",
        "sha256": "b02d34e08d0cdc8c96be7f1543011caf90adc9a5841fe42edcc4513226484262"
      },
      {
        "path": "docs/optical-detector.md",
        "sha256": "af5fb3ca260222146ae669decbecb567ff9e63fae6fe2e59ab64479504af9017"
      },
      {
        "path": "docs/photometry-calibration.md",
        "sha256": "98bc0a8ea98ef1e04468b86dde217954bc75997026e30bbcc164267bfbb6125e"
      }
    ],
    "compiler_flags": [
      "-std=c++20",
      "-O3",
      "-DNDEBUG",
      "-Wall",
      "-Wextra",
      "-Wpedantic",
      "-fno-fast-math",
      "-ffp-contract=off"
    ],
    "verification": "complete source/header/build/artifact maps at initial admission, precompile, prenative, prereference and finally; exact original before-each-numerical boundaries; required artifact/source/tool hash checks at consumption boundaries; actual compiled consumer identity before/after"
  },
  "reference_runtime": {
    "requested_executable": "/home/szymon/Projects/irreducible/evidence/project-review/science/recovery-abundance-growth-20261002/reference-venv/bin/python",
    "resolved_executable": "/home/szymon/.local/share/mise/installs/python/3.14.8/bin/python3.14",
    "executable_sha256": "815b1275bf87e7595fe1cdf3f2efc6ae58e27b842262a7b530de8610ac2760ba",
    "implementation": "CPython",
    "version": "3.14.8 (main, Oct  1 2026, 20:55:04) [GCC 16.2.1 20260810]",
    "prefix": "/home/szymon/Projects/irreducible/evidence/project-review/science/recovery-abundance-growth-20261002/reference-venv",
    "base_prefix": "/home/szymon/.local/share/mise/installs/python/latest",
    "pyvenv_cfg_sha256": "11a0241468a7d4cb44419cd7c7a87cd6149c3e90e38b75be6e1a201279621e43",
    "mpmath_version": "1.3.0",
    "mpmath_backend": "python",
    "mpmath_init_sha256": "b241584d2c1fc0304b0a1015ea923749d7b0800411dd406dcab7c82bf25d9fe8",
    "gmpy2_loaded": false,
    "original_seal_sha256": "ae4e395559ac113c6d7ecbbe165128509d350e4b5fa5d9656251e19b0f096745",
    "legacy_package_inventory_sha256": "982dc09f2ddf1206cc9cd4467c7be1ce4238a612489b43bf5f5d196a8b2495a8",
    "legacy_inventory_hash_encoding": "original seal reference_runtime.source_inventory absolute-path-to-sha map; sorted keys compact JSON UTF-8, no trailing LF",
    "expanded_inventory_rule": "complete metadata/scientific-code module imports before first actual map, no scientific calls; independently equal complete imported_modules and mapped_libraries sets plus inventory/config before/after; no numerical prewarm or installation"
  },
  "native_dependencies_proposed": {
    "enabled": true,
    "scope": "ELF loader resolved dependencies under scrubbed environment, not actual native process mapped inventory",
    "loader_path": "/usr/lib/ld-linux-x86-64.so.2",
    "loader_sha256": "f5e11cc62f8f2c24dff532982559f4303dba3b17e7450551372a88c8e7ea8757",
    "boundaries": [
      "before native",
      "terminal finally"
    ],
    "freeze_decision": "parent approved in principle; final literal invocation/environment/schema reviewed before executable request freeze"
  },
  "resources": {
    "jobs": 1,
    "threads": 1,
    "compile_cpu_seconds": 180,
    "compile_wall_seconds": 180,
    "native_cpu_seconds": 180,
    "native_wall_seconds": 180,
    "reference_cpu_seconds": 180,
    "reference_wall_seconds": 180,
    "structural_child_cpu_seconds": 180,
    "structural_child_wall_seconds": 180,
    "git_metadata_cpu_seconds": 30,
    "git_metadata_wall_seconds": 30,
    "git_metadata_regular_file_limit_bytes": 1048576,
    "git_archive_regular_file_limit_bytes": 4194304,
    "git_archive_stdout_bytes": 1048576,
    "git_archive_independent_output_bytes": 4194304,
    "address_space_bytes_each_child": 1073741824,
    "stdout_bytes_each_child": 1048576,
    "stderr_bytes_each_child": 1048576,
    "reference_json_bytes": 1048576,
    "compiled_executable_bytes": 67108864,
    "compile_regular_file_limit_bytes": 67108864,
    "native_reference_regular_file_limit_bytes": 1048576,
    "source_tree_bytes": 2097152,
    "source_archive_bytes": 4194304,
    "source_file_bytes": 1048576,
    "attempt_store_bytes": 268435456,
    "receipt_file_bytes": 1048576,
    "terminal_fallback_bytes": 65536,
    "terminal_fallback_artifact_count": 128,
    "terminal_artifact_relative_path_bytes": 256,
    "terminal_fallback_message_bytes": 2048,
    "runtime_inventory_files": 4096,
    "runtime_inventory_total_bytes_scanned": 67108864,
    "maximum_children": 28,
    "child_count_scope": "DIRECT controller launches; compiler descendants are one serial build job under process-group wall/store kill; RLIMIT_CPU is per-process, not aggregate descendant CPU",
    "git_metadata_children": 19,
    "non_git_children": 9,
    "maximum_aggregate_child_wall_seconds": 2190,
    "aggregate_wall_scope": "sum9*180+19*30 direct-child watchdog budgets, excluding bounded in-process hashing/serialization; not a whole-controller wall assertion",
    "cpu_measurement_scope": "actual wait4 rusage of reaped child as attributed by OS; do not claim an enforced aggregate compiler-tree CPU180 limit",
    "retry_count": 0,
    "schedule": [
      "repro-initial-head",
      "repro-initial-branch",
      "repro-initial-status",
      "repro-initial-tree",
      "repro-source-archive",
      "sdk-initial-head",
      "sdk-initial-status",
      "compiler-version",
      "discovery",
      "runtime-before",
      "sdk-precompile-head",
      "sdk-precompile-status",
      "compile",
      "loader-before",
      "sdk-prenative-head",
      "sdk-prenative-status",
      "native",
      "sdk-prereference-head",
      "sdk-prereference-status",
      "reference",
      "loader-after-finally",
      "runtime-after-finally",
      "sdk-final-head",
      "sdk-final-status",
      "repro-final-head",
      "repro-final-branch",
      "repro-final-status",
      "repro-final-tree"
    ],
    "git_fixed_global_options": [
      "-c",
      "core.fsmonitor=false",
      "-c",
      "core.preloadIndex=false",
      "-c",
      "index.threads=1",
      "--no-pager"
    ],
    "git_fixed_environment": {
      "GIT_OPTIONAL_LOCKS": "0",
      "GIT_CONFIG_NOSYSTEM": "1"
    },
    "metadata_ledger_rule": "literal controller operations described in PROPOSAL.md, never execute commands from JSON; one archive and no per-file git show/other helper launches",
    "environment": {
      "MKL_NUM_THREADS": "1",
      "NUMEXPR_NUM_THREADS": "1",
      "OMP_NUM_THREADS": "1",
      "OPENBLAS_NUM_THREADS": "1",
      "BLIS_NUM_THREADS": "1",
      "VECLIB_MAXIMUM_THREADS": "1",
      "PYTHONDONTWRITEBYTECODE": "1"
    },
    "environment_policy": "fixed PATH/tool paths, scrub LD_PRELOAD/LD_LIBRARY_PATH/PYTHONPATH/PYTHONHOME and unrelated numerical overrides; private compile TMPDIR inside attempt; actual effective environment recorded",
    "native_payload_policy_bytes": 1048576,
    "native_and_reference_math_policy": "unchanged original exact contract"
  },
  "transport": {
    "output_root": "results/temporal-shared-optical-control",
    "attempt_slug_pattern": "[A-Za-z0-9][A-Za-z0-9_-]{0,63}",
    "fresh_only": true,
    "reference_output_name": "reference.json",
    "source_layout": {
      "committed_source_directory": "source",
      "packet_directory": "source/experiments/temporal-shared-optical-control",
      "immediate_request": "request.json",
      "immediate_contract": "contract.snapshot.json",
      "reference_script": "source/experiments/temporal-shared-optical-control/reference.py",
      "contract_relative_to": "request parent",
      "source_ports_relative_to": "approved source/packet directory"
    },
    "maximum_json_depth": 32,
    "maximum_json_nodes": 65536,
    "maximum_string_bytes": 8192,
    "maximum_numeric_token_bytes": 256,
    "maximum_decimal_adjusted_exponent_absolute": 1000,
    "native_float_parser": "finite Decimal from exact emitted JSON number tokens; no binary64 intermediate",
    "reference_scalar_parser": "strict finite decimal strings retaining original60 digits",
    "decimal_string_pattern": "-?(0|[1-9][0-9]*)(\\.[0-9]+)?([eE][+-]?[0-9]+)?",
    "decimal_regex_encoding_rule": "JSON source contains two backslashes before dot, decoded regex exactly one backslash then dot (UTF8hex5c2e); decoded grammar is authoritative",
    "duplicate_keys": "reject",
    "nonfinite_constants": "reject",
    "source_admission": "clean committed repository HEAD and exact current/snapshot source blobs",
    "capture": "read raw stdout/stderr/reference JSON once bounded, hash and parse same bytes; preserve all outcomes; truncated overflow prefixes never impersonate complete streams",
    "terminal": "always rehash captured streams/files and all source/SDK/runtime/input/executable identities; raw drift revokes acceptance without replacing original parsed objects/checks",
    "terminal_size_strategy": "sealed artifacts retain full runtime/source/SDK maps, scientific objects and169checks; compact terminal binds artifact identities; oversize/serialization refuses via distinct closed64KiBfallback before creating record once",
    "earned_payload_rule": "retain a completed scientific payload if only later identity/transport fails; statusrefused and comparisonacceptance unavailable; null only never-earned fullpayload with retained partialprefixes",
    "seal": "exclusive new terminal record; regular output files0444/executable0555, no prior receipt edits"
  },
  "qualification": {
    "execution": "unperformed",
    "numerical": "new committed attempt and original169 comparisons required",
    "inference": "blocked",
    "interpretation": "synthetic-conditional-only",
    "physical_qualification": false,
    "observational_qualification": false,
    "unmet_gates": "all original contract unmet_gates retained; no measured template/calibration/rawdetector/population/covariance qualification"
  },
  "freeze_sequence": {
    "before_implementation": [
      "parent accepts embedded closed output/runtime/interface/model/domain/allocation contract",
      "parent accepts layout/path/resource/direct-child ledger and loader policy"
    ],
    "after_code_before_execution": [
      "I/O-only reference fork and dedicated controller/verifier source completed and deliberately reviewed",
      "all null port hashes replaced",
      "final executable request hash bound by experiment metadata",
      "clean committed checkpoint",
      "parent explicit single-job grant before tests/compile/science"
    ]
  },
  "shape_contract": {
    "schema_version": 1,
    "id": "temporal-shared-optical-closed-shape-proposal/v2",
    "status": "source-only proposal; requires parent contract/request approval before implementation",
    "scope": "transport/type/domain/order proposal only; original scientific arithmetic and budgets unchanged",
    "provenance": {
      "admission_sha256": "1ee77abd937be8342f3be946c17453a162f75e9cd315a9a53b05f8fe2db401cc",
      "original_contract_sha256": "e3a1ceb57d96f599cf82d938462d03feef5d424ad9d2c70f82316d63c34ba6d7",
      "original_consumer_sha256": "46e17d76d0b03390048fdc66ba97ca7012d8eda215beb9ab2aeea96a92905ab8",
      "original_reference_sha256": "f912a6f3b0b351b91ece9d26af959200618776908d1842bb64caa63ed87ebe0e",
      "original_controller_sha256": "4f76eca4f9ec4f01ead57978a7031f6dbf28e07eaf14c7e7e286424a1412e535",
      "original_runtime_subset_map_sha256": "982dc09f2ddf1206cc9cd4467c7be1ce4238a612489b43bf5f5d196a8b2495a8",
      "derivation": "Actual preserved consumer emitters and reference return/exception definitions; original sealed records used only for shape/identity/failure semantics. No computed acceptance values reused as constants or expected answers.",
      "previous_proposal_sha256": "30f3ba5add7849a1ebc2d6d223a0ddbc001258ad0b8bdd8be645a40e39a9df36"
    },
    "notation": {
      "object_rule": "Every named object below is closed; every listed field is required unless explicitly optional. Nullable means exactly JSON null or the named type. Unknown, duplicate or missing keys refuse transport admission.",
      "array_rule": "array(T,min,max) is ordered, bounded and element-typed. Complete routes require exact stated axis order; partial arrays are prefixes, never silently compressed or renumbered.",
      "bool": "type(value) is bool; never accept integers0/1",
      "uint": "type(value) is int and not bool, nonnegative, bounded as stated; integral Decimal or floating tokens do not become integer identity/work fields",
      "number": "finite JSON integer or decimal number token, decoded directly to int/Decimal; bool and decoded float objects refused. Preserve exact emitted decimal digits, including21-digit native longdouble serialization; no binary64 cast for comparisons.",
      "nonnegative": "number>=0, finite",
      "probability": "number in[0,1], finite",
      "decimal_string": "ASCII finite decimal string length<=256; grammar -?(0|[1-9][0-9]*)(\\.[0-9]+)?([eE][+-]?[0-9]+)?; adjusted exponent absolute value<=1000; parse directly as Decimal without float intermediary",
      "decimal_error": "decimal_string>=0; positivity/range relations checked by context",
      "sha256": "exact64 lowercase hexadecimal ASCII characters",
      "path": "absolute path string, UTF8<=4096bytes, no NUL; input/output/source admission refuses symlinks/escaping paths. Pinned readonly runtime requested aliases may be symlinks only where explicitly identified as requested paths; record and check their exact resolved target. Inventory filepaths are resolved regular files.",
      "status": [
        "ok",
        "invalid_input",
        "nonfinite_input",
        "overflow",
        "work_limit",
        "outside_domain",
        "singular",
        "not_positive_definite",
        "conditioning_budget_exceeded"
      ],
      "unknown_status": "The original emitter fallback 'unknown' is a transport refusal; raw bytes and failure diagnostics are retained.",
      "bounds": {
        "json_bytes": 1048576,
        "stderr_bytes": 1048576,
        "stdout_bytes": 1048576,
        "maximum_depth": 32,
        "maximum_total_items": 65536,
        "maximum_string_utf8_bytes": 8192,
        "maximum_numeric_token_characters": 256,
        "identifier_ascii_bytes": 256,
        "message_utf8_bytes": 4096
      }
    },
    "axes": {
      "synthetic_identity": "synthetic-fixed-bilinear-source-two-exposures-two-bands-one-common-two-state-optical-law",
      "grid": {
        "index": 0,
        "id": "synthetic-bilinear-ux16-v1",
        "source_origin": "original analytic synthetic control",
        "requested_mask": 4,
        "requested_group": "transmitted_photons",
        "omitted_groups": [
          "mean_flux",
          "energy"
        ]
      },
      "states": [
        {
          "index": 0,
          "id": "S0",
          "relative_mass": 1,
          "normalized_mass": "1/4"
        },
        {
          "index": 1,
          "id": "S1",
          "relative_mass": 3,
          "normalized_mass": "3/4"
        }
      ],
      "physical_bands": [
        {
          "index": 0,
          "id": "blue",
          "observed_wavelength_knots_metre": [
            2,
            3
          ]
        },
        {
          "index": 1,
          "id": "red",
          "observed_wavelength_knots_metre": [
            3,
            4
          ]
        }
      ],
      "prepared_band_order": [
        "S0/blue",
        "S0/red",
        "S1/blue",
        "S1/red"
      ],
      "exposures": [
        {
          "index": 0,
          "id": "E0",
          "observer_interval_second": [
            9,
            11
          ],
          "full_duration_second": 2,
          "covered_duration_second": 2,
          "expected_coverage": "full",
          "covered_fraction": "1"
        },
        {
          "index": 1,
          "id": "E1",
          "observer_interval_second": [
            12,
            15
          ],
          "full_duration_second": 3,
          "covered_duration_second": 2,
          "expected_coverage": "partial",
          "covered_fraction": "2/3"
        }
      ],
      "channels": [
        "E0/blue",
        "E0/red",
        "E1/blue",
        "E1/red"
      ],
      "primary_rows": [
        "S0/E0/blue",
        "S0/E0/red",
        "S0/E1/blue",
        "S0/E1/red",
        "S1/E0/blue",
        "S1/E0/red",
        "S1/E1/blue",
        "S1/E1/red"
      ],
      "row_mapping": "i=4*state+2*exposure+band; native band_index=2*state+band; detector state_index=i//4 and channel_index=i%4",
      "sigma_controls_electrons": [
        0,
        0.75
      ],
      "reference_record_order": [
        "discrete-all-detected",
        "discrete-selected-all-detected",
        "discrete-censored",
        "continuous-all-detected",
        "continuous-selected-all-detected",
        "continuous-censored"
      ],
      "per_sigma_record_indices": [
        0,
        1,
        2
      ],
      "detector_batch_slots": [
        "nominal",
        "lower-photon-sensitivity",
        "upper-photon-sensitivity"
      ],
      "detector_row_slots": [
        "joint-all-detected",
        "selected-control-joint-numerator",
        "joint-censored",
        "threshold-nondetection-probe0",
        "threshold-nondetection-probe1",
        "threshold-nondetection-probe2"
      ],
      "temporal_control_ids": [
        "no-time-overlap",
        "zero-duration",
        "invalid-redshift",
        "required-state-row-refusal",
        "temporal-work-refusal",
        "temporal-row-admission-refusal"
      ],
      "detector_control_ids": [
        "detector-QE-domain",
        "detector-work-refusal",
        "selected-nondetection",
        "zero-noise-threshold-boundary"
      ],
      "reducer_control_ids": [
        "required-QE-refusal",
        "required-work-refusal",
        "selected-censored-refusal"
      ],
      "record_measure": {
        "sigma0_record0": "joint discrete counting probability",
        "sigma0_record1": "selected-only discrete counting probability given same ALL-four-detected event",
        "sigma0_record2": "joint three detected counts times one nondetection probability",
        "sigma075_record0": "joint numerical ADU^-4 density",
        "sigma075_record1": "selected-only numerical ADU^-4 density given same ALL-four-detected event",
        "sigma075_record2": "joint numerical ADU^-3 density times one nondetection probability"
      },
      "selection_event": "ALL four channels detected, Y>=1.5ADU; zero-noise K>=1; one shared-state mixture denominator",
      "unit_policy": "photon expectation=count, lambda=expected electrons, durations=observer seconds, continuous density=caller ADU coordinate; reference log/difference measures follow record map. Native payload lacks these labels: immutable request/controller record binds them without modifying consumer bytes.",
      "alternative_records": "Never multiply record controls together as independent measurements; conditionally independent readouts share one optical state across both exposures.",
      "preservation": "State IDs/masses, rows and requested/omitted masks retained on failures; no dropping, reordering or renormalizing."
    },
    "native_types": {
      "Outcome": {
        "availability": "uint enum0=omitted,1=available,2=failed",
        "status": "status",
        "value": "nullable nonnegative number"
      },
      "TemporalRow": {
        "index": "uint[0,7]",
        "grid_index": "uint[0,0]",
        "band_index": "uint[0,4];4 only declared required-state-row-refusal",
        "source_epoch_second": "number exactly10",
        "observer_lower_second": "number",
        "observer_upper_second": "number",
        "admission_status": "status",
        "coverage": "enum unassessed,no_overlap,partial,full",
        "observer_duration_second": "nonnegative",
        "covered_observer_second": "nonnegative",
        "covered_fraction": "probability",
        "segment_work": "uint[0,512]",
        "photons": "Outcome",
        "mean_flux": "Outcome",
        "energy": "Outcome"
      },
      "TemporalBatch": {
        "status": "status",
        "segment_work": "uint[0,512]",
        "required_rows_admitted": "bool",
        "rows": "array(TemporalRow,0,8)"
      },
      "DetectorRow": {
        "status": "status",
        "detected": "bool",
        "threshold_adu": "number exactly1.5",
        "electron_count": "nullable uint[0,4294967295]",
        "measured_adu": "nullable number",
        "zero_probability": "bool",
        "log_value": "nullable number",
        "detection_probability": "nullable probability",
        "log_error": "nonnegative"
      },
      "DetectorBatch": {
        "status": "status",
        "photons": "nonnegative",
        "QE": "number; allow declared refusedQE1.01, never impose[0,1] on a retained failed-control source",
        "full_exposure_second": "number exactly2 or3 by channel/control",
        "sigma_electrons": "number exactly0 or0.75 by sigma/control",
        "poisson_terms": "uint[0,256]",
        "omitted_tail": "probability",
        "rows": "array(DetectorRow,0,6); primary successfulbatch6, single-detector controls0..1"
      },
      "DetectorAttempt": {
        "state_index": "uint[0,1]",
        "channel_index": "uint[0,3]",
        "lambda_electrons": "nullable nonnegative",
        "batches": "array(nullable DetectorBatch,3,3)"
      },
      "ConditionalState": {
        "state_id": "enum S0,S1 by slot",
        "mass": "number exactly0.25 or0.75 by slot",
        "status": "status",
        "log_record": "number",
        "log_event": "number",
        "record_error": "nonnegative",
        "event_error": "nonnegative"
      },
      "ReducerRecord": {
        "record_index": "uint[0,2]",
        "conditional": "array(ConditionalState,2,2)",
        "selected_only": "bool",
        "status": "status",
        "log_joint": "nullable number",
        "log_event": "nullable number",
        "log_value": "nullable number",
        "joint_log_error": "nonnegative",
        "event_log_error": "nonnegative",
        "value_log_error": "nonnegative",
        "product_per_exposure_log_record": "nullable number",
        "product_per_exposure_log_event": "nullable number",
        "record_relative_difference": "nullable number",
        "event_relative_difference": "nullable number"
      },
      "DetectorControl": {
        "sigma_electrons": "number exactly0 or0.75 by slot",
        "attempts": "array(DetectorAttempt,8,8)",
        "records": "array(ReducerRecord,3,3)"
      },
      "TemporalControl": {
        "id": "exacttemporalcontrolID by slot",
        "joint_withheld": "bool",
        "native": "TemporalBatch"
      },
      "SingleDetectorControl": {
        "id": "exactdetectorcontrolID by slot",
        "native": "DetectorBatch"
      },
      "ExternalReducerControl": {
        "id": "exactreducercontrolID by slot",
        "retained_state_masses": "array(number,2,2) exactly[0.25,0.75]",
        "replacement_batch": "nullable DetectorBatch",
        "reducer": "ReducerRecord"
      },
      "NativePayload": {
        "schema": "exact external-temporal-optical-sdk-control/v1",
        "model_id": "exact finite_bilinear_rest_spectral_time_zero_outside_full_observer_exposure_mean",
        "detector_model_id": "exact DETECTOR/Poisson-arrivals-fixed-QE-Gaussian-read/v1",
        "selection_id": "exact DETECTOR/threshold-joint-or-selected-censoring/v1",
        "sampled_empirical_sensitivity": "number; source-bound longdouble literal3e-12, exact emitted token preserved",
        "temporal_photon_endpoint_sensitivity": "number; source-bound longdouble literal3.2e-12, exact emitted token preserved",
        "prepared_status": "status",
        "grid_count": "uint[0,1]",
        "band_count": "uint[0,4]",
        "retained_temporal_payload_bytes": "nullable uint[0,1048576]",
        "prepare_seconds": "nonnegative",
        "temporal_evaluate_seconds": "nonnegative",
        "primary": "TemporalBatch",
        "detector_controls": "array(DetectorControl,0,2)",
        "primary_poisson_terms": "uint[0,12288]",
        "actual_temporal_controls": "array(TemporalControl,6,6)",
        "actual_detector_controls": "array(SingleDetectorControl,4,4)",
        "actual_external_reducer_controls": "array(ExternalReducerControl,0,3)",
        "all_poisson_terms": "uint[0,13824]",
        "all_temporal_segment_work": "uint[0,3584]",
        "detector_six_row_payload_bound_bytes": "nullable uint[0,1048576]",
        "total_native_seconds": "nonnegative"
      }
    },
    "native_semantic_rules": [
      "Primary required_rows_admitted exactly recomputes existing emitter rule: batch statusok, eight rows, each admissionok and photons availability1/statusok/finite nonnegative. This flag alone never establishes numerical acceptance.",
      "Outcome available1 requires statusok and finite nonnegative value; omitted0 and failed2 require null value. Requested mean_flux and energy remain omitted0/statusok/value=null on every row because mask4 requests photons only.",
      "Row identity and exposure bounds are exact from source/request; failed temporal rows may retain unassessed coverage and zero defaults. Do not apply successfulduration/coverage relations to withheld failed values.",
      "Native detector_controls has exactly2 groups and external_reducer_controls exactly3 iff primary is admitted; otherwise both are empty. Temporal6 and single-detector4 control families are emitted even when primary fails.",
      "Detector batch slot0 is mandatory; endpoint slots1/2 null only for exact nominal photons0. Positive nominal expectation requires both complete endpoint attempts, including their actual outside_domain/refusal records. Never infer skip from a missing/failed row.",
      "Nominal/endpoint six-row source: rows0/1 detected; row2 nondetected only channel3; rows3/4/5 threshold-only nondetection. sigma0 detected rows hold electron_count and nullmeasured_adu; sigma075 detected rows hold measured_adu and nullcount; nondetection holds bothnull.",
      "Detector failed status may coexist with zero rows or failed individual rows under batchstatusok. Structuralzero is statusok/zero_probabilitytrue/log_valuenull; it is not a numerical refusal. Failedstatus finite-log availability follows actual SDK; required reducers cannot consume a failedrow.",
      "Reducer statusok requires all seven aggregate/comparator scalar values finite; refusedstatus requires all sevennull. Conditional state logs/errors remain emitted even on refusal; status binds their validity and zero placeholders never become earned complete conditional data.",
      "Selected record1 divides its joint numerator by the same ALL-four-detected event exactly once. Selected censored record2 refusesinvalid_input and withholds aggregate/comparator groups. Replacement batch isnull only selected-censored-refusal; requiredQE/work controls retain the actual failedreplacement batch.",
      "All work sums recompute over actually retained nonnull batches/controls, without using historical measured work counts as expected constants. EachbatchPoissonterms<=256; primary48batches bound12288, all54batches bound13824; seven temporal calls each<=512 yield3584 totalcap.",
      "Exactly representable source values such as0.25/0.75/1.5 and integerduration retain exact numeric equality. QE source literal1.01 is binary64 and its full emitted decimal token is preserved; admit it only in declared failed-control contexts instead of asserting exactdecimal1.01. Longdouble sensitivity literal source identities are bound by exact unchanged compiled source and SDK headers, never rounded to a binary64 parser target.",
      "The unchanged C++ emitter can print a non-JSON inf/nan for an unexpected nonfinite diagnostic. Strict parsing must refuse such output and preserve complete bounded raw bytes/execution failure; do not sanitize it into a successful typed payload or alter the original consumer."
    ],
    "reference_types": {
      "FrequencyRow": {
        "polynomial": "decimal_string",
        "polynomial_error": "decimal_error",
        "frequency32": "decimal_string",
        "frequency64": "decimal_string",
        "frequency32_error": "decimal_error",
        "frequency64_error": "decimal_error"
      },
      "DiscreteDetectorReferenceRow": {
        "count_probability": "decimal_string",
        "count_probability_error": "decimal_error",
        "below": "decimal_string",
        "below_error": "decimal_error"
      },
      "ContinuousDetectorReferenceRow": {
        "density32_per_adu": "decimal_string",
        "density64_per_adu": "decimal_string",
        "density_error_per_adu": "decimal_error",
        "below32": "decimal_string",
        "below64": "decimal_string",
        "below_error": "decimal_error",
        "density_tail_per_adu": "decimal_error",
        "cdf_tail": "decimal_error"
      },
      "DetectorReferenceGroup": {
        "sigma_electrons": "number exactly0 or0.75 by slot",
        "rows": "array(DiscreteDetectorReferenceRow orContinuousDetectorReferenceRow accordingtosigma,8,8)"
      },
      "ReferenceRecord": {
        "sigma_index": "uint[0,1]",
        "record_index": "uint[0,2]",
        "conditional_record": "array(decimal_string,2,2)",
        "conditional_event": "array(decimal_string,2,2)",
        "log_joint": "decimal_string",
        "joint_log_error": "decimal_error",
        "log_event": "decimal_string",
        "event_log_error": "decimal_error",
        "log_value": "decimal_string",
        "value_log_error": "decimal_error",
        "product_per_exposure_log_record": "decimal_string",
        "product_per_exposure_log_event": "decimal_string",
        "product_per_exposure_record_log_error": "decimal_error",
        "product_per_exposure_event_log_error": "decimal_error",
        "record_relative_difference": "decimal_string",
        "event_relative_difference": "decimal_string",
        "record_signed_difference": "decimal_string",
        "record_signed_difference_error": "decimal_error",
        "event_signed_difference": "decimal_string",
        "event_signed_difference_error": "decimal_error",
        "record_relative_difference_error": "decimal_error",
        "event_relative_difference_error": "decimal_error",
        "positive_joint": "decimal_string",
        "positive_joint_error": "decimal_error",
        "positive_event": "decimal_string",
        "positive_event_error": "decimal_error"
      },
      "ReferencePayload": {
        "schema": "exact original-high-precision-temporal-optical-reference/v1",
        "mpmath_version": "exact1.3.0",
        "decimal_digits": "uint exactly60",
        "shared_ancestry": "array(ASCIIidentifier,2,2) exactly['declared equations and exact SI constants','mpmath arithmetic/GL node generator; no native kernels']",
        "frequency_seconds": "nonnegative",
        "elapsed_seconds": "nonnegative",
        "integrand_evaluations": "closed{time_frequency:uint[0,30720],characteristic_function:uint[0,73728]}; complete route exactly30720/73728 from original loops",
        "photons": "array(decimal_string,8,8)",
        "photon_errors": "array(decimal_error,8,8)",
        "frequency_rows": "array(FrequencyRow,8,8)",
        "lambda_electrons": "array(decimal_string,8,8)",
        "detector_rows": "array(DetectorReferenceGroup,2,2)",
        "records": "array(ReferenceRecord,6,6)"
      },
      "PositiveMixtureFailure": {
        "stage": "exact positive-mixture-error-gate",
        "value": "decimal_string",
        "error": "decimal_error",
        "conditional_values": "array(decimal_string,2,2)",
        "conditional_errors": "array(decimal_error,2,2)"
      },
      "ConditionalPositiveFailure": {
        "stage": "exact conditional-positive-reference-gate",
        "sigma_index": "uint[0,1]",
        "channel_index": "uint[0,7]",
        "lambda": "decimal_string",
        "density": "decimal_string",
        "density_error": "decimal_error",
        "nondetection": "decimal_string",
        "nondetection_error": "decimal_error",
        "coarse_refined_attempt": "DiscreteDetectorReferenceRow orContinuousDetectorReferenceRow accordingtosigma"
      },
      "TransportFailureDetails": {
        "operation": "boundedASCIIidentifier",
        "path": "nullable path",
        "reason": "boundedUTF8message"
      },
      "ReferenceError": {
        "kind": "enum admission,runtime,reference,io,resource",
        "stage": "boundedASCIIidentifier",
        "message": "boundedUTF8message",
        "message_truncated": "bool",
        "message_sha256": "sha256 of originaluntruncated UTF8message",
        "details": "nullable PositiveMixtureFailure orConditionalPositiveFailure orTransportFailureDetails"
      },
      "ReferencePartial": {
        "completed_groups": "array(enum frequency,photons,lambda,detector-sigma0,detector-sigma075,records,0,6); original group completion order",
        "photon_prefix": "array(decimal_string,0,8)",
        "photon_error_prefix": "array(decimal_error,0,8)",
        "frequency_prefix": "array(FrequencyRow,0,8)",
        "lambda_prefix": "array(decimal_string,0,8)",
        "detector_prefix": "array(DetectorReferenceGroup,0,2)",
        "active_sigma_index": "nullable uint[0,1]",
        "active_sigma_rows": "array(DiscreteDetectorReferenceRow orContinuousDetectorReferenceRow accordingtoactive_sigma_index,0,8)",
        "record_prefix": "array(ReferenceRecord,0,6)",
        "control_failure": "nullable PositiveMixtureFailure orConditionalPositiveFailure"
      },
      "ReferenceIdentities": {
        "request_sha256": "sha256",
        "original_contract_sha256": "sha256 exact originalv2",
        "admission_sha256": "sha256 exactpreserved admission",
        "original_reference_sha256": "sha256 exactf912 ancestry",
        "reference_script_sha256": "sha256 actualnewforksource",
        "current_source_sha256": "closed map of consumer.cpp,reference.py,controller.py,verify_sdk.py,contract.snapshot.json tosha256; actual admitted snapshot bytes"
      },
      "ReferenceWork": {
        "elapsed_seconds": "nonnegative",
        "time_frequency_evaluations": "uint[0,30720]",
        "characteristic_function_evaluations": "uint[0,73728]",
        "cpu_user_seconds": "nonnegative",
        "cpu_system_seconds": "nonnegative",
        "maximum_rss_kib": "uint[0,4294967295]"
      },
      "ReferenceQualification": {
        "scope": "exact empirical-synthetic-fixed-temporal-shared-optical-control/v1",
        "native_outputs_consumed": "bool exactlyfalse",
        "physical_qualification": "bool exactlyfalse",
        "observational_qualification": "bool exactlyfalse",
        "original_input_certificate": "bool exactlyfalse",
        "limits": "array(boundedUTF8message,1,16)"
      },
      "ReferenceEnvelope": {
        "schema_version": "uint exactly1",
        "interface_id": "exact temporal-shared-optical-reference-transport/v1",
        "status": "enum completed,refused",
        "error": "nullable ReferenceError",
        "identities": "nullable ReferenceIdentities untiladmissioncompleted",
        "axes": "exact immutable axes snapshot excluding documentaryunitnotes",
        "settings": "exact approvedrequest scientific/resource settings snapshot",
        "runtime_before": "nullable RuntimeFingerprint",
        "runtime_after": "nullable RuntimeFingerprint",
        "payload": "nullable ReferencePayload",
        "partial": "ReferencePartial",
        "work": "ReferenceWork",
        "qualification": "ReferenceQualification"
      }
    },
    "reference_semantic_rules": [
      "Complete payload/scientific data field names/order and arithmetic remain originalsource f912; the new transport envelope never supplies new scientific arguments or fitted values.",
      "Completed means scientific route emitted and source/runtime/transport identities unchanged; it is not controller comparison acceptance. Controller retains all original169 comparison semantics and original earned-error charging.",
      "Completed requires payloadnonnull, errornull, all full axes/array orders, runtimebefore/afternonnull equal, actual source/request identities equal. Refused preserves earned prefixes/details and actual work; complete payload may remain present if only after-seal identity fails but is unavailable for numerical acceptance.",
      "Partial photon/error/frequency prefixes have equal lengths. Full completed groups remain exact fullaxis; active_sigma_rows are an ordered prefix and control_failure retains failing row coarse/refined values. Records are earned in the original interleaved order: sigma0 rows then records0..2, sigma075 rows then records0..2. The records completed-group marker requires all six records, not the first three. Do not fabricate missing fine values, zeros or completed groups.",
      "All highprecision numbers stay60-digit scalarstrings; signed comparator differences may be negative, continuous density may exceed1, and errors are nonnegative. Complete positive densities/probabilities and mixtures must satisfy the original positive/error gate; failure details allow nonpositive failed values to remain diagnostic.",
      "Event probability is dimensionless and in(0,1]; continuous likelihood is ADU density, not constrained<=1. Positiveeventerror<positiveevent and positivejointerror<positivejoint are required. Selected value is same logjoint-logevent only record1.",
      "No independent event-dependence lower threshold is invented. Original numerator witness1e-3 and every event/comparator diagnostic retain separate roles and errors.",
      "Request SHA is actual admitted bytes echoed and independently matched by controller. Do not hardcode new requestSHA inside a reference that the request itself hashes: avoid cyclic source/request identity. Frozen interface/domain/budgets precede implementation; exact implemented script hashes and executable request follow source review. The approved closed axes/settings and shape/domain definitions are self-contained in public request bytes; ignored proposal hashes are documentary ancestry only.",
      "On killed/timeout process or unparseable/nonfinite exception payload, retain bounded raw output/partialfile and execution failure; absence of a typed final envelope cannot be relabeled as successful reference.",
      "No original acceptance values, historical fitted output, native photons/logs, comparison receipts or immutable attempts are consumed by reference arithmetic. Sealed old records are documentary lineage only."
    ],
    "runtime_types": {
      "FileIdentity": {
        "path": "path absolute resolved",
        "bytes": "uint[0,67108864]",
        "sha256": "sha256"
      },
      "PythonIdentity": {
        "version": "boundedUTF8message exactpinned3.14.8seal",
        "implementation": "exactCPython",
        "requested_executable": "path",
        "resolved_executable": "path",
        "executable_bytes": "uint[0,67108864]",
        "executable_sha256": "sha256 exact815b1275bf87e7595fe1cdf3f2efc6ae58e27b842262a7b530de8610ac2760ba",
        "prefix": "path",
        "base_prefix": "path",
        "pyvenv_cfg": "FileIdentity exactcfg11a0241468a7d4cb44419cd7c7a87cd6149c3e90e38b75be6e1a201279621e43",
        "dont_write_bytecode": "bool exactlytrue"
      },
      "MpmathIdentity": {
        "version": "exact1.3.0",
        "module_origin": "path",
        "module_root": "path",
        "backend": "exactpython",
        "decimal_digits": "uint exactly60",
        "binary_precision_bits": "uint actualmp context atdps60",
        "gmpy2_loaded": "bool exactlyfalse",
        "original_subset_map_sha256": "sha256 exact982dc09f2ddf1206cc9cd4467c7be1ce4238a612489b43bf5f5d196a8b2495a8"
      },
      "ImportedModuleIdentity": {
        "name": "boundedASCIIidentifier exactsys.moduleskey",
        "kind": "enum built-in,frozen,source,extension,namespace",
        "origin": "nullable path; source/extension require resolvedregular origin ininventory, built-in/frozen/namespace require null",
        "namespace_paths": "array(path,0,32); nonemptyonlynamespace, sortedunique"
      },
      "RuntimeFingerprint": {
        "schema_version": "uint exactly1",
        "interface_id": "exact temporal-shared-optical-reference-runtime/v1",
        "python": "PythonIdentity",
        "mpmath": "MpmathIdentity",
        "thread_environment": "closed OMP_NUM_THREADS,OPENBLAS_NUM_THREADS,MKL_NUM_THREADS,NUMEXPR_NUM_THREADS,BLIS_NUM_THREADS,VECLIB_MAXIMUM_THREADS,PYTHONDONTWRITEBYTECODE all exactstring'1'",
        "actual_thread_count": "uint exactly1 from/proc/self/status",
        "imported_modules": "array(ImportedModuleIdentity,1,4096) sorteduniquename; complete actualsys.modules set after metadata/helperimports",
        "mapped_libraries": "array(FileIdentity,0,256) sorteduniquepath; complete actual regular native .so mappings subset ofinventory",
        "inventory": "array(FileIdentity,1,4096) sorteduniquepath; totalbytes<=67108864",
        "inventory_sha256": "sha256 canonicalcompactUTF8 list sortedobjectkeys",
        "numeric_configuration": "closed{arithmetic_backend:exactmpmath-libmp-python,decimal_digits:uint60,binary_precision_bits:uintactual,active_native_numeric_backends:emptyarray,gmpy2_loaded:false,byteorder:enumlittle,big}"
      }
    },
    "runtime_semantic_rules": [
      "Seal original complete selected mpmath*.py/license/dist-info metadata subset map against its preserved canonical hash, not merely __init__.py/version. Add every relevant packagefile, importedstdlib source/extension, executable/cfg, currentreference source and mappedregular.so to expanded complete inventory.",
      "Record actual moduleorigin, BACKENDpython and actualprocessThreads1. Environment strings state policy only. Refuse a loaded gmpy2 or native numerical backend rather than silently changing arithmetic ancestry.",
      "Finish metadata/backend discovery before inventory. All mapped regular.so files from/proc/self/maps are independently hashed; no claim that native loaderresolution proves an actual processmapping.",
      "Fingerprint before/after source/science on success and failure. Runtime getter calls noGL nodegeneration, source integrals, Fourier evaluation or likelihood mixture; no package installation, SDK/library/runtime edits or bytecode writes.",
      "Expanded inventory canonical digest is SHA256(json.dumps(inventory,sort_keys=True,separators=(',',':'),allow_nan=False).encode('utf8')); paths are sorted unique resolved absolute strings. Every fileidentity hash is independently checked by controller.",
      "Separate complete imported_modules and mapped_libraries sets must also be exactly equal before/after and externalpreflight/terminal. Unionfileinventory equality alone cannot detect a new imported mpmath module whose file was already in the complete package inventory. Unsupported moduleorigin, disappeared/deleted mapping or unreadable library is an explicit runtime refusal, not an omitted identity.",
      "The controller records pinned native ELF loader --list output before native and in terminal finally as structural dependency resolution only. This does not establish actual native /proc mapping. Reference fingerprint derives its own actual mapped-library set directly from its current process."
    ],
    "unchanged_scientific_settings": {
      "decimal_digits": 60,
      "arithmetic_reservation": "1e-50*(1+abs(quantity))",
      "frequency_and_Fourier_GL_orders": [
        32,
        64
      ],
      "Fourier_panel_width_inverse_electron": "1/4",
      "Fourier_panels": 96,
      "Fourier_cutoff_inverse_electron": 24,
      "sigma_positive_electrons": "3/4",
      "state_masses": [
        "1/4",
        "3/4"
      ],
      "endpoint_relative_empirical_allocation": "3.2e-12",
      "native_photon_comparison": "1e-300+2e-12*abs(polynomial)",
      "frequency_comparison_and_refinement": "1e-300+2e-13*abs(polynomial)",
      "joint_event_selected_log_comparison": "2e-10*(1+abs(reference_log))",
      "reference_log_refinement_max_fraction": "0.05",
      "native_external_aggregate_log_diagnostic": "5e-12+1e-8*(1+abs(native_log))",
      "same_selection_identity": "1e-17 controllercheck",
      "numerator_dependence_witness": "1e-3, syntheticcriterion only",
      "event_dependence_minimum": null
    },
    "proposed_resource_limits": {
      "jobs": 1,
      "threads": 1,
      "compile_cpu_seconds": 180,
      "compile_wall_seconds": 180,
      "native_cpu_seconds": 180,
      "native_wall_seconds": 180,
      "reference_cpu_seconds": 180,
      "reference_wall_seconds": 180,
      "runtime_preflight_cpu_seconds": 180,
      "runtime_preflight_wall_seconds": 180,
      "address_limit_bytes": 1073741824,
      "stdout_limit_bytes_each": 1048576,
      "stderr_limit_bytes_each": 1048576,
      "reference_json_limit_bytes": 1048576,
      "source_tree_bytes": 2097152,
      "source_archive_bytes": 4194304,
      "source_file_bytes": 1048576,
      "attempt_store_bytes": 268435456
    },
    "gates": {
      "current": "proposal only; parent must approve closed contract and immutable request before any implementation",
      "compute": "zero allocated jobs; no imports, syntax checks, tests, builds, native/reference execution or scientific recomputation",
      "preserve": "All48 admitted files/original failedattempts/sealedreceipts remain byte-identical; public/SDK/runtime/package changes forbidden in this subtask",
      "interpretation": "Synthetic numerical control only; physical and observational source/calibration/detector/population/selection/covariance/NEXT17 gates remain open"
    }
  }
}

'''


def native_admission(value, request):
    types = request['shape_contract']['native_types']; axes = request['shape_contract']['axes']
    def obj(item, name): return closed(item, types[name], name)
    def outcome(item):
        obj(item, 'Outcome'); uint(item['availability'], 2); status(item['status'])
        if item['availability'] == 1:
            exact(item['status'], 'ok'); number(item['value'], 0)
        elif item['value'] is not None:
            raise ValueError('unavailable outcome has a value')
        if item['availability'] == 0: exact(item['status'], 'ok')
        if item['availability'] == 2 and item['status'] == 'ok': raise ValueError('failed outcome statusok')
    def temporal(item, control=None):
        obj(item, 'TemporalBatch'); status(item['status']); uint(item['segment_work'], 512)
        boolean(item['required_rows_admitted']); array(item['rows'], maximum=8)
        for i, row in enumerate(item['rows']):
            obj(row, 'TemporalRow'); exact(row['index'], i); exact(row['grid_index'], 0)
            wanted_band = 2*(i//4)+(i%2)
            if control == 'required-state-row-refusal' and i == 7: wanted_band = 4
            exact(row['band_index'], wanted_band); number(row['source_epoch_second']); exact(number(row['source_epoch_second']), Decimal(10))
            e = (i%4)//2; lo, hi = ((9,11),(12,15))[e]
            if control == 'no-time-overlap' and i == 0: lo, hi = 20,22
            if control == 'zero-duration' and i == 0: hi = lo
            exact(number(row['observer_lower_second']), Decimal(lo)); exact(number(row['observer_upper_second']), Decimal(hi))
            status(row['admission_status']); uint(row['segment_work'], 512)
            if row['coverage'] not in ('unassessed','no_overlap','partial','full'): raise ValueError('coverage enum')
            for key in ('observer_duration_second','covered_observer_second'): number(row[key],0)
            number(row['covered_fraction'],0,1)
            for key in ('photons','mean_flux','energy'): outcome(row[key])
            for key in ('mean_flux','energy'): exact(row[key]['availability'],0)
            if row['admission_status'] == 'ok':
                if row['covered_observer_second'] > row['observer_duration_second']: raise ValueError('coverage duration')
            elif row['photons']['availability'] == 1: raise ValueError('failed row has available photons')
        admitted = item['status']=='ok' and len(item['rows'])==8 and all(
            row['admission_status']=='ok' and row['photons']['availability']==1 and row['photons']['status']=='ok'
            and row['photons']['value'] is not None for row in item['rows'])
        exact(item['required_rows_admitted'], admitted)
        exact(item['segment_work'], sum(row['segment_work'] for row in item['rows']))
        return admitted
    def batch(item, sigma=None, channel=None, single=False, bad_qe=False):
        obj(item,'DetectorBatch'); status(item['status']); number(item['photons'],0)
        number(item['QE'],0,Decimal('1.1'))
        if bad_qe:
            if not Decimal(1) < number(item['QE']) < Decimal('1.1'): raise ValueError('declared badQE control')
        elif channel is not None: exact(number(item['QE']), Decimal('0.5' if channel%2 else '0.75'))
        else: exact(number(item['QE']),Decimal('0.75'))
        if channel is not None: exact(number(item['full_exposure_second']),Decimal(2 if channel<2 else 3))
        else: exact(number(item['full_exposure_second']),Decimal(2))
        if sigma is not None: exact(number(item['sigma_electrons']),Decimal(sigma))
        else: number(item['sigma_electrons'],0,Decimal('0.75'))
        uint(item['poisson_terms'],256); number(item['omitted_tail'],0,1)
        array(item['rows'],maximum=1 if single else 6)
        for i,row in enumerate(item['rows']):
            obj(row,'DetectorRow'); status(row['status']); boolean(row['detected']); boolean(row['zero_probability'])
            exact(number(row['threshold_adu']),Decimal('1.5'))
            if row['electron_count'] is not None: uint(row['electron_count'])
            if row['measured_adu'] is not None: number(row['measured_adu'])
            for key in ('log_value','detection_probability'):
                if row[key] is not None: number(row[key],0 if key=='detection_probability' else None,1 if key=='detection_probability' else None)
            number(row['log_error'],0)
            if not single:
                wanted = i < 3 and not (i==2 and channel==3)
                exact(row['detected'],wanted)
                if wanted:
                    if sigma==0:
                        exact(row['electron_count'],(2,3,3,2)[channel]); exact(row['measured_adu'],None)
                    else:
                        exact(number(row['measured_adu']),Decimal(3 if channel in (0,3) else '3.5')); exact(row['electron_count'],None)
                else: exact(row['electron_count'],None); exact(row['measured_adu'],None)
            if row['zero_probability'] and (row['status']!='ok' or row['log_value'] is not None): raise ValueError('structuralzero contradiction')
            if row['status']!='ok' and row['log_value'] is not None: raise ValueError('failed row finite log')
            if row['status']=='ok' and not row['zero_probability'] and row['log_value'] is None: raise ValueError('successful detector log unavailable')
        return item['poisson_terms']
    def reducer(item):
        obj(item,'ReducerRecord'); uint(item['record_index'],2); boolean(item['selected_only']); status(item['status'])
        array(item['conditional'],2)
        for i,row in enumerate(item['conditional']):
            obj(row,'ConditionalState'); exact(row['state_id'],'S'+str(i)); exact(number(row['mass']),Decimal('0.25' if i==0 else '0.75'))
            status(row['status'])
            for key in ('log_record','log_event'): number(row[key])
            for key in ('record_error','event_error'): number(row[key],0)
        if item['status']=='ok' and any(row['status']!='ok' for row in item['conditional']):
            raise ValueError('successful reducer has refused required conditional state')
        for key in ('joint_log_error','event_log_error','value_log_error'): number(item[key],0)
        for key in ('log_joint','log_event','log_value','product_per_exposure_log_record','product_per_exposure_log_event','record_relative_difference','event_relative_difference'):
            if item['status']=='ok': number(item[key])
            elif item[key] is not None: raise ValueError('refused reducer aggregate available')
    obj(value,'NativePayload')
    for key,wanted in [('schema','external-temporal-optical-sdk-control/v1'),('model_id','finite_bilinear_rest_spectral_time_zero_outside_full_observer_exposure_mean'),('detector_model_id','DETECTOR/Poisson-arrivals-fixed-QE-Gaussian-read/v1'),('selection_id','DETECTOR/threshold-joint-or-selected-censoring/v1')]: exact(value[key],wanted)
    for key in ('sampled_empirical_sensitivity','temporal_photon_endpoint_sensitivity'): number(value[key],0)
    status(value['prepared_status']); uint(value['grid_count'],1); uint(value['band_count'],4)
    for key in ('prepare_seconds','temporal_evaluate_seconds','total_native_seconds'): number(value[key],0)
    for key in ('retained_temporal_payload_bytes','detector_six_row_payload_bound_bytes'):
        if value[key] is not None: uint(value[key],MIB)
    admitted=temporal(value['primary']); array(value['detector_controls'],2 if admitted else 0)
    terms=0
    for noise,group in enumerate(value['detector_controls']):
        obj(group,'DetectorControl'); sigma=Decimal(0 if noise==0 else '0.75'); exact(number(group['sigma_electrons']),sigma)
        array(group['attempts'],8); array(group['records'],3)
        for i,attempt in enumerate(group['attempts']):
            obj(attempt,'DetectorAttempt'); exact(attempt['state_index'],i//4); exact(attempt['channel_index'],i%4)
            if attempt['lambda_electrons'] is not None: number(attempt['lambda_electrons'],0)
            array(attempt['batches'],3)
            nominal=value['primary']['rows'][i]['photons']['value']
            for j,b in enumerate(attempt['batches']):
                if b is None:
                    if j==0 or nominal!=0: raise ValueError('required endpoint missing')
                    continue
                terms+=batch(b,sigma,i%4)
                n=number(b['photons'],0)
                if j==0 and n!=nominal: raise ValueError('nominal photon ancestry')
                if j==1 and not 0<n<nominal: raise ValueError('lower endpoint order')
                if j==2 and not n>nominal: raise ValueError('upper endpoint order')
        for i,row in enumerate(group['records']): reducer(row); exact(row['record_index'],i); exact(row['selected_only'],i==1)
    uint(value['primary_poisson_terms'],12288); exact(value['primary_poisson_terms'],terms)
    array(value['actual_temporal_controls'],6)
    temporal_work=value['primary']['segment_work']
    for i,control in enumerate(value['actual_temporal_controls']):
        obj(control,'TemporalControl'); exact(control['id'],axes['temporal_control_ids'][i]); boolean(control['joint_withheld'])
        ok=temporal(control['native'],control['id']); exact(control['joint_withheld'],not ok)
        temporal_work+=control['native']['segment_work']
    array(value['actual_detector_controls'],4)
    for i,control in enumerate(value['actual_detector_controls']):
        obj(control,'SingleDetectorControl'); exact(control['id'],axes['detector_control_ids'][i])
        terms+=batch(control['native'],Decimal(0 if i==3 else '0.75'),single=True,bad_qe=i==0)
    array(value['actual_external_reducer_controls'],3 if admitted else 0)
    for i,control in enumerate(value['actual_external_reducer_controls']):
        obj(control,'ExternalReducerControl'); exact(control['id'],axes['reducer_control_ids'][i])
        array(control['retained_state_masses'],2)
        exact([number(v) for v in control['retained_state_masses']],[Decimal('0.25'),Decimal('0.75')])
        if i<2:
            if control['replacement_batch'] is None: raise ValueError('replacement missing')
            terms+=batch(control['replacement_batch'],Decimal('0.75'),3,bad_qe=i==0)
        else: exact(control['replacement_batch'],None)
        reducer(control['reducer'])
    uint(value['all_poisson_terms'],13824); exact(value['all_poisson_terms'],terms)
    uint(value['all_temporal_segment_work'],3584); exact(value['all_temporal_segment_work'],temporal_work)
    complete = admitted and value['prepared_status']=='ok'
    for group in value['detector_controls']:
        complete = complete and all(row['status']=='ok' for row in group['records'])
        for attempt in group['attempts']:
            complete = complete and all(b is not None and b['status']=='ok' and len(b['rows'])==6
                and all(row['status']=='ok' for row in b['rows']) for b in attempt['batches'])
    return complete


def runtime_admission(value, request):
    types=request['shape_contract']['runtime_types']; expected=request['reference_runtime']
    closed(value,types['RuntimeFingerprint']); exact(value['schema_version'],1); exact(value['interface_id'],'temporal-shared-optical-reference-runtime/v1')
    closed(value['python'],types['PythonIdentity']); closed(value['mpmath'],types['MpmathIdentity'])
    for key in ('version','implementation','requested_executable','resolved_executable','prefix','base_prefix'):
        exact(value['python'][key],expected[key])
    exact(value['python']['executable_sha256'],expected['executable_sha256']); exact(value['python']['dont_write_bytecode'],True)
    uint(value['python']['executable_bytes'],67108864)
    for key,wanted in [('version','1.3.0'),('backend','python'),('decimal_digits',60),('gmpy2_loaded',False),('original_subset_map_sha256',expected['legacy_package_inventory_sha256'])]: exact(value['mpmath'][key],wanted)
    uint(value['mpmath']['binary_precision_bits'],4096); exact(value['actual_thread_count'],1)
    exact(value['thread_environment'],request['resources']['environment'])
    inventory={}; total=0
    for item in array(value['inventory'],maximum=4096):
        closed(item,types['FileIdentity']); uint(item['bytes'],67108864); sha(item['sha256'])
        p=Path(item['path'])
        if not p.is_absolute() or str(p.resolve())!=str(p) or str(p) in inventory: raise ValueError('runtime path/duplicate')
        inventory[str(p)]=item; total+=item['bytes']
    if not inventory or total>67108864: raise ValueError('runtime inventory quota')
    exact(list(inventory),sorted(inventory))
    consumed=0
    for path,item in inventory.items():
        _,actual=file_identity(Path(path),67108864-consumed,item['sha256']); exact(actual['bytes'],item['bytes'])
        consumed+=actual['bytes']
        if consumed>67108864: raise ValueError('runtime actual inventory quota')
    canonical=json.dumps(value['inventory'],sort_keys=True,separators=(',',':'),allow_nan=False).encode()
    exact(value['inventory_sha256'],digest(canonical))
    executable=str(Path(expected['resolved_executable']).resolve())
    if executable not in inventory or inventory[executable]['sha256']!=expected['executable_sha256']: raise ValueError('runtime executable inventory')
    exact(value['python']['executable_bytes'],inventory[executable]['bytes'])
    cfg=value['python']['pyvenv_cfg']; closed(cfg,types['FileIdentity']); exact(cfg['sha256'],expected['pyvenv_cfg_sha256'])
    if cfg['path'] not in inventory or inventory[cfg['path']]!=cfg: raise ValueError('runtime cfg inventory')
    names=[]
    for module in array(value['imported_modules'],maximum=4096):
        closed(module,types['ImportedModuleIdentity']); names.append(module['name'])
        if type(module['name']) is not str or not module['name'].isascii() or len(module['name'])>256: raise ValueError('module name')
        if module['kind'] not in ('built-in','frozen','source','extension','namespace'): raise ValueError('module kind')
        if module['kind'] in ('source','extension'):
            if module['origin'] not in inventory: raise ValueError('imported module missing')
        elif module['origin'] is not None: raise ValueError('module origin')
        array(module['namespace_paths'],maximum=32)
        if module['namespace_paths']!=sorted(set(module['namespace_paths'])): raise ValueError('namespace paths')
        if module['kind']!='namespace' and module['namespace_paths']: raise ValueError('nonnamespace paths')
        for name in module['namespace_paths']:
            path=Path(name)
            if not path.is_absolute() or str(path.resolve())!=name or not path.is_dir(): raise ValueError('namespace directory identity')
    if names!=sorted(set(names)) or not names: raise ValueError('complete module order/uniqueness')
    mapped=[]
    for item in array(value['mapped_libraries'],maximum=256):
        closed(item,types['FileIdentity'])
        if item['path'] not in inventory or inventory[item['path']]!=item: raise ValueError('mapped library absent')
        mapped.append(item['path'])
    if mapped!=sorted(set(mapped)): raise ValueError('mapped order/uniqueness')
    configuration=value['numeric_configuration']
    closed(configuration,['arithmetic_backend','decimal_digits','binary_precision_bits','active_native_numeric_backends','gmpy2_loaded','byteorder'])
    exact(configuration['arithmetic_backend'],'mpmath-libmp-python'); exact(configuration['decimal_digits'],60)
    exact(configuration['binary_precision_bits'],value['mpmath']['binary_precision_bits'])
    exact(configuration['active_native_numeric_backends'],[]); exact(configuration['gmpy2_loaded'],False)
    if configuration['byteorder'] not in ('little','big'): raise ValueError('byteorder')
    root=Path(value['mpmath']['module_root']); origin=Path(value['mpmath']['module_origin'])
    if str(origin) not in inventory or inventory[str(origin)]['sha256']!=expected['mpmath_init_sha256']: raise ValueError('mpmath origin')
    if not root.is_absolute() or root.resolve()!=root or origin!=root/'__init__.py': raise ValueError('mpmath root identity')
    enumerated=0
    def bounded_files(directory):
        nonlocal enumerated
        stack=[directory]
        while stack:
            current=stack.pop(); walk_directory(current)
            with os.scandir(current) as entries:
                for entry in entries:
                    enumerated+=1
                    if enumerated>4096: raise ValueError('runtime package enumeration count')
                    if entry.is_symlink(): raise ValueError('runtime package link')
                    if entry.is_dir(follow_symlinks=False): stack.append(Path(entry.path))
                    elif entry.is_file(follow_symlinks=False): yield Path(entry.path)
                    else: raise ValueError('runtime package nonregular entry')
    source_paths=set()
    for path in bounded_files(root):
        if path.suffix=='.py': source_paths.add(str(path.resolve(strict=True)))
    if not source_paths.issubset(inventory): raise ValueError('incomplete mpmath sources')
    # Independently reconstruct the original pure-Python package/license subset.
    selected={Path(path) for path in source_paths}; distribution=root.parent/'mpmath-1.3.0.dist-info'
    walk_directory(distribution)
    for path in bounded_files(distribution):
        if str(path.resolve(strict=True)) not in inventory: raise ValueError('incomplete mpmath distribution metadata')
        if str(path).endswith(('LICENSE','METADATA','RECORD','WHEEL','top_level.txt')): selected.add(path)
    legacy={str(path):inventory[str(path.resolve())]['sha256'] for path in sorted(selected)}
    exact(digest(json.dumps(legacy,sort_keys=True,separators=(',',':')).encode()),expected['legacy_package_inventory_sha256'],'complete legacy runtime subset')
    return value


def discovery_admission(value, request):
    closed(value,['schema_version','product','executable','version','abi_version','build','capabilities','abi_schema','commands','interface_policy','scientific_qualifications'])
    exact(value['schema_version'],2); exact(value['product'],'Irreducible'); exact(value['executable'],'irred')
    root=Path(request['sdk']['artifacts_root'])
    manifest=parse(read(root/request['sdk']['manifest']['path']))
    exact(value['build'],manifest,'compiled discovery manifest')
    exact(value['version'],manifest['source_version'])
    abi=parse(read(Path(request['sdk']['source_root'])/'schema/abi.json'))
    exact(value['abi_version'],abi['abi_version']); exact(value['abi_schema'],abi)
    exact(value['commands'],['describe --json','version --json','run REQUEST STORE [--assurance numerical_contract|qualified]','stream STORE [LIMITS_JSON]'])
    exact(value['interface_policy'],'one current ABI revision; no compatibility aliases'); exact(value['scientific_qualifications'],[])
    ids=[]
    for item in array(value['capabilities'],maximum=64):
        if type(item) is not dict or not {'id','implementation','scientific','qualification'}.issubset(item): raise ValueError('discovery capability shape')
        if type(item['id']) is not str: raise ValueError('discovery capability id')
        exact(item['implementation'],'implemented'); boolean(item['scientific']); exact(item['qualification'],'unqualified'); ids.append(item['id'])
    if len(ids)!=len(set(ids)): raise ValueError('duplicate discovery capability')
    return value


def reference_admission(value, request, request_sha, source_hashes, runtime_before):
    """Admit the independently emitted transport object; never import reference math."""
    types=request['shape_contract']['reference_types']
    def obj(item,name): return closed(item,types[name],name)
    def scalars(item,name):
        obj(item,name)
        for key in types[name]: decimal_string(item[key], 'error' in key or key.endswith('_tail'))
    def frequency(item): scalars(item,'FrequencyRow')
    def detector_row(item,sigma,diagnostic=False):
        scalars(item,'DiscreteDetectorReferenceRow' if sigma==0 else 'ContinuousDetectorReferenceRow')
        if diagnostic: return
        if sigma==0:
            p=decimal_string(item['count_probability']); e=decimal_string(item['count_probability_error'],True)
            b=decimal_string(item['below']); be=decimal_string(item['below_error'],True)
            if not 0<p<=1 or not 0<b<=1 or e>=p or be>=b or not 1-b>be: raise ValueError('discrete positive reference gate')
        else:
            p=decimal_string(item['density64_per_adu']); e=decimal_string(item['density_error_per_adu'],True)
            b=decimal_string(item['below64']); be=decimal_string(item['below_error'],True)
            if not p>0 or not 0<b<=1 or e>=p or be>=b or not 1-b>be: raise ValueError('continuous positive reference gate')
            if not decimal_string(item['density32_per_adu'])>0 or not 0<decimal_string(item['below32'])<=1: raise ValueError('coarse detector domain')
    def group(item,i):
        obj(item,'DetectorReferenceGroup'); exact(number(item['sigma_electrons']),Decimal(0 if i==0 else '0.75'))
        for row in array(item['rows'],8): detector_row(row,i)
    def record(item,index):
        obj(item,'ReferenceRecord'); exact(item['sigma_index'],index//3); exact(item['record_index'],index%3)
        for key in ('conditional_record','conditional_event'):
            for text in array(item[key],2):
                n=decimal_string(text)
                if n<=0 or (key=='conditional_event' and n>1): raise ValueError('conditional probability domain')
        for key in types['ReferenceRecord']:
            if key in ('sigma_index','record_index','conditional_record','conditional_event'): continue
            decimal_string(item[key], 'error' in key)
        for key in ('joint','event'):
            n=decimal_string(item['positive_'+key]); e=decimal_string(item['positive_'+key+'_error'],True)
            if n<=0 or e>=n or (key=='event' and n>1): raise ValueError('positive mixture reference gate')
    def failure(item):
        if item is None: return
        if item.get('stage')=='positive-mixture-error-gate':
            obj(item,'PositiveMixtureFailure'); decimal_string(item['value']); decimal_string(item['error'],True)
            for text in array(item['conditional_values'],2): decimal_string(text)
            for text in array(item['conditional_errors'],2): decimal_string(text,True)
        elif item.get('stage')=='conditional-positive-reference-gate':
            obj(item,'ConditionalPositiveFailure'); uint(item['sigma_index'],1); uint(item['channel_index'],7)
            for key in ('lambda','density','nondetection'): decimal_string(item[key])
            for key in ('density_error','nondetection_error'): decimal_string(item[key],True)
            detector_row(item['coarse_refined_attempt'],item['sigma_index'],True)
        else: raise ValueError('scientific failure shape')
    obj(value,'ReferenceEnvelope'); exact(value['schema_version'],1); exact(value['interface_id'],'temporal-shared-optical-reference-transport/v1')
    if value['status'] not in ('completed','refused'): raise ValueError('reference status')
    axes={key:item for key,item in request['shape_contract']['axes'].items() if key!='unit_policy'}
    exact(value['axes'],axes); exact(value['settings'],{'scientific':request['shape_contract']['unchanged_scientific_settings'],'resources':request['resources']})
    obj(value['qualification'],'ReferenceQualification')
    exact(value['qualification']['scope'],'empirical-synthetic-fixed-temporal-shared-optical-control/v1')
    for key in ('native_outputs_consumed','physical_qualification','observational_qualification','original_input_certificate'): exact(value['qualification'][key],False)
    limits=array(value['qualification']['limits'],maximum=16)
    if not limits or any(type(item) is not str or len(item.encode())>8192 for item in limits): raise ValueError('qualification limits')
    obj(value['work'],'ReferenceWork')
    for key in ('elapsed_seconds','cpu_user_seconds','cpu_system_seconds'): number(value['work'][key],0)
    uint(value['work']['maximum_rss_kib']); uint(value['work']['time_frequency_evaluations'],30720); uint(value['work']['characteristic_function_evaluations'],73728)
    identities=value['identities']
    if identities is not None:
        obj(identities,'ReferenceIdentities')
        expected={'request_sha256':request_sha,'original_contract_sha256':CONTRACT_SHA,'admission_sha256':request['origin']['admission_sha256'],
                  'original_reference_sha256':request['origin']['original_sources'][2]['sha256'],'reference_script_sha256':source_hashes['reference.py'],'current_source_sha256':source_hashes}
        exact(identities,expected)
    for key in ('runtime_before','runtime_after'):
        if value[key] is not None: runtime_admission(value[key],request)
    if value['runtime_before'] is not None: exact(value['runtime_before'],runtime_before,'reference admitted runtime')
    partial=value['partial']; obj(partial,'ReferencePartial')
    order=['frequency','photons','lambda','detector-sigma0','detector-sigma075','records']
    groups=array(partial['completed_groups'],maximum=6)
    exact(groups,order[:len(groups)],'reference completed group order')
    for text in array(partial['photon_prefix'],maximum=8):
        if decimal_string(text)<=0: raise ValueError('partial photons')
    for text in array(partial['photon_error_prefix'],maximum=8): decimal_string(text,True)
    exact(len(partial['photon_prefix']),len(partial['photon_error_prefix']))
    for row in array(partial['frequency_prefix'],maximum=8): frequency(row)
    exact(len(partial['photon_prefix']),len(partial['frequency_prefix']),'photon/frequency earned prefix')
    for text in array(partial['lambda_prefix'],maximum=8):
        if decimal_string(text)<=0: raise ValueError('partial lambda')
    for i,item in enumerate(array(partial['detector_prefix'],maximum=2)): group(item,i)
    active=partial['active_sigma_index']
    if active is not None: uint(active,1)
    active_rows=array(partial['active_sigma_rows'],maximum=8)
    if active is None and active_rows: raise ValueError('active sigma absent')
    for row in active_rows: detector_row(row,active)
    for i,item in enumerate(array(partial['record_prefix'],maximum=6)): record(item,i)
    failure(partial['control_failure'])
    # Bind the unchanged observer's phases, including sigma0 records before sigma075.
    count=len(groups); detector_count=len(partial['detector_prefix']); record_count=len(partial['record_prefix'])
    if count==1: raise ValueError('frequency/photon completion is one observer operation')
    if count>=2:
        exact(len(partial['photon_prefix']),8,'completed photon/frequency prefix')
    exact(len(partial['lambda_prefix']),8 if count>=3 else 0,'lambda phase')
    exact(detector_count,2 if count>=5 else 1 if count==4 else 0,'completed sigma phase')
    if detector_count==0 and record_count!=0: raise ValueError('records precede completed sigma0')
    if detector_count==1 and record_count>3: raise ValueError('sigma075 records precede sigma075 completion')
    if detector_count==2 and record_count<3: raise ValueError('sigma075 precedes all sigma0 records')
    if (record_count==6)!=(count==6): raise ValueError('records completion marker differs')
    if active is not None:
        exact(count,3 if active==0 else 4,'active sigma phase')
        exact(record_count,0 if active==0 else 3,'active sigma preceding records')
    if partial['control_failure'] is not None:
        if record_count==6: raise ValueError('completed records retain scientific failure')
        problem=partial['control_failure']
        if problem['stage']=='conditional-positive-reference-gate':
            exact(active,problem['sigma_index'],'failed conditional active sigma')
            exact(len(active_rows),problem['channel_index'],'failed conditional earned row prefix')
        else:
            if active is not None or detector_count==0: raise ValueError('positive-mixture failure outside record phase')
    for label,field,n in [('frequency','frequency_prefix',8),('photons','photon_prefix',8),('lambda','lambda_prefix',8),('detector-sigma0','detector_prefix',1),('detector-sigma075','detector_prefix',2),('records','record_prefix',6)]:
        if label in groups and len(partial[field])<n: raise ValueError('completed group missing earned rows')
    payload=value['payload']
    if payload is not None:
        obj(payload,'ReferencePayload'); exact(payload['schema'],'original-high-precision-temporal-optical-reference/v1'); exact(payload['mpmath_version'],'1.3.0'); exact(payload['decimal_digits'],60)
        exact(payload['shared_ancestry'],['declared equations and exact SI constants','mpmath arithmetic/GL node generator; no native kernels'])
        for key in ('frequency_seconds','elapsed_seconds'): number(payload[key],0)
        exact(payload['integrand_evaluations'],{'time_frequency':30720,'characteristic_function':73728})
        for field in ('photons','lambda_electrons'):
            for text in array(payload[field],8):
                if decimal_string(text)<=0: raise ValueError('positive reference scalar')
        for text in array(payload['photon_errors'],8): decimal_string(text,True)
        for row in array(payload['frequency_rows'],8): frequency(row)
        for i,item in enumerate(array(payload['detector_rows'],2)): group(item,i)
        for i,item in enumerate(array(payload['records'],6)): record(item,i)
        exact(groups,order)
        for field,other in [('photon_prefix','photons'),('photon_error_prefix','photon_errors'),('frequency_prefix','frequency_rows'),('lambda_prefix','lambda_electrons'),('detector_prefix','detector_rows'),('record_prefix','records')]: exact(partial[field],payload[other])
        exact(value['work']['time_frequency_evaluations'],30720); exact(value['work']['characteristic_function_evaluations'],73728)
    error=value['error']
    if error is not None:
        obj(error,'ReferenceError')
        if error['kind'] not in ('admission','runtime','reference','io','resource'): raise ValueError('reference error kind')
        for key in ('stage','message'):
            if type(error[key]) is not str or len(error[key].encode())>8192: raise ValueError('error text')
        boolean(error['message_truncated']); sha(error['message_sha256'])
        if not error['message_truncated']: exact(digest(error['message'].encode()),error['message_sha256'])
        details=error['details']
        if details is not None:
            if details.get('stage') in ('positive-mixture-error-gate','conditional-positive-reference-gate'): failure(details)
            else:
                obj(details,'TransportFailureDetails')
                for key in ('operation','reason'):
                    if type(details[key]) is not str: raise ValueError('transport error text')
                if details['path'] is not None and type(details['path']) is not str: raise ValueError('transport error path')
    if value['status']=='completed':
        if payload is None or identities is None or error is not None or value['runtime_before'] is None or value['runtime_after'] is None: raise ValueError('completed reference missing earned evidence')
        if active is not None or active_rows or partial['control_failure'] is not None: raise ValueError('completed reference retains active/failed scientific phase')
        exact(value['runtime_before'],value['runtime_after'],'reference runtime before/after')
        return True
    if error is None: raise ValueError('refused reference has no cause')
    return False
