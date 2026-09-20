from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, subprocess, sys
import numpy, scipy, pandas

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'runs/mapping'
OUT.mkdir(parents=True, exist_ok=True)

def file_hash(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        while block := f.read(1024 * 1024):
            h.update(block)
    return h.hexdigest()

def manifest(name, config, inputs, outputs, extra=None):
    files = list(dict.fromkeys([str(Path(p).relative_to(ROOT)) if Path(p).is_absolute() else str(p) for p in inputs]))
    result = dict(name=name, utc=datetime.now(timezone.utc).isoformat(), config=config,
        python=sys.version, versions={'numpy': numpy.__version__, 'scipy': scipy.__version__, 'pandas': pandas.__version__},
        code_revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        input_sha256={p: file_hash(ROOT/p) for p in files},
        output_sha256={str(Path(p).relative_to(ROOT)): file_hash(p) for p in outputs}, extra=extra or {})
    (OUT / f'{name}-manifest.json').write_text(json.dumps(result, indent=2)+'\n')

