"""Explicit, numerically exact reassociation of the modern lensing likelihood."""
import hashlib
import json
from pathlib import Path
from modern_run import configuration as original_configuration
from target_identity import identify as original_identify, digest, ROOT
from fast_lensing import use_fast_lensing


def configuration(*args, **kwargs):
    return use_fast_lensing(original_configuration(*args, **kwargs))


def identify(info, sample_file, model_file=None):
    record = original_identify(info, sample_file, model_file)
    for path in [Path(__file__), Path(__file__).parent.parent/'external_probes/fast_lensing.py']:
        record['source_sha256'][str(path.relative_to(ROOT))] = digest(path)
    record.pop('identity')
    record['identity'] = hashlib.sha256(json.dumps(record, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    return record
