"""Paths and immutable identities for the unified survey-selection audit."""
from pathlib import Path
import hashlib
ROOT=Path(__file__).resolve().parents[4]
HERE=Path(__file__).resolve().parent
WORK=ROOT/'.work/unified-cosmology/survey-selection'
RESULTS=ROOT/'studies/unified_cosmology/results/survey_selection'
for p in [WORK,RESULTS]:p.mkdir(parents=True,exist_ok=True)
def sha(path,algorithm='sha256'):
    h=hashlib.new(algorithm)
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(2**20),b''):h.update(b)
    return h.hexdigest()
