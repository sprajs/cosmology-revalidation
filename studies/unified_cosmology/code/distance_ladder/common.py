"""Paths and immutable local evidence for the separate released distance ladder."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[4]
HERE=Path(__file__).resolve().parent
WORK=ROOT/'.work/unified-cosmology/distance-ladder'
OUT=ROOT/'studies/unified_cosmology/results/distance_ladder'
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def relative(path):return str(Path(path).resolve().relative_to(ROOT))
def write(path,value):
 path.parent.mkdir(parents=True,exist_ok=True)
 path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
