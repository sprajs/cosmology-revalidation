"""Paths and hashing for the survey-physics experiment; no downloaded code is vendored."""

import hashlib
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
HERE = Path(__file__).resolve().parent
WORK = ROOT / ".work/survey-physics"
RESULTS = ROOT / "studies/host_ages/results/survey_physics"
ARCHIVE = Path(
    os.environ.get(
        "SURVEY_PHYSICS_ARCHIVE",
        str(
            Path.home()
            / ".local/share/cosmology-revalidation/archive/17487bf659fcbdeeea072221492bac14b04a0a85"
        ),
    )
)
DATA = WORK / "SNDATA_ROOT"


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(2**20), b""):
            h.update(block)
    return h.hexdigest()


def native_env():
    env = os.environ.copy()
    env.update(
        SNANA_DIR=str(WORK / "SNANA"),
        SNDATA_ROOT=str(DATA),
        LD_LIBRARY_PATH=os.environ.get(
            "SURVEY_PHYSICS_LIBRARY_PATH",
            str(ARCHIVE / "phase2/official/build/sysroot/usr/lib"),
        ),
        OPENBLAS_NUM_THREADS="1",
        OMP_NUM_THREADS="1",
        MKL_NUM_THREADS="1",
    )
    return env
