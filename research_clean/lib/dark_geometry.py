"""Reference-only detector masks and exact circle/pixel area primitives."""

import math
import numpy as np
from astropy.io import fits
from pathlib import Path

N = 1024
ACTIVE = slice(5, 1019)


def map_bpixtab(path: Path) -> tuple[np.ndarray, dict]:
    """Exact IR doDQIIR/DQINormal map: one-index PIX, then raw LTV=+5."""
    mask = np.zeros((N, N), dtype=np.uint16)
    rows_used = 0
    clamped = 0
    with fits.open(path, memmap=True) as hdus:
        table = hdus["BPIX"]
        if table.header["SIZAXIS1"] != N or table.header["SIZAXIS2"] != N:
            raise ValueError("BPIXTAB detector dimension mismatch")
        for row in table.data:
            amp = str(row["CCDAMP"]).strip()
            chip = int(row["CCDCHIP"])
            gain = float(row["CCDGAIN"])
            if (
                amp not in ("N/A", "ABCD")
                or chip not in (-999, 1)
                or gain not in (-999.0, 2.5)
            ):
                continue
            x = int(row["PIX1"]) - 1 + 5  # raw SCI LTV1=5, LTM1_1=1
            y = int(row["PIX2"]) - 1 + 5
            length = int(row["LENGTH"])
            axis = int(row["AXIS"])
            flag_float = float(row["VALUE"])
            flag = int(flag_float)
            if (
                flag_float != flag
                or not 0 <= flag <= 65535
                or length <= 0
                or axis not in (1, 2)
            ):
                raise ValueError("invalid BPIXTAB row")
            rows_used += 1
            if axis == 1:
                if y < 0 or y >= N or x + length - 1 < 0 or x >= N:
                    clamped += 1
                    continue
                lo, hi = max(x, 0), min(x + length, N)
                clamped += int(lo != x or hi != x + length)
                mask[y, lo:hi] |= np.uint16(flag)
            else:
                if x < 0 or x >= N or y + length - 1 < 0 or y >= N:
                    clamped += 1
                    continue
                lo, hi = max(y, 0), min(y + length, N)
                clamped += int(lo != y or hi != y + length)
                mask[lo:hi, x] |= np.uint16(flag)
    return mask, {
        "rows_used": rows_used,
        "out_of_frame_or_clamped_rows": clamped,
        "nonzero_pixels_all_raw": int(np.count_nonzero(mask)),
        "nonzero_pixels_active": int(np.count_nonzero(mask[ACTIVE, ACTIVE])),
    }


def quadrant_area(x: float, y: float, radius: float) -> float:
    """Area of circle in [0,x]×[0,y] for x,y>=0."""
    x, y = min(max(x, 0.0), radius), min(max(y, 0.0), radius)
    if x == 0 or y == 0:
        return 0.0
    if x * x + y * y <= radius * radius:
        return x * y
    x0 = min(x, math.sqrt(max(radius * radius - y * y, 0.0)))

    def integral(t: float) -> float:
        return 0.5 * (
            t * math.sqrt(max(radius * radius - t * t, 0.0))
            + radius * radius * math.asin(min(1.0, t / radius))
        )

    return y * x0 + integral(x) - integral(x0)


def signed_area(x: float, y: float, radius: float) -> float:
    return (
        math.copysign(1.0, x)
        * math.copysign(1.0, y)
        * quadrant_area(abs(x), abs(y), radius)
    )


def pixel_circle_fraction(x: int, y: int, radius: float) -> float:
    x0, x1, y0, y1 = x - 0.5, x + 0.5, y - 0.5, y + 0.5
    return (
        signed_area(x1, y1, radius)
        - signed_area(x0, y1, radius)
        - signed_area(x1, y0, radius)
        + signed_area(x0, y0, radius)
    )
