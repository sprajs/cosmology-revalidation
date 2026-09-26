#!/usr/bin/env python3
"""Independent 1-D quadrature check of fixed circle-square fractions; synthetic only."""
import json
import math
from pathlib import Path

from scipy.integrate import quad

from fixed_mask import R0, R1, R2, pixel_circle_fraction

HERE = Path(__file__).resolve().parent
CASES = ((R0, (0, 0), (3, 0), (2, 2), (3, 1)),
         (R1, (0, 0), (9, 0), (7, 6), (9, 3)),
         (R2, (0, 0), (15, 0), (12, 10), (15, 4)))


def quadrature(x: int, y: int, radius: float) -> float:
    def vertical_length(u: float) -> float:
        if abs(u) >= radius:
            return 0.0
        edge = math.sqrt(max(radius * radius - u * u, 0.0))
        return max(0.0, min(y + 0.5, edge) - max(y - 0.5, -edge))
    points = [z for z in (0.0, -radius, radius,
                          -math.sqrt(max(radius * radius - (y - 0.5) ** 2, 0.0)),
                          math.sqrt(max(radius * radius - (y - 0.5) ** 2, 0.0)),
                          -math.sqrt(max(radius * radius - (y + 0.5) ** 2, 0.0)),
                          math.sqrt(max(radius * radius - (y + 0.5) ** 2, 0.0)))
              if x - 0.5 < z < x + 0.5]
    return quad(vertical_length, x - 0.5, x + 0.5,
                points=sorted(set(points)), epsabs=1e-12, epsrel=1e-12, limit=100)[0]


rows = []
for radius, *cells in CASES:
    for x, y in cells:
        analytic = pixel_circle_fraction(x, y, radius)
        numeric = quadrature(x, y, radius)
        rows.append({"radius": radius, "x": x, "y": y,
                     "analytic": analytic, "independent_quadrature": numeric,
                     "absolute_gap": abs(analytic - numeric)})
maximum = max(row["absolute_gap"] for row in rows)
result = {"pass": maximum <= 1e-9, "maximum_absolute_gap": maximum,
          "cases": rows, "real_SCI_read": False}
(HERE / "quadrature-check.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps({"pass": result["pass"], "maximum_absolute_gap": maximum, "cases": len(rows)}))
if not result["pass"]:
    raise SystemExit(1)
