"""Published age-table extraction and explicitly conditional regression."""

import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.linalg import cho_factor, cho_solve


def read_age(p, label):
    rows = []
    for line in p.read_text().splitlines()[1:]:
        c = [v.strip().replace("\\", "").strip() for v in line.split("&")]
        if len(c) != 5:
            continue
        hr = c[3].split("(")
        rows.append(
            [
                str(int(c[0])),
                label,
                float(c[1]),
                float(c[2]),
                float(hr[0]),
                float(hr[-1].strip(") ")),
                float(c[4]),
            ]
        )
    return pd.DataFrame(
        rows,
        columns=[
            "CID",
            "age_source",
            "age",
            "age_err",
            "hr_C25",
            "hr_original",
            "hr_C25_err",
        ],
    )


def mu(z, zhel=None, om=0.3):
    z = np.asarray(z)
    zhel = z if zhel is None else np.asarray(zhel)
    d = np.array(
        [
            quad(
                lambda t: 1 / np.sqrt(om * (1 + t) ** 3 + 1 - om), 0, zz, epsabs=1e-11
            )[0]
            for zz in z
        ]
    )
    return 5 * np.log10((1 + zhel) * 299792.458 / 70 * d) + 25


def wls(x, y, e, X=None, cov=None):
    if X is None:
        X = np.column_stack([np.ones(len(x)), x])
    if cov is None:
        wx = X / e[:, None]
        wy = y / e
        c = np.linalg.inv(wx.T @ wx)
        b = c @ wx.T @ wy
        chi = float(np.sum(((y - X @ b) / e) ** 2))
    else:
        cf = cho_factor(cov)
        c = np.linalg.inv(X.T @ cho_solve(cf, X))
        b = c @ X.T @ cho_solve(cf, y)
        chi = float((y - X @ b) @ cho_solve(cf, y - X @ b))
    return b, c, chi
