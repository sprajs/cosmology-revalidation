"""CLASS file and unit adapters. Cosmological calculations stay in CLASS."""
import bisect
import math
from pathlib import Path
import re

BASE_NS = "0.9660499"
REDSHIFTS = (0.0, 0.38, 0.51, 0.61)


def class_parameters(ns=BASE_NS):
    """The supplied CLASS point with explicit defaults and requested products."""
    if ns not in (BASE_NS, "0.9610499", "0.9710499"):
        raise ValueError("only the frozen anchor and ns +/-0.005 are supported")
    return {
        "H0": "67.32117", "omega_b": "0.02238280", "omega_cdm": "0.1201075",
        "T_cmb": "2.7255", "N_ur": "2.046", "N_ncdm": "1", "m_ncdm": "0.06",
        "T_ncdm": "0.7137658555036082", "deg_ncdm": "1", "ksi_ncdm": "0",
        "Omega_k": "0", "Omega_fld": "0", "Omega_scf": "0",
        "YHe": "0.2454006", "recombination": "HyRec", "reio_parametrization": "reio_camb",
        "tau_reio": "0.05430842", "reionization_exponent": "1.5", "reionization_width": "0.5",
        "helium_fullreio_redshift": "3.5", "helium_fullreio_width": "0.5",
        "P_k_ini type": "analytic_Pk", "A_s": "2.100549e-09", "n_s": ns,
        "alpha_s": "0", "k_pivot": "0.05",
        "modes": "s", "ic": "ad", "gauge": "synchronous",
        "output": "tCl,pCl,lCl,mPk", "lensing": "yes", "non linear": "halofit",
        "l_max_scalars": "2508", "P_k_max_1/Mpc": "1", "z_pk": "0,0.38,0.51,0.61",
        "format": "class", "headers": "yes", "write_background": "yes",
        "write_thermodynamics": "yes", "write_parameters": "yes", "write_warnings": "yes",
        "overwrite_root": "yes", "input_verbose": "1", "background_verbose": "1",
        "thermodynamics_verbose": "1", "perturbations_verbose": "1", "transfer_verbose": "1",
        "primordial_verbose": "1", "harmonic_verbose": "1", "fourier_verbose": "1",
        "lensing_verbose": "1", "output_verbose": "1",
    }


def render_ini(parameters, output_root):
    if parameters != class_parameters(parameters.get("n_s")):
        raise ValueError("parameter dictionary differs from the frozen physical point")
    root = str(Path(output_root).absolute())
    if any(c in root for c in "\n\r#=") or len(root.encode()) > 400:
        raise ValueError("invalid CLASS output root")
    lines = ["# Illustrative supplied CLASS point; no likelihood best-fit claim"]
    lines += [f"{key} = {value}" for key, value in parameters.items()]
    lines.append("root = " + root)
    return ("\n".join(lines) + "\n").encode("ascii")


def parse_table(raw, max_rows=100000):
    """Keep emitted order and numbered column titles; refuse incomplete rows."""
    if type(raw) is not bytes:
        raise ValueError("table input must be immutable bytes")
    headers, columns, rows = [], None, []
    for line in raw.decode("ascii").splitlines():
        if not line.strip():
            continue
        if len(line) > 65536:
            raise ValueError("table line too long")
        if line.lstrip().startswith("#"):
            headers.append(line)
            text = line.lstrip()[1:]
            markers = list(re.finditer(r"(?:^|\s)([0-9]+):", text))
            if markers:
                labels = [text[m.end():markers[i + 1].start() if i + 1 < len(markers)
                               else len(text)].strip() for i, m in enumerate(markers)]
                if [int(m.group(1)) for m in markers] != list(range(1, len(markers) + 1)):
                    raise ValueError("nonconsecutive column numbering")
                if columns is not None or len(set(labels)) != len(labels) or len(labels) > 96:
                    raise ValueError("duplicate or over-bound column titles")
                columns = labels
            continue
        if columns is None or len(rows) >= max_rows:
            raise ValueError("missing header or table row bound")
        fields = line.split()
        if len(fields) != len(columns):
            raise ValueError("table column mismatch")
        row = [float(x) for x in fields]
        if not all(math.isfinite(x) for x in row):
            raise ValueError("nonfinite table value")
        rows.append(row)
    if columns is None or not rows:
        raise ValueError("empty CLASS table")
    return {"headers": headers, "columns": columns, "rows": rows}


def read_table(path, max_bytes, max_rows=100000):
    # Controller supplies its already bounded, identity-checked file reader.
    with Path(path).open("rb") as stream:
        raw = stream.read(max_bytes + 1)
    if len(raw) > max_bytes:
        raise ValueError("table byte bound")
    return parse_table(raw, max_rows)


def cmb_spectra(table, tcmb=2.7255, lmax=2508):
    if not any("dimensionless total lensed [l(l+1)/2pi]" in x for x in table["headers"]):
        raise ValueError("expected dimensionless CLASS-format lensed Dl")
    columns = table["columns"]
    if columns[:4] != ["l", "TT", "EE", "TE"]:
        raise ValueError("unexpected CLASS spectrum order")
    rows = table["rows"]
    if [x[0] for x in rows] != list(range(2, lmax + 1)):
        raise ValueError("missing, duplicate or reordered multipoles")
    if not (math.isfinite(tcmb) and tcmb > 0):
        raise ValueError("invalid CMB temperature")
    ell = [int(x[0]) for x in rows]
    result = {"ell": ell, "raw_units": "dimensionless Dl=l(l+1)Cl/(2pi)",
              "Tcmb_K": tcmb, "lensed": True, "spectra": {}}
    for name in ("TT", "EE", "TE"):
        dl = [row[columns.index(name)] for row in rows]
        if name != "TE" and any(value < 0 for value in dl):
            raise ValueError("negative CMB auto-spectrum")
        cl = [value * 2 * math.pi / (l * (l + 1)) for l, value in zip(ell, dl)]
        scale = (tcmb * 1e6) ** 2
        result["spectra"][name.lower()] = {"Dl_dimensionless": dl, "Cl_dimensionless": cl,
                                         "Dl_uK2": [x * scale for x in dl],
                                         "Cl_uK2": [x * scale for x in cl]}
    if "phiphi" in columns:
        index = columns.index("phiphi")
        result["phiphi_Cl_dimensionless"] = [row[index] * 2 * math.pi / (l * (l + 1))
                                             for l, row in zip(ell, rows)]
    return result


def matter_spectrum(table, h, expected_z):
    if table["columns"] != ["k (h/Mpc)", "P (Mpc/h)^3"] or not (math.isfinite(h) and h > 0):
        raise ValueError("unexpected matter spectrum units")
    markers = [re.search(r"at redshift z=([^\s]+)", x) for x in table["headers"]]
    redshifts = [float(x.group(1)) for x in markers if x is not None]
    if redshifts != [expected_z]:
        raise ValueError("matter redshift mapping differs")
    rows = table["rows"]
    if any(k <= 0 or p <= 0 for k, p in rows) or any(rows[i][0] >= rows[i + 1][0]
                                                   for i in range(len(rows) - 1)):
        raise ValueError("nonpositive or unordered matter spectrum")
    return {"z": expected_z, "k_1_Mpc": [x[0] * h for x in rows],
            "P_Mpc3": [x[1] / h ** 3 for x in rows],
            "raw_k_unit": "h/Mpc", "raw_P_unit": "(Mpc/h)^3", "h": h}


def _background_grid(table):
    columns, rows = table["columns"], table["rows"]
    names = ("z", "H [1/Mpc]", "comov. dist.", "ang.diam.dist.", "lum. dist.")
    if any(name not in columns for name in names):
        raise ValueError("missing background unit/axis")
    indices = [columns.index(name) for name in names]
    selected = [[row[i] for i in indices] for row in rows]
    if all(selected[i][0] > selected[i + 1][0] for i in range(len(selected) - 1)):
        selected.reverse()
    if any(selected[i][0] >= selected[i + 1][0] for i in range(len(selected) - 1)):
        raise ValueError("background redshift axis is not monotonic")
    return selected


def _background_at_grid(selected, grid, z):
    if not (math.isfinite(z) and z >= 0):
        raise ValueError("invalid target redshift")
    index = bisect.bisect_left(grid, z)
    if index < len(grid) and grid[index] == z:
        lower = upper = selected[index]
        values = lower[1:]
    elif index == 0 or index == len(grid):
        raise ValueError("background target outside emitted grid")
    else:
        lower, upper = selected[index - 1:index + 1]
        weight = (z - lower[0]) / (upper[0] - lower[0])
        values = [a + weight * (b - a) for a, b in zip(lower[1:], upper[1:])]
    hubble, dm, da, dl = values
    if hubble <= 0 or min(dm, da, dl) < 0:
        raise ValueError("invalid emitted background domain")
    return {"z": z, "H_1_Mpc": hubble, "H_km_s_Mpc": hubble * 299792.458,
            "chi_Mpc": dm, "DM_Mpc": dm, "DM_scope": "flat-LCDM: transverse equals radial",
            "DA_Mpc": da, "DL_Mpc": dl,
            "bracket_z": [lower[0], upper[0]], "interpolation": "linear-on-emitted-background-table",
            "interpolation_error_bound": None}


def background_at(table, z):
    """Interpolate emitted background data; no new expansion calculation."""
    selected = _background_grid(table)
    return _background_at_grid(selected, [x[0] for x in selected], z)


def background_at_many(table, redshifts):
    selected = _background_grid(table)
    grid = [x[0] for x in selected]
    return [_background_at_grid(selected, grid, z) for z in redshifts]


def predicted_drag(logbytes):
    text = logbytes.decode("ascii")
    pattern = (r"baryon drag stops at z = ([0-9]+\.[0-9]{6})\s*"
               r"corresponding to conformal time = ([0-9]+\.[0-9]{6}) Mpc\s*"
               r"with comoving sound horizon rs = ([0-9]+\.[0-9]{6}) Mpc")
    matches = re.findall(pattern, text)
    if len(matches) != 1:
        raise ValueError("missing or ambiguous CLASS drag report")
    z, conformal, ruler = map(float, matches[0])
    if z <= 0 or ruler <= 0:
        raise ValueError("nonpositive predicted drag")
    return {"z_drag": z, "r_drag_Mpc": ruler, "conformal_time_Mpc": conformal,
            "source": "CLASS-thermodynamics-verbose",
            "rounding_half_width_z": 5e-7, "rounding_half_width_Mpc": 5e-7,
            "internal_unrounded_value_available": False}


def precision_difference(first, second):
    """Empirical same-solver difference, with a stable scale for signed TE."""
    if first["ell"] != second["ell"]:
        raise ValueError("precision comparison multipoles differ")
    result = {}
    for name in ("tt", "ee", "te"):
        a, b = first["spectra"][name]["Dl_uK2"], second["spectra"][name]["Dl_uK2"]
        scale = max(max(abs(x) for x in a), max(abs(x) for x in b))
        absolute = max(abs(x - y) for x, y in zip(a, b))
        result[name] = {"max_abs_Dl_uK2": absolute, "max_abs_over_global_spectrum_scale":
                        absolute / scale if scale else 0.0,
                        "scale_Dl_uK2": scale, "certified_error_bound": None}
    return result


def product_difference(first, second):
    result = {"cmb": precision_difference(first["cmb"], second["cmb"]),
              "drag": {key: second["drag"][key] - first["drag"][key]
                       for key in ("z_drag", "r_drag_Mpc")},
              "background": [], "linear_matter": [], "certified_error_bound": None}
    for a, b in zip(first["background"], second["background"]):
        if a["z"] != b["z"]:
            raise ValueError("background comparison redshifts differ")
        result["background"].append({"z": a["z"], "differences": {
            key: b[key] - a[key] for key in ("H_km_s_Mpc", "DM_Mpc", "DA_Mpc", "DL_Mpc")}})
    for a, b in zip(first["linear_matter"], second["linear_matter"]):
        if a["z"] != b["z"]:
            raise ValueError("matter comparison redshifts differ")
        values = []
        for k, power in zip(a["k_1_Mpc"], a["P_Mpc3"]):
            index = bisect.bisect_left(b["k_1_Mpc"], k)
            if index < len(b["k_1_Mpc"]) and b["k_1_Mpc"][index] == k:
                other = b["P_Mpc3"][index]
            elif index == 0 or index == len(b["k_1_Mpc"]):
                continue
            else:
                lo, hi = b["k_1_Mpc"][index - 1:index + 1]
                x, y = b["P_Mpc3"][index - 1:index + 1]
                other = x + (k - lo) * (y - x) / (hi - lo)
            values.append(abs(other - power) / abs(power))
        if not values:
            raise ValueError("no shared matter grid domain")
        result["linear_matter"].append({"z": a["z"], "max_abs_fractional_difference": max(values),
                                       "compared_anchor_rows": len(values),
                                       "excluded_anchor_rows": len(a["k_1_Mpc"]) - len(values),
                                       "grid_alignment": "linear-P interpolation on second grid; no extrapolation",
                                       "interpolation_error_bound": None})
    return result
