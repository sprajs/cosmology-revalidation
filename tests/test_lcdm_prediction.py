"""Transport fixtures for the external CLASS prediction reference."""
import importlib.util
import hashlib
import json
from pathlib import Path
import math
import sys
import tempfile
import unittest
from unittest.mock import patch

FULL_FOLDER = Path(__file__).resolve().parents[1] / "experiments/lcdm-reference"


def full_helper(name):
    """Load an adapter without installing theory or likelihood runtimes."""
    helper_spec = importlib.util.spec_from_file_location(
        "full_lcdm_fixture_" + name, FULL_FOLDER / (name + ".py"))
    helper = importlib.util.module_from_spec(helper_spec)
    helper_spec.loader.exec_module(helper)
    return helper



class FullLCDMTheoryMapping(unittest.TestCase):
    """Unit/axis fixtures, with no cosmological solver or accuracy claim."""

    def setUp(self):
        with patch.dict(sys.modules, {"classy": None, "numpy": None}):
            self.helper = full_helper("theory")

    def cmb(self):
        return self.helper.parse_table(
            b"# dimensionless total lensed [l(l+1)/2pi] C_l's\n"
            b"# 1:l 2:TT 3:EE 4:TE 5:BB 6:phiphi\n"
            b"2 3 6 -9 0 12\n3 6 12 0 0 24\n")

    def test_exact_parameter_mapping_and_response_cases(self):
        base = self.helper.class_parameters()
        self.assertEqual({key: base[key] for key in (
            "H0", "omega_b", "omega_cdm", "N_ur", "N_ncdm", "m_ncdm", "T_ncdm", "deg_ncdm")},
            {"H0": "67.32117", "omega_b": "0.02238280", "omega_cdm": "0.1201075",
             "N_ur": "2.046", "N_ncdm": "1", "m_ncdm": "0.06",
             "T_ncdm": "0.7137658555036082", "deg_ncdm": "1"})
        self.assertEqual(base["l_max_scalars"], "2508")
        self.assertEqual(base["z_pk"], "0,0.38,0.51,0.61")
        self.assertEqual(base["non linear"], "halofit")
        self.assertEqual(base["format"], "class")
        self.assertNotIn("z_drag", base)
        for ns in ("0.9610499", "0.9710499"):
            changed = self.helper.class_parameters(ns)
            self.assertEqual({k for k in base if base[k] != changed[k]}, {"n_s"})
        with self.assertRaises(ValueError):
            self.helper.class_parameters("0.96")

    def test_ini_refuses_changed_physics_and_line_injection(self):
        base = self.helper.class_parameters()
        rendered = self.helper.render_ini(base, "/tmp/fixture_class_")
        self.assertIn(b"omega_b = 0.02238280\n", rendered)
        self.assertIn(b"root = /tmp/fixture_class_\n", rendered)
        changed = dict(base, N_ur="3.044")
        with self.assertRaises(ValueError):
            self.helper.render_ini(changed, "/tmp/fixture_class_")
        with self.assertRaises(ValueError):
            self.helper.render_ini(base, "/tmp/fixture\nN_ur=0")

    def test_cmb_class_dl_to_cl_and_temperature_scaling(self):
        result = self.helper.cmb_spectra(self.cmb(), tcmb=2., lmax=3)
        self.assertEqual(result["ell"], [2, 3])
        self.assertAlmostEqual(result["spectra"]["tt"]["Cl_dimensionless"][0], math.pi)
        self.assertEqual(result["spectra"]["tt"]["Dl_uK2"], [12e12, 24e12])
        self.assertLess(result["spectra"]["te"]["Cl_uK2"][0], 0.)
        # Potential spectra carry no temperature-unit multiplication.
        self.assertAlmostEqual(result["phiphi_Cl_dimensionless"][0], 4 * math.pi)
        self.assertTrue(result["lensed"])

    def test_cmb_refuses_camb_format_reordered_or_missing_multipoles(self):
        table = self.cmb()
        for change in (
            lambda t: t.update(headers=["# total lensed (units: [microK]^2)"]),
            lambda t: t.update(columns=["l", "TT", "EE", "BB", "TE", "phiphi"]),
            lambda t: t["rows"].reverse(),
            lambda t: t["rows"].pop(),
            lambda t: t["rows"][0].__setitem__(1, -1.)):
            bad = {"headers": list(table["headers"]), "columns": list(table["columns"]),
                   "rows": [list(row) for row in table["rows"]]}
            change(bad)
            with self.assertRaises(ValueError):
                self.helper.cmb_spectra(bad, lmax=3)

    def test_table_rejects_nonfinite_duplicate_titles_and_row_mismatch(self):
        for raw in (b"# 1:x 2:y\n1 nan\n", b"# 1:x 2:x\n1 2\n",
                    b"# 1:x 3:y\n1 2\n", b"# 1:x 2:y\n1\n",
                    b"1 2\n", b"# 1:x 2:y\n"):
            with self.assertRaises(ValueError):
                self.helper.parse_table(raw)
        with self.assertRaises(ValueError):
            self.helper.parse_table(b"# 1:x\n1\n2\n", max_rows=1)

    def test_pk_h_units_and_redshift_are_not_interchanged(self):
        table = self.helper.parse_table(
            b"# Matter power spectrum P(k) at redshift z=0.38\n"
            b"# 1:k (h/Mpc) 2:P (Mpc/h)^3\n2 8\n4 1\n")
        result = self.helper.matter_spectrum(table, .5, .38)
        self.assertEqual(result["k_1_Mpc"], [1., 2.])
        self.assertEqual(result["P_Mpc3"], [64., 8.])
        with self.assertRaises(ValueError):
            self.helper.matter_spectrum(table, .5, .51)
        table["rows"].reverse()
        with self.assertRaises(ValueError):
            self.helper.matter_spectrum(table, .5, .38)

    def test_background_hubble_units_brackets_and_no_extrapolation(self):
        table = self.helper.parse_table(
            b"# 1:z 2:H [1/Mpc] 3:comov. dist. 4:ang.diam.dist. 5:lum. dist.\n"
            b"1 0.0002 200 100 400\n0 0.0001 0 0 0\n")
        exact = self.helper.background_at(table, 0.)
        self.assertEqual(exact["H_km_s_Mpc"], .0001 * 299792.458)
        self.assertEqual(exact["bracket_z"], [0., 0.])
        between = self.helper.background_at(table, .5)
        self.assertEqual(between["DM_Mpc"], 100.)
        self.assertEqual(between["bracket_z"], [0., 1.])
        self.assertIsNone(between["interpolation_error_bound"])
        with self.assertRaises(ValueError):
            self.helper.background_at(table, 1.01)

    def test_drag_is_predicted_drag_with_explicit_printing_precision(self):
        report = (b" -> recombination (maximum of visibility function) at z = 1100.000000\n"
                  b"    with comoving sound horizon = 144.000000 Mpc\n"
                  b" -> baryon drag stops at z = 1060.000000\n"
                  b"    corresponding to conformal time = 280.000000 Mpc\n"
                  b"    with comoving sound horizon rs = 147.000000 Mpc\n")
        result = self.helper.predicted_drag(report)
        self.assertEqual(result["r_drag_Mpc"], 147.)
        self.assertEqual(result["z_drag"], 1060.)
        self.assertEqual(result["rounding_half_width_Mpc"], 5e-7)
        self.assertFalse(result["internal_unrounded_value_available"])
        with self.assertRaises(ValueError):
            self.helper.predicted_drag(report + report)

    def test_precision_difference_handles_signed_te_without_certifying_error(self):
        first = self.helper.cmb_spectra(self.cmb(), tcmb=1e-6, lmax=3)
        changed = self.cmb()
        changed["rows"][0][3] = -8.
        second = self.helper.cmb_spectra(changed, tcmb=1e-6, lmax=3)
        result = self.helper.precision_difference(first, second)
        self.assertEqual(result["tt"]["max_abs_Dl_uK2"], 0.)
        self.assertEqual(result["te"]["max_abs_Dl_uK2"], 1.)
        self.assertEqual(result["te"]["max_abs_over_global_spectrum_scale"], 1 / 9)
        self.assertIsNone(result["te"]["certified_error_bound"])



class FullLCDMRunAdmission(unittest.TestCase):
    """Small input/fresh-output fixtures; no child process is launched."""

    def setUp(self):
        theory_helper = full_helper("theory")
        with patch.dict(sys.modules, {"theory": theory_helper}):
            self.helper = full_helper("run")
        # Public source metadata only; runtime/array inputs are not needed.
        self.config = json.loads((FULL_FOLDER / "reference.json").read_text())

    def test_runtime_design_refuses_type_aliases_physics_and_case_reordering(self):
        self.helper.validate_reference(self.config)
        for change in (
            lambda q: q["cases"][0].update(precision=0),
            lambda q: q["cases"][0].update(id="precision", precision=True),
            lambda q: q["parameters"].update(N_ur="3.044"),
            lambda q: q["limits"].update(case_cpu_seconds=True),
            lambda q: q["class"]["supplied_point"].update(sha256="0" * 64)):
            bad = json.loads(json.dumps(self.config))
            change(bad)
            with self.assertRaises(ValueError):
                self.helper.validate_reference(bad)

    def test_consumed_file_hash_and_fresh_output_refuse_changed_existing_bytes(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "fixture.txt"
            raw = b"small non-scientific fixture\n"
            pin = self.helper.write_new(path, raw)
            self.assertEqual(pin["sha256"], hashlib.sha256(raw).hexdigest())
            consumed, _ = self.helper.file_bytes(path, 1024, pin)
            self.assertEqual(consumed, raw)
            with self.assertRaises(FileExistsError):
                self.helper.write_new(path, b"replacement")
            path.chmod(0o600)
            path.write_bytes(raw + b"changed")
            with self.assertRaises(ValueError):
                self.helper.file_bytes(path, 1024, pin)
