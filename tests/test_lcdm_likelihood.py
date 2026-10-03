"""Official-likelihood transport controls with stdlib mocks, never native calls."""
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import struct
import sys
import tempfile
import time
import types
import unittest
from unittest.mock import patch


SOURCE = Path(__file__).resolve().parents[1] / "experiments/lcdm-reference/likelihood.py"
METADATA = ("free_calib str A_planck\nlmin int 2\nnell int 28\n"
            "nstepsEE int 3000\nlmax int 29\nstepEE float 0.0001\n"
            "lkl_type str simall\npipeid str simall_EE_BB_TE\nunit int 1\n").encode()


class NativeFixtureError(ValueError):
    record = {"operation": "clik_compute", "errors": [{"code": -1234, "text": "fixture"}]}


class LikelihoodTransport(unittest.TestCase):
    def setUp(self):
        spec = importlib.util.spec_from_file_location("likelihood_fixture", SOURCE)
        self.helper = importlib.util.module_from_spec(spec)
        with patch.dict(sys.modules, {"clik": None, "numpy": None}):
            spec.loader.exec_module(self.helper)
        self.folder = tempfile.TemporaryDirectory()
        self.addCleanup(self.folder.cleanup)
        self.root = Path(self.folder.name)
        self.calls, self.constructed = [], []
        self.results = {"commander_TT": [-3.], "simall_EE": [-2.], "plik_lite_TTTEEE": [-5.]}
        self.axes = {"commander_TT": (29, -1, -1, -1, -1, -1),
                     "simall_EE": (-1, 29, -1, -1, -1, -1),
                     "plik_lite_TTTEEE": (3, 3, -1, 3, -1, -1)}
        for key in self.results:
            (self.root / self.helper.PRODUCTS[key]).mkdir(parents=True)
        self.metadata_path = self.root / self.helper.PRODUCTS["simall_EE"] / "clik/lkl_0/_mdb"
        self.metadata_path.parent.mkdir(parents=True)
        self.metadata_path.write_bytes(METADATA)

    def config(self):
        return {"backend": "clik", "highl": "plik_lite_TTTEEE", "plc_root": str(self.root),
                "nuisance": {"A_planck": 1.}, "calibration_prior": {"convention": "relative_penalty"}}

    def owner(self):
        fixture = self

        class Object:
            def __init__(self, key):
                self.key = key

            def get_lmax(self):
                return fixture.axes[self.key]

            def get_extra_parameter_names(self):
                return ("A_planck",)

            def __call__(self, vector):
                fixture.calls.append((self.key, list(vector)))
                result = fixture.results[self.key]
                if isinstance(result, Exception):
                    raise result
                return result

        def construct(path):
            key = next(k for k in self.results if path == str(self.root / self.helper.PRODUCTS[k]))
            self.constructed.append(key)
            return Object(key)

        clik = types.SimpleNamespace(version=lambda: "fixture-only", clik=construct)
        with patch.dict(sys.modules, {"clik": clik}):
            return self.helper.prepare(self.config())

    def spectra(self):
        return {"TT": [0., 0.] + [.5] * 28,
                "EE": [0., 0.] + [.0001] * 28,
                "TE": [0., 0.] + [-.01] * 28}

    def evaluate(self, owner, spectra):
        return owner.evaluate(spectra, {"A_planck": 1.},
                              self.helper.calibration_prior(1., "relative_penalty"))

    def test_packing_is_cl_tt_ee_te_then_runtime_nuisance_order(self):
        vector = self.helper.build_clik_vector(
            {"TT": [10., 11., 12.], "EE": [20., 21.], "TE": [-4., -5., -6.], "EB": [60.]},
            {"second": 8., "first": 7.}, (2, 1, -1, 2, -1, 0), ("second", "first"))
        self.assertEqual(vector, [10., 11., 12., 20., 21., -4., -5., -6., 60., 8., 7.])
        for invalid in (True, math.nan, math.inf, -math.inf):
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                self.helper.build_clik_vector({"TT": [0., 0., invalid]}, {},
                                              (2, -1, -1, -1, -1, -1), ())
        for axes in ((True, -1, -1, -1, -1, -1), (2., -1, -1, -1, -1, -1)):
            with self.assertRaises(ValueError):
                self.helper.build_clik_vector({"TT": [0., 0., 1.]}, {}, axes, ())

    def test_metadata_identity_and_native_float32_step_are_retained(self):
        owner = self.owner()
        support = owner.contracts["simall_EE"]["support_metadata"]
        self.assertEqual(support["stepEE_native_float32"], 9.999999747378752e-05)
        self.assertEqual(support["metadata_identity"], {
            "path": str(self.metadata_path), "bytes": 155,
            "sha256": hashlib.sha256(METADATA).hexdigest()})
        self.assertEqual(support["nstepsEE"], 3000)

    def test_malformed_changed_duplicate_oversized_metadata_refuses_before_simall_init(self):
        for raw in (METADATA.replace(b"nell int 28", b"nell int 27"),
                    METADATA + b"unit int 1\n", METADATA + b"unexpected int 1\n",
                    METADATA.replace(b"stepEE float", b"stepEE str"), b"x" * 4097):
            with self.subTest(raw=raw[:30]):
                self.metadata_path.write_bytes(raw)
                self.constructed.clear()
                with self.assertRaises(self.helper.LikelihoodRefusal) as caught:
                    self.owner()
                self.assertEqual(self.constructed, ["commander_TT"])
                failed = caught.exception.record["components"][-1]
                self.assertEqual(failed["id"], "simall_EE")
                self.assertIsNone(failed["native_error"])
                self.assertEqual(failed["support_metadata_error"]["phase"], "simall-support-metadata")
                if len(raw) <= 4096:
                    self.assertEqual(failed["support_metadata_error"]["metadata_identity"]["sha256"],
                                     hashlib.sha256(raw).hexdigest())
                self.assertEqual(self.calls, [])

    def test_metadata_symlink_and_fifo_refuse_without_following_or_blocking(self):
        target = self.root / "target"
        target.write_bytes(METADATA)
        self.metadata_path.unlink()
        self.metadata_path.symlink_to(target)
        with self.assertRaises(self.helper.LikelihoodRefusal):
            self.owner()
        self.metadata_path.unlink()
        import os
        os.mkfifo(self.metadata_path)
        with self.assertRaises(self.helper.LikelihoodRefusal):
            self.owner()
        self.assertEqual(target.read_bytes(), METADATA)
        self.assertEqual(self.calls, [])

    def test_first_route_fixed_calibration_refuses_before_initialization_or_compute(self):
        for nuisance in ({}, {"A_planck": True}, {"A_planck": 1.01},
                         {"A_planck": math.nan}, {"A_planck": 1., "extra": 0.}):
            config = self.config(); config["nuisance"] = nuisance
            with patch.object(self.helper, "PlanckPrimary") as constructor:
                with self.assertRaises(ValueError):
                    self.helper.prepare(config)
                constructor.assert_not_called()
        owner = self.owner()
        with self.assertRaises(self.helper.LikelihoodRefusal) as caught:
            owner.evaluate(self.spectra(), {"A_planck": 1.01},
                           self.helper.calibration_prior(1.01, "relative_penalty"))
        self.assertEqual(self.calls, [])
        self.assertTrue(all(r["status"] == "not_started"
                            for r in caught.exception.record["component_attempts"]))
        config = self.config(); config["calibration_prior"] = {"convention": "normalized_density"}
        with self.assertRaises(ValueError):
            self.helper.validate_config(config)

    def test_negative_nonfinite_and_overflow_ee_refuse_before_simall_compute(self):
        owner = self.owner()
        for value in (-.0001, math.nan, math.inf, -math.inf, 1e308):
            with self.subTest(value=value):
                spectra = self.spectra(); spectra["EE"][2] = value
                self.calls.clear()
                with self.assertRaises(self.helper.LikelihoodRefusal) as caught:
                    self.evaluate(owner, spectra)
                self.assertEqual([k for k, _ in self.calls], ["commander_TT"])
                record = caught.exception.record
                self.assertEqual(record["components"], {"commander_TT": -3.})
                self.assertIsNone(record["combined"])
                failed = record["component_attempts"][1]
                self.assertEqual(failed["support_check"]["ell"], 2)
                self.assertEqual(failed["support_check"]["metadata_identity"],
                                 owner.contracts["simall_EE"]["support_metadata"]["metadata_identity"])
                self.assertIsNone(failed["input_vector"])
                self.assertIsNone(failed["native_error"])
                self.assertEqual(record["component_attempts"][2]["status"], "not_started")
                json.dumps(record, allow_nan=False)

    def test_exact_native_index_equality_3000_is_rejected(self):
        owner = self.owner()
        step = owner.contracts["simall_EE"]["support_metadata"]["stepEE_native_float32"]
        witness = None
        # Search a bounded neighbourhood for an exactly representable C-operation
        # boundary witness, rather than assuming inverted floating arithmetic.
        for ell in range(2, 30):
            cl = 3000. * step * 2. * math.pi / (ell * (ell + 1))
            for direction in (-math.inf, math.inf):
                value = cl
                for _ in range(16):
                    index = value * ell * (ell + 1) / 2. / math.pi / step
                    if index == 3000.:
                        witness = (ell, value); break
                    value = math.nextafter(value, direction)
                if witness:
                    break
            if witness:
                break
        self.assertIsNotNone(witness, "bounded exact equality fixture must exist")
        ell, value = witness
        spectra = self.spectra(); spectra["EE"][ell] = value
        with self.assertRaises(self.helper.LikelihoodRefusal) as caught:
            self.evaluate(owner, spectra)
        self.assertEqual([k for k, _ in self.calls], ["commander_TT"])
        guard = caught.exception.record["component_attempts"][1]["support_check"]
        self.assertEqual(guard["index_argument"]["raw_value"], 3000.)

    def test_zero_and_inside_support_return_actual_three_terms_without_rescaling(self):
        owner = self.owner()
        spectra = self.spectra(); spectra["EE"][2] = 0.
        result = self.evaluate(owner, spectra)
        self.assertEqual(result["loglike"], -10.)
        self.assertEqual(result["logtarget"], -10.)
        self.assertIsNone(result["logposterior"])
        self.assertEqual([k for k, _ in self.calls], list(owner.order))
        self.assertEqual(self.calls[1][1], spectra["EE"] + [1.])
        self.assertEqual(self.calls[2][1], spectra["TT"][:4] + spectra["EE"][:4]
                         + spectra["TE"][:4] + [1.])
        guard = result["component_attempts"][1]["support_check"]
        self.assertEqual(guard["ells_checked"], 28)
        self.assertEqual(guard["minimum_index_argument"], 0.)
        self.assertEqual(result["calibration_prior"]["applications"], 1)

    def test_simall_runtime_axes_mismatch_retains_metadata_prefix(self):
        self.axes["simall_EE"] = (-1, 30, -1, -1, -1, -1)
        with self.assertRaises(self.helper.LikelihoodRefusal) as caught:
            self.owner()
        failed = caught.exception.record["components"][-1]
        self.assertEqual(failed["support_metadata"]["metadata_identity"]["bytes"], 155)
        self.assertEqual(failed["status"], "refused")
        self.assertEqual(self.calls, [])

    def test_late_native_error_sentinel_and_nonfinite_return_preserve_earned_prefix(self):
        owner = self.owner()
        for result in (NativeFixtureError("fixture native error"), [-1e30], [math.nan]):
            self.results["plik_lite_TTTEEE"] = result
            with self.subTest(result=result), self.assertRaises(self.helper.LikelihoodRefusal) as caught:
                self.evaluate(owner, self.spectra())
            record = caught.exception.record
            self.assertEqual(record["components"], {"commander_TT": -3., "simall_EE": -2.})
            self.assertIsNone(record["combined"])
            failed = record["component_attempts"][2]
            self.assertEqual(failed["input_vector"]["scalars"], 13)
            if isinstance(result, Exception):
                self.assertEqual(failed["native_error"], NativeFixtureError.record)
                self.assertIsNone(failed["raw_observation"])
            else:
                self.assertEqual(failed["raw_observation"]["binary64_le_hex"],
                                 struct.pack("<d", result[0]).hex())
            json.dumps(record, allow_nan=False)

    def test_declared_prior_constant_is_once_and_never_a_posterior_alias(self):
        relative = self.helper.calibration_prior(1., "relative_penalty")
        normalized = self.helper.calibration_prior(1., "normalized_density")
        a = self.helper.combine_terms({"T": -3., "E": -2., "high": -5.}, relative)
        b = self.helper.combine_terms({"T": -3., "E": -2., "high": -5.}, normalized)
        constant = -math.log(.0025 * math.sqrt(2. * math.pi))
        self.assertEqual(a["loglike"], b["loglike"])
        self.assertEqual(a["logprior_term"], 0.)
        self.assertAlmostEqual(b["logtarget"] - a["logtarget"], constant)
        self.assertIsNone(a["logposterior"])
        self.assertIsNone(b["posterior_normalization"])
        self.assertIn("unverified", relative["source"]["A_planck_to_y_P_mapping"])
        altered = dict(relative); altered["applications"] = 3
        with self.assertRaises(ValueError):
            self.helper.combine_terms({"T": -3.}, altered)

    def test_bridge_keeps_cl_padding_and_refusal_file_has_partial_components(self):
        owner = self.owner()
        spectra = self.spectra()
        cmb = {"lensed": True, "ell": list(range(2, 30)), "spectra": {
            name.lower(): {"Cl_uK2": spectra[name][2:], "Dl_uK2": [999.] * 28}
            for name in ("TT", "EE", "TE")}}
        success = self.root / "success"; success.mkdir()
        result = self.helper.evaluate(cmb, self.config(), success, owner=owner)
        self.assertEqual(result["loglike"], -10.)
        self.assertEqual(self.calls[0][1][:2], [0., 0.])
        self.assertEqual(self.calls[2][1][8:12], spectra["TE"][:4])
        refused = self.root / "refused"; refused.mkdir()
        cmb["spectra"]["ee"]["Cl_uK2"][0] = -.0001
        self.calls.clear()
        with self.assertRaises(self.helper.LikelihoodRefusal):
            self.helper.evaluate(cmb, self.config(), refused, owner=owner)
        saved = json.loads((refused / "likelihood-failure.json").read_text())
        self.assertEqual(saved["components"], {"commander_TT": -3.})
        self.assertEqual([k for k, _ in self.calls], ["commander_TT"])
        self.assertFalse((refused / "likelihood.json").exists())


class ReleasedProductCustody(unittest.TestCase):
    def setUp(self):
        source = SOURCE.with_name("score.py")
        spec = importlib.util.spec_from_file_location("score_fixture", source)
        self.scorer = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.scorer)
        self.folder = tempfile.TemporaryDirectory()
        self.addCleanup(self.folder.cleanup)
        self.root = Path(self.folder.name)
        self.roots, self.pins = [], []
        for _, relative in self.scorer.PRODUCTS:
            root = self.root / relative
            path = root / "clik/lkl_0/_mdb"
            path.parent.mkdir(parents=True)
            path.write_bytes(b"fixture metadata\n")
            self.roots.append(root)
            self.pins.append({"path": str(path), "bytes": path.stat().st_size,
                              "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})

    def tree(self):
        return self.scorer.product_tree(self.root, self.pins, time.monotonic() + 5.)

    def test_exact_tree_closes_files_and_nested_directories(self):
        result = self.tree()
        self.assertEqual(result["files"], sorted(pin["path"] for pin in self.pins))
        self.assertEqual(result["directories"], sorted(str(path) for root in self.roots
                         for path in (root, root / "clik", root / "clik/lkl_0")))

    def test_unlisted_native_key_file_or_empty_subtree_refuses(self):
        extra = self.roots[0] / "clik/default"
        extra.write_bytes(b"unlisted optional native key")
        with self.assertRaisesRegex(ValueError, "undeclared or nonregular"):
            self.tree()
        extra.unlink()
        extra.mkdir()
        with self.assertRaisesRegex(ValueError, "undeclared released directory"):
            self.tree()

    def test_missing_duplicate_outside_or_expired_inventory_refuses(self):
        path = Path(self.pins[0]["path"])
        path.unlink()
        with self.assertRaisesRegex(ValueError, "inventory differs"):
            self.tree()
        path.write_bytes(b"fixture metadata\n")
        for pin in (self.pins[0], {"path": str(self.root / "outside")},
                    {"path": str(self.roots[0] / "../outside")}):
            with self.subTest(pin=pin), self.assertRaises(ValueError):
                self.scorer.product_tree(self.root, self.pins + [pin], time.monotonic() + 5.)
        with self.assertRaisesRegex(ValueError, "deadline"):
            self.scorer.product_tree(self.root, self.pins, time.monotonic() - 1.)

    def test_symlink_fifo_and_symlink_root_refuse_without_opening(self):
        extra = self.roots[0] / "clik/default"
        extra.symlink_to(Path(self.pins[0]["path"]))
        with self.assertRaises(ValueError):
            self.tree()
        extra.unlink()
        os.mkfifo(extra)
        with self.assertRaises(ValueError):
            self.tree()
        extra.unlink()
        root = self.roots[0]
        saved = root.with_name(root.name + "-saved")
        root.rename(saved)
        root.symlink_to(saved, target_is_directory=True)
        with self.assertRaisesRegex(ValueError, "directory type"):
            self.tree()

    def test_selfcheck_criterion_cannot_be_loosened_or_silently_replaced(self):
        raw = b"".join(("Checking likelihood '" + str(self.root / relative) +
                        "' on test data. got -1 expected -1 (diff 1e-9)\n").encode()
                       for _, relative in self.scorer.PRODUCTS)
        result = self.scorer.selfchecks(raw, self.root, {"max_abs_printed_difference": 1e-6})
        self.assertTrue(all(row["passed"] for row in result))
        for criterion in ({"max_abs_printed_difference": 1.},
                          {"max_abs_printed_difference": 1e-5},
                          {"max_abs_printed_difference": 0.},
                          {"max_abs_printed_difference": math.nan}, {}):
            with self.subTest(criterion=criterion), self.assertRaisesRegex(ValueError, "unchanged"):
                self.scorer.selfchecks(raw, self.root, criterion)


if __name__ == "__main__":
    unittest.main()
