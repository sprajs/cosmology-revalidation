"""Adversarial experiment admission/receipt controls; not numerical qualification."""
import copy
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
from decimal import Decimal

FOLDER=Path(__file__).resolve().parents[1]/"experiments/lcdm-baseline"
sys.path.insert(0,str(FOLDER))
spec=importlib.util.spec_from_file_location("lcdm_controller",FOLDER/"controller.py")
controller=importlib.util.module_from_spec(spec)
spec.loader.exec_module(controller)


class BaselineTests(unittest.TestCase):
    def setUp(self):
        self.q=controller.load(FOLDER/"request.json")

    def test_published_candidate_admission_and_full_model_refusal(self):
        packet,request,identities=controller.read_packet(FOLDER)
        self.assertEqual(packet["origin"]["revision"],"3e763bad692a13767e34e99af167623c3de4d0d8")
        self.assertEqual(identities["request"],"5d67253917a35031b69b85272fe381da8c946007da5615300a3dcf3381f12bed")
        self.assertEqual(request,FOLDER/"request.json")
        with tempfile.TemporaryDirectory() as t:
            folder=Path(t)/"lcdm-baseline";folder.mkdir()
            for name in ("experiment.json","request.json","candidate.json","README.md"):
                (folder/name).write_bytes((FOLDER/name).read_bytes())
            candidate=controller.load(folder/"candidate.json")
            candidate["readiness"]="blocked"
            candidate["blockers"]=["Full Planck massive-neutrino/thermal-drag/CMB closures unavailable"]
            (folder/"candidate.json").write_text(json.dumps(candidate))
            packet["origin"]["sha256"]=controller.sha256(folder/"candidate.json")
            (folder/"experiment.json").write_text(json.dumps(packet))
            with self.assertRaisesRegex(ValueError,"not ready"):
                controller.read_packet(folder)

    def test_supported_typed_request_and_full_model_rejection(self):
        controller.validate_request(self.q)
        for field in ("m_nu_ev","theta_mc","A_s","n_s","tau","omega_c"):
            q=copy.deepcopy(self.q);q["model"][field]=.06
            with self.assertRaises(ValueError):controller.validate_request(q)

    def test_bad_types_and_physical_subsets(self):
        for name,value in (("h0_km_s_mpc",True),("omega_gamma",float("nan")),("omega_b",1),("omega_r",-1),("z_drag",-1)):
            q=copy.deepcopy(self.q);q["model"][name]=value
            with self.assertRaises(ValueError):controller.validate_request(q)

    def test_axis_quotas_source_and_underbudget_refused(self):
        for change in (lambda q:q["rows"].pop(),lambda q:q["rows"].reverse(),
                       lambda q:q["rows"][0].update(observable="untyped"),
                       lambda q:q["budgets"].update(density_absolute=1e-4),
                       lambda q:q["resources"].update(reference_panels=[1,2]),
                       lambda q:q["model"].update(drag_origin="\nunknown")):
            q=copy.deepcopy(self.q);change(q)
            with self.assertRaises(ValueError):controller.validate_request(q)

    def test_packet_request_revision_timeout_and_assurance_must_agree(self):
        packet=controller.load(FOLDER/"experiment.json")
        controller.validate_admission(packet,self.q)
        for change in (lambda p:p["execution"].update(engine_revision="0"*40),
                       lambda p:p["execution"].update(timeout_seconds=1),
                       lambda p:p["execution"].update(assurance="qualified"),
                       lambda p:p["execution"].update(interface="cli")):
            bad=copy.deepcopy(packet);change(bad)
            with self.assertRaisesRegex(ValueError,"identity differs"):
                controller.validate_admission(bad,self.q)

    def test_changed_admitted_packet_or_request_retained_before_compile(self):
        for changed in ("experiment.json","request.json"):
            with tempfile.TemporaryDirectory() as t:
                root=Path(t);folder=root/"experiments/lcdm-baseline";folder.mkdir(parents=True)
                for name in ("experiment.json","request.json"):
                    (folder/name).write_bytes((FOLDER/name).read_bytes())
                admitted={"packet":controller.sha256(folder/"experiment.json"),"request":controller.sha256(folder/"request.json")}
                packet=controller.load(folder/"experiment.json")
                def tampered_admission(_):
                    (folder/changed).write_text((folder/changed).read_text()+" ")
                    return packet,folder/"request.json",admitted
                with patch.object(controller,"ROOT",root),patch.object(controller,"FOLDER",folder),patch.object(controller,"git",return_value="test"),patch.object(controller,"read_packet",side_effect=tampered_admission):
                    self.assertEqual(controller.execute(root,root,"mutation"),1)
                store=root/"results/lcdm-baseline/mutation"
                result=json.loads((store/"run.json").read_text())
                self.assertIn("changed after admission",result["error"])
                self.assertFalse((store/"compile.out").exists())

    def test_unique_exact_input_roles_and_paths_required(self):
        packet=controller.load(FOLDER/"experiment.json")
        for inputs in ([packet["inputs"][1]]*2, [packet["inputs"][0]]*2,
                       list(reversed(packet["inputs"])), packet["inputs"][:1]):
            bad=copy.deepcopy(packet);bad["inputs"]=inputs
            with patch.object(controller,"verify_inputs"):
                with self.assertRaisesRegex(ValueError,"exact unique"):
                    controller.scientific_inputs(bad,self.q)
        bad=copy.deepcopy(packet);bad["inputs"][0]["role"]="synthetic_control"
        with patch.object(controller,"verify_inputs"):
            with self.assertRaisesRegex(ValueError,"exact unique"):
                controller.scientific_inputs(bad,self.q)

    def test_consumed_bytes_are_hash_bound_not_only_prior_path_check(self):
        packet=controller.load(FOLDER/"experiment.json")
        class Stats:st_size=packet["inputs"][0]["bytes"]
        with patch.object(controller,"verify_inputs"),patch.object(Path,"stat",return_value=Stats()),patch.object(Path,"read_bytes",return_value=b"temporary substituted bytes"):
            with self.assertRaisesRegex(ValueError,"Consumed input"):
                controller.scientific_inputs(packet,self.q)

    def test_sdk_bad_hash_and_dirty_source_refused(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t);(p/"build-manifest.json").write_text('{}')
            with patch.object(controller,"git",side_effect=[self.q["engine_identity"]["revision"],""]):
                with self.assertRaisesRegex(ValueError,"manifest hash"):
                    controller.fingerprint(p,p,self.q)
            with patch.object(controller,"git",side_effect=[self.q["engine_identity"]["revision"]," M changed"]):
                with self.assertRaisesRegex(ValueError,"clean"):
                    controller.fingerprint(p,p,self.q)

    def test_deleted_even_unused_sdk_header_is_refused(self):
        with tempfile.TemporaryDirectory() as t:
            sdk=Path(t);folder=sdk/"include/irred";folder.mkdir(parents=True)
            (folder/"used.hpp").write_text("used")
            build={"sources":{"cpp/include/irred/used.hpp":controller.sha256(folder/"used.hpp"),
                              "cpp/include/irred/unused.hpp":"0"*64}}
            with self.assertRaisesRegex(ValueError,"inventory"):
                controller.verify_headers(sdk,build)

    def output(self):
        density={"status":"ok","predictions":[1.]*13,"quadratic":0.,"log_determinant":0.,"normalization":0.,"log_density":0.,"projection_estimate":1e-12,"callbacks":1}
        background=[{"status":"ok","z":z,"E":1.,"DM":0.,"DL":0.} for z in self.q["redshifts"]]
        return {"schema_version":1,"scientific_ids":dict(controller.SCIENTIFIC_IDS),"arithmetic":{"density_arithmetic_id":"F02/longdouble-cpu/v1","double_mantissa_bits":53,"long_double_mantissa_bits":64,"long_double_max_exponent":16384,"round_to_nearest":True},"producer_policy":dict(controller.PRODUCER),"density":[copy.deepcopy(density),copy.deepcopy(density)],"background":[copy.deepcopy(background),copy.deepcopy(background)],"callbacks":2}

    def test_status_missing_axis_and_nonfinite_outputs_refused(self):
        controller.check_outputs(self.output(),self.q)
        for change in (lambda o:o["density"][0].update(status="failed"),
                       lambda o:o["density"][0]["predictions"].pop(),
                       lambda o:o["density"][0].update(projection_estimate=None),
                       lambda o:o["background"][0][0].update(z=99),
                       lambda o:o["density"][0].update(log_density=float("inf")),
                       lambda o:o["scientific_ids"].update(physical_model="other-model"),
                       lambda o:o["arithmetic"].update(round_to_nearest=False)):
            o=self.output();change(o)
            with self.assertRaises((ValueError,TypeError)):controller.check_outputs(o,self.q)

    def test_independent_comparison_and_refinement_failure(self):
        fine={"background":[{"z":Decimal.from_float(float(z)),"E":Decimal(1),"DM":Decimal(0),"DL":Decimal(0)} for z in self.q["redshifts"]],"predictions":[Decimal(1)]*13,"quadratic":Decimal(0),"log_determinant":Decimal(0),"normalization":Decimal(0),"log_density":Decimal(0)}
        for change in (lambda r:r["predictions"].pop(),lambda r:r["background"].pop(),
                       lambda r:r["predictions"].append(Decimal(1)),
                       lambda r:r["background"][0].update(E=Decimal("NaN"))):
            bad=copy.deepcopy(fine);change(bad)
            with self.assertRaisesRegex(ValueError,"reference"):
                controller.comparisons(self.output(),bad,fine,self.q)
        bad=copy.deepcopy(fine);bad["predictions"][0]=Decimal(2)
        with self.assertRaisesRegex(ValueError,"refinement"):
            controller.comparisons(self.output(),bad,fine,self.q)
        o=self.output();o["density"][0]["predictions"][0]=2
        with self.assertRaisesRegex(ValueError,"Scientific comparison"):
            controller.comparisons(o,fine,fine,self.q)

    def test_dimensionless_twin_E_does_not_receive_Mpc_additive_slack(self):
        fine={"background":[{"z":Decimal.from_float(float(z)),"E":Decimal(1),"DM":Decimal(0),"DL":Decimal(0)} for z in self.q["redshifts"]],"predictions":[Decimal(1)]*13,"quadratic":Decimal(0),"log_determinant":Decimal(0),"normalization":Decimal(0),"log_density":Decimal(0)}
        output=self.output();output["background"][1][0]["E"]=1+5e-10
        with self.assertRaisesRegex(ValueError,r"H0_background\[0\].E"):
            controller.comparisons(output,fine,fine,self.q)

    def test_fresh_failed_attempt_immutable_and_never_overwritten(self):
        with tempfile.TemporaryDirectory() as t:
            with patch.object(controller,"ROOT",Path(t)),patch.object(controller,"git",return_value="test"),patch.object(controller,"read_packet",side_effect=ValueError("bad input hash")):
                self.assertEqual(controller.execute(Path(t),Path(t),"failure"),1)
                record=Path(t)/"results/lcdm-baseline/failure/run.json"
                result=json.loads(record.read_text())
                self.assertEqual(result["execution"],"failed")
                self.assertEqual(result["qualification"]["numerical"],"not_assessed")
                self.assertIn("bad input hash",result["error"])
                self.assertEqual(record.stat().st_mode&0o222,0)
                before=record.read_bytes()
                with self.assertRaises(FileExistsError):controller.execute(Path(t),Path(t),"failure")
                self.assertEqual(before,record.read_bytes())

    def test_committed_snapshot_preserves_bytes_and_rejects_dirty_replacement(self):
        import io,tarfile
        for changed in (False,True):
            with tempfile.TemporaryDirectory() as t:
                root=Path(t);store=root/"results";store.mkdir()
                (root/"contract.py").write_bytes(b"changed" if changed else b"committed source")
                def archive(command,*args,**kwargs):
                    with tarfile.open(store/"source.tar","w") as output:
                        item=tarfile.TarInfo("contract.py");item.size=len(b"committed source")
                        output.addfile(item,io.BytesIO(b"committed source"))
                with patch.object(controller,"ROOT",root),patch.object(controller,"bounded",side_effect=archive):
                    if changed:
                        with self.assertRaisesRegex(ValueError,"immutable committed"):
                            controller.snapshot_source(store,"reviewed-revision")
                    else:
                        destination=controller.snapshot_source(store,"reviewed-revision")
                        self.assertEqual((destination/"contract.py").read_bytes(),b"committed source")
                        self.assertEqual((destination/"contract.py").stat().st_mode&0o222,0)

    def test_timeout_and_output_resource_caps_preserve_failure(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)
            with self.assertRaises(Exception):controller.bounded([sys.executable,"-c","import time;time.sleep(2)"],p,"timeout",.05)
            self.assertTrue((p/"timeout.err").exists())
            with self.assertRaises(ValueError):controller.bounded([sys.executable,"-c","print('x'*100000)"],p,"overflow",2,64)
            self.assertLessEqual((p/"overflow.out").stat().st_size,64)
