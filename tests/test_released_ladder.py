"""Meaningful adversarial lineage/measure/budget checks; no released assets in CI."""
import copy
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location("released_ladder",ROOT/"experiments/released-ladder/controller.py")
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)

class ReleasedLineageTests(unittest.TestCase):
    def setUp(self):self.manifest=module.load(ROOT/"experiments/released-ladder/lineage.json")

    def test_order_and_zero_width_measure_are_required(self):
        module.validate_lineage(self.manifest)
        changed=copy.deepcopy(self.manifest);changed["axes"][0],changed["axes"][1]=changed["axes"][1],changed["axes"][0]
        with self.assertRaisesRegex(ValueError,"ordered full47"):module.validate_lineage(changed)
        changed=copy.deepcopy(self.manifest);changed["released_mcmc_prior"]["fixed_coordinates"]=[]
        with self.assertRaisesRegex(ValueError,"zero-width"):module.validate_lineage(changed)
        changed=copy.deepcopy(self.manifest);changed["axes"][44]["physical_identity"]="SMC"
        with self.assertRaisesRegex(ValueError,"no admitted"):module.validate_lineage(changed)

    def test_unreviewed_budgets_source_target_and_perturbations_rejected(self):
        for field in module.FROZEN_BUDGETS:
            changed=copy.deepcopy(self.manifest);changed["target"]["reference_budgets"][field]*=100
            with self.assertRaisesRegex(ValueError,"comparison budgets"):module.validate_lineage(changed)
        changed=copy.deepcopy(self.manifest);changed["target"]["paper_rounded_coordinate"]["half_last_printed_digit"]*=100
        with self.assertRaisesRegex(ValueError,"rounded-coordinate"):module.validate_lineage(changed)
        changed=copy.deepcopy(self.manifest);changed["constraint_rows"]["sensitivity_control"]["delta_y"]*=2
        with self.assertRaisesRegex(ValueError,"sensitivity rows"):module.validate_lineage(changed)
        changed=copy.deepcopy(self.manifest);changed["constraint_rows"]["indices"].reverse()
        with self.assertRaisesRegex(ValueError,"sensitivity rows"):module.validate_lineage(changed)

    def test_oversized_source_refused_before_read(self):
        with tempfile.TemporaryDirectory() as td:
            Path(td,"file").write_bytes(b"abcd")
            source={"sources":[{"name":"file","bytes":3,"sha256":"0"*64}]}
            with patch.object(Path,"open",side_effect=AssertionError("must not allocate/read oversized source")):
                with self.assertRaisesRegex(ValueError,"byte count differs"):module.admitted_sources(Path(td),source)

    def test_failed_child_retains_after_integrity_and_original_failure(self):
        with tempfile.TemporaryDirectory() as td:
            args=SimpleNamespace(name="failure",sources=Path(td),engine_source=Path(td),sdk=Path(td),reference_python=None)
            with patch.object(module,"ROOT",Path(td)),patch.object(module,"admitted_sources",return_value={"test":b"abc"}),patch.object(module,"source_audit",return_value={}),patch.object(module,"fingerprint",side_effect=[{"sdk":"before"},{"sdk":"changed"}]),patch.object(module,"child",side_effect=ValueError("compiler failure")),patch.object(module.subprocess,"check_output",return_value="diagnostic"):
                self.assertEqual(module.execute(args),1)
            record=module.load(Path(td)/"results/released-ladder/failure/record.json")
            self.assertEqual(record["error"],"compiler failure")
            self.assertEqual(record["sdk_after"],{"sdk":"changed"})
            self.assertEqual(record["integrity_errors"],["SDK changed"])
            self.assertEqual(record["source_hashes_after"],record["source_hashes"])
            self.assertEqual((Path(td)/"results/released-ladder/failure/record.json").stat().st_mode & 0o222,0)

    def test_malformed_or_changed_source_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            source={"sources":[{"name":"file","bytes":3,"sha256":"0"*64}]}
            Path(td,"file").write_bytes(b"abc")
            with self.assertRaisesRegex(ValueError,"hash differ"):module.admitted_sources(Path(td),source)
        with self.assertRaisesRegex(ValueError,"primary profile"):
            module.fits_image(b"SIMPLE  = F".ljust(2880,b" "),(1,))

    def test_source_join_ambiguity_and_host_mismatch_rejected(self):
        # One source row per each of37 hosts; complete original axis-index mapping.
        n,p=37,47;x=[0.]*(n*p);table=[];hosts=[]
        for i in range(n):
            x[i*p+i]=1;x[i*p+41]=.25;x[i*p+43]=i
            hosts.append(f"h{i}")
            table.append(dict(host=f"h{i}",id=str(i),line=i+1,period_lower=.24,period_upper=.26,metal=i,metal_half=.005))
        result=module.host_join(x,table,hosts,n=n,p=p,stop=n)
        self.assertEqual([r["host"] for r in result],hosts)
        wrong=hosts.copy();wrong[0]="guessed"
        with self.assertRaisesRegex(ValueError,"unresolved/changed"):module.host_join(x,table,wrong,n=n,p=p,stop=n)
        ambiguous=table+[dict(table[0],host="another")]
        with self.assertRaisesRegex(ValueError,"unresolved/changed"):module.host_join(x,ambiguous,hosts,n=n,p=p,stop=n)

    def test_frozen_coefficient_budget_not_relaxed(self):
        native={"coefficients":[0.]*47,"quadratic":0.,"variance46":1.,"sensitivities":[]}
        ref={"rank":47,"algorithms":{name:{"coefficients":[0.]*47,"quadratic":0.,"variance46":1.} for name in ("LAPACK_gesdd_SVD","LAPACK_pivoted_QR")},"sensitivities":[]}
        native["coefficients"][46]=2e-8
        with self.assertRaisesRegex(ValueError,"frozen comparison failed"):module.compare(native,ref,self.manifest)

    def test_printed_period_intervals_use_actual_precision(self):
        rows=module.photometry_rows("A & 0 & 0 & 1 & 10.00 & 0 & 0 & 0 & 0 & -0.10 & HST")
        self.assertEqual(rows[0]["metal_half"],.005)
        self.assertLess(rows[0]["period_lower"],0.)
        self.assertGreater(rows[0]["period_upper"],0.)
        self.assertGreater(module.binary32_half_ulp(.25),0.)

if __name__=="__main__":unittest.main()
