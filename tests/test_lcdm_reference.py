"""Reference-only analytic controls independent of reference self-refinement."""
from decimal import Decimal as D, localcontext
import importlib.util
from pathlib import Path
import unittest

FOLDER=Path(__file__).resolve().parents[1]/"experiments/lcdm-baseline"
spec=importlib.util.spec_from_file_location("lcdm_reference_controls",FOLDER/"reference.py")
reference=importlib.util.module_from_spec(spec)
spec.loader.exec_module(reference)


class ReferenceControls(unittest.TestCase):
    def test_gl8_polynomial_moments_through_degree15(self):
        with localcontext() as c:
            c.prec=60
            rule=reference.nodes()
            for degree in range(16):
                actual=reference.quadrature(lambda x:x**degree,D(1),1,rule)
                self.assertLess(abs(actual-D(1)/(degree+1)),D("1e-52"))

    def test_dense_ldlt_known_correlated_quadratic_and_determinant(self):
        with localcontext() as c:
            c.prec=60
            result=reference.ldlt([[D(4),D(1)],[D(1),D(9)]],[D(1),D(2)])
            self.assertLess(abs(result["quadratic"]-D(".6")),D("1e-55"))
            self.assertLess(abs(result["log_determinant"]-D(35).ln()),D("1e-55"))
            self.assertLess(abs(reference.pi()-D("3.141592653589793238462643383279502884197169399375105820974944")),D("1e-54"))

    def test_radiation_only_distance_and_supplied_ruler_analytic_limit(self):
        q={"model":{"h0_km_s_mpc":70.,"omega_m":0.,"omega_r":1.,"omega_b":0.,"omega_gamma":1.,"z_drag":9.,"drag_origin":"radiation-only analytic control"},"redshifts":[0.,.5,1.,2.],"rows":[{"z":.5,"observable":"DM_over_rs"} for _ in range(13)]}
        covariance=[[float(i==j) for j in range(13)] for i in range(13)]
        result=reference.evaluate(q,[1.]*13,covariance,128)
        with localcontext() as c:
            c.prec=60
            scale=D("299792.458")/70
            for row in result["background"]:
                z=row["z"]
                self.assertLess(abs(row["E"]-(1+z)**2),D("1e-50"))
                self.assertLess(abs(row["DM"]-scale*z/(1+z)),D("1e-45"))
                self.assertLess(abs(row["DL"]-scale*z),D("1e-45"))
            ratio=D(3).sqrt()*10/3
            for value in result["predictions"]:
                self.assertLess(abs(value-ratio),D("1e-48"))
