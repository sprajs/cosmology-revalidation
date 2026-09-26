#!/usr/bin/env python3
"""Independent invariants and numerical probes for matrix_audit outputs."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from astropy.cosmology import Flatw0waCDM
import matrix_audit as audit


def main():
    manifest = json.loads((audit.OUT/"manifest.json").read_text())
    assert manifest["script_sha256"] == audit.sha(Path(audit.__file__))
    for path, sha in {**manifest["inputs"], **manifest["outputs"]}.items():
        assert hashlib.sha256((audit.ROOT/path).read_bytes()).hexdigest() == sha, path
    result = {"manifest_input_and_output_hashes_match":True,
              "audit_manifest_sha256":audit.sha(audit.OUT/"manifest.json"),
              "verification_script_sha256":audit.sha(Path(__file__)),"releases":{}}
    grouping_path=audit.ROOT/'runs/salt_dust_audit/matrix_grouping/results.json'
    grouping=json.loads(grouping_path.read_text())
    for path,sha in grouping['input_sha256'].items():
        assert audit.sha(audit.ROOT/path)==sha,path
    result['independent_grouping_result_sha256']=audit.sha(grouping_path)
    for release,cfg in audit.RELEASES.items():
        out = audit.OUT/release
        base = audit.ROOT/cfg["root"]
        d = pd.read_csv(out/"rows_and_correction_jacobian.csv",dtype={"CID":str})
        modes = pd.read_csv(out/"systematic_spectral_modes.csv.gz",dtype={"CID":str})
        op = pd.read_csv(out/"cosmology_response_operator.csv",dtype={"CID":str})
        summary = json.loads((out/"summary.json").read_text())
        g=summary['grouping_adjusted_additivity']
        assert g['calibration_minus_CALSPEC_positive_rank']==9
        assert g['calibration_minus_CALSPEC_negative_rank']==0
        assert abs(g['relative_stat_whitened_frobenius']-grouping[release]['residual_relative_total_sys_stat_whitened_frobenius'])<1e-12
        if release=='Dovekie':
            assert g['relative_stat_whitened_frobenius']<1e-6
        else:
            rem=pd.read_csv(out/'unassigned_grouping_remainder_modes.csv',dtype={'CID':str})
            assert rem[['CID','IDSURVEY']].equals(d[['CID','IDSURVEY']])
            assert rem.mode_0based.nunique()==1 and (rem.eigenvalue_sign==1).all()
            assert g['remainder_resolved_positive_rank_at_1e_minus5_total_max']==1
        assert np.array_equal(d.row_index_0based,np.arange(len(d)))
        assert d[["CID","IDSURVEY"]].equals(op[["CID","IDSURVEY"]])
        assert summary["stat"]["min_eigenvalue_mag2"] > 0
        assert summary["total"]["min_eigenvalue_mag2"] > 0
        assert (modes.loc[modes.component=="VPEC","eigenvalue_sign"] < 0).any()
        # Both signal and manifest preserve source ordering; every mode carries IDs.
        for (component,k),g in modes.groupby(["component","mode_0based"],sort=False):
            assert g[["CID","IDSURVEY"]].reset_index(drop=True).equals(d[["CID","IDSURVEY"]])

        # Cross-check our quadrature and parameter derivatives against Astropy,
        # using the release's (1+zHEL) observer-frame luminosity-distance factor.
        ix = np.unique(np.linspace(0,len(d)-1,7,dtype=int).tolist()+[int(np.argmax(d.zHD))])
        z,zh = d.zHD.to_numpy()[ix],d.zHEL.to_numpy()[ix]
        def astropy_mu(Om0=.3,w0=-1,wa=0):
            cosmo = Flatw0waCDM(H0=70,Om0=Om0,w0=w0,wa=wa,Tcmb0=0)
            return cosmo.distmod(z).value+5*np.log10((1+zh)/(1+z))
        astropy_error = np.max(abs(astropy_mu()-audit.mu(z,zh)))
        assert astropy_error < 1e-9
        derivative_errors = {}
        for p,col,v in [("Om0","dmu_dom",.3),("w0","dmu_dw0",-1.),("wa","dmu_dwa",0.)]:
            derivative = (astropy_mu(**{p:v+1e-4})-astropy_mu(**{p:v-1e-4}))/2e-4
            derivative_errors[p] = float(np.max(abs(derivative-d[col].to_numpy()[ix])))
            assert derivative_errors[p] < 1e-6

        # Replay a small actual nonlinear cosmological fit around the stated
        # fiducial to check sign and normalization of A for a beta perturbation.
        p = np.diag(1/d.stat_sigma_mag.to_numpy()**2)
        one = np.ones(len(d)); po = p@one; den = one@po
        baseline = audit.mu(d.zHD.to_numpy(),d.zHEL.to_numpy())
        perturbation = d.dmu_dbeta.to_numpy()
        eps = 1e-3
        def fit(shift):
            y = baseline+shift*perturbation
            def objective(om):
                residual = y-audit.mu(d.zHD.to_numpy(),d.zHEL.to_numpy(),om=om)
                residual -= (po@residual)/den
                return float(residual@p@residual)
            return minimize_scalar(objective,bounds=(.29,.31),method="bounded",options={"xatol":1e-14}).x
        nonlinear_response = (fit(eps)-fit(-eps))/(2*eps)
        linear_response = op.stat_flat_LCDM__Omega_m.to_numpy()@perturbation
        response_difference = nonlinear_response-linear_response
        assert abs(response_difference) < max(1e-5,abs(linear_response)*1e-3)

        # Directly reconstruct two released MW matrices from signed spectral
        # modes. This verifies unpacking, subtraction, row order and units.
        suffix = ".txt.gz" if release=="original" else ".npz"
        cstat = np.diag(d.stat_sigma_mag.to_numpy()**2)
        sig = d.stat_sigma_mag.to_numpy()
        reconstruction = {}
        for component in (["MWEBV","MWCOLORLAW"] if release=="original" else ["MWEBV","COLORLAW"]):
            cs = audit.load_matrix(base/"SingleSYS_CovMatrix"/(component+suffix),len(d))
            if release=="Dovekie":
                cs-=cstat
            g = modes.loc[modes.component==component]
            pivot = g.pivot(index="row_index_0based",columns="mode_0based",values="loading_mag_arbitrary_sign")
            signs = g.groupby("mode_0based").eigenvalue_sign.first().to_numpy()
            b = pivot.to_numpy()
            residual = (cs-(b*signs)@b.T)/sig[:,None]/sig[None,:]
            relative_error = np.linalg.norm(residual)/np.linalg.norm(cs/sig[:,None]/sig[None,:])
            reconstruction[component] = float(relative_error)
            assert relative_error < 1e-4
        result["releases"][release] = {
            "grouping_agrees_with_independent_calculation":True,
            "all_mode_rows_match_HD":True,"VPEC_negative_modes_retained":True,
            "astropy_distance_max_error_mag":float(astropy_error),
            "astropy_derivative_max_errors":derivative_errors,
            "finite_nonlinear_beta_response_Omega_m":float(nonlinear_response),
            "linear_beta_response_Omega_m":float(linear_response),
            "beta_response_difference":float(response_difference),
            "MW_stat_whitened_factor_reconstruction_relative_error":reconstruction}
        if release == "original":
            vpec = audit.load_matrix(base/"SingleSYS_CovMatrix/VPEC.txt.gz",len(d))
            diagonal = vpec.diagonal()
            violation = vpec*vpec-np.outer(diagonal,diagonal)
            i,j = np.unravel_index(np.argmax(violation),violation.shape)
            assert violation[i,j] > 0
            result["releases"][release]["VPEC_serialization_and_PSD_witness"] = {
                "maximum_distance_from_1e_minus5_mag2_grid":float(np.max(abs(vpec*1e5-np.round(vpec*1e5)))/1e5),
                "row_i":int(i),"CID_i":str(d.CID.iloc[i]),"row_j":int(j),"CID_j":str(d.CID.iloc[j]),
                "Cii_mag2":float(vpec[i,i]),"Cij_mag2":float(vpec[i,j]),"Cjj_mag2":float(vpec[j,j]),
                "Cauchy_Schwarz_violation_mag4":float(violation[i,j]),
                "interpretation":"This released single component is not exactly PSD; compatible with coarse serialization, not proof of a physical defect or invalid total covariance."}
    (audit.OUT/"verification.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__ == "__main__":
    main()
