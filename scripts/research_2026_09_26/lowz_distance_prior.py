"""Exact low-z distance-shape check, not a retraining or observed-bias estimate."""
from pathlib import Path
import hashlib
import json
import tarfile

import numpy as np
import pandas as pd
from scipy.integrate import quad
from astropy.cosmology import FlatLambdaCDM

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "runs/research_2026_09_26/lowz_distance_prior"
ARCHIVE = ROOT / "sources/updates/2026-09-20-ztf/originals/official-bayesn--08eef9e54188f601506ef9f7f47fc65f9f52a946.tar.gz"


def dimensionless_luminosity_distance(z, omega):
    return (1 + z) * quad(lambda x: 1 / np.sqrt(omega * (1 + x)**3 + 1 - omega),
                          0, z, epsabs=1e-13, epsrel=1e-13)[0]


def q_distance(z, q):
    return (1 + z) * (np.log1p(z) if q == 0 else -np.expm1(-q * np.log1p(z)) / q)


def main():
    OUT.mkdir(exist_ok=False)
    protocol = {
        "status": "Analytic follow-up specified before these numeric outcomes; not a statistical test.",
        "question": "Does an informative distance prior below z=.08 depend exactly only on H0?",
        "reference": "Flat LCDM Omega_m=.3, plus .28 archived code-default sensitivity; same arbitrary H0.",
        "alternatives": "Flat constant q=0 and q=.5 over the low-redshift interval.",
        "grid": [.005, .01, .02, .04, .06, .079, .08],
        "anchor": "Subtract alternative-minus-reference difference at z=.02 to remove one gray intercept.",
        "scale_comparison": "5/ln(10)*(150/300000)/z, peculiar-velocity-only scale, not total prior width or significance.",
        "boundary": "z=.08 is a limiting comparison; paper's informative branch is strictly below .08.",
    }
    (OUT / "protocol.json").write_text(json.dumps(protocol, indent=2) + "\n")
    with tarfile.open(ARCHIVE) as archive:
        name = "bayesn-bayesn-08eef9e/bayesn/bayesn_model.py"
        source = archive.extractfile(name).read()
    text = source.decode()
    assert 'fiducial_cosmology={"H0": 73.24, "Om0": 0.28}' in text
    assert '150 / 3e5' in text and 'Ds_err = jnp.sqrt(muhat_err * muhat_err + sigma0 * sigma0)' in text
    selected = []
    for i, line in enumerate(text.splitlines(), 1):
        if any(s in line for s in ['fiducial_cosmology={', '150 / 3e5', 'muhat_err = 5 /',
                                   'Ds_err = jnp.sqrt(muhat_err * muhat_err + sigma0 * sigma0)']):
            selected.append({"archive_member": name, "line": i, "source": line})
    rows = []
    verification = []
    for omega in [.3, .28]:
        independent = FlatLambdaCDM(H0=70., Om0=omega, Tcmb0=0)
        anchor = dimensionless_luminosity_distance(.02, omega)
        for q in [0., .5]:
            anchor_delta = 5 * np.log10(q_distance(.02, q) / anchor)
            for z in protocol["grid"]:
                fiducial = dimensionless_luminosity_distance(z, omega)
                alt = q_distance(z, q)
                # c/H0 cancels in the ratio; verify the fiducial with Astropy.
                other = independent.luminosity_distance(z).value / independent.hubble_distance.value
                error = float(abs(other - fiducial) / fiducial)
                assert error < 1e-10
                if q == .5:
                    analytic = 2 * (1 + z - np.sqrt(1 + z))
                    assert np.allclose(alt, analytic, rtol=1e-12, atol=1e-14)
                delta = 5 * np.log10(alt / fiducial)
                rows.append({"reference_Omega_m": omega, "alternative_q": q, "z": z,
                             "delta_mu_same_H0_mag": float(delta),
                             "delta_mu_after_z002_anchor_mag": float(delta - anchor_delta),
                             "peculiar_velocity_only_scale_mag": float(5 / np.log(10) * .0005 / z)})
                verification.append(error)
    table = pd.DataFrame(rows)
    table.to_csv(OUT / "distance-shape.csv", index=False)
    (OUT / "source-checks.json").write_text(json.dumps(selected, indent=2) + "\n")
    result = {"question_answer": "No: low-redshift cosmological shape dependence is suppressed by z but not zero.",
              "primary_near_threshold": table[(table.reference_Omega_m == .3) & (table.z == .079)].to_dict("records"),
              "source_default_near_threshold": table[(table.reference_Omega_m == .28) & (table.z == .079)].to_dict("records"),
              "max_relative_Astropy_verification_error": max(verification),
              "interpretation": "A fixed gray calibration cannot absorb all redshift-dependent prior changes. This only motivates a training-prior sensitivity.",
              "limitations": ["No G26 training rerun, matched training-redshift sample or inferred calibration/dust response.",
                              "The 150km/s value is only a scale comparison. Spectroscopic errors and latent intrinsic scatter also matter; inspected code combines sigma0 with this term in a training Ds prior.",
                              "Archived default code does not establish the exact G26 training execution or hybrid-prior implementation.",
                              "q=0 and q=.5 are declared contrasting kinematic histories, not estimates of the observed universe."]}
    (OUT / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
    manifest = {"inputs_sha256": {str(Path(__file__).relative_to(ROOT)): sha(__file__),
                                 str(ARCHIVE.relative_to(ROOT)): sha(ARCHIVE)},
                "inspected_source_sha256": hashlib.sha256(source).hexdigest(),
                "outputs_sha256": {str(p.relative_to(ROOT)): sha(p) for p in OUT.iterdir() if p.is_file()}}
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
