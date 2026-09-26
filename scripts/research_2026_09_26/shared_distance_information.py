"""Design-only shared-mode distance information, without using residual values."""
from pathlib import Path
import hashlib
import json

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "runs/research_2026_09_26/astra_design"
OUT = ROOT / "runs/research_2026_09_26/shared_distance_information"


def main():
    OUT.mkdir(exist_ok=False)
    protocol = {"scope": "Posterior variance geometry in the same finite Gaussian approximation; not a fitted mean or physical bias bound.",
                "target": "Unchanged high255-minus-low255 validation zHEL contrast from shared_distance_response.",
                "inputs_allowed": "Only joint_F information matrices and saved contrast response; no joint_u or residual arrays.",
                "calculation": "Var=h.T (I+F)^-1 h, with F from discovery alone versus all43+1020 objects; original and isotropic priors retained."}
    (OUT / "protocol.json").write_text(json.dumps(protocol, indent=2) + "\n")
    inputs = {}

    def read(path, key):
        inputs[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
        with np.load(path, allow_pickle=False) as data:
            return data[key]

    records = []
    bandmap = np.array([[1., 0, 0], [-1, -1, -1], [0, 1, 0], [0, 0, 1]])
    iso = np.linalg.cholesky(2 * .02**2 * np.linalg.inv(bandmap.T @ bandmap))
    for n in [9, 12]:
        if n == 9:
            d_path, v_path = BASE / "shared43/projected-modes.npz", BASE / "shared1020/projected-modes.npz"
        else:
            d_path, v_path = BASE / "expanded12/discovery-projected-modes.npz", BASE / "expanded12/validation-projected-modes.npz"
        fd = read(d_path, "joint_F").sum(0)
        fv = read(v_path, "joint_F").sum(0)
        # These raw Gram matrices have already-scaled calibration columns and
        # unscaled observer columns. The saved h instead includes .02 observer.
        h_scaled = read(ROOT / f"runs/research_2026_09_26/shared_distance_response_{n}/responses.npz", "contrast")
        h_raw = h_scaled.copy()
        h_raw[n:] /= .02
        transforms = {"systematics_only": np.eye(n + 3)[:, :n],
                      "inherited_observer": np.diag(np.r_[np.ones(n), [.02] * 3])}
        isotropic = np.eye(n + 3)
        isotropic[n:, n:] = iso
        transforms["isotropic_observer"] = isotropic
        for name, transform in transforms.items():
            h = h_raw @ transform
            f_discovery = transform.T @ fd @ transform
            f_all = transform.T @ (fd + fv) @ transform
            variance_d = float(h @ np.linalg.solve(np.eye(len(h)) + f_discovery, h))
            variance_all = float(h @ np.linalg.solve(np.eye(len(h)) + f_all, h))
            assert 0 <= variance_all <= variance_d + 1e-12 <= h @ h + 1e-10
            eigen, vectors = np.linalg.eigh(f_all)
            zero = eigen < max(eigen) * 1e-12
            null_variance = float(np.sum((h @ vectors[:, zero])**2))
            records.append({"calibration_modes": n, "observer_family": name,
                            "prior_sd_mag": float(np.linalg.norm(h)),
                            "discovery43_conditional_sd_mag": float(np.sqrt(variance_d)),
                            "all1063_conditional_sd_mag": float(np.sqrt(variance_all)),
                            "all1063_fraction_of_prior_variance": variance_all / float(h @ h),
                            "latent_dimensions": len(h), "numerical_information_rank": int((~zero).sum()),
                            "target_variance_along_numerical_null": null_variance})
    result = {"scope": protocol["scope"], "records": records,
              "qualification": "Uses all design information but no observed residual values or new fitted amplitudes. Selection, unknown modes, intrinsic populations and arbitrary gray evolution remain outside this calculation."}
    (OUT / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    inputs[str(Path(__file__).relative_to(ROOT))] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    outputs = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.iterdir() if p.is_file()}
    (OUT / "manifest.json").write_text(json.dumps({"inputs_sha256": inputs, "outputs_sha256": outputs}, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
