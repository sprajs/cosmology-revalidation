# Sources and attribution

Scientific inputs retain their upstream ownership and terms. This repository supplies no new license for third-party data, model assets or software. Exact file identities and retrieval URLs are in `provenance/inputs.json`; package versions and dependencies are fixed in `uv.lock`.

| Resource | Frozen source |
|---|---|
| Compact Pantheon+ data | [Cobaya sn_data](https://github.com/CobayaSampler/sn_data/tree/61d96434cafc2770928322c38e5a750e686368ae), commit `61d96434cafc2770928322c38e5a750e686368ae` |
| DESI DR2 BAO data | [Cobaya bao_data](https://github.com/CobayaSampler/bao_data/tree/bb0c1c9009dc76d1391300e169e8df38fd1096db), commit `bb0c1c9009dc76d1391300e169e8df38fd1096db` |
| DES release models | [DES-SN5YR](https://github.com/des-science/DES-SN5YR/tree/e3493cb3b9505fc1f3f392d887364b50dae21439), release 1.3 commit |
| RAISIN released products | [RAISIN_DataRelease](https://github.com/djones1040/RAISIN_DataRelease/tree/a383c4bd03c9fbfd64bf5bfda525aec38d32b39c), commit `a383c4bd03c9fbfd64bf5bfda525aec38d32b39c` |
| Signed RAISIN precursor photometry | [RAISIN_cosmo](https://github.com/djones1040/RAISIN_cosmo/tree/b888214a5cbae38ac0bf886488ce7734f5b77e87), commit `b888214a5cbae38ac0bf886488ce7734f5b77e87` |
| Coherent 2021 simulation tables | [RAISIN_cosmo](https://github.com/djones1040/RAISIN_cosmo/tree/aaa709ead7a7339d56a2a4604d78991329d9368f), commit `aaa709ead7a7339d56a2a4604d78991329d9368f` |
| CSP physical throughput | [Pantheon+ calibration release](https://github.com/PantheonPlusSH0ES/DataRelease/tree/7fc680548d1ea9fe6ee983a89d8b5635bb5784a8), with retained upstream filter README |
| CSP DR3 photometry | [CSP data archive](https://csp.obs.carnegiescience.edu/data); exact TGZ bytes hashed |
| Published host-age tables | Chung et al. supplementary archive, DOI [10.1093/mnras/staf497](https://doi.org/10.1093/mnras/staf497); frozen ZIP and extracted tables |
| Current HST exposures | [MAST](https://mast.stsci.edu/); per-file product URLs and SHA-256 in the manifest |
| HST PAM and flat references | STScI WFC3 pixel-area map and [CRDS](https://hst-crds.stsci.edu/); exact URLs/filenames in the manifest |
| HST dark RAWs and bad-pixel/gain references | Fourteen original MAST RAWs and three CRDS tables; exact product URLs and hashes in the manifest |
| Native dark FLTs | Two local CALWF3 3.7.3 first-eight-read products; [native execution provenance](../provenance/hst-dark-native.json) |

The population kernels retain the explicitly implemented B13/MD14 cosmic star-formation and C14/W26 delay-law parameterizations from the earlier investigation. Their code constants and assumptions are documented in `methods/populations.md`; they are empirical functions, not laws derived from first principles.

DES objective exports, fit-summary arrays, projected calibration matrices and original frozen HST design arrays were generated in the earlier investigation. They are labelled as local derivatives. `provenance/code-origin.json` maps retained or adapted numerical code to its original file hash and symbols. Rewritten interfaces and concise methods form this repository; the preserved source snapshot remains the record for historical decisions and unresolved experimental branches.
