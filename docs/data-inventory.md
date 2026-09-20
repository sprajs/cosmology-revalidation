# Data inventory

Checked 20 September 2026. **Local** means acquired bytes, not a completed scientific replication. **Public but not recovered** describes a retrieval failure; it does not mean private. **Not located** describes the limits of this search. Exact versions, hashes and failed attempts are retained in `catalog/`.

## Supernova and host data

| Product | Local location | Contents and limitations |
|---|---|---|
| Chung Paper I supplementary tables | `data/host_ages/chung2025/` | Original ZIP plus tables: 199 G11 rows and 102 R19 rows. Ages, uncertainties, SN IDs and HRs. These are table rows, not 301 proven independent SNe. Full age posterior chains are not supplied here. |
| Gupta 2011 | `data/host_ages/gupta2011/` | CDS ReadMe and table2.dat; historical host photometry/properties, original ages and SN measurements. [CDS J/ApJ/740/92](https://cdsarc.cds.unistra.fr/ftp/J/ApJ/740/92/) |
| Rose 2019 | `data/host_ages/rose2019/` | CDS ReadMe and tables 1, 2, 3, 7; original global/local host photometry/properties and comparison material. [CDS J/ApJ/874/32](https://cdsarc.cds.unistra.fr/ftp/J/ApJ/874/32/) |
| Rose age-fitting code and auxiliary data | `sources/repos/benjaminrose__mc-age/` | Pinned source snapshot. Its historical implementation is not Chung's exact updated FSPS environment. The SFD dust-map submodule was acquired separately at its exact recorded commit, under `sources/repos/kbarbary__sfddata@7a5fe7fadf086561ba4748756e59a4c51d0ec632/`; it has not been wired into an execution environment. |
| Pantheon+ modern release | `sources/repos/PantheonPlusSH0ES__DataRelease/` | Light curves, calibration, light-curve fits, distance table, covariance and likelihood material. The main distance table has 1,701 rows and 1,543 distinct CIDs. Includes SH0ES products, which must not silently become part of an uncalibrated SN-only analysis. |
| Pantheon+ snapshot cited by W26 | `sources/repos/PantheonPlusSH0ES__DataRelease@7fc6805/` | Full exact commit named by W26, preserved alongside the current release. W26's matched-row selection and paper-specific regression script were not located. |
| Original DES-SN5YR, tag 1.3 | `sources/repos/des-science__DES-SN5YR@1.3/` | Full repository archive and extracted files. HD contains 1,829 rows: 1,635 DES and 194 low-z. Original distances, covariance, simulations, SALT3, classification and chains. |
| DES-Dovekie | `sources/repos/des-science__DES-SN5YR/` | Commit c9a4fca, also the DES commit cited by W26. HD contains 1,820 rows: 1,623 DES and 197 low-z. Distances, packed precision matrices, metadata, calibration models, classification and published chains are local. 250 simulated light-curve files remain LFS pointers after retrieval failures; see below. |
| Full DES scene-modelled photometry | The two DES repositories' `0_DATA/` trees | HEAD/PHOT FITS pairs are present. Header inspection finds **19,706 DES candidates**, 185 Foundation and 342 other low-z entries at this earlier stage. These counts are not the final cosmology sample. |
| DES difference-image photometry | [Zenodo 14332953](https://zenodo.org/records/14332953); `data/des-diffimg/` | README and both HEAD/PHOT FITS files acquired; all three match Zenodo's MD5 checksums. Public 31,636-candidate, 18,106,216-photometry-row release. The release explicitly says these photometry products are **not suitable for cosmological inference** and were not used for DES parameter estimation. |
| Dovekie calibration inputs/software | `sources/repos/bap37__Dovekie/` | Pinned calibration repository, including stellar photometry/filter and calibration products present in the public tree. It is not a reconstruction of every upstream survey exposure. |
| Galaxy/SN population model | `sources/repos/wisemanp__des_sn_hosts/` | Model scripts, configurations and model-data products. Its 250 MB LFS host library was recovered and SHA-256 checked. Two machine-specific symlinks were preserved in the archive but not materialized. Not a guaranteed exact W26 environment. |
| Sah et al. analysis products | `sources/repos/Shin107__Anisotropy-in-Pantheon-Plus/` | Public Python/notebook scripts, Pantheon data/covariance and mean/median age correction tables. Some upstream `.pkl` files are zero bytes; not executed or treated as completed chains. |
| Supporting SN likelihood data | `sources/repos/CobayaSampler__sn_data/` | Compact data/covariance products for Pantheon+, DESY5, DES-Dovekie, Union3 and historical datasets. Retained as alternative upstream packaging, not an instruction to mix releases. |
| GSWLC-2 galaxy catalogues | `data/dust/GSWLC-X2.dat.gz`, `GSWLC-M2.dat.gz` and column-description PDF | Both the combined best-depth catalogue (659,229 rows) and the medium-depth parent catalogue used by Salim are acquired. Columns include stellar mass, SFR, FUV/B/V attenuation and flags. Missing values are −99. These are not automatically the exact quality-selected 230,000-galaxy attenuation-curve sample. No explicit per-galaxy `R_V` column is listed in the schema. [Author release](https://salims.pages.iu.edu/gswlc/) |
| YONSEI residuals used by Son | **Not local** | Son says these were obtained privately. Kim 2019 says the full catalogue is available from its corresponding author on request. The published paper and some tables are available; the precise privately transmitted residual file is not identified. |
| TITAN 6,983 host histories/posteriors | **Not local** | Murakami's paper uses TITAN DR1 host data and describes forthcoming catalogue/selection releases. The [official DR1 page](https://titan-snia.github.io/dr1.html) currently has no download. Public ATLAS measurements alone are not the fitted host histories used in this paper. |

Chung table 1 includes HR values already adjusted for redshift, with original G11 HRs in parentheses. Preserve both; a parser must not silently discard the original values. Table 2 uses a different historical HR source. The July correction concerns Figure 1, not a documented replacement of these downloaded supplementary tables.

## BAO, CMB and existing chains

| Product | Local location | Interpretation |
|---|---|---|
| DESI DR2 BAO vector and covariance | `sources/repos/CobayaSampler__bao_data/desi_bao_dr2/` | Compact published cosmology inputs, including the 13-observable combined vector/covariance. Full survey spectra and galaxy catalogues are not required at this reproduction level and have not been mirrored. |
| Official DESI DR2 reference chains | `data/bao/desi-dr2-reference/` | Three flat-CPL combinations: BAO+CMB, BAO+CMB+Pantheon+, BAO+CMB+DESY5. Input/updated YAML, four chain files per combination, marginal summaries and ancillary output acquired. **All 30 files match the official SHA-256 list.** These are unmodified DESI reference results, not Son's age-corrected chains. |
| Full DESI archive catalogue | `data/bao/desi-dr2-chains-README.md`, `data/bao/desi-dr2-chains.sha256sum` | Describes the broader public model/likelihood grid; only the three directly relevant combinations are mirrored here. [Official directory](https://data.desi.lbl.gov/public/papers/y3/bao-cosmo-params/) |
| ACT DR6 + Planck lensing likelihood | `sources/repos/ACTCollaboration__act_dr6_lenslike/`, `data/cmb/ACT_dr6_likelihood_v1.2.tgz` | Code plus NASA LAMBDA v1.2 likelihood data. Supports the ACT+Planck lensing combination cited by the DESI baseline. This is likelihood-level data, not sky maps. |
| Planck PR4 high-ell CamSpec | `data/cmb/CamSpec_NPIPE.zip` | Public native likelihood data linked by Cobaya, release v1. |
| Planck PR3 low-ell | `data/cmb/Planck_PR3_baseline_R3.00.tar.gz` | ESA product 151902, the baseline likelihood bundle referenced by Cobaya's Commander/SimAll definitions. Neither the likelihood library nor its compilation toolchain has been installed. |
| Likelihood and simulation source | `sources/repos/CobayaSampler__cobaya/`, `sources/repos/RickKessler__SNANA/` | Current pinned source and configuration references. They are not a claim about the historical software versions used by every paper. CAMB/CosmoSIS/FSPS/SNANA execution environments remain future work. |

Son Section 4.1 specifies **Planck PR3 Commander and SimAll, PR4 NPIPE CamSpec, and Planck PR4 + ACT DR6 lensing**. The saved DESI YAML confirms this combination. Do not replace it with an unspecified “Planck 2018” likelihood, a compressed prior, or the newer Planck/ACT/SPT combination in DES-Dovekie without explicitly defining a separate experiment. [Son](https://arxiv.org/abs/2510.13121), [DESI archive documentation](https://data.desi.lbl.gov/public/papers/y3/bao-cosmo-params/README.md)

## Version and format traps found during acquisition

1. **DES-Dovekie `STATONLY.npz` and `STAT+SYS.npz` store packed upper-triangular inverse covariance**, not an ordinary full covariance. Follow the unpacking convention in the release likelihood. No matrix inversion or scientific likelihood evaluation has been done here.
2. DES-Dovekie's file named `DES-Dovekie_HD.csv` is actually a SNANA whitespace table with comment lines, `VARNAMES:` and `SN:` prefixes. The old DES HD is ordinary CSV. Several README/default-path names retain older labels; never infer a file format from its suffix alone.
3. The Dovekie metadata rows are not in the same order as the HD/matrix. Use explicit ID joins and preserve HD order. The classification README also retains a 1,635-object count; the actual Dovekie distance file has 1,623 DES rows.
4. Original DES covariance conventions differ from Dovekie's. Preserve their respective likelihoods and statistical-error treatment. Pantheon+ uses its own ordering and repeated measurements of some SNe; rows and unique explosions are different quantities.
5. The original DES Zenodo DOI points to a **1.2 ZIP**, whereas the repository recommends **tag 1.3** for the old analysis. Its metadata and checksum are saved; the redundant 1.2 ZIP is not mirrored. Tag 1.3 is not asserted to be Son's exact unspecified snapshot.
6. GitHub source archives do not necessarily include Git LFS objects or submodule contents. The hydration log distinguishes recovered data from pointers. DES's 250 current mock-file pointers describe about 400 MB; direct object URLs returned 404, and the LFS batch endpoint returned 403 in this environment. This does not establish that the data are private or lost.
7. None of these datasets should be multiplied together as independent likelihoods merely because they have different names. Shared SNe, calibrations, low-z anchors and lensing information require an overlap/covariance audit.

## Exact snapshot ledger

Full 40-character hashes, original archive paths, materialized-file counts, symlinks and LFS state are in [repository_inventory.json](../catalog/repository_inventory.json). W26's explicit paper references are Pantheon+ **7fc6805** and DES **c9a4fca**. All other current snapshots are collection-time pins, not undocumented claims of author-version equivalence.
