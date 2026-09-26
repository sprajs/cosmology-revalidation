# Pinned CALWF3 RAW-to-FLT closure feasibility

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../README.md) for current execution status.

This is a header and runtime inventory for the two already acquired science RAWs, `icxoi1bcq` and `icxoi4hgq`. No CALWF3 build or execution, reference download, or science-pixel read was performed. The [machine-readable inventory](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_feasibility/inventory.json) freezes RAW/FLT hashes, their header switches and reference names; the [HEAD protocol](../specifications/experiments/hst/calwf3_native_feasibility/head-protocol.json) and [responses](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_feasibility/head-results.json) record eight official CRDS metadata requests in the declared 60-second limit.

## Execution closure and blockers

Both current archive FLTs report `CAL_VER=3.7.3 (Jan-07-2026)`, matching the pinned HSTCAL source [version header](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_variance_source/source-3.7.3/wf3version.h). Their CRDS contexts differ (`hst_1339.pmap` search, `hst_1337.pmap` template), but *all eleven explicit reference filenames per RAW exactly match its corresponding archive FLT*. The RAW switches for DQI, zero signal/offset, dark, bias level, nonlinearity, flat, cosmic ray, unit and photometric calibration say `PERFORM`; corresponding FLT switches say `COMPLETE`. `RPTCORR=OMIT`, and `DRIZCORR=PERFORM` remains pending in the FLT: a RAW-to-FLT replay is separate from AstroDrizzle.

There is no local `calwf3.e` or `wf3ir.e`, no CMake, and no full pinned source checkout. The preserved source directory contains twelve selected files, not the buildable tree. The pinned [root CMake file](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_feasibility/source/CMakeLists.txt) requires CMake 3.11+, C99, pkg-config and CFITSIO and adds `ctegen2`, `cvos`, `hstio`, `lib`, `tables`, and `pkg` subdirectories. `/usr/bin/gcc`, `make`, `pkg-config`, CFITSIO 4.7.0 and its shared library/header are present. `gfortran` is absent; the [upstream install guide](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_feasibility/source/INSTALL.md) lists it in Linux package recipes, while the inspected top-level CMake activates C only. Whether any nested target needs Fortran remains untested. No `crds`, `wfc3tools`, `drizzlepac`, or `stenv` Python package is present in the project environment. The native C executable can use explicit `iref$` filenames without those Python wrappers once it exists. The [upstream HSTCAL README](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_feasibility/source/README.md) explicitly requires `iref` to identify the WFC3 reference directory.

Six of fourteen distinct named references are already local and hashed in the inventory (31,680,000 bytes total): `t2c16200i_ccd.fits` and `q911321mi_osc.fits` in `quality_reference_audit/reference_files/`; `4ac1921li_pfl.fits`, `4ac1822ni_dfl.fits`, and `4ac1831ei_dfl.fits` in `calwf3_variance_source/`; and `w3m18525i_idc.fits` in `pixel_acquisition/pam_validation/`. Eight are missing:

| Reference | Use | Official HEAD bytes |
|---|---|---:|
| `3562029fi_bpx.fits` | search bad-pixel table | 912,960 |
| `3562018mi_bpx.fits` | template bad-pixel table | 964,800 |
| `u6a1748ri_crr.fits` | cosmic-ray table, common | 14,400 |
| `3562025ji_drk.fits` | search dark | 302,607,360 |
| `3562012ji_drk.fits` | template dark | 302,607,360 |
| `a2412448i_lin.fits` | nonlinearity file, common | 79,902,720 |
| `8ch15233i_imp.fits` | photometric calibration table, common | 83,520 |
| `3562021pi_mdz.fits` | drizzle table named in header | 40,320 |

All eight official CRDS `HEAD https://hst-crds.stsci.edu/unchecked_get/references/hst/<name>` requests returned HTTP 200 with these `Content-Length`s. The total missing transfer/storage is **687,133,440 bytes (655.3 MiB)** if every header-named file is staged. `MDRIZTAB` is likely downstream of FLT production, but preserving the exact header and staging all named references avoids silently changing calibration behavior; the narrow native dependency could only be established by source trace or an authorized run. The two RAWs (50,843,520 bytes) and two comparison FLTs (33,223,680 bytes) are already local. Isolated replay would additionally create two IMAs and two FLTs plus build artifacts; their sizes are not established by this read-only inventory. Reserve at least ~1.5–2 GB free scratch **as an engineering allowance, not a measured requirement**, above the 655 MiB missing-reference transfer. The full pinned source archive/build size remains unknown.

## Minimal future isolated replay route

Acquire the full HSTCAL tree at exact commit `6a1147d7bdccb7e2a7b73af276f4597c6420fbb8` and a project-local CMake >=3.11, then configure/install into new private directories, leaving system packages and existing environments untouched. The pinned upstream [installation instructions](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_feasibility/source/INSTALL.md) support a private prefix and `-DENABLE_OPENMP=OFF`. A concrete proposed build is:

```sh
cmake -S "$PINNED_HSTCAL_SOURCE" -B "$PRIVATE_BUILD" \
  -DCMAKE_INSTALL_PREFIX="$PRIVATE_PREFIX" -DENABLE_OPENMP=OFF \
  -DCMAKE_C_COMPILER=/usr/bin/gcc
cmake --build "$PRIVATE_BUILD" -j2
cmake --install "$PRIVATE_BUILD"
```

Here `cmake` denotes a separately supplied project-local binary; none is installed now. First verify the resulting binary's version and dependencies. Then make an isolated `iref/` containing the fourteen exact basenames (symlinks to immutable validated references are sufficient if the native reader follows them), copy each RAW into a separate scratch working directory, and invoke one process at a time:

```sh
cd "$SEARCH_WORK"
env iref="$PRIVATE_IREF/" OMP_NUM_THREADS=1 \
  "$PRIVATE_PREFIX/bin/calwf3.e" icxoi1bcq_raw.fits
cd "$TEMPLATE_WORK"
env iref="$PRIVATE_IREF/" OMP_NUM_THREADS=1 \
  "$PRIVATE_PREFIX/bin/calwf3.e" icxoi4hgq_raw.fits
```

This is a **proposal, not an executed or validated command**. The [WFC3 handbook pipeline description](https://hst-docs.stsci.edu/wfc3dhb/chapter-3-wfc3-data-calibration/3-1-the-calwf3-data-processing-pipeline) documents `calwf3.e [options] input`; its [manual recalibration guide](https://hst-docs.stsci.edu/wfc3dhb/chapter-3-wfc3-data-calibration/3-5-manual-recalibration-of-wfc3-data) requires `iref` for the named reference files. An authorized closure would verify binary/source version, each staged reference's identity, processing switches, CALWF3 trailer, raw/IMA/FLT provenance and per-HDU SCI/ERR/DQ agreement against the current archive FLT. Version-string and filename agreement alone do not establish numerical replay or historical 2016–2017 processing equivalence.

## Quadrant/gain coordinate closure for later RAW diagnostics

The pinned source [doir.c:115](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_variance_source/source-3.7.3/doir.c#L115) calls `GetGrp`; the exact-commit [getgrp.c:46–48](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_feasibility/source/getgrp.c#L46) reads RAW SCI `LTV1=LTV2=5` into `offsetx/y`. [noiscalc.c:106–122](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_variance_source/source-3.7.3/noiscalc.c#L106) reads `trimx/y[0]=5` and sets the boundary to `AMPX + offsetx − trimx[0] = 512 + 5 − 5 = 512`, likewise y. Thus zero-index RAW active pixels 5–1018 divide at i/j=512; after trimming five pixels, the 1014×1014 FLT boundary is at zero-index x/y=507. RAW SCI headers are 1024×1024 with LTV=5; FLT SCI headers are 1014×1014 with LTV=0. [noiscalc.c:125–170](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_variance_source/source-3.7.3/noiscalc.c#L125) maps lower-left B, lower-right C, upper-left A and upper-right D. The matched CCDTAB `CCDAMP=ABCD`, `CCDGAIN=2.5`, `AMPX=AMPY=512` row gives A/B/C/D gain 2.34/2.37/2.31/2.38 electrons per DN. The source/header boundary is now explicit, but RAW DN-to-electron propagation remains a distinct conditional analysis until CALWF3 replay closes.
