# Handoff: physical-equation audit

Saved 21 September 2026 at the user's request to finalize within five minutes and preserve progress for later agents.

## User's objective and delivered work

The user asked to gather and write the actual underlying physical equations for luminosity, extinction and the rest of this project, establish them from physical foundations, correct errors and evaluate whether assumptions are justified. The deliverable is [README.md](README.md) and four detailed LaTeX-in-Markdown chapters, an equation index, finding register, four executable validation scripts and their new results. This is an equation/implementation audit, not a replacement empirical cosmology analysis.

Start with the README and `findings.json`; then read the relevant chapter. Stable equation IDs L, D, O, C and P allow later tests to cite exact assumptions. `runs/physics_audit/manifest.json` records hashes and commands. The index records chapter/section/line for each numbered equation.

## Completed changes

- Added luminosity equations from hydrodynamics/nuclear energy and Ni–Co decay through deposition, expansion work, diffusion, Arnett limits, spectral transfer and the empirical SALT interface.
- Added extinction equations from Maxwell/material response and grain cross sections through transfer, dust screens, CCM/O'Donnell/F99, broad-band response, geometry and dust-map interpretation.
- Added the emitted-spectrum → redshift → instrument photons → calibrated flux/magnitude chain, K corrections, detector/noise, standardization and selection boundaries.
- Added GR/FLRW/acceleration, CPL/q bins, cosmic time, stellar-population/SFH/DTD weighting, BAO/CMB and directional-model boundaries.
- Fixed `scripts/cosmology/core.py`: explicit piecewise q-bin E(z); reject outside bin domain in E, q and integral.
- Fixed `scripts/phase2/hierarchy/core.py`: q-bin distances return NaN outside valid domain under JIT; `scripts/phase2/hierarchy/run.py` validates data before sampling.
- Fixed `scripts/phase2/independent_flux/engine.py`: evaluated filters must have full nonzero response inside spectral support; unused unsupported filters remain preparable.
- Fixed `scripts/standardization/dust_identifiability.py:slab`: transparent limit, stable thin limit, invalid optical-depth rejection.
- Clarified acceleration density/pressure units in `docs/first-principles.md`; added navigation from that file and root README.

## Verification

All four new check scripts passed on the saved implementations. Logs and JSON are under `runs/physics_audit/`. The manifest contains exact rerun commands; normal reruns take seconds with existing local environments/data. Cosmology checks use root `.venv` and spawn the separate hierarchy JAX interpreter. Extinction checks freshly compile bounded vendor routines in temporary storage.

Important evidence: 12 saved reference objects / 528 flux epochs unchanged exactly; positive slab experiment grid unchanged exactly; inspected Pantheon/hierarchy redshifts inside their respective 2.5/1.3 bounds. No existing cosmological posterior or historical output was overwritten. Before-edit snapshots for the root's edited kernels and a pre-edit flux fixture are in `runs/physics_audit/before/` and `flux-before.json`; the originally clean tracked slab function is also available in Git HEAD.

## Scientific questions still open, not silently resolved

1. Extrapolated historical and spline F99 curves permit negative passive-screen extinction in the saved simulated low-RV population. Independent raw-header census: 71,946 rows; 4,601 negative historical A8000 and 4,791 negative spline A8000. These are simulation truth rows, not observed objects or selected cosmology weights. A supported nonnegative dust family/population must be chosen, refitted and propagated through selection/BBC before a distance or cosmology consequence can be claimed. No arbitrary clamp, R_V truncation or vendor rewrite was made.
2. SALT is an empirical surface, not a closed explosion luminosity calculation. The new nuclear/diffusion equations are established with regime limits and checked synthetically; they were not installed as a purported first-principles prediction of the data.
3. Some native SALT nodes permit negative spectra for part of the tested x1 interval: 703 of 27,511 nodes, not 703 observations. Actual passband/epoch/parameter impact remains to be traced before changing learned surfaces.
4. Host-age/progenitor-delay mapping, mean versus median, population transport, dust/intrinsic-colour separation, selection and calibration require empirical tests. Correct formulas do not establish their assumed distributions.
5. Imported CMB chains and external full SNANA/Cobaya/Planck/ACT code are not comprehensively audited. See explicit coverage limits rather than treating the chapter as a certificate for every downloaded repository.

The finding register provides concrete next tests. These unresolved questions are scientific model choices, not unfinished local fixes that can safely be filled in by guessing. No new cosmology sampling or full production simulation rerun is needed just to finish the equation record.

## Workspace precautions for a future agent

The tree already contained extensive untracked phase-2 and dust-audit work when this task began. Other phase-2 files appeared while work proceeded, suggesting concurrent workspace activity. Do not revert, clean or stage the entire tree. Limit changes to the named audit files and clearly identified kernels; check current hashes before resuming. No commit, push, upload or external communication was performed.

The original primary data, vendor sources and older results remain preserved. Git alone is not the full local data bundle. Do not reinterpret an old manifest as describing a newly modified script; new science runs need fresh provenance. No persistent memory files were changed.
