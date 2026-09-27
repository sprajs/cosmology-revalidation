# Nebular-line overlap test

This tests whether measured Hα/Hβ information changes incremental age–brightness prediction on the existing spectroscopic ZTF/TITAN cohort. Gas attenuation, stellar Av and supernova extinction remain distinct quantities.

Run from the repository root after restoring the [galaxy observation inputs](../galaxy_validation/README.md):

~~~bash
.venv/bin/python studies/host_ages/code/nebular_dust/acquire.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python studies/host_ages/code/nebular_dust/analyze.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python studies/host_ages/code/nebular_dust/validate.py
~~~

The existing repository environment supplies NumPy, pandas, SciPy, Astropy and iminuit. No new environment or numerical catalogue download is required. The acquisition script verifies inherited numerical input hashes and restores pinned public measurement-method pages if missing; changed upstream bytes fail the identity check.

The frozen design records exploratory timing, common cohorts, models and limitations. The primary ratio comparison uses detected Balmer lines and first-order propagated errors. A broad signed contrast retains weak/negative measurements, with its nonlinear uncertainty limitation measured separately. Star-forming classification, stronger line cuts, historical error scaling, line covariance, marginal age errors and a conditional Hβ continuum correction are explicit sensitivities.

The analysis retains the existing SALT covariance likelihood and physical-host folds. Each full fit is checked using an independently expanded objective and SLSQP. Compact records go to [the results directory](../../results/nebular_dust); generated cohort and prediction tables go to the ignored .work/nebular-dust directory.

The validator checks source identities, exact integer identifiers, fitted/predicted records, analytic gradients, signed-flux units and nonlinear propagation. The [scientific note](../../notes/nebular-dust-results.md) explains the findings and inference limits. This branch cannot produce a certified unblinded or causal cosmological correction.
