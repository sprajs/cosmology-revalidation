# Infrared source separation

This experiment combines the actual DES galaxy positions and public empirical SPIRE beams. It measures conditional linear information about separating a host from neighbouring sources. It does not estimate a host dust luminosity or append a supernova correction.

First recover the DES deep catalogues and SPIRE images using the [host photometry instructions](../../notes/host-transport-results.md). That branch creates the isolated `.work/host-transport/.venv` environment used here. From the repository root:

```bash
.work/host-transport/.venv/bin/python studies/host_ages/code/infrared_resolution/acquire.py
OPENBLAS_NUM_THREADS=1 .work/host-transport/.venv/bin/python studies/host_ages/code/infrared_resolution/analyze.py
OPENBLAS_NUM_THREADS=1 .work/host-transport/.venv/bin/python studies/host_ages/code/infrared_resolution/validate.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python studies/host_ages/code/infrared_resolution/figure.py
```

`acquire.py` downloads and hashes three official ESA Neptune beam images and their source README. It checks any existing acquisition identities instead of accepting silently changed inputs. Large downloads and generated object-level tables remain under ignored `.work/infrared-resolution/`.

The figure uses the main repository environment, which supplies Matplotlib, and records the hashes of its data, source and exported PNG.

The analysis selects all 265 eligible eight-band hosts with valid measurements in all three SPIRE bands, independently of their flux sign. It reconstructs optical neighbours, actual image sampling and masks, fits empirical beam centres, and computes nearest-pair and multi-source information under independent equal-variance pixel errors. No infrared image flux enters the fitted response.

The validator independently reconstructs neighbour geometry and linear projections, checks Gaussian integrals, tests conditional noise recovery, and measures sensitivity to source radius, fit radius, beam support, pixel averaging and numerical rank. The empirical beam profile remains a shared input. These checks validate the declared calculation, not the assumed image covariance or completeness of the optical source list.

See the [scientific result and equations](../../notes/infrared-resolution-results.md), [exploratory design](design.json), [validation design and amendment](validation-design.json), [compact information result](../../results/infrared_resolution/information.json) and [validation record](../../results/infrared_resolution/validation.json).
