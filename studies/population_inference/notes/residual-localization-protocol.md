# Fixed-vector wavelength and phase localization

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../README.md) for current execution status.

This is a post-discovery mechanism diagnostic, registered before reading its
partitioned scores. It neither provides an independent detection significance
nor selects a correction. Keep the exact43 discovery observer vector unchanged
and use all 1,020 validation objects, the nominal exported means/covariance and
the archived four-coordinate nuisance Jacobians, with the exact amplitude
column. Existing all-epoch results are already known.

The primary wavelength partition asks whether the pattern persists using
passbands with at most 1% of their positive photon throughput below rest-frame
3500 Angstrom. Compute the fraction from lambda*T(lambda), on a maximum
5-Angstrom observer grid, and the measured zHEL. This is a throughput criterion,
not an SED-weighted or photon-counted fraction. Declare thresholds 3000 and 4000
Angstrom as sensitivities. Evaluate each complement separately, and preserve
objects with zero or insufficient surviving information in the ledger.

Additional descriptive partitions are pre-peak and post-peak rest phase,
using the original fitted peak time and split phase<=0 versus >0. They are
not independent tests or prespecified population priors. No cuts depend on
observed flux, residual size, information gain or preferred sign.

First run design-only: record exact selected epoch indices, ranks and frozen
vector information before loading observed flux. Use marginal covariance
submatrices and recompute the weighted nuisance projection for each partition.
Do not subset a precision matrix or an already projected residual. Then open
the observed flux and compute M=prediction dot residual, I=prediction squared,
and G=M-I/2. Report M/I only where I>0 and label it descriptive. Retain all
rank-deficient objects, using their numerical rank with threshold 1e-10 of the
largest singular value; report their contribution separately. Reconstruct the
all-epoch aggregate and compare to the saved primary result as a gate.

Partitions are not additive after separate nuisance projection. Covariance
and tangent remain fixed at the original all-epoch fit; this is a linear
compensation diagnostic, not full nonlinear refitting after data exclusion.
Correlated calibration/training modes, population selection and cosmology
are unchanged. A persistence result rules out confinement to these excluded
epochs; it cannot rule out general spectral evolution or dust.

Independent motivation includes the tentative UV evolution studied by
[Wang et al., 2512.25064v2](https://arxiv.org/abs/2512.25064v2).
That study's model and calibration caveats are not used as an empirical prior
or as evidence that this DES pattern has the same cause.

Design-only amendment, before partitioned scores: the UV-exposed 3000/3500
complements have information below 3e-30 because their band-offset direction
is absorbed into each object's amplitude. Thus report no descriptive amplitude
when aggregate I<=1e-10, rather than dividing numerical roundoff by roundoff.
The 4000 complement has I=.00214 and remains explicitly labelled as almost
uninformative. The original native fit already imposed a representative
rest-wavelength range 3500–8000 Angstrom; the present throughput-tail cut tests
the additional blue tails of those accepted passbands.
