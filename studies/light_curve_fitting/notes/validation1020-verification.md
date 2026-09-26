# Independent arithmetic check of the frozen 1,020-object validation

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../README.md) for current execution status.

The saved conditional validation scores pass an independent recomputation from
the frozen coefficients, per-object scores, and six-coordinate sufficient
statistics. This checks the saved arithmetic; it does not refit the light curves
or include cross-object calibration and SALT3 uncertainty. The reproducible
check is [`validation1020_verify.py`](../code/validation1020_verify.py),
with its separately saved [verification JSON](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/validation1020_independent_verification/verification.json).
The original analysis and its artifacts were left unchanged.

For the published mask, the recomputed frozen observer score gain is **27.2332138551**
from matched filter **192.4012835777** and information **330.3361394452**.
The frozen rest/phase gain is **−0.6923843631**. Retaining only the first epoch
of each exact duplicate group changes those gains to **28.2205246050** and **1.6871389643**.
Every per-object score agrees with its saved value within **8.9×10⁻¹⁶**.

I enumerated all **1,024** ten-field sign vectors using field sums rebuilt from
object scores. Exactly **one** assignment meets or exceeds the observed
observer, rest/phase, paired observer-minus-rest, or shared-posterior observer
statistic in either duplicate arm. This reproduces each reported conditional
field-sign tail of **1/1,024**. Recounting exact `(MJD, band, flux, quoted error)`
tuples in all 1,020 saved object exports finds **37 objects**, **59 groups**, and
**147 extra epochs**. The first-exact-duplicate arm retains **39,459** of the
published **39,606** epochs.

The shared three-coefficient result was recalculated as one Gaussian integral
over the *summed* validation sufficient statistics. With frozen discovery mean
`c`, covariance `V`, and validation totals `u=Σuᵢ`, `F=ΣFᵢ`, define
`A=V⁻¹+F` and `b=V⁻¹c+u`. The dense integral is
`(bᵀA⁻¹b − cᵀV⁻¹c − log det V − log det A)/2`.
It yields observer gains **59.2291792505** (published) and **59.7246592784**
(duplicate sensitivity), and frozen rest/phase gains **30.5925837705** and
**31.4919109849**. The independently recomputed 1,024 signed integral values
agree with the saved array within **1.8×10⁻¹²**. This calculation integrates
the single discovery posterior once; it does not multiply 1,020 separately
marginalized object likelihoods.

The observer coefficient file matches the hash frozen in the validation input
manifest, the result's input hash, and the reported coefficient mean/covariance.
Its mean also equals the saved 43-object discovery ridge posterior solution.
The rest/phase coefficient file matches its pre-score amendment hash. The
executed analysis source matches the result's source hash; inspection shows
the primary score uses those loaded coefficient arrays, while its reported
transfer amplitude is calculated only after the fixed score. No coefficient
estimated from these 1,020 validation objects replaced the frozen prediction.

All Gaussian Z values and field-sign tails in this validation remain **conditional
on fixed nominal calibration, SALT3, per-object covariance, and selection**.
They are not full-systematics significance measures or evidence of an observed
calibration bias or cosmological correction.
