# Why the first signed-light-curve distance shift is not an identified correction

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../README.md) for current execution status.

The first predeclared object, DES16C1cim, now has a stable restricted native
profile over shape and empirical extinction at fixed header peak and fixed
covariance. Restoring nine measured negative epochs to 65 positives changes
the global best-fit distance parameter by **−0.16495 mag**. However, the
signed fit retains a high-shape solution only **0.11574 in Q** above its
minimum, with a distance response of **+0.01417 mag**. The large point shift
therefore switches between weakly distinguished template branches.

The [native report](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/astra_design/raisin_signed_refit/fixed-c-profile/global-profile/report.md)
preserves the frozen specification, numerical/domain gates, all branches,
alternative covariance anchor and limitations. The covariance-anchor
sensitivity gives −0.16564 mag for the global-minimum shift, but does not
resolve the competing template branches. This is a conditional diagnostic,
not a measured physical dust correction or a revised cosmological distance.

[Conditional native profiles](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_profile_decomposition/conditional-shape-profile.png)

The root's separately recorded [descriptive decomposition](../specifications/experiments/raisin_profile_decomposition/protocol.json)
was specified **after** seeing the first global shift. It reuses saved native
means, with no new fit or model evaluation. At the positive-only optimum's
shape and extinction, adding the negative values and refitting amplitude
changes the distance by **+0.012454 mag**. Allowing shape and extinction to
move then changes it by **−0.177403 mag**. This is an exact path decomposition
of the fitted response; it is not a unique causal split between physical
extinction and intrinsic luminosity.

For the common covariance C, partition rows into positives A and restored
negatives N. The Gaussian quadratic decomposes exactly:

    Q_B = r_Aᵀ C_AA⁻¹ r_A
          + (r_N − C_NA C_AA⁻¹ r_A)ᵀ C_N|A⁻¹
            (r_N − C_NA C_AA⁻¹ r_A),
    C_N|A = C_NN − C_NA C_AA⁻¹ C_AN.

| Parameter point | Positive marginal Q | Added-negative conditional Q | Full signed Q |
|---|---:|---:|---:|
| Positive-only optimum | 56.69454 | 37.56380 | 94.25833 |
| Same shape/extinction, amplitude refit on signed rows | 56.86929 | 37.21316 | 94.08245 |
| Signed global optimum | 58.85709 | 35.10937 | 93.96646 |

The global switch improves the full signed Q by only **0.29188** relative
to the positive-only parameter point. A cost of 2.16255 in the positive
rows is offset by a gain of 2.45443 in the negative rows. The block identity
closes to `1.42e-14`. These are comparisons of parameter points within a
specified metric; differences between separately minimized Q values on
different row sets are not a likelihood-ratio significance.

The shape profile's broad range of nearly equal Q demonstrates weak distance
identification within this conditional model. The plotted Q=1 guide and any
reported profile ranges are not calibrated confidence intervals. A valid
bias estimate still needs signed/censored selection-aware likelihood or
matched recovery simulations, justified measurement/population assumptions,
timing propagation and full-cohort uncertainty. The underlying Gaussian
metric applied to positive-only data is itself not the likelihood conditioned
on that sign selection.

[Numerical result](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_profile_decomposition/result.json),
[restored epoch ledger](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_profile_decomposition/negative-epochs.csv),
[decomposition source](../code/raisin_profile_decomposition.py),
[figure source](../code/raisin_profile_figure.py), and
[exportable PDF](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_profile_decomposition/conditional-shape-profile.pdf)
preserve the calculation and display inputs. The completed
[independent global-array review](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_profile_solver_review/global-profile-review/README.md)
reproduces all 86,813 mean vectors/four metrics, native stream replies and
numerical gates, and this four-point decomposition to `4.27e-14`.
Allowing exact amplitude freedom within the sampled Q≤Qmin+1 set widens the
signed distance range to 42.47139–42.86564 mag; it is still a descriptive set,
not a calibrated interval. The Q≤Qmin+9 sets hit the declared shape-domain
edges, so broad envelopes remain domain-truncated. No full research direction
has been closed by this single-object numerical result.
