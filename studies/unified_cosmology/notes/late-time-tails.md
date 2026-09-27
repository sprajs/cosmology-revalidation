# Late-time CPL tails and prior sensitivity

The Dovekie supernova distances and DESI DR2 BAO permit a secondary, low-matter
ridge in the stated late-time CPL model. This matters for the uncertainty on
dark-energy evolution even though both the main solution and the low-matter
ridge imply present acceleration. The calculation uses the released corrected
supernova distances and full covariance. It does not establish that those
corrections capture every population effect.

The baseline is flat matter plus CPL dark energy over the observed redshift
range, with a free scale $H_0r_d$. It omits radiation and does not impose an
early-time cut on $w_0+w_a$. Its uniform priors are
$0.01<\Omega_m<0.99$, $5000<H_0r_d<15000$,
$-3<w_0<1$, and $-5<w_a<3$. There is no CMB information or empirical host-age
correction in this target. These bounds remain assumptions.

At fixed $\Omega_m$, every BAO prediction is proportional to
$y=(H_0r_d)^{-1}$. Thus its contribution is
$\chi^2_{\rm BAO}=Ay^2-2By+C$, with the constrained minimum at
$y=\operatorname{clip}(B/A,1/15000,1/5000)$. The common supernova intercept
removes the supernova dependence on this scale. We profile the remaining
$(w_0,w_a)$ coordinates from three starts at each of sixteen matter densities.

The [profile](../results/inference/late-tail-profile.json) gives:

| Fixed $\Omega_m$ | Best $w_0$ | Best $w_a$ | $\Delta\chi^2$ from best profile grid point |
|---:|---:|---:|---:|
| 0.010001 | −0.6908 | 0.8871 | 2.915 |
| 0.10 | −0.7776 | 0.9271 | 4.898 |
| 0.20 | −0.8840 | 0.8469 | 6.870 |
| 0.31416 | −0.8434 | −0.5357 | 0 |

The best grid point has $\chi^2=1638.11985$. The lower-matter ridge reaches
the declared lower prior bound, so its posterior mass cannot be described
without that bound. A profile height is not posterior mass or a Bayes factor.
The three optimization starts agree within $9\times10^{-11}$ in
$\chi^2$; independent direct supernova distance quadrature at four ridge
positions agrees within $2\times10^{-12}$. Those checks establish numerical
closure at the tested positions, not global uniqueness of the optimum.

The first two independent nested runs used 1000 live points and a random-walk
kernel. Both met their individual call-limit and weighted effective-size
checks, but [failed the declared cross-seed check](../results/inference/late-nested-comparison.json):
the largest marginal CDF difference was 0.0644, exceeding 0.05. Their inferred
$P(w_0+w_a>0)$ differed, 3.07% versus 1.37%, and their lower 2.5% matter-density
quantiles were 0.157 versus 0.251. Their similar evidence estimates did not
resolve the tail. These initial records are retained.

A [declared refinement](../results/inference/late-nested-refinement-design.json)
uses 5000 live points and random-direction slice sampling in two new seeds,
with the same physical target and numerical comparison gates. Both completed,
using 2.20 and 2.21 million likelihood evaluations. Their weighted effective
sizes were 27,460 and 27,750. The [refined comparison](../results/inference/late-nested-refined-comparison.json)
passes the declared gates: maximum marginal CDF difference 0.0218, and
unnormalized log-evidence difference 0.076 against a combined reported
numerical error of 0.079. Evidence here only checks two computations of the
same target; no model comparison follows.

The median matter densities are 0.3131 and 0.3126, the median $w_0$ values
−0.8440 and −0.8443, and the median $q_0$ values −0.3753 and −0.3776.
The corresponding central 95% $q_0$ intervals are [−0.565, −0.203] and
[−0.566, −0.206]. These support present acceleration within the stated
distance model and priors.

The coarse comparison gate does **not** establish precise tail probabilities.
The two runs assign 2.87% and 3.90% to $w_0+w_a>0$, and their lower 2.5%
matter-density quantiles remain 0.157 and 0.108. This between-run variation
exceeds the saved prior-volume jitter errors, which omit sampling and
missed-mode uncertainty. The separate runs are therefore retained; tiny
weighted nonacceleration probabilities are not converted into sigma claims.

The [conditional-cut design](../results/inference/late-tail-design.json)
specifies separate, post-hoc summaries conditional on $\Omega_m>0.1$
and on $w_0+w_a\leq0$. Their [computed effects](../results/inference/late-tail-conditional.json)
are visible in the following ranges across the two seeds:

| Conditioning assumption | Retained baseline probability | Lower 2.5% $\Omega_m$ quantile | Upper 97.5% $w_a$ quantile |
|---|---:|---:|---:|
| Original priors | 100% | 0.108–0.157 | 0.867–0.899 |
| $\Omega_m>0.1$ | 97.73–98.35% | 0.231–0.251 | 0.637–0.750 |
| $w_0+w_a\leq0$ | 96.10–97.13% | 0.264–0.266 | 0.485–0.515 |

Applying both cuts has nearly the same result as the second cut. These are
sensitivity calculations, not baseline restrictions. In particular, the
$w_0+w_a$ inequality does not reproduce a measured CMB likelihood. The tail
responds appreciably to these assumptions while the present-acceleration
conclusion remains unchanged for this conditional model.

Reproduce the profile with the pinned nested-sampler environment:

```bash
uv venv --python .venv/bin/python .work/unified-cosmology/inference/.nested-venv
uv pip install --python .work/unified-cosmology/inference/.nested-venv/bin/python \
  -r studies/unified_cosmology/code/inference/nested-requirements-lock.txt
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .work/unified-cosmology/inference/.nested-venv/bin/python \
  studies/unified_cosmology/code/inference/late_tail.py
```

The source and results pin the physical likelihood, normalized data, full BAO
covariance, and design hashes. Original and refined nested sampling are
implemented in [late_nested.py](../code/inference/late_nested.py) and
[late_nested_refine.py](../code/inference/late_nested_refine.py). Generated
posterior arrays and checkpoints stay in the ignored working directory.
The [independent numerical audit](../code/inference/late_tail_validate.py)
checks identities, conditional-summary arithmetic, the weighted-CDF
comparison, and acceleration/jerk against finite derivatives of $E(z)$.
