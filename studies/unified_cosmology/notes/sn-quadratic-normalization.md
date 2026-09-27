# Supernova residuals and the published fit score

The released 1,820-row Dovekie distance likelihood reproduces the published SN-only flat-ΛCDM score when the same normalization is used. Minimizing over matter density gives Ωm = 0.3303168. The residual quadratic is **1631.42054**; integrating the flat magnitude offset adds **8.85267**, giving **1640.27320**, which rounds to the paper's **1640.3**. Its Table 10 reports Ωm = 0.330 ± 0.015. [Dovekie, version 3, Table 10](https://arxiv.org/html/2511.07517v3#S10.T10)

For residual vector r, released covariance C, and a vector of ones u, the quantities are

$$
A=u^TC^{-1}u,\qquad
Q=r^TC^{-1}r-\frac{(u^TC^{-1}r)^2}{A},\qquad
S=Q+\log\!\left(\frac{A}{2\pi}\right).
$$

The [pinned author likelihood](https://github.com/des-science/DES-SN5YR/blob/c9a4fcafc4cbd19bd750dee47fc76194a45c181f/5_COSMOLOGY/Dovekie_cosmosis_likelihood.py) returns −S/2. The additive term follows from integrating one common magnitude with a flat measure. This is the author's score convention; it is not an absolute Bayesian evidence with a normalized magnitude prior and all Gaussian determinant constants.

An independent adaptive distance integral agrees with 256-node Gaussian quadrature to 7.2×10⁻¹⁵ mag over five checked geometries. LU and Cholesky evaluations agree within 1.2×10⁻¹² in Q, and extraction of the author's numerical function confirms the offset. The checked model contains matter and Λ, omits radiation, and absorbs the arbitrary H0 = 70 normalization into the magnitude offset. [Executable check](../code/inference/sn_quadratic_normalization.py) · [Result and input hashes](../results/inference/sn-quadratic-normalization.json)

At a fixed geometry, Q has an N−1 Gaussian reference distribution **if** the released covariance is the true covariance of that Gaussian generating model: one unpenalized constant has been projected out. After fitting Ωm to the same observations, N−2 is a local linear or asymptotic reference, not an exact finite-sample result for an arbitrary nonlinear fit. This calculation does not report a fitted-model tail probability. A quadratic evaluated over the joint CMB–BAO–SN posterior also uses different geometries from this SN-only optimum, so it should not be compared directly with the published score.

There is a further distinction from the original supernova analysis: the paper reports **1,684 effective data points**, obtained by summing BBC BEAMS probabilities. That mixture-model count is not the dimension of the released 1,820-row Gaussian vector. Neither count can silently replace the other in a reference distribution. A small tail probability under the frozen Gaussian replication model is not, by itself, a p-value for the full BEAMS generator. [Dovekie, version 3, §10.5 and Appendix C](https://arxiv.org/html/2511.07517v3)

Checking the exact released cohort reveals a naming distinction. The 1,623 DES classifier probabilities `PROB_SNNV19`, plus 197 spectroscopic low-redshift objects, sum to **1684.2649**, which rounds to the printed count. The released BBC posterior diagnostic `PROBIA_BEAMS` instead sums to **1714.22219**; its per-row agreement with metadata `1-PROBCC_BEAMS` is within 5×10⁻⁶, the printed rounding. Identifiers, survey, redshifts, distances and normalized covariance order match exactly. This supplies a plausible interpretation of the printed count, not proof of the author's counting implementation or an implementation error. [Repeatable count check](../code/inference/dovekie_probability_counts.py) · [Results and hashes](../results/inference/dovekie-probability-counts.json)

This accounting agreement resolves the score comparison. It does not explain the low residual dispersion, establish a covariance pathology, or validate all selection and contamination assumptions. No covariance rescaling or change to the cosmological likelihood was made.
