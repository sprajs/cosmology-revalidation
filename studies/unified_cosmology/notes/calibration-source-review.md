# What the public calibration covariance does and does not establish

The released Pantheon+SH0ES statistical covariance contains shared terms between distinct supernovae in the same Cepheid host. For four such hosts, those terms are approximately twice the independently reconstructed variance of the host distance. This pattern is independently verified, but its production-code origin is unresolved. It is not yet an established error in the release.

| Host | Distinct supernovae | Released shared STATONLY covariance (mag²) | Reconstructed SN-free host variance (mag²) | Ratio |
|---|---|---:|---:|---:|
| NGC 1448 | 2001el, 2021pit | 0.00265025 | 0.001343615 | 1.97248 |
| NGC 3147 | 1997bq, 2008fv, 2021hpr | 0.05463299 | 0.027345159 | 1.99790 |
| NGC 5468 | 1999cp, 2002cr | 0.01095859 | 0.005497042 | 1.99354 |
| NGC 5643 | 2013aa, 2017cbv | 0.00538959 | 0.002722283 | 1.97981 |

The check preserves every released light-curve row and compares different physical supernovae, excluding duplicate observations of the same event. The [source and calculation record](../results/distance_ladder/calibration-source-review.json) stores the exact row pairs, matrix/input hashes, source URLs and short supporting quotations. The host variances come from the separately validated [Cepheid-only factor](distance-ladder.md), not from fitting the supernova magnitudes.

Brout et al. write the combined residual covariance as the sum of an SN covariance and one Cepheid covariance. Their statistical prescription also distinguishes repeated light curves of the **same** supernova from different events. The paper supplies the likelihood equation, but not the complete numerical host-to-light-curve embedding that generated these particular release files. [Brout et al., equations 5 and 14–15](https://arxiv.org/html/2202.04077v2#S2.SS3)

Riess et al. describe the 37-host covariance as derived from the first two rungs without using SN measurements. Their Table 6 distinguishes approximate distances excluding an individual host's SN measurements from distances excluding all SN measurements. That supports comparing the release with an SN-free host factor, but does not prove that our reconstructed factor is the exact numerical version used to assemble the public Pantheon+ covariance. [Riess et al., Table 6, Figure 16 and §5.2](https://arxiv.org/pdf/2112.04510v3)

Scolnic et al. find no evidence for a shared intrinsic distance residual among different SNe in the same host, while treating repeated observations of one SN as sharing its intrinsic scatter. The former statement concerns SN intrinsic scatter; it does not make a common Cepheid distance independent for each sibling. [Scolnic et al., §§3.1–3.2](https://arxiv.org/pdf/2112.03863v2)

For a single uncertain host distance (u), independent of the SN measurements, ordinary covariance propagation gives

$$
\operatorname{Cov}(m_i-u,m_j-u)
=\operatorname{Cov}(m_i,m_j)+\operatorname{Var}(u).
$$

Thus, having two siblings does not by itself require two copies of the host variance. An additional term, a different input version, or a construction issue could produce a different released value. The numerical pattern alone cannot select among them.

The inspected [pipeline configuration](https://github.com/PantheonPlusSH0ES/DataRelease/blob/c447f0fea703fcd0fff57de5000947b5ca81286b/Pantheon%2B_Data/7_PIPPIN/PPLUS.yml#L4126) refers to an external duplicate-scatter covariance and several velocity covariances. Its active blocks do not specify the later Cepheid embedding, and the named duplicate file and embedding program are absent from the pinned release tree. Both public covariance file-path histories contain a single upload on 5 July 2022, so no later tracked edit at these paths explains the contrast. This does not establish the history of the unpublished construction inputs.

The [release documentation](https://github.com/PantheonPlusSH0ES/DataRelease/blob/c447f0fea703fcd0fff57de5000947b5ca81286b/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/README) and an [author clarification](https://github.com/PantheonPlusSH0ES/DataRelease/issues/6#issuecomment-1385104075) confirm that Cepheid uncertainty is already present in the covariance used for absolute calibration. They do not resolve this sibling pattern. Consequently, adding a separate Cepheid likelihood to that released covariance would reuse calibration information; subtracting our host factor would not establish an exact SN-only covariance. The released joint likelihood remains usable as its own conditional model. Reconstructing a different joint model requires the exact component covariance and embedding inputs, or a complete reconstruction from the underlying measurements.
