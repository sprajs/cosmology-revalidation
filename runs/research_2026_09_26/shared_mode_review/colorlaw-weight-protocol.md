# Frozen COLORLAW weight-provenance sensitivity

Frozen after the expanded12 result and before executing this weight variant.
This is a targeted retrospective sensitivity requested by the parent, not a
new independent confirmatory analysis or evidence of the historical execution.

The local DES methods paper 2401.02945v2, section 6.3, explicitly gives the
COLORLAW covariance weight W_S^2=1/3. The archived fitopts.yml gives amplitude
0.3, and the inspected covariance code applies that amplitude before forming
the outer product, giving weight0.09. Preserve that original result.

For this one sensitivity, read the completed expanded12 discovery/validation
projected arrays. Column11 (zero-based) is the already0.3-weighted COLORLAW
response. Multiply it by sqrt(1/3)/0.3 in BOTH cohorts, once. Leave all other
mode columns, all data, object/epoch membership, nuisance projections, nominal
covariance, mode identities, and original observer coefficients unchanged.
Keep both the inherited observer prior and the already declared equal-total-
variance isotropic prior. Do not select whichever gives the better score.

Recompute each model's joint latent posterior from the43 discovery objects
with independent standard-normal prior coordinates, then integrate that single
shared posterior once over the1020 validation objects. Recomputing these
posteriors is necessary after changing the prior-scaled design; no validation
coefficient is estimated or inserted. Compute the12-coordinate null and the
two15-coordinate alternatives. Independently check low-rank joint-minus-
discovery scores against rank-aware SVD conditional predictive densities.
Recompute the same original frozen observer direction's centered matched
statistic and full deterministic shift score using the new null mean and
covariance. Report differences from all original expanded12 scores.

No native refits, new mode construction, prior tuning, covariance subtraction,
or changed physical correction follows. Neither the paper nor config/source
alone establishes the original executed weight. The question is only whether
this particular discrepancy materially changes the present finite-mode clue.
