# Populations, prediction and selection

Light-curve width improves held-out prediction in the selected DES sample. A selection-complete population inference is not yet established. Nominal SNNV19 reconstruction gives RMS probability error 0.1301 against a target of 0.005, and 89.54% membership agreement at the 0.999 threshold against a target of 98%. This prevents treating that reconstruction as the released selection model. Synthetic recovery and altered-population experiments probe the consequences rather than repairing this gap by assumption.

## Code and scope

Hierarchical likelihoods, conditional predictors, classifier preprocessing and comparisons, forward SNANA simulations, paired predictive scores, population-delay integrals and recovery tests.

The [source directory](code/) retains the original calculations. See the [execution guide](../README.md#execution) before running them: historical imports and relative paths are assembled into an isolated workspace, and some studies require additional data or native software. Source preservation does not certify that an incomplete experiment succeeded.

## Supporting findings

- [Classification reconstruction status](notes/classification.md)
- [Nominal DES SNNV19 reconstruction: failed release reproduction](notes/classifier.md)
- [Common reconstructed-classifier sensitivity in P21/G10 pilots](notes/common-classifier-residual.md)
- [Selected-sample conditional sampler repair (26 September 2026)](notes/conditional-repair.md)
- [Forward population and selection checks](notes/forward-discrimination.md)
- [Registered forward-comparison robustness continuation (26 September 2026)](notes/forward-robustness.md)
- [Independent hierarchical population model — specification before data comparisons](notes/generative-model.md)
- [Independent inference audit: preflight](notes/inference-audit.md)
- [Phase 2 literature and selection audit](notes/literature-selection.md)
- [Independent host-to-progenitor mapping investigation](notes/mapping-physics.md)
- [Coupled P21 noise-variant mechanism control](notes/noise-variant-control.md)
- [Synthetic recovery of a one-factor censored amplitude](notes/onefactor-recovery.md)
- [Fixed-vector wavelength and phase localization](notes/residual-localization-protocol.md)
- [Where the frozen residual pattern is measured](notes/residual-localization.md)
- [Disjoint same-catalogue P21/G10 residual holdout](notes/simulation-holdout-control.md)
- [P21/G10 native epoch-residual control](notes/simulation-residual-control.md)
- [Existing simulation epoch-residual control: feasibility inventory](notes/simulation-residual-feasibility.md)
