# Independently specified spectral alternatives: acquisition audit

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../README.md) for current execution status.

The broad SED construction establishes ambiguity but supplies no empirical
population prior. Two published alternatives provide more constrained next
tests. Neither has been silently substituted into the current analysis.

**SALT3+.** The study adds a measured light-curve component beyond stretch and
colour, with coherent colour/secondary-maximum variation. It reports a
Hubble-residual association when that component is omitted, while finding no
evidence of bias in the cosmologies it examined. This motivates testing the
published component, rather than inferring that our arbitrary smooth family is
present. Its ZTF plus earlier SALT training and shared calibration must be
accounted for when assessing independence.
[Kenworthy et al., 2502.09713v1](https://arxiv.org/abs/2502.09713v1).

The primary PDF was acquired and hashed. `pdftotext` reported an xref repair
warning but produced the 18-page text; no figure was digitized. The current
public SALTShaker tree was pinned at
`5183c2dcd8f0a6fba61e5a2137bcf186120ef43f`. It includes a
`SALT3PlusSource` implementation expecting `salt3_template_2.dat`, but its
recursive file listing did not identify the paper's trained M2 surface.
That source clips its combined spectral flux to zero, so adopting it requires
an explicit implementation comparison; copying only a derivative would not
reproduce its mean model. No fetched code was installed or executed.
[Pinned implementation](https://github.com/djones1040/SALTShaker/blob/5183c2dcd8f0a6fba61e5a2137bcf186120ef43f/saltshaker/util/salt3plussncosmo.py).

**SALT3-UV.** A newer study incorporates well-calibrated HST UV spectra and
reports tentative low/high-redshift spectral differences. Its authors also
identify auxiliary-photometry calibration limitations and problematic model
support at the shortest wavelengths. These are reasons to test the exact
released alternatives and observational constraints, not to treat its reported
difference as a DES correction or as a prior validating our SED amplitudes.
[Wang et al., 2512.25064v2](https://arxiv.org/abs/2512.25064v2).

The primary full HTML/text is preserved. Its identified Zenodo link is a
**sncosmo software citation**, not the spectral dataset. Neither the inspected
paper links, SALTShaker tree, nor the author's inspected public repository
listing yielded the exact low/high-z trained surfaces. This is a bounded
acquisition gap, not a claim that the assets do not exist. No model was
reconstructed by digitizing a plot or replacing missing inputs with a generic
SALT template. The new [throughput localization](../../population_inference/notes/residual-localization.md)
tests the narrower UV-confinement hypothesis using existing exact inputs.

The same read-only refresh checked the public calibration repositories.
DES-SN5YR HEAD remains `c9a4fcafc4cbd19bd750dee47fc76194a45c181f`;
Dovekie HEAD is `67b8b28c9f9033732524c71e271dddeab3cec097`.
The latter's live `dovekie.py` and `DOVEKIE_DEFS.yml` match the inspected local
bytes exactly. The previously documented missing execution mapping has not
been resolved by a changed version of those two files.

The acquisition folder (`sources/updates/2026-09-26-spectral-models`; historical local reference)
retains primary pages, URLs, GitHub metadata, bounded search results, code and
hashes. No contact with authors or other external communication was made.
