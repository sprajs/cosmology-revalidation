# Literal native law identifiers: pre-fit clarification

The frozen protocol phrase “generation-consistent current option-99” means the literal namelist assignment `OPT_MWCOLORLAW=-99`. This is the historical approximate F99 implementation exposed under a negative option in the modern native executable. The implementation control is `OPT_MWCOLORLAW=99`, the current exact F99 law used in the real residual analysis. This resolves typography before simulation fits or residual scores; it does not change either intended hypothesis.

Sol independently traced the actual historical host-screen path: `snlc_sim.c:25508–25535` passes GENLC.RV/AV to genmag_SALT2; `genmag_SALT2.c:2813–2817` uses host attenuation when both are positive; `genmag_SEDtools.c:2401–2484` evaluates GALextinct at rest-frame wavelength with MWXT_SEDMODEL.OPT_COLORLAW. Thus P21 positive AV/RV invokes the passive screen and G10 -9 sentinels bypass it. The frozen optical-grid support labels remain unchanged.
