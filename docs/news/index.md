# Changelog

## DPMirt 0.1.0

*Initial release*

### Models

- Rasch, 2PL, and 3PL item response theory models.
- Parametric (Normal) and semiparametric (Dirichlet Process Mixture)
  latent trait priors.
- IRT and slope-intercept (SI) parameterizations for 2PL and 3PL.
- Three identification strategies: `constrained_item` (Rasch default),
  `constrained_ability`, and `unconstrained` with post-hoc rescaling
  (2PL/3PL default).

### Workflow

- [`dpmirt()`](https://joonho112.github.io/DPMirt/reference/dpmirt.md) —
  all-in-one model fitting (spec → compile → sample → rescale).
- Modular pipeline:
  [`dpmirt_spec()`](https://joonho112.github.io/DPMirt/reference/dpmirt_spec.md)
  →
  [`dpmirt_compile()`](https://joonho112.github.io/DPMirt/reference/dpmirt_compile.md)
  →
  [`dpmirt_sample()`](https://joonho112.github.io/DPMirt/reference/dpmirt_sample.md)
  →
  [`dpmirt_rescale()`](https://joonho112.github.io/DPMirt/reference/dpmirt_rescale.md).
- [`dpmirt_resume()`](https://joonho112.github.io/DPMirt/reference/dpmirt_resume.md)
  — continue sampling from a compiled model without recompilation.

### Posterior Summaries

- Three estimators via
  [`dpmirt_estimates()`](https://joonho112.github.io/DPMirt/reference/dpmirt_estimates.md):
  - **PM** (Posterior Mean) — optimal individual-level MSE.
  - **CB** (Constrained Bayes; Ghosh, 1992) — moment-matched to the
    posterior predictive distribution.
  - **GR** (Triple-Goal; Shen & Louis, 1998) — simultaneous estimation,
    ranking, and distributional recovery.
- [`dpmirt_draws()`](https://joonho112.github.io/DPMirt/reference/dpmirt_draws.md)
  — extract posterior samples in matrix or long format.
- [`dpmirt_loss()`](https://joonho112.github.io/DPMirt/reference/dpmirt_loss.md)
  — evaluate MSEL, MSELR, KS, and custom loss metrics.

### Diagnostics

- [`dpmirt_diagnostics()`](https://joonho112.github.io/DPMirt/reference/dpmirt_diagnostics.md)
  — ESS, R-hat, trace summaries, and WAIC.
- [`dpmirt_compare()`](https://joonho112.github.io/DPMirt/reference/dpmirt_compare.md)
  — WAIC-based model comparison across fits.
- [`dpmirt_dp_density()`](https://joonho112.github.io/DPMirt/reference/dpmirt_dp_density.md)
  — posterior density estimation from the DP mixture.

### Simulation

- [`dpmirt_simulate()`](https://joonho112.github.io/DPMirt/reference/dpmirt_simulate.md)
  — generate IRT data with 12 latent distribution shapes (normal,
  bimodal, skewed, heavy-tailed, etc.).
- IRTsimrel integration for reliability-targeted simulation via EQC
  calibration.

### Prior Elicitation

- [`dpmirt_alpha_prior()`](https://joonho112.github.io/DPMirt/reference/dpmirt_alpha_prior.md)
  — principled concentration-parameter selection via DPprior, with
  graceful fallback to Gamma(1, 3).

### Visualization

- S3 [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods
  for `dpmirt_fit`, `dpmirt_estimates`, and `dpmirt_sim` objects.
- 12 standalone ggplot2 functions:
  [`dpmirt_plot_trace()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_trace.md),
  [`dpmirt_plot_density()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_density.md),
  [`dpmirt_plot_caterpillar()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_caterpillar.md),
  [`dpmirt_plot_items()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_items.md),
  [`dpmirt_plot_icc()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_icc.md),
  [`dpmirt_plot_info()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_info.md),
  [`dpmirt_plot_dp_density()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_dp_density.md),
  [`dpmirt_plot_clusters()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_clusters.md),
  [`dpmirt_plot_wright_map()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_wright_map.md),
  [`dpmirt_plot_parameter_trace()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_parameter_trace.md),
  [`dpmirt_plot_density_compare()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_density_compare.md),
  [`dpmirt_plot_pp_check()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_pp_check.md).

### Documentation

- Eight vignettes covering quick start, models and workflow, posterior
  summaries, prior elicitation, simulation studies, mathematical
  foundations, and NIMBLE internals.
- 31 manual pages with full roxygen2 documentation.
