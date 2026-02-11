# DPMirt 0.1.0

*Initial release*

## Models

* Rasch, 2PL, and 3PL item response theory models.
* Parametric (Normal) and semiparametric (Dirichlet Process Mixture) latent
  trait priors.
* IRT and slope-intercept (SI) parameterizations for 2PL and 3PL.
* Three identification strategies: `constrained_item` (Rasch default),
  `constrained_ability`, and `unconstrained` with post-hoc rescaling
  (2PL/3PL default).

## Workflow

* `dpmirt()` — all-in-one model fitting (spec → compile → sample → rescale).
* Modular pipeline: `dpmirt_spec()` → `dpmirt_compile()` → `dpmirt_sample()`
  → `dpmirt_rescale()`.
* `dpmirt_resume()` — continue sampling from a compiled model without
  recompilation.

## Posterior Summaries

* Three estimators via `dpmirt_estimates()`:
  - **PM** (Posterior Mean) — optimal individual-level MSE.
  - **CB** (Constrained Bayes; Ghosh, 1992) — moment-matched to the
    posterior predictive distribution.
  - **GR** (Triple-Goal; Shen & Louis, 1998) — simultaneous estimation,
    ranking, and distributional recovery.
* `dpmirt_draws()` — extract posterior samples in matrix or long format.
* `dpmirt_loss()` — evaluate MSEL, MSELR, KS, and custom loss metrics.

## Diagnostics

* `dpmirt_diagnostics()` — ESS, R-hat, trace summaries, and WAIC.
* `dpmirt_compare()` — WAIC-based model comparison across fits.
* `dpmirt_dp_density()` — posterior density estimation from the DP mixture.

## Simulation

* `dpmirt_simulate()` — generate IRT data with 12 latent distribution
  shapes (normal, bimodal, skewed, heavy-tailed, etc.).
* IRTsimrel integration for reliability-targeted simulation via EQC
  calibration.

## Prior Elicitation

* `dpmirt_alpha_prior()` — principled concentration-parameter selection
  via DPprior, with graceful fallback to Gamma(1, 3).

## Visualization

* S3 `plot()` methods for `dpmirt_fit`, `dpmirt_estimates`, and
  `dpmirt_sim` objects.
* 12 standalone ggplot2 functions: `dpmirt_plot_trace()`,
  `dpmirt_plot_density()`, `dpmirt_plot_caterpillar()`,
  `dpmirt_plot_items()`, `dpmirt_plot_icc()`, `dpmirt_plot_info()`,
  `dpmirt_plot_dp_density()`, `dpmirt_plot_clusters()`,
  `dpmirt_plot_wright_map()`, `dpmirt_plot_parameter_trace()`,
  `dpmirt_plot_density_compare()`, `dpmirt_plot_pp_check()`.

## Documentation

* Eight vignettes covering quick start, models and workflow, posterior
  summaries, prior elicitation, simulation studies, mathematical
  foundations, and NIMBLE internals.
* 31 manual pages with full roxygen2 documentation.
