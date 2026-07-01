# DPMirt 0.2.0

DPMirt 0.2.0 is the overhaul release. It expands the package from the
initial development release into a release-candidate implementation with a
clearer model contract, richer posterior-draw metadata, broader diagnostics,
updated vignettes, and a rebuilt pkgdown site.

## Models

* Rasch, 2PL, and 3PL item response theory models.
* Parametric (Normal) and semiparametric (Dirichlet Process Mixture) latent
  trait priors.
* IRT and slope-intercept (SI) parameterizations for 2PL and 3PL.
* Three identification strategies: `constrained_item` (Rasch default),
  `constrained_ability`, and `unconstrained` with post-hoc rescaling
  (2PL/3PL default).
* 3PL models include item guessing parameters (`delta`), with posterior
  summaries, rescaled draw extraction, and plots where applicable.
* Some identification combinations are intentionally unavailable:
  `constrained_item` is rejected for 3PL models, and `constrained_ability`
  is rejected for DPM models.

## Workflow

* `dpmirt()` — all-in-one model fitting (spec → compile → sample → rescale).
* Modular pipeline: `dpmirt_spec()` → `dpmirt_compile()` → `dpmirt_sample()`
  → `dpmirt_rescale()`.
* `dpmirt_resume()` — continue sampling without recompilation while the
  compiled NIMBLE object is live in the current R session. Saved RDS objects
  cannot restore compiled pointers.
* `sampler_config` accepts an advanced function hook for NIMBLE sampler
  customization. List-based sampler schemas are reserved.
* `item_priors` is reserved; use `NULL` or `list()` for DPMirt's fixed item
  priors.

## Posterior Summaries

* Three estimators via `dpmirt_estimates()`:
  - **PM** (Posterior Mean) — optimal individual-level MSE.
  - **CB** (Constrained Bayes; Ghosh, 1992) — moment-matched to the
    posterior predictive distribution.
  - **GR** (Triple-Goal; Shen & Louis, 1998) — simultaneous estimation,
    ranking, and distributional recovery.
* `dpmirt_draws()` — extract posterior samples in matrix or long format.
  Current extraction is for rescaled draws only; raw draw extraction is
  reserved.
* `dpmirt_loss()` — evaluate MSEL, MSELR, KS, and custom loss metrics.

## Diagnostics

* `dpmirt_fit` objects carry `schema_version`, `chain_info`, `draw_index`,
  and `run_history` metadata.
* `dpmirt_diagnostics()` — ESS, optional chain-aware R-hat, trace summaries,
  WAIC provenance, timing, and DPM cluster diagnostics.
* `dpmirt_compare()` — WAIC-based model comparison across fits, including
  `waic_aggregation` provenance.
* `dpmirt_dp_density()` — posterior density estimation from the DP mixture.
  This can be computationally expensive because it reconstructs posterior DP
  measure samples through NIMBLE's public `modelValues` and
  `getSamplesDPmeasure()` workflow.
* `save_draws = FALSE` and `save_path` are reserved; DPMirt currently returns
  in-memory draw-retaining fit objects.

## Simulation

* `dpmirt_simulate()` — generate IRT data with Rasch, 2PL, and 3PL response
  models. 3PL fallback simulation includes Beta-distributed guessing
  parameters (`delta`) and does not use `target_rho`.
* IRTsimrel integration for Rasch/2PL reliability-targeted simulation via
  EQC calibration by default, with SAC calibration for MSEM targets; 3PL uses
  DPMirt's internal fallback simulator.

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
* Precomputed vignette fixtures follow an 11 `.rds` file compact fixture
  contract in `inst/extdata`, plus an `inst/extdata/README.md` note. Fit
  fixtures now ship as xz-compressed, evenly thinned
  250-draw `dpmirt_fit` objects with single-chain `chain-aware-v1` metadata
  and plot-ready DP-density summaries; session-bound compiled/raw NIMBLE
  storage is not included.
* 31 manual pages generated with roxygen2 and including the current reserved
  interfaces and diagnostic contracts.

## Compatibility Notes

* New and upgraded `dpmirt_fit` objects use the `chain-aware-v1` draw-storage
  schema. Fits created by older development versions may not contain
  `schema_version`, `chain_info`, `draw_index`, or `run_history` metadata.
* R-hat is reported only when at least two labeled chain/run streams are
  available. Single-chain fits should be assessed with ESS, trace summaries,
  posterior predictive checks, and substantive model diagnostics.
* `dpmirt_resume()` requires a live compiled NIMBLE object in the current R
  session. Saved RDS files cannot restore compiled external pointers.
* `dpmirt_draws()` currently returns retained rescaled posterior draws.
  Raw/unrescaled draw extraction and disk-backed draw storage are reserved for
  a future release.
* DP density recomputation no longer uses private NIMBLE namespace helpers.
  It loads retained DP draws through public `modelValues` accessors before
  calling `getSamplesDPmeasure()`.
* Transformed-scale DP density output for 2PL/3PL IRT and SI fits should be
  interpreted as a diagnostic summary. Full scale/Jacobian-adjusted density
  reconstruction remains a documented future enhancement.

# DPMirt 0.1.0

* Initial development release.
