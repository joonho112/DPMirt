# Fit a Bayesian IRT model with DPMirt

The main entry point for DPMirt. Handles the full pipeline from data to
fitted model: specification, compilation, MCMC sampling, rescaling, and
diagnostics.

## Usage

``` r
dpmirt(
  data,
  model = c("rasch", "2pl", "3pl"),
  prior = c("normal", "dpm"),
  parameterization = c("irt", "si"),
  identification = NULL,
  niter = 10000L,
  nburnin = 2000L,
  thin = 1L,
  thin2 = NULL,
  nchains = 1L,
  seed = NULL,
  alpha_prior = NULL,
  mu_K = NULL,
  confidence = "medium",
  base_measure = list(s2_mu = 2, nu1 = 2.01, nu2 = 1.01),
  M = 50L,
  item_priors = list(),
  rescale = TRUE,
  compute_waic = TRUE,
  compute_dp_density = TRUE,
  save_draws = TRUE,
  save_path = NULL,
  sampler_config = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- data:

  A matrix or data.frame of binary (0/1) responses. Persons in rows,
  items in columns.

- model:

  Character. IRT model type: `"rasch"`, `"2pl"`, or `"3pl"`.

- prior:

  Character. Latent trait prior: `"normal"` or `"dpm"`.

- parameterization:

  Character. `"irt"` or `"si"`. Default `"irt"`.

- identification:

  Character or NULL. Identification strategy. If NULL, uses
  model-specific default.

- niter:

  Integer. Total MCMC iterations. Default 10000.

- nburnin:

  Integer. Burn-in iterations. Default 2000.

- thin:

  Integer. Thinning for main monitors. Default 1.

- thin2:

  Integer or NULL. Thinning for eta. NULL = same as thin.

- nchains:

  Integer. Number of chains. Default 1.

- seed:

  Integer or NULL. Random seed.

- alpha_prior:

  Alpha hyperprior for DPM. NULL (default Gamma(1,3) or auto-elicit if
  `mu_K` is set), a numeric vector `c(a, b)` for Gamma(a, b), or the
  output of
  [`dpmirt_alpha_prior`](https://joonho112.github.io/DPMirt/reference/dpmirt_alpha_prior.md).

- mu_K:

  Numeric or NULL. Expected number of clusters for automatic alpha prior
  elicitation via
  [`dpmirt_alpha_prior`](https://joonho112.github.io/DPMirt/reference/dpmirt_alpha_prior.md).
  Only used when `prior = "dpm"` and `alpha_prior = NULL`. Requires the
  DPprior package; falls back to Gamma(1, 3) if not installed.

- confidence:

  Character. DPprior confidence level for alpha elicitation: `"low"`,
  `"medium"` (default), or `"high"`. Only used when `mu_K` is set.

- base_measure:

  List of DPM base measure hyperparameters.

- M:

  Integer. Max clusters for CRP. Default 50.

- item_priors:

  List of custom item priors.

- rescale:

  Logical. Apply post-hoc rescaling. Default TRUE.

- compute_waic:

  Logical. Compute WAIC. Default TRUE.

- compute_dp_density:

  Logical. Compute DP density. Default TRUE for DPM.

- save_draws:

  Logical. Save full theta posterior. Default TRUE.

- save_path:

  Character or NULL. Path for disk-backed storage.

- sampler_config:

  Optional custom sampler configuration.

- verbose:

  Logical. Print progress. Default TRUE.

- ...:

  Additional arguments.

## Value

A `dpmirt_fit` S3 object.

## References

Paganin, S., Paciorek, C. J., Wehrhahn, C., Rodriguez, A., Rabe-Hesketh,
S., & de Valpine, P. (2023). Computational strategies and estimation
performance with Bayesian semiparametric item response theory models.
*Journal of Educational and Behavioral Statistics, 48*(2), 147–188.

## See also

[`dpmirt_spec`](https://joonho112.github.io/DPMirt/reference/dpmirt_spec.md),
[`dpmirt_compile`](https://joonho112.github.io/DPMirt/reference/dpmirt_compile.md),
[`dpmirt_sample`](https://joonho112.github.io/DPMirt/reference/dpmirt_sample.md),
[`dpmirt_estimates`](https://joonho112.github.io/DPMirt/reference/dpmirt_estimates.md),
[`dpmirt_simulate`](https://joonho112.github.io/DPMirt/reference/dpmirt_simulate.md)

Other model fitting:
[`dpmirt_compile()`](https://joonho112.github.io/DPMirt/reference/dpmirt_compile.md),
[`dpmirt_fit-methods`](https://joonho112.github.io/DPMirt/reference/dpmirt_fit-methods.md),
[`dpmirt_resume()`](https://joonho112.github.io/DPMirt/reference/dpmirt_resume.md),
[`dpmirt_sample()`](https://joonho112.github.io/DPMirt/reference/dpmirt_sample.md),
[`dpmirt_spec()`](https://joonho112.github.io/DPMirt/reference/dpmirt_spec.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# --- Simulate test data ---
sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)

# --- Rasch model with Normal prior ---
fit_normal <- dpmirt(sim$response, model = "rasch", prior = "normal",
                     niter = 5000, nburnin = 1000, seed = 123)
summary(fit_normal)
coef(fit_normal)

# --- Rasch model with DPM prior ---
fit_dpm <- dpmirt(sim$response, model = "rasch", prior = "dpm",
                  niter = 10000, nburnin = 3000, seed = 123)

# --- Compare Normal vs DPM ---
dpmirt_compare(fit_normal, fit_dpm)
} # }
```
