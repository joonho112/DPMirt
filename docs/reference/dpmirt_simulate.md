# Simulate IRT response data

Generates binary response data under known IRT parameters. When
IRTsimrel is available, uses IRTsimrel calibration to target a marginal
reliability. The returned `reliability` field is the empirical KR-20
reliability of the generated response matrix; the calibration-level
achieved marginal reliability is returned as `achieved_rho`. Otherwise
falls back to Paganin-style simulation without reliability targeting.

## Usage

``` r
dpmirt_simulate(
  n_persons,
  n_items,
  model = c("rasch", "2pl", "3pl"),
  target_rho = 0.8,
  latent_shape = "normal",
  item_source = "parametric",
  reliability_metric = c("info", "msem"),
  latent_params = list(),
  item_params = list(),
  M = 10000L,
  seed = NULL,
  use_irtsimrel = TRUE,
  verbose = FALSE,
  ...
)

# S3 method for class 'dpmirt_sim'
print(x, ...)
```

## Arguments

- n_persons:

  Integer. Number of persons.

- n_items:

  Integer. Number of items.

- model:

  Character. IRT model type: `"rasch"`, `"2pl"`, or `"3pl"`. IRTsimrel
  reliability targeting is available for Rasch and 2PL; 3PL simulations
  use DPMirt's fallback generator.

- target_rho:

  Numeric in (0, 1). Target marginal reliability. Only used when
  IRTsimrel is available and `model` is `"rasch"` or `"2pl"`. Default
  0.8.

- latent_shape:

  Character. Shape of latent ability distribution. When using IRTsimrel,
  supports all 12 shapes: `"normal"`, `"bimodal"`, `"trimodal"`,
  `"multimodal"`, `"skew_pos"`, `"skew_neg"`, `"heavy_tail"`,
  `"light_tail"`, `"uniform"`, `"floor"`, `"ceiling"`, `"custom"`. For
  `"custom"`, supply
  `latent_params = list(mixture_spec = list(weights = ..., means = ..., sds = ...))`;
  custom shapes require IRTsimrel with `model = "rasch"` or
  `model = "2pl"`. Fallback supports: `"normal"`, `"bimodal"`,
  `"skewed"`.

- item_source:

  Character. Source for item parameters in IRTsimrel. One of
  `"parametric"` (default), `"irw"`, `"hierarchical"`, `"custom"`.

- reliability_metric:

  Character. Reliability metric for IRTsimrel calibration. `"info"`
  (default, average-information, recommended) uses EQC calibration;
  `"msem"` uses IRTsimrel's SAC calibration.

- latent_params:

  List. Additional parameters passed to
  [`IRTsimrel::sim_latentG()`](https://joonho112.github.io/IRTsimrel/reference/sim_latentG.html)
  (e.g., `list(shape_params = list(delta = 0.8))`). If
  `latent_shape = "custom"`, this must include `mixture_spec`.

- item_params:

  List. Additional parameters passed to
  [`IRTsimrel::sim_item_params()`](https://joonho112.github.io/IRTsimrel/reference/sim_item_params.html)
  (e.g., `list(discrimination_params = list(rho = -0.3))`).

- M:

  Integer. Quadrature sample size for EQC, and pre-calibration
  quadrature size for SAC. Default 10000.

- seed:

  Integer or NULL. Random seed for reproducibility.

- use_irtsimrel:

  Logical. If TRUE (default), attempt to use IRTsimrel for
  reliability-targeted simulation. Falls back to internal simulation if
  IRTsimrel is not installed.

- verbose:

  Logical. Print progress messages. Default FALSE.

- ...:

  Additional arguments (currently unused).

- x:

  A `dpmirt_sim` object.

## Value

A `dpmirt_sim` S3 object containing:

- response:

  N x I binary response matrix

- theta:

  True person abilities (length N)

- beta:

  True item difficulties (length I)

- lambda:

  True discriminations (length I for 2PL/3PL, NULL for Rasch)

- delta:

  True guessing parameters (length I for 3PL, NULL otherwise)

- n_persons, n_items, model:

  Simulation settings

- reliability:

  Empirical KR-20 reliability of the response matrix

- achieved_rho:

  IRTsimrel calibration achieved marginal reliability (NULL if fallback)

- target_rho:

  Requested target reliability (NULL if fallback)

- latent_shape:

  Distribution shape used

- eqc_result:

  IRTsimrel calibration result (NULL if fallback). The field name is
  retained for backward compatibility; the object class distinguishes
  EQC and SAC results.

- method:

  Character: "irtsimrel" or "fallback"

## See also

[`dpmirt`](https://joonho112.github.io/DPMirt/reference/dpmirt.md),
[`dpmirt_loss`](https://joonho112.github.io/DPMirt/reference/dpmirt_loss.md),
[`plot.dpmirt_sim`](https://joonho112.github.io/DPMirt/reference/plot.dpmirt_sim.md)

Other simulation:
[`dpmirt_loss()`](https://joonho112.github.io/DPMirt/reference/dpmirt_loss.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Simple fallback simulation
sim <- dpmirt_simulate(200, 25, model = "rasch", seed = 42)

# With IRTsimrel (reliability-targeted)
sim <- dpmirt_simulate(200, 25, model = "rasch",
                       target_rho = 0.85,
                       latent_shape = "bimodal",
                       seed = 42)

# Full simulation study workflow
sim <- dpmirt_simulate(200, 25, model = "rasch",
                       target_rho = 0.8,
                       latent_shape = "bimodal")
fit <- dpmirt(sim$response, model = "rasch", prior = "dpm")
est <- dpmirt_estimates(fit)
dpmirt_loss(est, true_theta = sim$theta, true_beta = sim$beta)
} # }
```
