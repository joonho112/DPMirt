# Create a DPMirt model specification

Constructs a complete specification for a Bayesian IRT model, including
NIMBLE model code, constants, data, initial values, and monitor
configuration. This is the first step in the step-by-step workflow.

## Usage

``` r
dpmirt_spec(
  data,
  model = c("rasch", "2pl", "3pl"),
  prior = c("normal", "dpm"),
  parameterization = c("irt", "si"),
  identification = NULL,
  alpha_prior = NULL,
  base_measure = list(s2_mu = 2, nu1 = 2.01, nu2 = 1.01),
  item_priors = list(),
  M = 50L,
  data_format = c("auto", "matrix", "long"),
  ...
)

# S3 method for class 'dpmirt_spec'
print(x, ...)
```

## Arguments

- data:

  A matrix or data.frame of binary (0/1) responses with persons in rows
  and items in columns. Can also be a long-format data.frame with
  columns for person, item, and response.

- model:

  Character. IRT model type: `"rasch"`, `"2pl"`, or `"3pl"`.

- prior:

  Character. Latent trait prior: `"normal"` (parametric) or `"dpm"`
  (Dirichlet Process Mixture).

- parameterization:

  Character. `"irt"` for standard IRT or `"si"` for slope-intercept.
  Only relevant for 2PL/3PL.

- identification:

  Character or NULL. Identification strategy: `"constrained_item"`,
  `"constrained_ability"`, or `"unconstrained"`. If NULL, uses
  model-specific default (constrained_item for Rasch, unconstrained for
  2PL/3PL).

- alpha_prior:

  Alpha hyperprior specification. NULL for default (auto-elicit or
  Gamma(1,3)), a numeric vector c(a, b) for Gamma(a, b), or a
  DPprior_fit object. Only used when prior = "dpm".

- base_measure:

  List with DPM base measure hyperparameters: s2_mu, nu1, nu2. Defaults
  from Paganin et al. (2023).

- item_priors:

  List of custom item priors (advanced use).

- M:

  Integer. Maximum number of clusters for CRP truncation. Only used when
  prior = "dpm".

- data_format:

  Character. `"auto"`, `"matrix"`, or `"long"`.

- ...:

  Additional arguments (currently unused).

- x:

  A `dpmirt_spec` object.

## Value

A `dpmirt_spec` S3 object containing:

- code:

  A `nimbleCode` object.

- constants:

  List of constants (N, I, M, etc.).

- data:

  List with the response data.

- inits:

  List of initial values.

- monitors:

  Character vector of parameters to track.

- monitors2:

  Character vector of parameters for thinned monitoring.

- config:

  List of all model configuration options.

## See also

[`dpmirt`](https://joonho112.github.io/DPMirt/reference/dpmirt.md),
[`dpmirt_compile`](https://joonho112.github.io/DPMirt/reference/dpmirt_compile.md)

Other model fitting:
[`dpmirt()`](https://joonho112.github.io/DPMirt/reference/dpmirt.md),
[`dpmirt_compile()`](https://joonho112.github.io/DPMirt/reference/dpmirt_compile.md),
[`dpmirt_fit-methods`](https://joonho112.github.io/DPMirt/reference/dpmirt_fit-methods.md),
[`dpmirt_resume()`](https://joonho112.github.io/DPMirt/reference/dpmirt_resume.md),
[`dpmirt_sample()`](https://joonho112.github.io/DPMirt/reference/dpmirt_sample.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)

# Rasch-Normal specification
spec <- dpmirt_spec(sim$response, model = "rasch", prior = "normal")
print(spec)

# Rasch-DPM with custom alpha prior
spec_dpm <- dpmirt_spec(sim$response, model = "rasch", prior = "dpm",
                        alpha_prior = c(1, 3))

# 2PL specification
sim2 <- dpmirt_simulate(300, 25, model = "2pl", seed = 42)
spec_2pl <- dpmirt_spec(sim2$response, model = "2pl", prior = "normal",
                        parameterization = "irt")
} # }
```
