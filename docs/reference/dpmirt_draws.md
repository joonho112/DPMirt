# Extract posterior draws from a DPMirt fit

Provides convenient access to posterior samples in matrix or long
format.

## Usage

``` r
dpmirt_draws(
  fit,
  vars = c("theta", "beta", "lambda", "delta"),
  format = c("matrix", "long"),
  use_rescaled = TRUE
)
```

## Arguments

- fit:

  A `dpmirt_fit` object.

- vars:

  Character vector of variables to extract. Options: `"theta"`,
  `"beta"`, `"lambda"`, `"delta"`.

- format:

  Character. Output format: `"matrix"` or `"long"`.

- use_rescaled:

  Logical. If TRUE (default), return rescaled samples.

## Value

A matrix (niter x N or niter x I) or long data.frame.

## Details

This function provides direct access to the posterior MCMC samples
stored in a `dpmirt_fit` object. The matrix format (default) is
efficient for computation, while the long format is convenient for
visualization with ggplot2.

## See also

[`dpmirt`](https://joonho112.github.io/DPMirt/reference/dpmirt.md),
[`dpmirt_estimates`](https://joonho112.github.io/DPMirt/reference/dpmirt_estimates.md)

Other estimation:
[`dpmirt_estimates()`](https://joonho112.github.io/DPMirt/reference/dpmirt_estimates.md),
[`dpmirt_rescale()`](https://joonho112.github.io/DPMirt/reference/dpmirt_rescale.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
fit <- dpmirt(sim$response, model = "rasch", prior = "normal",
              niter = 5000, nburnin = 1000, seed = 123)

# Extract theta draws as matrix (niter x N)
theta <- dpmirt_draws(fit, vars = "theta")
dim(theta$theta)

# Extract beta draws in long format
beta_long <- dpmirt_draws(fit, vars = "beta", format = "long")
head(beta_long$beta)
} # }
```
