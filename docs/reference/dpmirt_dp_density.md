# Compute posterior density of the DP mixture

Samples from the posterior Dirichlet Process mixing distribution using
NIMBLE's
[`getSamplesDPmeasure()`](https://rdrr.io/pkg/nimble/man/getSamplesDPmeasure.html)
and evaluates the resulting mixture density on a grid. The density is
computed by summing weighted Normal components from the DP base measure.

## Usage

``` r
dpmirt_dp_density(
  fit,
  grid = seq(-6, 6, length.out = 500),
  credible_interval = 0.95,
  apply_rescaling = TRUE,
  verbose = TRUE,
  ...
)

# S3 method for class 'dpmirt_dp_density'
print(x, ...)
```

## Arguments

- fit:

  A `dpmirt_fit` object with `prior = "dpm"`.

- grid:

  Numeric vector. Grid points for density evaluation. Default:
  `seq(-6, 6, length.out = 500)`.

- credible_interval:

  Numeric. Width of the pointwise credible band. Default: 0.95 (i.e., 95
  percent band).

- apply_rescaling:

  Logical. If TRUE (default), shift the grid by the
  retained-draw-specific location shift from post-hoc rescaling. Only
  changes output for unconstrained Rasch/location-shift fits.

- verbose:

  Logical. Print progress messages. Default TRUE.

- ...:

  Additional arguments (currently unused).

- x:

  A `dpmirt_dp_density` object.

## Value

A list of class `dpmirt_dp_density` containing:

- grid:

  Numeric vector of evaluation points.

- density_mean:

  Numeric vector of posterior mean densities.

- density_lower:

  Numeric vector of lower credible band.

- density_upper:

  Numeric vector of upper credible band.

- density_samples:

  Matrix (`n_retained_draws` x `length(grid)`) of per-retained-draw
  densities (for custom summaries).

- dp_samples:

  List from
  [`getSamplesDPmeasure()`](https://rdrr.io/pkg/nimble/man/getSamplesDPmeasure.html)
  – each element is a matrix with columns (weights, means, variances).

- ci_level:

  The credible interval level used.

## Details

This function follows Paganin et al.'s (2023) workflow:

1.  Extract posterior samples for DP parameters (alpha, zi, muTilde,
    s2Tilde) from the fitted model.

2.  Reconstruct a NIMBLE model and MCMC with monitors set to only DP
    parameters.

3.  Populate the MCMC's sample storage with the posterior samples using
    NIMBLE modelValues accessors.

4.  Call
    [`getSamplesDPmeasure()`](https://rdrr.io/pkg/nimble/man/getSamplesDPmeasure.html)
    to compute stick-breaking weights and atoms for each posterior draw.

5.  Evaluate the mixture density
    `f(x|Gs) = sum_k w_k * phi(x; mu_k, s2_k)` for each posterior sample
    s and grid point x.

For Rasch models with unconstrained identification, a location shift
(mean(beta) per iteration) is applied so the density is on the rescaled
theta scale. For 2PL/3PL IRT and SI parameterizations, full
transformed-scale density reconstruction would also require scale and
Jacobian adjustments. The current implementation applies the
Rasch/location-shift contract, so DP densities for transformed-scale
2PL/3PL fits should be treated as diagnostic summaries rather than
definitive latent-density estimates.

## See also

[`dpmirt`](https://joonho112.github.io/DPMirt/reference/dpmirt.md),
[`dpmirt_plot_dp_density`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_dp_density.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- dpmirt_simulate(200, 20, model = "rasch",
                       latent_shape = "bimodal", seed = 42)
fit <- dpmirt(sim$response, model = "rasch", prior = "dpm",
              niter = 10000, nburnin = 3000, seed = 123)

# Compute DP density on default grid
dpd <- dpmirt_dp_density(fit)
print(dpd)

# Custom grid
dpd2 <- dpmirt_dp_density(fit, grid = seq(-4, 4, length.out = 200))
} # }
```
