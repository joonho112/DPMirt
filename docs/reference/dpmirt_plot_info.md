# Plot Test Information Function

Displays the test information function (TIF), which shows the precision
of measurement across the ability range. Higher information indicates
more precise measurement. Requires ggplot2.

## Usage

``` r
dpmirt_plot_info(
  fit,
  theta_range = c(-4, 4),
  n_points = 201,
  show_items = FALSE,
  show_density = TRUE,
  ...
)
```

## Arguments

- fit:

  A `dpmirt_fit` object from
  [`dpmirt`](https://joonho112.github.io/DPMirt/reference/dpmirt.md).

- theta_range:

  Numeric vector of length 2. Default: `c(-4, 4)`.

- n_points:

  Integer. Number of grid points. Default: 201.

- show_items:

  Logical. Show individual item information curves. Default FALSE.

- show_density:

  Logical. Overlay person ability density (scaled). Default TRUE.

- ...:

  Currently unused.

## Value

A `ggplot` object.

## Details

Fisher information for each item is:

**Rasch**: \\I_j(\theta) = P(1-P)\\

**2PL**: \\I_j(\theta) = \lambda^2 P(1-P)\\

**3PL**: \\I_j(\theta) = \lambda^2 \frac{(P - \delta)^2}{(1 - \delta)^2
P(1 - P)}\\

The test information function is the sum: \\I(\theta) = \sum_j
I_j(\theta)\\.

## See also

[`plot.dpmirt_fit`](https://joonho112.github.io/DPMirt/reference/plot.dpmirt_fit.md),
[`dpmirt_plot_icc`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_icc.md)

Other visualization:
[`dpmirt_plot_caterpillar()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_caterpillar.md),
[`dpmirt_plot_clusters()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_clusters.md),
[`dpmirt_plot_density()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_density.md),
[`dpmirt_plot_density_compare()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_density_compare.md),
[`dpmirt_plot_dp_density()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_dp_density.md),
[`dpmirt_plot_icc()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_icc.md),
[`dpmirt_plot_items()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_items.md),
[`dpmirt_plot_parameter_trace()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_parameter_trace.md),
[`dpmirt_plot_pp_check()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_pp_check.md),
[`dpmirt_plot_trace()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_trace.md),
[`dpmirt_plot_wright_map()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_wright_map.md),
[`plot.dpmirt_estimates()`](https://joonho112.github.io/DPMirt/reference/plot.dpmirt_estimates.md),
[`plot.dpmirt_fit()`](https://joonho112.github.io/DPMirt/reference/plot.dpmirt_fit.md),
[`plot.dpmirt_sim()`](https://joonho112.github.io/DPMirt/reference/plot.dpmirt_sim.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
fit <- dpmirt(sim$response, model = "rasch", prior = "normal",
              niter = 5000, nburnin = 1000, seed = 123)
dpmirt_plot_info(fit)
} # }
```
