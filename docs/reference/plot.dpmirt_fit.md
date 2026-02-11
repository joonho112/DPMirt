# Plot a DPMirt fit object

Produces visualizations for fitted IRT models. Supports 12 plot types
with automatic selection between base R and ggplot2 backends.

## Usage

``` r
# S3 method for class 'dpmirt_fit'
plot(
  x,
  type = c("density", "items", "trace", "clusters", "dp_density", "icc", "wright_map",
    "parameter_trace", "caterpillar", "density_compare", "info", "pp_check"),
  engine = c("auto", "base", "ggplot2"),
  ...
)
```

## Arguments

- x:

  A `dpmirt_fit` object.

- type:

  Character. Plot type. One of:

  `"density"`

  :   Kernel density of posterior mean theta.

  `"items"`

  :   Item difficulty estimates with error bars.

  `"trace"`

  :   Log-likelihood MCMC trace.

  `"clusters"`

  :   Cluster count trace and histogram (DPM only).

  `"dp_density"`

  :   DP mixture density with credible band (DPM only).

  `"icc"`

  :   Item Characteristic Curves.

  `"wright_map"`

  :   Person-Item map (Wright map).

  `"parameter_trace"`

  :   Individual parameter MCMC traces.

  `"caterpillar"`

  :   Sorted estimates with credible intervals.

  `"density_compare"`

  :   Posterior density vs reference overlay.

  `"info"`

  :   Test Information Function.

  `"pp_check"`

  :   Posterior predictive check.

- engine:

  Character. Plotting backend: `"auto"` (default) uses ggplot2 if
  available, `"base"` forces base R, `"ggplot2"` requires ggplot2.

- ...:

  Additional arguments passed to the specific plotting function.

## Value

Invisibly returns the plot object (ggplot) or NULL (base R).

## See also

[`dpmirt_plot_density`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_density.md),
[`dpmirt_plot_items`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_items.md),
[`dpmirt_plot_trace`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_trace.md),
[`dpmirt_plot_icc`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_icc.md)

Other visualization:
[`dpmirt_plot_caterpillar()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_caterpillar.md),
[`dpmirt_plot_clusters()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_clusters.md),
[`dpmirt_plot_density()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_density.md),
[`dpmirt_plot_density_compare()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_density_compare.md),
[`dpmirt_plot_dp_density()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_dp_density.md),
[`dpmirt_plot_icc()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_icc.md),
[`dpmirt_plot_info()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_info.md),
[`dpmirt_plot_items()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_items.md),
[`dpmirt_plot_parameter_trace()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_parameter_trace.md),
[`dpmirt_plot_pp_check()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_pp_check.md),
[`dpmirt_plot_trace()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_trace.md),
[`dpmirt_plot_wright_map()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_wright_map.md),
[`plot.dpmirt_estimates()`](https://joonho112.github.io/DPMirt/reference/plot.dpmirt_estimates.md),
[`plot.dpmirt_sim()`](https://joonho112.github.io/DPMirt/reference/plot.dpmirt_sim.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
fit <- dpmirt(sim$response, model = "rasch", prior = "normal",
              niter = 5000, nburnin = 1000, seed = 123)

# Theta density
plot(fit, type = "density")

# Item difficulty estimates
plot(fit, type = "items")

# MCMC trace
plot(fit, type = "trace")

# Force base R backend
plot(fit, type = "density", engine = "base")
} # }
```
