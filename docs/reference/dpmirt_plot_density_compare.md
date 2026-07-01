# Plot Posterior Density vs Reference Distribution

Overlays the estimated posterior mean density with the N(0,1) reference
distribution and, when available, the fitted DP mixture density.
Requires ggplot2.

## Usage

``` r
dpmirt_plot_density_compare(fit, reference = c("normal"), ...)
```

## Arguments

- fit:

  A `dpmirt_fit` object from
  [`dpmirt`](https://joonho112.github.io/DPMirt/reference/dpmirt.md).

- reference:

  Character. Currently only `"normal"` is supported, which overlays
  N(0,1).

- ...:

  Currently unused.

## Value

A `ggplot` object.

## See also

[`plot.dpmirt_fit`](https://joonho112.github.io/DPMirt/reference/plot.dpmirt_fit.md),
[`dpmirt_plot_density`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_density.md)

Other visualization:
[`dpmirt_plot_caterpillar()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_caterpillar.md),
[`dpmirt_plot_clusters()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_clusters.md),
[`dpmirt_plot_density()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_density.md),
[`dpmirt_plot_dp_density()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_dp_density.md),
[`dpmirt_plot_icc()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_icc.md),
[`dpmirt_plot_info()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_info.md),
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

# Compare to N(0,1)
dpmirt_plot_density_compare(fit)

} # }
```
