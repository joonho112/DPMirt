# Plot dpmirt_estimates object

Plot dpmirt_estimates object

## Usage

``` r
# S3 method for class 'dpmirt_estimates'
plot(x, type = c("estimates", "shrinkage"), param = c("theta", "beta"), ...)
```

## Arguments

- x:

  A `dpmirt_estimates` object.

- type:

  Character. Plot type: `"estimates"` (caterpillar of PM/CB/GR) or
  `"shrinkage"` (PM vs CB scatter).

- param:

  Character. Which parameter: `"theta"` or `"beta"`.

- ...:

  Additional graphical parameters.

## See also

[`dpmirt_estimates`](https://joonho112.github.io/DPMirt/reference/dpmirt_estimates.md),
[`dpmirt_plot_caterpillar`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_caterpillar.md)

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
[`plot.dpmirt_fit()`](https://joonho112.github.io/DPMirt/reference/plot.dpmirt_fit.md),
[`plot.dpmirt_sim()`](https://joonho112.github.io/DPMirt/reference/plot.dpmirt_sim.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
fit <- dpmirt(sim$response, model = "rasch", prior = "normal",
              niter = 5000, nburnin = 1000, seed = 123)
est <- dpmirt_estimates(fit)

# PM/CB/GR caterpillar plot
plot(est, type = "estimates")

# Shrinkage plot (PM vs CB)
plot(est, type = "shrinkage")
} # }
```
