# Plot Item Characteristic Curves

Displays the Item Characteristic Curves (ICCs) for all items, showing
the probability of a correct response as a function of ability. Requires
ggplot2.

## Usage

``` r
dpmirt_plot_icc(fit, items = NULL, theta_range = c(-4, 4), n_points = 201, ...)
```

## Arguments

- fit:

  A `dpmirt_fit` object from
  [`dpmirt`](https://joonho112.github.io/DPMirt/reference/dpmirt.md).

- items:

  Integer vector of item indices to plot. Default: all items (up to 10).

- theta_range:

  Numeric vector of length 2 for the ability axis range. Default:
  `c(-4, 4)`.

- n_points:

  Integer. Number of grid points. Default: 201.

- ...:

  Currently unused.

## Value

A `ggplot` object.

## Details

The ICC gives the probability of a correct response at each ability
level:

**Rasch**: \\P(\theta) = \mathrm{logistic}(\theta - \beta)\\

**2PL**: \\P(\theta) = \mathrm{logistic}(\lambda(\theta - \beta))\\

**3PL**: \\P(\theta) = \delta + (1 - \delta) \cdot
\mathrm{logistic}(\lambda(\theta - \beta))\\

## See also

[`plot.dpmirt_fit`](https://joonho112.github.io/DPMirt/reference/plot.dpmirt_fit.md),
[`dpmirt_plot_info`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_info.md)

Other visualization:
[`dpmirt_plot_caterpillar()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_caterpillar.md),
[`dpmirt_plot_clusters()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_clusters.md),
[`dpmirt_plot_density()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_density.md),
[`dpmirt_plot_density_compare()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_density_compare.md),
[`dpmirt_plot_dp_density()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_dp_density.md),
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
dpmirt_plot_icc(fit)
dpmirt_plot_icc(fit, items = 1:5)
} # }
```
