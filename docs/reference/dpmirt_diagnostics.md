# Compute MCMC Diagnostics for a DPMirt Fit

Returns a structured list of diagnostic information including effective
sample sizes (ESS), WAIC, log-likelihood trace, timing, and DPM-specific
cluster diagnostics (number of clusters, alpha posterior).

## Usage

``` r
dpmirt_diagnostics(fit)

# S3 method for class 'dpmirt_diagnostics'
print(x, ...)
```

## Arguments

- fit:

  A `dpmirt_fit` object from
  [`dpmirt`](https://joonho112.github.io/DPMirt/reference/dpmirt.md).

- x:

  A `dpmirt_diagnostics` object.

- ...:

  Additional arguments (currently unused).

## Value

A `dpmirt_diagnostics` S3 object containing:

- ess:

  List of ESS vectors for items and persons.

- waic:

  WAIC value (if computed).

- loglik_trace:

  Log-likelihood trace vector.

- n_clusters:

  Posterior cluster counts (DPM only).

- alpha_summary:

  Alpha posterior summary (DPM only).

- compilation_time, sampling_time, total_time:

  Timing information.

## Details

Effective sample size (ESS) measures the number of effectively
independent draws from the posterior. Low ESS (\< 100) suggests poor
mixing and the need for longer chains or different samplers. For DPM
models, the cluster count trace is a key diagnostic — stable oscillation
indicates convergence, while monotonic trends suggest the chain has not
yet mixed.

## See also

[`dpmirt`](https://joonho112.github.io/DPMirt/reference/dpmirt.md),
[`dpmirt_compare`](https://joonho112.github.io/DPMirt/reference/dpmirt_compare.md),
[`plot.dpmirt_fit`](https://joonho112.github.io/DPMirt/reference/plot.dpmirt_fit.md)

Other diagnostics:
[`dpmirt_compare()`](https://joonho112.github.io/DPMirt/reference/dpmirt_compare.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
fit <- dpmirt(sim$response, model = "rasch", prior = "dpm",
              niter = 5000, nburnin = 1000, seed = 123)
diag <- dpmirt_diagnostics(fit)
print(diag)
} # }
```
