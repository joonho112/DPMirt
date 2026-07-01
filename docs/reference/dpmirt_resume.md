# Resume MCMC sampling from a previous run

Continues MCMC sampling without recompilation, using the NIMBLE
`Cmcmc$run(niter, reset = FALSE)` pattern.

## Usage

``` r
dpmirt_resume(fit_or_compiled, niter_more, reset = FALSE, verbose = TRUE, ...)
```

## Arguments

- fit_or_compiled:

  A `dpmirt_samples`, `dpmirt_fit`, or `dpmirt_compiled` object.

- niter_more:

  Integer. Number of additional iterations.

- reset:

  Logical. FALSE to continue from current state (default), TRUE to
  restart.

- verbose:

  Logical. Print progress messages.

- ...:

  Additional arguments.

## Value

A `dpmirt_samples` object with extended samples. Pass the returned
object to
[`dpmirt_rescale`](https://joonho112.github.io/DPMirt/reference/dpmirt_rescale.md)
before using fit methods such as
[`summary()`](https://rdrr.io/r/base/summary.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html), or
[`dpmirt_estimates()`](https://joonho112.github.io/DPMirt/reference/dpmirt_estimates.md).

## Details

Resume requires the original live compiled NIMBLE object; saved RDS
objects cannot restore compiled external pointers. Multi-chain fitted
objects cannot be continued with `reset = FALSE`; use `reset = TRUE` for
a fresh single-chain restart or rerun
[`dpmirt`](https://joonho112.github.io/DPMirt/reference/dpmirt.md).

## See also

[`dpmirt_sample`](https://joonho112.github.io/DPMirt/reference/dpmirt_sample.md),
[`dpmirt`](https://joonho112.github.io/DPMirt/reference/dpmirt.md)

Other model fitting:
[`dpmirt()`](https://joonho112.github.io/DPMirt/reference/dpmirt.md),
[`dpmirt_compile()`](https://joonho112.github.io/DPMirt/reference/dpmirt_compile.md),
[`dpmirt_fit-methods`](https://joonho112.github.io/DPMirt/reference/dpmirt_fit-methods.md),
[`dpmirt_sample()`](https://joonho112.github.io/DPMirt/reference/dpmirt_sample.md),
[`dpmirt_spec()`](https://joonho112.github.io/DPMirt/reference/dpmirt_spec.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Continue from a previous fit for more iterations
fit <- dpmirt(sim$response, model = "rasch", prior = "normal",
              niter = 5000, nburnin = 1000)
resumed <- dpmirt_resume(fit, niter_more = 5000)
fit2 <- dpmirt_rescale(resumed)
summary(fit2)
} # }
```
