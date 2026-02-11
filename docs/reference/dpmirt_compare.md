# Compare DPMirt models using information criteria

Compare DPMirt models using information criteria

## Usage

``` r
dpmirt_compare(..., criterion = "waic")
```

## Arguments

- ...:

  Two or more `dpmirt_fit` objects.

- criterion:

  Character. Comparison criterion. Default "waic".

## Value

A data.frame ranking models by the criterion.

## See also

[`dpmirt`](https://joonho112.github.io/DPMirt/reference/dpmirt.md),
[`dpmirt_diagnostics`](https://joonho112.github.io/DPMirt/reference/dpmirt_diagnostics.md)

Other diagnostics:
[`dpmirt_diagnostics()`](https://joonho112.github.io/DPMirt/reference/dpmirt_diagnostics.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
fit1 <- dpmirt(sim$response, model = "rasch", prior = "normal",
               niter = 5000, nburnin = 1000, seed = 123)
fit2 <- dpmirt(sim$response, model = "rasch", prior = "dpm",
               niter = 5000, nburnin = 1000, seed = 123)

# Compare via WAIC
comp <- dpmirt_compare(fit1, fit2)
print(comp)
} # }
```
