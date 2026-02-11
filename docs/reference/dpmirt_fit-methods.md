# Methods for dpmirt_fit Objects

Print, summarize, and extract coefficients from a fitted DPMirt model.

## Usage

``` r
# S3 method for class 'dpmirt_fit'
print(x, ...)

# S3 method for class 'dpmirt_fit'
summary(object, ...)

# S3 method for class 'dpmirt_fit'
coef(object, type = c("items", "persons"), ...)
```

## Arguments

- x, object:

  A `dpmirt_fit` object from
  [`dpmirt`](https://joonho112.github.io/DPMirt/reference/dpmirt.md).

- ...:

  Additional arguments (currently unused).

- type:

  For `coef`: character, `"items"` (default) or `"persons"`.

## Value

`print` and `summary` return the input invisibly. `coef` returns a
data.frame of posterior mean estimates with one row per item or person.

## Details

`print` displays a concise one-screen summary: model type, data
dimensions, MCMC settings, WAIC, minimum ESS, and DPM cluster summary.

`summary` displays a comprehensive report including item parameter
posterior means and SDs, person ability range, DPM diagnostics, and
timing.

`coef` extracts posterior mean point estimates for items
(`type = "items"`) or persons (`type = "persons"`).

## See also

[`dpmirt`](https://joonho112.github.io/DPMirt/reference/dpmirt.md) for
model fitting,
[`dpmirt_estimates`](https://joonho112.github.io/DPMirt/reference/dpmirt_estimates.md)
for PM/CB/GR estimators,
[`dpmirt_diagnostics`](https://joonho112.github.io/DPMirt/reference/dpmirt_diagnostics.md)
for detailed MCMC diagnostics

Other model fitting:
[`dpmirt()`](https://joonho112.github.io/DPMirt/reference/dpmirt.md),
[`dpmirt_compile()`](https://joonho112.github.io/DPMirt/reference/dpmirt_compile.md),
[`dpmirt_resume()`](https://joonho112.github.io/DPMirt/reference/dpmirt_resume.md),
[`dpmirt_sample()`](https://joonho112.github.io/DPMirt/reference/dpmirt_sample.md),
[`dpmirt_spec()`](https://joonho112.github.io/DPMirt/reference/dpmirt_spec.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
fit <- dpmirt(sim$response, model = "rasch", prior = "normal",
              niter = 5000, nburnin = 1000, seed = 123)

# Print concise summary
print(fit)

# Detailed summary with item parameters
summary(fit)

# Extract item difficulty estimates
coef(fit, type = "items")

# Extract person ability estimates
coef(fit, type = "persons")
} # }
```
