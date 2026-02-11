# Run MCMC sampling on a compiled DPMirt model

Executes MCMC sampling using the compiled model. This is the lightweight
step that can be called repeatedly from the same compiled model
(compile-once, sample-many pattern).

## Usage

``` r
dpmirt_sample(
  compiled,
  niter = 10000L,
  nburnin = 2000L,
  thin = 1L,
  thin2 = NULL,
  seed = NULL,
  reset = TRUE,
  verbose = TRUE,
  ...
)

# S3 method for class 'dpmirt_samples'
print(x, ...)
```

## Arguments

- compiled:

  A `dpmirt_compiled` object from
  [`dpmirt_compile`](https://joonho112.github.io/DPMirt/reference/dpmirt_compile.md).

- niter:

  Integer. Total number of MCMC iterations.

- nburnin:

  Integer. Number of burn-in iterations to discard.

- thin:

  Integer. Thinning interval for main monitors.

- thin2:

  Integer or NULL. Thinning interval for monitors2 (eta/theta). If NULL,
  uses the same value as `thin`.

- seed:

  Integer or NULL. Random seed for reproducibility.

- reset:

  Logical. If TRUE (default), reset the MCMC state before sampling. Set
  to FALSE for chain continuation.

- verbose:

  Logical. Print progress messages.

- ...:

  Additional arguments (currently unused).

- x:

  A `dpmirt_samples` object.

## Value

A `dpmirt_samples` S3 object containing:

- samples:

  Matrix of posterior samples from main monitors.

- samples2:

  Matrix of posterior samples from thinned monitors (eta).

- waic:

  WAIC value if computed, otherwise NULL.

- sampling_time:

  Time taken for sampling.

- mcmc_control:

  List of MCMC settings used.

- model_config:

  Reference to model configuration.

- compiled:

  Reference to compiled object (for resume).

## See also

[`dpmirt_compile`](https://joonho112.github.io/DPMirt/reference/dpmirt_compile.md),
[`dpmirt_resume`](https://joonho112.github.io/DPMirt/reference/dpmirt_resume.md),
[`dpmirt_rescale`](https://joonho112.github.io/DPMirt/reference/dpmirt_rescale.md)

Other model fitting:
[`dpmirt()`](https://joonho112.github.io/DPMirt/reference/dpmirt.md),
[`dpmirt_compile()`](https://joonho112.github.io/DPMirt/reference/dpmirt_compile.md),
[`dpmirt_fit-methods`](https://joonho112.github.io/DPMirt/reference/dpmirt_fit-methods.md),
[`dpmirt_resume()`](https://joonho112.github.io/DPMirt/reference/dpmirt_resume.md),
[`dpmirt_spec()`](https://joonho112.github.io/DPMirt/reference/dpmirt_spec.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
spec <- dpmirt_spec(sim$response, model = "rasch", prior = "normal")
compiled <- dpmirt_compile(spec)

# Run MCMC sampling
samples <- dpmirt_sample(compiled, niter = 5000, nburnin = 1000, seed = 123)
print(samples)
} # }
```
