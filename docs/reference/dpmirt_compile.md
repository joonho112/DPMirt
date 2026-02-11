# Compile a DPMirt model specification

Takes a `dpmirt_spec` object and compiles the NIMBLE model and MCMC
engine. This is the expensive step (~30-120 seconds) that only needs to
be done once per model specification.

## Usage

``` r
dpmirt_compile(
  spec,
  sampler_config = NULL,
  use_centered_sampler = "auto",
  enable_waic = TRUE,
  enable_logprob_monitor = TRUE,
  verbose = TRUE,
  ...
)

# S3 method for class 'dpmirt_compiled'
print(x, ...)
```

## Arguments

- spec:

  A `dpmirt_spec` object from
  [`dpmirt_spec`](https://joonho112.github.io/DPMirt/reference/dpmirt_spec.md).

- sampler_config:

  Optional custom MCMC sampler configuration.

- use_centered_sampler:

  Character or logical. Whether to use the centered sampler for SI
  parameterization. "auto" enables it when appropriate (SI param +
  2PL/3PL). Only relevant for Phase 3+.

- enable_waic:

  Logical. Whether to enable WAIC computation.

- enable_logprob_monitor:

  Logical. Whether to add log-probability monitoring samplers.

- verbose:

  Logical. Print progress messages.

- ...:

  Additional arguments (currently unused).

- x:

  A `dpmirt_compiled` object.

## Value

A `dpmirt_compiled` S3 object containing the compiled model, compiled
MCMC, specification reference, and session signature.

## Details

The compiled NIMBLE objects contain external C++ pointers that **cannot
be serialized** across R sessions. The compile-once pattern is therefore
a within-session optimization. For cross-session workflows, save the
`dpmirt_spec` object and recompile in the new session.

## See also

[`dpmirt_spec`](https://joonho112.github.io/DPMirt/reference/dpmirt_spec.md),
[`dpmirt_sample`](https://joonho112.github.io/DPMirt/reference/dpmirt_sample.md),
[`dpmirt`](https://joonho112.github.io/DPMirt/reference/dpmirt.md)

Other model fitting:
[`dpmirt()`](https://joonho112.github.io/DPMirt/reference/dpmirt.md),
[`dpmirt_fit-methods`](https://joonho112.github.io/DPMirt/reference/dpmirt_fit-methods.md),
[`dpmirt_resume()`](https://joonho112.github.io/DPMirt/reference/dpmirt_resume.md),
[`dpmirt_sample()`](https://joonho112.github.io/DPMirt/reference/dpmirt_sample.md),
[`dpmirt_spec()`](https://joonho112.github.io/DPMirt/reference/dpmirt_spec.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
spec <- dpmirt_spec(sim$response, model = "rasch", prior = "normal")

# Compile (takes 30-120 seconds)
compiled <- dpmirt_compile(spec)
print(compiled)

# Compile-once, sample-many pattern
samples1 <- dpmirt_sample(compiled, niter = 5000, nburnin = 1000, seed = 1)
samples2 <- dpmirt_sample(compiled, niter = 5000, nburnin = 1000, seed = 2)
} # }
```
