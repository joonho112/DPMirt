# Rescale posterior samples for identification

Applies post-hoc rescaling to MCMC posterior samples to resolve
identification indeterminacy. For unconstrained models, this centers
item difficulties and (for 2PL/3PL) normalizes discriminations.

## Usage

``` r
dpmirt_rescale(samples_obj, rescale = TRUE)
```

## Arguments

- samples_obj:

  A `dpmirt_samples` object.

- rescale:

  Logical. If TRUE (default), apply rescaling. If FALSE, return samples
  as-is.

## Value

A `dpmirt_fit` S3 object containing rescaled posterior samples,
diagnostics (ESS, WAIC, cluster info), and model configuration. This
object can be passed directly to
[`dpmirt_estimates`](https://joonho112.github.io/DPMirt/reference/dpmirt_estimates.md),
[`dpmirt_resume`](https://joonho112.github.io/DPMirt/reference/dpmirt_resume.md),
and other downstream functions.

## Details

Post-hoc rescaling resolves the identification indeterminacy inherent in
unconstrained IRT models:

**Rasch** (location only): For each MCMC iteration \\s\\:
\$\$\beta^\*\_i = \beta_i - \bar{\beta}\$\$ \$\$\theta^\*\_j =
\theta_j - \bar{\beta}\$\$

**2PL/3PL IRT parameterization** (location + scale): Let \\c_s =
\bar{\beta}\_s\\ and \\d_s = (\prod_i \lambda_i)^{-1/I}\\:
\$\$\beta^\*\_i = (\beta_i - c_s) / d_s\$\$ \$\$\lambda^\*\_i =
\lambda_i \cdot d_s\$\$ \$\$\theta^\*\_j = (\theta_j - c_s) / d_s\$\$

After rescaling: \\\bar{\beta}^\* = 0\\ and \\(\prod_i
\lambda^\*\_i)^{1/I} = 1\\.

## References

Paganin, S., Paciorek, C. J., Wehrhahn, C., Rodriguez, A., Rabe-Hesketh,
S., & de Valpine, P. (2023). Computational strategies and estimation
performance with Bayesian semiparametric item response theory models.
*Journal of Educational and Behavioral Statistics, 48*(2), 147–188.

## See also

[`dpmirt_sample`](https://joonho112.github.io/DPMirt/reference/dpmirt_sample.md),
[`dpmirt_estimates`](https://joonho112.github.io/DPMirt/reference/dpmirt_estimates.md)

Other estimation:
[`dpmirt_draws()`](https://joonho112.github.io/DPMirt/reference/dpmirt_draws.md),
[`dpmirt_estimates()`](https://joonho112.github.io/DPMirt/reference/dpmirt_estimates.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
spec <- dpmirt_spec(sim$response, model = "rasch", prior = "normal",
                    identification = "unconstrained")
compiled <- dpmirt_compile(spec)
samples <- dpmirt_sample(compiled, niter = 5000, nburnin = 1000)

# Apply rescaling
rescaled <- dpmirt_rescale(samples)
} # }
```
