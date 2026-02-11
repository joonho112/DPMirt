# Elicit a Gamma Hyperprior for the DP Concentration Parameter

Uses the DPprior package to calibrate a Gamma(a, b) hyperprior for the
Dirichlet Process concentration parameter \\\alpha\\, based on the
expected number of clusters \\\mu_K\\ and a confidence level. Falls back
to Paganin's default Gamma(1, 3) if DPprior is not installed.

## Usage

``` r
dpmirt_alpha_prior(N, mu_K = NULL, confidence = "medium", ...)
```

## Arguments

- N:

  Integer. Sample size (number of persons).

- mu_K:

  Numeric or NULL. Expected number of clusters. If NULL, defaults to
  `max(3, ceiling(log(N)))`.

- confidence:

  Character. DPprior confidence level: `"low"`, `"medium"`, or `"high"`.

- ...:

  Additional arguments passed to
  [`DPprior::DPprior_fit`](https://rdrr.io/pkg/DPprior/man/DPprior_fit.html).

## Value

Named numeric vector `c(a = ..., b = ...)` for Gamma(a, b).

## Details

In the CRP representation used by DPMirt, the concentration parameter
\\\alpha\\ controls the expected number of clusters. Larger \\\alpha\\
favors more clusters (closer to nonparametric), while smaller \\\alpha\\
concentrates mass on fewer clusters (closer to parametric).

The DPprior package (Lee, 2026) provides a principled elicitation
framework: given the sample size \\N\\ and an expected cluster count
\\\mu_K\\, it calibrates Gamma(a, b) so that \\E\[K \| \alpha\] \approx
\mu_K\\ with the specified confidence level.

The Paganin et al. (2023) default is Gamma(1, 3), which implies
\\E\[\alpha\] = 1/3\\ — a mildly informative prior favoring few
clusters.

## References

Lee, J. (2026). Design-conditional prior elicitation for Dirichlet
process mixtures. *arXiv:2602.06301*.

Paganin, S., Paciorek, C. J., Wehrhahn, C., Rodriguez, A., Rabe-Hesketh,
S., & de Valpine, P. (2023). Computational strategies and estimation
performance with Bayesian semiparametric item response theory models.
*Journal of Educational and Behavioral Statistics, 48*(2), 147–188.

## See also

[`dpmirt`](https://joonho112.github.io/DPMirt/reference/dpmirt.md),
[`dpmirt_spec`](https://joonho112.github.io/DPMirt/reference/dpmirt_spec.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Default: uses Gamma(1, 3) if DPprior not installed
alpha <- dpmirt_alpha_prior(N = 200)

# Specify expected clusters and confidence
alpha <- dpmirt_alpha_prior(N = 500, mu_K = 5, confidence = "medium")

# Use in model fitting
sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
fit <- dpmirt(sim$response, model = "rasch", prior = "dpm",
              alpha_prior = alpha, niter = 5000, nburnin = 1000)
} # }
```
