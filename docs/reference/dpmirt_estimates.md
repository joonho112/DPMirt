# Compute posterior estimates using PM, CB, and GR methods

Given a `dpmirt_fit` object, computes person-specific and item-specific
point estimates using multiple posterior summary methods:

- **PM**: Posterior Mean — optimal for individual MSE (Goal 1)

- **CB**: Constrained Bayes (Ghosh, 1992) — optimal for EDF estimation
  (Goal 3)

- **GR**: Triple-Goal (Shen & Louis, 1998) — optimal for ranking + EDF
  (Goals 2 & 3)

## Usage

``` r
dpmirt_estimates(
  fit,
  methods = c("pm", "cb", "gr"),
  alpha = 0.05,
  quantile_type = 7,
  stop_if_ties = FALSE,
  include_items = TRUE
)

# S3 method for class 'dpmirt_estimates'
print(x, ...)
```

## Arguments

- fit:

  A `dpmirt_fit` object from
  [`dpmirt`](https://joonho112.github.io/DPMirt/reference/dpmirt.md).

- methods:

  Character vector of methods to compute. Default `c("pm", "cb", "gr")`.

- alpha:

  Significance level for credible intervals. Default 0.05.

- quantile_type:

  Integer 1-9 for [`quantile()`](https://rdrr.io/r/stats/quantile.html)
  type parameter. Default 7 (R default).

- stop_if_ties:

  Logical. If TRUE, stop when ties detected in posterior mean ranks.
  Default FALSE.

- include_items:

  Logical. If TRUE, apply CB/GR to beta as well. Default TRUE.

- x:

  A `dpmirt_estimates` object.

- ...:

  Additional arguments (currently unused).

## Value

A `dpmirt_estimates` S3 object.

## Details

The three estimators target different inferential goals (Shen & Louis,
1998):

**Goal 1 — Individual estimation**: Minimize individual mean squared
error. The posterior mean (PM) is optimal: \$\$\hat{\theta}^{PM}\_j =
E\[\theta_j \| y\]\$\$

**Goal 2 — Ranking**: Correctly rank individuals. The GR estimator uses
quantiles of the marginal posterior predictive distribution evaluated at
the posterior mean rank of each individual.

**Goal 3 — Distribution estimation**: Recover the empirical distribution
function (EDF). The constrained Bayes (CB) estimator rescales the PM to
match the correct first two moments: \$\$\hat{\theta}^{CB}\_j =
\bar{\theta}^{PM} + \sqrt{1 +
\frac{\bar{\lambda}}{\mathrm{Var}(\theta^{PM})}} \cdot
(\hat{\theta}^{PM}\_j - \bar{\theta}^{PM})\$\$

where \\\bar{\lambda} = \frac{1}{K}\sum_k \lambda_k\\ is the mean
posterior variance.

## References

Ghosh, M. (1992). Constrained Bayes estimation with applications. *JASA,
87*(418), 533–540.

Shen, W., & Louis, T. A. (1998). Triple-goal estimates in two-stage
hierarchical models. *JRSS-B, 60*(2), 455–471.

## See also

[`dpmirt`](https://joonho112.github.io/DPMirt/reference/dpmirt.md),
[`dpmirt_draws`](https://joonho112.github.io/DPMirt/reference/dpmirt_draws.md),
[`dpmirt_loss`](https://joonho112.github.io/DPMirt/reference/dpmirt_loss.md),
[`plot.dpmirt_estimates`](https://joonho112.github.io/DPMirt/reference/plot.dpmirt_estimates.md)

Other estimation:
[`dpmirt_draws()`](https://joonho112.github.io/DPMirt/reference/dpmirt_draws.md),
[`dpmirt_rescale()`](https://joonho112.github.io/DPMirt/reference/dpmirt_rescale.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
fit <- dpmirt(sim$response, model = "rasch", prior = "normal",
              niter = 5000, nburnin = 1000, seed = 123)

# Compute all three estimators
est <- dpmirt_estimates(fit)
print(est)

# Access person estimates
head(est$theta)

# Access item estimates
head(est$beta)

# PM only (faster)
est_pm <- dpmirt_estimates(fit, methods = "pm")
} # }
```
