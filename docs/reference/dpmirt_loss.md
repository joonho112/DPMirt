# Evaluate estimator performance using loss functions

Computes loss metrics comparing estimated parameters to true values.
Designed for simulation studies where true parameter values are known.

## Usage

``` r
dpmirt_loss(
  estimates,
  true_theta,
  true_beta = NULL,
  true_lambda = NULL,
  metrics = c("msel", "mselr", "ks"),
  custom_loss = NULL
)
```

## Arguments

- estimates:

  A `dpmirt_estimates` object.

- true_theta:

  Numeric vector of true person abilities.

- true_beta:

  Numeric vector of true item difficulties (optional). For SI
  parameterization, this should be the true gamma (intercept) values.

- true_lambda:

  Numeric vector of true item discriminations (optional). Only relevant
  for 2PL/3PL models.

- metrics:

  Character vector of metrics to compute. Options: `"msel"` (mean
  squared error loss), `"mselr"` (MSE of ranks), `"ks"`
  (Kolmogorov-Smirnov statistic).

- custom_loss:

  Optional custom loss function with signature
  `function(estimate, true)`.

## Value

A data.frame with loss values for each method x parameter combination.

## Details

Three built-in loss metrics measure different aspects of estimation
quality:

**MSEL** (Mean Squared Error Loss): Measures individual-level accuracy.
\$\$MSEL = \frac{1}{K} \sum\_{k=1}^{K} (\hat{\theta}\_k -
\theta_k)^2\$\$

**MSELR** (MSE of Ranks): Measures ranking accuracy using normalized
ranks. \$\$MSELR = \frac{1}{K} \sum\_{k=1}^{K}
\left(\frac{R(\hat{\theta}\_k)}{K} - \frac{R(\theta_k)}{K}\right)^2\$\$

**KS** (Kolmogorov-Smirnov): Measures distributional accuracy. \$\$KS =
\sup_t \|F\_{\hat{\theta}}(t) - F\_{\theta}(t)\|\$\$

Custom loss functions can be supplied via `custom_loss`; they must
accept two vectors (estimates, true values) and return a scalar.

## References

Shen, W., & Louis, T. A. (1998). Triple-goal estimates in two-stage
hierarchical models. *JRSS-B, 60*(2), 455–471.

## See also

[`dpmirt_estimates`](https://joonho112.github.io/DPMirt/reference/dpmirt_estimates.md),
[`dpmirt_simulate`](https://joonho112.github.io/DPMirt/reference/dpmirt_simulate.md)

Other simulation:
[`dpmirt_simulate()`](https://joonho112.github.io/DPMirt/reference/dpmirt_simulate.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
fit <- dpmirt(sim$response, model = "rasch", prior = "normal",
              niter = 5000, nburnin = 1000, seed = 123)
est <- dpmirt_estimates(fit)

# Evaluate against true values
loss <- dpmirt_loss(est, true_theta = sim$theta, true_beta = sim$beta)
print(loss)

# With custom loss function
mae <- function(est, true) mean(abs(est - true))
loss2 <- dpmirt_loss(est, true_theta = sim$theta, custom_loss = mae)
} # }
```
