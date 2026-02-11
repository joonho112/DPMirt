# Quick Start: Your First IRT Model in 5 Minutes

## Overview

This vignette walks you through a complete DPMirt workflow from start to
finish. By the end you will know how to:

1.  **Simulate** IRT response data with
    [`dpmirt_simulate()`](https://joonho112.github.io/DPMirt/reference/dpmirt_simulate.md).
2.  **Fit** a Rasch model under both a Normal prior and a Dirichlet
    Process Mixture (DPM) prior.
3.  **Visualize** posterior densities, item parameters, and MCMC
    diagnostics.
4.  **Extract** person ability estimates using posterior mean (PM),
    constrained Bayes (CB), and triple-goal (GR) estimators.

No NIMBLE compilation is required for the simulation or estimation steps
shown here — we load pre-computed MCMC results so you can follow along
instantly.

## Simulate Data

[`dpmirt_simulate()`](https://joonho112.github.io/DPMirt/reference/dpmirt_simulate.md)
generates binary response data under a known IRT model. When the
optional **IRTsimrel** package is installed, it uses Empirical
Quadrature Calibration (EQC) to hit a target marginal reliability;
otherwise it falls back to a Paganin-style simulation.

``` r

sim <- dpmirt_simulate(
  n_persons    = 200,
  n_items      = 25,
  model        = "rasch",
  target_rho   = 0.8,
  latent_shape = "normal",
  seed         = 42
)
#> Note: Target rho* = 0.800 is near the achievable maximum (0.824) for this configuration.

str(sim, max.level = 1)
#> List of 13
#>  $ response    : num [1:200, 1:25] 0 0 0 1 1 0 1 0 1 0 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ theta       : num [1:200] -2.5215 0.0563 -0.6466 1.135 -0.3096 ...
#>  $ beta        : num [1:25] 0.131 -1.295 0.384 -0.636 0.507 ...
#>  $ lambda      : NULL
#>  $ delta       : NULL
#>  $ n_persons   : num 200
#>  $ n_items     : num 25
#>  $ model       : chr "rasch"
#>  $ reliability : num 0.807
#>  $ target_rho  : num 0.8
#>  $ latent_shape: chr "normal"
#>  $ eqc_result  :List of 16
#>   ..- attr(*, "class")= chr [1:2] "eqc_result" "list"
#>  $ method      : chr "irtsimrel"
#>  - attr(*, "class")= chr "dpmirt_sim"
```

The returned `dpmirt_sim` object contains the binary response matrix
(`sim$response`), the true person abilities (`sim$theta`), and the true
item difficulties (`sim$beta`). You can inspect the simulation with
[`print()`](https://rdrr.io/r/base/print.html):

``` r

sim
#> DPMirt Simulated Data
#> =====================
#> Model:         RASCH 
#> Persons:       200 
#> Items:         25 
#> Distribution:  normal 
#> Method:        irtsimrel 
#> Target rho:    0.8 
#> KR-20:         0.807 
#> EQC c*:        0.9082 
#> EQC rho:       0.8
```

Because
[`dpmirt_simulate()`](https://joonho112.github.io/DPMirt/reference/dpmirt_simulate.md)
is pure R, it runs instantly — no NIMBLE compilation needed.

## Fit a Rasch Model with a Normal Prior

The main entry point is
[`dpmirt()`](https://joonho112.github.io/DPMirt/reference/dpmirt.md).
Below is the call you would run in an interactive session. NIMBLE
compiles the model on first use, which takes roughly 1–2 minutes; after
that, additional sampling is fast.

``` r

# This is the code you would run (takes ~2 minutes for compilation):
fit <- dpmirt(
  sim$response,
  model   = "rasch",
  prior   = "normal",
  niter   = 10000,
  nburnin = 2000,
  seed    = 100
)
```

For this vignette we load a pre-computed result instead:

``` r

fit <- readRDS(find_extdata("vignette_fit_rasch_normal.rds"))
```

## Understanding the Output

A `dpmirt_fit` object stores posterior samples, diagnostics, and model
configuration. The [`print()`](https://rdrr.io/r/base/print.html) method
gives a compact overview:

``` r

print(fit)
#> DPMirt Model Fit
#> ================
#> Model:            RASCH 
#> Prior:            normal 
#> Identification:   constrained_item 
#> Persons (N):      200 
#> Items (I):        25 
#> MCMC:            10000 iterations (2000 burnin, thin=1)
#> WAIC:             6079.03 
#> Total time:       55.4 sec 
#> Min ESS (items):  1714 
#> Min ESS (theta):  1372
```

Key fields:

- **Model / Prior / Identification** — what was specified.
- **MCMC** — iteration count, burn-in, and thinning.
- **WAIC** — Watanabe–Akaike information criterion for model comparison.
- **Min ESS** — minimum effective sample size across items and persons;
  a quick convergence check.
- **Total time** — wall-clock time for the full pipeline.

For a richer view, call
[`summary()`](https://rdrr.io/r/base/summary.html):

``` r

summary(fit)
#> DPMirt Model Summary
#> ====================
#> 
#> Model Configuration:
#>   Model:            RASCH 
#>   Prior:            normal 
#>   Identification:   constrained_item 
#>   Rescaled:         TRUE 
#> 
#> Data:
#>   Persons (N): 200 
#>   Items (I):   25 
#> 
#> MCMC Settings:
#>   Iterations:  10000 
#>   Burn-in:     2000 
#>   Thinning:    1 
#>   Chains:      1 
#> 
#> Timing:
#>   Compilation:  18.1 sec 
#>   Sampling:     36.6 sec 
#>   Total:        55.4 sec 
#> 
#> Item Difficulty (beta) Summary:
#>            Mean    SD
#> beta[1]   0.143 0.151
#> beta[2]  -1.022 0.161
#> beta[3]   0.535 0.158
#> beta[4]  -0.688 0.154
#> beta[5]   0.191 0.152
#> beta[6]   0.501 0.154
#> beta[7]  -0.347 0.149
#> beta[8]  -0.378 0.153
#> beta[9]   0.482 0.156
#> beta[10] -0.843 0.158
#> beta[11]  0.760 0.157
#> beta[12]  0.335 0.155
#> beta[13]  0.147 0.152
#> beta[14] -0.252 0.148
#> beta[15]  0.031 0.148
#> beta[16] -0.038 0.157
#> beta[17] -0.186 0.151
#> beta[18]  0.383 0.151
#> beta[19]  0.999 0.163
#> beta[20]  0.897 0.165
#> beta[21] -0.113 0.153
#> beta[22] -0.692 0.154
#> beta[23]  0.678 0.155
#> beta[24] -0.180 0.151
#> beta[25] -1.342 0.170
#> 
#> Person Ability (theta) Summary:
#>   Range: [ -2.052 ,  1.768 ]
#>   Mean:   -0.073 
#>   SD:     0.81 
#> 
#> Model Comparison:
#>   WAIC:  6079.03
```

The summary adds item-by-item parameter estimates (posterior mean and
SD), a distributional summary of person abilities, and — for DPM models
— cluster and concentration-parameter diagnostics.

## Visualizations

`plot(fit, type = ...)` dispatches to 12 plot types. If **ggplot2** is
installed the package uses it automatically; otherwise base R graphics
are produced.

### Posterior Density of Theta

``` r

plot(fit, type = "density")
```

![Kernel density of the posterior mean theta under a Normal
prior.](quick-start_files/figure-html/plot-density-1.png)

Kernel density of the posterior mean theta under a Normal prior.

This shows the kernel density of the $`N = 200`$ posterior-mean person
abilities. Under a Normal prior the density is smooth and unimodal by
construction.

### Item Difficulty Estimates

``` r

plot(fit, type = "items")
```

![Item difficulty estimates with +/- 2 posterior SD error
bars.](quick-start_files/figure-html/plot-items-1.png)

Item difficulty estimates with +/- 2 posterior SD error bars.

Each point is the posterior mean of $`\beta_j`$; error bars span
$`\pm 2`$ posterior standard deviations. Items are ordered from easiest
(most negative) to hardest.

### MCMC Trace

``` r

plot(fit, type = "trace")
```

![Log-likelihood trace plot. A stationary trace suggests
convergence.](quick-start_files/figure-html/plot-trace-1.png)

Log-likelihood trace plot. A stationary trace suggests convergence.

A stationary log-likelihood trace with no visible trend or drift is a
first-pass convergence check. For formal diagnostics see the **Models
and Workflow** vignette.

## Adding DPM Flexibility

The Normal prior assumes latent abilities follow a Gaussian
distribution. When that assumption is questionable — bimodal
populations, floor/ceiling effects, skew — a Dirichlet Process Mixture
(DPM) prior lets the data speak.

``` r

# This is the code you would run:
fit_dpm <- dpmirt(
  sim$response,
  model   = "rasch",
  prior   = "dpm",
  niter   = 10000,
  nburnin = 2000,
  seed    = 200
)
```

Again, we load a pre-computed result:

``` r

fit_dpm <- readRDS(find_extdata("vignette_fit_rasch_dpm.rds"))
```

### Posterior Density Comparison

``` r

plot(fit_dpm, type = "density")
```

![Posterior mean theta density under the DPM
prior.](quick-start_files/figure-html/plot-dpm-density-1.png)

Posterior mean theta density under the DPM prior.

The overlay below places both posteriors on the same axes for a direct
comparison. Because the true latent distribution here is Normal, the two
densities are nearly identical — the DPM prior adapts to the data and
does not distort the estimates.

``` r

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  pm_normal <- colMeans(fit$theta_samp)
  pm_dpm    <- colMeans(fit_dpm$theta_samp)

  df_overlay <- data.frame(
    value  = c(pm_normal, pm_dpm),
    Prior  = factor(rep(c("Normal", "DPM"), each = length(pm_normal)),
                    levels = c("Normal", "DPM"))
  )

  ggplot(df_overlay, aes(x = value, colour = Prior, fill = Prior)) +
    geom_density(alpha = 0.25, linewidth = 0.9) +
    scale_colour_manual(values = c(Normal = pal$parametric,
                                    DPM   = pal$semiparametric)) +
    scale_fill_manual(values = c(Normal = pal$parametric,
                                  DPM   = pal$semiparametric)) +
    labs(title = "Normal vs. DPM Prior (Normal Population)",
         x = expression(theta), y = "Density") +
    theme_bw() +
    theme(legend.position = "top")
}
#> Warning: package 'ggplot2' was built under R version 4.5.2
```

![Normal-prior vs. DPM-prior posterior mean densities on the same data.
When the true distribution is Normal, both priors recover it equally
well.](quick-start_files/figure-html/density-overlay-1.png)

Normal-prior vs. DPM-prior posterior mean densities on the same data.
When the true distribution is Normal, both priors recover it equally
well.

### Cluster Diagnostics

``` r

plot(fit_dpm, type = "clusters")
```

![Number of active clusters across MCMC iterations (left) and its
posterior distribution
(right).](quick-start_files/figure-html/plot-clusters-1.png)

Number of active clusters across MCMC iterations (left) and its
posterior distribution (right).

The left panel traces the number of active clusters at each post-burn-in
iteration; the right panel shows the posterior distribution of cluster
counts. A single dominant mode suggests the data are well-described by
that many latent subgroups.

### DP Mixture Density

``` r

plot(fit_dpm, type = "dp_density")
```

![DP mixture posterior mean density with 95 percent pointwise credible
band. Dashed line: N(0,1)
reference.](quick-start_files/figure-html/plot-dp-density-1.png)

DP mixture posterior mean density with 95 percent pointwise credible
band. Dashed line: N(0,1) reference.

The solid curve is the posterior mean of the DP mixture density
evaluated on a fine grid; the shaded ribbon is a 95% pointwise credible
band. The dashed line shows the standard Normal for reference.

## Seeing the DPM Advantage

The example above used a truly Normal population, so both priors
performed equally well. But what happens when the population departs
from normality? The DPMirt package ships with pre-computed results for a
**bimodal** population — a 50/50 mixture of two groups centered at
$`\theta = -1.5`$ and $`\theta = 1.5`$ — that reveals the DPM prior’s
key advantage.

``` r

sim_bm   <- readRDS(find_extdata("vignette_sim_bimodal.rds"))
fit_n_bm <- readRDS(find_extdata("vignette_fit_rasch_normal_bimodal.rds"))
fit_d_bm <- readRDS(find_extdata("vignette_fit_rasch_dpm_bimodal.rds"))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  true_theta <- sim_bm$theta
  pm_normal  <- colMeans(fit_n_bm$theta_samp)
  pm_dpm     <- colMeans(fit_d_bm$theta_samp)

  df_bm <- data.frame(
    value  = c(true_theta, pm_normal, pm_dpm),
    source = factor(rep(c("True", "Normal Prior (PM)", "DPM Prior (PM)"),
                        each = length(true_theta)),
                    levels = c("True", "Normal Prior (PM)", "DPM Prior (PM)"))
  )

  ggplot(df_bm, aes(x = value, fill = source, colour = source)) +
    geom_density(alpha = 0.30, linewidth = 0.8) +
    scale_fill_manual(values = c("True"              = pal$reference,
                                  "Normal Prior (PM)" = pal$parametric,
                                  "DPM Prior (PM)"    = pal$semiparametric)) +
    scale_colour_manual(values = c("True"              = pal$reference,
                                    "Normal Prior (PM)" = pal$parametric,
                                    "DPM Prior (PM)"    = pal$semiparametric)) +
    labs(title = "Normal vs. DPM Prior: Recovering a Bimodal Population",
         subtitle = "Rasch model, 25 items, 200 persons",
         x = expression(theta), y = "Density") +
    theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank())
}
```

![When the true population is bimodal, the Normal prior forces a
unimodal fit (blue), masking the two-group structure. The DPM prior
(orange) recovers both
modes.](quick-start_files/figure-html/bimodal-comparison-1.png)

When the true population is bimodal, the Normal prior forces a unimodal
fit (blue), masking the two-group structure. The DPM prior (orange)
recovers both modes.

> **Key insight.** The Normal prior has no mechanism to represent two
> modes, so it forces all estimates toward a single central peak. The
> DPM prior adapts its shape to the data, preserving the bimodal
> structure. This difference is most consequential for classification
> decisions (e.g., identifying students for intervention) and for
> reporting on the shape of the population distribution. See
> [`vignette("posterior-summaries")`](https://joonho112.github.io/DPMirt/articles/posterior-summaries.md)
> for a detailed comparison of estimators that exploit this flexibility.

## Extracting Estimates

[`dpmirt_estimates()`](https://joonho112.github.io/DPMirt/reference/dpmirt_estimates.md)
computes person and item point estimates using three complementary
posterior summary methods:

| Method | Full name                        | Optimizes                   |
|:------:|:---------------------------------|:----------------------------|
| **PM** | Posterior Mean                   | Individual MSE (Goal 1)     |
| **CB** | Constrained Bayes (Ghosh, 1992)  | EDF estimation (Goal 3)     |
| **GR** | Triple-Goal (Shen & Louis, 1998) | Ranking + EDF (Goals 2 & 3) |

``` r

est <- dpmirt_estimates(fit_dpm, methods = c("pm", "cb", "gr"))
```

The `theta` element is a data frame with one row per person:

``` r

head(est$theta, 10)
#>            theta_pm theta_psd    theta_cb     theta_gr      rbar rhat
#> eta[1]  -1.64208358 0.4850620 -1.82658631 -1.898414370  14.10438    6
#> eta[2]   0.08246350 0.3785971  0.10105187  0.130777459 110.04200  114
#> eta[3]   0.05810186 0.3859230  0.07382128  0.001450727 108.21700  103
#> eta[4]   1.15463248 0.4241983  1.29948472  1.393735077 179.26688  192
#> eta[5]  -0.68190514 0.3877922 -0.75333256 -0.665180404  53.89862   52
#> eta[6]   0.81760681 0.3950633  0.92276925  0.811079297 162.99625  167
#> eta[7]   1.72015915 0.4758415  1.93161057  2.206872410 194.11588  200
#> eta[8]  -0.06517854 0.3837235 -0.06397722 -0.057278283  98.39575   98
#> eta[9]   0.66101983 0.3824489  0.74774183  0.688329335 153.44163  159
#> eta[10] -0.06306564 0.3877453 -0.06161549 -0.033409063  98.65662  100
#>         theta_lower theta_upper
#> eta[1]  -2.71809760 -0.76137677
#> eta[2]  -0.63931121  0.81955388
#> eta[3]  -0.68900612  0.82161938
#> eta[4]   0.34563696  2.00481098
#> eta[5]  -1.44006283  0.07806641
#> eta[6]   0.05942028  1.60027832
#> eta[7]   0.84623291  2.75717002
#> eta[8]  -0.84149547  0.68498198
#> eta[9]  -0.06723822  1.41774065
#> eta[10] -0.82119382  0.69405528
```

> **Which estimator should you use?** It depends on your inferential
> goal. If you need the best point prediction for each individual, use
> **PM**. If you need the set of estimates to reproduce the shape of the
> ability distribution (e.g., for group-level reporting), use **CB**. If
> you need both accurate rankings *and* distributional fidelity, use
> **GR**. See the **Posterior Summaries** vignette for an in-depth
> comparison.

Item estimates are available in `est$beta`:

``` r

head(est$beta, 10)
#>             beta_pm  beta_psd    beta_cb     beta_gr     rbar rhat beta_lower
#> beta[1]   0.1393281 0.1481244  0.1437448  0.09951606 14.54713   14 -0.1448619
#> beta[2]  -1.0281008 0.1593125 -1.0606917 -1.05388296  2.26200    2 -1.3531175
#> beta[3]   0.5266330 0.1558563  0.5433273  0.59943417 19.94138   21  0.2277232
#> beta[4]  -0.6922001 0.1558382 -0.7141429 -0.73611032  4.36675    4 -0.9959212
#> beta[5]   0.1946660 0.1501926  0.2008369  0.23378009 15.36788   16 -0.1017568
#> beta[6]   0.5082932 0.1550549  0.5244062  0.52121043 19.70375   20  0.2010165
#> beta[7]  -0.3503717 0.1511989 -0.3614786 -0.36061570  7.43475    7 -0.6562346
#> beta[8]  -0.3770295 0.1534740 -0.3889813 -0.46590990  7.14050    6 -0.6858021
#> beta[9]   0.4868289 0.1543031  0.5022614  0.44842581 19.43700   19  0.1870308
#> beta[10] -0.8442989 0.1608205 -0.8710633 -0.87137757  3.31225    3 -1.1662977
#>           beta_upper
#> beta[1]   0.44597696
#> beta[2]  -0.72312083
#> beta[3]   0.84615054
#> beta[4]  -0.39020070
#> beta[5]   0.48423148
#> beta[6]   0.80991762
#> beta[7]  -0.06704874
#> beta[8]  -0.08586970
#> beta[9]   0.78795591
#> beta[10] -0.52819222
```

## What’s Next?

You now have a working end-to-end pipeline. Depending on your goal, the
following vignettes go deeper:

| Step | What to read | Why |
|:---|:---|:---|
| Understand models | *Models and Workflow* | Step-by-step control over specification, compilation, and sampling; 2PL and 3PL models |
| Choose an estimator | *Posterior Summaries* | When PM vs CB vs GR; shrinkage diagnostics and loss evaluation |
| Set DPM priors | *Prior Elicitation* | Principled choice of the concentration parameter $`\alpha`$ via **DPprior** |

------------------------------------------------------------------------

**Happy modeling!** You are now equipped to fit your own Rasch models
with both parametric (Normal) and semiparametric (DPM) priors. For
questions, bugs, or suggestions, visit the [DPMirt
repository](https://github.com/joonho112/DPMirt).
