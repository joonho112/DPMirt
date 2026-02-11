# The Complete Guide to Models and Workflows

## 1. Overview

DPMirt provides two ways to fit Bayesian IRT models.

**One-step workflow.** Call
[`dpmirt()`](https://joonho112.github.io/DPMirt/reference/dpmirt.md) and
get a fully processed `dpmirt_fit` object in one function call. Best for
standard analyses where the defaults are appropriate.

**Step-by-step workflow.** Walk through the pipeline explicitly:
[`dpmirt_spec()`](https://joonho112.github.io/DPMirt/reference/dpmirt_spec.md)
$`\rightarrow`$[`dpmirt_compile()`](https://joonho112.github.io/DPMirt/reference/dpmirt_compile.md)
$`\rightarrow`$[`dpmirt_sample()`](https://joonho112.github.io/DPMirt/reference/dpmirt_sample.md)
$`\rightarrow`$[`dpmirt_rescale()`](https://joonho112.github.io/DPMirt/reference/dpmirt_rescale.md)
$`\rightarrow`$[`dpmirt_estimates()`](https://joonho112.github.io/DPMirt/reference/dpmirt_estimates.md).
Best when you need to inspect intermediate objects, run multiple chains
from a single compiled model (compile-once, sample-many), or modify the
MCMC configuration before sampling.

Both workflows produce the same statistical results.

## 2. The One-Step Workflow

### 2.1 Key Arguments

| Argument | Default | Description |
|:---|:---|:---|
| data | — | N x I binary response matrix (persons in rows, items in columns) |
| model | “rasch” | IRT model: “rasch”, “2pl”, or “3pl” |
| prior | “normal” | Latent trait prior: “normal” or “dpm” |
| parameterization | “irt” | “irt” or “si” (slope-intercept; 2PL/3PL only) |
| identification | NULL | Identification strategy; NULL selects model default |
| niter | 10000 | Total MCMC iterations |
| nburnin | 2000 | Burn-in iterations to discard |
| thin | 1 | Thinning interval for item parameters |
| thin2 | NULL | Thinning for theta (NULL = same as thin) |
| nchains | 1 | Number of independent chains |
| seed | NULL | Random seed for reproducibility |
| alpha_prior | NULL | DPM concentration hyperprior: NULL, c(a,b), or DPprior object |
| base_measure | list(s2_mu=2, …) | DPM base measure hyperparameters (Paganin defaults) |
| M | 50 | CRP truncation level (max clusters) |
| rescale | TRUE | Apply post-hoc identification rescaling |
| compute_waic | TRUE | Compute WAIC for model comparison |
| save_draws | TRUE | Store full theta posterior (needed for CB/GR estimators) |

### 2.2 Examples

Each call below takes roughly 1–3 minutes (dominated by NIMBLE
compilation on first use). We show the code and load pre-computed
results.

**Rasch – Normal prior**

``` r

fit_rasch_n <- dpmirt(
  sim$response, model = "rasch", prior = "normal",
  niter = 10000, nburnin = 2000, seed = 100
)
```

``` r

fit_rasch_n <- readRDS(find_extdata("vignette_fit_rasch_normal.rds"))
```

**Rasch – DPM prior**

``` r

fit_rasch_dpm <- dpmirt(
  sim$response, model = "rasch", prior = "dpm",
  niter = 10000, nburnin = 2000, seed = 200
)
```

``` r

fit_rasch_dpm <- readRDS(find_extdata("vignette_fit_rasch_dpm.rds"))
```

**2PL – DPM prior**

``` r

fit_2pl <- dpmirt(
  sim$response, model = "2pl", prior = "dpm",
  niter = 15000, nburnin = 5000, seed = 300
)
```

``` r

fit_2pl <- readRDS(find_extdata("vignette_fit_2pl_dpm.rds"))
```

The returned `dpmirt_fit` object gives a compact overview when printed,
and [`summary()`](https://rdrr.io/r/base/summary.html) adds item-level
estimates and DPM cluster diagnostics:

``` r

fit_rasch_n
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

``` r

summary(fit_rasch_dpm)
#> DPMirt Model Summary
#> ====================
#> 
#> Model Configuration:
#>   Model:            RASCH 
#>   Prior:            dpm 
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
#>   Compilation:  21.3 sec 
#>   Sampling:     37.0 sec 
#>   Total:        1.5 min 
#> 
#> Item Difficulty (beta) Summary:
#>            Mean    SD
#> beta[1]   0.139 0.148
#> beta[2]  -1.028 0.159
#> beta[3]   0.527 0.156
#> beta[4]  -0.692 0.156
#> beta[5]   0.195 0.150
#> beta[6]   0.508 0.155
#> beta[7]  -0.350 0.151
#> beta[8]  -0.377 0.153
#> beta[9]   0.487 0.154
#> beta[10] -0.844 0.161
#> beta[11]  0.767 0.159
#> beta[12]  0.333 0.150
#> beta[13]  0.140 0.152
#> beta[14] -0.251 0.152
#> beta[15]  0.029 0.145
#> beta[16] -0.039 0.147
#> beta[17] -0.178 0.152
#> beta[18]  0.382 0.155
#> beta[19]  1.003 0.166
#> beta[20]  0.897 0.164
#> beta[21] -0.116 0.150
#> beta[22] -0.688 0.153
#> beta[23]  0.688 0.156
#> beta[24] -0.182 0.151
#> beta[25] -1.350 0.173
#> 
#> Person Ability (theta) Summary:
#>   Range: [ -2.127 ,  1.72 ]
#>   Mean:   -0.075 
#>   SD:     0.814 
#> 
#> Model Comparison:
#>   WAIC:  6078.17 
#> 
#> DPM Diagnostics:
#>   Alpha (concentration):
#>     Posterior mean:  0.271 
#>     95% CI: [ 0.007 ,  0.982 ]
#>   Number of clusters:
#>     Posterior mean:  2.3 
#>     Mode:           1 
#>     Range: [ 1 ,  13 ]
#>   DP density:    computed (500 grid points)
```

## 3. The Step-by-Step Workflow

### 3.1 Specification

[`dpmirt_spec()`](https://joonho112.github.io/DPMirt/reference/dpmirt_spec.md)
translates your model request into a complete NIMBLE specification:
programmatically generated model code, constants, data, initial values,
and monitor lists.

``` r

spec <- dpmirt_spec(
  data = sim$response, model = "rasch",
  prior = "dpm", identification = "constrained_item"
)
```

Inspect the returned `dpmirt_spec` object:

``` r

spec$code       # nimbleCode object
str(spec$constants)
spec$monitors   # beta, alpha, zi, muTilde, s2Tilde, logprob nodes
spec$monitors2  # eta (thinned separately via thin2)
```

Because the spec is pure R data (no C++ pointers), it can be safely
saved to disk with [`saveRDS()`](https://rdrr.io/r/base/readRDS.html)
and loaded in a new session.

### 3.2 Compilation

[`dpmirt_compile()`](https://joonho112.github.io/DPMirt/reference/dpmirt_compile.md)
feeds the spec into NIMBLE’s compiler, producing C++ objects that run
the sampler at native speed.

``` r

compiled <- dpmirt_compile(spec)
```

> **Timing warning.** Compilation typically takes 30–120 seconds. It
> only needs to be done **once** per model specification.

> **Session-bound caveat.** The compiled object contains external C++
> pointers and **cannot** be serialized with
> [`saveRDS()`](https://rdrr.io/r/base/readRDS.html). If you restart R,
> re-run
> [`dpmirt_compile()`](https://joonho112.github.io/DPMirt/reference/dpmirt_compile.md)
> on the saved spec.

### 3.3 Sampling

[`dpmirt_sample()`](https://joonho112.github.io/DPMirt/reference/dpmirt_sample.md)
runs MCMC on the compiled model. Because compilation is already done,
this step is fast and repeatable.

``` r

samples <- dpmirt_sample(
  compiled, niter = 10000, nburnin = 2000,
  thin = 1, thin2 = 1, seed = 42
)
```

**Compile-once, sample-many.** Draw multiple independent chains from the
same compiled model without paying the compilation cost each time:

``` r

chains <- lapply(1:4, function(i) {
  dpmirt_sample(compiled, niter = 10000, nburnin = 2000, seed = i)
})
```

For convenience,
[`dpmirt()`](https://joonho112.github.io/DPMirt/reference/dpmirt.md)
also supports multi-chain runs directly via the `nchains` argument.

### 3.4 Rescaling

[`dpmirt_rescale()`](https://joonho112.github.io/DPMirt/reference/dpmirt_rescale.md)
applies post-hoc identification rescaling to the raw posterior samples
and returns a `dpmirt_fit` object that can be passed directly to
[`dpmirt_estimates()`](https://joonho112.github.io/DPMirt/reference/dpmirt_estimates.md),
[`dpmirt_resume()`](https://joonho112.github.io/DPMirt/reference/dpmirt_resume.md),
and other downstream functions.

| Model | Indeterminacy | Rescaling |
|:---|:---|:---|
| Rasch | Location only | Center beta at 0; shift theta accordingly |
| 2PL/3PL (IRT) | Location + scale | Center beta, normalize lambda to geom-mean 1, rescale theta |
| 2PL/3PL (SI) | Location + scale | Center gamma via weighted sum, normalize lambda, rescale theta |

For the Rasch model:

``` math
\beta_i^* = \beta_i - \bar{\beta}, \qquad
\theta_j^* = \theta_j - \bar{\beta}
```

For 2PL/3PL (IRT parameterization):

``` math
s = \left(\prod_{i=1}^{I} \lambda_i\right)^{-1/I}, \quad
\beta_i^* = \frac{\beta_i - \bar{\beta}}{s}, \quad
\lambda_i^* = \lambda_i \cdot s, \quad
\theta_j^* = \frac{\theta_j - \bar{\beta}}{s}
```

``` r

fit <- dpmirt_rescale(samples)

# fit is a dpmirt_fit: pass directly to estimates or resume
est <- dpmirt_estimates(fit)
```

> If the model uses `"constrained_item"` or `"constrained_ability"`
> identification, rescaling is a no-op — parameters are already
> identified through in-model constraints.

### 3.5 Estimates

[`dpmirt_estimates()`](https://joonho112.github.io/DPMirt/reference/dpmirt_estimates.md)
computes person and item point estimates using three posterior summary
methods.

| Method | Name                             | Optimizes                   |
|:-------|:---------------------------------|:----------------------------|
| PM     | Posterior Mean                   | Individual MSE (Goal 1)     |
| CB     | Constrained Bayes (Ghosh, 1992)  | EDF estimation (Goal 3)     |
| GR     | Triple-Goal (Shen & Louis, 1998) | Ranking + EDF (Goals 2 & 3) |

``` r

est <- dpmirt_estimates(fit_rasch_dpm, methods = c("pm", "cb", "gr"))
head(est$theta, 8)
#>           theta_pm theta_psd    theta_cb     theta_gr      rbar rhat
#> eta[1] -1.64208358 0.4850620 -1.82658631 -1.898414370  14.10438    6
#> eta[2]  0.08246350 0.3785971  0.10105187  0.130777459 110.04200  114
#> eta[3]  0.05810186 0.3859230  0.07382128  0.001450727 108.21700  103
#> eta[4]  1.15463248 0.4241983  1.29948472  1.393735077 179.26688  192
#> eta[5] -0.68190514 0.3877922 -0.75333256 -0.665180404  53.89862   52
#> eta[6]  0.81760681 0.3950633  0.92276925  0.811079297 162.99625  167
#> eta[7]  1.72015915 0.4758415  1.93161057  2.206872410 194.11588  200
#> eta[8] -0.06517854 0.3837235 -0.06397722 -0.057278283  98.39575   98
#>        theta_lower theta_upper
#> eta[1] -2.71809760 -0.76137677
#> eta[2] -0.63931121  0.81955388
#> eta[3] -0.68900612  0.82161938
#> eta[4]  0.34563696  2.00481098
#> eta[5] -1.44006283  0.07806641
#> eta[6]  0.05942028  1.60027832
#> eta[7]  0.84623291  2.75717002
#> eta[8] -0.84149547  0.68498198
```

``` r

head(est$beta, 8)
#>            beta_pm  beta_psd    beta_cb     beta_gr     rbar rhat beta_lower
#> beta[1]  0.1393281 0.1481244  0.1437448  0.09951606 14.54713   14 -0.1448619
#> beta[2] -1.0281008 0.1593125 -1.0606917 -1.05388296  2.26200    2 -1.3531175
#> beta[3]  0.5266330 0.1558563  0.5433273  0.59943417 19.94138   21  0.2277232
#> beta[4] -0.6922001 0.1558382 -0.7141429 -0.73611032  4.36675    4 -0.9959212
#> beta[5]  0.1946660 0.1501926  0.2008369  0.23378009 15.36788   16 -0.1017568
#> beta[6]  0.5082932 0.1550549  0.5244062  0.52121043 19.70375   20  0.2010165
#> beta[7] -0.3503717 0.1511989 -0.3614786 -0.36061570  7.43475    7 -0.6562346
#> beta[8] -0.3770295 0.1534740 -0.3889813 -0.46590990  7.14050    6 -0.6858021
#>          beta_upper
#> beta[1]  0.44597696
#> beta[2] -0.72312083
#> beta[3]  0.84615054
#> beta[4] -0.39020070
#> beta[5]  0.48423148
#> beta[6]  0.80991762
#> beta[7] -0.06704874
#> beta[8] -0.08586970
```

## 4. Model Specifications

All models use a logistic link for
$`\pi_{ij} = P(Y_{ij} = 1 \mid \theta_j, \text{item}_i)`$.

### 4.1 Rasch Model

``` math
\text{logit}(\pi_{ij}) = \theta_j - \beta_i
```

The Rasch model has a single item parameter — difficulty $`\beta_i`$.
The response probability depends only on the difference between person
ability and item difficulty.

**Identification.** Location indeterminacy only (adding a constant $`c`$
to all $`\theta`$ and $`\beta`$ leaves the likelihood unchanged).
Default: `"constrained_item"` (mean-centers $`\beta`$ during MCMC).

**When to use.** Most parsimonious IRT model. Assumes equal
discrimination and no guessing. Stable item estimates even with modest
sample sizes.

``` r

fit_rasch <- dpmirt(data, model = "rasch", prior = "normal")
```

### 4.2 Two-Parameter Logistic (2PL)

Adds a discrimination parameter $`\lambda_i`$. Two equivalent
parameterizations are available.

**IRT parameterization:**
``` math
\text{logit}(\pi_{ij}) = \lambda_i(\theta_j - \beta_i)
```

**Slope–intercept (SI) parameterization:**
``` math
\text{logit}(\pi_{ij}) = \lambda_i \theta_j + \gamma_i
```

where $`\gamma_i = -\lambda_i \beta_i`$ is an intercept.

| Feature | IRT | SI |
|:---|:---|:---|
| Formula | $`\lambda_i(\theta_j - \beta_i)`$ | $`\lambda_i \theta_j + \gamma_i`$ |
| Item location | $`\beta_i`$ (difficulty) | $`\gamma_i`$ (intercept) |
| Interpretation | $`\beta_i`$ on theta scale | Baseline log-odds at $`\theta = 0`$ |
| Posterior correlation | Higher ($`\lambda`$–$`\beta`$) | Lower ($`\lambda`$–$`\gamma`$) |
| Centered sampler | Not applicable | Auto-enabled |
| Default identification | unconstrained | unconstrained |

The SI parameterization reduces posterior correlations between
discrimination and intercept, improving MCMC mixing. DPMirt
automatically enables a centered sampler for SI (Paganin et al., 2023).

**Identification.** Both location and scale indeterminacy. Default:
`"unconstrained"` + post-hoc rescaling.

``` r

fit_2pl_irt <- dpmirt(data, model = "2pl", prior = "dpm")
fit_2pl_si  <- dpmirt(data, model = "2pl", prior = "dpm",
                       parameterization = "si")
```

### 4.3 Three-Parameter Logistic (3PL)

Adds a lower-asymptote (guessing) parameter $`\delta_i \in (0, 1)`$:

``` math
\pi_{ij} = \delta_i + (1 - \delta_i)\,\text{logistic}\!\big(\lambda_i(\theta_j - \beta_i)\big)
```

DPMirt places a $`\text{Beta}(4, 12)`$ prior on each $`\delta_i`$ (prior
mean $`\approx 0.25`$, concentrating mass in 0.1–0.4).

**Practical recommendations:**

- Use longer chains (20,000+ iterations, 5,000+ burn-in).
- Most appropriate for multiple-choice items with genuine guessing.
- Both IRT and SI parameterizations are supported.

``` r

fit_3pl <- dpmirt(
  data, model = "3pl", prior = "dpm",
  niter = 20000, nburnin = 5000, seed = 400
)
```

## 5. Priors: Parametric vs Semiparametric

### 5.1 Normal (Parametric) Prior

``` math
\theta_j \sim N(\mu, \sigma^2), \qquad
\mu \sim N(0, 3), \qquad
\sigma^2 \sim \text{Inv-Gamma}(2.01, 1.01)
```

Standard assumption in most IRT software. Works well when the ability
distribution is approximately unimodal and symmetric.

### 5.2 DPM (Semiparametric) Prior

The DPM prior uses the Chinese Restaurant Process (CRP) representation:

``` math
z_j \mid \alpha \sim \text{CRP}(\alpha), \qquad
\theta_j \mid z_j \sim N\!\big(\tilde{\mu}_{z_j},\; \tilde{\sigma}^2_{z_j}\big)
```

Each cluster $`m`$ draws from the base measure $`G_0`$:

``` math
\tilde{\mu}_m \sim N(0, \sigma^2_\mu), \qquad
\tilde{\sigma}^2_m \sim \text{Inv-Gamma}(\nu_1, \nu_2)
```

The concentration parameter $`\alpha`$ controls the expected number of
clusters.

### 5.3 Hyperparameter Defaults

| Parameter | Default | Description |
|:---|:---|:---|
| $`\sigma^2_\mu`$ | 2 | Base measure: cluster mean prior variance |
| $`\nu_1`$ | 2.01 | Base measure: Inv-Gamma shape for cluster variance |
| $`\nu_2`$ | 1.01 | Base measure: Inv-Gamma rate for cluster variance |
| $`M`$ | 50 | CRP truncation level (maximum clusters) |
| $`\alpha`$ | Gamma(1, 3) | Concentration parameter prior (Paganin default) |

These can be modified via `base_measure`:

``` r

fit <- dpmirt(data, model = "rasch", prior = "dpm",
              base_measure = list(s2_mu = 4, nu1 = 2.5, nu2 = 1.5),
              M = 100)
```

### 5.4 When to Choose the DPM Prior

- **Non-normal populations.** Heterogeneous test-taker groups (e.g.,
  mixing advanced and introductory students).
- **Multimodal distributions.** Floor/ceiling effects or distinct
  subgroups.
- **Exploratory analysis.** When you have little prior information about
  the ability distribution shape.
- **Shrinkage calibration.** The DPM prior produces less severe
  shrinkage when the true distribution is non-Normal.

When the true distribution is Normal, the DPM converges to a single
cluster and reproduces the Normal-prior results with only a modest
computation-time increase.

### 5.5 Normal vs. DPM: When Normality Fails

The guidelines above are easier to appreciate with a concrete example.
The pre-computed bimodal simulation ships with the package: 200 persons
drawn from a two-group mixture ($`\theta \sim 0.5\,N(-1.5, 0.5) +
0.5\,N(1.5, 0.5)`$) assessed on 25 Rasch items
($`\bar{w} \approx 0.8`$).

``` r

sim_bm    <- readRDS(find_extdata("vignette_sim_bimodal.rds"))
fit_bm_n  <- readRDS(find_extdata("vignette_fit_rasch_normal_bimodal.rds"))
fit_bm_d  <- readRDS(find_extdata("vignette_fit_rasch_dpm_bimodal.rds"))
```

WAIC provides a model-level comparison:

``` r

dpmirt_compare(fit_bm_n, fit_bm_d)
#>          model     waic delta_waic
#> 1 RASCH-normal 4009.778   0.000000
#> 2    RASCH-dpm 4010.544   0.766108
```

In this example the two WAIC values are nearly identical
($`|\Delta\text{WAIC}| < 2`$), which is expected: WAIC measures
predictive accuracy for the *binary responses*, and both priors predict
item responses equally well. The DPM’s real advantage shows in the
*recovered ability distribution* — a difference that WAIC does not
capture but the density plot below makes visually obvious.

``` r

true_theta <- sim_bm$theta
pm_normal  <- colMeans(fit_bm_n$theta_samp)
pm_dpm     <- colMeans(fit_bm_d$theta_samp)

df_bm <- data.frame(
  value  = c(true_theta, pm_normal, pm_dpm),
  source = factor(rep(c("True", "Normal (PM)", "DPM (PM)"),
                      each = length(true_theta)),
                  levels = c("True", "Normal (PM)", "DPM (PM)"))
)

ggplot(df_bm, aes(x = value, fill = source, colour = source)) +
  geom_density(alpha = 0.25, linewidth = 0.8) +
  scale_fill_manual(values = c("True"         = pal$reference,
                                "Normal (PM)" = pal$parametric,
                                "DPM (PM)"    = pal$semiparametric)) +
  scale_colour_manual(values = c("True"         = pal$reference,
                                  "Normal (PM)" = pal$parametric,
                                  "DPM (PM)"    = pal$semiparametric)) +
  labs(title = "Recovering a Bimodal Population",
       x = expression(theta), y = "Density") +
  theme_bw() +
  theme(legend.position = "top", legend.title = element_blank())
```

![True bimodal density (gray) versus posterior mean estimates under
Normal (blue) and DPM (orange) priors. The Normal prior forces a
unimodal fit; the DPM prior recovers both
modes.](models-and-workflow_files/figure-html/bimodal-density-1.png)

True bimodal density (gray) versus posterior mean estimates under Normal
(blue) and DPM (orange) priors. The Normal prior forces a unimodal fit;
the DPM prior recovers both modes.

> **Take-away.** When the true distribution departs from normality, the
> DPM prior produces posterior mean estimates that track the true
> distributional shape far more accurately than the Normal prior. WAIC
> may not detect this advantage because it evaluates item-response
> prediction, not ability-distribution recovery. For practical guidance
> on choosing the posterior summary that best serves your inferential
> goal, see
> [`vignette("posterior-summaries")`](https://joonho112.github.io/DPMirt/articles/posterior-summaries.md).

**Practical recommendation.** For applied researchers who suspect
non-normality, fit both priors and compare them using (1) WAIC for
item-response prediction and (2) a density plot of posterior mean
estimates for distributional recovery. When the population is truly
non-normal and reliability is at least moderate ($`\bar{w} \ge 0.7`$),
the density comparison will typically reveal the DPM’s advantage even
when WAIC values are similar.

## 6. Identification Strategies

| Strategy | Rasch | 2PL/3PL | How It Works |
|:---|:--:|:--:|:---|
| constrained_ability | Yes | Yes | Fix $`\theta \sim N(0, 1)`$; no hyperparameters |
| constrained_item | Yes (default) | Yes | Mean-center $`\beta`$ during MCMC (+ geom-mean $`\lambda = 1`$) |
| unconstrained | Yes | Yes (default) | No constraints; post-hoc rescaling via dpmirt_rescale() |

> **Paganin’s finding.** The unconstrained + post-hoc rescaling approach
> yields the most efficient MCMC sampler for DPM-IRT models.
> `constrained_ability` is **incompatible** with the DPM prior (DPMirt
> will raise an error) since the DPM’s purpose is to learn the ability
> distribution from data.

``` r

# Rasch with unconstrained + post-hoc rescaling
fit <- dpmirt(data, model = "rasch", prior = "normal",
              identification = "unconstrained")

# 2PL with in-model constrained_item centering
fit <- dpmirt(data, model = "2pl", prior = "normal",
              identification = "constrained_item")
```

**How to choose.** For DPM models, use the defaults
(`"constrained_item"` for Rasch, `"unconstrained"` for 2PL/3PL). For
Normal-prior models all three are valid; `"constrained_ability"` gives
the strongest $`N(0,1)`$ shrinkage.

## 7. Chain Continuation

When an initial run needs more iterations,
[`dpmirt_resume()`](https://joonho112.github.io/DPMirt/reference/dpmirt_resume.md)
continues the sampler from its current state without recompilation.

``` r

# Initial fit
fit <- dpmirt(data, model = "rasch", prior = "dpm",
              niter = 5000, nburnin = 1000, seed = 42)

# Trace looks non-stationary...
plot(fit, type = "trace")

# Resume: add 10,000 more iterations (no burn-in needed)
resumed <- dpmirt_resume(fit, niter_more = 10000)
```

With the step-by-step workflow:

``` r

samples <- dpmirt_sample(compiled, niter = 5000, nburnin = 1000)
resumed <- dpmirt_resume(compiled, niter_more = 10000)
```

> [`dpmirt_resume()`](https://joonho112.github.io/DPMirt/reference/dpmirt_resume.md)
> requires the compiled model’s C++ pointers to be alive (same R
> session). If pointers are expired, recompile.

The returned `dpmirt_samples` object can then be passed through
[`dpmirt_rescale()`](https://joonho112.github.io/DPMirt/reference/dpmirt_rescale.md)
(which returns a `dpmirt_fit`) and
[`dpmirt_estimates()`](https://joonho112.github.io/DPMirt/reference/dpmirt_estimates.md).

## 8. Model Comparison with WAIC

Every `dpmirt_fit` stores a WAIC (Watanabe–Akaike Information Criterion)
value, computed automatically during sampling. **Lower WAIC indicates
better out-of-sample prediction.**

``` r

comparison <- dpmirt_compare(fit_rasch_n, fit_rasch_dpm)
comparison
```

``` r

comparison <- dpmirt_compare(fit_rasch_n, fit_rasch_dpm)
knitr::kable(comparison, digits = 2, row.names = FALSE)
```

| model        |    waic | delta_waic |
|:-------------|--------:|-----------:|
| RASCH-dpm    | 6078.17 |       0.00 |
| RASCH-normal | 6079.03 |       0.86 |

The `delta_waic` column shows the difference from the best model. You
can compare any number of fits:

``` r

dpmirt_compare(fit_rasch_n, fit_rasch_dpm, fit_2pl)
```

> **Interpreting WAIC differences.**
>
> - $`|\Delta\text{WAIC}| < 2`$: Models essentially equivalent.
> - $`2 < |\Delta\text{WAIC}| < 10`$: Moderate preference.
> - $`|\Delta\text{WAIC}| > 10`$: Strong preference.

## 9. Data Formats

### Matrix Format (Default)

Binary matrix with persons in rows and items in columns. Missing
responses coded as `NA`.

``` r

head(sim$response[, 1:6])
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    1    0    1    1    0    1
#> [2,]    0    1    0    1    1    0
```

### Long Format

A data.frame with columns for person, item, and response. DPMirt detects
the format automatically and converts internally.

``` r

head(response_long)
#>   person item response
#> 1      1    1        1
#> 2      1    2        0
```

``` r

fit <- dpmirt(response_long, model = "rasch", prior = "normal")
```

> If the long-to-matrix conversion produces a very large or sparse
> matrix, verify dimensions with
> [`dpmirt_spec()`](https://joonho112.github.io/DPMirt/reference/dpmirt_spec.md)
> before running the full pipeline.

## 10. Visualizations

DPMirt provides 12 plot types via `plot(fit, type = ...)`. When
**ggplot2** is installed it is used automatically; otherwise base R
graphics are produced.

| Type            | Description                                         |
|:----------------|:----------------------------------------------------|
| density         | Kernel density of posterior mean theta              |
| items           | Item parameter estimates with posterior intervals   |
| trace           | Log-likelihood trace (convergence diagnostic)       |
| clusters        | Active cluster count trace and posterior (DPM only) |
| dp_density      | DP mixture density with credible band (DPM only)    |
| icc             | Item Characteristic Curves                          |
| wright_map      | Wright map: items and persons on common logit scale |
| parameter_trace | Trace plots for individual parameters               |
| caterpillar     | Caterpillar plot of sorted estimates with CIs       |
| density_compare | Overlaid densities from two fits                    |
| info            | Test and item information functions (2PL/3PL)       |
| pp_check        | Posterior predictive check                          |

``` r

plot(fit_rasch_dpm, type = "density")
```

![Posterior density of person
abilities.](models-and-workflow_files/figure-html/plot-density-guide-1.png)

Posterior density of person abilities.

``` r

plot(fit_rasch_dpm, type = "items")
```

![Item difficulty estimates with posterior
intervals.](models-and-workflow_files/figure-html/plot-items-guide-1.png)

Item difficulty estimates with posterior intervals.

``` r

plot(fit_rasch_dpm, type = "trace")
```

![Log-likelihood
trace.](models-and-workflow_files/figure-html/plot-trace-guide-1.png)

Log-likelihood trace.

``` r

plot(fit_rasch_dpm, type = "icc", items = 1:6)
```

![Item Characteristic Curves for items
1--6.](models-and-workflow_files/figure-html/plot-icc-guide-1.png)

Item Characteristic Curves for items 1–6.

``` r

plot(fit_rasch_dpm, type = "wright_map")
```

![Wright
map.](models-and-workflow_files/figure-html/plot-wright-guide-1.png)

Wright map.

## 11. Summary and What’s Next?

This vignette covered the full model and workflow space in DPMirt:

- **Two workflow modes**: one-step
  [`dpmirt()`](https://joonho112.github.io/DPMirt/reference/dpmirt.md)
  and the step-by-step pipeline for compile-once sampling.
- **Three IRT models**: Rasch, 2PL, and 3PL, each with Normal or DPM
  priors.
- **Two parameterizations** for 2PL/3PL: IRT and slope–intercept.
- **Three identification strategies**: constrained_ability,
  constrained_item, and unconstrained + post-hoc rescaling.
- **Chain continuation** with
  [`dpmirt_resume()`](https://joonho112.github.io/DPMirt/reference/dpmirt_resume.md).
- **Model comparison** with WAIC via
  [`dpmirt_compare()`](https://joonho112.github.io/DPMirt/reference/dpmirt_compare.md).
- **Twelve built-in plot types** for posterior visualization.

| Topic | Vignette | What |
|:---|:---|:---|
| Choosing an estimator | *Posterior Summaries* | In-depth PM vs CB vs GR comparison; shrinkage diagnostics; loss evaluation |
| Setting DPM priors | *Prior Elicitation* | Principled choice of the concentration parameter via DPprior |
| NIMBLE internals | *NIMBLE Internals* | Customizing samplers, NIMBLE compilation model, advanced tuning |
