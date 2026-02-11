# Under the Hood: NIMBLE Backend and Advanced Usage

## Overview

DPMirt is built on top of the NIMBLE probabilistic programming framework
(de Valpine et al., 2017). This vignette is for users who want to:

- Understand **how** DPMirt generates NIMBLE model code.
- Inspect and modify the **custom NIMBLE components** (distributions and
  samplers).
- Leverage the **compile-once, sample-many** pattern for efficient
  multi-chain or exploratory workflows.
- Understand the **post-hoc rescaling** that resolves identification
  indeterminacy.
- Reconstruct the **DP mixture density** from posterior samples.

Most users will never need this level of detail — the
[`dpmirt()`](https://joonho112.github.io/DPMirt/reference/dpmirt.md)
function handles everything automatically. But if you want to write
custom samplers, hook into the compiled MCMC, or adapt DPMirt for a new
model, this is where to start.

## Programmatic Code Generation

### Why Programmatic?

Paganin et al. (2023) provide over 20 separate NIMBLE code files, one
for each model–prior–identification combination. DPMirt instead
generates NIMBLE code programmatically via
[`dpmirt_spec()`](https://joonho112.github.io/DPMirt/reference/dpmirt_spec.md),
dispatching on three configuration axes:

| Axis           | Options                                              |
|:---------------|:-----------------------------------------------------|
| Model          | rasch, 2pl, 3pl                                      |
| Prior          | normal, dpm                                          |
| Identification | constrained_item, constrained_ability, unconstrained |

This yields $`3 \times 2 \times 3 = 18`$ potential combinations (not all
valid), all generated from a single specification call.

### Inspecting a Specification

The
[`dpmirt_spec()`](https://joonho112.github.io/DPMirt/reference/dpmirt_spec.md)
function is lightweight — it builds the NIMBLE code, constants, data,
initial values, and monitor configuration without compiling anything.
This is fast (under 1 second) and produces a fully inspectable object:

``` r

# Generate some example data
sim <- dpmirt_simulate(200, 25, model = "rasch", latent_shape = "bimodal",
                       seed = 42)

# Create a specification (no compilation)
spec <- dpmirt_spec(
  sim$response,
  model  = "rasch",
  prior  = "dpm",
  alpha_prior = c(a = 1, b = 3),
  base_measure = list(s2_mu = 2, nu1 = 2.01, nu2 = 1.01)
)

print(spec)
#> DPMirt Model Specification
#> ==========================
#> Model:            RASCH 
#> Prior:            dpm 
#> Identification:   constrained_item 
#> Persons (N):      200 
#> Items (I):        25 
#> Max clusters (M): 50 
#> Alpha prior:      Gamma(1, 3)
#> Monitors:         beta, alpha, zi, muTilde, s2Tilde, myLogProbAll, myLogProbSome, myLogLik 
#> Monitors2:        eta
```

### Inspecting the Generated Code

The `code` element is a `nimbleCode` object that you can print or
manipulate:

``` r

spec$code
#> {
#>     for (j in 1:N) {
#>         for (i in 1:I) {
#>             y[j, i] ~ dbern(pi[j, i])
#>             logit(pi[j, i]) <- eta[j] - beta[i]
#>         }
#>     }
#>     for (i in 1:I) {
#>         beta.tmp[i] ~ dnorm(0, var = sigma2_beta)
#>     }
#>     beta[1:I] <- beta.tmp[1:I] - mean(beta.tmp[1:I])
#>     zi[1:N] ~ dCRP(alpha, size = N)
#>     alpha ~ dgamma(a, b)
#>     for (j in 1:N) {
#>         eta[j] ~ dnorm(mu_j[j], var = s2_j[j])
#>         mu_j[j] <- muTilde[zi[j]]
#>         s2_j[j] <- s2Tilde[zi[j]]
#>     }
#>     for (m in 1:M) {
#>         muTilde[m] ~ dnorm(0, var = s2_mu)
#>         s2Tilde[m] ~ dinvgamma(nu1, nu2)
#>     }
#>     myLogProbAll ~ dnorm(0, 1)
#>     myLogProbSome ~ dnorm(0, 1)
#>     myLogLik ~ dnorm(0, 1)
#> }
```

### Inspecting Constants and Monitors

``` r

# Constants passed to the NIMBLE model
str(spec$constants)
#> List of 9
#>  $ N          : int 200
#>  $ I          : int 25
#>  $ sigma2_beta: num 3
#>  $ M          : int 50
#>  $ a          : Named num 1
#>   ..- attr(*, "names")= chr "a"
#>  $ b          : Named num 3
#>   ..- attr(*, "names")= chr "b"
#>  $ s2_mu      : num 2
#>  $ nu1        : num 2.01
#>  $ nu2        : num 1.01
```

``` r

# What gets tracked during MCMC
cat("Primary monitors:", paste(spec$monitors, collapse = ", "), "\n")
#> Primary monitors: beta, alpha, zi, muTilde, s2Tilde, myLogProbAll, myLogProbSome, myLogLik
cat("Thinned monitors:", paste(spec$monitors2, collapse = ", "), "\n")
#> Thinned monitors: eta
```

### Annotated Walkthrough: Rasch-DPM-Unconstrained

The generated code for a Rasch model with DPM prior and unconstrained
identification has this structure (shown here as pseudocode for
clarity):

    nimbleCode({
      # --- Likelihood ---
      for (i in 1:I) {
        for (j in 1:N) {
          y[j, i] ~ dbern(pi[j, i])
          logit(pi[j, i]) <- eta[j] - beta[i]      # Rasch ICC
        }
      }

      # --- Item parameters: free (rescaled post-hoc) ---
      for (i in 1:I) {
        beta[i] ~ dnorm(0, var = sigma2_beta)        # Prior on difficulties
      }

      # --- DPM prior for abilities via CRP ---
      zi[1:N] ~ dCRP(alpha, size = N)                # Cluster assignments
      alpha ~ dgamma(a, b)                           # Concentration parameter

      for (j in 1:N) {
        eta[j] ~ dnorm(mu_j[j], var = s2_j[j])      # Person ability
        mu_j[j]  <- muTilde[zi[j]]                   # Cluster mean
        s2_j[j]  <- s2Tilde[zi[j]]                   # Cluster variance
      }

      for (m in 1:M) {
        muTilde[m]  ~ dnorm(0, var = s2_mu)          # Base measure: mean
        s2Tilde[m]  ~ dinvgamma(nu1, nu2)            # Base measure: variance
      }

      # --- Log-probability monitoring nodes ---
      myLogProbAll  ~ dnorm(0, 1)                    # Replaced by logProb_summer
      myLogProbSome ~ dnorm(0, 1)
      myLogLik      ~ dnorm(0, 1)
    })

Key design choices:

- **CRP representation**: `dCRP(alpha, size = N)` uses the Chinese
  Restaurant Process to assign each person to a cluster. NIMBLE handles
  the stick-breaking conversion internally.
- **Truncation**: `M = 50` clusters (configurable) are pre-allocated.
  The CRP can use fewer but never more.
- **Dummy monitoring nodes**: `myLogProbAll`, `myLogProbSome`, and
  `myLogLik` are declared as dummy `dnorm(0,1)` nodes. Their default
  samplers are replaced by `logProb_summer` during MCMC configuration
  (see below).

## Custom NIMBLE Components

DPMirt registers three custom NIMBLE components, all adapted from
Paganin et al. (2023) and implemented using a lazy-initialization
pattern to avoid executing
[`nimbleFunction()`](https://rdrr.io/pkg/nimble/man/nimbleFunction.html)
at package load time.

| Component | Type | Purpose | When_Used |
|:---|:---|:---|:---|
| dBernoulliVector | Distribution | Vectorized Bernoulli likelihood for constrained_item identification | constrained_item with 2PL/3PL only |
| logProb_summer | Sampler | Log-probability monitoring: tracks logProb(all), logProb(params), logLik(data) | Always (all models, all priors) |
| sampler_centered | Sampler | Joint adaptive MH for correlated (log_lambda, gamma) in SI parameterization | 2PL/3PL with SI parameterization only |

Custom NIMBLE components registered by DPMirt. {.table}

### dBernoulliVector

For `constrained_item` identification in 2PL/3PL models, the likelihood
for each person’s entire response vector is evaluated jointly:

    y[j, 1:I] ~ dBernoulliVector(prob = pi[j, 1:I])

This is necessary because the constraint
`beta[1:I] <- beta.tmp[1:I] - mean(beta.tmp[1:I])` creates a
deterministic dependency across all items, requiring a vectorized
likelihood.

The density function sums element-wise Bernoulli log-probabilities:

``` r

# Simplified version of the registered nimbleFunction
dBernoulliVector <- function(x, prob, log = 0) {
  logProb <- sum(dbinom(x, size = 1, prob = prob, log = TRUE))
  if (log) return(logProb) else return(exp(logProb))
}
```

### logProb_summer

This pseudo-sampler does not propose new values — it computes and stores
the log-probability of a set of nodes at each MCMC iteration. Three
instances are configured:

| Target node | What it tracks | nodeList |
|:---|:---|:---|
| `myLogProbAll` | Total model log-probability | All nodes |
| `myLogProbSome` | Parameter log-probability | `c("beta", "eta")` or `c("gamma", "lambda", "eta")` |
| `myLogLik` | Data log-likelihood | `"y"` |

The `myLogLik` trace is used for: - Visual convergence checking (trace
plots). - WAIC computation (when NIMBLE’s built-in WAIC is not
available).

### sampler_centered

The centered sampler is critical for efficient mixing in the
slope-intercept (SI) parameterization of 2PL and 3PL models. It
addresses the strong posterior correlation between $`\log\lambda_i`$ and
$`\gamma_i`$ by making a joint proposal:

1.  **Propose** $`\log\lambda'_i = \log\lambda_i + \epsilon`$, where
    $`\epsilon \sim N(0, \text{scale}^2)`$.
2.  **Center** $`\gamma'_i = \gamma_i + \bar\eta \cdot
    (\exp(\log\lambda_i) - \exp(\log\lambda'_i))`$.
3.  **Accept/reject** via Metropolis-Hastings.

The centering mean $`\bar\eta`$ is updated each iteration to
`mean(eta)`. Adaptive scaling targets an acceptance rate of 0.44
(optimal for univariate random-walk MH).

> **When is it used?** Only for 2PL/3PL models with SI parameterization
> and `identification = "unconstrained"` or `"constrained_ability"`. For
> the IRT parameterization, NIMBLE’s default samplers are sufficient
> because the correlation structure is different.

The
[`dpmirt_compile()`](https://joonho112.github.io/DPMirt/reference/dpmirt_compile.md)
function automatically enables the centered sampler when
`use_centered_sampler = "auto"` (the default).

## Compile-Once, Sample-Many Pattern

NIMBLE compilation translates the model and MCMC configuration into C++
code, compiles it, and loads the resulting shared library. This is the
most expensive step in the pipeline.

DPMirt separates compilation from sampling, allowing you to compile once
and then run multiple chains, extend runs, or experiment with different
MCMC lengths — all without recompiling.

### Architecture

    dpmirt_spec()  ──>  dpmirt_compile()  ──>  dpmirt_sample()
      [< 1 sec]        [30-120 sec]           [10-60 sec per call]
                             │
                             ├── dpmirt_sample(seed = 1)
                             ├── dpmirt_sample(seed = 2)
                             ├── dpmirt_sample(seed = 3)
                             └── dpmirt_sample(seed = 4)

### Timing Table

| Operation          | Typical_Time | Bottleneck           |
|:-------------------|:-------------|:---------------------|
| dpmirt_spec()      | \< 1 sec     | Code generation      |
| dpmirt_compile()   | 30-120 sec   | C++ compilation      |
| dpmirt_sample(10K) | 10-60 sec    | MCMC iterations      |
| dpmirt_rescale()   | \< 1 sec     | Matrix operations    |
| dpmirt_estimates() | \< 1 sec     | Quantile computation |

Typical timing for each pipeline stage (Rasch, N=200, I=25). {.table}

### Step-by-Step Workflow

``` r

# Step 1: Specification (fast)
spec <- dpmirt_spec(sim$response, model = "rasch", prior = "dpm")

# Step 2: Compilation (expensive --- do this once)
compiled <- dpmirt_compile(spec, verbose = TRUE)

# Step 3: Multiple chains (fast --- reuse compiled object)
chains <- lapply(1:4, function(i) {
  dpmirt_sample(compiled, niter = 10000, nburnin = 2000, seed = i)
})

# Step 4: Chain continuation (extend without recompiling)
# Option A: low-level via dpmirt_sample with reset = FALSE
extended <- dpmirt_sample(compiled, niter = 5000, nburnin = 0,
                           reset = FALSE)

# Option B: higher-level via dpmirt_resume (preferred)
# extended <- dpmirt_resume(compiled, niter = 5000)
```

### Custom Sampler Configuration

For advanced users,
[`dpmirt_compile()`](https://joonho112.github.io/DPMirt/reference/dpmirt_compile.md)
accepts a custom sampler configuration function:

``` r

# Define a custom sampler configuration function
my_config <- function(conf, model, spec) {
  # Remove all default samplers for eta
  for (j in seq_len(spec$config$N)) {
    conf$removeSamplers(paste0("eta[", j, "]"))
  }

  # Add slice samplers instead of the default RW
  for (j in seq_len(spec$config$N)) {
    conf$addSampler(
      target = paste0("eta[", j, "]"),
      type   = "slice"
    )
  }

  conf
}

compiled <- dpmirt_compile(spec, sampler_config = my_config)
```

> **Serialization caveat:** Compiled NIMBLE objects contain external C++
> pointers that **cannot be serialized** across R sessions. If you call
> `saveRDS(compiled, ...)` and then `readRDS(...)` in a new session, the
> pointers will be dead and all MCMC operations will fail. Instead, save
> the `dpmirt_spec` object and recompile:
>
> ``` r
>
> # CORRECT: Save spec, recompile in new session
> saveRDS(spec, "my_spec.rds")
> # ... in new R session ...
> spec <- readRDS("my_spec.rds")
> compiled <- dpmirt_compile(spec)
> ```

### Accessing NIMBLE Objects Directly

The compiled object provides direct access to the underlying NIMBLE
objects for advanced manipulation:

``` r

# The compiled NIMBLE model (C object)
compiled$Cmodel

# The compiled MCMC (C object)
compiled$Cmcmc

# List all node names
compiled$Cmodel$getNodeNames()

# Get current values of a parameter
compiled$Cmodel$eta

# Calculate log-probability of the current model state
compiled$Cmodel$calculate()

# Inspect sampler configuration
compiled$mcmc_conf$printSamplers()
```

## Post-hoc Rescaling Implementation

IRT models have identification indeterminacy: without constraints, the
likelihood is invariant to certain transformations of the parameters.
DPMirt’s default approach is to leave the model unconstrained during
MCMC and apply rescaling to the posterior samples afterwards.

### Rasch: Location Shift Only

The Rasch model has only location indeterminacy. For each MCMC iteration
$`s`$, the rescaling is:

``` math
\beta^*_{i,s} = \beta_{i,s} - \bar{\beta}_s, \qquad
\theta^*_{p,s} = \theta_{p,s} - \bar{\beta}_s
```

where $`\bar{\beta}_s = \frac{1}{I}\sum_{i=1}^I \beta_{i,s}`$.

**Invariant check:** After rescaling, $`|\bar{\beta}^*_s| < \epsilon`$
for every iteration.

### 2PL IRT: Location + Scale

The 2PL model has both location and scale indeterminacy. For each
iteration $`s`$:

``` math
\text{location}_s = \bar{\beta}_s, \qquad
\text{scale}_s = \left(\prod_{i=1}^I \lambda_{i,s}\right)^{-1/I}
```

``` math
\beta^*_{i,s} = \frac{\beta_{i,s} - \text{location}_s}{\text{scale}_s}, \qquad
\lambda^*_{i,s} = \lambda_{i,s} \cdot \text{scale}_s, \qquad
\theta^*_{p,s} = \frac{\theta_{p,s} - \text{location}_s}{\text{scale}_s}
```

**Invariant checks:** - $`|\bar{\beta}^*_s| < \epsilon`$ -
$`|\text{geom\_mean}(\lambda^*_s) - 1| < \epsilon`$

### 2PL SI: Different Sign Convention

For the slope-intercept parameterization
($`\text{logit}(p) = \lambda\theta + \gamma`$), the location shift has
opposite sign because $`\gamma = -\lambda\beta`$:

``` math
\text{location}_s = \frac{\sum_i \gamma_{i,s}}{\sum_i \lambda_{i,s}}
```

``` math
\gamma^*_{i,s} = \gamma_{i,s} - \lambda_{i,s} \cdot \text{location}_s, \qquad
\theta^*_{p,s} = \frac{\theta_{p,s} + \text{location}_s}{\text{scale}_s}
```

Note the $`+`$ sign for theta rescaling (vs. $`-`$ in IRT
parameterization).

### Why Post-hoc?

Post-hoc rescaling has two advantages over in-MCMC constraints:

1.  **Simpler MCMC**: The unconstrained model has a simpler posterior
    geometry, avoiding the need for constrained samplers.
2.  **Flexibility**: Different identification conventions can be applied
    to the same posterior samples without re-running the MCMC.

The disadvantage is that the raw (unrescaled) samples may have poor
mixing if the location/scale dimensions are poorly identified. In
practice, DPMirt’s default priors provide sufficient regularization.

## DP Density Reconstruction

For DPM models, DPMirt can reconstruct the posterior predictive density
of the latent trait. This follows Paganin et al.’s
[`getSamplesDPmeasure()`](https://rdrr.io/pkg/nimble/man/getSamplesDPmeasure.html)
approach.

### Algorithm

For each post-burn-in MCMC iteration $`s`$:

1.  **Extract** the stick-breaking weights
    $`w_{s,1}, \ldots, w_{s,K_s}`$ and atoms
    $`(\tilde\mu_{s,k}, \tilde\sigma^2_{s,k})`$ from the CRP posterior
    via NIMBLE’s
    [`getSamplesDPmeasure()`](https://rdrr.io/pkg/nimble/man/getSamplesDPmeasure.html).

2.  **Evaluate** the finite mixture density on a grid:

``` math
f_s(x) = \sum_{k=1}^{K_s} w_{s,k} \cdot \phi\!\left(x;\, \tilde\mu_{s,k},\, \tilde\sigma^2_{s,k}\right)
```

where $`\phi(x; \mu, \sigma^2)`$ is the Normal density.

3.  **Average** over iterations for the posterior predictive mean:

``` math
\hat{f}(x) = \frac{1}{S}\sum_{s=1}^S f_s(x)
```

4.  **Pointwise credible bands**: Take quantiles across iterations at
    each grid point.

### Rescaling Adjustment

For unconstrained models, the DP density is on the raw (unrescaled)
scale. To evaluate it on the rescaled theta scale, the grid is shifted
by the iteration-specific location shift:

``` math
\hat{f}_{\text{rescaled}}(x) = f_s(x + \text{location}_s)
```

The Jacobian is 1 for a location shift, so no density adjustment is
needed.

### Using `dpmirt_dp_density()`

``` r

# Compute DP density from a fitted DPM model
dp_dens <- dpmirt_dp_density(
  fit,
  grid               = seq(-6, 6, length.out = 500),
  credible_interval  = 0.95,
  apply_rescaling    = TRUE
)

# Access the results
dp_dens$grid            # Evaluation points
dp_dens$density_mean    # Posterior mean density
dp_dens$density_lower   # 2.5th percentile (lower band)
dp_dens$density_upper   # 97.5th percentile (upper band)
dp_dens$density_samples # Full matrix: niter x n_grid
```

### Visualization

``` r

# Plot with credible band
plot(dp_dens$grid, dp_dens$density_mean, type = "l", lwd = 2,
     col = pal$semiparametric,
     xlab = expression(theta), ylab = "Density",
     main = "DP Mixture Posterior Density")

polygon(c(dp_dens$grid, rev(dp_dens$grid)),
        c(dp_dens$density_lower, rev(dp_dens$density_upper)),
        col = adjustcolor(pal$semiparametric, alpha.f = 0.2),
        border = NA)

# Reference: N(0,1)
curve(dnorm, add = TRUE, lty = 2, col = pal$reference, lwd = 1.5)

legend("topright",
       legend = c("DP posterior mean", "95% credible band", "N(0,1)"),
       col = c(pal$semiparametric, pal$semiparametric, pal$reference),
       lwd = c(2, NA, 1.5), lty = c(1, NA, 2),
       fill = c(NA, adjustcolor(pal$semiparametric, 0.2), NA),
       border = c(NA, NA, NA), bty = "n")
```

> **Note on computation cost:**
> [`dpmirt_dp_density()`](https://joonho112.github.io/DPMirt/reference/dpmirt_dp_density.md)
> internally rebuilds and recompiles a NIMBLE model to call
> [`getSamplesDPmeasure()`](https://rdrr.io/pkg/nimble/man/getSamplesDPmeasure.html).
> This adds 30–60 seconds. The
> [`dpmirt()`](https://joonho112.github.io/DPMirt/reference/dpmirt.md)
> function runs this automatically when `compute_dp_density = TRUE` (the
> default for DPM models). If you want to skip it, set
> `compute_dp_density = FALSE` and call `dpmirt_dp_density(fit)` later.

## MCMC Pipeline Internals

The following diagram shows the complete data flow from raw data to the
final `dpmirt_fit` object:

    data (N x I matrix)
      |
      v
    dpmirt_spec()
      |-- code:       nimbleCode object
      |-- constants:  N, I, M, hyperparams
      |-- data:       list(y = ...)
      |-- inits:      smart init from k-means
      |-- monitors:   c("beta", "alpha", "zi", ...)
      |-- monitors2:  c("eta")
      v
    dpmirt_compile()
      |-- nimbleModel()      --> Rmodel
      |-- configureMCMC()    --> add logProb_summer, centered sampler
      |-- buildMCMC()        --> Rmcmc
      |-- compileNimble()    --> Cmodel, Cmcmc
      v
    dpmirt_sample()
      |-- runMCMC()          --> samples (niter x params), samples2 (niter x N)
      v
    dpmirt_rescale()
      |-- .rescale_rasch()   or .rescale_irt() or .rescale_si()
      |-- theta_samp, beta_samp, lambda_samp (all rescaled)
      v
    dpmirt_fit object
      |-- theta_samp, beta_samp, ...
      |-- ess, waic, loglik_trace
      |-- cluster_info (DPM only)
      |-- dp_density (DPM only, optional)

### Initial Values: K-means Strategy

DPMirt uses a k-means clustering strategy (following Paganin et al.) for
initializing DPM cluster assignments:

1.  Compute standardized sum scores for each person.
2.  Run k-means with $`k = \min(5, \lfloor N/4 \rfloor)`$ clusters.
3.  Use the cluster assignments as initial values for `zi[1:N]`.

This provides a reasonable starting point that reduces the burn-in
needed for the CRP to find a good configuration.

### Monitor Thinning

DPMirt supports differential thinning:

- **monitors** (beta, lambda, alpha, zi, muTilde, s2Tilde, logprob):
  Thinned at rate `thin` (default 1 = every iteration).
- **monitors2** (eta/theta): Thinned at rate `thin2` (default = `thin`).

Since theta has $`N`$ dimensions (one per person), it can dominate
memory for large samples. Setting `thin2 > thin` reduces storage:

``` r

# Thin theta more aggressively than item parameters
fit <- dpmirt(
  data,
  model = "rasch",
  prior = "dpm",
  niter = 20000,
  thin  = 1,     # Keep every iteration for items
  thin2 = 5      # Keep every 5th iteration for theta
)
```

## Extending DPMirt

### Writing a Custom Sampler

To replace or supplement DPMirt’s default samplers:

``` r

# Example: Block sampler for all item difficulties
my_block_beta_config <- function(conf, model, spec) {
  I <- spec$config$I

  # Remove individual beta samplers
  for (i in seq_len(I)) {
    conf$removeSamplers(paste0("beta[", i, "]"))
  }

  # Add a single block (AF_slice) sampler
  conf$addSampler(
    target = paste0("beta[1:", I, "]"),
    type   = "AF_slice"
  )

  conf
}

compiled <- dpmirt_compile(spec, sampler_config = my_block_beta_config)
```

### Adding a New Model

To add a new model type (e.g., a multidimensional IRT model), you would:

1.  Add a new code builder function in `R/model_spec.R` (e.g.,
    `.build_mirt_code()`).
2.  Register it in the `.build_nimble_code()` dispatcher.
3.  Update `.generate_inits()` for the new parameter structure.
4.  Update `.setup_monitors()` for new parameters.
5.  Add rescaling logic in `R/rescale.R` if needed.

### Extracting Raw NIMBLE Objects

For maximum flexibility, you can extract the underlying NIMBLE objects
from a `dpmirt_compiled` object:

``` r

# Get the compiled model
Cmodel <- compiled$Cmodel

# Query model structure
Cmodel$getNodeNames()
Cmodel$getDependencies("beta[1]")
Cmodel$getLogProb("beta")

# Get the compiled MCMC
Cmcmc <- compiled$Cmcmc

# Run MCMC manually
Cmcmc$run(5000, reset = TRUE)

# Extract samples directly
samples <- as.matrix(Cmcmc$mvSamples)
samples2 <- as.matrix(Cmcmc$mvSamples2)
```

## Troubleshooting

### Common Issues

| Symptom | Likely Cause | Solution |
|:---|:---|:---|
| “Compiled model is no longer valid” | Serialized/loaded compiled object | Recompile from spec |
| Very slow compilation (\> 5 min) | Large N or I | Expected; reduce M for DPM |
| Poor mixing for alpha | Tight Gamma prior | Use wider prior or DP-diffuse |
| All persons in one cluster | alpha too small | Increase mu_K or use “low” confidence |
| WAIC computation failed | NIMBLE version mismatch | Update nimble; use compute_waic=FALSE |

### Checking Sampler Configuration

``` r

# After compilation, inspect what samplers are assigned
compiled$mcmc_conf$printSamplers()

# Count samplers by type
sampler_types <- sapply(
  seq_len(length(compiled$mcmc_conf$getSamplers())),
  function(i) compiled$mcmc_conf$getSamplers()[[i]]$name
)
table(sampler_types)
```

## What’s Next?

| Vignette            | Why read it                                          |
|:--------------------|:-----------------------------------------------------|
| *Quick Start*       | End-to-end pipeline without worrying about internals |
| *Prior Elicitation* | Principled $`\alpha`$ choice via DPprior             |
| *Simulation Study*  | Full factorial comparison of methods and priors      |

## References

de Valpine, P., Turek, D., Paciorek, C. J., Anderson-Bergman, C., Lang,
D. T., & Bodik, R. (2017). Programming with models: Writing statistical
algorithms for general model structures with NIMBLE. *Journal of
Computational and Graphical Statistics*, 26(2), 403–413.

Paganin, S., Paciorek, C. J., Wehrhahn, C., Rodríguez, A., Rabe-Hesketh,
S., & de Valpine, P. (2023). Computational strategies and estimation
performance with Bayesian semiparametric item response theory models.
*Journal of Educational and Behavioral Statistics*, 48(2), 147–188.
<https://doi.org/10.3102/10769986221136105>

NIMBLE Development Team. (2024). NIMBLE: MCMC, particle filtering, and
programmable hierarchical modeling. R package version 1.2.1.
<https://r-nimble.org>
