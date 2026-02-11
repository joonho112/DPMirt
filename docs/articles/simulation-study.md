# Simulation Study: Evaluating Prior Models and Posterior Summaries

## Overview

Simulation studies are the primary tool for evaluating IRT methods
because they provide access to the ground truth — the true person
abilities and item parameters that generated the data. DPMirt provides
an integrated simulation framework that lets you:

1.  **Simulate** data under known conditions using IRTsimrel’s EQC
    calibration (or a built-in fallback).
2.  **Fit** models under competing prior specifications (Normal
    vs. DPM).
3.  **Estimate** person abilities using three posterior summary methods
    (PM, CB, GR).
4.  **Evaluate** performance using multiple loss functions (MSEL, MSELR,
    KS).

This vignette walks through a complete simulation study design, executes
a single-condition deep dive, interprets the results, and provides
templates for scaling to a full factorial study.

## Simulation Design

The APM manuscript (Lee & Wind) uses a $`3 \times 2 \times 3`$ factorial
design crossing three factors:

### Factor 1: Latent Distribution Shape ($`G`$)

| Shape | Description | Generating mechanism |
|:---|:---|:---|
| Normal | Standard Gaussian | $`\theta \sim N(0, 1)`$ |
| Bimodal | Two-group mixture | $`0.5 \cdot N(-1.5, 0.5) + 0.5 \cdot N(1.5, 0.5)`$ |
| Skew | Right-skewed | Shifted exponential: $`\text{Exp}(1) - 1`$ |

The **Normal** condition is the parametric model’s “home turf” where DPM
offers little advantage. The **Bimodal** and **Skew** conditions
represent departures where flexible priors should excel.

### Factor 2: Sample Size ($`N`$)

| $`N`$ | Interpretation                                        |
|:------|:------------------------------------------------------|
| 50    | Small sample: high uncertainty, strong shrinkage      |
| 200   | Moderate sample: reliable estimation for most methods |

### Factor 3: Marginal Reliability ($`\bar{w}`$)

Reliability $`\bar{w}`$ is the most influential factor in the
simulation. It controls the information per person and therefore the
degree of posterior shrinkage.

| w-bar | ~Items (Rasch) | Shrinkage | Separation |
|:-----:|:--------------:|:---------:|:----------:|
| 0.50  |      ~10       |  Severe   |    1.0     |
| 0.70  |      ~23       | Moderate  |    1.5     |
| 0.90  |      ~61       |   Mild    |    3.0     |

Reliability levels and their implications for IRT estimation. {.table}

> **Why reliability dominates:** At $`\bar{w} = 0.50`$, each person’s
> posterior is heavily shrunk toward the prior — there is little room
> for any method to differ from another. At $`\bar{w} = 0.90`$,
> posteriors are tightly concentrated around data-driven estimates,
> amplifying the impact of prior misspecification.

### Full Design Matrix

The full factorial design has $`3 \times 2 \times 3 = 18`$ conditions.
Each condition is replicated (e.g., 100 times) for Monte Carlo
stability.

| Condition | G       |  N  | w_bar |
|:---------:|:--------|:---:|:-----:|
|     1     | Normal  | 50  |  0.5  |
|     2     | Bimodal | 50  |  0.5  |
|     3     | Skew    | 50  |  0.5  |
|     4     | Normal  | 200 |  0.5  |
|     5     | Bimodal | 200 |  0.5  |
|     6     | Skew    | 200 |  0.5  |
|     7     | Normal  | 50  |  0.7  |
|     8     | Bimodal | 50  |  0.7  |
|     9     | Skew    | 50  |  0.7  |
|    10     | Normal  | 200 |  0.7  |
|    11     | Bimodal | 200 |  0.7  |
|    12     | Skew    | 200 |  0.7  |
|    13     | Normal  | 50  |  0.9  |
|    14     | Bimodal | 50  |  0.9  |
|    15     | Skew    | 50  |  0.9  |
|    16     | Normal  | 200 |  0.9  |
|    17     | Bimodal | 200 |  0.9  |
|    18     | Skew    | 200 |  0.9  |

Full 3 x 2 x 3 factorial design (18 conditions). {.table}

### IRTsimrel Integration

DPMirt uses the **IRTsimrel** package (when available) to achieve
precise reliability targeting via Empirical Quadrature Calibration
(EQC). The key function is
[`dpmirt_simulate()`](https://joonho112.github.io/DPMirt/reference/dpmirt_simulate.md):

``` r

# Simulate bimodal data with target reliability 0.8
sim <- dpmirt_simulate(
  n_persons    = 200,
  n_items      = 25,
  model        = "rasch",
  target_rho   = 0.8,
  latent_shape = "bimodal",
  seed         = 42
)

sim
#> DPMirt Simulated Data
#> =====================
#> Model:         RASCH 
#> Persons:       200 
#> Items:         25 
#> Distribution:  bimodal 
#> Method:        irtsimrel 
#> Target rho:    0.8 
#> KR-20:         0.81 
#> EQC c*:        0.935 
#> EQC rho:       0.8
```

The returned `dpmirt_sim` object reports:

- **Method**: “irtsimrel” (EQC-calibrated) or “fallback” (Paganin-style)
- **KR-20**: Empirical reliability of the generated data
- **EQC c**\*: The calibration constant that maps target reliability to
  the number of items

When IRTsimrel is not installed,
[`dpmirt_simulate()`](https://joonho112.github.io/DPMirt/reference/dpmirt_simulate.md)
falls back to a simpler simulation with evenly spaced item difficulties
and no reliability targeting.

### Visualizing the Simulated Data

``` r

theta_dens <- density(sim$theta)
plot(theta_dens, main = "True theta distribution (bimodal)",
     xlab = expression(theta), lwd = 2)
rug(sim$theta, col = adjustcolor("black", alpha.f = 0.2))
```

![Density of true theta values from the bimodal
simulation.](simulation-study_files/figure-html/sim-theta-density-1.png)

Density of true theta values from the bimodal simulation.

``` r

sum_scores <- rowSums(sim$response)
hist(sum_scores, breaks = 20, main = "Sum score distribution",
     xlab = "Sum score", col = adjustcolor(pal$parametric, alpha.f = 0.5),
     border = "white")
```

![Distribution of sum scores from the simulated response
matrix.](simulation-study_files/figure-html/sim-sumscores-1.png)

Distribution of sum scores from the simulated response matrix.

## Analysis Pipeline

The full analysis pipeline for a single simulation condition consists of
six steps. Below we show the complete code you would run in a live
session, then load pre-computed results for the remaining analysis.

### Step 1: Simulate Data

``` r

sim <- dpmirt_simulate(
  n_persons    = 200,
  n_items      = 25,
  model        = "rasch",
  target_rho   = 0.8,
  latent_shape = "bimodal",
  seed         = 42
)
```

### Step 2: Fit Models

``` r

# --- Fit Normal prior model (~2 min compilation + sampling) ---
fit_normal <- dpmirt(
  sim$response,
  model   = "rasch",
  prior   = "normal",
  niter   = 10000,
  nburnin = 2000,
  seed    = 100
)

# --- Fit DPM prior model (~2 min compilation + sampling) ---
fit_dpm <- dpmirt(
  sim$response,
  model      = "rasch",
  prior      = "dpm",
  mu_K       = 5,
  confidence = "medium",
  niter      = 10000,
  nburnin    = 2000,
  seed       = 200
)
```

### Step 3: Extract Estimates

``` r

# Compute PM, CB, GR for both models
est_normal <- dpmirt_estimates(fit_normal, methods = c("pm", "cb", "gr"))
est_dpm    <- dpmirt_estimates(fit_dpm,    methods = c("pm", "cb", "gr"))
```

### Step 4: Evaluate Losses

``` r

# Compare to true theta
loss_normal <- dpmirt_loss(est_normal, true_theta = sim$theta,
                            true_beta = sim$beta)
loss_dpm    <- dpmirt_loss(est_dpm,    true_theta = sim$theta,
                            true_beta = sim$beta)
```

### Step 5: Combine Results

``` r

loss_normal$prior <- "Normal"
loss_dpm$prior    <- "DPM"
combined_loss <- rbind(loss_normal, loss_dpm)
```

### Loading Pre-Computed Results

For this vignette, all fitting results are pre-computed:

``` r

# Load pre-computed simulation data
sim_data <- readRDS(find_extdata("vignette_sim_bimodal.rds"))

# Load estimates comparison (all 6 method-prior combos)
est_comparison <- readRDS(find_extdata("vignette_estimates_comparison.rds"))

# Load loss results
loss_results <- readRDS(find_extdata("vignette_loss_results.rds"))
```

## Single Condition Deep Dive

We focus on the most informative condition: **Bimodal**, $`N = 200`$,
$`\bar{w} \approx 0.8`$. This is where the DPM prior should show the
clearest advantage over the Normal prior, because the true latent
distribution violates the Gaussian assumption and the data are
informative enough for the posterior to reflect that violation.

### Posterior Density Overlay

``` r

theta_true <- est_comparison$true_theta
N <- length(theta_true)

df_dens <- data.frame(
  value = c(theta_true,
            est_comparison$normal$theta$theta_pm,
            est_comparison$normal$theta$theta_cb,
            est_comparison$normal$theta$theta_gr,
            est_comparison$dpm$theta$theta_pm,
            est_comparison$dpm$theta$theta_cb,
            est_comparison$dpm$theta$theta_gr),
  Method = factor(rep(c("True", "Normal-PM", "Normal-CB", "Normal-GR",
                         "DPM-PM", "DPM-CB", "DPM-GR"), each = N),
                  levels = c("True", "Normal-PM", "Normal-CB", "Normal-GR",
                             "DPM-PM", "DPM-CB", "DPM-GR"))
)

method_colors <- c(
  "True"      = "black",
  "Normal-PM" = pal$normal_pm, "Normal-CB" = pal$normal_cb,
  "Normal-GR" = pal$normal_gr,
  "DPM-PM"    = pal$dpm_pm, "DPM-CB" = pal$dpm_cb,
  "DPM-GR"    = pal$dpm_gr
)
method_lty <- c("True" = "solid",
                "Normal-PM" = "solid", "Normal-CB" = "dashed",
                "Normal-GR" = "dotted",
                "DPM-PM"    = "solid", "DPM-CB"    = "dashed",
                "DPM-GR"    = "dotted")

ggplot(df_dens, aes(x = value, colour = Method, linetype = Method)) +
  geom_density(linewidth = 0.9, fill = NA) +
  scale_colour_manual(values = method_colors) +
  scale_linetype_manual(values = method_lty) +
  labs(title = "Posterior density comparison: 6 methods vs. truth",
       subtitle = "Bimodal population, 25 items, 200 persons",
       x = expression(theta), y = "Density") +
  theme_bw() +
  theme(legend.position = "right")
```

![Density estimates from six method-prior combinations overlaid with the
true latent density. DPM-based methods (warm colors) better capture the
bimodal
shape.](simulation-study_files/figure-html/density-overlay-1.png)

Density estimates from six method-prior combinations overlaid with the
true latent density. DPM-based methods (warm colors) better capture the
bimodal shape.

> **What to look for:** Under the Normal prior, all three estimators
> (PM, CB, GR) tend to “fill in” the valley between the two modes
> because the Gaussian assumption forces a unimodal shape. The DPM-based
> estimators, especially CB and GR, better preserve the bimodal
> structure.

### Shrinkage Comparison

``` r

par(mfrow = c(1, 2))

# Normal-PM
plot(theta_true, est_comparison$normal$theta$theta_pm,
     pch = 16, cex = 0.6, col = adjustcolor(pal$normal_pm, 0.6),
     xlab = expression(theta[true]), ylab = expression(hat(theta)[PM]),
     main = "Normal-PM", asp = 1)
abline(0, 1, col = pal$reference, lwd = 1.5, lty = 2)

# DPM-PM
plot(theta_true, est_comparison$dpm$theta$theta_pm,
     pch = 16, cex = 0.6, col = adjustcolor(pal$dpm_pm, 0.6),
     xlab = expression(theta[true]), ylab = expression(hat(theta)[PM]),
     main = "DPM-PM", asp = 1)
abline(0, 1, col = pal$reference, lwd = 1.5, lty = 2)
```

![Posterior mean (PM) estimates vs. true theta. Normal-PM shows stronger
shrinkage toward zero; DPM-PM preserves extreme values
better.](simulation-study_files/figure-html/shrinkage-plot-1.png)

Posterior mean (PM) estimates vs. true theta. Normal-PM shows stronger
shrinkage toward zero; DPM-PM preserves extreme values better.

``` r


par(mfrow = c(1, 1))
```

### Loss Function Results

Three complementary loss functions capture different aspects of
estimator quality:

| Loss | Formula | Measures |
|:---|:---|:---|
| **MSEL** | $`\frac{1}{N}\sum_p(\hat{\theta}_p - \theta_p)^2`$ | Individual accuracy (Goal 1) |
| **MSELR** | $`\frac{1}{N}\sum_p(R(\hat{\theta}_p)/N - R(\theta_p)/N)^2`$ | Ranking accuracy (Goal 2) |
| **KS** | $`\max_t\lvert F_{\hat{\theta}}(t) - F_\theta(t)\rvert`$ | Distributional fidelity (Goal 3) |

``` r

# loss_results is a list with $normal and $dpm data.frames,
# each with columns: parameter, method, msel, mselr, ks

# Combine into one table with a prior column
loss_normal <- loss_results$normal
loss_normal$prior <- "Normal"
loss_dpm <- loss_results$dpm
loss_dpm$prior <- "DPM"
loss_combined <- rbind(loss_normal, loss_dpm)

# Filter to theta and format
loss_display <- loss_combined[loss_combined$parameter == "theta", ]
loss_display <- loss_display[order(loss_display$prior, loss_display$method), ]

loss_display$Label <- paste0(loss_display$prior, "-",
                              toupper(loss_display$method))
loss_display <- loss_display[, c("Label", "msel", "mselr", "ks")]
names(loss_display) <- c("Method", "MSEL", "MSELR", "KS")
loss_display$MSEL  <- round(loss_display$MSEL, 4)
loss_display$MSELR <- round(loss_display$MSELR, 5)
loss_display$KS    <- round(loss_display$KS, 4)

knitr::kable(loss_display, row.names = FALSE,
             caption = "Loss function results: Bimodal, N=200, w-bar=0.8.",
             align = "lccc")
```

| Method    |  MSEL  |  MSELR  |  KS   |
|:----------|:------:|:-------:|:-----:|
| DPM-CB    | 0.3029 | 0.02499 | 0.080 |
| DPM-GR    | 0.3136 | 0.02487 | 0.055 |
| DPM-PM    | 0.2544 | 0.02499 | 0.095 |
| Normal-CB | 0.3010 | 0.02460 | 0.090 |
| Normal-GR | 0.3057 | 0.02439 | 0.060 |
| Normal-PM | 0.2544 | 0.02460 | 0.095 |

Loss function results: Bimodal, N=200, w-bar=0.8. {.table}

### Interpreting the Loss Table

``` r

# Identify best method for each loss
best_msel  <- loss_display$Method[which.min(loss_display$MSEL)]
best_mselr <- loss_display$Method[which.min(loss_display$MSELR)]
best_ks    <- loss_display$Method[which.min(loss_display$KS)]

cat("Best for MSEL  (individual accuracy):", best_msel, "\n")
#> Best for MSEL  (individual accuracy): DPM-PM
cat("Best for MSELR (ranking accuracy):   ", best_mselr, "\n")
#> Best for MSELR (ranking accuracy):    Normal-GR
cat("Best for KS    (distributional):     ", best_ks, "\n")
#> Best for KS    (distributional):      DPM-GR
```

### Visual Loss Comparison

``` r

method_colors_vec <- c(
  "Normal-CB" = pal$normal_cb, "Normal-GR" = pal$normal_gr,
  "Normal-PM" = pal$normal_pm,
  "DPM-CB"    = pal$dpm_cb, "DPM-GR" = pal$dpm_gr,
  "DPM-PM"    = pal$dpm_pm
)

loss_long <- data.frame(
  Method = rep(loss_display$Method, 3),
  Metric = rep(c("MSEL", "MSELR", "KS"), each = nrow(loss_display)),
  Loss   = c(loss_display$MSEL, loss_display$MSELR, loss_display$KS)
)
loss_long$Method <- factor(loss_long$Method, levels = loss_display$Method)
loss_long$Metric <- factor(loss_long$Metric, levels = c("MSEL", "MSELR", "KS"))

ggplot(loss_long, aes(x = Method, y = Loss, fill = Method)) +
  geom_col(width = 0.7) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = method_colors_vec, guide = "none") +
  labs(title = "Loss comparison: 6 method-prior combinations",
       y = "Loss (lower is better)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        strip.text = element_text(face = "bold", size = 11))
```

![Bar chart of loss values across six methods. Lower is better for all
three metrics. PM wins on MSEL (individual accuracy); GR wins on MSELR
(ranking) and KS (distributional
fidelity).](simulation-study_files/figure-html/loss-barplot-1.png)

Bar chart of loss values across six methods. Lower is better for all
three metrics. PM wins on MSEL (individual accuracy); GR wins on MSELR
(ranking) and KS (distributional fidelity).

### Low-Reliability Contrast

The deep dive above used a 25-item test ($`\bar{w} \approx 0.8`$) where
the data are informative enough for the DPM prior to shine. How do the
patterns change when reliability is much lower?

``` r

comp_hi <- est_comparison
comp_lo <- readRDS(find_extdata("vignette_estimates_comparison_lowrel.rds"))

build_panel <- function(comp_obj, label) {
  N <- length(comp_obj$true_theta)
  data.frame(
    value = c(comp_obj$true_theta,
              comp_obj$dpm$theta$theta_pm,
              comp_obj$dpm$theta$theta_cb,
              comp_obj$dpm$theta$theta_gr),
    Method = factor(rep(c("True", "DPM-PM", "DPM-CB", "DPM-GR"), each = N),
                    levels = c("True", "DPM-PM", "DPM-CB", "DPM-GR")),
    Reliability = label
  )
}

df_rel <- rbind(
  build_panel(comp_lo, "\u03C1 \u2248 0.5  (10 items)"),
  build_panel(comp_hi, "\u03C1 \u2248 0.8  (25 items)")
)
df_rel$Reliability <- factor(df_rel$Reliability,
  levels = c("\u03C1 \u2248 0.5  (10 items)",
             "\u03C1 \u2248 0.8  (25 items)"))

rel_colors <- c("True" = "gray40", "DPM-PM" = pal$dpm_pm,
                "DPM-CB" = pal$dpm_cb, "DPM-GR" = pal$dpm_gr)

ggplot(df_rel, aes(x = value, colour = Method, fill = Method)) +
  geom_density(alpha = 0.12, linewidth = 0.9) +
  facet_wrap(~ Reliability, ncol = 2) +
  scale_colour_manual(values = rel_colors) +
  scale_fill_manual(values = rel_colors) +
  labs(title = "Reliability determines when estimator choice matters",
       x = expression(theta), y = "Density") +
  theme_bw() +
  theme(legend.position = "top",
        strip.text = element_text(face = "bold", size = 11))
```

![Posterior density comparison at two reliability levels (DPM prior
only). Left (\$\rho \approx 0.5\$, 10 items): severe shrinkage causes PM
to collapse the bimodal distribution into a single narrow peak. CB and
GR resist compression, preserving the two-group structure. Right (\$\rho
\approx 0.8\$, 25 items): shrinkage is mild; all three estimators
recover the bimodal shape, with GR tracking the truth most
closely.](simulation-study_files/figure-html/lowrel-comparison-1.png)

Posterior density comparison at two reliability levels (DPM prior only).
Left ($`\rho \approx 0.5`$, 10 items): severe shrinkage causes PM to
collapse the bimodal distribution into a single narrow peak. CB and GR
resist compression, preserving the two-group structure. Right
($`\rho \approx 0.8`$, 25 items): shrinkage is mild; all three
estimators recover the bimodal shape, with GR tracking the truth most
closely.

| Reliability | Method |   MSEL |    MSELR |    KS |
|:------------|:-------|-------:|---------:|------:|
| ρ ≈ 0.5     | PM     | 0.5561 | 0.054310 | 0.285 |
| ρ ≈ 0.5     | CB     | 0.5234 | 0.054310 | 0.220 |
| ρ ≈ 0.5     | GR     | 0.5338 | 0.054209 | 0.165 |
| ρ ≈ 0.8     | PM     | 0.2544 | 0.024992 | 0.095 |
| ρ ≈ 0.8     | CB     | 0.3029 | 0.024992 | 0.080 |
| ρ ≈ 0.8     | GR     | 0.3136 | 0.024868 | 0.055 |

Loss values at two reliability levels. At low reliability the gap
between PM and GR on KS widens substantially. {.table}

> **Key observation.** At $`\bar{w} \approx 0.5`$, the choice of
> estimator (PM vs. CB vs. GR) matters more than the choice of prior,
> because every person’s posterior is heavily shrunk toward the prior
> center. CB and GR resist this compression, producing estimates whose
> distribution more closely tracks the truth — even though their
> individual-level MSE (MSEL) is slightly higher than PM’s.

### Dispersion Recovery at Two Reliability Levels

| Reliability | SD(True) | SD(PM) | SD(CB) | SD(GR) | PM/True | GR/True |
|:------------|---------:|-------:|-------:|-------:|--------:|--------:|
| ρ ≈ 0.5     |    0.991 |  0.493 |  0.700 |  0.696 |    0.50 |    0.70 |
| ρ ≈ 0.8     |    0.955 |  0.948 |  1.075 |  1.073 |    0.99 |    1.12 |

Standard deviation of estimates vs. truth. The PM/True ratio quantifies
shrinkage severity; GR/True shows how well the triple-goal estimator
recovers the true dispersion. {.table}

> **Pattern.** At low reliability, PM compresses the ability
> distribution to roughly half its true spread (PM/True
> $`\approx 0.50`$), effectively destroying the bimodal structure. The
> GR estimator retains about 70% of the true dispersion (GR/True
> $`\approx 0.70`$), substantially reducing the information loss. At
> high reliability, all estimators converge toward the true spread — the
> shrinkage effect becomes negligible.

## Key Patterns

The simulation study from Lee & Wind reveals four main findings:

### Finding 1: Reliability Dominates Over Sample Size

Across all conditions, the marginal reliability $`\bar{w}`$ is a
stronger predictor of estimator performance than sample size $`N`$. At
$`\bar{w} = 0.50`$, all six methods produce similar results because the
data are too weak to distinguish methods. At $`\bar{w} = 0.90`$, method
differences are amplified.

> **Implication:** Invest in longer tests (more items) rather than
> larger samples if your goal is to discriminate between Normal and DPM
> priors.

### Finding 2: DPM Advantage Under Non-Normality + High Reliability

The DPM prior shows its largest advantage when:

- The true latent distribution is non-normal (bimodal or skewed),
  **and**
- Reliability is at least moderate ($`\bar{w} \ge 0.70`$).

Under normality, the DPM prior performs comparably to the Normal prior —
the DP adapts to the true distribution and does not distort it.

### Finding 3: PM Best for Individual MSE; GR Best for KS

The three posterior summary methods optimize different goals:

| Method | Strength | When to use |
|:---|:---|:---|
| **PM** | Lowest individual MSE (MSEL) | Point predictions for each person |
| **CB** | Better distributional match than PM | Group-level summary statistics |
| **GR** | Lowest KS distance | Percentile reports, distributional inference |

PM is the traditional choice and remains best for Goal 1 (minimizing
per-person squared error). However, PM always over-shrinks toward the
prior center, producing estimates whose empirical distribution is too
narrow.

CB and GR “de-shrink” the estimates to better match the true
distribution, at a small cost to individual MSE.

### Finding 4: No One-Size-Fits-All

No single method–prior combination dominates across all conditions and
loss functions simultaneously. The optimal choice depends on:

- The shape of the true latent distribution (or your best guess).
- The reliability of the instrument.
- The inferential goal (individual prediction vs. distributional
  reporting).

> **Recommendation:** For routine use, fit both Normal and DPM models,
> compare via WAIC, and select the estimator (PM, CB, or GR) that
> matches your reporting goals.

## Scaling to a Full Study

The single-condition pipeline above can be wrapped in a loop for a full
Monte Carlo study. Below is a template that you would run in a
non-interactive session (e.g., on a computing cluster). This code is
shown but **not executed** in the vignette:

``` r

# ============================================================
# Full simulation study template
# ============================================================
# WARNING: This runs 18 conditions x 100 replications x 2 models
# Estimated time: ~72 hours on a single core
# ============================================================

library(DPMirt)

# --- Design ---
conditions <- expand.grid(
  latent_shape = c("normal", "bimodal", "skew_pos"),
  n_persons    = c(50, 200),
  target_rho   = c(0.5, 0.7, 0.9),
  stringsAsFactors = FALSE
)

n_reps  <- 100
n_items <- 25

# --- Storage ---
all_results <- vector("list", nrow(conditions) * n_reps)
result_idx  <- 0

for (cond in seq_len(nrow(conditions))) {
  cfg <- conditions[cond, ]
  cat("Condition", cond, "/", nrow(conditions), ":",
      cfg$latent_shape, "N =", cfg$n_persons,
      "rho =", cfg$target_rho, "\n")

  for (rep in seq_len(n_reps)) {
    result_idx <- result_idx + 1
    rep_seed   <- cond * 1000 + rep

    # --- Step 1: Simulate ---
    sim <- dpmirt_simulate(
      n_persons    = cfg$n_persons,
      n_items      = n_items,
      model        = "rasch",
      target_rho   = cfg$target_rho,
      latent_shape = cfg$latent_shape,
      seed         = rep_seed
    )

    # --- Step 2: Fit Normal ---
    fit_n <- dpmirt(
      sim$response,
      model   = "rasch",
      prior   = "normal",
      niter   = 10000,
      nburnin = 2000,
      seed    = rep_seed + 1,
      verbose = FALSE
    )

    # --- Step 3: Fit DPM ---
    fit_d <- dpmirt(
      sim$response,
      model      = "rasch",
      prior      = "dpm",
      mu_K       = 5,
      confidence = "medium",
      niter      = 10000,
      nburnin    = 2000,
      seed       = rep_seed + 2,
      verbose    = FALSE
    )

    # --- Step 4: Estimate ---
    est_n <- dpmirt_estimates(fit_n, methods = c("pm", "cb", "gr"))
    est_d <- dpmirt_estimates(fit_d, methods = c("pm", "cb", "gr"))

    # --- Step 5: Evaluate ---
    loss_n <- dpmirt_loss(est_n, true_theta = sim$theta,
                           true_beta = sim$beta)
    loss_d <- dpmirt_loss(est_d, true_theta = sim$theta,
                           true_beta = sim$beta)

    loss_n$prior <- "Normal"
    loss_d$prior <- "DPM"

    # --- Step 6: Store ---
    all_results[[result_idx]] <- cbind(
      rbind(loss_n, loss_d),
      condition    = cond,
      replication  = rep,
      latent_shape = cfg$latent_shape,
      n_persons    = cfg$n_persons,
      target_rho   = cfg$target_rho,
      waic_normal  = fit_n$waic,
      waic_dpm     = fit_d$waic,
      reliability  = sim$reliability
    )
  }
}

# --- Combine ---
results_df <- do.call(rbind, all_results)

# --- Save ---
saveRDS(results_df, "simulation_results.rds")
```

### Parallelization

For faster execution, parallelize across conditions using
[`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html) or a
cluster scheduler:

``` r

library(parallel)

run_one_condition <- function(cond, n_reps = 100) {
  cfg <- conditions[cond, ]
  # ... same inner loop as above ...
}

# Run on 4 cores (adjust to your machine)
results_list <- mclapply(
  seq_len(nrow(conditions)),
  run_one_condition,
  mc.cores = 4
)
```

### Aggregation and Reporting

Once all replications are complete, aggregate results by condition:

``` r

library(dplyr)

summary_table <- results_df %>%
  filter(parameter == "theta") %>%
  group_by(latent_shape, n_persons, target_rho, prior, method) %>%
  summarise(
    mean_msel  = mean(msel),
    se_msel    = sd(msel) / sqrt(n()),
    mean_mselr = mean(mselr),
    mean_ks    = mean(ks),
    .groups    = "drop"
  )
```

## Extending the Study

The DPMirt simulation framework supports several extensions beyond the
basic Rasch study:

### 2PL Models

``` r

sim_2pl <- dpmirt_simulate(
  n_persons    = 200,
  n_items      = 25,
  model        = "2pl",
  target_rho   = 0.8,
  latent_shape = "bimodal",
  seed         = 42
)

fit_2pl <- dpmirt(
  sim_2pl$response,
  model = "2pl",
  prior = "dpm",
  niter = 15000,
  nburnin = 3000
)
```

### Custom Latent Distributions

With IRTsimrel installed, you have access to 12 distribution shapes:

``` r

# All supported shapes (when IRTsimrel is available)
shapes <- c("normal", "bimodal", "trimodal", "multimodal",
            "skew_pos", "skew_neg", "heavy_tail", "light_tail",
            "uniform", "floor", "ceiling", "custom")

for (shape in shapes) {
  sim <- dpmirt_simulate(200, 25, latent_shape = shape, seed = 42)
  cat(shape, ": KR-20 =", round(sim$reliability, 3), "\n")
}
```

### Custom Loss Functions

Extend evaluation with your own loss function:

``` r

# Bias (signed error)
bias_loss <- function(estimate, true) mean(estimate - true)

# Coverage of 95% CI
coverage_loss <- function(estimate, true) {
  # Requires access to full posterior...
  # Use dpmirt_draws() for this
}

loss_custom <- dpmirt_loss(
  est_dpm,
  true_theta  = sim$theta,
  custom_loss = bias_loss
)
```

## Diagnostic Checks for Simulation Validity

Before interpreting results, verify that the simulation and MCMC
converged properly.

### Check 1: Achieved Reliability

``` r

cat("Target reliability:", 0.8, "\n")
#> Target reliability: 0.8
cat("Achieved KR-20:   ", round(sim$reliability, 4), "\n")
#> Achieved KR-20:    0.8099
cat("Simulation method: ", sim$method, "\n")
#> Simulation method:  irtsimrel
```

### Check 2: MCMC Convergence

For pre-computed fits, you would check:

``` r

# Minimum ESS across all parameters
min(fit$ess$items)   # Should be > 400 for 8000 post-burnin samples
min(fit$ess$theta)

# Trace plot
plot(fit, type = "trace")

# WAIC comparison
cat("WAIC Normal:", fit_normal$waic, "\n")
cat("WAIC DPM:   ", fit_dpm$waic, "\n")
```

### Check 3: Posterior Predictive Check

``` r

# Compare observed vs. predicted sum score distribution
plot(fit, type = "ppc")
```

## Summary Table of Functions

| Function | Purpose | Phase |
|:---|:---|:---|
| [`dpmirt_simulate()`](https://joonho112.github.io/DPMirt/reference/dpmirt_simulate.md) | Generate IRT data with target reliability | Simulation |
| [`dpmirt()`](https://joonho112.github.io/DPMirt/reference/dpmirt.md) | Fit model (one-step) | Fitting |
| [`dpmirt_estimates()`](https://joonho112.github.io/DPMirt/reference/dpmirt_estimates.md) | Compute PM, CB, GR estimates | Estimation |
| [`dpmirt_loss()`](https://joonho112.github.io/DPMirt/reference/dpmirt_loss.md) | Evaluate MSEL, MSELR, KS | Evaluation |
| [`dpmirt_draws()`](https://joonho112.github.io/DPMirt/reference/dpmirt_draws.md) | Extract raw posterior samples | Post-hoc analysis |
| [`dpmirt_dp_density()`](https://joonho112.github.io/DPMirt/reference/dpmirt_dp_density.md) | Compute DP mixture density | DPM diagnostics |

## What’s Next?

| Vignette | Why read it |
|:---|:---|
| *Quick Start* | Simpler walkthrough of the basic pipeline |
| *Prior Elicitation* | Principled $`\alpha`$ choice via DPprior |
| *Under the Hood* | NIMBLE internals, compile-once pattern, custom samplers |

## References

Lee, J. & Wind, S. Targeting toward inferential goals in Bayesian Rasch
models for estimating person-specific latent traits. *OSF Preprint*.
<https://doi.org/10.31219/osf.io/qrw4n>

Lee, J. (2025). Reliability-targeted simulation of item response data:
Solving the inverse design problem. arXiv preprint arXiv:2512.16012.
<https://arxiv.org/abs/2512.16012>

Paganin, S., Paciorek, C. J., Wehrhahn, C., Rodríguez, A., Rabe-Hesketh,
S., & de Valpine, P. (2023). Computational strategies and estimation
performance with Bayesian semiparametric item response theory models.
*Journal of Educational and Behavioral Statistics*, 48(2), 147–188.
<https://doi.org/10.3102/10769986221136105>

Ghosh, M. (1992). Constrained Bayes estimation with applications.
*JASA*, 87(418), 533–540.

Shen, W., & Louis, T. A. (1998). Triple-goal estimates in two-stage
hierarchical models. *JRSS-B*, 60(2), 455–471.
