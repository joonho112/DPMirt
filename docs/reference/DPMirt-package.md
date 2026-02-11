# DPMirt: Bayesian Semiparametric IRT Models Using DPM Priors

Fits Bayesian semiparametric Item Response Theory (IRT) models using
Dirichlet Process Mixture (DPM) priors via NIMBLE. Supports Rasch, 2PL,
and 3PL models with parametric (Normal) or semiparametric (Dirichlet
Process Mixture) priors on the latent ability distribution, and provides
triple-goal posterior summaries (PM, CB, GR) for simultaneous estimation
and ranking of person abilities.

## Details

The DPMirt package supports:

- Three IRT models: Rasch, 2PL, and 3PL

- Two latent trait priors: parametric (Normal) and semiparametric (DPM)

- Three identification strategies: constrained_item,
  constrained_ability, unconstrained

- Three posterior summary methods: PM, CB (Ghosh 1992), GR (Shen & Louis
  1998)

- Compile-once, sample-many MCMC workflow

- Principled alpha hyperprior elicitation via DPprior

- Reliability-targeted simulation via IRTsimrel

The backbone NIMBLE model code is adapted from Paganin et al. (2023).
MCMC compilation, sampling, and model management are handled by NIMBLE's
C++ infrastructure.

## Typical Workflow

A standard DPMirt analysis proceeds in five steps:

1.  **Simulate or load data**: Use
    [`dpmirt_simulate`](https://joonho112.github.io/DPMirt/reference/dpmirt_simulate.md)
    or provide a binary response matrix.

2.  **Fit the model**: Call
    [`dpmirt`](https://joonho112.github.io/DPMirt/reference/dpmirt.md)
    for one-step fitting, or use the step-by-step pipeline
    [`dpmirt_spec`](https://joonho112.github.io/DPMirt/reference/dpmirt_spec.md)
    \\\rightarrow\\
    [`dpmirt_compile`](https://joonho112.github.io/DPMirt/reference/dpmirt_compile.md)
    \\\rightarrow\\
    [`dpmirt_sample`](https://joonho112.github.io/DPMirt/reference/dpmirt_sample.md)
    \\\rightarrow\\
    [`dpmirt_rescale`](https://joonho112.github.io/DPMirt/reference/dpmirt_rescale.md).

3.  **Summarize**: Compute triple-goal estimates with
    [`dpmirt_estimates`](https://joonho112.github.io/DPMirt/reference/dpmirt_estimates.md)
    and extract posterior draws with
    [`dpmirt_draws`](https://joonho112.github.io/DPMirt/reference/dpmirt_draws.md).

4.  **Diagnose**: Evaluate convergence and model comparison via
    [`dpmirt_diagnostics`](https://joonho112.github.io/DPMirt/reference/dpmirt_diagnostics.md)
    and
    [`dpmirt_compare`](https://joonho112.github.io/DPMirt/reference/dpmirt_compare.md).

5.  **Visualize**: Use
    [`plot`](https://joonho112.github.io/DPMirt/reference/plot.dpmirt_fit.md)
    and the `dpmirt_plot_*` family for publication-quality figures.

## Model Fitting

- [`dpmirt`](https://joonho112.github.io/DPMirt/reference/dpmirt.md):

  One-step model fitting (specification + compilation + sampling +
  rescaling)

- [`dpmirt_spec`](https://joonho112.github.io/DPMirt/reference/dpmirt_spec.md):

  Create model specification

- [`dpmirt_compile`](https://joonho112.github.io/DPMirt/reference/dpmirt_compile.md):

  Compile NIMBLE model and MCMC

- [`dpmirt_sample`](https://joonho112.github.io/DPMirt/reference/dpmirt_sample.md):

  Run MCMC sampling

- [`dpmirt_resume`](https://joonho112.github.io/DPMirt/reference/dpmirt_resume.md):

  Continue sampling from a fitted model

## Estimation and Rescaling

- [`dpmirt_rescale`](https://joonho112.github.io/DPMirt/reference/dpmirt_rescale.md):

  Post-hoc identification rescaling (Rasch, IRT, SI)

- [`dpmirt_estimates`](https://joonho112.github.io/DPMirt/reference/dpmirt_estimates.md):

  Compute PM, CB, and GR triple-goal posterior summaries

- [`dpmirt_draws`](https://joonho112.github.io/DPMirt/reference/dpmirt_draws.md):

  Extract posterior draws as matrix or long-format data frame

## Diagnostics and Model Comparison

- [`dpmirt_diagnostics`](https://joonho112.github.io/DPMirt/reference/dpmirt_diagnostics.md):

  MCMC convergence diagnostics (ESS, Rhat, trace summaries)

- [`dpmirt_compare`](https://joonho112.github.io/DPMirt/reference/dpmirt_compare.md):

  WAIC-based model comparison

## Simulation

- [`dpmirt_simulate`](https://joonho112.github.io/DPMirt/reference/dpmirt_simulate.md):

  Simulate IRT data with flexible latent distributions and reliability
  targeting

- [`dpmirt_loss`](https://joonho112.github.io/DPMirt/reference/dpmirt_loss.md):

  Evaluate estimator loss (MSEL, MSELR, KS, custom)

## Prior Specification

- [`dpmirt_alpha_prior`](https://joonho112.github.io/DPMirt/reference/dpmirt_alpha_prior.md):

  Principled DPM concentration parameter elicitation

## DP Density Estimation

- [`dpmirt_dp_density`](https://joonho112.github.io/DPMirt/reference/dpmirt_dp_density.md):

  Posterior density estimation from the Dirichlet Process mixture

## Visualization

S3 plot methods:

- [`plot(fit)`](https://joonho112.github.io/DPMirt/reference/plot.dpmirt_fit.md):

  Trace plots, density plots, and caterpillar plots for fitted models

- [`plot(estimates)`](https://joonho112.github.io/DPMirt/reference/plot.dpmirt_estimates.md):

  Shrinkage, rank, and comparison plots for PM/CB/GR estimates

- [`plot(sim)`](https://joonho112.github.io/DPMirt/reference/plot.dpmirt_sim.md):

  ICC and distribution plots for simulated data

Standalone ggplot2 functions (require ggplot2):
[`dpmirt_plot_trace`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_trace.md),
[`dpmirt_plot_density`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_density.md),
[`dpmirt_plot_caterpillar`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_caterpillar.md),
[`dpmirt_plot_items`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_items.md),
[`dpmirt_plot_icc`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_icc.md),
[`dpmirt_plot_info`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_info.md),
[`dpmirt_plot_dp_density`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_dp_density.md),
[`dpmirt_plot_clusters`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_clusters.md),
[`dpmirt_plot_wright_map`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_wright_map.md),
[`dpmirt_plot_parameter_trace`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_parameter_trace.md),
[`dpmirt_plot_density_compare`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_density_compare.md),
[`dpmirt_plot_pp_check`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_pp_check.md)

## References

Paganin, S., Paciorek, C. J., Wehrhahn, C., Rodriguez, A., Rabe-Hesketh,
S., & de Valpine, P. (2023). Computational strategies and estimation
performance with Bayesian semiparametric item response theory models.
*Journal of Educational and Behavioral Statistics, 48*(2), 147–188.

Ghosh, M. (1992). Constrained Bayes estimation with applications.
*Journal of the American Statistical Association, 87*(418), 533–540.

Shen, W., & Louis, T. A. (1998). Triple-goal estimates in two-stage
hierarchical models. *Journal of the Royal Statistical Society: Series
B, 60*(2), 455–471.

## See also

Useful links:

- <https://github.com/joonho112/DPMirt>

- Report bugs at <https://github.com/joonho112/DPMirt/issues>

## Author

**Maintainer**: JoonHo Lee <jlee296@ua.edu>
