# Package index

## Package

- [`DPMirt-package`](https://joonho112.github.io/DPMirt/reference/DPMirt-package.md)
  [`DPMirt`](https://joonho112.github.io/DPMirt/reference/DPMirt-package.md)
  : DPMirt: Bayesian Semiparametric IRT Models Using DPM Priors

## Main Fitting

All-in-one model fitting function

- [`dpmirt()`](https://joonho112.github.io/DPMirt/reference/dpmirt.md) :
  Fit a Bayesian IRT model with DPMirt

## Step-by-Step Workflow

Modular pipeline for fine-grained control

- [`dpmirt_spec()`](https://joonho112.github.io/DPMirt/reference/dpmirt_spec.md)
  [`print(`*`<dpmirt_spec>`*`)`](https://joonho112.github.io/DPMirt/reference/dpmirt_spec.md)
  : Create a DPMirt model specification
- [`dpmirt_compile()`](https://joonho112.github.io/DPMirt/reference/dpmirt_compile.md)
  [`print(`*`<dpmirt_compiled>`*`)`](https://joonho112.github.io/DPMirt/reference/dpmirt_compile.md)
  : Compile a DPMirt model specification
- [`dpmirt_sample()`](https://joonho112.github.io/DPMirt/reference/dpmirt_sample.md)
  [`print(`*`<dpmirt_samples>`*`)`](https://joonho112.github.io/DPMirt/reference/dpmirt_sample.md)
  : Run MCMC sampling on a compiled DPMirt model
- [`dpmirt_rescale()`](https://joonho112.github.io/DPMirt/reference/dpmirt_rescale.md)
  : Rescale posterior samples for identification
- [`dpmirt_resume()`](https://joonho112.github.io/DPMirt/reference/dpmirt_resume.md)
  : Resume MCMC sampling from a previous run

## Posterior Estimation

Triple-goal estimation (PM, CB, GR) and raw posterior access

- [`dpmirt_estimates()`](https://joonho112.github.io/DPMirt/reference/dpmirt_estimates.md)
  [`print(`*`<dpmirt_estimates>`*`)`](https://joonho112.github.io/DPMirt/reference/dpmirt_estimates.md)
  : Compute posterior estimates using PM, CB, and GR methods
- [`dpmirt_draws()`](https://joonho112.github.io/DPMirt/reference/dpmirt_draws.md)
  : Extract posterior draws from a DPMirt fit

## Model Assessment

Diagnostics, comparison, and loss evaluation

- [`dpmirt_diagnostics()`](https://joonho112.github.io/DPMirt/reference/dpmirt_diagnostics.md)
  [`print(`*`<dpmirt_diagnostics>`*`)`](https://joonho112.github.io/DPMirt/reference/dpmirt_diagnostics.md)
  : Compute MCMC Diagnostics for a DPMirt Fit
- [`dpmirt_compare()`](https://joonho112.github.io/DPMirt/reference/dpmirt_compare.md)
  : Compare DPMirt models using information criteria
- [`dpmirt_loss()`](https://joonho112.github.io/DPMirt/reference/dpmirt_loss.md)
  : Evaluate estimator performance using loss functions

## DPM-Specific Analysis

Dirichlet Process Mixture density and prior elicitation

- [`dpmirt_dp_density()`](https://joonho112.github.io/DPMirt/reference/dpmirt_dp_density.md)
  [`print(`*`<dpmirt_dp_density>`*`)`](https://joonho112.github.io/DPMirt/reference/dpmirt_dp_density.md)
  : Compute posterior density of the DP mixture
- [`dpmirt_alpha_prior()`](https://joonho112.github.io/DPMirt/reference/dpmirt_alpha_prior.md)
  : Elicit a Gamma Hyperprior for the DP Concentration Parameter

## Simulation

Data generation with reliability targeting

- [`dpmirt_simulate()`](https://joonho112.github.io/DPMirt/reference/dpmirt_simulate.md)
  [`print(`*`<dpmirt_sim>`*`)`](https://joonho112.github.io/DPMirt/reference/dpmirt_simulate.md)
  : Simulate IRT response data

## Visualization

Plotting functions for model results (ggplot2-based)

- [`dpmirt_plot_density()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_density.md)
  : Plot Posterior Mean Density
- [`dpmirt_plot_items()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_items.md)
  : Plot Item Parameter Estimates
- [`dpmirt_plot_trace()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_trace.md)
  : Plot Log-Likelihood MCMC Trace
- [`dpmirt_plot_parameter_trace()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_parameter_trace.md)
  : Plot Individual Parameter MCMC Traces
- [`dpmirt_plot_caterpillar()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_caterpillar.md)
  : Plot Caterpillar (Forest) Plot
- [`dpmirt_plot_icc()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_icc.md)
  : Plot Item Characteristic Curves
- [`dpmirt_plot_info()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_info.md)
  : Plot Test Information Function
- [`dpmirt_plot_wright_map()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_wright_map.md)
  : Plot Person-Item Map (Wright Map)
- [`dpmirt_plot_clusters()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_clusters.md)
  : Plot Cluster Count Diagnostics
- [`dpmirt_plot_dp_density()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_dp_density.md)
  : Plot DP Mixture Density
- [`dpmirt_plot_density_compare()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_density_compare.md)
  : Plot Posterior Density vs Reference Distribution
- [`dpmirt_plot_pp_check()`](https://joonho112.github.io/DPMirt/reference/dpmirt_plot_pp_check.md)
  : Plot Posterior Predictive Check

## S3 Methods

Print, summary, coef, and plot methods for DPMirt objects

- [`print(`*`<dpmirt_fit>`*`)`](https://joonho112.github.io/DPMirt/reference/dpmirt_fit-methods.md)
  [`summary(`*`<dpmirt_fit>`*`)`](https://joonho112.github.io/DPMirt/reference/dpmirt_fit-methods.md)
  [`coef(`*`<dpmirt_fit>`*`)`](https://joonho112.github.io/DPMirt/reference/dpmirt_fit-methods.md)
  : Methods for dpmirt_fit Objects
- [`plot(`*`<dpmirt_fit>`*`)`](https://joonho112.github.io/DPMirt/reference/plot.dpmirt_fit.md)
  : Plot a DPMirt fit object
- [`plot(`*`<dpmirt_estimates>`*`)`](https://joonho112.github.io/DPMirt/reference/plot.dpmirt_estimates.md)
  : Plot dpmirt_estimates object
- [`plot(`*`<dpmirt_sim>`*`)`](https://joonho112.github.io/DPMirt/reference/plot.dpmirt_sim.md)
  : Plot dpmirt_sim object
