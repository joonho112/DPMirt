# Articles

### All vignettes

- [DPMirt: Bayesian Semiparametric Item Response
  Theory](https://joonho112.github.io/DPMirt/articles/introduction.md):

  A gateway introduction to the DPMirt package for fitting Bayesian
  semiparametric IRT models using Dirichlet Process Mixture priors, with
  emphasis on the three inferential goals that motivate posterior
  summary selection and the complete model family supported by the
  package.

- [The Complete Guide to Models and
  Workflows](https://joonho112.github.io/DPMirt/articles/models-and-workflow.md):

  A comprehensive guide to every model, prior, parameterization, and
  identification strategy available in DPMirt. Covers both the one-step
  dpmirt() wrapper and the step-by-step pipeline, with practical advice
  on when to use each workflow.

- [Under the Hood: NIMBLE Backend and Advanced
  Usage](https://joonho112.github.io/DPMirt/articles/nimble-internals.md):

  A deep dive into DPMirt’s NIMBLE implementation — programmatic code
  generation, custom samplers and distributions, the compile-once
  sample-many architecture, post-hoc rescaling, and DP density
  reconstruction — for users who want to understand or extend the
  package internals.

- [Posterior Summary Methods: Matching Estimators to Inferential
  Goals](https://joonho112.github.io/DPMirt/articles/posterior-summaries.md):

  How to choose among posterior mean (PM), constrained Bayes (CB), and
  triple-goal (GR) estimators depending on whether the goal is
  individual scoring, ranking, or distribution recovery. Tied closely to
  the APM manuscript (Lee & Wind).

- [Principled Prior Elicitation for the DP Concentration
  Parameter](https://joonho112.github.io/DPMirt/articles/prior-elicitation.md):

  Learn why the DP concentration parameter alpha matters, how the
  Paganin default compares to principled elicitation via DPprior, and
  how to integrate both approaches into your DPMirt workflow — including
  sensitivity analysis and custom base measure tuning.

- [Quick Start: Your First IRT Model in 5
  Minutes](https://joonho112.github.io/DPMirt/articles/quick-start.md):

  A hands-on introduction to DPMirt: simulate data, fit Rasch models
  with Normal and DPM priors, visualize results, and extract person
  ability estimates — all in under five minutes of reading.

- [Simulation Study: Evaluating Prior Models and Posterior
  Summaries](https://joonho112.github.io/DPMirt/articles/simulation-study.md):

  A comprehensive guide to designing, executing, and analyzing
  simulation studies with DPMirt — covering factorial design, IRTsimrel
  integration, the full analysis pipeline, loss function evaluation, and
  key findings from the APM manuscript.

- [Mathematical Foundations: IRT Models and Dirichlet Process Mixture
  Priors](https://joonho112.github.io/DPMirt/articles/theory-irt-dpm.md):

  A comprehensive treatment of the mathematical foundations underlying
  the DPMirt package, covering Item Response Theory models (Rasch, 2PL,
  3PL), identification and post-hoc rescaling, Dirichlet Process Mixture
  priors, the concentration parameter, and posterior summary theory (PM,
  CB, GR).
