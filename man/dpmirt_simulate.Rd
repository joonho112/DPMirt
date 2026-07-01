% Please edit documentation in R/simulate.R
\name{dpmirt_simulate}
\alias{dpmirt_simulate}
\alias{print.dpmirt_sim}
\title{Simulate IRT response data}
\usage{
dpmirt_simulate(
  n_persons,
  n_items,
  model = c("rasch", "2pl", "3pl"),
  target_rho = 0.8,
  latent_shape = "normal",
  item_source = "parametric",
  reliability_metric = c("info", "msem"),
  latent_params = list(),
  item_params = list(),
  M = 10000L,
  seed = NULL,
  use_irtsimrel = TRUE,
  verbose = FALSE,
  ...
)

\method{print}{dpmirt_sim}(x, ...)
}
\arguments{
\item{n_persons}{Integer. Number of persons.}

\item{n_items}{Integer. Number of items.}

\item{model}{Character. IRT model type: \code{"rasch"}, \code{"2pl"},
or \code{"3pl"}. IRTsimrel reliability targeting is available for Rasch
and 2PL; 3PL simulations use DPMirt's fallback generator.}

\item{target_rho}{Numeric in (0, 1). Target marginal reliability. Only used
when IRTsimrel is available and \code{model} is \code{"rasch"} or
\code{"2pl"}. Default 0.8.}

\item{latent_shape}{Character. Shape of latent ability distribution.
When using IRTsimrel, supports all 12 shapes: \code{"normal"},
\code{"bimodal"}, \code{"trimodal"}, \code{"multimodal"},
\code{"skew_pos"}, \code{"skew_neg"}, \code{"heavy_tail"},
\code{"light_tail"}, \code{"uniform"}, \code{"floor"},
\code{"ceiling"}, \code{"custom"}.
For \code{"custom"}, supply
\code{latent_params = list(mixture_spec = list(weights = ..., means = ..., sds = ...))};
custom shapes require IRTsimrel with \code{model = "rasch"} or
\code{model = "2pl"}.
Fallback supports: \code{"normal"}, \code{"bimodal"}, \code{"skewed"}.}

\item{item_source}{Character. Source for item parameters in IRTsimrel.
One of \code{"parametric"} (default), \code{"irw"},
\code{"hierarchical"}, \code{"custom"}.}

\item{reliability_metric}{Character. Reliability metric for IRTsimrel
calibration. \code{"info"} (default, average-information, recommended)
uses EQC calibration; \code{"msem"} uses IRTsimrel's SAC calibration.}

\item{latent_params}{List. Additional parameters passed to
\code{IRTsimrel::sim_latentG()} (e.g.,
\code{list(shape_params = list(delta = 0.8))}). If
\code{latent_shape = "custom"}, this must include
\code{mixture_spec}.}

\item{item_params}{List. Additional parameters passed to
\code{IRTsimrel::sim_item_params()} (e.g.,
\code{list(discrimination_params = list(rho = -0.3))}).}

\item{M}{Integer. Quadrature sample size for EQC, and pre-calibration
quadrature size for SAC. Default 10000.}

\item{seed}{Integer or NULL. Random seed for reproducibility.}

\item{use_irtsimrel}{Logical. If TRUE (default), attempt to use
IRTsimrel for reliability-targeted simulation. Falls back to internal
simulation if IRTsimrel is not installed.}

\item{verbose}{Logical. Print progress messages. Default FALSE.}

\item{...}{Additional arguments (currently unused).}

\item{x}{A \code{dpmirt_sim} object.}
}
\value{
A \code{dpmirt_sim} S3 object containing:
\describe{
\item{response}{N x I binary response matrix}
\item{theta}{True person abilities (length N)}
\item{beta}{True item difficulties (length I)}
\item{lambda}{True discriminations (length I for 2PL/3PL, NULL for Rasch)}
\item{delta}{True guessing parameters (length I for 3PL, NULL otherwise)}
\item{n_persons, n_items, model}{Simulation settings}
\item{reliability}{Empirical KR-20 reliability of the response matrix}
\item{achieved_rho}{IRTsimrel calibration achieved marginal reliability
(NULL if fallback)}
\item{target_rho}{Requested target reliability (NULL if fallback)}
\item{latent_shape}{Distribution shape used}
\item{eqc_result}{IRTsimrel calibration result (NULL if fallback).
The field name is retained for backward compatibility; the object class
distinguishes EQC and SAC results.}
\item{method}{Character: "irtsimrel" or "fallback"}
}
}
\description{
Generates binary response data under known IRT parameters.
When IRTsimrel is available, uses IRTsimrel calibration to target a
marginal reliability. The returned \code{reliability} field is the
empirical KR-20 reliability of the generated response matrix; the
calibration-level achieved marginal reliability is returned as
\code{achieved_rho}. Otherwise falls back to Paganin-style simulation
without reliability targeting.
}
\examples{
\dontrun{
# Simple fallback simulation
sim <- dpmirt_simulate(200, 25, model = "rasch", seed = 42)

# With IRTsimrel (reliability-targeted)
sim <- dpmirt_simulate(200, 25, model = "rasch",
                       target_rho = 0.85,
                       latent_shape = "bimodal",
                       seed = 42)

# Full simulation study workflow
sim <- dpmirt_simulate(200, 25, model = "rasch",
                       target_rho = 0.8,
                       latent_shape = "bimodal")
fit <- dpmirt(sim$response, model = "rasch", prior = "dpm")
est <- dpmirt_estimates(fit)
dpmirt_loss(est, true_theta = sim$theta, true_beta = sim$beta)
}

}
\seealso{
\code{\link{dpmirt}}, \code{\link{dpmirt_loss}},
\code{\link{plot.dpmirt_sim}}

Other simulation:
\code{\link[=dpmirt_loss]{dpmirt_loss()}}
}
\concept{simulation}
