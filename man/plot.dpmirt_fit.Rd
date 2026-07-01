% Please edit documentation in R/plot.R
\name{plot.dpmirt_fit}
\alias{plot.dpmirt_fit}
\title{Plot a DPMirt fit object}
\usage{
\method{plot}{dpmirt_fit}(
  x,
  type = c("density", "items", "trace", "clusters", "dp_density", "icc", "wright_map",
    "parameter_trace", "caterpillar", "density_compare", "info", "pp_check"),
  engine = c("auto", "base", "ggplot2"),
  ...
)
}
\arguments{
\item{x}{A \code{dpmirt_fit} object.}

\item{type}{Character. Plot type. One of:
\describe{
\item{\code{"density"}}{Kernel density of posterior mean theta.}
\item{\code{"items"}}{Item difficulty estimates with error bars.}
\item{\code{"trace"}}{Log-likelihood MCMC trace.}
\item{\code{"clusters"}}{Cluster count trace and histogram (DPM only).}
\item{\code{"dp_density"}}{DP mixture density with credible band
(DPM only; see \code{\link{dpmirt_dp_density}} for the
transformed-scale caveat).}
\item{\code{"icc"}}{Item Characteristic Curves.}
\item{\code{"wright_map"}}{Person-Item map (Wright map).}
\item{\code{"parameter_trace"}}{Individual parameter MCMC traces.}
\item{\code{"caterpillar"}}{Sorted estimates with credible intervals.}
\item{\code{"density_compare"}}{Posterior density vs reference overlay.}
\item{\code{"info"}}{Test Information Function.}
\item{\code{"pp_check"}}{Posterior predictive check.}
}}

\item{engine}{Character. Plotting backend: \code{"auto"} (default) uses
ggplot2 if available, \code{"base"} forces base R, \code{"ggplot2"}
requires ggplot2.}

\item{...}{Additional arguments passed to the specific plotting function.}
}
\value{
Invisibly returns a ggplot object when the ggplot2 backend is used;
invisibly returns \code{NULL} when the base R backend is used.
}
\description{
Produces visualizations for fitted IRT models. Supports 12 plot types
with automatic selection between base R and ggplot2 backends.
}
\details{
For \code{type = "dp_density"}, Rasch/location-shift settings are
the most directly interpretable. For transformed-scale 2PL/3PL IRT and SI
fits, use the plot as a diagnostic summary because full scale/Jacobian
density adjustment is not yet implemented.
}
\examples{
\dontrun{
sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
fit <- dpmirt(sim$response, model = "rasch", prior = "normal",
              niter = 5000, nburnin = 1000, seed = 123)

# Theta density
plot(fit, type = "density")

# Item difficulty estimates
plot(fit, type = "items")

# MCMC trace
plot(fit, type = "trace")

# Force base R backend
plot(fit, type = "density", engine = "base")
}

}
\seealso{
\code{\link{dpmirt_plot_density}}, \code{\link{dpmirt_plot_items}},
\code{\link{dpmirt_plot_trace}}, \code{\link{dpmirt_plot_icc}}

Other visualization:
\code{\link[=dpmirt_plot_caterpillar]{dpmirt_plot_caterpillar()}},
\code{\link[=dpmirt_plot_clusters]{dpmirt_plot_clusters()}},
\code{\link[=dpmirt_plot_density]{dpmirt_plot_density()}},
\code{\link[=dpmirt_plot_density_compare]{dpmirt_plot_density_compare()}},
\code{\link[=dpmirt_plot_dp_density]{dpmirt_plot_dp_density()}},
\code{\link[=dpmirt_plot_icc]{dpmirt_plot_icc()}},
\code{\link[=dpmirt_plot_info]{dpmirt_plot_info()}},
\code{\link[=dpmirt_plot_items]{dpmirt_plot_items()}},
\code{\link[=dpmirt_plot_parameter_trace]{dpmirt_plot_parameter_trace()}},
\code{\link[=dpmirt_plot_pp_check]{dpmirt_plot_pp_check()}},
\code{\link[=dpmirt_plot_trace]{dpmirt_plot_trace()}},
\code{\link[=dpmirt_plot_wright_map]{dpmirt_plot_wright_map()}},
\code{\link[=plot.dpmirt_estimates]{plot.dpmirt_estimates()}},
\code{\link[=plot.dpmirt_sim]{plot.dpmirt_sim()}}
}
\concept{visualization}
