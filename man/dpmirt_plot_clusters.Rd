% Please edit documentation in R/plot_gg.R
\name{dpmirt_plot_clusters}
\alias{dpmirt_plot_clusters}
\title{Plot Cluster Count Diagnostics}
\usage{
dpmirt_plot_clusters(fit, ...)
}
\arguments{
\item{fit}{A \code{dpmirt_fit} object with \code{prior = "dpm"}.}

\item{...}{Currently unused.}
}
\value{
A \code{ggplot} object.
}
\description{
For DPM models, displays the trace and histogram of the number of active
clusters across MCMC iterations. Stable oscillation is a useful first-pass
visual mixing check. Requires ggplot2.
}
\examples{
\dontrun{
sim <- dpmirt_simulate(200, 20, model = "rasch",
                       latent_shape = "bimodal", seed = 42)
fit <- dpmirt(sim$response, model = "rasch", prior = "dpm",
              niter = 10000, nburnin = 3000, seed = 123)
dpmirt_plot_clusters(fit)
}

}
\seealso{
\code{\link{plot.dpmirt_fit}}, \code{\link{dpmirt_plot_dp_density}}

Other visualization:
\code{\link[=dpmirt_plot_caterpillar]{dpmirt_plot_caterpillar()}},
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
\code{\link[=plot.dpmirt_fit]{plot.dpmirt_fit()}},
\code{\link[=plot.dpmirt_sim]{plot.dpmirt_sim()}}
}
\concept{visualization}
