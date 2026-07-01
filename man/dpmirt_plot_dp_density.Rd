% Please edit documentation in R/plot_gg.R
\name{dpmirt_plot_dp_density}
\alias{dpmirt_plot_dp_density}
\title{Plot DP Mixture Density}
\usage{
dpmirt_plot_dp_density(fit, ...)
}
\arguments{
\item{fit}{A \code{dpmirt_fit} object with \code{prior = "dpm"} and
a computed DP density.}

\item{...}{Currently unused.}
}
\value{
A \code{ggplot} object.
}
\description{
Displays the posterior mean density of the Dirichlet Process mixture
with pointwise credible bands. Requires a pre-computed
\code{dpmirt_dp_density} object (computed automatically by
\code{\link{dpmirt}} or manually via \code{\link{dpmirt_dp_density}}).
Requires ggplot2.
}
\details{
For 2PL/3PL IRT and SI parameterizations, transformed-scale DP
density reconstruction currently uses the same location-shift contract as
\code{\link{dpmirt_dp_density}}. Treat these plots as diagnostic summaries;
Rasch/location-shift settings are the most directly interpretable.
}
\examples{
\dontrun{
sim <- dpmirt_simulate(200, 20, model = "rasch",
                       latent_shape = "bimodal", seed = 42)
fit <- dpmirt(sim$response, model = "rasch", prior = "dpm",
              niter = 10000, nburnin = 3000, seed = 123)
dpmirt_plot_dp_density(fit)
}

}
\seealso{
\code{\link{dpmirt_dp_density}}, \code{\link{plot.dpmirt_fit}}

Other visualization:
\code{\link[=dpmirt_plot_caterpillar]{dpmirt_plot_caterpillar()}},
\code{\link[=dpmirt_plot_clusters]{dpmirt_plot_clusters()}},
\code{\link[=dpmirt_plot_density]{dpmirt_plot_density()}},
\code{\link[=dpmirt_plot_density_compare]{dpmirt_plot_density_compare()}},
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
