% Please edit documentation in R/plot.R
\name{plot.dpmirt_estimates}
\alias{plot.dpmirt_estimates}
\title{Plot dpmirt_estimates object}
\usage{
\method{plot}{dpmirt_estimates}(x, type = c("estimates", "shrinkage"), param = c("theta", "beta"), ...)
}
\arguments{
\item{x}{A \code{dpmirt_estimates} object.}

\item{type}{Character. Plot type: \code{"estimates"} (caterpillar of
PM/CB/GR) or \code{"shrinkage"} (PM vs CB scatter).}

\item{param}{Character. Which parameter: \code{"theta"} or \code{"beta"}.}

\item{...}{Additional graphical parameters.}
}
\value{
Invisibly returns the input \code{dpmirt_estimates} object.
}
\description{
Plot dpmirt_estimates object
}
\examples{
\dontrun{
sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
fit <- dpmirt(sim$response, model = "rasch", prior = "normal",
              niter = 5000, nburnin = 1000, seed = 123)
est <- dpmirt_estimates(fit)

# PM/CB/GR caterpillar plot
plot(est, type = "estimates")

# Shrinkage plot (PM vs CB)
plot(est, type = "shrinkage")
}

}
\seealso{
\code{\link{dpmirt_estimates}}, \code{\link{dpmirt_plot_caterpillar}}

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
\code{\link[=plot.dpmirt_fit]{plot.dpmirt_fit()}},
\code{\link[=plot.dpmirt_sim]{plot.dpmirt_sim()}}
}
\concept{visualization}
