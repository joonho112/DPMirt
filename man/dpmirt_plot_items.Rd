% Please edit documentation in R/plot_gg.R
\name{dpmirt_plot_items}
\alias{dpmirt_plot_items}
\title{Plot Item Parameter Estimates}
\usage{
dpmirt_plot_items(fit, ...)
}
\arguments{
\item{fit}{A \code{dpmirt_fit} object from \code{\link{dpmirt}}.}

\item{...}{Currently unused.}
}
\value{
A \code{ggplot} object.
}
\description{
Displays item difficulty (beta) estimates with posterior credible intervals.
For 2PL/3PL models, also shows discrimination (lambda) estimates.
Requires ggplot2.
}
\examples{
\dontrun{
sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
fit <- dpmirt(sim$response, model = "rasch", prior = "normal",
              niter = 5000, nburnin = 1000, seed = 123)
dpmirt_plot_items(fit)
}

}
\seealso{
\code{\link{plot.dpmirt_fit}}, \code{\link{dpmirt_plot_caterpillar}}

Other visualization:
\code{\link[=dpmirt_plot_caterpillar]{dpmirt_plot_caterpillar()}},
\code{\link[=dpmirt_plot_clusters]{dpmirt_plot_clusters()}},
\code{\link[=dpmirt_plot_density]{dpmirt_plot_density()}},
\code{\link[=dpmirt_plot_density_compare]{dpmirt_plot_density_compare()}},
\code{\link[=dpmirt_plot_dp_density]{dpmirt_plot_dp_density()}},
\code{\link[=dpmirt_plot_icc]{dpmirt_plot_icc()}},
\code{\link[=dpmirt_plot_info]{dpmirt_plot_info()}},
\code{\link[=dpmirt_plot_parameter_trace]{dpmirt_plot_parameter_trace()}},
\code{\link[=dpmirt_plot_pp_check]{dpmirt_plot_pp_check()}},
\code{\link[=dpmirt_plot_trace]{dpmirt_plot_trace()}},
\code{\link[=dpmirt_plot_wright_map]{dpmirt_plot_wright_map()}},
\code{\link[=plot.dpmirt_estimates]{plot.dpmirt_estimates()}},
\code{\link[=plot.dpmirt_fit]{plot.dpmirt_fit()}},
\code{\link[=plot.dpmirt_sim]{plot.dpmirt_sim()}}
}
\concept{visualization}
