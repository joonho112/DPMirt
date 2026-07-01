% Please edit documentation in R/plot_gg.R
\name{dpmirt_plot_wright_map}
\alias{dpmirt_plot_wright_map}
\title{Plot Person-Item Map (Wright Map)}
\usage{
dpmirt_plot_wright_map(fit, bins = 30, item_labels = NULL, ...)
}
\arguments{
\item{fit}{A \code{dpmirt_fit} object from \code{\link{dpmirt}}.}

\item{bins}{Integer. Number of histogram bins. Default: 30.}

\item{item_labels}{Optional character vector of item labels.}

\item{...}{Currently unused.}
}
\value{
A \code{ggplot} object.
}
\description{
Displays a Wright map showing the joint distribution of person abilities
and item difficulties on a common logit scale. Useful for identifying
gaps in item coverage. Requires ggplot2.
}
\examples{
\dontrun{
sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
fit <- dpmirt(sim$response, model = "rasch", prior = "normal",
              niter = 5000, nburnin = 1000, seed = 123)
dpmirt_plot_wright_map(fit)
}

}
\seealso{
\code{\link{plot.dpmirt_fit}}, \code{\link{dpmirt_plot_items}}

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
\code{\link[=plot.dpmirt_estimates]{plot.dpmirt_estimates()}},
\code{\link[=plot.dpmirt_fit]{plot.dpmirt_fit()}},
\code{\link[=plot.dpmirt_sim]{plot.dpmirt_sim()}}
}
\concept{visualization}
