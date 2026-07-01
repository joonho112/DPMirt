% Please edit documentation in R/plot_gg.R
\name{dpmirt_plot_caterpillar}
\alias{dpmirt_plot_caterpillar}
\title{Plot Caterpillar (Forest) Plot}
\usage{
dpmirt_plot_caterpillar(
  fit,
  param = c("beta", "theta", "lambda", "delta"),
  sort = TRUE,
  ci_level = 0.95,
  max_show = 50,
  ...
)
}
\arguments{
\item{fit}{A \code{dpmirt_fit} object from \code{\link{dpmirt}}.}

\item{param}{Character. Parameter type: \code{"beta"} (default),
\code{"theta"}, \code{"lambda"}, or \code{"delta"}.}

\item{sort}{Logical. Sort by estimate magnitude. Default TRUE.}

\item{ci_level}{Numeric. Credible interval level. Default 0.95.}

\item{max_show}{Integer. Maximum number of parameters to display. Default 50.}

\item{...}{Currently unused.}
}
\value{
A \code{ggplot} object.
}
\description{
Displays sorted point estimates with credible intervals for items or
persons. Useful for identifying extreme values and assessing uncertainty.
Requires ggplot2.
}
\examples{
\dontrun{
sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
fit <- dpmirt(sim$response, model = "rasch", prior = "normal",
              niter = 5000, nburnin = 1000, seed = 123)
dpmirt_plot_caterpillar(fit, param = "beta")
}

}
\seealso{
\code{\link{plot.dpmirt_fit}}, \code{\link{dpmirt_plot_items}}

Other visualization:
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
\code{\link[=plot.dpmirt_fit]{plot.dpmirt_fit()}},
\code{\link[=plot.dpmirt_sim]{plot.dpmirt_sim()}}
}
\concept{visualization}
