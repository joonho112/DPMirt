% Please edit documentation in R/plot.R
\name{plot.dpmirt_sim}
\alias{plot.dpmirt_sim}
\title{Plot dpmirt_sim object}
\usage{
\method{plot}{dpmirt_sim}(x, type = c("parameters", "response"), ...)
}
\arguments{
\item{x}{A \code{dpmirt_sim} object.}

\item{type}{Character. Plot type: \code{"parameters"} (histograms of true
parameters) or \code{"response"} (heatmap of response matrix).}

\item{...}{Additional graphical parameters.}
}
\value{
Invisibly returns the input \code{dpmirt_sim} object.
}
\description{
Plot dpmirt_sim object
}
\examples{
\dontrun{
sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)

# True parameter histograms
plot(sim, type = "parameters")

# Response matrix heatmap
plot(sim, type = "response")
}

}
\seealso{
\code{\link{dpmirt_simulate}}

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
\code{\link[=plot.dpmirt_fit]{plot.dpmirt_fit()}}
}
\concept{visualization}
