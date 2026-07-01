% Please edit documentation in R/draws.R
\name{dpmirt_draws}
\alias{dpmirt_draws}
\title{Extract posterior draws from a DPMirt fit}
\usage{
dpmirt_draws(
  fit,
  vars = c("theta", "beta", "lambda", "delta"),
  format = c("matrix", "long"),
  use_rescaled = TRUE
)
}
\arguments{
\item{fit}{A \code{dpmirt_fit} object.}

\item{vars}{Character vector of variables to extract. Options:
\code{"theta"}, \code{"beta"}, \code{"lambda"}, \code{"delta"}.}

\item{format}{Character. Output format: \code{"matrix"} or \code{"long"}.}

\item{use_rescaled}{Reserved. Only \code{TRUE} is currently supported.}
}
\value{
If one variable is requested in matrix format, a matrix
(\code{n_retained_draws} x N or \code{n_retained_draws} x I). If multiple
variables are requested in matrix format, a named list of matrices. In long
format, a data.frame stacking all requested variables and chain/draw
metadata, including retained MCMC iteration labels after \code{thin} or
\code{thin2}.
}
\description{
Provides convenient access to posterior samples in matrix or long format.
}
\details{
This function provides direct access to the posterior MCMC samples
stored in a \code{dpmirt_fit} object. The matrix format (default) is
efficient for computation, while the long format is convenient for
visualization with ggplot2 and includes chain metadata when available.
Raw draw extraction
(\code{use_rescaled = FALSE}) is reserved until the raw/rescaled draw
schema is finalized.
}
\examples{
\dontrun{
sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
fit <- dpmirt(sim$response, model = "rasch", prior = "normal",
              niter = 5000, nburnin = 1000, seed = 123)

# Extract theta draws as matrix (retained draws x N)
theta <- dpmirt_draws(fit, vars = "theta")
dim(theta)

# Extract beta draws in long format
beta_long <- dpmirt_draws(fit, vars = "beta", format = "long")
head(beta_long)
}

}
\seealso{
\code{\link{dpmirt}}, \code{\link{dpmirt_estimates}}

Other estimation:
\code{\link[=dpmirt_estimates]{dpmirt_estimates()}},
\code{\link[=dpmirt_rescale]{dpmirt_rescale()}}
}
\concept{estimation}
