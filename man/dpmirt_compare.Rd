% Please edit documentation in R/diagnostics.R
\name{dpmirt_compare}
\alias{dpmirt_compare}
\title{Compare DPMirt models using information criteria}
\usage{
dpmirt_compare(..., criterion = "waic")
}
\arguments{
\item{...}{Two or more \code{dpmirt_fit} objects.}

\item{criterion}{Character. Comparison criterion. Currently only
\code{"waic"} is supported.}
}
\value{
A data.frame ranking models by the criterion with columns
\code{model}, \code{waic}, \code{delta_waic}, and
\code{waic_aggregation}. Fits with unavailable WAIC are kept in the
table with \code{NA} deltas; all-missing WAIC values raise an error.
}
\description{
Compare DPMirt models using information criteria
}
\details{
For chain-labeled multi-run fits, \code{waic_aggregation} records whether
the top-level WAIC is a mean of per-run WAIC values. Such values are
retained for provenance and backward compatibility but should not be
interpreted as pooled posterior WAIC.
}
\examples{
\dontrun{
sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
fit1 <- dpmirt(sim$response, model = "rasch", prior = "normal",
               niter = 5000, nburnin = 1000, seed = 123)
fit2 <- dpmirt(sim$response, model = "rasch", prior = "dpm",
               niter = 5000, nburnin = 1000, seed = 123)

# Compare via WAIC
comp <- dpmirt_compare(fit1, fit2)
print(comp)
}

}
\seealso{
\code{\link{dpmirt}}, \code{\link{dpmirt_diagnostics}}

Other diagnostics:
\code{\link[=dpmirt_diagnostics]{dpmirt_diagnostics()}}
}
\concept{diagnostics}
