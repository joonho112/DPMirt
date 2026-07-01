% Please edit documentation in R/diagnostics.R
\name{dpmirt_diagnostics}
\alias{dpmirt_diagnostics}
\alias{print.dpmirt_diagnostics}
\title{Compute MCMC Diagnostics for a DPMirt Fit}
\usage{
dpmirt_diagnostics(fit)

\method{print}{dpmirt_diagnostics}(x, ...)
}
\arguments{
\item{fit}{A \code{dpmirt_fit} object from \code{\link{dpmirt}}.}

\item{x}{A \code{dpmirt_diagnostics} object.}

\item{...}{Additional arguments (currently unused).}
}
\value{
A \code{dpmirt_diagnostics} S3 object containing:
\describe{
\item{ess}{List of ESS vectors for items and persons.}
\item{ess_by_chain}{Per-chain ESS data frames when chain metadata exists.}
\item{ess_min_items, ess_min_theta}{Minimum finite item/person ESS.}
\item{rhat}{R-hat vectors by parameter family when computable.}
\item{rhat_max}{Maximum finite R-hat value, or \code{NA}.}
\item{waic}{WAIC value (if computed).}
\item{waic_by_chain}{Per-chain WAIC provenance when available.}
\item{waic_aggregation}{How top-level WAIC was aggregated.}
\item{loglik_trace}{Log-likelihood trace vector.}
\item{loglik_by_chain}{Log-likelihood trace with chain metadata.}
\item{chain_info}{Stored chain/run metadata when available.}
\item{n_clusters}{Posterior cluster counts (DPM only).}
\item{n_clusters_summary}{Summary of posterior cluster counts (DPM only).}
\item{alpha_summary}{Alpha posterior summary (DPM only).}
\item{compilation_time, sampling_time, total_time}{Timing information.}
}
}
\description{
Returns a structured list of diagnostic information including effective
sample sizes (ESS), optional chain-aware R-hat, WAIC, log-likelihood trace,
timing, and DPM-specific cluster diagnostics (number of clusters, alpha
posterior).
}
\details{
Effective sample size (ESS) measures the number of effectively independent
draws from the posterior. Low ESS (< 100) suggests poor mixing and the
need for longer chains or different samplers. R-hat is reported only when
\code{draw_index} contains at least two labeled chains with retained draws.
Sequential resume segments are tracked in \code{run_history} but are not
treated as independent chains. For DPM models, the cluster count trace is a
first-pass visual check: stable oscillation is reassuring, while monotonic
trends suggest the run has not yet mixed.
}
\examples{
\dontrun{
sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
fit <- dpmirt(sim$response, model = "rasch", prior = "dpm",
              niter = 5000, nburnin = 1000, seed = 123)
diag <- dpmirt_diagnostics(fit)
print(diag)
}

}
\seealso{
\code{\link{dpmirt}}, \code{\link{dpmirt_compare}},
\code{\link{plot.dpmirt_fit}}

Other diagnostics:
\code{\link[=dpmirt_compare]{dpmirt_compare()}}
}
\concept{diagnostics}
