% Please edit documentation in R/rescale.R
\name{dpmirt_rescale}
\alias{dpmirt_rescale}
\title{Rescale posterior samples for identification}
\usage{
dpmirt_rescale(samples_obj, rescale = TRUE)
}
\arguments{
\item{samples_obj}{A \code{dpmirt_samples} object.}

\item{rescale}{Logical. If TRUE (default), apply rescaling.
If FALSE, wrap untransformed posterior samples in a \code{dpmirt_fit}
object using the same downstream fit schema.}
}
\value{
A \code{dpmirt_fit} S3 object containing rescaled posterior
samples, chain/run metadata, diagnostics (ESS, WAIC, cluster info), and
model configuration. This object can be passed directly to
\code{\link{dpmirt_estimates}}, \code{\link{dpmirt_resume}}, and
other downstream functions.
}
\description{
Applies post-hoc rescaling to MCMC posterior samples to resolve
identification indeterminacy. For unconstrained models, this
centers item difficulties and (for 2PL/3PL) normalizes
discriminations.
}
\details{
Post-hoc rescaling resolves the identification indeterminacy inherent
in unconstrained IRT models:

\strong{Rasch} (location only): For each MCMC iteration \eqn{s}:
\deqn{\beta^*_i = \beta_i - \bar{\beta}}
\deqn{\theta^*_j = \theta_j - \bar{\beta}}

\strong{2PL/3PL IRT parameterization} (location + scale): Let
\eqn{c_s = \bar{\beta}_s} and \eqn{d_s = (\prod_i \lambda_i)^{-1/I}}:
\deqn{\beta^*_i = (\beta_i - c_s) / d_s}
\deqn{\lambda^*_i = \lambda_i \cdot d_s}
\deqn{\theta^*_j = (\theta_j - c_s) / d_s}

\strong{2PL/3PL slope-intercept parameterization}: Let
\eqn{c_s = \sum_i \gamma_i / \sum_i \lambda_i} and
\eqn{d_s = (\prod_i \lambda_i)^{-1/I}}:
\deqn{\gamma^*_i = \gamma_i - \lambda_i c_s}
\deqn{\lambda^*_i = \lambda_i \cdot d_s}
\deqn{\theta^*_j = (\theta_j + c_s) / d_s}

Constrained models and calls with \code{rescale = FALSE} use the same
\code{dpmirt_fit} wrapper but preserve the raw sampled scale.

After rescaling: \eqn{\bar{\beta}^* = 0} and
\eqn{(\prod_i \lambda^*_i)^{1/I} = 1}.
}
\examples{
\dontrun{
sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
spec <- dpmirt_spec(sim$response, model = "rasch", prior = "normal",
                    identification = "unconstrained")
compiled <- dpmirt_compile(spec)
samples <- dpmirt_sample(compiled, niter = 5000, nburnin = 1000)

# Apply rescaling
rescaled <- dpmirt_rescale(samples)
}

}
\references{
Paganin, S., Paciorek, C. J., Wehrhahn, C., Rodriguez, A.,
Rabe-Hesketh, S., & de Valpine, P. (2023). Computational strategies
and estimation performance with Bayesian semiparametric item response
theory models. \emph{Journal of Educational and Behavioral Statistics,
48}(2), 147--188.
}
\seealso{
\code{\link{dpmirt_sample}}, \code{\link{dpmirt_estimates}}

Other estimation:
\code{\link[=dpmirt_draws]{dpmirt_draws()}},
\code{\link[=dpmirt_estimates]{dpmirt_estimates()}}
}
\concept{estimation}
