% Please edit documentation in R/sample.R
\name{dpmirt_resume}
\alias{dpmirt_resume}
\title{Resume MCMC sampling from a previous run}
\usage{
dpmirt_resume(fit_or_compiled, niter_more, reset = FALSE, verbose = TRUE, ...)
}
\arguments{
\item{fit_or_compiled}{A \code{dpmirt_samples}, \code{dpmirt_fit}, or
\code{dpmirt_compiled} object.}

\item{niter_more}{Integer. Number of additional iterations.}

\item{reset}{Logical. FALSE to continue from current state (default),
TRUE to restart.}

\item{verbose}{Logical. Print progress messages.}

\item{...}{Additional arguments.}
}
\value{
A \code{dpmirt_samples} object with extended samples. Pass the
returned object to \code{\link{dpmirt_rescale}} before using fit methods
such as \code{summary()}, \code{plot()}, or \code{dpmirt_estimates()}.
}
\description{
Continues MCMC sampling without recompilation, using the NIMBLE
\code{Cmcmc$run(niter, reset = FALSE)} pattern.
}
\details{
Resume requires the original live compiled NIMBLE object; saved RDS objects
cannot restore compiled external pointers. Multi-chain fitted objects cannot
be continued with \code{reset = FALSE}; use \code{reset = TRUE} for a fresh
single-chain restart or rerun \code{\link{dpmirt}}.
}
\examples{
\dontrun{
# Continue from a previous fit for more iterations
fit <- dpmirt(sim$response, model = "rasch", prior = "normal",
              niter = 5000, nburnin = 1000)
resumed <- dpmirt_resume(fit, niter_more = 5000)
fit2 <- dpmirt_rescale(resumed)
summary(fit2)
}

}
\seealso{
\code{\link{dpmirt_sample}}, \code{\link{dpmirt}}

Other model fitting:
\code{\link[=dpmirt]{dpmirt()}},
\code{\link[=dpmirt_compile]{dpmirt_compile()}},
\code{\link{dpmirt_fit-methods}},
\code{\link[=dpmirt_sample]{dpmirt_sample()}},
\code{\link[=dpmirt_spec]{dpmirt_spec()}}
}
\concept{model fitting}
