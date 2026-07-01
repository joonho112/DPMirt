% Please edit documentation in R/sample.R
\name{dpmirt_sample}
\alias{dpmirt_sample}
\alias{print.dpmirt_samples}
\title{Run MCMC sampling on a compiled DPMirt model}
\usage{
dpmirt_sample(
  compiled,
  niter = 10000L,
  nburnin = 2000L,
  thin = 1L,
  thin2 = NULL,
  seed = NULL,
  reset = TRUE,
  verbose = TRUE,
  inits = NULL,
  init_seed = NULL,
  init_strategy = NULL,
  ...
)

\method{print}{dpmirt_samples}(x, ...)
}
\arguments{
\item{compiled}{A \code{dpmirt_compiled} object from
\code{\link{dpmirt_compile}}.}

\item{niter}{Integer. Total number of MCMC iterations.}

\item{nburnin}{Integer. Number of burn-in iterations to discard.}

\item{thin}{Integer. Thinning interval for main monitors.}

\item{thin2}{Integer or NULL. Thinning interval for monitors2 (eta/theta).
If NULL, uses the same value as \code{thin}.}

\item{seed}{Integer or NULL. Random seed for reproducibility.}

\item{reset}{Logical. If TRUE (default), reset sample and WAIC storage before
sampling. Set to FALSE to continue from the compiled object's current
sampler state and append to existing sample storage.}

\item{verbose}{Logical. Print progress messages.}

\item{inits}{Optional list of initial values to apply before sampling. This
is mainly used internally for chain-specific starts and requires
\code{reset = TRUE}.}

\item{init_seed}{Optional integer seed used to generate \code{inits}, stored
as provenance.}

\item{init_strategy}{Optional character label describing how initial values
were chosen, stored as provenance.}

\item{...}{Additional arguments (currently unused).}

\item{x}{A \code{dpmirt_samples} object.}
}
\value{
A \code{dpmirt_samples} S3 object containing:
\describe{
\item{samples}{Matrix of posterior samples from main monitors.}
\item{samples2}{Matrix of posterior samples from thinned monitors (eta).}
\item{waic}{WAIC value if computed, otherwise NULL.}
\item{sampling_time}{Time taken for sampling.}
\item{mcmc_control}{List of MCMC settings used.}
\item{model_config}{Reference to model configuration.}
\item{compiled}{Reference to compiled object (for resume).}
\item{schema_version}{Draw-storage schema version.}
\item{chain_info}{Data frame describing retained rows by chain.}
\item{draw_index}{Main/theta row maps with chain and MCMC iteration IDs.}
\item{run_history}{Data frame describing sampling runs.}
}
}
\description{
Executes MCMC sampling using the compiled model. This is the lightweight
step that can be called repeatedly from the same compiled model
(compile-once, sample-many pattern).
}
\examples{
\dontrun{
sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
spec <- dpmirt_spec(sim$response, model = "rasch", prior = "normal")
compiled <- dpmirt_compile(spec)

# Run MCMC sampling
samples <- dpmirt_sample(compiled, niter = 5000, nburnin = 1000, seed = 123)
print(samples)
}

}
\seealso{
\code{\link{dpmirt_compile}}, \code{\link{dpmirt_resume}},
\code{\link{dpmirt_rescale}}

Other model fitting:
\code{\link[=dpmirt]{dpmirt()}},
\code{\link[=dpmirt_compile]{dpmirt_compile()}},
\code{\link{dpmirt_fit-methods}},
\code{\link[=dpmirt_resume]{dpmirt_resume()}},
\code{\link[=dpmirt_spec]{dpmirt_spec()}}
}
\concept{model fitting}
