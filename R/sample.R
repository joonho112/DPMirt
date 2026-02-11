# ============================================================================
# Module 3: Sampling
# ============================================================================
# Purpose: Run MCMC on compiled model. Lightweight and repeatable.
# Blueprint: Section 6 - Module 3
# ============================================================================

#' Run MCMC sampling on a compiled DPMirt model
#'
#' Executes MCMC sampling using the compiled model. This is the lightweight
#' step that can be called repeatedly from the same compiled model
#' (compile-once, sample-many pattern).
#'
#' @param compiled A \code{dpmirt_compiled} object from
#'   \code{\link{dpmirt_compile}}.
#' @param niter Integer. Total number of MCMC iterations.
#' @param nburnin Integer. Number of burn-in iterations to discard.
#' @param thin Integer. Thinning interval for main monitors.
#' @param thin2 Integer or NULL. Thinning interval for monitors2 (eta/theta).
#'   If NULL, uses the same value as \code{thin}.
#' @param seed Integer or NULL. Random seed for reproducibility.
#' @param reset Logical. If TRUE (default), reset the MCMC state before
#'   sampling. Set to FALSE for chain continuation.
#' @param verbose Logical. Print progress messages.
#' @param ... Additional arguments.
#'
#' @return A \code{dpmirt_samples} S3 object containing:
#' \describe{
#'   \item{samples}{Matrix of posterior samples from main monitors.}
#'   \item{samples2}{Matrix of posterior samples from thinned monitors (eta).}
#'   \item{waic}{WAIC value if computed, otherwise NULL.}
#'   \item{sampling_time}{Time taken for sampling.}
#'   \item{mcmc_control}{List of MCMC settings used.}
#'   \item{model_config}{Reference to model configuration.}
#'   \item{compiled}{Reference to compiled object (for resume).}
#' }
#'
#' @examples
#' \dontrun{
#' sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
#' spec <- dpmirt_spec(sim$response, model = "rasch", prior = "normal")
#' compiled <- dpmirt_compile(spec)
#'
#' # Run MCMC sampling
#' samples <- dpmirt_sample(compiled, niter = 5000, nburnin = 1000, seed = 123)
#' print(samples)
#' }
#'
#' @family model fitting
#' @seealso \code{\link{dpmirt_compile}}, \code{\link{dpmirt_resume}},
#'   \code{\link{dpmirt_rescale}}
#'
#' @export
dpmirt_sample <- function(compiled,
                          niter = 10000L,
                          nburnin = 2000L,
                          thin = 1L,
                          thin2 = NULL,
                          seed = NULL,
                          reset = TRUE,
                          verbose = TRUE,
                          ...) {

  # --- Validate input ---
  if (!inherits(compiled, "dpmirt_compiled")) {
    stop("Input must be a dpmirt_compiled object from dpmirt_compile().",
         call. = FALSE)
  }

  # Validate compiled object is still alive (C++ pointers valid)
  .validate_compiled_alive(compiled)

  # --- Set seed ---
  .set_seed(seed)

  # --- Resolve thin2 ---
  if (is.null(thin2)) {
    thin2 <- thin
  }

  # --- Run MCMC ---
  .vmsg("Running MCMC (niter=", niter, ", nburnin=", nburnin,
        ", thin=", thin, ", thin2=", thin2, ")...", verbose = verbose)

  timer_start <- .start_timer()

  # Use NIMBLE's runMCMC for clean interface
  mcmc_output <- tryCatch(
    runMCMC(
      compiled$Cmcmc,
      niter   = niter,
      nburnin = nburnin,
      thin    = thin,
      thin2   = thin2,
      setSeed = if (!is.null(seed)) seed else FALSE,
      WAIC    = TRUE
    ),
    error = function(e) {
      # If WAIC=TRUE fails (not configured), retry without
      if (grepl("WAIC", conditionMessage(e), ignore.case = TRUE)) {
        .vmsg("  Note: WAIC not available, running without it.",
              verbose = verbose)
        runMCMC(
          compiled$Cmcmc,
          niter   = niter,
          nburnin = nburnin,
          thin    = thin,
          thin2   = thin2,
          setSeed = if (!is.null(seed)) seed else FALSE,
          WAIC    = FALSE
        )
      } else {
        stop(e)
      }
    }
  )

  sampling_time <- .elapsed_time(timer_start)
  .vmsg("  Sampling complete in ", .format_time(sampling_time),
        verbose = verbose)

  # --- Extract samples ---
  # runMCMC returns different formats depending on configuration
  samples_result <- .extract_mcmc_output(mcmc_output, compiled)

  # --- Compute WAIC if enabled ---
  waic_val <- tryCatch(
    compiled$Cmcmc$getWAIC()$WAIC,
    error = function(e) {
      if (verbose) message("  Note: WAIC computation not available.")
      NULL
    }
  )

  # --- Build samples object ---
  result <- structure(
    list(
      samples       = samples_result$samples,
      samples2      = samples_result$samples2,
      waic          = waic_val,
      sampling_time = sampling_time,
      mcmc_control  = list(
        niter   = niter,
        nburnin = nburnin,
        thin    = thin,
        thin2   = thin2,
        seed    = seed,
        reset   = reset
      ),
      model_config = compiled$spec$config,
      compiled     = compiled
    ),
    class = "dpmirt_samples"
  )

  result
}


#' Resume MCMC sampling from a previous run
#'
#' Continues MCMC sampling without recompilation, using the NIMBLE
#' \code{Cmcmc$run(niter, reset = FALSE)} pattern.
#'
#' @param fit_or_compiled A \code{dpmirt_samples}, \code{dpmirt_fit}, or
#'   \code{dpmirt_compiled} object.
#' @param niter_more Integer. Number of additional iterations.
#' @param reset Logical. FALSE to continue from current state (default),
#'   TRUE to restart.
#' @param verbose Logical. Print progress messages.
#' @param ... Additional arguments.
#'
#' @return A \code{dpmirt_samples} object with extended samples.
#'
#' @examples
#' \dontrun{
#' # Continue from a previous fit for more iterations
#' fit <- dpmirt(sim$response, model = "rasch", prior = "normal",
#'               niter = 5000, nburnin = 1000)
#' fit2 <- dpmirt_resume(fit, niter_more = 5000)
#' summary(fit2)
#' }
#'
#' @family model fitting
#' @seealso \code{\link{dpmirt_sample}}, \code{\link{dpmirt}}
#'
#' @export
dpmirt_resume <- function(fit_or_compiled,
                          niter_more,
                          reset = FALSE,
                          verbose = TRUE,
                          ...) {

  # Extract compiled object
  if (inherits(fit_or_compiled, "dpmirt_samples")) {
    compiled <- fit_or_compiled$compiled
  } else if (inherits(fit_or_compiled, "dpmirt_fit")) {
    compiled <- fit_or_compiled$compiled
  } else if (inherits(fit_or_compiled, "dpmirt_compiled")) {
    compiled <- fit_or_compiled
  } else {
    stop("Input must be a dpmirt_samples, dpmirt_fit, or dpmirt_compiled object.",
         call. = FALSE)
  }

  if (is.null(compiled)) {
    stop("No compiled model reference found. Cannot resume. ",
         "Recompile using dpmirt_compile().", call. = FALSE)
  }

  .validate_compiled_alive(compiled)

  .vmsg("Resuming MCMC for ", niter_more, " additional iterations ",
        "(reset=", reset, ")...", verbose = verbose)

  timer_start <- .start_timer()

  # Run additional iterations
  compiled$Cmcmc$run(niter_more, reset = reset)

  sampling_time <- .elapsed_time(timer_start)
  .vmsg("  Resume complete in ", .format_time(sampling_time),
        verbose = verbose)

  # Extract samples from the extended run
  samples  <- as.matrix(compiled$Cmcmc$mvSamples)
  samples2 <- as.matrix(compiled$Cmcmc$mvSamples2)

  # WAIC
  waic_val <- tryCatch(
    compiled$Cmcmc$getWAIC()$WAIC,
    error = function(e) NULL
  )

  result <- structure(
    list(
      samples       = samples,
      samples2      = samples2,
      waic          = waic_val,
      sampling_time = sampling_time,
      mcmc_control  = list(
        niter_more = niter_more,
        reset      = reset,
        resumed    = TRUE
      ),
      model_config = compiled$spec$config,
      compiled     = compiled
    ),
    class = "dpmirt_samples"
  )

  result
}


# ============================================================================
# Internal Helpers
# ============================================================================

#' Validate that a compiled object's C++ pointers are still alive
#' @noRd
.validate_compiled_alive <- function(compiled) {
  tryCatch({
    # Try accessing the model — will error if pointers are dead
    compiled$Cmodel$getNodeNames()
    invisible(TRUE)
  }, error = function(e) {
    stop("Compiled model is no longer valid (C++ pointers expired). ",
         "This can happen if you saved/loaded the object across R sessions. ",
         "Please recompile using dpmirt_compile().",
         call. = FALSE)
  })
}


#' Extract samples from runMCMC output
#'
#' Handles the various output formats from NIMBLE's runMCMC.
#'
#' @noRd
.extract_mcmc_output <- function(mcmc_output, compiled) {

  # runMCMC can return a list with $samples and $samples2
  # or just a matrix (if no monitors2)
  if (is.list(mcmc_output) && !is.null(mcmc_output$samples)) {
    samples  <- mcmc_output$samples
    samples2 <- mcmc_output$samples2
  } else if (is.matrix(mcmc_output)) {
    samples  <- mcmc_output
    samples2 <- NULL
  } else {
    # Try to extract from compiled MCMC directly
    samples  <- as.matrix(compiled$Cmcmc$mvSamples)
    samples2 <- tryCatch(
      as.matrix(compiled$Cmcmc$mvSamples2),
      error = function(e) NULL
    )
  }

  list(samples = samples, samples2 = samples2)
}


# ============================================================================
# S3 Methods for dpmirt_samples
# ============================================================================

#' @rdname dpmirt_sample
#' @param x A \code{dpmirt_samples} object.
#' @param ... Additional arguments (currently unused).
#' @export
print.dpmirt_samples <- function(x, ...) {
  cat("DPMirt MCMC Samples\n")
  cat("===================\n")
  cat("Model:          ", toupper(x$model_config$model), "\n")
  cat("Prior:          ", x$model_config$prior, "\n")

  if (!is.null(x$samples)) {
    cat("Main samples:   ", nrow(x$samples), " iterations x ",
        ncol(x$samples), " parameters\n", sep = "")
  }
  if (!is.null(x$samples2)) {
    cat("Eta samples:    ", nrow(x$samples2), " iterations x ",
        ncol(x$samples2), " persons\n", sep = "")
  }
  if (!is.null(x$waic)) {
    cat("WAIC:           ", round(x$waic, 2), "\n")
  }
  cat("Sampling time:  ", .format_time(x$sampling_time), "\n")
  invisible(x)
}
