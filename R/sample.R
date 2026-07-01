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
#' @param reset Logical. If TRUE (default), reset sample and WAIC storage before
#'   sampling. Set to FALSE to continue from the compiled object's current
#'   sampler state and append to existing sample storage.
#' @param verbose Logical. Print progress messages.
#' @param inits Optional list of initial values to apply before sampling. This
#'   is mainly used internally for chain-specific starts and requires
#'   \code{reset = TRUE}.
#' @param init_seed Optional integer seed used to generate \code{inits}, stored
#'   as provenance.
#' @param init_strategy Optional character label describing how initial values
#'   were chosen, stored as provenance.
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
#'   \item{schema_version}{Draw-storage schema version.}
#'   \item{chain_info}{Data frame describing retained rows by chain.}
#'   \item{draw_index}{Main/theta row maps with chain and MCMC iteration IDs.}
#'   \item{run_history}{Data frame describing sampling runs.}
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
                          inits = NULL,
                          init_seed = NULL,
                          init_strategy = NULL,
                          ...) {

  # --- Validate input ---
  if (!inherits(compiled, "dpmirt_compiled")) {
    stop("Input must be a dpmirt_compiled object from dpmirt_compile().",
         call. = FALSE)
  }

  control <- .validate_mcmc_control_args(
    niter = niter,
    nburnin = nburnin,
    thin = thin,
    thin2 = thin2
  )
  niter <- control$niter
  nburnin <- control$nburnin
  thin <- control$thin
  thin2 <- control$thin2
  if (is.null(init_strategy)) {
    init_strategy <- if (is.null(inits)) "compiled_spec" else "supplied"
  } else {
    init_strategy <- as.character(init_strategy)[[1]]
  }

  # Validate compiled object is still alive (C++ pointers valid)
  .validate_compiled_alive(compiled)

  # --- Set seed ---
  .set_seed(seed)

  # --- Run MCMC ---
  .vmsg("Running MCMC (niter=", niter, ", nburnin=", nburnin,
        ", thin=", thin, ", thin2=", thin2, ")...", verbose = verbose)

  timer_start <- .start_timer()

  .run_compiled_mcmc(
    compiled = compiled,
    niter = niter,
    nburnin = nburnin,
    thin = thin,
    thin2 = thin2,
    reset = reset,
    resetWAIC = reset,
    inits = inits
  )

  sampling_time <- .elapsed_time(timer_start)
  .vmsg("  Sampling complete in ", .format_time(sampling_time),
        verbose = verbose)

  # --- Extract samples ---
  samples_result <- .extract_compiled_samples(compiled)

  # --- Compute WAIC if enabled ---
  waic_val <- .get_compiled_waic(compiled, verbose = verbose)

  # --- Build samples object ---
  niter_total_retained <- if (isTRUE(reset)) {
    NULL
  } else {
    nrow(samples_result$samples)
  }
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
        reset   = reset,
        resumed = !isTRUE(reset),
        init_seed = init_seed,
        init_strategy = init_strategy,
        niter_requested = niter,
        niter_total_retained = niter_total_retained,
        n_draws_main = nrow(samples_result$samples),
        n_draws_theta = .n_draw_rows(samples_result$samples2),
        waic_enabled = .compiled_waic_enabled(compiled)
      ),
      model_config = compiled$spec$config,
      compiled     = compiled
    ),
    class = "dpmirt_samples"
  )

  .add_chain_metadata(result, chain_id = 1L, seed = seed)
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
#' @return A \code{dpmirt_samples} object with extended samples. Pass the
#'   returned object to \code{\link{dpmirt_rescale}} before using fit methods
#'   such as \code{summary()}, \code{plot()}, or \code{dpmirt_estimates()}.
#'
#' @details
#' Resume requires the original live compiled NIMBLE object; saved RDS objects
#' cannot restore compiled external pointers. Multi-chain fitted objects cannot
#' be continued with \code{reset = FALSE}; use \code{reset = TRUE} for a fresh
#' single-chain restart or rerun \code{\link{dpmirt}}.
#'
#' @examples
#' \dontrun{
#' # Continue from a previous fit for more iterations
#' fit <- dpmirt(sim$response, model = "rasch", prior = "normal",
#'               niter = 5000, nburnin = 1000)
#' resumed <- dpmirt_resume(fit, niter_more = 5000)
#' fit2 <- dpmirt_rescale(resumed)
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

  # Extract compiled object and previous controls
  previous_control <- list()
  if (inherits(fit_or_compiled, "dpmirt_samples")) {
    compiled <- fit_or_compiled$compiled
    previous_control <- fit_or_compiled$mcmc_control
  } else if (inherits(fit_or_compiled, "dpmirt_fit")) {
    compiled <- fit_or_compiled$compiled
    previous_control <- fit_or_compiled$config
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

  resume_nchains <- .resume_input_nchains(fit_or_compiled, previous_control)
  if (!isTRUE(reset) && resume_nchains > 1L) {
    stop(
      "Multi-chain fitted objects cannot be resumed with reset = FALSE ",
      "because DPMirt stores one live compiled sampler, not one sampler per ",
      "chain. Re-run dpmirt() or use reset = TRUE for a fresh single-chain ",
      "restart.",
      call. = FALSE
    )
  }

  .validate_compiled_alive(compiled)
  .validate_resume_state(fit_or_compiled, compiled, reset = reset)

  thin <- .control_get(previous_control, "thin", 1L)
  thin2 <- .control_get(previous_control, "thin2", thin)
  control <- .validate_mcmc_control_args(
    niter = niter_more,
    nburnin = 0L,
    thin = thin,
    thin2 = thin2
  )
  niter_more <- control$niter
  thin <- control$thin
  thin2 <- control$thin2

  .vmsg("Resuming MCMC for ", niter_more, " additional iterations ",
        "(reset=", reset, ")...", verbose = verbose)

  reset_inits <- if (isTRUE(reset)) {
    .resume_reset_inits(compiled, previous_control)
  } else {
    NULL
  }
  reset_init_seed <- if (isTRUE(reset) && !is.null(reset_inits)) {
    .control_get(previous_control, "seed", NULL)
  } else {
    NULL
  }
  reset_init_strategy <- if (isTRUE(reset) && !is.null(reset_inits)) {
    if (is.null(reset_init_seed)) "compiled_spec_reset" else "chain_seeded_reset"
  } else {
    .control_get(previous_control, "init_strategy", NULL)
  }

  timer_start <- .start_timer()

  # Run additional iterations
  .run_compiled_mcmc(
    compiled = compiled,
    niter = niter_more,
    nburnin = 0L,
    thin = thin,
    thin2 = thin2,
    reset = reset,
    resetWAIC = reset,
    inits = reset_inits
  )

  sampling_time <- .elapsed_time(timer_start)
  .vmsg("  Resume complete in ", .format_time(sampling_time),
        verbose = verbose)

  # Extract samples from the extended run
  samples_result <- .extract_compiled_samples(compiled)
  samples <- samples_result$samples
  samples2 <- samples_result$samples2

  # WAIC
  waic_val <- .get_compiled_waic(compiled, verbose = verbose)
  previous_mcmc_end <- if (isTRUE(reset)) {
    .control_get(previous_control, "niter", NULL)
  } else {
    .previous_mcmc_end(fit_or_compiled, previous_control)
  }
  cumulative_niter <- if (isTRUE(reset)) {
    niter_more
  } else {
    previous_mcmc_end + niter_more
  }

  result <- structure(
    list(
      samples       = samples,
      samples2      = samples2,
      waic          = waic_val,
      sampling_time = sampling_time,
      mcmc_control  = list(
        niter = cumulative_niter,
        niter_more = niter_more,
        niter_requested = niter_more,
        niter_total_retained = nrow(samples),
        n_draws_total_retained = nrow(samples),
        nburnin = 0L,
        thin = thin,
        thin2 = thin2,
        seed = .control_get(previous_control, "seed", NULL),
        reset      = reset,
        resumed    = !isTRUE(reset),
        init_seed = reset_init_seed,
        init_strategy = reset_init_strategy,
        previous_niter = previous_mcmc_end,
        previous_n_draws_main = .control_get(
          previous_control, "n_draws_main", NULL
        ),
        n_draws_main = nrow(samples),
        n_draws_theta = .n_draw_rows(samples2),
        waic_enabled = .compiled_waic_enabled(compiled)
      ),
      model_config = compiled$spec$config,
      compiled     = compiled
    ),
    class = "dpmirt_samples"
  )

  .add_resume_metadata(
    result,
    previous_obj = fit_or_compiled,
    previous_control = previous_control,
    reset = reset,
    chain_id = 1L
  )
}


# ============================================================================
# Internal Helpers
# ============================================================================

#' Validate MCMC iteration, thinning, and chain controls
#' @noRd
.validate_mcmc_control_args <- function(niter,
                                        nburnin = 0L,
                                        thin = 1L,
                                        thin2 = NULL,
                                        nchains = NULL) {
  if (!.is_single_whole_number(niter) || niter <= 0) {
    stop("niter must be a single positive whole number.", call. = FALSE)
  }
  if (!.is_single_whole_number(nburnin) || nburnin < 0) {
    stop("nburnin must be a single non-negative whole number.",
         call. = FALSE)
  }
  if (!.is_single_whole_number(thin) || thin <= 0) {
    stop("thin must be a single positive whole number.", call. = FALSE)
  }
  if (is.null(thin2)) {
    thin2 <- thin
  }
  if (!.is_single_whole_number(thin2) || thin2 <= 0) {
    stop("thin2 must be NULL or a single positive whole number.",
         call. = FALSE)
  }
  if (!is.null(nchains) &&
      (!.is_single_whole_number(nchains) || nchains <= 0)) {
    stop("nchains must be a single positive whole number.", call. = FALSE)
  }

  niter <- as.integer(niter)
  nburnin <- as.integer(nburnin)
  thin <- as.integer(thin)
  thin2 <- as.integer(thin2)

  if (niter <= nburnin) {
    stop("niter must be greater than nburnin.", call. = FALSE)
  }
  if (floor((niter - nburnin) / thin) < 1L) {
    stop(
      "MCMC settings retain zero main draws; use niter, nburnin, and thin ",
      "so floor((niter - nburnin) / thin) >= 1.",
      call. = FALSE
    )
  }
  if (floor((niter - nburnin) / thin2) < 1L) {
    stop(
      "MCMC settings retain zero theta draws; use niter, nburnin, and thin2 ",
      "so floor((niter - nburnin) / thin2) >= 1.",
      call. = FALSE
    )
  }

  out <- list(
    niter = niter,
    nburnin = nburnin,
    thin = thin,
    thin2 = thin2
  )
  if (!is.null(nchains)) {
    out$nchains <- as.integer(nchains)
  }
  out
}


#' @noRd
.is_single_whole_number <- function(x) {
  is.numeric(x) &&
    length(x) == 1L &&
    !is.na(x) &&
    is.finite(x) &&
    x <= .Machine$integer.max &&
    abs(x - round(x)) < sqrt(.Machine$double.eps)
}


#' Validate that a compiled object's C++ pointers are still alive
#' @noRd
.validate_compiled_alive <- function(compiled) {
  tryCatch({
    # Try accessing the model; this errors if pointers are dead.
    compiled$Cmodel$getNodeNames()
    if (is.null(compiled$Cmcmc) || is.null(compiled$Cmcmc$run)) {
      stop("compiled MCMC pointer is missing")
    }
    invisible(TRUE)
  }, error = function(e) {
    .stop_compiled_expired()
  })
}


#' @noRd
.stop_compiled_expired <- function() {
  stop(
    "Compiled model is no longer valid (C++ pointers expired). ",
    "This can happen if you saved/loaded the object across R sessions. ",
    "Please recompile using dpmirt_compile().",
    call. = FALSE
  )
}


#' @noRd
.is_compiled_pointer_error <- function(e) {
  msg <- conditionMessage(e)
  grepl(
    "C symbol|not in load table|external pointer|NULL value|bad_weak_ptr|expired",
    msg,
    ignore.case = TRUE
  )
}


#' Run a compiled NIMBLE MCMC with explicit reset controls
#' @noRd
.run_compiled_mcmc <- function(compiled,
                               niter,
                               nburnin = 0L,
                               thin = 1L,
                               thin2 = thin,
                               reset = TRUE,
                               resetWAIC = reset,
                               inits = NULL) {
  control <- .validate_mcmc_control_args(
    niter = niter,
    nburnin = nburnin,
    thin = thin,
    thin2 = thin2
  )
  niter <- control$niter
  nburnin <- control$nburnin
  thin <- control$thin
  thin2 <- control$thin2

  if (!is.null(inits)) {
    if (!isTRUE(reset)) {
      stop("inits can only be supplied when reset = TRUE.", call. = FALSE)
    }
    .set_compiled_inits(compiled, inits)
  }

  run_args <- list(niter = as.integer(niter))
  run_formals <- tryCatch(
    names(formals(compiled$Cmcmc$run)),
    error = function(e) character(0)
  )

  add_if_supported <- function(name, value) {
    if (length(run_formals) == 0L || name %in% run_formals) {
      run_args[[name]] <<- value
    }
  }

  add_if_supported("nburnin", as.integer(nburnin))
  add_if_supported("thin", as.integer(thin))
  add_if_supported("thin2", as.integer(thin2))
  add_if_supported("reset", isTRUE(reset))
  add_if_supported("resetMV", isTRUE(reset))
  add_if_supported("resetWAIC", isTRUE(resetWAIC))
  add_if_supported("initializeModel", isTRUE(reset))

  tryCatch(
    do.call(compiled$Cmcmc$run, run_args),
    error = function(e) {
      if (.is_compiled_pointer_error(e)) {
        .stop_compiled_expired()
      }
      stop(e)
    }
  )

  invisible(TRUE)
}


#' Apply initial values to a compiled NIMBLE model
#' @noRd
.set_compiled_inits <- function(compiled, inits) {
  if (!is.list(inits)) {
    stop("inits must be a list of NIMBLE initial values.", call. = FALSE)
  }

  set_inits <- .compiled_set_inits_function(compiled)
  if (is.null(set_inits)) {
    stop("The compiled NIMBLE model does not expose setInits().",
         call. = FALSE)
  }

  tryCatch(
    set_inits(inits),
    error = function(e) {
      if (.is_compiled_pointer_error(e)) {
        .stop_compiled_expired()
      }
      stop("Failed to set initial values on the compiled NIMBLE model: ",
           conditionMessage(e), call. = FALSE)
    }
  )

  invisible(TRUE)
}


#' Locate the compiled model setInits function across NIMBLE object layouts
#' @noRd
.compiled_set_inits_function <- function(compiled) {
  candidates <- list(
    compiled$Cmodel,
    compiled$Cmcmc$Robject$model$CobjectInterface,
    compiled$Cmcmc$Robject$model
  )

  for (candidate in candidates) {
    set_inits <- tryCatch(
      {
        if (is.null(candidate)) {
          NULL
        } else {
          candidate$setInits
        }
      },
      error = function(e) NULL
    )
    if (is.function(set_inits)) {
      return(set_inits)
    }
  }

  NULL
}


#' Extract samples from the compiled NIMBLE MCMC object
#' @noRd
.extract_compiled_samples <- function(compiled) {
  tryCatch({
    samples <- as.matrix(compiled$Cmcmc$mvSamples)
    samples2 <- tryCatch(
      as.matrix(compiled$Cmcmc$mvSamples2),
      error = function(e) NULL
    )
    list(samples = samples, samples2 = samples2)
  }, error = function(e) {
    if (.is_compiled_pointer_error(e)) {
      .stop_compiled_expired()
    }
    stop(e)
  })
}


#' @noRd
.compiled_waic_enabled <- function(compiled) {
  if (is.null(compiled$waic_enabled)) {
    return(TRUE)
  }
  isTRUE(compiled$waic_enabled)
}


#' @noRd
.normalize_waic <- function(waic) {
  if (is.null(waic) || length(waic) == 0L) {
    return(NULL)
  }
  waic <- suppressWarnings(as.numeric(waic[[1]]))
  if (!is.finite(waic)) {
    return(NULL)
  }
  waic
}


#' @noRd
.get_compiled_waic <- function(compiled, verbose = FALSE) {
  if (!.compiled_waic_enabled(compiled)) {
    return(NULL)
  }

  out <- tryCatch(
    compiled$Cmcmc$getWAIC()$WAIC,
    error = function(e) {
      if (.is_compiled_pointer_error(e)) {
        .stop_compiled_expired()
      }
      if (verbose) {
        message("  Note: WAIC computation not available.")
      }
      NULL
    }
  )
  .normalize_waic(out)
}


#' @noRd
.control_get <- function(x, name, default = NULL) {
  if (is.null(x) || !(name %in% names(x))) {
    return(default)
  }
  x[[name]]
}


#' @noRd
.resume_input_nchains <- function(fit_or_compiled, previous_control = list()) {
  vals <- integer(0)

  control_n <- suppressWarnings(
    as.integer(.control_get(previous_control, "nchains", NA_integer_))
  )
  if (!is.na(control_n)) {
    vals <- c(vals, control_n)
  }

  chain_info <- fit_or_compiled$chain_info
  if (!is.null(chain_info) && !is.null(chain_info$chain_id)) {
    vals <- c(vals, length(unique(chain_info$chain_id)))
  }

  draw_index <- fit_or_compiled$draw_index
  if (!is.null(draw_index$main) && !is.null(draw_index$main$chain_id)) {
    vals <- c(vals, length(unique(draw_index$main$chain_id)))
  }

  vals <- vals[!is.na(vals)]
  if (length(vals) == 0L) {
    return(1L)
  }
  max(vals)
}


#' @noRd
.resume_reset_inits <- function(compiled, previous_control = list()) {
  if (is.null(compiled$spec)) {
    return(NULL)
  }

  seed <- .control_get(previous_control, "seed", NULL)
  if (!is.null(seed)) {
    return(.generate_chain_inits(compiled$spec, chain_id = 1L, seed = seed))
  }

  compiled$spec$inits
}


#' Validate that resume input still matches the live compiled sample state
#' @noRd
.validate_resume_state <- function(fit_or_compiled, compiled, reset = FALSE) {
  if (isTRUE(reset) || inherits(fit_or_compiled, "dpmirt_compiled")) {
    return(invisible(TRUE))
  }

  stored_samples <- NULL
  stored_samples2 <- NULL
  if (inherits(fit_or_compiled, "dpmirt_samples")) {
    stored_samples <- fit_or_compiled$samples
    stored_samples2 <- fit_or_compiled$samples2
  } else if (inherits(fit_or_compiled, "dpmirt_fit")) {
    stored_samples <- fit_or_compiled$samples_raw
    stored_samples2 <- fit_or_compiled$samples2_raw
  }

  if (is.null(stored_samples)) {
    return(invisible(TRUE))
  }

  live <- .extract_compiled_samples(compiled)
  if (!.samples_match(live$samples, stored_samples) ||
      !.samples_match(live$samples2, stored_samples2)) {
    stop(
      "Cannot resume because the live compiled MCMC state no longer matches ",
      "the supplied sample object. Use the most recent dpmirt_samples object, ",
      "set reset = TRUE, or recompile and sample again.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


.samples_match <- function(live, stored) {
  if (is.null(stored)) {
    return(TRUE)
  }
  if (is.null(live)) {
    return(FALSE)
  }
  nrow(live) == nrow(stored) &&
    isTRUE(all.equal(live, stored, check.attributes = FALSE))
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
# Chain-Aware Metadata Helpers
# ============================================================================

.schema_version <- function() "chain-aware-v1"


.scalar_or_na <- function(x, na_value) {
  if (is.null(x) || length(x) == 0) {
    return(na_value)
  }
  x[[1]]
}


.scalar_int_or_na <- function(x) {
  val <- .scalar_or_na(x, NA_integer_)
  if (is.na(val)) NA_integer_ else as.integer(val)
}


.scalar_num_or_na <- function(x) {
  val <- .scalar_or_na(x, NA_real_)
  if (is.na(val)) NA_real_ else as.numeric(val)
}


.scalar_lgl_or_na <- function(x) {
  val <- .scalar_or_na(x, NA)
  if (is.na(val)) NA else as.logical(val)
}


.scalar_chr_or_na <- function(x) {
  val <- .scalar_or_na(x, NA_character_)
  if (is.na(val)) NA_character_ else as.character(val)
}


.n_draw_rows <- function(x) {
  if (is.null(x)) {
    return(0L)
  }
  as.integer(nrow(x))
}


.row_start <- function(n_draws, offset = 0L) {
  if (n_draws > 0L) as.integer(offset + 1L) else NA_integer_
}


.row_end <- function(n_draws, offset = 0L) {
  if (n_draws > 0L) as.integer(offset + n_draws) else NA_integer_
}


.retained_mcmc_iterations <- function(n_draws, nburnin, thin) {
  if (n_draws <= 0L) {
    return(integer(0))
  }
  nburnin <- .scalar_int_or_na(nburnin)
  thin <- .scalar_int_or_na(thin)
  if (is.na(nburnin)) {
    nburnin <- 0L
  }
  if (is.na(thin) || thin <= 0L) {
    thin <- 1L
  }
  as.integer(nburnin + seq.int(from = thin, by = thin, length.out = n_draws))
}


.make_draw_index <- function(n_draws,
                             chain_id,
                             nburnin,
                             thin,
                             offset = 0L,
                             iteration_start = NULL) {
  if (n_draws <= 0L) {
    return(data.frame(
      row = integer(0),
      chain_id = integer(0),
      within_chain_draw = integer(0),
      mcmc_iteration = integer(0)
    ))
  }

  mcmc_iteration <- if (is.null(iteration_start)) {
    .retained_mcmc_iterations(n_draws, nburnin, thin)
  } else {
    as.integer(seq.int(
      from = as.integer(iteration_start),
      by = as.integer(thin),
      length.out = n_draws
    ))
  }

  data.frame(
    row = as.integer(offset + seq_len(n_draws)),
    chain_id = rep(as.integer(chain_id), n_draws),
    within_chain_draw = seq_len(n_draws),
    mcmc_iteration = mcmc_iteration
  )
}


.first_or_na_int <- function(x) {
  if (length(x) == 0L) NA_integer_ else as.integer(x[[1]])
}


.last_or_na_int <- function(x) {
  if (length(x) == 0L) NA_integer_ else as.integer(x[[length(x)]])
}


.max_int_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_integer_ else as.integer(max(x))
}


.previous_mcmc_end <- function(previous_obj, previous_control = NULL) {
  candidates <- c(
    previous_obj$run_history$mcmc_end_main,
    previous_obj$run_history$mcmc_end_theta,
    previous_obj$draw_index$main$mcmc_iteration,
    previous_obj$draw_index$theta$mcmc_iteration
  )
  out <- .max_int_or_na(candidates)
  if (!is.na(out)) {
    return(out)
  }

  out <- .scalar_int_or_na(.control_get(previous_control, "niter", NULL))
  if (!is.na(out)) {
    return(out)
  }

  0L
}


.run_history_row <- function(run_id,
                             segment_id,
                             chain_id,
                             seed,
                             reset,
                             resumed,
                             niter,
                             nburnin,
                             thin,
                             thin2,
                             init_seed = NULL,
                             init_strategy = NULL,
                             main_index,
                             theta_index) {
  data.frame(
    run_id = as.integer(run_id),
    chain_id = as.integer(chain_id),
    segment_id = as.integer(segment_id),
    seed = .scalar_int_or_na(seed),
    reset = as.logical(reset),
    resumed = as.logical(resumed),
    niter = .scalar_int_or_na(niter),
    nburnin = .scalar_int_or_na(nburnin),
    thin = .scalar_int_or_na(thin),
    thin2 = .scalar_int_or_na(thin2),
    init_seed = .scalar_int_or_na(init_seed),
    init_strategy = .scalar_chr_or_na(init_strategy),
    n_draws_main = nrow(main_index),
    n_draws_theta = nrow(theta_index),
    row_start_main = .first_or_na_int(main_index$row),
    row_end_main = .last_or_na_int(main_index$row),
    row_start_theta = .first_or_na_int(theta_index$row),
    row_end_theta = .last_or_na_int(theta_index$row),
    mcmc_start_main = .first_or_na_int(main_index$mcmc_iteration),
    mcmc_end_main = .last_or_na_int(main_index$mcmc_iteration),
    mcmc_start_theta = .first_or_na_int(theta_index$mcmc_iteration),
    mcmc_end_theta = .last_or_na_int(theta_index$mcmc_iteration)
  )
}


.bind_run_history <- function(...) {
  pieces <- list(...)
  pieces <- pieces[!vapply(pieces, is.null, logical(1))]
  if (length(pieces) == 0L) {
    return(NULL)
  }
  all_names <- unique(unlist(lapply(pieces, names), use.names = FALSE))
  pieces <- lapply(pieces, function(x) {
    missing <- setdiff(all_names, names(x))
    for (nm in missing) {
      x[[nm]] <- NA
    }
    x[, all_names, drop = FALSE]
  })
  do.call(rbind, pieces)
}


.sample_chain_metadata <- function(samples_obj,
                                   chain_id = 1L,
                                   seed = NULL,
                                   resumed = FALSE,
                                   main_offset = 0L,
                                   theta_offset = 0L) {
  ctrl <- samples_obj$mcmc_control
  if (is.null(ctrl)) {
    ctrl <- list()
  }

  if (is.null(seed)) {
    seed <- .control_get(ctrl, "seed", NULL)
  }

  thin <- .scalar_int_or_na(.control_get(ctrl, "thin", NULL))
  if (is.na(thin)) {
    thin <- 1L
  }
  thin2 <- .scalar_int_or_na(.control_get(ctrl, "thin2", NULL))
  if (is.na(thin2)) {
    thin2 <- thin
  }

  n_main <- .n_draw_rows(samples_obj$samples)
  n_theta <- .n_draw_rows(samples_obj$samples2)

  chain_info <- data.frame(
    chain_id = as.integer(chain_id),
    seed = .scalar_int_or_na(seed),
    niter = .scalar_int_or_na(.control_get(ctrl, "niter", NULL)),
    nburnin = .scalar_int_or_na(.control_get(ctrl, "nburnin", NULL)),
    thin = thin,
    thin2 = thin2,
    reset = .scalar_lgl_or_na(.control_get(ctrl, "reset", NULL)),
    resumed = isTRUE(resumed) || isTRUE(.control_get(ctrl, "resumed", FALSE)),
    init_seed = .scalar_int_or_na(.control_get(ctrl, "init_seed", NULL)),
    init_strategy = .scalar_chr_or_na(
      .control_get(ctrl, "init_strategy", NULL)
    ),
    n_draws_main = n_main,
    n_draws_theta = n_theta,
    row_start_main = .row_start(n_main, main_offset),
    row_end_main = .row_end(n_main, main_offset),
    row_start_theta = .row_start(n_theta, theta_offset),
    row_end_theta = .row_end(n_theta, theta_offset),
    sampling_time = .scalar_num_or_na(samples_obj$sampling_time),
    waic = .scalar_num_or_na(samples_obj$waic)
  )

  draw_index <- list(
    main = .make_draw_index(
      n_main, chain_id, .control_get(ctrl, "nburnin", NULL), thin,
      offset = main_offset
    ),
    theta = .make_draw_index(
      n_theta, chain_id, .control_get(ctrl, "nburnin", NULL), thin2,
      offset = theta_offset
    )
  )

  run_history <- .run_history_row(
    run_id = 1L,
    segment_id = 1L,
    chain_id = chain_id,
    seed = seed,
    reset = chain_info$reset,
    resumed = chain_info$resumed,
    niter = chain_info$niter,
    nburnin = chain_info$nburnin,
    thin = chain_info$thin,
    thin2 = chain_info$thin2,
    init_seed = chain_info$init_seed,
    init_strategy = chain_info$init_strategy,
    main_index = draw_index$main,
    theta_index = draw_index$theta
  )

  list(
    chain_info = chain_info,
    draw_index = draw_index,
    run_history = run_history
  )
}


.add_chain_metadata <- function(samples_obj,
                                chain_id = 1L,
                                seed = NULL,
                                resumed = FALSE) {
  meta <- .sample_chain_metadata(
    samples_obj,
    chain_id = chain_id,
    seed = seed,
    resumed = resumed
  )
  samples_obj$schema_version <- .schema_version()
  samples_obj$chain_info <- meta$chain_info
  samples_obj$draw_index <- meta$draw_index
  samples_obj$run_history <- meta$run_history
  samples_obj
}


.resume_previous_draw_index <- function(previous_obj,
                                        previous_control,
                                        family = c("main", "theta"),
                                        chain_id = 1L) {
  family <- match.arg(family)

  if (!is.null(previous_obj$draw_index) &&
      !is.null(previous_obj$draw_index[[family]])) {
    return(previous_obj$draw_index[[family]])
  }

  if (identical(family, "main")) {
    n_draws <- .control_get(previous_control, "n_draws_main", NULL)
    samples <- if (!is.null(previous_obj$samples)) {
      previous_obj$samples
    } else {
      previous_obj$samples_raw
    }
    if (is.null(n_draws)) {
      n_draws <- .n_draw_rows(samples)
    }
    thin <- .control_get(previous_control, "thin", 1L)
  } else {
    n_draws <- .control_get(previous_control, "n_draws_theta", NULL)
    samples <- if (!is.null(previous_obj$samples2)) {
      previous_obj$samples2
    } else {
      previous_obj$samples2_raw
    }
    if (is.null(n_draws)) {
      n_draws <- .n_draw_rows(samples)
    }
    thin <- .control_get(
      previous_control, "thin2",
      .control_get(previous_control, "thin", 1L)
    )
  }

  .make_draw_index(
    n_draws = as.integer(n_draws),
    chain_id = chain_id,
    nburnin = .control_get(previous_control, "nburnin", 0L),
    thin = thin
  )
}


.add_resume_metadata <- function(samples_obj,
                                 previous_obj,
                                 previous_control,
                                 reset = FALSE,
                                 chain_id = 1L) {
  if (isTRUE(reset) || inherits(previous_obj, "dpmirt_compiled")) {
    return(.add_chain_metadata(samples_obj, chain_id = chain_id,
                               resumed = !isTRUE(reset)))
  }

  previous_main <- .resume_previous_draw_index(
    previous_obj, previous_control, family = "main", chain_id = chain_id
  )
  previous_theta <- .resume_previous_draw_index(
    previous_obj, previous_control, family = "theta", chain_id = chain_id
  )

  n_main_total <- .n_draw_rows(samples_obj$samples)
  n_theta_total <- .n_draw_rows(samples_obj$samples2)
  n_main_previous <- nrow(previous_main)
  n_theta_previous <- nrow(previous_theta)
  n_main_new <- n_main_total - n_main_previous
  n_theta_new <- n_theta_total - n_theta_previous

  if (n_main_new < 0L || n_theta_new < 0L) {
    return(.add_chain_metadata(samples_obj, chain_id = chain_id,
                               resumed = TRUE))
  }

  thin <- .control_get(samples_obj$mcmc_control, "thin", 1L)
  thin2 <- .control_get(samples_obj$mcmc_control, "thin2", thin)
  previous_niter <- .previous_mcmc_end(previous_obj, previous_control)

  appended_main <- .make_draw_index(
    n_draws = n_main_new,
    chain_id = chain_id,
    nburnin = 0L,
    thin = thin,
    offset = n_main_previous,
    iteration_start = previous_niter + thin
  )
  appended_theta <- .make_draw_index(
    n_draws = n_theta_new,
    chain_id = chain_id,
    nburnin = 0L,
    thin = thin2,
    offset = n_theta_previous,
    iteration_start = previous_niter + thin2
  )

  base_meta <- .sample_chain_metadata(
    samples_obj,
    chain_id = chain_id,
    resumed = TRUE
  )

  previous_history <- previous_obj$run_history
  if (is.null(previous_history)) {
    previous_history <- .run_history_row(
      run_id = 1L,
      segment_id = 1L,
      chain_id = chain_id,
      seed = .control_get(previous_control, "seed", NULL),
      reset = .control_get(previous_control, "reset", TRUE),
      resumed = .control_get(previous_control, "resumed", FALSE),
      niter = previous_niter,
      nburnin = .control_get(previous_control, "nburnin", NULL),
      thin = .control_get(previous_control, "thin", 1L),
      thin2 = .control_get(
        previous_control, "thin2",
        .control_get(previous_control, "thin", 1L)
      ),
      init_seed = .control_get(previous_control, "init_seed", NULL),
      init_strategy = .control_get(previous_control, "init_strategy", NULL),
      main_index = previous_main,
      theta_index = previous_theta
    )
  }

  next_run_id <- max(
    .max_int_or_na(previous_history$run_id),
    0L,
    na.rm = TRUE
  ) + 1L
  next_segment_id <- max(
    .max_int_or_na(previous_history$segment_id),
    0L,
    na.rm = TRUE
  ) + 1L

  appended_history <- .run_history_row(
    run_id = next_run_id,
    segment_id = next_segment_id,
    chain_id = chain_id,
    seed = .control_get(samples_obj$mcmc_control, "seed", NULL),
    reset = FALSE,
    resumed = TRUE,
    niter = .control_get(samples_obj$mcmc_control, "niter_more", NULL),
    nburnin = 0L,
    thin = thin,
    thin2 = thin2,
    init_seed = .control_get(samples_obj$mcmc_control, "init_seed", NULL),
    init_strategy = .control_get(
      samples_obj$mcmc_control, "init_strategy", NULL
    ),
    main_index = appended_main,
    theta_index = appended_theta
  )

  samples_obj$schema_version <- .schema_version()
  samples_obj$chain_info <- base_meta$chain_info
  samples_obj$draw_index <- list(
    main = rbind(previous_main, appended_main),
    theta = rbind(previous_theta, appended_theta)
  )
  samples_obj$run_history <- .bind_run_history(
    previous_history,
    appended_history
  )

  samples_obj
}


.ensure_sample_chain_metadata <- function(samples_obj,
                                          fallback_chain_id = 1L,
                                          main_offset = 0L,
                                          theta_offset = 0L) {
  if (!is.null(samples_obj$schema_version) &&
      identical(samples_obj$schema_version, .schema_version()) &&
      !is.null(samples_obj$chain_info) &&
      !is.null(samples_obj$draw_index)) {
    if (main_offset == 0L && theta_offset == 0L) {
      return(samples_obj)
    }
  }

  meta <- .sample_chain_metadata(
    samples_obj,
    chain_id = fallback_chain_id,
    main_offset = main_offset,
    theta_offset = theta_offset
  )
  samples_obj$schema_version <- .schema_version()
  samples_obj$chain_info <- meta$chain_info
  samples_obj$draw_index <- meta$draw_index
  samples_obj$run_history <- meta$run_history
  samples_obj
}


.combine_chain_metadata <- function(all_samples) {
  chain_info <- list()
  main_index <- list()
  theta_index <- list()
  run_history <- list()
  main_offset <- 0L
  theta_offset <- 0L

  for (i in seq_along(all_samples)) {
    chain_id <- i
    if (!is.null(all_samples[[i]]$chain_info$chain_id)) {
      chain_id <- all_samples[[i]]$chain_info$chain_id[[1]]
    }

    sample_i <- .ensure_sample_chain_metadata(
      all_samples[[i]],
      fallback_chain_id = chain_id,
      main_offset = main_offset,
      theta_offset = theta_offset
    )
    meta <- .sample_chain_metadata(
      sample_i,
      chain_id = chain_id,
      seed = sample_i$mcmc_control$seed,
      main_offset = main_offset,
      theta_offset = theta_offset
    )

    chain_info[[i]] <- meta$chain_info
    main_index[[i]] <- meta$draw_index$main
    theta_index[[i]] <- meta$draw_index$theta
    run_history[[i]] <- meta$run_history
    run_history[[i]]$run_id <- i

    main_offset <- main_offset + meta$chain_info$n_draws_main
    theta_offset <- theta_offset + meta$chain_info$n_draws_theta
  }

  list(
    chain_info = do.call(rbind, chain_info),
    draw_index = list(
      main = do.call(rbind, main_index),
      theta = do.call(rbind, theta_index)
    ),
    run_history = do.call(rbind, run_history)
  )
}


.build_chain_diagnostics <- function(ess,
                                     waic,
                                     loglik_trace,
                                     cluster_info,
                                     chain_info,
                                     draw_index,
                                     draws = list()) {
  waic_by_chain <- if (!is.null(chain_info)) {
    chain_info[, c("chain_id", "waic"), drop = FALSE]
  } else {
    data.frame(chain_id = integer(0), waic = numeric(0))
  }

  loglik_by_chain <- NULL
  if (!is.null(loglik_trace) && !is.null(draw_index$main)) {
    n <- min(length(loglik_trace), nrow(draw_index$main))
    loglik_by_chain <- data.frame(
      chain_id = draw_index$main$chain_id[seq_len(n)],
      row = draw_index$main$row[seq_len(n)],
      within_chain_draw = draw_index$main$within_chain_draw[seq_len(n)],
      mcmc_iteration = draw_index$main$mcmc_iteration[seq_len(n)],
      value = as.numeric(loglik_trace[seq_len(n)])
    )
  }

  cluster_by_chain <- NULL
  if (!is.null(cluster_info$n_clusters) && !is.null(draw_index$main)) {
    n <- min(length(cluster_info$n_clusters), nrow(draw_index$main))
    cluster_by_chain <- data.frame(
      chain_id = draw_index$main$chain_id[seq_len(n)],
      row = draw_index$main$row[seq_len(n)],
      within_chain_draw = draw_index$main$within_chain_draw[seq_len(n)],
      mcmc_iteration = draw_index$main$mcmc_iteration[seq_len(n)],
      n_clusters = as.integer(cluster_info$n_clusters[seq_len(n)])
    )
  }

  ess_by_chain <- .compute_ess_by_chain(draws, draw_index)
  rhat <- .compute_rhat_by_family(draws, draw_index)

  list(
    ess = list(
      pooled = ess,
      by_chain = ess_by_chain
    ),
    rhat = rhat,
    waic = list(
      value = waic,
      by_chain = waic_by_chain,
      aggregation = if (!is.null(chain_info) && nrow(chain_info) > 1L) {
        "mean_of_chain_waic"
      } else {
        "single_chain"
      }
    ),
    loglik = list(
      trace = loglik_trace,
      by_chain = loglik_by_chain
    ),
    clusters = list(
      trace = cluster_info$n_clusters,
      by_chain = cluster_by_chain
    )
  )
}


.draw_index_for_family <- function(family, draw_index) {
  if (is.null(draw_index)) {
    return(NULL)
  }
  if (identical(family, "theta")) {
    draw_index$theta
  } else {
    draw_index$main
  }
}


.draw_param_names <- function(mat) {
  cn <- colnames(mat)
  if (!is.null(cn)) {
    return(cn)
  }
  paste0("V", seq_len(ncol(mat)))
}


.split_matrix_by_chain <- function(mat, index) {
  if (is.null(mat) || is.null(index) || nrow(mat) != nrow(index) ||
      nrow(mat) == 0L || ncol(mat) == 0L) {
    return(NULL)
  }

  split(seq_len(nrow(index)), index$chain_id)
}


.compute_ess_by_chain <- function(draws, draw_index) {
  result <- list()

  for (family in names(draws)) {
    mat <- draws[[family]]
    index <- .draw_index_for_family(family, draw_index)
    split_rows <- .split_matrix_by_chain(mat, index)
    if (is.null(split_rows)) {
      next
    }

    family_rows <- lapply(names(split_rows), function(chain_id) {
      rows <- split_rows[[chain_id]]
      ess <- tryCatch(
        coda::effectiveSize(coda::as.mcmc(mat[rows, , drop = FALSE])),
        error = function(e) rep(NA_real_, ncol(mat))
      )
      data.frame(
        chain_id = as.integer(chain_id),
        parameter = .draw_param_names(mat),
        ess = as.numeric(ess),
        stringsAsFactors = FALSE
      )
    })
    result[[family]] <- do.call(rbind, family_rows)
  }

  if (length(result) == 0L) {
    NULL
  } else {
    result
  }
}


.compute_rhat_by_family <- function(draws, draw_index) {
  result <- list()

  for (family in names(draws)) {
    mat <- draws[[family]]
    index <- .draw_index_for_family(family, draw_index)
    rhat <- .compute_rhat_matrix(mat, index)
    if (!is.null(rhat)) {
      result[[family]] <- rhat
    }
  }

  if (length(result) == 0L) {
    NULL
  } else {
    result
  }
}


.compute_rhat_matrix <- function(mat, index) {
  split_rows <- .split_matrix_by_chain(mat, index)
  if (is.null(split_rows) || length(split_rows) < 2L) {
    return(NULL)
  }

  split_rows <- split_rows[vapply(split_rows, length, integer(1)) >= 2L]
  if (length(split_rows) < 2L) {
    out <- rep(NA_real_, ncol(mat))
    names(out) <- .draw_param_names(mat)
    return(out)
  }

  out <- tryCatch({
    mcmc_list <- coda::mcmc.list(lapply(split_rows, function(rows) {
      coda::as.mcmc(mat[rows, , drop = FALSE])
    }))
    psrf <- coda::gelman.diag(
      mcmc_list,
      autoburnin = FALSE,
      multivariate = FALSE
    )$psrf
    vals <- psrf[, "Point est."]
    names(vals) <- .draw_param_names(mat)
    vals
  }, error = function(e) {
    vals <- rep(NA_real_, ncol(mat))
    names(vals) <- .draw_param_names(mat)
    vals
  })

  out
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
