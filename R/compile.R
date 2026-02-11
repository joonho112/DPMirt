# ============================================================================
# Module 2: Compilation
# ============================================================================
# Purpose: Compile NIMBLE model and MCMC, separating this expensive step
#          from fast re-sampling.
# Blueprint: Section 6 - Module 2, Section 7.3
# ============================================================================

#' Compile a DPMirt model specification
#'
#' Takes a \code{dpmirt_spec} object and compiles the NIMBLE model and MCMC
#' engine. This is the expensive step (~30-120 seconds) that only needs to be
#' done once per model specification.
#'
#' @param spec A \code{dpmirt_spec} object from \code{\link{dpmirt_spec}}.
#' @param sampler_config Optional custom MCMC sampler configuration.
#' @param use_centered_sampler Character or logical. Whether to use the
#'   centered sampler for SI parameterization. "auto" enables it when
#'   appropriate (SI param + 2PL/3PL). Only relevant for Phase 3+.
#' @param enable_waic Logical. Whether to enable WAIC computation.
#' @param enable_logprob_monitor Logical. Whether to add log-probability
#'   monitoring samplers.
#' @param verbose Logical. Print progress messages.
#' @param ... Additional arguments passed to NIMBLE.
#'
#' @return A \code{dpmirt_compiled} S3 object containing the compiled model,
#'   compiled MCMC, specification reference, and session signature.
#'
#' @details
#' The compiled NIMBLE objects contain external C++ pointers that
#' \strong{cannot be serialized} across R sessions. The compile-once pattern
#' is therefore a within-session optimization. For cross-session workflows,
#' save the \code{dpmirt_spec} object and recompile in the new session.
#'
#' @examples
#' \dontrun{
#' sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
#' spec <- dpmirt_spec(sim$response, model = "rasch", prior = "normal")
#'
#' # Compile (takes 30-120 seconds)
#' compiled <- dpmirt_compile(spec)
#' print(compiled)
#'
#' # Compile-once, sample-many pattern
#' samples1 <- dpmirt_sample(compiled, niter = 5000, nburnin = 1000, seed = 1)
#' samples2 <- dpmirt_sample(compiled, niter = 5000, nburnin = 1000, seed = 2)
#' }
#'
#' @family model fitting
#' @seealso \code{\link{dpmirt_spec}}, \code{\link{dpmirt_sample}},
#'   \code{\link{dpmirt}}
#'
#' @export
dpmirt_compile <- function(spec,
                           sampler_config = NULL,
                           use_centered_sampler = "auto",
                           enable_waic = TRUE,
                           enable_logprob_monitor = TRUE,
                           verbose = TRUE,
                           ...) {

  # --- Validate input ---
  if (!inherits(spec, "dpmirt_spec")) {
    stop("Input must be a dpmirt_spec object from dpmirt_spec().",
         call. = FALSE)
  }

  timer_start <- .start_timer()

  # --- Step 0: Register custom distributions if needed ---
  if (spec$config$identification == "constrained_item" &&
      spec$config$model %in% c("2pl", "3pl")) {
    .vmsg("Registering dBernoulliVector...", verbose = verbose)
    .register_dBernoulliVector()
  }

  # --- Step 1: Build NIMBLE model ---
  .vmsg("Building NIMBLE model...", verbose = verbose)

  nimble_model <- nimbleModel(
    code      = spec$code,
    constants = spec$constants,
    data      = spec$data,
    inits     = spec$inits,
    name      = paste0("DPMirt_", spec$config$model, "_", spec$config$prior)
  )

  .vmsg("  Model built successfully.", verbose = verbose)

  # --- Step 2: Configure MCMC ---
  .vmsg("Configuring MCMC...", verbose = verbose)

  mcmc_conf <- .configure_mcmc(
    nimble_model     = nimble_model,
    spec             = spec,
    sampler_config   = sampler_config,
    use_centered_sampler = use_centered_sampler,
    enable_waic      = enable_waic,
    enable_logprob_monitor = enable_logprob_monitor,
    verbose          = verbose
  )

  .vmsg("  MCMC configured.", verbose = verbose)

  # --- Step 3: Build MCMC ---
  nimble_mcmc <- buildMCMC(mcmc_conf)
  .vmsg("  MCMC built.", verbose = verbose)

  # --- Step 4: Compile model and MCMC ---
  .vmsg("Compiling model (this may take 30-120 seconds)...", verbose = verbose)

  Cmodel <- compileNimble(nimble_model)
  Cmcmc  <- compileNimble(nimble_mcmc, project = nimble_model)

  compilation_time <- .elapsed_time(timer_start)
  .vmsg("  Compilation complete in ", .format_time(compilation_time),
        verbose = verbose)

  # --- Step 5: Build session signature ---
  sig <- .session_signature(spec)

  # --- Construct and return compiled object ---
  compiled <- structure(
    list(
      Cmodel           = Cmodel,
      Cmcmc            = Cmcmc,
      spec             = spec,
      mcmc_conf        = mcmc_conf,
      compilation_time = compilation_time,
      session_signature = sig,
      nimble_version   = as.character(packageVersion("nimble"))
    ),
    class = "dpmirt_compiled"
  )

  compiled
}


# ============================================================================
# MCMC Configuration
# ============================================================================

#' Configure MCMC with custom samplers and monitors
#'
#' Internal function that sets up:
#' 1. Standard MCMC configuration with monitors
#' 2. logProb_summer samplers for probability monitoring
#' 3. Centered samplers for SI parameterization (Phase 3+)
#' 4. WAIC computation settings
#' 5. Differential thinning for person parameters
#'
#' @noRd
.configure_mcmc <- function(nimble_model, spec, sampler_config,
                             use_centered_sampler, enable_waic,
                             enable_logprob_monitor, verbose) {

  # --- Basic configuration with monitors ---
  conf <- configureMCMC(
    nimble_model,
    monitors  = spec$monitors,
    enableWAIC = enable_waic
  )

  # --- Add monitors2 for thinned parameters ---
  for (m2 in spec$monitors2) {
    conf$addMonitors2(m2)
  }

  # --- Log-probability monitoring (Paganin's logProb_summer pattern) ---
  if (enable_logprob_monitor) {
    .add_logprob_monitors(conf, nimble_model, spec)
  }

  # --- Centered sampler for SI parameterization ---
  .add_centered_sampler(conf, nimble_model, spec, use_centered_sampler,
                         verbose)

  # --- WAIC grouping (by person) ---
  if (enable_waic) {
    # NIMBLE's WAIC uses conditional WAIC by default
    # For IRT models, group by person
    # This is automatically handled by NIMBLE when data are
    # structured as y[j, i] with j = person index
  }

  conf
}


#' Add log-probability monitoring samplers
#'
#' Replaces the default samplers for the dummy monitoring nodes
#' (myLogProbAll, myLogProbSome, myLogLik) with custom logProb_summer
#' samplers that track log-probabilities during MCMC.
#'
#' Adapted from Paganin's monitorLogProb.R
#'
#' @noRd
.add_logprob_monitors <- function(conf, nimble_model, spec) {

  # Remove default samplers for monitoring nodes
  conf$removeSampler("myLogProbAll")
  conf$removeSampler("myLogProbSome")
  conf$removeSampler("myLogLik")

  # Get the logProb_summer sampler via lazy init (R/samplers.R).
  lps <- .get_logProb_summer()

  # Add logProb_summer for all parameters
  conf$addSampler(
    target  = "myLogProbAll",
    type    = lps
  )

  # Build "some_nodes" list based on model type
  # Paganin: nodeList = c("beta", "lambda", "eta") for IRT
  #          nodeList = c("gamma", "lambda", "eta") for SI
  is_2pl_or_3pl <- spec$config$model %in% c("2pl", "3pl")
  is_si <- !is.null(spec$config$parameterization) &&
    spec$config$parameterization == "si"

  if (is_si && is_2pl_or_3pl) {
    some_nodes <- c("gamma", "lambda", "eta")
  } else if (is_2pl_or_3pl) {
    some_nodes <- c("beta", "lambda", "eta")
  } else {
    some_nodes <- c("beta", "eta")
  }

  conf$addSampler(
    target  = "myLogProbSome",
    type    = lps,
    control = list(nodeList = some_nodes)
  )

  # Add logProb_summer for data likelihood
  conf$addSampler(
    target  = "myLogLik",
    type    = lps,
    control = list(nodeList = "y")
  )
}


#' Add centered sampler for SI parameterization
#'
#' For 2PL/3PL SI parameterization with unconstrained or constrained_ability
#' identification, replaces default samplers for log_lambda with the centered
#' sampler that jointly proposes (log_lambda\[i\], gamma\[i\]).
#'
#' @noRd
.add_centered_sampler <- function(conf, nimble_model, spec,
                                   use_centered_sampler, verbose) {

  model_type <- spec$config$model
  param_type <- spec$config$parameterization
  id_type    <- spec$config$identification

  # Auto-resolve: only use centered for SI + 2PL/3PL + (unconstrained|constrained_ability)
  if (identical(use_centered_sampler, "auto")) {
    use_it <- (model_type %in% c("2pl", "3pl")) &&
      (param_type == "si") &&
      (id_type %in% c("unconstrained", "constrained_ability"))
  } else {
    use_it <- isTRUE(use_centered_sampler)
  }

  if (!use_it) return(invisible(NULL))

  .vmsg("  Adding centered sampler for SI parameterization...",
        verbose = verbose)

  # Get the centered sampler components
  centered <- .get_centered_sampler()

  # Determine target node names
  I <- spec$config$I

  # Remove default samplers for log_lambda (or logLambda.tmp)
  # and add centered sampler
  if (id_type == "constrained_item") {
    # For constrained_item, nodes are logLambda.tmp and gamma.tmp
    # Centered sampler pairs: (logLambda.tmp[i], gamma.tmp[i])
    for (i in seq_len(I)) {
      conf$removeSamplers(paste0("logLambda.tmp[", i, "]"))
    }
    conf$addSampler(
      type    = centered$multi,
      target  = c("logLambda.tmp", "gamma.tmp"),
      control = list(nodesToCenter = "eta", scale = 0.1, adaptive = TRUE)
    )
  } else {
    # For unconstrained/constrained_ability, nodes are log_lambda and gamma
    # But NIMBLE uses: log(lambda[i]) ~ dnorm(...) → the stochastic node
    # is named "log_lambda[i]" by NIMBLE internally
    # We need to remove the sampler on "log_lambda"
    for (i in seq_len(I)) {
      tryCatch(
        conf$removeSamplers(paste0("log_lambda[", i, "]")),
        error = function(e) {
          # NIMBLE might name it differently; try the lifted node
          tryCatch(
            conf$removeSamplers(paste0("lifted_log_lambda_oBi_cB_L", i)),
            error = function(e2) NULL
          )
        }
      )
    }
    conf$addSampler(
      type    = centered$multi,
      target  = c("log_lambda", "gamma"),
      control = list(nodesToCenter = "eta", scale = 0.1, adaptive = TRUE)
    )
  }

  invisible(NULL)
}


# ============================================================================
# S3 Methods for dpmirt_compiled
# ============================================================================

#' @rdname dpmirt_compile
#' @param x A \code{dpmirt_compiled} object.
#' @param ... Additional arguments (currently unused).
#' @export
print.dpmirt_compiled <- function(x, ...) {
  cat("DPMirt Compiled Model\n")
  cat("=====================\n")
  cat("Model:           ", toupper(x$spec$config$model), "\n")
  cat("Prior:           ", x$spec$config$prior, "\n")
  cat("Identification:  ", x$spec$config$identification, "\n")
  cat("Persons (N):     ", x$spec$config$N, "\n")
  cat("Items (I):       ", x$spec$config$I, "\n")
  cat("Compilation time:", .format_time(x$compilation_time), "\n")
  cat("NIMBLE version:  ", x$nimble_version, "\n")
  cat("\nNOTE: This object contains compiled C++ pointers and\n")
  cat("      cannot be saved/loaded across R sessions.\n")
  invisible(x)
}
