# ============================================================================
# Main Entry Point: dpmirt()
# ============================================================================
# Purpose: One-step model fitting that orchestrates the full pipeline:
#          spec -> compile -> sample -> rescale -> fit object
# Blueprint: Section 5.1
# ============================================================================

#' Fit a Bayesian IRT model with DPMirt
#'
#' The main entry point for DPMirt. Handles the full pipeline from data
#' to fitted model: specification, compilation, MCMC sampling, rescaling,
#' and diagnostics.
#'
#' @param data A matrix or data.frame of binary (0/1) responses. Persons in
#'   rows, items in columns.
#' @param model Character. IRT model type: \code{"rasch"}, \code{"2pl"},
#'   or \code{"3pl"}.
#' @param prior Character. Latent trait prior: \code{"normal"} or \code{"dpm"}.
#' @param parameterization Character. \code{"irt"} or \code{"si"}. Default
#'   \code{"irt"}.
#' @param identification Character or NULL. Identification strategy. If NULL,
#'   uses model-specific default.
#' @param niter Integer. Total MCMC iterations. Default 10000.
#' @param nburnin Integer. Burn-in iterations. Default 2000.
#' @param thin Integer. Thinning for main monitors. Default 1.
#' @param thin2 Integer or NULL. Thinning for eta. NULL = same as thin.
#' @param nchains Integer. Number of chains. Default 1.
#' @param seed Integer or NULL. Random seed.
#' @param alpha_prior Alpha hyperprior for DPM. NULL (default Gamma(1,3) or
#'   auto-elicit if \code{mu_K} is set), a numeric vector \code{c(a, b)} for
#'   Gamma(a, b), or the output of \code{\link{dpmirt_alpha_prior}}.
#' @param mu_K Numeric or NULL. Expected number of clusters for automatic
#'   alpha prior elicitation via \code{\link{dpmirt_alpha_prior}}. Only used
#'   when \code{prior = "dpm"} and \code{alpha_prior = NULL}. Requires
#'   the DPprior package; falls back to Gamma(1, 3) if not installed.
#' @param confidence Character. DPprior confidence level for alpha elicitation:
#'   \code{"low"}, \code{"medium"} (default), or \code{"high"}. Only used when
#'   \code{mu_K} is set.
#' @param base_measure List of DPM base measure hyperparameters.
#' @param M Integer. Max clusters for CRP. Default 50.
#' @param item_priors List of custom item priors.
#' @param rescale Logical. Apply post-hoc rescaling. Default TRUE.
#' @param compute_waic Logical. Compute WAIC. Default TRUE.
#' @param compute_dp_density Logical. Compute DP density. Default TRUE for DPM.
#' @param save_draws Logical. Save full theta posterior. Default TRUE.
#' @param save_path Character or NULL. Path for disk-backed storage.
#' @param sampler_config Optional custom sampler configuration.
#' @param verbose Logical. Print progress. Default TRUE.
#' @param ... Additional arguments.
#'
#' @return A \code{dpmirt_fit} S3 object.
#'
#' @examples
#' \dontrun{
#' # --- Simulate test data ---
#' sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
#'
#' # --- Rasch model with Normal prior ---
#' fit_normal <- dpmirt(sim$response, model = "rasch", prior = "normal",
#'                      niter = 5000, nburnin = 1000, seed = 123)
#' summary(fit_normal)
#' coef(fit_normal)
#'
#' # --- Rasch model with DPM prior ---
#' fit_dpm <- dpmirt(sim$response, model = "rasch", prior = "dpm",
#'                   niter = 10000, nburnin = 3000, seed = 123)
#'
#' # --- Compare Normal vs DPM ---
#' dpmirt_compare(fit_normal, fit_dpm)
#' }
#'
#' @references
#' Paganin, S., Paciorek, C. J., Wehrhahn, C., Rodriguez, A.,
#' Rabe-Hesketh, S., & de Valpine, P. (2023). Computational strategies
#' and estimation performance with Bayesian semiparametric item response
#' theory models. \emph{Journal of Educational and Behavioral Statistics,
#' 48}(2), 147--188.
#'
#' @family model fitting
#' @seealso \code{\link{dpmirt_spec}}, \code{\link{dpmirt_compile}},
#'   \code{\link{dpmirt_sample}}, \code{\link{dpmirt_estimates}},
#'   \code{\link{dpmirt_simulate}}
#'
#' @export
dpmirt <- function(data,
                   model = c("rasch", "2pl", "3pl"),
                   prior = c("normal", "dpm"),
                   parameterization = c("irt", "si"),
                   identification = NULL,

                   # MCMC settings
                   niter = 10000L,
                   nburnin = 2000L,
                   thin = 1L,
                   thin2 = NULL,
                   nchains = 1L,
                   seed = NULL,

                   # DPM settings
                   alpha_prior = NULL,
                   mu_K = NULL,
                   confidence = "medium",
                   base_measure = list(s2_mu = 2, nu1 = 2.01, nu2 = 1.01),
                   M = 50L,

                   # Item prior settings
                   item_priors = list(),

                   # Output control
                   rescale = TRUE,
                   compute_waic = TRUE,
                   compute_dp_density = TRUE,
                   save_draws = TRUE,
                   save_path = NULL,

                   # Advanced
                   sampler_config = NULL,
                   verbose = TRUE,
                   ...) {

  total_timer <- .start_timer()

  model <- match.arg(model)
  prior <- match.arg(prior)
  parameterization <- match.arg(parameterization)

  .vmsg("=== DPMirt: Fitting ", toupper(model), " model with ",
        prior, " prior ===", verbose = verbose)

  # --- Step 0.5: Auto-elicit alpha prior if mu_K is provided ---
  if (prior == "dpm" && is.null(alpha_prior) && !is.null(mu_K)) {
    .vmsg("  Eliciting alpha prior via DPprior (mu_K = ", mu_K,
          ", confidence = \"", confidence, "\")...", verbose = verbose)
    alpha_prior <- dpmirt_alpha_prior(
      N          = nrow(as.matrix(data)),
      mu_K       = mu_K,
      confidence = confidence
    )
  }

  # --- Step 1: Model Specification ---
  .vmsg("\n[1/4] Creating model specification...", verbose = verbose)
  spec <- dpmirt_spec(
    data             = data,
    model            = model,
    prior            = prior,
    parameterization = parameterization,
    identification   = identification,
    alpha_prior      = alpha_prior,
    base_measure     = base_measure,
    item_priors      = item_priors,
    M                = M
  )

  .vmsg("  ", spec$config$N, " persons x ", spec$config$I, " items",
        verbose = verbose)
  .vmsg("  Identification: ", spec$config$identification, verbose = verbose)

  # --- Step 2: Compilation ---
  .vmsg("\n[2/4] Compiling NIMBLE model...", verbose = verbose)
  compiled <- dpmirt_compile(
    spec                   = spec,
    sampler_config         = sampler_config,
    enable_waic            = compute_waic,
    enable_logprob_monitor = TRUE,
    verbose                = verbose
  )

  # --- Step 3: Sampling ---
  .vmsg("\n[3/4] Running MCMC...", verbose = verbose)

  if (nchains == 1) {
    samples_obj <- dpmirt_sample(
      compiled = compiled,
      niter    = niter,
      nburnin  = nburnin,
      thin     = thin,
      thin2    = thin2,
      seed     = seed,
      verbose  = verbose
    )
    all_samples <- list(samples_obj)
  } else {
    # Multiple chains
    all_samples <- vector("list", nchains)
    for (ch in seq_len(nchains)) {
      ch_seed <- if (!is.null(seed)) seed + ch - 1 else NULL
      .vmsg("  Chain ", ch, "/", nchains, "...", verbose = verbose)
      all_samples[[ch]] <- dpmirt_sample(
        compiled = compiled,
        niter    = niter,
        nburnin  = nburnin,
        thin     = thin,
        thin2    = thin2,
        seed     = ch_seed,
        verbose  = FALSE
      )
    }
  }

  # --- Step 4: Rescale and build fit object ---
  .vmsg("\n[4/4] Post-processing (rescaling, diagnostics)...",
        verbose = verbose)

  # Combine chains: rescale each chain, then row-bind
  if (nchains == 1) {
    primary_samples <- all_samples[[1]]
    rescaled <- dpmirt_rescale(primary_samples, rescale = rescale)
    waic_val <- primary_samples$waic
    loglik_trace <- .extract_loglik_trace(primary_samples$samples)
    raw_samples <- primary_samples$samples
    sampling_time_total <- primary_samples$sampling_time
  } else {
    .vmsg("  Combining ", nchains, " chains...", verbose = verbose)
    combined <- .combine_chains(all_samples, rescale = rescale)
    rescaled        <- combined$rescaled
    waic_val        <- combined$waic
    loglik_trace    <- combined$loglik_trace
    raw_samples     <- combined$raw_samples
    sampling_time_total <- combined$sampling_time
  }

  # --- Compute diagnostics ---
  ess <- .compute_ess(rescaled)

  # --- DPM-specific: cluster diagnostics ---
  cluster_info <- NULL
  if (prior == "dpm" && !is.null(raw_samples)) {
    .vmsg("  Computing cluster diagnostics...", verbose = verbose)
    cluster_info <- .extract_cluster_info(raw_samples, spec$config$N)
  }

  # --- Build fit object ---
  total_time <- .elapsed_time(total_timer)

  fit <- structure(
    list(
      # Core results
      theta_samp     = rescaled$theta_samp,
      beta_samp      = rescaled$beta_samp,
      lambda_samp    = rescaled$lambda_samp,
      delta_samp     = rescaled$delta_samp,
      other_samp     = rescaled$other_samp,

      # Rescaling info
      scale_shift    = rescaled$scale_shift,
      location_shift = rescaled$location_shift,

      # DPM-specific
      dp_density     = NULL,   # Computed below if requested
      cluster_info   = cluster_info,

      # Model comparison
      waic           = waic_val,

      # Diagnostics
      ess            = ess,
      loglik_trace   = loglik_trace,

      # Timings
      compilation_time = compiled$compilation_time,
      sampling_time    = sampling_time_total,
      total_time       = total_time,

      # Compiled reference (for resume)
      compiled       = compiled,

      # Configuration
      config         = list(
        model            = spec$config$model,
        prior            = spec$config$prior,
        parameterization = spec$config$parameterization,
        identification   = spec$config$identification,
        niter            = niter,
        nburnin          = nburnin,
        thin             = thin,
        thin2            = thin2,
        nchains          = nchains,
        alpha_prior      = spec$config$alpha_prior,
        base_measure     = base_measure,
        M                = M,
        N                = spec$config$N,
        I                = spec$config$I,
        seed             = seed,
        rescale          = rescale
      )
    ),
    class = "dpmirt_fit"
  )

  # --- DPM-specific: DP density (optional, can be expensive) ---
  if (prior == "dpm" && compute_dp_density) {
    .vmsg("  Computing DP density...", verbose = verbose)
    fit$dp_density <- tryCatch(
      dpmirt_dp_density(fit, verbose = verbose),
      error = function(e) {
        .vmsg("  Note: DP density computation failed: ", conditionMessage(e),
              verbose = verbose)
        .vmsg("  You can compute it manually with dpmirt_dp_density(fit).",
              verbose = verbose)
        NULL
      }
    )
  }

  # Update total time to include DP density computation
  fit$total_time <- .elapsed_time(total_timer)

  .vmsg("\n=== Done in ", .format_time(fit$total_time), " ===",
        verbose = verbose)

  fit
}


# ============================================================================
# Internal Helpers for dpmirt()
# ============================================================================

#' Compute effective sample sizes for key parameters
#' @noRd
.compute_ess <- function(rescaled) {
  ess <- list()

  if (!is.null(rescaled$beta_samp)) {
    ess$items <- tryCatch(
      coda::effectiveSize(coda::as.mcmc(rescaled$beta_samp)),
      error = function(e) rep(NA_real_, ncol(rescaled$beta_samp))
    )
  }

  if (!is.null(rescaled$lambda_samp)) {
    ess$lambda <- tryCatch(
      coda::effectiveSize(coda::as.mcmc(rescaled$lambda_samp)),
      error = function(e) rep(NA_real_, ncol(rescaled$lambda_samp))
    )
  }

  if (!is.null(rescaled$delta_samp)) {
    ess$delta <- tryCatch(
      coda::effectiveSize(coda::as.mcmc(rescaled$delta_samp)),
      error = function(e) rep(NA_real_, ncol(rescaled$delta_samp))
    )
  }

  if (!is.null(rescaled$theta_samp)) {
    ess$theta <- tryCatch(
      coda::effectiveSize(coda::as.mcmc(rescaled$theta_samp)),
      error = function(e) rep(NA_real_, ncol(rescaled$theta_samp))
    )
  }

  ess
}


#' Extract log-likelihood trace from samples
#' @noRd
.extract_loglik_trace <- function(samples) {
  ll_col <- grep("^myLogLik$", colnames(samples))
  if (length(ll_col) == 1) {
    return(samples[, ll_col])
  }
  NULL
}


#' Combine multiple MCMC chains
#'
#' Rescales each chain independently, then row-binds the rescaled samples.
#' Independent rescaling is required because the posterior from each chain
#' may have different centering before identification constraints are applied.
#'
#' @param all_samples List of \code{dpmirt_samples} objects.
#' @param rescale Logical. Apply rescaling.
#' @return A list with combined rescaled results and aggregated diagnostics.
#' @noRd
.combine_chains <- function(all_samples, rescale = TRUE) {

  nchains <- length(all_samples)

  # Rescale each chain independently
  rescaled_list <- lapply(all_samples, dpmirt_rescale, rescale = rescale)

  # Row-bind rescaled sample matrices
  combined_theta  <- do.call(rbind, lapply(rescaled_list, `[[`, "theta_samp"))
  combined_beta   <- do.call(rbind, lapply(rescaled_list, `[[`, "beta_samp"))

  combined_lambda <- NULL
  if (!is.null(rescaled_list[[1]]$lambda_samp)) {
    combined_lambda <- do.call(rbind,
                               lapply(rescaled_list, `[[`, "lambda_samp"))
  }

  combined_delta <- NULL
  if (!is.null(rescaled_list[[1]]$delta_samp)) {
    combined_delta <- do.call(rbind,
                              lapply(rescaled_list, `[[`, "delta_samp"))
  }

  combined_other <- NULL
  if (!is.null(rescaled_list[[1]]$other_samp)) {
    combined_other <- do.call(rbind,
                              lapply(rescaled_list, `[[`, "other_samp"))
  }

  combined_scale    <- unlist(lapply(rescaled_list, `[[`, "scale_shift"))
  combined_location <- unlist(lapply(rescaled_list, `[[`, "location_shift"))

  # Aggregate WAIC: use mean across chains (each chain provides one WAIC)
  waic_vals <- sapply(all_samples, function(s) s$waic)
  waic_val  <- if (all(is.null(waic_vals))) NULL else mean(unlist(waic_vals), na.rm = TRUE)

  # Concatenate log-likelihood traces
  loglik_traces <- lapply(all_samples, function(s) {
    .extract_loglik_trace(s$samples)
  })
  loglik_trace <- if (all(sapply(loglik_traces, is.null))) {
    NULL
  } else {
    unlist(loglik_traces)
  }

  # Row-bind raw samples for cluster diagnostics
  raw_samples <- do.call(rbind, lapply(all_samples, `[[`, "samples"))

  # Total sampling time across chains
  sampling_time <- sum(sapply(all_samples, function(s) s$sampling_time))

  list(
    rescaled = list(
      theta_samp     = combined_theta,
      beta_samp      = combined_beta,
      lambda_samp    = combined_lambda,
      delta_samp     = combined_delta,
      other_samp     = combined_other,
      scale_shift    = combined_scale,
      location_shift = combined_location,
      samples_raw    = raw_samples,
      samples2_raw   = do.call(rbind,
                               lapply(rescaled_list, `[[`, "samples2_raw"))
    ),
    waic          = waic_val,
    loglik_trace  = loglik_trace,
    raw_samples   = raw_samples,
    sampling_time = sampling_time
  )
}
