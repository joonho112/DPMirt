# ============================================================================
# Module 7: DP Measure Sampling & Density Computation
# ============================================================================
# Purpose: Sample from the posterior predictive distribution under the DPM
#          prior and compute the posterior density on a grid.
# Blueprint: Section 6 - Module 7
#
# Adapted from Paganin et al. (2023):
#   - 3_simulateFromDPmeasure.R  (getSamplesDPmeasure workflow)
#   - 4_computeQuantitesForFigures.R  (density evaluation from DP samples)
#
# Workflow:
# 1. Extract DP-related MCMC samples (alpha, zi, muTilde, s2Tilde)
# 2. Reconstruct model + compiled MCMC, populate with posterior samples
# 3. Call NIMBLE's getSamplesDPmeasure() to get stick-breaking weights + atoms
# 4. Evaluate the finite mixture density on a grid
# 5. Apply rescaling (location shift for Rasch) to the grid
# 6. Return density with pointwise credible intervals
# ============================================================================


#' Compute posterior density of the DP mixture
#'
#' Samples from the posterior Dirichlet Process mixing distribution using
#' NIMBLE's \code{getSamplesDPmeasure()} and evaluates the resulting
#' mixture density on a grid. The density is computed by summing
#' weighted Normal components from the DP base measure.
#'
#' @param fit A \code{dpmirt_fit} object with \code{prior = "dpm"}.
#' @param grid Numeric vector. Grid points for density evaluation.
#'   Default: \code{seq(-6, 6, length.out = 500)}.
#' @param credible_interval Numeric. Width of the pointwise credible band.
#'   Default: 0.95 (i.e., 95 percent band).
#' @param apply_rescaling Logical. If TRUE (default), shift the grid by the
#'   iteration-specific location shift from post-hoc rescaling. Only
#'   relevant for unconstrained models.
#' @param verbose Logical. Print progress messages. Default TRUE.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list of class \code{dpmirt_dp_density} containing:
#' \describe{
#'   \item{grid}{Numeric vector of evaluation points.}
#'   \item{density_mean}{Numeric vector of posterior mean densities.}
#'   \item{density_lower}{Numeric vector of lower credible band.}
#'   \item{density_upper}{Numeric vector of upper credible band.}
#'   \item{density_samples}{Matrix (niter x length(grid)) of per-iteration
#'     densities (for custom summaries).}
#'   \item{dp_samples}{List from \code{getSamplesDPmeasure()} -- each element
#'     is a matrix with columns (weights, means, variances).}
#'   \item{ci_level}{The credible interval level used.}
#' }
#'
#' @details
#' This function follows Paganin et al.'s (2023) workflow:
#'
#' \enumerate{
#'   \item Extract posterior samples for DP parameters (alpha, zi, muTilde,
#'     s2Tilde) from the fitted model.
#'   \item Reconstruct a NIMBLE model and compiled MCMC with monitors set
#'     to only DP parameters.
#'   \item Populate the compiled MCMC's sample storage with the posterior
#'     samples using \code{nimble:::matrix2mv()}.
#'   \item Call \code{getSamplesDPmeasure()} to compute stick-breaking
#'     weights and atoms for each posterior draw.
#'   \item Evaluate the mixture density
#'     `f(x|Gs) = sum_k w_k * phi(x; mu_k, s2_k)`
#'     for each posterior sample s and grid point x.
#' }
#'
#' For Rasch models with unconstrained identification, a location shift
#' (mean(beta) per iteration) is applied so the density is on the
#' rescaled theta scale.
#'
#' @examples
#' \dontrun{
#' sim <- dpmirt_simulate(200, 20, model = "rasch",
#'                        latent_shape = "bimodal", seed = 42)
#' fit <- dpmirt(sim$response, model = "rasch", prior = "dpm",
#'               niter = 10000, nburnin = 3000, seed = 123)
#'
#' # Compute DP density on default grid
#' dpd <- dpmirt_dp_density(fit)
#' print(dpd)
#'
#' # Custom grid
#' dpd2 <- dpmirt_dp_density(fit, grid = seq(-4, 4, length.out = 200))
#' }
#'
#' @family DP density
#' @seealso \code{\link{dpmirt}}, \code{\link{dpmirt_plot_dp_density}}
#'
#' @export
dpmirt_dp_density <- function(fit,
                              grid = seq(-6, 6, length.out = 500),
                              credible_interval = 0.95,
                              apply_rescaling = TRUE,
                              verbose = TRUE,
                              ...) {

  # --- Validate input ---
  if (!inherits(fit, "dpmirt_fit")) {
    stop("Input must be a dpmirt_fit object.", call. = FALSE)
  }

  if (fit$config$prior != "dpm") {
    stop("DP density is only available for DPM models.", call. = FALSE)
  }

  .vmsg("Computing DP density...", verbose = verbose)
  timer_start <- .start_timer()

  # --- Step 1: Extract DP-related posterior samples ---
  .vmsg("  Extracting DP posterior samples...", verbose = verbose)

  dp_samples_raw <- .extract_dp_samples(fit)

  # --- Step 2: Call getSamplesDPmeasure via NIMBLE ---
  .vmsg("  Sampling from DP measure (getSamplesDPmeasure)...",
        verbose = verbose)

  dp_measure <- .get_dp_measure_samples(fit, dp_samples_raw, verbose)

  # --- Step 3: Evaluate density on grid ---
  .vmsg("  Evaluating density on grid (", length(grid), " points)...",
        verbose = verbose)

  density_mat <- .evaluate_dp_density(
    dp_measure   = dp_measure,
    grid         = grid,
    fit          = fit,
    apply_rescaling = apply_rescaling
  )

  # --- Step 4: Compute summary statistics ---
  alpha_lower <- (1 - credible_interval) / 2
  alpha_upper <- 1 - alpha_lower

  density_mean  <- colMeans(density_mat)
  density_lower <- apply(density_mat, 2, quantile, probs = alpha_lower)
  density_upper <- apply(density_mat, 2, quantile, probs = alpha_upper)
  density_median <- apply(density_mat, 2, median)

  elapsed <- .elapsed_time(timer_start)
  .vmsg("  DP density computed in ", .format_time(elapsed), verbose = verbose)

  # --- Return ---
  result <- structure(
    list(
      grid            = grid,
      density_mean    = density_mean,
      density_median  = density_median,
      density_lower   = density_lower,
      density_upper   = density_upper,
      density_samples = density_mat,
      dp_samples      = dp_measure,
      ci_level        = credible_interval,
      computation_time = elapsed
    ),
    class = "dpmirt_dp_density"
  )

  result
}


# ============================================================================
# Internal: Extract DP Samples from Fit Object
# ============================================================================

#' Extract DP-related columns from posterior samples
#'
#' Extracts alpha, zi\[1:N\], muTilde\[1:M\], s2Tilde\[1:M\] from
#' the raw samples matrix stored in the fit object.
#'
#' @noRd
.extract_dp_samples <- function(fit) {

  # Raw samples are stored in rescaled$samples_raw via the fit construction
  # We need to get them from the other_samp or reconstruct from compiled
  raw_samples <- NULL

  # Try to get raw samples from the stored rescaling output
  # The fit object stores other_samp which contains all non-beta, non-logprob columns
  # For DPM: alpha, zi[1:N], muTilde[1:M], s2Tilde[1:M] should be in other_samp
  if (!is.null(fit$other_samp)) {
    raw_samples <- fit$other_samp
  }

  if (is.null(raw_samples)) {
    stop("Cannot find raw posterior samples in fit object. ",
         "DP density requires monitors for alpha, zi, muTilde, s2Tilde.",
         call. = FALSE)
  }

  # Verify required columns exist
  cnames <- colnames(raw_samples)

  alpha_col <- grep("^alpha$", cnames)
  zi_cols   <- grep("^zi\\[", cnames)
  mu_cols   <- grep("^muTilde\\[", cnames)
  s2_cols   <- grep("^s2Tilde\\[", cnames)

  if (length(alpha_col) == 0) {
    stop("alpha not found in posterior samples. ",
         "Ensure alpha is in monitors.", call. = FALSE)
  }
  if (length(zi_cols) == 0) {
    stop("zi[...] not found in posterior samples. ",
         "Ensure zi is in monitors.", call. = FALSE)
  }
  if (length(mu_cols) == 0) {
    stop("muTilde[...] not found in posterior samples. ",
         "Ensure muTilde is in monitors.", call. = FALSE)
  }
  if (length(s2_cols) == 0) {
    stop("s2Tilde[...] not found in posterior samples. ",
         "Ensure s2Tilde is in monitors.", call. = FALSE)
  }

  # Extract the DP-only columns (same column order as Paganin)
  dp_col_indices <- c(alpha_col, zi_cols, mu_cols, s2_cols)
  dp_samples <- raw_samples[, dp_col_indices, drop = FALSE]

  dp_samples
}


# ============================================================================
# Internal: Get DP Measure Samples via NIMBLE
# ============================================================================

#' Reconstruct NIMBLE model and call getSamplesDPmeasure
#'
#' Follows Paganin's pattern from 3_simulateFromDPmeasure.R:
#' 1. Rebuild NIMBLE model from spec
#' 2. Configure MCMC with DP-only monitors
#' 3. Compile model and MCMC
#' 4. Populate compiled MCMC with posterior samples via matrix2mv()
#' 5. Call getSamplesDPmeasure()
#'
#' @noRd
.get_dp_measure_samples <- function(fit, dp_samples, verbose) {

  # Get the spec from the compiled object
  spec <- fit$compiled$spec

  if (is.null(spec)) {
    stop("Cannot find model specification in fit object. ",
         "DP density requires the compiled model reference.",
         call. = FALSE)
  }

  # --- Rebuild NIMBLE model ---
  .vmsg("    Rebuilding NIMBLE model for DP sampling...", verbose = verbose)

  model <- nimbleModel(
    code      = spec$code,
    constants = spec$constants,
    data      = spec$data,
    inits     = spec$inits,
    calculate = FALSE
  )

  # --- Configure MCMC with DP-only monitors ---
  conf <- configureMCMC(model)
  conf$monitors <- c("alpha", "zi", "muTilde", "s2Tilde")

  mcmc <- buildMCMC(conf)

  # --- Compile ---
  .vmsg("    Compiling for DP sampling...", verbose = verbose)
  cmodel <- compileNimble(model)
  cmcmc  <- compileNimble(mcmc, project = model)

  # --- Populate compiled MCMC with posterior samples ---
  .vmsg("    Loading posterior samples into compiled MCMC...",
        verbose = verbose)

  # Use NIMBLE's internal matrix2mv to load samples
  # This is the same approach used by Paganin
  nimble:::matrix2mv(dp_samples, cmcmc$mvSamples)

  # --- Call getSamplesDPmeasure ---
  .vmsg("    Calling getSamplesDPmeasure()...", verbose = verbose)

  dp_measure <- getSamplesDPmeasure(cmcmc)

  dp_measure
}


# ============================================================================
# Internal: Evaluate Density on Grid
# ============================================================================

#' Evaluate DP mixture density on a grid
#'
#' For each posterior sample s, computes:
#'   `f_s(x) = sum_k w_{s,k} * dnorm(x, mu_{s,k}, sqrt(s2_{s,k}))`
#'
#' With optional rescaling (location shift) for unconstrained models.
#'
#' @noRd
.evaluate_dp_density <- function(dp_measure, grid, fit,
                                 apply_rescaling = TRUE) {

  n_samples <- length(dp_measure)
  n_grid    <- length(grid)

  density_mat <- matrix(0, nrow = n_samples, ncol = n_grid)

  # Get location shift if rescaling is needed
  has_location_shift <- apply_rescaling &&
    fit$config$identification == "unconstrained" &&
    !is.null(fit$location_shift)

  for (s in seq_len(n_samples)) {
    # dp_measure[[s]] is a matrix with columns: [weights, means, variances]
    dp_s <- dp_measure[[s]]

    if (is.null(dp_s) || nrow(dp_s) == 0) next

    weights  <- dp_s[, 1]
    means    <- dp_s[, 2]
    variances <- dp_s[, 3]

    # Apply rescaling: shift grid to raw (unrescaled) scale
    if (has_location_shift) {
      # For Rasch: rescaled_theta = raw_theta - location_shift
      # So to evaluate density of rescaled_theta, we need density of
      # (grid + location_shift) under the raw DP measure
      # The Jacobian is 1 for location shift
      loc_s <- fit$location_shift[min(s, length(fit$location_shift))]
      eval_grid <- grid + loc_s
    } else {
      eval_grid <- grid
    }

    # Evaluate mixture density at each grid point
    for (g in seq_len(n_grid)) {
      density_mat[s, g] <- sum(
        weights * dnorm(eval_grid[g], mean = means, sd = sqrt(variances))
      )
    }
  }

  density_mat
}


# ============================================================================
# Convenience: Compute DP Percentiles
# ============================================================================

#' Compute posterior percentile function from DP mixture
#'
#' For each posterior sample, computes the CDF of the DP mixture at
#' given theta values. Useful for percentile-based interpretation.
#'
#' @param dp_density A \code{dpmirt_dp_density} object.
#' @param theta_values Numeric vector. Theta values at which to
#'   evaluate the CDF.
#'
#' @return A matrix (n_samples x length(theta_values)) of CDF values.
#'
#' @noRd
.dp_percentile <- function(dp_density, theta_values) {

  dp_measure <- dp_density$dp_samples
  n_samples <- length(dp_measure)
  n_theta   <- length(theta_values)

  perc_mat <- matrix(0, nrow = n_samples, ncol = n_theta)

  for (s in seq_len(n_samples)) {
    dp_s <- dp_measure[[s]]
    if (is.null(dp_s) || nrow(dp_s) == 0) next

    weights   <- dp_s[, 1]
    means     <- dp_s[, 2]
    variances <- dp_s[, 3]

    for (t in seq_len(n_theta)) {
      perc_mat[s, t] <- sum(
        weights * pnorm(theta_values[t], mean = means, sd = sqrt(variances))
      )
    }
  }

  perc_mat
}


# ============================================================================
# S3 Methods
# ============================================================================

#' @rdname dpmirt_dp_density
#' @param x A \code{dpmirt_dp_density} object.
#' @param ... Additional arguments (currently unused).
#' @export
print.dpmirt_dp_density <- function(x, ...) {
  cat("DPMirt DP Density\n")
  cat("=================\n")
  cat("Grid points:      ", length(x$grid), "\n")
  cat("Grid range:       [", min(x$grid), ", ", max(x$grid), "]\n", sep = "")
  cat("Posterior samples: ", nrow(x$density_samples), "\n")
  cat("Credible interval:", x$ci_level * 100, "%\n")
  cat("Computation time: ", .format_time(x$computation_time), "\n")

  # Summary of density
  cat("\nDensity summary:\n")
  cat("  Peak density:   ", round(max(x$density_mean), 4), "\n")
  cat("  Peak location:  ", round(x$grid[which.max(x$density_mean)], 3), "\n")

  invisible(x)
}
