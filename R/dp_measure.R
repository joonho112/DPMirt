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
# 2. Reconstruct model + MCMC, populate modelValues with posterior samples
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
#'   retained-draw-specific location shift from post-hoc rescaling. Only
#'   changes output for unconstrained Rasch/location-shift fits.
#' @param verbose Logical. Print progress messages. Default TRUE.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list of class \code{dpmirt_dp_density} containing:
#' \describe{
#'   \item{grid}{Numeric vector of evaluation points.}
#'   \item{density_mean}{Numeric vector of posterior mean densities.}
#'   \item{density_lower}{Numeric vector of lower credible band.}
#'   \item{density_upper}{Numeric vector of upper credible band.}
#'   \item{density_samples}{Matrix
#'     (\code{n_retained_draws} x \code{length(grid)}) of per-retained-draw
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
#'   \item Reconstruct a NIMBLE model and MCMC with monitors set
#'     to only DP parameters.
#'   \item Populate the MCMC's sample storage with the posterior samples using
#'     NIMBLE modelValues accessors.
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
#' For 2PL/3PL IRT and SI parameterizations, full transformed-scale density
#' reconstruction would also require scale and Jacobian adjustments. The
#' current implementation applies the Rasch/location-shift contract, so DP
#' densities for transformed-scale 2PL/3PL fits should be treated as
#' diagnostic summaries rather than definitive latent-density estimates.
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

  extra <- list(...)
  if (length(extra) > 0L) {
    stop("Unused argument(s): ", paste(names(extra), collapse = ", "),
         call. = FALSE)
  }

  # --- Validate input ---
  if (!inherits(fit, "dpmirt_fit")) {
    stop("Input must be a dpmirt_fit object.", call. = FALSE)
  }

  if (is.null(fit$config$prior) || length(fit$config$prior) != 1L ||
      !identical(fit$config$prior, "dpm")) {
    stop("DP density is only available for DPM models.", call. = FALSE)
  }

  grid <- .validate_dp_density_grid(grid)
  credible_interval <- .validate_dp_credible_interval(credible_interval)
  apply_rescaling <- .validate_dp_apply_rescaling(apply_rescaling)

  .vmsg("Computing DP density...", verbose = verbose)
  timer_start <- .start_timer()

  # --- Step 1: Extract DP-related posterior samples ---
  .vmsg("  Extracting DP posterior samples...", verbose = verbose)

  dp_samples_raw <- .extract_dp_samples(fit)

  # --- Step 2: Call getSamplesDPmeasure via NIMBLE ---
  .vmsg("  Sampling from DP measure (getSamplesDPmeasure)...",
        verbose = verbose)

  dp_measure <- .get_dp_measure_samples(fit, dp_samples_raw, verbose)

  if (length(dp_measure) != nrow(dp_samples_raw)) {
    stop(
      "getSamplesDPmeasure returned ", length(dp_measure),
      " draw(s), but DP posterior samples have ", nrow(dp_samples_raw),
      " row(s).",
      call. = FALSE
    )
  }

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

.validate_dp_density_grid <- function(grid) {
  if (!is.numeric(grid) || length(grid) < 2L ||
      any(!is.finite(grid))) {
    stop(
      "grid must be a finite numeric vector with length at least 2.",
      call. = FALSE
    )
  }
  as.numeric(grid)
}


.validate_dp_credible_interval <- function(credible_interval) {
  if (!is.numeric(credible_interval) ||
      length(credible_interval) != 1L ||
      !is.finite(credible_interval) ||
      credible_interval <= 0 ||
      credible_interval >= 1) {
    stop(
      "credible_interval must be a single numeric value between 0 and 1.",
      call. = FALSE
    )
  }
  as.numeric(credible_interval)
}


.validate_dp_apply_rescaling <- function(apply_rescaling) {
  if (!is.logical(apply_rescaling) ||
      length(apply_rescaling) != 1L ||
      is.na(apply_rescaling)) {
    stop("apply_rescaling must be TRUE or FALSE.", call. = FALSE)
  }
  apply_rescaling
}

#' Extract DP-related columns from posterior samples
#'
#' Extracts alpha, zi\[1:N\], muTilde\[1:M\], s2Tilde\[1:M\] from
#' the raw samples matrix stored in the fit object.
#'
#' @noRd
.extract_dp_samples <- function(fit) {

  candidates <- list(
    other_samp = fit$other_samp,
    samples_raw = fit$samples_raw
  )
  candidates <- candidates[!vapply(candidates, is.null, logical(1))]

  if (length(candidates) == 0L) {
    stop("Cannot find raw posterior samples in fit object. ",
         "DP density requires monitors for alpha, zi, muTilde, s2Tilde.",
         call. = FALSE)
  }

  errors <- character(0)
  for (nm in names(candidates)) {
    out <- tryCatch(
      .extract_dp_samples_from_matrix(candidates[[nm]], fit),
      error = function(e) {
        errors <<- c(errors, paste0(nm, ": ", conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(out)) {
      return(out)
    }
  }

  if (length(errors) > 0L) {
    stop(paste(errors, collapse = "; "), call. = FALSE)
  }

  stop("Cannot find raw posterior samples in fit object. ",
       "DP density requires monitors for alpha, zi, muTilde, s2Tilde.",
       call. = FALSE)
}


.extract_dp_samples_from_matrix <- function(raw_samples, fit) {
  if (is.null(raw_samples)) {
    stop("Cannot find raw posterior samples in fit object. ",
         "DP density requires monitors for alpha, zi, muTilde, s2Tilde.",
         call. = FALSE)
  }

  if (is.null(dim(raw_samples)) || length(dim(raw_samples)) != 2L) {
    stop("DP posterior samples must be a matrix-like object.", call. = FALSE)
  }

  if (nrow(raw_samples) == 0L) {
    stop("DP posterior samples must contain at least one draw.", call. = FALSE)
  }

  # Verify required columns exist
  cnames <- colnames(raw_samples)

  if (is.null(cnames)) {
    stop("DP posterior samples must have column names.", call. = FALSE)
  }

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

  N <- fit$config$N
  M <- fit$config$M

  if (!is.null(N) && length(N) == 1L && is.finite(N)) {
    if (length(zi_cols) != as.integer(N)) {
      stop(
        "Expected ", as.integer(N), " zi[...] columns, found ",
        length(zi_cols), ".",
        call. = FALSE
      )
    }
  }

  if (!is.null(M) && length(M) == 1L && is.finite(M)) {
    if (length(mu_cols) != as.integer(M)) {
      stop(
        "Expected ", as.integer(M), " muTilde[...] columns, found ",
        length(mu_cols), ".",
        call. = FALSE
      )
    }
    if (length(s2_cols) != as.integer(M)) {
      stop(
        "Expected ", as.integer(M), " s2Tilde[...] columns, found ",
        length(s2_cols), ".",
        call. = FALSE
      )
    }
  } else if (length(mu_cols) != length(s2_cols)) {
    stop(
      "muTilde[...] and s2Tilde[...] column counts must match.",
      call. = FALSE
    )
  }

  # Extract the DP-only columns (same column order as Paganin)
  dp_col_indices <- c(alpha_col, zi_cols, mu_cols, s2_cols)
  dp_samples <- raw_samples[, dp_col_indices, drop = FALSE]

  if (!is.numeric(dp_samples)) {
    stop("DP posterior samples must be numeric.", call. = FALSE)
  }

  if (any(!is.finite(dp_samples))) {
    stop("DP posterior samples must contain finite numeric values.",
         call. = FALSE)
  }

  dp_samples
}


# ============================================================================
# Internal: Get DP Measure Samples via NIMBLE
# ============================================================================

#' Load a posterior sample matrix into NIMBLE modelValues storage
#'
#' This replaces the older private matrix-to-modelValues path with public
#' NIMBLE modelValues accessors.
#'
#' @noRd
.matrix_to_model_values <- function(mat, mv) {
  if (is.null(dim(mat)) || length(dim(mat)) != 2L) {
    stop("Posterior samples must be a matrix-like object.", call. = FALSE)
  }
  if (!is.numeric(mat)) {
    stop("Posterior samples must be numeric.", call. = FALSE)
  }
  if (is.null(colnames(mat))) {
    stop("Posterior samples must have column names.", call. = FALSE)
  }
  if (nrow(mat) == 0L) {
    stop("Posterior samples must contain at least one draw.", call. = FALSE)
  }
  if (is.null(mv$varNames) || is.null(mv$sizes)) {
    stop("NIMBLE modelValues object is missing variable metadata.",
         call. = FALSE)
  }

  n_draws <- nrow(mat)
  nimble::resize(mv, n_draws)

  base_names <- sub("\\[.*\\]$", "", colnames(mat))

  for (var_name in unique(base_names)) {
    if (!var_name %in% mv$varNames) {
      stop("Posterior sample variable '", var_name,
           "' is not present in NIMBLE modelValues storage.",
           call. = FALSE)
    }

    var_cols <- base_names == var_name
    expected_size <- prod(mv$sizes[[var_name]])
    if (sum(var_cols) != expected_size) {
      stop("Posterior sample dimensions for variable '", var_name,
           "' do not match NIMBLE modelValues storage: expected ",
           expected_size, " column(s), found ", sum(var_cols), ".",
           call. = FALSE)
    }

    var_mat <- as.matrix(mat[, var_cols, drop = FALSE])
    values <- lapply(seq_len(n_draws), function(i) as.vector(var_mat[i, ]))
    mv[var_name, seq_len(n_draws)] <- values
  }

  invisible(mv)
}


#' Reconstruct NIMBLE model and call getSamplesDPmeasure
#'
#' Follows Paganin's pattern from 3_simulateFromDPmeasure.R:
#' 1. Rebuild NIMBLE model from spec
#' 2. Configure MCMC with DP-only monitors
#' 3. Populate MCMC modelValues storage with posterior samples
#' 4. Call getSamplesDPmeasure()
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
  dp_monitors <- c("alpha", "zi", "muTilde", "s2Tilde")
  conf <- configureMCMC(model, monitors = dp_monitors, print = verbose)

  mcmc <- buildMCMC(conf)

  # --- Populate MCMC with posterior samples ---
  .vmsg("    Loading posterior samples into MCMC storage...",
        verbose = verbose)

  .matrix_to_model_values(dp_samples, mcmc$mvSamples)

  # --- Call getSamplesDPmeasure ---
  .vmsg("    Calling getSamplesDPmeasure()...", verbose = verbose)

  dp_measure <- getSamplesDPmeasure(mcmc, progressBar = verbose)

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

  grid <- .validate_dp_density_grid(grid)

  if (!is.list(dp_measure) || length(dp_measure) == 0L) {
    stop("dp_measure must be a non-empty list of DP draw matrices.",
         call. = FALSE)
  }

  n_samples <- length(dp_measure)
  n_grid    <- length(grid)

  density_mat <- matrix(0, nrow = n_samples, ncol = n_grid)

  # Get location shift if rescaling is needed
  has_location_shift <- apply_rescaling &&
    fit$config$identification == "unconstrained" &&
    isTRUE(fit$config$model == "rasch")

  if (has_location_shift &&
      (is.null(fit$location_shift) ||
       !is.numeric(fit$location_shift) ||
       length(fit$location_shift) != n_samples ||
       any(!is.finite(fit$location_shift)))) {
    stop(
      "location_shift must be finite and have one value per DP draw.",
      call. = FALSE
    )
  }

  for (s in seq_len(n_samples)) {
    # dp_measure[[s]] is a matrix with columns: [weights, means, variances]
    dp_s <- .validate_dp_measure_draw(dp_measure[[s]], draw_id = s)

    weights  <- dp_s[, 1]
    means    <- dp_s[, 2]
    variances <- dp_s[, 3]

    # Apply rescaling: shift grid to raw (unrescaled) scale
    if (has_location_shift) {
      # For Rasch: rescaled_theta = raw_theta - location_shift
      # So to evaluate density of rescaled_theta, we need density of
      # (grid + location_shift) under the raw DP measure
      # The Jacobian is 1 for location shift
      loc_s <- fit$location_shift[s]
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


.validate_dp_measure_draw <- function(dp_s, draw_id = NA_integer_) {
  draw_label <- if (is.na(draw_id)) "DP measure draw" else
    paste0("DP measure draw ", draw_id)

  if (is.null(dp_s) || is.null(dim(dp_s)) || length(dim(dp_s)) != 2L) {
    stop(draw_label, " must be a matrix-like object.", call. = FALSE)
  }

  if (nrow(dp_s) == 0L) {
    stop(draw_label, " must contain at least one mixture component.",
         call. = FALSE)
  }

  if (ncol(dp_s) < 3L) {
    stop(
      draw_label,
      " must have at least three columns: weights, means, variances.",
      call. = FALSE
    )
  }

  dp_s <- as.matrix(dp_s[, 1:3, drop = FALSE])

  if (!is.numeric(dp_s) || any(!is.finite(dp_s))) {
    stop(draw_label, " must contain finite numeric values.", call. = FALSE)
  }

  weights <- dp_s[, 1]
  variances <- dp_s[, 3]

  if (any(weights < 0)) {
    stop(draw_label, " contains negative mixture weights.", call. = FALSE)
  }

  if (sum(weights) <= 0) {
    stop(draw_label, " must have positive total mixture weight.",
         call. = FALSE)
  }

  if (any(variances <= 0)) {
    stop(draw_label, " contains non-positive variances.", call. = FALSE)
  }

  dp_s
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
  if (!is.list(dp_measure) || length(dp_measure) == 0L) {
    stop("dp_density$dp_samples must be a non-empty list of DP draw matrices.",
         call. = FALSE)
  }
  theta_values <- .validate_dp_density_grid(theta_values)

  n_samples <- length(dp_measure)
  n_theta   <- length(theta_values)

  perc_mat <- matrix(0, nrow = n_samples, ncol = n_theta)

  for (s in seq_len(n_samples)) {
    dp_s <- .validate_dp_measure_draw(dp_measure[[s]], draw_id = s)

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
