# ============================================================================
# Module 8: Model Diagnostics
# ============================================================================
# Blueprint: Section 6 - Module 8
# ============================================================================

#' Compute MCMC Diagnostics for a DPMirt Fit
#'
#' Returns a structured list of diagnostic information including effective
#' sample sizes (ESS), optional chain-aware R-hat, WAIC, log-likelihood trace,
#' timing, and DPM-specific cluster diagnostics (number of clusters, alpha
#' posterior).
#'
#' @param fit A \code{dpmirt_fit} object from \code{\link{dpmirt}}.
#'
#' @return A \code{dpmirt_diagnostics} S3 object containing:
#' \describe{
#'   \item{ess}{List of ESS vectors for items and persons.}
#'   \item{ess_by_chain}{Per-chain ESS data frames when chain metadata exists.}
#'   \item{ess_min_items, ess_min_theta}{Minimum finite item/person ESS.}
#'   \item{rhat}{R-hat vectors by parameter family when computable.}
#'   \item{rhat_max}{Maximum finite R-hat value, or \code{NA}.}
#'   \item{waic}{WAIC value (if computed).}
#'   \item{waic_by_chain}{Per-chain WAIC provenance when available.}
#'   \item{waic_aggregation}{How top-level WAIC was aggregated.}
#'   \item{loglik_trace}{Log-likelihood trace vector.}
#'   \item{loglik_by_chain}{Log-likelihood trace with chain metadata.}
#'   \item{chain_info}{Stored chain/run metadata when available.}
#'   \item{n_clusters}{Posterior cluster counts (DPM only).}
#'   \item{n_clusters_summary}{Summary of posterior cluster counts (DPM only).}
#'   \item{alpha_summary}{Alpha posterior summary (DPM only).}
#'   \item{compilation_time, sampling_time, total_time}{Timing information.}
#' }
#'
#' @details
#' Effective sample size (ESS) measures the number of effectively independent
#' draws from the posterior. Low ESS (< 100) suggests poor mixing and the
#' need for longer chains or different samplers. R-hat is reported only when
#' \code{draw_index} contains at least two labeled chains with retained draws.
#' Sequential resume segments are tracked in \code{run_history} but are not
#' treated as independent chains. For DPM models, the cluster count trace is a
#' first-pass visual check: stable oscillation is reassuring, while monotonic
#' trends suggest the run has not yet mixed.
#'
#' @examples
#' \dontrun{
#' sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
#' fit <- dpmirt(sim$response, model = "rasch", prior = "dpm",
#'               niter = 5000, nburnin = 1000, seed = 123)
#' diag <- dpmirt_diagnostics(fit)
#' print(diag)
#' }
#'
#' @family diagnostics
#' @seealso \code{\link{dpmirt}}, \code{\link{dpmirt_compare}},
#'   \code{\link{plot.dpmirt_fit}}
#'
#' @export
dpmirt_diagnostics <- function(fit) {

  if (!inherits(fit, "dpmirt_fit")) {
    stop("Input must be a dpmirt_fit object.", call. = FALSE)
  }

  chain_diag <- fit$diagnostics
  diag <- list()

  # --- ESS ---
  diag$ess <- fit$ess
  diag$ess_by_chain <- if (!is.null(chain_diag$ess$by_chain)) {
    chain_diag$ess$by_chain
  } else {
    NULL
  }
  diag$ess_min_items <- .min_finite_or_na(fit$ess$items)
  diag$ess_min_theta <- .min_finite_or_na(fit$ess$theta)

  # --- R-hat ---
  diag$rhat <- if (!is.null(chain_diag$rhat)) chain_diag$rhat else NULL
  diag$rhat_max <- .max_rhat_or_na(diag$rhat)

  # --- WAIC ---
  diag$waic <- fit$waic
  diag$waic_by_chain <- if (!is.null(chain_diag$waic$by_chain)) {
    chain_diag$waic$by_chain
  } else {
    NULL
  }
  diag$waic_aggregation <- if (!is.null(chain_diag$waic$aggregation)) {
    chain_diag$waic$aggregation
  } else {
    NULL
  }

  # --- Log-likelihood trace ---
  diag$loglik_trace <- fit$loglik_trace
  diag$loglik_by_chain <- if (!is.null(chain_diag$loglik$by_chain)) {
    chain_diag$loglik$by_chain
  } else {
    NULL
  }
  diag$chain_info <- fit$chain_info

  # --- DPM-specific diagnostics ---
  if (fit$config$prior == "dpm") {
    ci <- fit$cluster_info
    if (!is.null(ci)) {
      diag$n_clusters <- ci$n_clusters
      diag$n_clusters_summary <- .summarize_n_clusters(ci$n_clusters)
      diag$alpha_summary <- ci$alpha_summary
    } else {
      diag$n_clusters <- NULL
      diag$n_clusters_summary <- NULL
      diag$alpha_summary <- NULL
    }
  }

  # --- Timing ---
  diag$compilation_time <- fit$compilation_time
  diag$sampling_time    <- fit$sampling_time
  diag$total_time       <- fit$total_time

  class(diag) <- "dpmirt_diagnostics"
  diag
}


.min_finite_or_na <- function(x) {
  if (is.null(x) || length(x) == 0L) {
    return(NA_real_)
  }
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    NA_real_
  } else {
    min(x)
  }
}


.max_rhat_or_na <- function(rhat) {
  if (is.null(rhat) || length(rhat) == 0L) {
    return(NA_real_)
  }
  vals <- unlist(rhat, use.names = FALSE)
  vals <- as.numeric(vals)
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0L) {
    NA_real_
  } else {
    max(vals)
  }
}


.format_diag_time <- function(seconds) {
  if (is.null(seconds) || length(seconds) == 0L || !is.finite(seconds[[1]])) {
    return("not available")
  }
  .format_time(seconds[[1]])
}


#' Summarize the posterior of number of clusters
#' @noRd
.summarize_n_clusters <- function(n_clusters) {
  if (is.null(n_clusters) || length(n_clusters) == 0) return(NULL)

  list(
    mean   = mean(n_clusters),
    median = median(n_clusters),
    sd     = sd(n_clusters),
    q025   = quantile(n_clusters, 0.025, names = FALSE),
    q975   = quantile(n_clusters, 0.975, names = FALSE),
    min    = min(n_clusters),
    max    = max(n_clusters),
    mode   = .numeric_mode(n_clusters)
  )
}


#' Compute mode of a numeric vector
#' @noRd
.numeric_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


#' @rdname dpmirt_diagnostics
#' @param x A \code{dpmirt_diagnostics} object.
#' @param ... Additional arguments (currently unused).
#' @export
print.dpmirt_diagnostics <- function(x, ...) {
  cat("DPMirt MCMC Diagnostics\n")
  cat("=======================\n\n")

  # ESS
  cat("Effective Sample Size (ESS):\n")
  cat("  Min ESS (items): ", .format_diag_number(x$ess_min_items), "\n")
  cat("  Min ESS (theta): ", .format_diag_number(x$ess_min_theta), "\n")

  if (!is.null(x$rhat) && is.finite(x$rhat_max)) {
    cat("\nR-hat:\n")
    cat("  Max R-hat:       ", round(x$rhat_max, 3), "\n")
  }

  # WAIC
  if (!is.null(x$waic) && length(x$waic) > 0L && is.finite(x$waic[[1]])) {
    cat("\nWAIC: ", round(x$waic, 2), "\n")
  }

  # DPM-specific
  if (!is.null(x$n_clusters_summary)) {
    s <- x$n_clusters_summary
    cat("\nNumber of Clusters (posterior):\n")
    cat("  Mean:   ", round(s$mean, 1), "\n")
    cat("  Median: ", s$median, "\n")
    cat("  Mode:   ", s$mode, "\n")
    cat("  95% CI: [", s$q025, ", ", s$q975, "]\n")
    cat("  Range:  [", s$min, ", ", s$max, "]\n")
  }

  if (!is.null(x$alpha_summary)) {
    a <- x$alpha_summary
    cat("\nDP Concentration (alpha, posterior):\n")
    cat("  Mean:   ", round(a$mean, 3), "\n")
    cat("  Median: ", round(a$median, 3), "\n")
    cat("  SD:     ", round(a$sd, 3), "\n")
    cat("  95% CI: [", round(a$q025, 3), ", ", round(a$q975, 3), "]\n")
  }

  # Timing
  cat("\nTiming:\n")
  cat("  Compilation: ", .format_diag_time(x$compilation_time), "\n")
  cat("  Sampling:    ", .format_diag_time(x$sampling_time), "\n")
  cat("  Total:       ", .format_diag_time(x$total_time), "\n")

  invisible(x)
}


.format_diag_number <- function(x) {
  if (is.null(x) || length(x) == 0L || !is.finite(x[[1]])) {
    "not available"
  } else {
    as.character(round(x[[1]], 1))
  }
}


#' Compare DPMirt models using information criteria
#'
#' @param ... Two or more \code{dpmirt_fit} objects.
#' @param criterion Character. Comparison criterion. Currently only
#'   \code{"waic"} is supported.
#' @return A data.frame ranking models by the criterion with columns
#'   \code{model}, \code{waic}, \code{delta_waic}, and
#'   \code{waic_aggregation}. Fits with unavailable WAIC are kept in the
#'   table with \code{NA} deltas; all-missing WAIC values raise an error.
#'
#' @details
#' For chain-labeled multi-run fits, \code{waic_aggregation} records whether
#' the top-level WAIC is a mean of per-run WAIC values. Such values are
#' retained for provenance and backward compatibility but should not be
#' interpreted as pooled posterior WAIC.
#'
#' @examples
#' \dontrun{
#' sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
#' fit1 <- dpmirt(sim$response, model = "rasch", prior = "normal",
#'                niter = 5000, nburnin = 1000, seed = 123)
#' fit2 <- dpmirt(sim$response, model = "rasch", prior = "dpm",
#'                niter = 5000, nburnin = 1000, seed = 123)
#'
#' # Compare via WAIC
#' comp <- dpmirt_compare(fit1, fit2)
#' print(comp)
#' }
#'
#' @family diagnostics
#' @seealso \code{\link{dpmirt}}, \code{\link{dpmirt_diagnostics}}
#'
#' @export
dpmirt_compare <- function(..., criterion = "waic") {
  fits <- list(...)

  if (length(fits) < 2L) {
    stop("At least two dpmirt_fit objects are required.", call. = FALSE)
  }

  # Validate
  for (i in seq_along(fits)) {
    if (!inherits(fits[[i]], "dpmirt_fit")) {
      stop("Argument ", i, " is not a dpmirt_fit object.", call. = FALSE)
    }
  }

  if (criterion != "waic") {
    stop("Only 'waic' criterion is currently supported.", call. = FALSE)
  }

  # Extract WAIC values
  waic_vals <- vapply(fits, function(f) {
    .normalize_fit_waic(f$waic)
  }, numeric(1), USE.NAMES = FALSE)

  if (all(is.na(waic_vals))) {
    stop(
      "No comparable WAIC values are available. Refit with compute_waic = TRUE ",
      "or compare fits with another validated criterion.",
      call. = FALSE
    )
  }

  if (any(is.na(waic_vals))) {
    warning(
      "Some fits have unavailable WAIC and will not receive a finite ",
      "delta_waic.",
      call. = FALSE
    )
  }

  waic_aggregation <- vapply(fits, function(f) {
    agg <- NULL
    if (!is.null(f$diagnostics) && !is.null(f$diagnostics$waic)) {
      agg <- f$diagnostics$waic$aggregation
    }
    if (is.null(agg) && !is.null(f$chain_info) && nrow(f$chain_info) > 1L) {
      agg <- "mean_of_chain_waic"
    }
    if (is.null(agg)) NA_character_ else as.character(agg[[1]])
  }, character(1), USE.NAMES = FALSE)

  if (any(waic_aggregation == "mean_of_chain_waic", na.rm = TRUE)) {
    warning(
      "At least one fit uses mean-of-run WAIC provenance. This is not a ",
      "pooled posterior WAIC; interpret model ranking cautiously.",
      call. = FALSE
    )
  }

  delta_waic <- rep(NA_real_, length(waic_vals))
  finite_waic <- is.finite(waic_vals)
  delta_waic[finite_waic] <- waic_vals[finite_waic] -
    min(waic_vals[finite_waic])

  # Build model labels
  labels <- vapply(fits, function(f) {
    paste0(toupper(f$config$model), "-", f$config$prior)
  }, character(1), USE.NAMES = FALSE)

  result <- data.frame(
    model = labels,
    waic  = waic_vals,
    delta_waic = delta_waic,
    waic_aggregation = waic_aggregation,
    stringsAsFactors = FALSE
  )
  result <- result[order(is.na(result$waic), result$waic), ]
  rownames(result) <- NULL

  result
}


.normalize_fit_waic <- function(waic) {
  if (is.null(waic) || length(waic) == 0L) {
    return(NA_real_)
  }
  waic <- suppressWarnings(as.numeric(waic[[1]]))
  if (!is.finite(waic)) NA_real_ else waic
}


# ============================================================================
# Internal: Extract Cluster Info from Posterior Samples
# ============================================================================

#' Extract DPM cluster diagnostics from posterior samples
#'
#' Computes the number of active clusters per MCMC iteration and
#' summarizes the posterior distribution of the DP concentration parameter.
#'
#' @param samples Matrix. Posterior samples with zi\[1:N\] and alpha columns.
#' @param N Integer. Number of persons.
#' @return A list with n_clusters (numeric vector) and alpha_summary.
#' @noRd
.extract_cluster_info <- function(samples, N) {

  cnames <- colnames(samples)

  # --- Number of active clusters per iteration ---
  zi_cols <- grep("^zi\\[", cnames)
  n_clusters <- NULL

  if (length(zi_cols) > 0) {
    zi_samp <- samples[, zi_cols, drop = FALSE]

    # Count unique cluster IDs per iteration
    n_clusters <- apply(zi_samp, 1, function(row) {
      length(unique(row))
    })
  }

  # --- Alpha posterior summary ---
  alpha_col <- grep("^alpha$", cnames)
  alpha_summary <- NULL

  if (length(alpha_col) > 0) {
    alpha_samp <- samples[, alpha_col]
    alpha_summary <- list(
      mean   = mean(alpha_samp),
      median = median(alpha_samp),
      sd     = sd(alpha_samp),
      q025   = quantile(alpha_samp, 0.025, names = FALSE),
      q975   = quantile(alpha_samp, 0.975, names = FALSE)
    )
  }

  list(
    n_clusters    = n_clusters,
    alpha_summary = alpha_summary
  )
}
