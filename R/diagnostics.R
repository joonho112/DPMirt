# ============================================================================
# Module 8: Model Diagnostics
# ============================================================================
# Blueprint: Section 6 - Module 8
# ============================================================================

#' Compute MCMC Diagnostics for a DPMirt Fit
#'
#' Returns a structured list of diagnostic information including effective
#' sample sizes (ESS), WAIC, log-likelihood trace, timing, and
#' DPM-specific cluster diagnostics (number of clusters, alpha posterior).
#'
#' @param fit A \code{dpmirt_fit} object from \code{\link{dpmirt}}.
#'
#' @return A \code{dpmirt_diagnostics} S3 object containing:
#' \describe{
#'   \item{ess}{List of ESS vectors for items and persons.}
#'   \item{waic}{WAIC value (if computed).}
#'   \item{loglik_trace}{Log-likelihood trace vector.}
#'   \item{n_clusters}{Posterior cluster counts (DPM only).}
#'   \item{alpha_summary}{Alpha posterior summary (DPM only).}
#'   \item{compilation_time, sampling_time, total_time}{Timing information.}
#' }
#'
#' @details
#' Effective sample size (ESS) measures the number of effectively independent
#' draws from the posterior. Low ESS (< 100) suggests poor mixing and the
#' need for longer chains or different samplers. For DPM models, the cluster
#' count trace is a key diagnostic — stable oscillation indicates convergence,
#' while monotonic trends suggest the chain has not yet mixed.
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

  diag <- list()

  # --- ESS ---
  diag$ess <- fit$ess
  diag$ess_min_items <- if (!is.null(fit$ess$items)) {
    min(fit$ess$items, na.rm = TRUE)
  } else {
    NA
  }
  diag$ess_min_theta <- if (!is.null(fit$ess$theta)) {
    min(fit$ess$theta, na.rm = TRUE)
  } else {
    NA
  }

  # --- WAIC ---
  diag$waic <- fit$waic

  # --- Log-likelihood trace ---
  diag$loglik_trace <- fit$loglik_trace

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
  cat("  Min ESS (items): ", round(x$ess_min_items, 1), "\n")
  cat("  Min ESS (theta): ", round(x$ess_min_theta, 1), "\n")

  # WAIC
  if (!is.null(x$waic)) {
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
  cat("  Compilation: ", .format_time(x$compilation_time), "\n")
  cat("  Sampling:    ", .format_time(x$sampling_time), "\n")
  cat("  Total:       ", .format_time(x$total_time), "\n")

  invisible(x)
}


#' Compare DPMirt models using information criteria
#'
#' @param ... Two or more \code{dpmirt_fit} objects.
#' @param criterion Character. Comparison criterion. Default "waic".
#' @return A data.frame ranking models by the criterion.
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
  waic_vals <- sapply(fits, function(f) {
    if (!is.null(f$waic)) f$waic else NA
  })

  # Build model labels
  labels <- sapply(fits, function(f) {
    paste0(toupper(f$config$model), "-", f$config$prior)
  })

  result <- data.frame(
    model = labels,
    waic  = waic_vals,
    stringsAsFactors = FALSE
  )
  result <- result[order(result$waic), ]
  result$delta_waic <- result$waic - min(result$waic, na.rm = TRUE)

  result
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
