# ============================================================================
# S3 Methods for DPMirt Objects
# ============================================================================
# Blueprint: Section 5.4
# ============================================================================

# ============================================================================
# dpmirt_fit methods
# ============================================================================

#' Methods for dpmirt_fit Objects
#'
#' Print, summarize, and extract coefficients from a fitted DPMirt model.
#'
#' \code{print} displays a concise one-screen summary: model type, data
#' dimensions, MCMC settings, WAIC, minimum ESS, and DPM cluster summary.
#'
#' \code{summary} displays a comprehensive report including item parameter
#' posterior means and SDs, person ability range, DPM diagnostics, and timing.
#'
#' \code{coef} extracts posterior mean point estimates for items
#' (\code{type = "items"}) or persons (\code{type = "persons"}).
#'
#' @param x,object A \code{dpmirt_fit} object from \code{\link{dpmirt}}.
#' @param type For \code{coef}: character, \code{"items"} (default)
#'   or \code{"persons"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return
#' \code{print} and \code{summary} return the input invisibly.
#' \code{coef} returns a data.frame of posterior mean estimates with
#' one row per item or person.
#'
#' @examples
#' \dontrun{
#' sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
#' fit <- dpmirt(sim$response, model = "rasch", prior = "normal",
#'               niter = 5000, nburnin = 1000, seed = 123)
#'
#' # Print concise summary
#' print(fit)
#'
#' # Detailed summary with item parameters
#' summary(fit)
#'
#' # Extract item difficulty estimates
#' coef(fit, type = "items")
#'
#' # Extract person ability estimates
#' coef(fit, type = "persons")
#' }
#'
#' @family model fitting
#' @seealso \code{\link{dpmirt}} for model fitting,
#'   \code{\link{dpmirt_estimates}} for PM/CB/GR estimators,
#'   \code{\link{dpmirt_diagnostics}} for detailed MCMC diagnostics
#'
#' @name dpmirt_fit-methods
NULL

#' @rdname dpmirt_fit-methods
#' @export
print.dpmirt_fit <- function(x, ...) {
  cat("DPMirt Model Fit\n")
  cat("================\n")
  cat("Model:           ", toupper(x$config$model), "\n")
  cat("Prior:           ", x$config$prior, "\n")
  cat("Identification:  ", x$config$identification, "\n")
  cat("Persons (N):     ", x$config$N, "\n")
  cat("Items (I):       ", x$config$I, "\n")
  cat("MCMC:            ", x$config$niter, " iterations (",
      x$config$nburnin, " burnin, thin=", x$config$thin, ")\n", sep = "")

  if (!is.null(x$waic)) {
    cat("WAIC:            ", round(x$waic, 2), "\n")
  }

  cat("Total time:      ", .format_time(x$total_time), "\n")

  # ESS summary
  if (!is.null(x$ess$items) && length(x$ess$items) > 0) {
    min_ess_items <- min(x$ess$items, na.rm = TRUE)
    cat("Min ESS (items): ", round(min_ess_items), "\n")
  }
  if (!is.null(x$ess$theta) && length(x$ess$theta) > 0) {
    min_ess_theta <- min(x$ess$theta, na.rm = TRUE)
    cat("Min ESS (theta): ", round(min_ess_theta), "\n")
  }

  # DPM summary
  if (x$config$prior == "dpm" && !is.null(x$cluster_info)) {
    if (!is.null(x$cluster_info$n_clusters)) {
      cat("Clusters (mean): ", round(mean(x$cluster_info$n_clusters), 1), "\n")
    }
    if (!is.null(x$cluster_info$alpha_summary)) {
      cat("Alpha (mean):    ", round(x$cluster_info$alpha_summary$mean, 3),
          "\n")
    }
  }

  invisible(x)
}


#' @rdname dpmirt_fit-methods
#' @export
summary.dpmirt_fit <- function(object, ...) {
  cat("DPMirt Model Summary\n")
  cat("====================\n\n")

  # --- Model Configuration ---
  cat("Model Configuration:\n")
  cat("  Model:           ", toupper(object$config$model), "\n")
  cat("  Prior:           ", object$config$prior, "\n")
  if (object$config$model != "rasch") {
    cat("  Parameterization:", object$config$parameterization, "\n")
  }
  cat("  Identification:  ", object$config$identification, "\n")
  cat("  Rescaled:        ", object$config$rescale, "\n")

  # --- Data Dimensions ---
  cat("\nData:\n")
  cat("  Persons (N):", object$config$N, "\n")
  cat("  Items (I):  ", object$config$I, "\n")

  # --- MCMC Settings ---
  cat("\nMCMC Settings:\n")
  cat("  Iterations: ", object$config$niter, "\n")
  cat("  Burn-in:    ", object$config$nburnin, "\n")
  cat("  Thinning:   ", object$config$thin, "\n")
  cat("  Chains:     ", object$config$nchains, "\n")

  # --- Timing ---
  cat("\nTiming:\n")
  cat("  Compilation: ", .format_time(object$compilation_time), "\n")
  cat("  Sampling:    ", .format_time(object$sampling_time), "\n")
  cat("  Total:       ", .format_time(object$total_time), "\n")

  # --- Item Parameters ---
  is_si <- !is.null(object$config$parameterization) &&
    object$config$parameterization == "si" &&
    object$config$model != "rasch"
  is_2pl <- object$config$model %in% c("2pl", "3pl")

  if (!is.null(object$beta_samp)) {
    param_label <- if (is_si) "Item Intercept (gamma)" else "Item Difficulty (beta)"
    param_sym   <- if (is_si) "gamma" else "beta"
    cat("\n", param_label, " Summary:\n", sep = "")
    beta_pm  <- colMeans(object$beta_samp)
    beta_psd <- apply(object$beta_samp, 2, sd)
    item_summary <- data.frame(
      Mean = round(beta_pm, 3),
      SD   = round(beta_psd, 3)
    )
    rownames(item_summary) <- paste0(param_sym, "[", seq_along(beta_pm), "]")

    # Add lambda for 2PL/3PL
    if (is_2pl && !is.null(object$lambda_samp)) {
      lambda_pm  <- colMeans(object$lambda_samp)
      lambda_psd <- apply(object$lambda_samp, 2, sd)
      item_summary$Lambda_Mean <- round(lambda_pm, 3)
      item_summary$Lambda_SD   <- round(lambda_psd, 3)
    }

    # Add delta for 3PL
    if (!is.null(object$delta_samp)) {
      delta_pm  <- colMeans(object$delta_samp)
      delta_psd <- apply(object$delta_samp, 2, sd)
      item_summary$Delta_Mean <- round(delta_pm, 3)
      item_summary$Delta_SD   <- round(delta_psd, 3)
    }

    print(item_summary)
  }

  # --- Theta Summary ---
  if (!is.null(object$theta_samp)) {
    cat("\nPerson Ability (theta) Summary:\n")
    theta_pm <- colMeans(object$theta_samp)
    cat("  Range: [", round(min(theta_pm), 3), ", ",
        round(max(theta_pm), 3), "]\n")
    cat("  Mean:  ", round(mean(theta_pm), 3), "\n")
    cat("  SD:    ", round(sd(theta_pm), 3), "\n")
  }

  # --- Model Comparison ---
  if (!is.null(object$waic)) {
    cat("\nModel Comparison:\n")
    cat("  WAIC: ", round(object$waic, 2), "\n")
  }

  # --- DPM-specific info ---
  if (object$config$prior == "dpm") {
    cat("\nDPM Diagnostics:\n")

    # Cluster info
    if (!is.null(object$cluster_info)) {
      ci <- object$cluster_info

      if (!is.null(ci$alpha_summary)) {
        cat("  Alpha (concentration):\n")
        cat("    Posterior mean: ", round(ci$alpha_summary$mean, 3), "\n")
        cat("    95% CI: [", round(ci$alpha_summary$q025, 3), ", ",
            round(ci$alpha_summary$q975, 3), "]\n")
      }

      if (!is.null(ci$n_clusters)) {
        n_cl <- ci$n_clusters
        cat("  Number of clusters:\n")
        cat("    Posterior mean: ", round(mean(n_cl), 1), "\n")
        cat("    Mode:          ", .numeric_mode(n_cl), "\n")
        cat("    Range: [", min(n_cl), ", ", max(n_cl), "]\n")
      }
    } else if (!is.null(object$other_samp)) {
      # Fallback: extract from other_samp
      alpha_col <- grep("^alpha$", colnames(object$other_samp))
      if (length(alpha_col) > 0) {
        alpha_samp <- object$other_samp[, alpha_col]
        cat("  Alpha (concentration):\n")
        cat("    Posterior mean: ", round(mean(alpha_samp), 3), "\n")
        cat("    95% CI: [", round(quantile(alpha_samp, 0.025), 3), ", ",
            round(quantile(alpha_samp, 0.975), 3), "]\n")
      }
    }

    # DP density
    if (!is.null(object$dp_density)) {
      cat("  DP density:    computed (", length(object$dp_density$grid),
          " grid points)\n", sep = "")
    } else {
      cat("  DP density:    not computed ",
          "(use dpmirt_dp_density(fit) to compute)\n")
    }
  }

  invisible(object)
}


#' @rdname dpmirt_fit-methods
#' @export
coef.dpmirt_fit <- function(object, type = c("items", "persons"), ...) {
  type <- match.arg(type)

  if (type == "items") {
    is_si <- !is.null(object$config$parameterization) &&
      object$config$parameterization == "si" &&
      object$config$model != "rasch"

    if (is_si) {
      result <- data.frame(gamma = colMeans(object$beta_samp))
    } else {
      result <- data.frame(beta = colMeans(object$beta_samp))
    }
    if (!is.null(object$lambda_samp)) {
      result$lambda <- colMeans(object$lambda_samp)
    }
    if (!is.null(object$delta_samp)) {
      result$delta <- colMeans(object$delta_samp)
    }
    rownames(result) <- paste0("item_", seq_len(nrow(result)))
    return(result)
  }

  if (type == "persons") {
    result <- data.frame(
      theta = colMeans(object$theta_samp)
    )
    rownames(result) <- paste0("person_", seq_len(nrow(result)))
    return(result)
  }
}
