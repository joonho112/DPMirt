# ============================================================================
# Module 5: Person/Item Estimates
# ============================================================================
# Purpose: Compute person-specific and item-specific point estimates using
#          multiple posterior summary methods (PM, CB, GR).
#
# Contains the internalized triple-goal algorithm, adapted from
# HETOP's triple_goal.R (Shen & Louis, 1998).
#
# Blueprint: Section 6 - Module 5, Appendix D, Appendix E
# ============================================================================

#' Compute posterior estimates using PM, CB, and GR methods
#'
#' Given a \code{dpmirt_fit} object, computes person-specific and
#' item-specific point estimates using multiple posterior summary methods:
#' \itemize{
#'   \item \strong{PM}: Posterior Mean — optimal for individual MSE (Goal 1)
#'   \item \strong{CB}: Constrained Bayes (Ghosh, 1992) — optimal for EDF
#'     estimation (Goal 3)
#'   \item \strong{GR}: Triple-Goal (Shen & Louis, 1998) — optimal for
#'     ranking + EDF (Goals 2 & 3)
#' }
#'
#' @param fit A \code{dpmirt_fit} object from \code{\link{dpmirt}}.
#' @param methods Character vector of methods to compute. Default
#'   \code{c("pm", "cb", "gr")}.
#' @param alpha Significance level for credible intervals. Default 0.05.
#' @param quantile_type Integer 1-9 for \code{quantile()} type parameter.
#'   Default 7 (R default).
#' @param stop_if_ties Logical. If TRUE, stop when ties detected in
#'   posterior mean ranks. Default FALSE.
#' @param include_items Logical. If TRUE, apply CB/GR to beta as well.
#'   Default TRUE.
#'
#' @return A \code{dpmirt_estimates} S3 object.
#'
#' @references
#' Ghosh, M. (1992). Constrained Bayes estimation with applications.
#' \emph{JASA, 87}(418), 533--540.
#'
#' Shen, W., & Louis, T. A. (1998). Triple-goal estimates in two-stage
#' hierarchical models. \emph{JRSS-B, 60}(2), 455--471.
#'
#' @details
#' The three estimators target different inferential goals (Shen & Louis, 1998):
#'
#' \strong{Goal 1 — Individual estimation}: Minimize individual mean squared
#' error. The posterior mean (PM) is optimal:
#' \deqn{\hat{\theta}^{PM}_j = E[\theta_j | y]}
#'
#' \strong{Goal 2 — Ranking}: Correctly rank individuals. The GR estimator
#' uses quantiles of the marginal posterior predictive distribution evaluated
#' at the posterior mean rank of each individual.
#'
#' \strong{Goal 3 — Distribution estimation}: Recover the empirical
#' distribution function (EDF). The constrained Bayes (CB) estimator rescales
#' the PM to match the correct first two moments:
#' \deqn{\hat{\theta}^{CB}_j = \bar{\theta}^{PM} + \sqrt{1 + \frac{\bar{\lambda}}{\mathrm{Var}(\theta^{PM})}} \cdot (\hat{\theta}^{PM}_j - \bar{\theta}^{PM})}
#'
#' where \eqn{\bar{\lambda} = \frac{1}{K}\sum_k \lambda_k} is the mean
#' posterior variance.
#'
#' @examples
#' \dontrun{
#' sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
#' fit <- dpmirt(sim$response, model = "rasch", prior = "normal",
#'               niter = 5000, nburnin = 1000, seed = 123)
#'
#' # Compute all three estimators
#' est <- dpmirt_estimates(fit)
#' print(est)
#'
#' # Access person estimates
#' head(est$theta)
#'
#' # Access item estimates
#' head(est$beta)
#'
#' # PM only (faster)
#' est_pm <- dpmirt_estimates(fit, methods = "pm")
#' }
#'
#' @family estimation
#' @seealso \code{\link{dpmirt}}, \code{\link{dpmirt_draws}},
#'   \code{\link{dpmirt_loss}}, \code{\link{plot.dpmirt_estimates}}
#'
#' @export
dpmirt_estimates <- function(fit,
                             methods = c("pm", "cb", "gr"),
                             alpha = 0.05,
                             quantile_type = 7,
                             stop_if_ties = FALSE,
                             include_items = TRUE) {

  if (!inherits(fit, "dpmirt_fit")) {
    stop("Input must be a dpmirt_fit object.", call. = FALSE)
  }

  methods <- match.arg(methods, c("pm", "cb", "gr"), several.ok = TRUE)

  # --- Person estimates ---
  theta_samp <- fit$theta_samp
  if (is.null(theta_samp) || ncol(theta_samp) == 0) {
    stop("No theta posterior samples found in fit object.", call. = FALSE)
  }

  # Apply .triple_goal for person parameters
  theta_result <- .triple_goal(
    s = theta_samp,
    quantile_type = quantile_type,
    stop_if_ties  = stop_if_ties
  )

  # Add credible intervals
  lower_q <- alpha / 2
  upper_q <- 1 - alpha / 2
  theta_result$theta_lower <- apply(theta_samp, 2, quantile,
                                     probs = lower_q)
  theta_result$theta_upper <- apply(theta_samp, 2, quantile,
                                     probs = upper_q)

  # --- Item estimates ---
  beta_samp <- fit$beta_samp
  beta_result <- NULL
  if (!is.null(beta_samp) && ncol(beta_samp) > 0) {
    if (include_items && any(methods %in% c("cb", "gr"))) {
      beta_result <- .triple_goal(
        s = beta_samp,
        quantile_type = quantile_type,
        stop_if_ties  = stop_if_ties
      )
      # Rename columns for beta
      names(beta_result) <- gsub("theta_", "beta_", names(beta_result))
    } else {
      # Just PM and intervals
      beta_result <- data.frame(
        beta_pm  = colMeans(beta_samp),
        beta_psd = apply(beta_samp, 2, sd)
      )
    }
    beta_result$beta_lower <- apply(beta_samp, 2, quantile,
                                     probs = lower_q)
    beta_result$beta_upper <- apply(beta_samp, 2, quantile,
                                     probs = upper_q)
  }

  # --- Lambda estimates (2PL/3PL only) ---
  lambda_result <- NULL
  if (!is.null(fit$lambda_samp)) {
    lambda_result <- data.frame(
      lambda_pm    = colMeans(fit$lambda_samp),
      lambda_psd   = apply(fit$lambda_samp, 2, sd),
      lambda_lower = apply(fit$lambda_samp, 2, quantile, probs = lower_q),
      lambda_upper = apply(fit$lambda_samp, 2, quantile, probs = upper_q)
    )
  }

  # --- Delta estimates (3PL only) ---
  delta_result <- NULL
  if (!is.null(fit$delta_samp)) {
    delta_result <- data.frame(
      delta_pm    = colMeans(fit$delta_samp),
      delta_psd   = apply(fit$delta_samp, 2, sd),
      delta_lower = apply(fit$delta_samp, 2, quantile, probs = lower_q),
      delta_upper = apply(fit$delta_samp, 2, quantile, probs = upper_q)
    )
  }

  # --- Quality flags ---
  quality_flags <- attr(theta_result, "quality_flags")
  if (is.null(quality_flags)) quality_flags <- list()

  # --- Build estimates object ---
  result <- structure(
    list(
      theta         = theta_result,
      beta          = beta_result,
      lambda        = lambda_result,
      delta         = delta_result,
      methods       = methods,
      alpha         = alpha,
      quality_flags = quality_flags
    ),
    class = "dpmirt_estimates"
  )

  result
}


# ============================================================================
# Internal Triple-Goal Implementation
# ============================================================================
# Adapted from HETOP-master/R/triple_goal.R
# Original: Shen & Louis (1998); implementation by HETOP team
# See Blueprint Appendix E for full algorithm specification.
# ============================================================================

#' Triple-goal estimator (internal)
#'
#' Computes PM, PSD, CB, and GR estimates from a (n_iter x K) matrix of
#' posterior samples.
#'
#' @param s Numeric matrix (n_iter x K) of posterior samples for K units.
#' @param quantile_type Integer 1-9. Quantile interpolation method. Default 7.
#' @param stop_if_ties Logical. If TRUE, stop on tied posterior mean ranks.
#' @param ties_method Character. Tie-breaking method for rank(). Default
#'   "random".
#'
#' @return A data.frame with columns: theta_pm, theta_psd, theta_cb,
#'   theta_gr, rbar, rhat. Attribute "quality_flags" contains stability
#'   indicators.
#'
#' @noRd
.triple_goal <- function(s,
                         quantile_type = 7,
                         stop_if_ties  = FALSE,
                         ties_method   = "random") {

  K <- ncol(s)

  # --- Step 1: Posterior Means (PM) ---
  theta_pm  <- colMeans(s)

  # --- Step 2: Grand mean ---
  etadot <- mean(theta_pm)

  # --- Step 3: Posterior Standard Deviations ---
  theta_psd <- apply(s, 2, sd)

  # --- Step 4: Posterior variances per unit ---
  lambda_k <- theta_psd^2

  # --- Step 5: Variance of posterior means ---
  var_pm <- var(theta_pm)

  # --- Quality flags ---
  quality_flags <- list()

  # --- Step 6: Constrained Bayes (CB) — Ghosh (1992) ---
  if (var_pm < 1e-6) {
    warning("Near-zero variance of posterior means; CB estimates set to PM.",
            call. = FALSE)
    theta_cb <- theta_pm
    quality_flags$cb_fallback <- TRUE
  } else {
    cb_factor <- sqrt(1 + mean(lambda_k) / var_pm)
    if (cb_factor > 5) {
      warning("CB scaling factor > 5; estimates may be unstable.",
              call. = FALSE)
      quality_flags$cb_extreme_factor <- cb_factor
    }
    theta_cb <- etadot + cb_factor * (theta_pm - etadot)
  }

  # --- Step 7: Posterior mean ranks ---
  # For each MCMC iteration, rank all K units; then average ranks
  rbar <- apply(
    t(apply(s, 1, rank, ties.method = "average")),
    2, mean
  )

  # --- Step 8: Check for ties in rbar ---
  if (stop_if_ties && any(duplicated(round(rbar, 10)))) {
    stop("Ties detected in posterior mean ranks (rbar). ",
         "Set stop_if_ties = FALSE to proceed with random tie-breaking.",
         call. = FALSE)
  }

  # --- Step 9: Integer ranks ---
  rhat <- rank(rbar, ties.method = ties_method)

  # --- Step 10: Triple-Goal (GR) — Shen & Louis (1998) ---
  # Use ECDF of all pooled posterior samples as G-hat (ISEL estimator)
  theta_gr <- quantile(
    c(s),
    probs = (2 * rhat - 1) / (2 * K),
    names = FALSE,
    type  = quantile_type
  )

  # --- Build result ---
  result <- data.frame(
    theta_pm  = theta_pm,
    theta_psd = theta_psd,
    theta_cb  = theta_cb,
    theta_gr  = theta_gr,
    rbar      = rbar,
    rhat      = as.integer(rhat)
  )
  attr(result, "quality_flags") <- quality_flags

  result
}


# ============================================================================
# S3 Methods for dpmirt_estimates
# ============================================================================

#' @rdname dpmirt_estimates
#' @param x A \code{dpmirt_estimates} object.
#' @param ... Additional arguments (currently unused).
#' @export
print.dpmirt_estimates <- function(x, ...) {
  cat("DPMirt Posterior Estimates\n")
  cat("=========================\n")
  cat("Methods:     ", paste(x$methods, collapse = ", "), "\n")
  cat("Persons:     ", nrow(x$theta), "\n")

  if (!is.null(x$beta)) {
    cat("Items:       ", nrow(x$beta), "\n")
  }
  if (!is.null(x$lambda)) {
    cat("Lambda:       included (2PL/3PL)\n")
  }
  cat("CI level:    ", (1 - x$alpha) * 100, "%\n")

  # Quality flags
  qf <- x$quality_flags
  if (length(qf) > 0) {
    cat("\nQuality Flags:\n")
    if (isTRUE(qf$cb_fallback)) {
      cat("  ! CB fell back to PM (near-zero variance)\n")
    }
    if (!is.null(qf$cb_extreme_factor)) {
      cat("  ! CB scaling factor =", round(qf$cb_extreme_factor, 2),
          "(> 5; unstable)\n")
    }
  }

  cat("\nTheta summary (first 6 persons):\n")
  theta_cols <- intersect(c("theta_pm", "theta_cb", "theta_gr"),
                           names(x$theta))
  print(head(x$theta[, theta_cols, drop = FALSE]), digits = 3)

  # Show beta/gamma summary
  if (!is.null(x$beta)) {
    beta_cols <- intersect(c("beta_pm", "beta_cb", "beta_gr"),
                            names(x$beta))
    if (length(beta_cols) > 0) {
      cat("\nItem parameter summary (first 6 items):\n")
      print(head(x$beta[, beta_cols, drop = FALSE]), digits = 3)
    }
  }

  # Show lambda summary for 2PL/3PL
  if (!is.null(x$lambda)) {
    cat("\nLambda summary (first 6 items):\n")
    lambda_cols <- intersect(c("lambda_pm", "lambda_psd", "lambda_lower",
                                "lambda_upper"), names(x$lambda))
    print(head(x$lambda[, lambda_cols, drop = FALSE]), digits = 3)
  }

  invisible(x)
}
