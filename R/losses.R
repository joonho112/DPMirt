# ============================================================================
# Module 6: Loss Functions
# ============================================================================
# Purpose: Evaluate estimator performance when true values are known.
# Blueprint: Section 6 - Module 6
# Implemented early for simulation study readiness.
# ============================================================================

#' Evaluate estimator performance using loss functions
#'
#' Computes loss metrics comparing estimated parameters to true values.
#' Designed for simulation studies where true parameter values are known.
#'
#' @param estimates A \code{dpmirt_estimates} object.
#' @param true_theta Numeric vector of true person abilities.
#' @param true_beta Numeric vector of true item difficulties (optional).
#'   For SI parameterization, this should be the true gamma (intercept) values.
#' @param true_lambda Numeric vector of true item discriminations (optional).
#'   Only relevant for 2PL/3PL models.
#' @param metrics Character vector of metrics to compute. Options:
#'   \code{"msel"} (mean squared error loss), \code{"mselr"} (MSE of ranks),
#'   \code{"ks"} (Kolmogorov-Smirnov statistic).
#' @param custom_loss Optional custom loss function with signature
#'   \code{function(estimate, true)}.
#'
#' @return A data.frame with loss values for each method x parameter
#'   combination.
#'
#' @details
#' Three built-in loss metrics measure different aspects of estimation quality:
#'
#' \strong{MSEL} (Mean Squared Error Loss): Measures individual-level accuracy.
#' \deqn{MSEL = \frac{1}{K} \sum_{k=1}^{K} (\hat{\theta}_k - \theta_k)^2}
#'
#' \strong{MSELR} (MSE of Ranks): Measures ranking accuracy using normalized ranks.
#' \deqn{MSELR = \frac{1}{K} \sum_{k=1}^{K} \left(\frac{R(\hat{\theta}_k)}{K} - \frac{R(\theta_k)}{K}\right)^2}
#'
#' \strong{KS} (Kolmogorov-Smirnov): Measures distributional accuracy.
#' \deqn{KS = \sup_t |F_{\hat{\theta}}(t) - F_{\theta}(t)|}
#'
#' Custom loss functions can be supplied via \code{custom_loss}; they must
#' accept two vectors (estimates, true values) and return a scalar.
#'
#' @references
#' Shen, W., & Louis, T. A. (1998). Triple-goal estimates in two-stage
#' hierarchical models. \emph{JRSS-B, 60}(2), 455--471.
#'
#' @examples
#' \dontrun{
#' sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
#' fit <- dpmirt(sim$response, model = "rasch", prior = "normal",
#'               niter = 5000, nburnin = 1000, seed = 123)
#' est <- dpmirt_estimates(fit)
#'
#' # Evaluate against true values
#' loss <- dpmirt_loss(est, true_theta = sim$theta, true_beta = sim$beta)
#' print(loss)
#'
#' # With custom loss function
#' mae <- function(est, true) mean(abs(est - true))
#' loss2 <- dpmirt_loss(est, true_theta = sim$theta, custom_loss = mae)
#' }
#'
#' @family simulation
#' @seealso \code{\link{dpmirt_estimates}}, \code{\link{dpmirt_simulate}}
#'
#' @export
dpmirt_loss <- function(estimates,
                        true_theta,
                        true_beta = NULL,
                        true_lambda = NULL,
                        metrics = c("msel", "mselr", "ks"),
                        custom_loss = NULL) {

  if (!inherits(estimates, "dpmirt_estimates")) {
    stop("Input must be a dpmirt_estimates object.", call. = FALSE)
  }

  metrics <- match.arg(metrics, c("msel", "mselr", "ks"), several.ok = TRUE)

  results <- list()

  # --- Person parameter losses ---
  theta_df <- estimates$theta

  methods_available <- c()
  if ("theta_pm" %in% names(theta_df)) methods_available <- c(methods_available, "pm")
  if ("theta_cb" %in% names(theta_df)) methods_available <- c(methods_available, "cb")
  if ("theta_gr" %in% names(theta_df)) methods_available <- c(methods_available, "gr")

  for (method in methods_available) {
    est_col <- paste0("theta_", method)
    est_vals <- theta_df[[est_col]]

    row <- data.frame(
      parameter = "theta",
      method    = method,
      stringsAsFactors = FALSE
    )

    if ("msel" %in% metrics) {
      row$msel <- mean((est_vals - true_theta)^2)
    }
    if ("mselr" %in% metrics) {
      N <- length(true_theta)
      rank_est  <- rank(est_vals) / N
      rank_true <- rank(true_theta) / N
      row$mselr <- mean((rank_est - rank_true)^2)
    }
    if ("ks" %in% metrics) {
      row$ks <- .ks_stat(est_vals, true_theta)
    }
    if (!is.null(custom_loss)) {
      row$custom <- custom_loss(est_vals, true_theta)
    }

    results[[length(results) + 1]] <- row
  }

  # --- Item parameter losses ---
  if (!is.null(true_beta) && !is.null(estimates$beta)) {
    beta_df <- estimates$beta
    beta_methods <- c()
    if ("beta_pm" %in% names(beta_df)) beta_methods <- c(beta_methods, "pm")
    if ("beta_cb" %in% names(beta_df)) beta_methods <- c(beta_methods, "cb")
    if ("beta_gr" %in% names(beta_df)) beta_methods <- c(beta_methods, "gr")

    for (method in beta_methods) {
      est_col <- paste0("beta_", method)
      est_vals <- beta_df[[est_col]]

      row <- data.frame(
        parameter = "beta",
        method    = method,
        stringsAsFactors = FALSE
      )

      if ("msel" %in% metrics) {
        row$msel <- mean((est_vals - true_beta)^2)
      }
      if ("mselr" %in% metrics) {
        I_val <- length(true_beta)
        rank_est  <- rank(est_vals) / I_val
        rank_true <- rank(true_beta) / I_val
        row$mselr <- mean((rank_est - rank_true)^2)
      }
      if ("ks" %in% metrics) {
        row$ks <- .ks_stat(est_vals, true_beta)
      }

      results[[length(results) + 1]] <- row
    }
  }

  # --- Lambda parameter losses (2PL/3PL) ---
  if (!is.null(true_lambda) && !is.null(estimates$lambda)) {
    lambda_df <- estimates$lambda

    # Lambda only has PM estimates (no CB/GR — not hierarchically exchangeable)
    if ("lambda_pm" %in% names(lambda_df)) {
      est_vals <- lambda_df$lambda_pm

      row <- data.frame(
        parameter = "lambda",
        method    = "pm",
        stringsAsFactors = FALSE
      )

      if ("msel" %in% metrics) {
        row$msel <- mean((est_vals - true_lambda)^2)
      }
      if ("mselr" %in% metrics) {
        I_val <- length(true_lambda)
        rank_est  <- rank(est_vals) / I_val
        rank_true <- rank(true_lambda) / I_val
        row$mselr <- mean((rank_est - rank_true)^2)
      }
      if ("ks" %in% metrics) {
        row$ks <- .ks_stat(est_vals, true_lambda)
      }
      if (!is.null(custom_loss)) {
        row$custom <- custom_loss(est_vals, true_lambda)
      }

      results[[length(results) + 1]] <- row
    }
  }

  do.call(rbind, results)
}


#' Kolmogorov-Smirnov statistic between two samples
#' @noRd
.ks_stat <- function(x, y) {
  # KS = max|F_x(t) - F_y(t)|
  n_x <- length(x)
  n_y <- length(y)
  combined <- sort(c(x, y))
  ecdf_x <- ecdf(x)
  ecdf_y <- ecdf(y)
  max(abs(ecdf_x(combined) - ecdf_y(combined)))
}
