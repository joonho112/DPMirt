# ============================================================================
# Module 10: Alpha Prior Elicitation
# ============================================================================
# Blueprint: Section 6 - Module 10, Section 8.1
# ============================================================================

#' Elicit a Gamma Hyperprior for the DP Concentration Parameter
#'
#' Uses the DPprior package to calibrate a Gamma(a, b) hyperprior for the
#' Dirichlet Process concentration parameter \eqn{\alpha}, based on the
#' expected number of clusters \eqn{\mu_K} and a confidence level. Falls
#' back to Paganin's default Gamma(1, 3) if DPprior is not installed.
#'
#' @param N Integer. Sample size (number of persons).
#' @param mu_K Numeric or NULL. Expected number of clusters. If NULL,
#'   defaults to \code{max(3, ceiling(log(N)))}.
#' @param confidence Character. DPprior confidence level:
#'   \code{"low"}, \code{"medium"}, or \code{"high"}.
#' @param ... Additional arguments passed to \code{DPprior::DPprior_fit}.
#'
#' @return Named numeric vector \code{c(a = ..., b = ...)} for Gamma(a, b).
#'
#' @details
#' In the CRP representation used by DPMirt, the concentration parameter
#' \eqn{\alpha} controls the expected number of clusters. Larger \eqn{\alpha}
#' favors more clusters (closer to nonparametric), while smaller \eqn{\alpha}
#' concentrates mass on fewer clusters (closer to parametric).
#'
#' The DPprior package (Lee, 2026) provides a principled elicitation framework:
#' given the sample size \eqn{N} and an expected cluster count \eqn{\mu_K},
#' it calibrates Gamma(a, b) so that \eqn{E[K | \alpha] \approx \mu_K} with
#' the specified confidence level.
#'
#' The Paganin et al. (2023) default is Gamma(1, 3), which implies
#' \eqn{E[\alpha] = 1/3} — a mildly informative prior favoring few clusters.
#'
#' @references
#' Lee, J. (2026). Design-conditional prior elicitation for Dirichlet process
#' mixtures. \emph{arXiv:2602.06301}.
#'
#' Paganin, S., Paciorek, C. J., Wehrhahn, C., Rodriguez, A.,
#' Rabe-Hesketh, S., & de Valpine, P. (2023). Computational strategies
#' and estimation performance with Bayesian semiparametric item response
#' theory models. \emph{Journal of Educational and Behavioral Statistics,
#' 48}(2), 147--188.
#'
#' @examples
#' \dontrun{
#' # Default: uses Gamma(1, 3) if DPprior not installed
#' alpha <- dpmirt_alpha_prior(N = 200)
#'
#' # Specify expected clusters and confidence
#' alpha <- dpmirt_alpha_prior(N = 500, mu_K = 5, confidence = "medium")
#'
#' # Use in model fitting
#' sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
#' fit <- dpmirt(sim$response, model = "rasch", prior = "dpm",
#'               alpha_prior = alpha, niter = 5000, nburnin = 1000)
#' }
#'
#' @family prior specification
#' @seealso \code{\link{dpmirt}}, \code{\link{dpmirt_spec}}
#'
#' @export
dpmirt_alpha_prior <- function(N,
                               mu_K = NULL,
                               confidence = "medium",
                               ...) {

  if (!requireNamespace("DPprior", quietly = TRUE)) {
    message("DPprior not installed. Using default: Gamma(1, 3) ",
            "[Paganin et al., 2023].\n",
            "Install DPprior for principled alpha elicitation: ",
            "remotes::install_github(\"joonho112/DPprior\")")
    return(c(a = 1, b = 3))
  }

  # Auto-set mu_K if not provided
  if (is.null(mu_K)) {
    mu_K <- max(3, ceiling(log(N)))
    message(sprintf("Using default mu_K = %d (based on N = %d)", mu_K, N))
  }

  # Call DPprior
  fit <- DPprior::DPprior_fit(J = N, mu_K = mu_K,
                               confidence = confidence, ...)

  result <- c(a = fit$a, b = fit$b)
  message(sprintf("Alpha prior: Gamma(%.2f, %.2f) [E[alpha]=%.2f]",
                  result["a"], result["b"],
                  result["a"] / result["b"]))

  result
}
