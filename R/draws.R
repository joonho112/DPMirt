# ============================================================================
# Module 12: Draws Extraction
# ============================================================================
# Blueprint: Section 6 - Module 12
# ============================================================================

#' Extract posterior draws from a DPMirt fit
#'
#' Provides convenient access to posterior samples in matrix or long format.
#'
#' @param fit A \code{dpmirt_fit} object.
#' @param vars Character vector of variables to extract. Options:
#'   \code{"theta"}, \code{"beta"}, \code{"lambda"}, \code{"delta"}.
#' @param format Character. Output format: \code{"matrix"} or \code{"long"}.
#' @param use_rescaled Logical. If TRUE (default), return rescaled samples.
#'
#' @return A matrix (niter x N or niter x I) or long data.frame.
#'
#' @details
#' This function provides direct access to the posterior MCMC samples
#' stored in a \code{dpmirt_fit} object. The matrix format (default) is
#' efficient for computation, while the long format is convenient for
#' visualization with ggplot2.
#'
#' @examples
#' \dontrun{
#' sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
#' fit <- dpmirt(sim$response, model = "rasch", prior = "normal",
#'               niter = 5000, nburnin = 1000, seed = 123)
#'
#' # Extract theta draws as matrix (niter x N)
#' theta <- dpmirt_draws(fit, vars = "theta")
#' dim(theta$theta)
#'
#' # Extract beta draws in long format
#' beta_long <- dpmirt_draws(fit, vars = "beta", format = "long")
#' head(beta_long$beta)
#' }
#'
#' @family estimation
#' @seealso \code{\link{dpmirt}}, \code{\link{dpmirt_estimates}}
#'
#' @export
dpmirt_draws <- function(fit,
                         vars = c("theta", "beta", "lambda", "delta"),
                         format = c("matrix", "long"),
                         use_rescaled = TRUE) {

  if (!inherits(fit, "dpmirt_fit")) {
    stop("Input must be a dpmirt_fit object.", call. = FALSE)
  }

  vars <- match.arg(vars, c("theta", "beta", "lambda", "delta"),
                    several.ok = TRUE)
  format <- match.arg(format)

  # Collect requested draws
  draws_list <- list()

  for (v in vars) {
    samp_name <- paste0(v, "_samp")
    samp <- fit[[samp_name]]

    if (is.null(samp)) {
      next
    }

    if (format == "matrix") {
      draws_list[[v]] <- samp
    } else {
      # Long format
      niter <- nrow(samp)
      n_units <- ncol(samp)
      long_df <- data.frame(
        iteration = rep(seq_len(niter), times = n_units),
        index     = rep(seq_len(n_units), each = niter),
        value     = as.vector(samp),
        variable  = v,
        stringsAsFactors = FALSE
      )
      draws_list[[v]] <- long_df
    }
  }

  if (length(draws_list) == 0) {
    stop("No samples found for requested variables: ",
         paste(vars, collapse = ", "), call. = FALSE)
  }

  if (format == "matrix") {
    if (length(draws_list) == 1) {
      return(draws_list[[1]])
    }
    return(draws_list)
  } else {
    return(do.call(rbind, draws_list))
  }
}
