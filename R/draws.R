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
#' @param use_rescaled Reserved. Only \code{TRUE} is currently supported.
#'
#' @return If one variable is requested in matrix format, a matrix
#'   (\code{n_retained_draws} x N or \code{n_retained_draws} x I). If multiple
#'   variables are requested in matrix format, a named list of matrices. In long
#'   format, a data.frame stacking all requested variables and chain/draw
#'   metadata, including retained MCMC iteration labels after \code{thin} or
#'   \code{thin2}.
#'
#' @details
#' This function provides direct access to the posterior MCMC samples
#' stored in a \code{dpmirt_fit} object. The matrix format (default) is
#' efficient for computation, while the long format is convenient for
#' visualization with ggplot2 and includes chain metadata when available.
#' Raw draw extraction
#' (\code{use_rescaled = FALSE}) is reserved until the raw/rescaled draw
#' schema is finalized.
#'
#' @examples
#' \dontrun{
#' sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
#' fit <- dpmirt(sim$response, model = "rasch", prior = "normal",
#'               niter = 5000, nburnin = 1000, seed = 123)
#'
#' # Extract theta draws as matrix (retained draws x N)
#' theta <- dpmirt_draws(fit, vars = "theta")
#' dim(theta)
#'
#' # Extract beta draws in long format
#' beta_long <- dpmirt_draws(fit, vars = "beta", format = "long")
#' head(beta_long)
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

  if (!isTRUE(use_rescaled)) {
    stop(
      "use_rescaled = FALSE is reserved; raw draw extraction is not ",
      "currently implemented. Use use_rescaled = TRUE.",
      call. = FALSE
    )
  }

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
      row_meta <- .draw_row_metadata(fit, v, niter)
      long_df <- data.frame(
        iteration = rep(seq_len(niter), times = n_units),
        index     = rep(seq_len(n_units), each = niter),
        value     = as.vector(samp),
        variable  = v,
        chain_id  = rep(row_meta$chain_id, times = n_units),
        within_chain_draw = rep(row_meta$within_chain_draw, times = n_units),
        mcmc_iteration = rep(row_meta$mcmc_iteration, times = n_units),
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


.draw_row_metadata <- function(fit, var, niter) {
  draw_index <- NULL
  if (!is.null(fit$draw_index)) {
    draw_index <- if (identical(var, "theta")) {
      fit$draw_index$theta
    } else {
      fit$draw_index$main
    }
  }

  if (is.null(draw_index) || nrow(draw_index) != niter) {
    return(data.frame(
      chain_id = rep(1L, niter),
      within_chain_draw = seq_len(niter),
      mcmc_iteration = seq_len(niter)
    ))
  }

  data.frame(
    chain_id = draw_index$chain_id,
    within_chain_draw = draw_index$within_chain_draw,
    mcmc_iteration = draw_index$mcmc_iteration
  )
}
