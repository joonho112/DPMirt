# ============================================================================
# Module 9: Simulation
# ============================================================================
# Blueprint: Section 6 - Module 9
# Phase 5: Full IRTsimrel integration with fallback simulation.
#
# Workflow:
#   1. If use_irtsimrel=TRUE and IRTsimrel is installed:
#      eqc_calibrate() -> simulate_response_data() -> dpmirt_sim
#   2. Otherwise: Paganin-style fallback (no reliability targeting)
# ============================================================================

#' Simulate IRT response data
#'
#' Generates binary response data under known IRT parameters.
#' When IRTsimrel is available, uses EQC (Empirical Quadrature Calibration)
#' to achieve a target marginal reliability. Otherwise falls back to
#' Paganin-style simulation without reliability targeting.
#'
#' @param n_persons Integer. Number of persons.
#' @param n_items Integer. Number of items.
#' @param model Character. IRT model type: \code{"rasch"} or \code{"2pl"}.
#' @param target_rho Numeric in (0, 1). Target marginal reliability. Only used
#'   when IRTsimrel is available. Default 0.8.
#' @param latent_shape Character. Shape of latent ability distribution.
#'   When using IRTsimrel, supports all 12 shapes: \code{"normal"},
#'   \code{"bimodal"}, \code{"trimodal"}, \code{"multimodal"},
#'   \code{"skew_pos"}, \code{"skew_neg"}, \code{"heavy_tail"},
#'   \code{"light_tail"}, \code{"uniform"}, \code{"floor"},
#'   \code{"ceiling"}, \code{"custom"}.
#'   Fallback supports: \code{"normal"}, \code{"bimodal"}, \code{"skewed"}.
#' @param item_source Character. Source for item parameters in IRTsimrel.
#'   One of \code{"irw"} (default), \code{"parametric"},
#'   \code{"hierarchical"}, \code{"custom"}.
#' @param reliability_metric Character. Reliability metric for EQC calibration.
#'   \code{"msem"} (default) or \code{"info"}.
#' @param latent_params List. Additional parameters passed to
#'   \code{IRTsimrel::sim_latentG()} (e.g.,
#'   \code{list(shape_params = list(delta = 0.8))}).
#' @param item_params List. Additional parameters passed to
#'   \code{IRTsimrel::sim_item_params()} (e.g.,
#'   \code{list(discrimination_params = list(rho = -0.3))}).
#' @param M Integer. Quadrature sample size for EQC. Default 10000.
#' @param seed Integer or NULL. Random seed for reproducibility.
#' @param use_irtsimrel Logical. If TRUE (default), attempt to use
#'   IRTsimrel for reliability-targeted simulation. Falls back to internal
#'   simulation if IRTsimrel is not installed.
#' @param verbose Logical. Print progress messages. Default FALSE.
#' @param ... Additional arguments (currently unused).
#'
#' @return A \code{dpmirt_sim} S3 object containing:
#' \describe{
#'   \item{response}{N x I binary response matrix}
#'   \item{theta}{True person abilities (length N)}
#'   \item{beta}{True item difficulties (length I)}
#'   \item{lambda}{True discriminations (length I for 2PL, NULL for Rasch)}
#'   \item{n_persons, n_items, model}{Simulation settings}
#'   \item{reliability}{Achieved marginal reliability}
#'   \item{target_rho}{Requested target reliability (NULL if fallback)}
#'   \item{latent_shape}{Distribution shape used}
#'   \item{eqc_result}{EQC calibration result (NULL if fallback)}
#'   \item{method}{Character: "irtsimrel" or "fallback"}
#' }
#'
#' @examples
#' \dontrun{
#' # Simple fallback simulation
#' sim <- dpmirt_simulate(200, 25, model = "rasch", seed = 42)
#'
#' # With IRTsimrel (reliability-targeted)
#' sim <- dpmirt_simulate(200, 25, model = "rasch",
#'                        target_rho = 0.85,
#'                        latent_shape = "bimodal",
#'                        seed = 42)
#'
#' # Full simulation study workflow
#' sim <- dpmirt_simulate(200, 25, model = "rasch",
#'                        target_rho = 0.8,
#'                        latent_shape = "bimodal")
#' fit <- dpmirt(sim$response, model = "rasch", prior = "dpm")
#' est <- dpmirt_estimates(fit)
#' dpmirt_loss(est, true_theta = sim$theta, true_beta = sim$beta)
#' }
#'
#' @family simulation
#' @seealso \code{\link{dpmirt}}, \code{\link{dpmirt_loss}},
#'   \code{\link{plot.dpmirt_sim}}
#' @export
dpmirt_simulate <- function(n_persons,
                            n_items,
                            model = c("rasch", "2pl", "3pl"),
                            target_rho = 0.8,
                            latent_shape = "normal",
                            item_source = "irw",
                            reliability_metric = c("msem", "info"),
                            latent_params = list(),
                            item_params = list(),
                            M = 10000L,
                            seed = NULL,
                            use_irtsimrel = TRUE,
                            verbose = FALSE,
                            ...) {

  model <- match.arg(model)
  reliability_metric <- match.arg(reliability_metric)
  .set_seed(seed)

  # Validate inputs
  stopifnot(is.numeric(n_persons), length(n_persons) == 1, n_persons >= 2)
  stopifnot(is.numeric(n_items), length(n_items) == 1, n_items >= 2)
  stopifnot(is.numeric(target_rho), length(target_rho) == 1,
            target_rho > 0, target_rho < 1)

  # --- Try IRTsimrel integration (Rasch/2PL only; 3PL uses fallback) ---
  if (model != "3pl" &&
      use_irtsimrel && requireNamespace("IRTsimrel", quietly = TRUE)) {
    result <- .simulate_irtsimrel(
      n_persons          = n_persons,
      n_items            = n_items,
      model              = model,
      target_rho         = target_rho,
      latent_shape       = latent_shape,
      item_source        = item_source,
      reliability_metric = reliability_metric,
      latent_params      = latent_params,
      item_params        = item_params,
      M                  = M,
      seed               = seed,
      verbose            = verbose
    )
    sim_method <- "irtsimrel"
  } else {
    # --- Fallback: Paganin-style simulation ---
    if (model == "3pl" && use_irtsimrel) {
      message("IRTsimrel does not support 3PL. Using fallback simulation.")
    } else if (use_irtsimrel &&
               !requireNamespace("IRTsimrel", quietly = TRUE)) {
      message("IRTsimrel not installed. Using fallback simulation. ",
              "Install with: remotes::install_github('itemresponsewarehouse/IRTsimrel')")
    }
    result <- .simulate_fallback(n_persons, n_items, model, latent_shape)
    sim_method <- "fallback"
  }

  # --- Build dpmirt_sim object ---
  sim <- structure(
    list(
      response     = result$y,
      theta        = result$theta,
      beta         = result$beta,
      lambda       = result$lambda,
      delta        = result$delta,
      n_persons    = n_persons,
      n_items      = n_items,
      model        = model,
      reliability  = result$reliability,
      target_rho   = if (sim_method == "irtsimrel") target_rho else NULL,
      latent_shape = latent_shape,
      eqc_result   = result$eqc_result,
      method       = sim_method
    ),
    class = "dpmirt_sim"
  )

  sim
}


# ============================================================================
# IRTsimrel Integration
# ============================================================================

#' Simulate data using IRTsimrel (EQC calibration)
#'
#' Internal function that orchestrates the IRTsimrel workflow:
#' eqc_calibrate() -> simulate_response_data()
#'
#' @noRd
.simulate_irtsimrel <- function(n_persons,
                                n_items,
                                model,
                                target_rho,
                                latent_shape,
                                item_source,
                                reliability_metric,
                                latent_params,
                                item_params,
                                M,
                                seed,
                                verbose) {

  if (verbose) message("  Using IRTsimrel for reliability-targeted simulation...")

  # --- Map DPMirt shape names to IRTsimrel ---
  irtsimrel_shape <- .map_latent_shape(latent_shape)

  # --- Step 1: EQC Calibration ---
  if (verbose) message("  Running EQC calibration (target rho = ", target_rho, ")...")

  eqc_args <- list(
    target_rho         = target_rho,
    n_items            = n_items,
    model              = model,
    latent_shape       = irtsimrel_shape,
    item_source        = item_source,
    reliability_metric = reliability_metric,
    M                  = M,
    seed               = seed,
    verbose            = verbose
  )

  # Pass through latent_params and item_params
  if (length(latent_params) > 0) {
    eqc_args$latent_params <- latent_params
  }
  if (length(item_params) > 0) {
    eqc_args$item_params <- item_params
  }

  eqc_result <- do.call(IRTsimrel::eqc_calibrate, eqc_args)

  if (verbose) {
    message("  EQC result: c* = ", round(eqc_result$c_star, 4),
            ", achieved rho = ", round(eqc_result$achieved_rho, 4))
  }

  # --- Step 2: Generate response data ---
  if (verbose) message("  Generating response data (N = ", n_persons, ")...")

  sim_args <- list(
    eqc_result    = eqc_result,
    n_persons     = n_persons,
    latent_shape  = irtsimrel_shape,
    seed          = if (!is.null(seed)) seed + 1000L else NULL
  )

  # Pass through latent_params for data generation
  if (length(latent_params) > 0) {
    sim_args$latent_params <- latent_params
  }

  sim_data <- do.call(IRTsimrel::simulate_response_data, sim_args)

  # --- Compute empirical reliability (KR-20) ---
  y <- sim_data$response_matrix
  storage.mode(y) <- "double"
  reliability <- .compute_kr20(y)

  if (verbose) {
    message("  Empirical KR-20 reliability: ", round(reliability, 4))
  }

  # Rasch: lambda should be NULL (all discriminations = 1 in Rasch)
  lambda_out <- if (model == "2pl") sim_data$lambda else NULL

  list(
    y           = y,
    theta       = sim_data$theta,
    beta        = sim_data$beta,
    lambda      = lambda_out,
    delta       = NULL,   # IRTsimrel does not support 3PL
    reliability = reliability,
    eqc_result  = eqc_result
  )
}


#' Map DPMirt latent shape names to IRTsimrel shape names
#'
#' Handles backward compatibility: DPMirt's "skewed" maps to
#' IRTsimrel's "skew_pos". All other names pass through directly.
#'
#' @param shape Character. DPMirt shape name.
#' @return Character. IRTsimrel shape name.
#' @noRd
.map_latent_shape <- function(shape) {
  # DPMirt's legacy "skewed" -> IRTsimrel's "skew_pos"
  if (shape == "skewed") {
    return("skew_pos")
  }

  # All IRTsimrel shapes are valid (validated by IRTsimrel itself)
  valid_irtsimrel <- c("normal", "bimodal", "trimodal", "multimodal",
                       "skew_pos", "skew_neg", "heavy_tail",
                       "light_tail", "uniform", "floor", "ceiling",
                       "custom")
  valid_fallback <- c("normal", "bimodal", "skewed")

  if (!(shape %in% c(valid_irtsimrel, valid_fallback))) {
    stop("Unknown latent_shape: '", shape, "'. ",
         "Valid shapes: ", paste(valid_irtsimrel, collapse = ", "),
         call. = FALSE)
  }

  shape
}


# ============================================================================
# Fallback Simulation (Paganin-style, no reliability targeting)
# ============================================================================

#' Fallback simulation without IRTsimrel
#'
#' Generates IRT data without reliability targeting.
#' Supports: normal, bimodal, skewed latent distributions.
#' Item difficulties: evenly spaced on \[-2, 2\].
#' Item discriminations (2PL): log-normal.
#'
#' @noRd
.simulate_fallback <- function(n_persons, n_items, model, latent_shape) {

  # --- Generate latent abilities ---
  theta <- switch(latent_shape,
    "normal"  = rnorm(n_persons, 0, 1),
    "bimodal" = {
      # 50/50 mixture of N(-1.5, 0.5) and N(1.5, 0.5)
      groups <- sample(1:2, n_persons, replace = TRUE)
      ifelse(groups == 1,
             rnorm(n_persons, -1.5, sqrt(0.5)),
             rnorm(n_persons,  1.5, sqrt(0.5)))
    },
    "skewed"  = {
      # Right-skewed: shifted exponential
      rexp(n_persons, rate = 1) - 1
    },
    rnorm(n_persons, 0, 1)
  )

  # --- Generate item parameters ---
  beta <- seq(-2, 2, length.out = n_items)  # Evenly spaced difficulties

  lambda <- NULL
  delta <- NULL
  if (model %in% c("2pl", "3pl")) {
    lambda <- exp(rnorm(n_items, 0.5, 0.3))  # Log-normal discrimination
  }
  if (model == "3pl") {
    delta <- rbeta(n_items, 4, 12)  # Guessing parameters ~ Beta(4, 12)
  }

  # --- Generate responses ---
  y <- .generate_responses(theta, beta, lambda, model, delta)

  # --- Compute empirical reliability (KR-20) ---
  reliability <- .compute_kr20(y)

  list(y = y, theta = theta, beta = beta, lambda = lambda,
       delta = delta, reliability = reliability, eqc_result = NULL)
}


# ============================================================================
# Shared Internal Helpers
# ============================================================================

#' Generate binary response matrix from IRT parameters
#'
#' Supports Rasch, 2PL, and 3PL models.
#' 3PL ICC: pi = delta + (1 - delta) * logistic(lambda * (theta - beta))
#'
#' @noRd
.generate_responses <- function(theta, beta, lambda, model, delta = NULL) {
  n_persons <- length(theta)
  n_items <- length(beta)

  y <- matrix(NA_real_, nrow = n_persons, ncol = n_items)
  for (i in seq_len(n_items)) {
    if (model == "rasch") {
      logit_p <- theta - beta[i]
    } else {
      logit_p <- lambda[i] * (theta - beta[i])
    }
    prob <- 1 / (1 + exp(-logit_p))

    # 3PL: apply guessing parameter
    if (model == "3pl" && !is.null(delta)) {
      prob <- delta[i] + (1 - delta[i]) * prob
    }

    y[, i] <- rbinom(n_persons, 1, prob)
  }

  storage.mode(y) <- "double"
  y
}


#' Compute KR-20 reliability
#' @noRd
.compute_kr20 <- function(y) {
  n_items <- ncol(y)
  p_vals <- colMeans(y)
  item_var <- p_vals * (1 - p_vals)
  total_var <- var(rowSums(y))

  if (total_var < .Machine$double.eps) {
    return(0)
  }

  (n_items / (n_items - 1)) * (1 - sum(item_var) / total_var)
}


# ============================================================================
# S3 Methods
# ============================================================================

#' @rdname dpmirt_simulate
#' @param x A \code{dpmirt_sim} object.
#' @param ... Additional arguments (currently unused).
#' @export
print.dpmirt_sim <- function(x, ...) {
  cat("DPMirt Simulated Data\n")
  cat("=====================\n")
  cat("Model:        ", toupper(x$model), "\n")
  cat("Persons:      ", x$n_persons, "\n")
  cat("Items:        ", x$n_items, "\n")
  cat("Distribution: ", x$latent_shape, "\n")
  cat("Method:       ", x$method, "\n")

  if (!is.null(x$target_rho)) {
    cat("Target rho:   ", round(x$target_rho, 3), "\n")
  }
  cat("KR-20:        ", round(x$reliability, 3), "\n")

  if (!is.null(x$eqc_result)) {
    cat("EQC c*:       ", round(x$eqc_result$c_star, 4), "\n")
    cat("EQC rho:      ", round(x$eqc_result$achieved_rho, 4), "\n")
  }

  invisible(x)
}
