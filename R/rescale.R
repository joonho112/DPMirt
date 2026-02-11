# ============================================================================
# Module 4: Rescaling
# ============================================================================
# Purpose: Apply post-hoc identification rescaling to posterior samples.
# Blueprint: Section 6 - Module 4, Section 7.4
#
# Adapted from Paganin et al. (2023) rescalingFunctions.R.
# For Rasch models, only location shift is needed (no scale ambiguity).
# For 2PL/3PL, both location and scale rescaling are applied.
# ============================================================================

#' Rescale posterior samples for identification
#'
#' Applies post-hoc rescaling to MCMC posterior samples to resolve
#' identification indeterminacy. For unconstrained models, this
#' centers item difficulties and (for 2PL/3PL) normalizes
#' discriminations.
#'
#' @param samples_obj A \code{dpmirt_samples} object.
#' @param rescale Logical. If TRUE (default), apply rescaling.
#'   If FALSE, return samples as-is.
#'
#' @return A \code{dpmirt_fit} S3 object containing rescaled posterior
#'   samples, diagnostics (ESS, WAIC, cluster info), and model
#'   configuration.  This object can be passed directly to
#'   \code{\link{dpmirt_estimates}}, \code{\link{dpmirt_resume}}, and
#'   other downstream functions.
#'
#' @details
#' Post-hoc rescaling resolves the identification indeterminacy inherent
#' in unconstrained IRT models:
#'
#' \strong{Rasch} (location only): For each MCMC iteration \eqn{s}:
#' \deqn{\beta^*_i = \beta_i - \bar{\beta}}
#' \deqn{\theta^*_j = \theta_j - \bar{\beta}}
#'
#' \strong{2PL/3PL IRT parameterization} (location + scale): Let
#' \eqn{c_s = \bar{\beta}_s} and \eqn{d_s = (\prod_i \lambda_i)^{-1/I}}:
#' \deqn{\beta^*_i = (\beta_i - c_s) / d_s}
#' \deqn{\lambda^*_i = \lambda_i \cdot d_s}
#' \deqn{\theta^*_j = (\theta_j - c_s) / d_s}
#'
#' After rescaling: \eqn{\bar{\beta}^* = 0} and
#' \eqn{(\prod_i \lambda^*_i)^{1/I} = 1}.
#'
#' @references
#' Paganin, S., Paciorek, C. J., Wehrhahn, C., Rodriguez, A.,
#' Rabe-Hesketh, S., & de Valpine, P. (2023). Computational strategies
#' and estimation performance with Bayesian semiparametric item response
#' theory models. \emph{Journal of Educational and Behavioral Statistics,
#' 48}(2), 147--188.
#'
#' @examples
#' \dontrun{
#' sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
#' spec <- dpmirt_spec(sim$response, model = "rasch", prior = "normal",
#'                     identification = "unconstrained")
#' compiled <- dpmirt_compile(spec)
#' samples <- dpmirt_sample(compiled, niter = 5000, nburnin = 1000)
#'
#' # Apply rescaling
#' rescaled <- dpmirt_rescale(samples)
#' }
#'
#' @family estimation
#' @seealso \code{\link{dpmirt_sample}}, \code{\link{dpmirt_estimates}}
#'
#' @export
dpmirt_rescale <- function(samples_obj, rescale = TRUE) {

  if (!inherits(samples_obj, "dpmirt_samples")) {
    stop("Input must be a dpmirt_samples object.", call. = FALSE)
  }

  config <- samples_obj$model_config

  # No rescaling needed for constrained models
  if (config$identification != "unconstrained" || !rescale) {
    rescaled <- .no_rescale(samples_obj)
  } else if (config$model == "rasch") {
    rescaled <- .rescale_rasch(samples_obj)
  } else if (config$parameterization == "si") {
    rescaled <- .rescale_si(samples_obj)
  } else {
    rescaled <- .rescale_irt(samples_obj)
  }

  # Wrap the plain-list result into a dpmirt_fit S3 object
  .build_fit_from_rescaled(rescaled, samples_obj)
}


# ============================================================================
# Rasch Rescaling (Location Shift Only)
# ============================================================================

#' Rescale Rasch model posterior samples
#'
#' For the Rasch model, only location indeterminacy exists.
#' Rescaling centers beta at zero and shifts eta accordingly.
#'
#' @noRd
.rescale_rasch <- function(samples_obj) {

  samples  <- samples_obj$samples
  samples2 <- samples_obj$samples2
  config   <- samples_obj$model_config

  I <- config$I
  N <- config$N

  # --- Extract beta samples ---
  beta_cols <- grep("^beta\\[", colnames(samples))
  beta_samp <- samples[, beta_cols, drop = FALSE]

  niter <- nrow(beta_samp)

  # --- Extract eta samples ---
  eta_cols <- grep("^eta\\[", colnames(samples2))
  eta_samp <- samples2[, eta_cols, drop = FALSE]

  # Handle potential dimension mismatch (different thinning)
  niter_eta <- nrow(eta_samp)

  # --- Compute location shift per iteration ---
  # location_s = mean(beta_s) for each MCMC iteration
  location_shift <- rowMeans(beta_samp)

  # --- Apply rescaling ---
  # beta*_s  = beta_s - location_s
  # theta*_s = theta_s - location_s
  beta_rescaled <- beta_samp - location_shift

  # Handle different thinning rates for eta
  if (niter_eta == niter) {
    eta_rescaled <- eta_samp - location_shift
    location_shift_eta <- location_shift
  } else {
    # If thin2 differs from thin, we need to align
    # Use interpolated or subsetted location shifts
    thin_ratio <- niter / niter_eta
    if (thin_ratio == round(thin_ratio)) {
      # Clean ratio: subsample location shifts
      idx <- seq(1, niter, by = round(thin_ratio))
      location_shift_eta <- location_shift[idx[1:niter_eta]]
    } else {
      # Approximate: use overall mean shift
      location_shift_eta <- rep(mean(location_shift), niter_eta)
    }
    eta_rescaled <- eta_samp - location_shift_eta
  }

  # --- Extract delta (3PL guessing, no rescaling needed) ---
  delta_samp <- .extract_delta_samp(samples)

  # --- Other parameters (pass through) ---
  skip_cols <- c(beta_cols, grep("^delta\\[", colnames(samples)),
                 grep("^myLog", colnames(samples)))
  other_cols <- setdiff(seq_len(ncol(samples)), skip_cols)
  other_samp <- if (length(other_cols) > 0) {
    samples[, other_cols, drop = FALSE]
  } else {
    NULL
  }

  list(
    theta_samp     = eta_rescaled,
    beta_samp      = beta_rescaled,
    lambda_samp    = NULL,      # Rasch: no discrimination
    delta_samp     = delta_samp,
    other_samp     = other_samp,
    scale_shift    = rep(1, niter),   # No scale shift for Rasch
    location_shift = location_shift,
    samples_raw    = samples,
    samples2_raw   = samples2
  )
}


# ============================================================================
# 2PL/3PL IRT Rescaling (Location + Scale)
# ============================================================================

#' Rescale 2PL/3PL model posterior samples (IRT parameterization)
#'
#' Adapted from Paganin's posteriorRescalingBeta().
#' For each MCMC iteration s:
#'   location_s = mean(beta_s)
#'   scale_s    = prod(lambda_s)^(-1/I)   (inverse geometric mean)
#'   beta*_s    = (beta_s - location_s) / scale_s
#'   lambda*_s  = lambda_s * scale_s
#'   theta*_s   = (theta_s - location_s) / scale_s
#'
#' Invariants: mean(beta*) approx 0, geom_mean(lambda*) approx 1
#'
#' @noRd
.rescale_irt <- function(samples_obj) {

  samples  <- samples_obj$samples
  samples2 <- samples_obj$samples2
  config   <- samples_obj$model_config

  I <- config$I
  N <- config$N

  # --- Extract lambda samples ---
  lambda_cols <- grep("^lambda\\[", colnames(samples))
  lambda_samp <- samples[, lambda_cols, drop = FALSE]

  # --- Extract beta samples ---
  beta_cols <- grep("^beta\\[", colnames(samples))
  beta_samp <- samples[, beta_cols, drop = FALSE]

  niter <- nrow(beta_samp)

  # --- Extract eta samples ---
  eta_cols <- grep("^eta\\[", colnames(samples2))
  eta_samp <- samples2[, eta_cols, drop = FALSE]
  niter_eta <- nrow(eta_samp)

  # --- Compute location and scale shifts per iteration ---
  location_shift <- rowMeans(beta_samp)
  scale_shift <- apply(lambda_samp, 1, function(x) prod(x)^(-1 / I))

  # --- Apply rescaling ---
  beta_rescaled <- (beta_samp - location_shift) / scale_shift
  lambda_rescaled <- lambda_samp * scale_shift

  # Handle different thinning rates for eta
  if (niter_eta == niter) {
    eta_rescaled <- (eta_samp - location_shift) / scale_shift
    location_shift_eta <- location_shift
    scale_shift_eta <- scale_shift
  } else {
    thin_ratio <- niter / niter_eta
    if (thin_ratio == round(thin_ratio)) {
      idx <- seq(1, niter, by = round(thin_ratio))
      location_shift_eta <- location_shift[idx[1:niter_eta]]
      scale_shift_eta    <- scale_shift[idx[1:niter_eta]]
    } else {
      location_shift_eta <- rep(mean(location_shift), niter_eta)
      scale_shift_eta    <- rep(mean(scale_shift), niter_eta)
    }
    eta_rescaled <- (eta_samp - location_shift_eta) / scale_shift_eta
  }

  # --- Extract delta (3PL guessing, no rescaling needed) ---
  delta_samp <- .extract_delta_samp(samples)

  # --- Other parameters (pass through) ---
  skip_cols <- c(beta_cols, lambda_cols,
                 grep("^delta\\[", colnames(samples)),
                 grep("^myLog", colnames(samples)))
  other_cols <- setdiff(seq_len(ncol(samples)), skip_cols)
  other_samp <- if (length(other_cols) > 0) {
    samples[, other_cols, drop = FALSE]
  } else {
    NULL
  }

  list(
    theta_samp     = eta_rescaled,
    beta_samp      = beta_rescaled,
    lambda_samp    = lambda_rescaled,
    delta_samp     = delta_samp,
    other_samp     = other_samp,
    scale_shift    = scale_shift,
    location_shift = location_shift,
    samples_raw    = samples,
    samples2_raw   = samples2
  )
}


#' Rescale 2PL/3PL model posterior samples (SI parameterization)
#'
#' Adapted from Paganin's posteriorRescalingGamma().
#' For each MCMC iteration s:
#'   location_s = sum(gamma_s) / sum(lambda_s)
#'   scale_s    = prod(lambda_s)^(-1/I)
#'   gamma*_s   = gamma_s - lambda_s * location_s
#'   lambda*_s  = lambda_s * scale_s
#'   theta*_s   = (theta_s + location_s) / scale_s
#'
#' Note the sign difference from IRT: theta gets +location (not -location)
#' because gamma = -lambda*beta, so the location shift has opposite sign.
#'
#' @noRd
.rescale_si <- function(samples_obj) {

  samples  <- samples_obj$samples
  samples2 <- samples_obj$samples2
  config   <- samples_obj$model_config

  I <- config$I
  N <- config$N

  # --- Extract lambda samples ---
  lambda_cols <- grep("^lambda\\[", colnames(samples))
  lambda_samp <- samples[, lambda_cols, drop = FALSE]

  # --- Extract gamma samples ---
  gamma_cols <- grep("^gamma\\[", colnames(samples))
  gamma_samp <- samples[, gamma_cols, drop = FALSE]

  niter <- nrow(gamma_samp)

  # --- Extract eta samples ---
  eta_cols <- grep("^eta\\[", colnames(samples2))
  eta_samp <- samples2[, eta_cols, drop = FALSE]
  niter_eta <- nrow(eta_samp)

  # --- Compute location and scale shifts ---
  # SI location: sum(gamma) / sum(lambda) per iteration
  location_shift <- sapply(seq_len(niter), function(s) {
    sum(gamma_samp[s, ]) / sum(lambda_samp[s, ])
  })
  scale_shift <- apply(lambda_samp, 1, function(x) prod(x)^(-1 / I))

  # --- Apply rescaling ---
  gamma_rescaled <- gamma_samp - lambda_samp * location_shift
  lambda_rescaled <- lambda_samp * scale_shift

  # Handle different thinning rates for eta
  if (niter_eta == niter) {
    # SI: theta* = (theta + location) / scale  [note: + not -]
    eta_rescaled <- (eta_samp + location_shift) / scale_shift
    location_shift_eta <- location_shift
    scale_shift_eta <- scale_shift
  } else {
    thin_ratio <- niter / niter_eta
    if (thin_ratio == round(thin_ratio)) {
      idx <- seq(1, niter, by = round(thin_ratio))
      location_shift_eta <- location_shift[idx[1:niter_eta]]
      scale_shift_eta    <- scale_shift[idx[1:niter_eta]]
    } else {
      location_shift_eta <- rep(mean(location_shift), niter_eta)
      scale_shift_eta    <- rep(mean(scale_shift), niter_eta)
    }
    eta_rescaled <- (eta_samp + location_shift_eta) / scale_shift_eta
  }

  # --- Extract delta (3PL guessing, no rescaling needed) ---
  delta_samp <- .extract_delta_samp(samples)

  # --- Other parameters (pass through) ---
  skip_cols <- c(gamma_cols, lambda_cols,
                 grep("^delta\\[", colnames(samples)),
                 grep("^myLog", colnames(samples)))
  other_cols <- setdiff(seq_len(ncol(samples)), skip_cols)
  other_samp <- if (length(other_cols) > 0) {
    samples[, other_cols, drop = FALSE]
  } else {
    NULL
  }

  list(
    theta_samp     = eta_rescaled,
    beta_samp      = gamma_rescaled,  # Store gamma as "beta_samp" for unified interface
    lambda_samp    = lambda_rescaled,
    delta_samp     = delta_samp,
    other_samp     = other_samp,
    scale_shift    = scale_shift,
    location_shift = location_shift,
    samples_raw    = samples,
    samples2_raw   = samples2
  )
}


# ============================================================================
# No-Rescale Passthrough
# ============================================================================

#' Pass through samples without rescaling
#' @noRd
.no_rescale <- function(samples_obj) {

  samples  <- samples_obj$samples
  samples2 <- samples_obj$samples2
  config   <- samples_obj$model_config

  is_si <- !is.null(config$parameterization) && config$parameterization == "si"

  # Extract beta or gamma (for SI)
  if (is_si) {
    # SI: look for gamma[...] or gamma.tmp[...]
    item_cols <- grep("^gamma\\[", colnames(samples))
    beta_samp <- if (length(item_cols) > 0) {
      samples[, item_cols, drop = FALSE]
    } else {
      gamma_tmp_cols <- grep("^gamma\\.tmp\\[", colnames(samples))
      if (length(gamma_tmp_cols) > 0) {
        tmp <- samples[, gamma_tmp_cols, drop = FALSE]
        tmp - rowMeans(tmp)
      } else {
        NULL
      }
    }
  } else {
    item_cols <- grep("^beta\\[", colnames(samples))
    beta_samp <- if (length(item_cols) > 0) {
      samples[, item_cols, drop = FALSE]
    } else {
      beta_tmp_cols <- grep("^beta\\.tmp\\[", colnames(samples))
      if (length(beta_tmp_cols) > 0) {
        tmp <- samples[, beta_tmp_cols, drop = FALSE]
        tmp - rowMeans(tmp)
      } else {
        NULL
      }
    }
  }

  # Extract lambda (for 2PL/3PL)
  lambda_cols <- grep("^lambda\\[", colnames(samples))
  lambda_samp <- if (length(lambda_cols) > 0) {
    samples[, lambda_cols, drop = FALSE]
  } else {
    NULL
  }

  # Extract eta
  eta_cols <- grep("^eta\\[", colnames(samples2))
  eta_samp <- if (!is.null(samples2) && length(eta_cols) > 0) {
    samples2[, eta_cols, drop = FALSE]
  } else {
    NULL
  }

  # Extract delta (3PL guessing, no rescaling needed)
  delta_samp <- .extract_delta_samp(samples)

  # Other parameters (exclude item params, lambda, delta, logprob, eta)
  skip_patterns <- c("^beta\\[", "^beta\\.tmp\\[", "^gamma\\[", "^gamma\\.tmp\\[",
                      "^lambda\\[", "^logLambda\\.tmp\\[", "^delta\\[", "^myLog")
  skip_cols <- unlist(lapply(skip_patterns, function(p) grep(p, colnames(samples))))
  skip_cols <- c(skip_cols, item_cols, lambda_cols)
  skip_cols <- unique(skip_cols)
  other_cols <- setdiff(seq_len(ncol(samples)), skip_cols)
  other_samp <- if (length(other_cols) > 0) {
    samples[, other_cols, drop = FALSE]
  } else {
    NULL
  }

  niter <- nrow(samples)

  list(
    theta_samp     = eta_samp,
    beta_samp      = beta_samp,
    lambda_samp    = lambda_samp,
    delta_samp     = delta_samp,
    other_samp     = other_samp,
    scale_shift    = rep(1, niter),
    location_shift = rep(0, niter),
    samples_raw    = samples,
    samples2_raw   = samples2
  )
}


# ============================================================================
# Delta Extraction Helper
# ============================================================================

#' Extract delta (guessing) samples from raw MCMC output
#'
#' Delta is scale-invariant and does not need rescaling.
#' Returns NULL if delta columns are not present (Rasch/2PL models).
#'
#' @noRd
.extract_delta_samp <- function(samples) {
  delta_cols <- grep("^delta\\[", colnames(samples))
  if (length(delta_cols) > 0) {
    samples[, delta_cols, drop = FALSE]
  } else {
    NULL
  }
}


# ============================================================================
# Build dpmirt_fit from Rescaled Samples
# ============================================================================

#' Assemble a dpmirt_fit S3 object from rescaled samples
#'
#' Takes the plain-list output of an internal rescaling function
#' (`.rescale_rasch`, `.rescale_irt`, `.rescale_si`, `.no_rescale`)
#' and the original `dpmirt_samples` input and constructs a full
#' `dpmirt_fit` object with diagnostics and configuration.
#'
#' @param rescaled Plain list from the internal rescaler.
#' @param samples_obj The original `dpmirt_samples` object.
#' @return A `dpmirt_fit` S3 object.
#' @noRd
.build_fit_from_rescaled <- function(rescaled, samples_obj) {

  config <- samples_obj$model_config

  # Compute diagnostics
  ess <- .compute_ess(rescaled)
  loglik_trace <- .extract_loglik_trace(rescaled$samples_raw)

  # Cluster diagnostics for DPM models
  cluster_info <- NULL
  if (!is.null(config$prior) && config$prior == "dpm") {
    cluster_info <- tryCatch(
      .extract_cluster_info(rescaled$samples_raw, config$N),
      error = function(e) NULL
    )
  }

  # Extract timings (may be NULL for mock objects)
  compilation_time <- NULL
  if (!is.null(samples_obj$compiled)) {
    compilation_time <- samples_obj$compiled$compilation_time
  }

  # Build the MCMC control info from samples_obj if available
  mcmc_ctrl <- samples_obj$mcmc_control
  niter   <- if (!is.null(mcmc_ctrl$niter))   mcmc_ctrl$niter   else nrow(rescaled$beta_samp)
  nburnin <- if (!is.null(mcmc_ctrl$nburnin)) mcmc_ctrl$nburnin else NULL
  thin    <- if (!is.null(mcmc_ctrl$thin))    mcmc_ctrl$thin    else 1L
  thin2   <- if (!is.null(mcmc_ctrl$thin2))   mcmc_ctrl$thin2   else 1L

  structure(
    list(
      # Core rescaled samples
      theta_samp     = rescaled$theta_samp,
      beta_samp      = rescaled$beta_samp,
      lambda_samp    = rescaled$lambda_samp,
      delta_samp     = rescaled$delta_samp,
      other_samp     = rescaled$other_samp,

      # Rescaling info
      scale_shift    = rescaled$scale_shift,
      location_shift = rescaled$location_shift,

      # Raw samples (needed by .combine_chains and dp_density)
      samples_raw    = rescaled$samples_raw,
      samples2_raw   = rescaled$samples2_raw,

      # DPM-specific
      dp_density     = NULL,
      cluster_info   = cluster_info,

      # Model comparison
      waic           = samples_obj$waic,

      # Diagnostics
      ess            = ess,
      loglik_trace   = loglik_trace,

      # Timings
      compilation_time = compilation_time,
      sampling_time    = samples_obj$sampling_time,
      total_time       = NULL,

      # Compiled reference (for resume)
      compiled       = samples_obj$compiled,

      # Configuration
      config         = list(
        model            = config$model,
        prior            = config$prior,
        parameterization = config$parameterization,
        identification   = config$identification,
        N                = config$N,
        I                = config$I,
        niter            = niter,
        nburnin          = nburnin,
        thin             = thin,
        thin2            = thin2,
        nchains          = 1L,
        rescale          = TRUE
      )
    ),
    class = "dpmirt_fit"
  )
}
