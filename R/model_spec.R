# ============================================================================
# Module 1: Model Specification
# ============================================================================
# Purpose: Build NIMBLE model code, constants, data, and initial values.
# Blueprint: Section 6 - Module 1
# ============================================================================

#' Create a DPMirt model specification
#'
#' Constructs a complete specification for a Bayesian IRT model, including
#' NIMBLE model code, constants, data, initial values, and monitor
#' configuration. This is the first step in the step-by-step workflow.
#'
#' @param data A matrix or data.frame of binary (0/1) responses with persons
#'   in rows and items in columns. Can also be a long-format data.frame with
#'   columns for person, item, and response.
#' @param model Character. IRT model type: \code{"rasch"}, \code{"2pl"},
#'   or \code{"3pl"}.
#' @param prior Character. Latent trait prior: \code{"normal"} (parametric)
#'   or \code{"dpm"} (Dirichlet Process Mixture).
#' @param parameterization Character. \code{"irt"} for standard IRT or
#'   \code{"si"} for slope-intercept. Only relevant for 2PL/3PL.
#' @param identification Character or NULL. Identification strategy:
#'   \code{"constrained_item"}, \code{"constrained_ability"},
#'   or \code{"unconstrained"}. If NULL, uses model-specific default
#'   (constrained_item for Rasch, unconstrained for 2PL/3PL).
#' @param alpha_prior Alpha hyperprior specification. NULL for default
#'   (auto-elicit or Gamma(1,3)), a numeric vector c(a, b) for
#'   Gamma(a, b), or a DPprior_fit object. Only used when prior = "dpm".
#' @param base_measure List with DPM base measure hyperparameters:
#'   s2_mu, nu1, nu2. Defaults from Paganin et al. (2023).
#' @param item_priors List of custom item priors (advanced use).
#' @param M Integer. Maximum number of clusters for CRP truncation.
#'   Only used when prior = "dpm".
#' @param data_format Character. \code{"auto"}, \code{"matrix"}, or
#'   \code{"long"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A \code{dpmirt_spec} S3 object containing:
#' \describe{
#'   \item{code}{A \code{nimbleCode} object.}
#'   \item{constants}{List of constants (N, I, M, etc.).}
#'   \item{data}{List with the response data.}
#'   \item{inits}{List of initial values.}
#'   \item{monitors}{Character vector of parameters to track.}
#'   \item{monitors2}{Character vector of parameters for thinned monitoring.}
#'   \item{config}{List of all model configuration options.}
#' }
#'
#' @examples
#' \dontrun{
#' sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
#'
#' # Rasch-Normal specification
#' spec <- dpmirt_spec(sim$response, model = "rasch", prior = "normal")
#' print(spec)
#'
#' # Rasch-DPM with custom alpha prior
#' spec_dpm <- dpmirt_spec(sim$response, model = "rasch", prior = "dpm",
#'                         alpha_prior = c(1, 3))
#'
#' # 2PL specification
#' sim2 <- dpmirt_simulate(300, 25, model = "2pl", seed = 42)
#' spec_2pl <- dpmirt_spec(sim2$response, model = "2pl", prior = "normal",
#'                         parameterization = "irt")
#' }
#'
#' @family model fitting
#' @seealso \code{\link{dpmirt}}, \code{\link{dpmirt_compile}}
#'
#' @export
dpmirt_spec <- function(data,
                        model = c("rasch", "2pl", "3pl"),
                        prior = c("normal", "dpm"),
                        parameterization = c("irt", "si"),
                        identification = NULL,
                        alpha_prior = NULL,
                        base_measure = list(s2_mu = 2,
                                            nu1 = 2.01,
                                            nu2 = 1.01),
                        item_priors = list(),
                        M = 50L,
                        data_format = c("auto", "matrix", "long"),
                        ...) {

  # --- Argument matching and validation ---
  model <- .validate_model(match.arg(model))
  prior <- .validate_prior(match.arg(prior))
  parameterization <- .validate_parameterization(match.arg(parameterization),
                                                  model)
  identification <- .resolve_identification(identification, model, prior)
  data_format <- match.arg(data_format)

  # Validate the full combination
  .validate_model_combination(model, prior, parameterization, identification)

  # --- Prepare data ---
  data_info <- .validate_data(data, data_format)
  N <- data_info$N
  I <- data_info$I
  y <- data_info$y

  # --- Resolve alpha prior (DPM only) ---
  alpha_ab <- NULL
  if (prior == "dpm") {
    alpha_ab <- .resolve_alpha_prior(alpha_prior, N)
    M <- as.integer(M)
  }

  # --- Build NIMBLE code ---
  nimble_code <- .build_nimble_code(model, prior, parameterization,
                                     identification)

  # --- Build constants ---
  constants <- .build_constants(N, I, M, prior, alpha_ab, base_measure,
                                 model, parameterization)

  # --- Build data list ---
  nimble_data <- list(y = y)

  # --- Generate initial values ---
  inits <- .generate_inits(y, model, prior, parameterization,
                            identification, N, I, M)

  # --- Configure monitors ---
  mon <- .setup_monitors(model, prior, parameterization, identification)

  # --- Build config record ---
  config <- list(
    model = model,
    prior = prior,
    parameterization = parameterization,
    identification = identification,
    alpha_prior = alpha_ab,
    base_measure = base_measure,
    M = M,
    N = N,
    I = I,
    item_priors = item_priors,
    data_format = data_info$data_format
  )

  # --- Construct and return spec object ---
  spec <- structure(
    list(
      code      = nimble_code,
      constants = constants,
      data      = nimble_data,
      inits     = inits,
      monitors  = mon$monitors,
      monitors2 = mon$monitors2,
      config    = config
    ),
    class = "dpmirt_spec"
  )

  spec
}


# ============================================================================
# NIMBLE Code Generation (Programmatic)
# ============================================================================

#' Build NIMBLE code for specified model configuration
#' @noRd
.build_nimble_code <- function(model, prior, parameterization,
                                identification) {
  code_builder <- switch(model,
    "rasch" = .build_rasch_code,
    "2pl"   = .build_2pl_code,
    "3pl"   = .build_3pl_code
  )
  code_builder(prior, parameterization, identification)
}


# --------------------------------------------------------------------------
# Rasch Model Code Builders
# --------------------------------------------------------------------------

#' Build NIMBLE code for Rasch model
#' @noRd
.build_rasch_code <- function(prior, parameterization, identification) {
  if (prior == "normal") {
    .build_rasch_normal_code(identification)
  } else {
    .build_rasch_dpm_code(identification)
  }
}


#' Rasch + Normal prior NIMBLE code
#' @noRd
.build_rasch_normal_code <- function(identification) {

  if (identification == "constrained_item") {
    # Rasch-Normal with constrained item (mean-center beta during MCMC)
    # Uses dBernoulliVector for vectorized likelihood
    code <- nimbleCode({
      # --- Likelihood ---
      for (j in 1:N) {
        for (i in 1:I) {
          y[j, i] ~ dbern(pi[j, i])
          logit(pi[j, i]) <- eta[j] - beta[i]
        }
      }

      # --- Item parameters: mean-centered ---
      for (i in 1:I) {
        beta.tmp[i] ~ dnorm(0, var = sigma2_beta)
      }
      beta[1:I] <- beta.tmp[1:I] - mean(beta.tmp[1:I])

      # --- Latent abilities: Normal prior ---
      for (j in 1:N) {
        eta[j] ~ dnorm(mu, var = s2.eta)
      }
      mu ~ dnorm(0, var = 3)
      s2.eta ~ dinvgamma(2.01, 1.01)

      # --- Log-probability monitoring nodes ---
      myLogProbAll  ~ dnorm(0, 1)
      myLogProbSome ~ dnorm(0, 1)
      myLogLik      ~ dnorm(0, 1)
    })

  } else if (identification == "constrained_ability") {
    # Rasch-Normal with constrained abilities: eta ~ N(0, 1)
    code <- nimbleCode({
      for (j in 1:N) {
        for (i in 1:I) {
          y[j, i] ~ dbern(pi[j, i])
          logit(pi[j, i]) <- eta[j] - beta[i]
        }
      }

      for (i in 1:I) {
        beta[i] ~ dnorm(0, var = sigma2_beta)
      }

      for (j in 1:N) {
        eta[j] ~ dnorm(0, 1)
      }

      myLogProbAll  ~ dnorm(0, 1)
      myLogProbSome ~ dnorm(0, 1)
      myLogLik      ~ dnorm(0, 1)
    })

  } else {
    # Rasch-Normal unconstrained (post-hoc rescaling)
    code <- nimbleCode({
      for (i in 1:I) {
        for (j in 1:N) {
          y[j, i] ~ dbern(pi[j, i])
          logit(pi[j, i]) <- eta[j] - beta[i]
        }
      }

      for (i in 1:I) {
        beta[i] ~ dnorm(0, var = sigma2_beta)
      }

      for (j in 1:N) {
        eta[j] ~ dnorm(mu, var = s2.eta)
      }
      mu ~ dnorm(0, var = 3)
      s2.eta ~ dinvgamma(2.01, 1.01)

      myLogProbAll  ~ dnorm(0, 1)
      myLogProbSome ~ dnorm(0, 1)
      myLogLik      ~ dnorm(0, 1)
    })
  }

  code
}


#' Rasch + DPM prior NIMBLE code
#' @noRd
.build_rasch_dpm_code <- function(identification) {

  if (identification == "constrained_item") {
    # Rasch-DPM with constrained item (mean-center beta)
    # Blueprint Section 3.2, first code block
    code <- nimbleCode({
      # --- Likelihood ---
      for (j in 1:N) {
        for (i in 1:I) {
          y[j, i] ~ dbern(pi[j, i])
          logit(pi[j, i]) <- eta[j] - beta[i]
        }
      }

      # --- Item parameters: mean-centered ---
      for (i in 1:I) {
        beta.tmp[i] ~ dnorm(0, var = sigma2_beta)
      }
      beta[1:I] <- beta.tmp[1:I] - mean(beta.tmp[1:I])

      # --- DPM prior for abilities via CRP ---
      zi[1:N] ~ dCRP(alpha, size = N)
      alpha ~ dgamma(a, b)

      for (j in 1:N) {
        eta[j] ~ dnorm(mu_j[j], var = s2_j[j])
        mu_j[j]  <- muTilde[zi[j]]
        s2_j[j]  <- s2Tilde[zi[j]]
      }

      for (m in 1:M) {
        muTilde[m]  ~ dnorm(0, var = s2_mu)
        s2Tilde[m]  ~ dinvgamma(nu1, nu2)
      }

      # --- Log-probability monitoring ---
      myLogProbAll  ~ dnorm(0, 1)
      myLogProbSome ~ dnorm(0, 1)
      myLogLik      ~ dnorm(0, 1)
    })

  } else {
    # Rasch-DPM unconstrained (post-hoc rescaling)
    # Blueprint Section 3.2, second code block
    code <- nimbleCode({
      for (i in 1:I) {
        for (j in 1:N) {
          y[j, i] ~ dbern(pi[j, i])
          logit(pi[j, i]) <- eta[j] - beta[i]
        }
      }

      for (i in 1:I) {
        beta[i] ~ dnorm(0, var = sigma2_beta)
      }

      zi[1:N] ~ dCRP(alpha, size = N)
      alpha ~ dgamma(a, b)

      for (j in 1:N) {
        eta[j] ~ dnorm(mu_j[j], var = s2_j[j])
        mu_j[j]  <- muTilde[zi[j]]
        s2_j[j]  <- s2Tilde[zi[j]]
      }

      for (m in 1:M) {
        muTilde[m]  ~ dnorm(0, var = s2_mu)
        s2Tilde[m]  ~ dinvgamma(nu1, nu2)
      }

      myLogProbAll  ~ dnorm(0, 1)
      myLogProbSome ~ dnorm(0, 1)
      myLogLik      ~ dnorm(0, 1)
    })
  }

  code
}


# --------------------------------------------------------------------------
# 2PL Model Code Builders
# --------------------------------------------------------------------------

#' Build NIMBLE code for 2PL model
#'
#' Dispatches to IRT or SI builders based on parameterization,
#' then to Normal or DPM builders based on prior.
#' @noRd
.build_2pl_code <- function(prior, parameterization, identification) {
  if (parameterization == "irt") {
    if (prior == "normal") {
      .build_2pl_irt_normal_code(identification)
    } else {
      .build_2pl_irt_dpm_code(identification)
    }
  } else {
    # SI parameterization
    if (prior == "normal") {
      .build_2pl_si_normal_code(identification)
    } else {
      .build_2pl_si_dpm_code(identification)
    }
  }
}


# --- 2PL IRT Parameterization: logit(p) = lambda * (eta - beta) ---

#' 2PL IRT + Normal prior
#' @noRd
.build_2pl_irt_normal_code <- function(identification) {

  if (identification == "unconstrained") {
    code <- nimbleCode({
      for (i in 1:I) {
        for (j in 1:N) {
          y[j, i] ~ dbern(pi[j, i])
          logit(pi[j, i]) <- lambda[i] * (eta[j] - beta[i])
        }
      }
      for (i in 1:I) {
        log(lambda[i]) ~ dnorm(0.5, var = 0.5)
        beta[i] ~ dnorm(0, var = sigma2_beta)
      }
      for (j in 1:N) {
        eta[j] ~ dnorm(mu, var = s2.eta)
      }
      mu ~ dnorm(0, var = 3)
      s2.eta ~ dinvgamma(2.01, 1.01)
      myLogProbAll  ~ dnorm(0, 1)
      myLogProbSome ~ dnorm(0, 1)
      myLogLik      ~ dnorm(0, 1)
    })

  } else if (identification == "constrained_item") {
    # Requires dBernoulliVector for vectorized likelihood
    code <- nimbleCode({
      for (j in 1:N) {
        y[j, 1:I] ~ dBernoulliVector(prob = pi[j, 1:I])
        logit(pi[j, 1:I]) <- lambda[1:I] * (eta[j] - beta[1:I])
      }
      for (i in 1:I) {
        beta.tmp[i] ~ dnorm(0, var = sigma2_beta)
        logLambda.tmp[i] ~ dnorm(0.5, var = 0.5)
      }
      log(lambda[1:I]) <- logLambda.tmp[1:I] - mean(logLambda.tmp[1:I])
      beta[1:I] <- beta.tmp[1:I] - mean(beta.tmp[1:I])
      for (j in 1:N) {
        eta[j] ~ dnorm(mu, var = s2.eta)
      }
      mu ~ dnorm(0, var = 3)
      s2.eta ~ dinvgamma(2.01, 1.01)
      myLogProbAll  ~ dnorm(0, 1)
      myLogProbSome ~ dnorm(0, 1)
      myLogLik      ~ dnorm(0, 1)
    })

  } else {
    # constrained_ability: eta ~ N(0, 1)
    code <- nimbleCode({
      for (i in 1:I) {
        for (j in 1:N) {
          y[j, i] ~ dbern(pi[j, i])
          logit(pi[j, i]) <- lambda[i] * (eta[j] - beta[i])
        }
      }
      for (i in 1:I) {
        log(lambda[i]) ~ dnorm(0.5, var = 0.5)
        beta[i] ~ dnorm(0, var = sigma2_beta)
      }
      for (j in 1:N) {
        eta[j] ~ dnorm(0, 1)
      }
      myLogProbAll  ~ dnorm(0, 1)
      myLogProbSome ~ dnorm(0, 1)
      myLogLik      ~ dnorm(0, 1)
    })
  }
  code
}


#' 2PL IRT + DPM prior
#' @noRd
.build_2pl_irt_dpm_code <- function(identification) {

  if (identification == "unconstrained") {
    code <- nimbleCode({
      for (i in 1:I) {
        for (j in 1:N) {
          y[j, i] ~ dbern(pi[j, i])
          logit(pi[j, i]) <- lambda[i] * (eta[j] - beta[i])
        }
      }
      for (i in 1:I) {
        log(lambda[i]) ~ dnorm(0.5, var = 0.5)
        beta[i] ~ dnorm(0, var = sigma2_beta)
      }
      zi[1:N] ~ dCRP(alpha, size = N)
      alpha ~ dgamma(a, b)
      for (j in 1:N) {
        eta[j] ~ dnorm(mu_j[j], var = s2_j[j])
        mu_j[j]  <- muTilde[zi[j]]
        s2_j[j]  <- s2Tilde[zi[j]]
      }
      for (m in 1:M) {
        muTilde[m]  ~ dnorm(0, var = s2_mu)
        s2Tilde[m]  ~ dinvgamma(nu1, nu2)
      }
      myLogProbAll  ~ dnorm(0, 1)
      myLogProbSome ~ dnorm(0, 1)
      myLogLik      ~ dnorm(0, 1)
    })

  } else {
    # constrained_item + DPM
    code <- nimbleCode({
      for (j in 1:N) {
        y[j, 1:I] ~ dBernoulliVector(prob = pi[j, 1:I])
        logit(pi[j, 1:I]) <- lambda[1:I] * (eta[j] - beta[1:I])
      }
      for (i in 1:I) {
        beta.tmp[i] ~ dnorm(0, var = sigma2_beta)
        logLambda.tmp[i] ~ dnorm(0.5, var = 0.5)
      }
      log(lambda[1:I]) <- logLambda.tmp[1:I] - mean(logLambda.tmp[1:I])
      beta[1:I] <- beta.tmp[1:I] - mean(beta.tmp[1:I])
      zi[1:N] ~ dCRP(alpha, size = N)
      alpha ~ dgamma(a, b)
      for (j in 1:N) {
        eta[j] ~ dnorm(mu_j[j], var = s2_j[j])
        mu_j[j]  <- muTilde[zi[j]]
        s2_j[j]  <- s2Tilde[zi[j]]
      }
      for (m in 1:M) {
        muTilde[m]  ~ dnorm(0, var = s2_mu)
        s2Tilde[m]  ~ dinvgamma(nu1, nu2)
      }
      myLogProbAll  ~ dnorm(0, 1)
      myLogProbSome ~ dnorm(0, 1)
      myLogLik      ~ dnorm(0, 1)
    })
  }
  code
}


# --- 2PL SI Parameterization: logit(p) = lambda * eta + gamma ---

#' 2PL SI + Normal prior
#' @noRd
.build_2pl_si_normal_code <- function(identification) {

  if (identification == "unconstrained") {
    code <- nimbleCode({
      for (i in 1:I) {
        for (j in 1:N) {
          y[j, i] ~ dbern(pi[j, i])
          logit(pi[j, i]) <- lambda[i] * eta[j] + gamma[i]
        }
      }
      for (i in 1:I) {
        log(lambda[i]) ~ dnorm(0.5, var = 0.5)
        gamma[i] ~ dnorm(0, var = sigma2_gamma)
      }
      for (j in 1:N) {
        eta[j] ~ dnorm(mu, var = s2.eta)
      }
      mu ~ dnorm(0, var = 3)
      s2.eta ~ dinvgamma(2.01, 1.01)
      myLogProbAll  ~ dnorm(0, 1)
      myLogProbSome ~ dnorm(0, 1)
      myLogLik      ~ dnorm(0, 1)
    })

  } else if (identification == "constrained_item") {
    code <- nimbleCode({
      for (j in 1:N) {
        y[j, 1:I] ~ dBernoulliVector(prob = pi[j, 1:I])
        logit(pi[j, 1:I]) <- lambda[1:I] * eta[j] + gamma[1:I]
      }
      for (i in 1:I) {
        gamma.tmp[i] ~ dnorm(0, var = sigma2_gamma)
        logLambda.tmp[i] ~ dnorm(0.5, var = 0.5)
      }
      log(lambda[1:I]) <- logLambda.tmp[1:I] - mean(logLambda.tmp[1:I])
      gamma[1:I] <- gamma.tmp[1:I] - mean(gamma.tmp[1:I])
      for (j in 1:N) {
        eta[j] ~ dnorm(mu, var = s2.eta)
      }
      mu ~ dnorm(0, var = 3)
      s2.eta ~ dinvgamma(2.01, 1.01)
      myLogProbAll  ~ dnorm(0, 1)
      myLogProbSome ~ dnorm(0, 1)
      myLogLik      ~ dnorm(0, 1)
    })

  } else {
    # constrained_ability: eta ~ N(0, 1)
    code <- nimbleCode({
      for (i in 1:I) {
        for (j in 1:N) {
          y[j, i] ~ dbern(pi[j, i])
          logit(pi[j, i]) <- lambda[i] * eta[j] + gamma[i]
        }
      }
      for (i in 1:I) {
        log(lambda[i]) ~ dnorm(0.5, var = 0.5)
        gamma[i] ~ dnorm(0, var = sigma2_gamma)
      }
      for (j in 1:N) {
        eta[j] ~ dnorm(0, 1)
      }
      myLogProbAll  ~ dnorm(0, 1)
      myLogProbSome ~ dnorm(0, 1)
      myLogLik      ~ dnorm(0, 1)
    })
  }
  code
}


#' 2PL SI + DPM prior
#' @noRd
.build_2pl_si_dpm_code <- function(identification) {

  if (identification == "unconstrained") {
    code <- nimbleCode({
      for (i in 1:I) {
        for (j in 1:N) {
          y[j, i] ~ dbern(pi[j, i])
          logit(pi[j, i]) <- lambda[i] * eta[j] + gamma[i]
        }
      }
      for (i in 1:I) {
        log(lambda[i]) ~ dnorm(0.5, var = 0.5)
        gamma[i] ~ dnorm(0, var = sigma2_gamma)
      }
      zi[1:N] ~ dCRP(alpha, size = N)
      alpha ~ dgamma(a, b)
      for (j in 1:N) {
        eta[j] ~ dnorm(mu_j[j], var = s2_j[j])
        mu_j[j]  <- muTilde[zi[j]]
        s2_j[j]  <- s2Tilde[zi[j]]
      }
      for (m in 1:M) {
        muTilde[m]  ~ dnorm(0, var = s2_mu)
        s2Tilde[m]  ~ dinvgamma(nu1, nu2)
      }
      myLogProbAll  ~ dnorm(0, 1)
      myLogProbSome ~ dnorm(0, 1)
      myLogLik      ~ dnorm(0, 1)
    })

  } else {
    # constrained_item + DPM + SI
    code <- nimbleCode({
      for (j in 1:N) {
        y[j, 1:I] ~ dBernoulliVector(prob = pi[j, 1:I])
        logit(pi[j, 1:I]) <- lambda[1:I] * eta[j] + gamma[1:I]
      }
      for (i in 1:I) {
        gamma.tmp[i] ~ dnorm(0, var = sigma2_gamma)
        logLambda.tmp[i] ~ dnorm(0.5, var = 0.5)
      }
      log(lambda[1:I]) <- logLambda.tmp[1:I] - mean(logLambda.tmp[1:I])
      gamma[1:I] <- gamma.tmp[1:I] - mean(gamma.tmp[1:I])
      zi[1:N] ~ dCRP(alpha, size = N)
      alpha ~ dgamma(a, b)
      for (j in 1:N) {
        eta[j] ~ dnorm(mu_j[j], var = s2_j[j])
        mu_j[j]  <- muTilde[zi[j]]
        s2_j[j]  <- s2Tilde[zi[j]]
      }
      for (m in 1:M) {
        muTilde[m]  ~ dnorm(0, var = s2_mu)
        s2Tilde[m]  ~ dinvgamma(nu1, nu2)
      }
      myLogProbAll  ~ dnorm(0, 1)
      myLogProbSome ~ dnorm(0, 1)
      myLogLik      ~ dnorm(0, 1)
    })
  }
  code
}


# --------------------------------------------------------------------------
# 3PL Model Code Builders (Phase 6)
# --------------------------------------------------------------------------
# 3PL ICC: pi = delta + (1 - delta) * logistic(linear_predictor)
# where delta ~ Beta(4, 12) is the guessing (lower asymptote) parameter.
# Same structure as 2PL, but with an auxiliary linearReg node.
# --------------------------------------------------------------------------

#' Build NIMBLE code for 3PL model
#'
#' Dispatches to IRT or SI builders based on parameterization,
#' then to Normal or DPM builders based on prior.
#' @noRd
.build_3pl_code <- function(prior, parameterization, identification) {
  if (parameterization == "irt") {
    if (prior == "normal") {
      .build_3pl_irt_normal_code(identification)
    } else {
      .build_3pl_irt_dpm_code(identification)
    }
  } else {
    # SI parameterization
    if (prior == "normal") {
      .build_3pl_si_normal_code(identification)
    } else {
      .build_3pl_si_dpm_code(identification)
    }
  }
}


# --- 3PL IRT Parameterization ---
# logit(linearReg) = lambda * (eta - beta); pi = delta + (1 - delta) * linearReg

#' 3PL IRT + Normal prior
#' @noRd
.build_3pl_irt_normal_code <- function(identification) {

  if (identification == "unconstrained") {
    code <- nimbleCode({
      for (i in 1:I) {
        for (j in 1:N) {
          y[j, i] ~ dbern(pi[j, i])
          pi[j, i] <- delta[i] + (1 - delta[i]) * linearReg[j, i]
          logit(linearReg[j, i]) <- lambda[i] * (eta[j] - beta[i])
        }
      }
      for (i in 1:I) {
        log(lambda[i]) ~ dnorm(0.5, var = 0.5)
        beta[i] ~ dnorm(0, var = sigma2_beta)
        delta[i] ~ dbeta(4, 12)
      }
      for (j in 1:N) {
        eta[j] ~ dnorm(mu, var = s2.eta)
      }
      mu ~ dnorm(0, var = 3)
      s2.eta ~ dinvgamma(2.01, 1.01)
      myLogProbAll  ~ dnorm(0, 1)
      myLogProbSome ~ dnorm(0, 1)
      myLogLik      ~ dnorm(0, 1)
    })

  } else {
    # constrained_ability: eta ~ N(0, 1)
    code <- nimbleCode({
      for (i in 1:I) {
        for (j in 1:N) {
          y[j, i] ~ dbern(pi[j, i])
          pi[j, i] <- delta[i] + (1 - delta[i]) * linearReg[j, i]
          logit(linearReg[j, i]) <- lambda[i] * (eta[j] - beta[i])
        }
      }
      for (i in 1:I) {
        log(lambda[i]) ~ dnorm(0.5, var = 0.5)
        beta[i] ~ dnorm(0, var = sigma2_beta)
        delta[i] ~ dbeta(4, 12)
      }
      for (j in 1:N) {
        eta[j] ~ dnorm(0, 1)
      }
      myLogProbAll  ~ dnorm(0, 1)
      myLogProbSome ~ dnorm(0, 1)
      myLogLik      ~ dnorm(0, 1)
    })
  }
  code
}


#' 3PL IRT + DPM prior
#' @noRd
.build_3pl_irt_dpm_code <- function(identification) {

  # Only unconstrained is valid for 3PL + DPM
  code <- nimbleCode({
    for (i in 1:I) {
      for (j in 1:N) {
        y[j, i] ~ dbern(pi[j, i])
        pi[j, i] <- delta[i] + (1 - delta[i]) * linearReg[j, i]
        logit(linearReg[j, i]) <- lambda[i] * (eta[j] - beta[i])
      }
    }
    for (i in 1:I) {
      log(lambda[i]) ~ dnorm(0.5, var = 0.5)
      beta[i] ~ dnorm(0, var = sigma2_beta)
      delta[i] ~ dbeta(4, 12)
    }
    zi[1:N] ~ dCRP(alpha, size = N)
    alpha ~ dgamma(a, b)
    for (j in 1:N) {
      eta[j] ~ dnorm(mu_j[j], var = s2_j[j])
      mu_j[j]  <- muTilde[zi[j]]
      s2_j[j]  <- s2Tilde[zi[j]]
    }
    for (m in 1:M) {
      muTilde[m]  ~ dnorm(0, var = s2_mu)
      s2Tilde[m]  ~ dinvgamma(nu1, nu2)
    }
    myLogProbAll  ~ dnorm(0, 1)
    myLogProbSome ~ dnorm(0, 1)
    myLogLik      ~ dnorm(0, 1)
  })
  code
}


# --- 3PL SI Parameterization ---
# logit(linearReg) = lambda * eta + gamma; pi = delta + (1 - delta) * linearReg

#' 3PL SI + Normal prior
#' @noRd
.build_3pl_si_normal_code <- function(identification) {

  if (identification == "unconstrained") {
    code <- nimbleCode({
      for (i in 1:I) {
        for (j in 1:N) {
          y[j, i] ~ dbern(pi[j, i])
          pi[j, i] <- delta[i] + (1 - delta[i]) * linearReg[j, i]
          logit(linearReg[j, i]) <- lambda[i] * eta[j] + gamma[i]
        }
      }
      for (i in 1:I) {
        log(lambda[i]) ~ dnorm(0.5, var = 0.5)
        gamma[i] ~ dnorm(0, var = sigma2_gamma)
        delta[i] ~ dbeta(4, 12)
      }
      for (j in 1:N) {
        eta[j] ~ dnorm(mu, var = s2.eta)
      }
      mu ~ dnorm(0, var = 3)
      s2.eta ~ dinvgamma(2.01, 1.01)
      myLogProbAll  ~ dnorm(0, 1)
      myLogProbSome ~ dnorm(0, 1)
      myLogLik      ~ dnorm(0, 1)
    })

  } else {
    # constrained_ability: eta ~ N(0, 1)
    code <- nimbleCode({
      for (i in 1:I) {
        for (j in 1:N) {
          y[j, i] ~ dbern(pi[j, i])
          pi[j, i] <- delta[i] + (1 - delta[i]) * linearReg[j, i]
          logit(linearReg[j, i]) <- lambda[i] * eta[j] + gamma[i]
        }
      }
      for (i in 1:I) {
        log(lambda[i]) ~ dnorm(0.5, var = 0.5)
        gamma[i] ~ dnorm(0, var = sigma2_gamma)
        delta[i] ~ dbeta(4, 12)
      }
      for (j in 1:N) {
        eta[j] ~ dnorm(0, 1)
      }
      myLogProbAll  ~ dnorm(0, 1)
      myLogProbSome ~ dnorm(0, 1)
      myLogLik      ~ dnorm(0, 1)
    })
  }
  code
}


#' 3PL SI + DPM prior
#' @noRd
.build_3pl_si_dpm_code <- function(identification) {

  # Only unconstrained is valid for 3PL + DPM
  code <- nimbleCode({
    for (i in 1:I) {
      for (j in 1:N) {
        y[j, i] ~ dbern(pi[j, i])
        pi[j, i] <- delta[i] + (1 - delta[i]) * linearReg[j, i]
        logit(linearReg[j, i]) <- lambda[i] * eta[j] + gamma[i]
      }
    }
    for (i in 1:I) {
      log(lambda[i]) ~ dnorm(0.5, var = 0.5)
      gamma[i] ~ dnorm(0, var = sigma2_gamma)
      delta[i] ~ dbeta(4, 12)
    }
    zi[1:N] ~ dCRP(alpha, size = N)
    alpha ~ dgamma(a, b)
    for (j in 1:N) {
      eta[j] ~ dnorm(mu_j[j], var = s2_j[j])
      mu_j[j]  <- muTilde[zi[j]]
      s2_j[j]  <- s2Tilde[zi[j]]
    }
    for (m in 1:M) {
      muTilde[m]  ~ dnorm(0, var = s2_mu)
      s2Tilde[m]  ~ dinvgamma(nu1, nu2)
    }
    myLogProbAll  ~ dnorm(0, 1)
    myLogProbSome ~ dnorm(0, 1)
    myLogLik      ~ dnorm(0, 1)
  })
  code
}


# ============================================================================
# Constants and Initialization
# ============================================================================

#' Build constants list for NIMBLE model
#' @noRd
.build_constants <- function(N, I, M, prior, alpha_ab, base_measure,
                              model = "rasch", parameterization = "irt") {
  constants <- list(
    N = N,
    I = I,
    sigma2_beta = 3    # Default beta prior variance (IRT)
  )

  # SI parameterization uses gamma instead of beta
  if (model != "rasch" && parameterization == "si") {
    constants$sigma2_gamma <- 3   # Gamma (intercept) prior variance
  }

  if (prior == "dpm") {
    constants$M     <- M
    constants$a     <- alpha_ab[1]
    constants$b     <- alpha_ab[2]
    constants$s2_mu <- base_measure$s2_mu
    constants$nu1   <- base_measure$nu1
    constants$nu2   <- base_measure$nu2
  }

  constants
}


#' Generate smart initial values
#'
#' Uses standardized sum scores for theta initialization and
#' k-means clustering for DPM cluster assignments (following Paganin).
#'
#' @noRd
.generate_inits <- function(y, model, prior, parameterization,
                             identification, N, I, M) {

  inits <- list()

  # --- Compute standardized sum scores ---
  scores <- rowSums(y, na.rm = TRUE)
  # Guard against zero-variance scores
  if (sd(scores) < .Machine$double.eps * 100) {
    std_scores <- rep(0, N)
  } else {
    std_scores <- (scores - mean(scores)) / sd(scores)
  }

  # --- Initialize person abilities ---
  inits$eta <- std_scores

  # --- Initialize item parameters ---
  is_2pl_or_3pl <- model %in% c("2pl", "3pl")

  if (is_2pl_or_3pl && parameterization == "si") {
    # SI parameterization: gamma (intercept) + lambda (discrimination)
    if (identification == "constrained_item") {
      inits[["gamma.tmp"]]     <- rnorm(I, 0, 1)
      inits[["logLambda.tmp"]] <- runif(I, -1, 1)
    } else {
      inits$gamma <- rnorm(I, 0, 1)
      log_lambda <- runif(I, -1, 1)
      inits$lambda <- exp(log_lambda)
    }
  } else if (is_2pl_or_3pl && parameterization == "irt") {
    # IRT parameterization: beta (difficulty) + lambda (discrimination)
    if (identification == "constrained_item") {
      inits[["beta.tmp"]]      <- rnorm(I, 0, 1)
      inits[["logLambda.tmp"]] <- runif(I, -1, 1)
    } else {
      inits$beta <- rnorm(I, 0, 1)
      log_lambda <- runif(I, -1, 1)
      inits$lambda <- exp(log_lambda)
    }
  } else {
    # Rasch: only beta (no discrimination)
    if (identification == "constrained_item") {
      inits[["beta.tmp"]] <- rnorm(I, 0, 1)
    } else {
      inits$beta <- rnorm(I, 0, 1)
    }
  }

  # --- 3PL: delta (guessing) initialization ---
  if (model == "3pl") {
    inits$delta <- rbeta(I, 4, 12)
  }

  # --- Normal prior hyperparameters ---
  if (prior == "normal" && identification != "constrained_ability") {
    inits$mu     <- 0
    inits$s2.eta <- 1
  }

  # --- DPM-specific initialization ---
  if (prior == "dpm") {
    # K-means for cluster initialization (Paganin pattern)
    n_init_clusters <- min(5, floor(N / 4))
    n_init_clusters <- max(2, n_init_clusters)

    km <- tryCatch(
      kmeans(std_scores, centers = n_init_clusters, nstart = 5),
      error = function(e) {
        # Fallback: assign sequentially
        list(cluster = rep(1:n_init_clusters, length.out = N))
      }
    )
    inits$zi <- km$cluster

    inits$alpha <- 1

    inits$muTilde  <- rnorm(M, 0, 1)
    inits$s2Tilde  <- 1 / rgamma(M, 2.01, 1.01)
  }

  # --- Log-probability monitoring node inits ---
  inits$myLogProbAll  <- 0
  inits$myLogProbSome <- 0
  inits$myLogLik      <- 0

  inits
}


#' Configure monitors for MCMC sampling
#'
#' Sets up primary monitors (sampled every iteration * thin) and
#' secondary monitors (sampled every iteration * thin2, for theta).
#'
#' @noRd
.setup_monitors <- function(model, prior, parameterization, identification) {

  is_2pl_or_3pl <- model %in% c("2pl", "3pl")

  # Primary monitors: item parameters + hyperparameters
  if (is_2pl_or_3pl && parameterization == "si") {
    monitors <- c("gamma", "lambda")
  } else if (is_2pl_or_3pl) {
    # IRT parameterization for 2PL/3PL
    monitors <- c("beta", "lambda")
  } else {
    # Rasch: only beta
    monitors <- c("beta")
  }

  # 3PL: add delta (guessing parameter) to monitors
  if (model == "3pl") {
    monitors <- c(monitors, "delta")
  }

  if (prior == "normal" && identification != "constrained_ability") {
    monitors <- c(monitors, "mu", "s2.eta")
  }

  if (prior == "dpm") {
    monitors <- c(monitors, "alpha", "zi", "muTilde", "s2Tilde")
  }

  # Log-probability monitors
  monitors <- c(monitors, "myLogProbAll", "myLogProbSome", "myLogLik")

  # Secondary monitors (thinned): person abilities
  monitors2 <- "eta"

  list(monitors = monitors, monitors2 = monitors2)
}


# ============================================================================
# Alpha Prior Resolution
# ============================================================================

#' Resolve alpha prior specification
#'
#' Accepts NULL (auto-elicit or default), numeric c(a, b), or DPprior_fit object.
#'
#' @noRd
.resolve_alpha_prior <- function(alpha_prior, N) {

  if (is.null(alpha_prior)) {
    # Default: Paganin Gamma(1, 3)
    return(c(a = 1, b = 3))
  }

  # Check if DPprior_fit object (class-based or duck-typed)
  if (inherits(alpha_prior, "DPprior_fit") ||
      (is.list(alpha_prior) && !is.null(alpha_prior$a) &&
       !is.null(alpha_prior$b))) {
    ap <- c(a = alpha_prior$a, b = alpha_prior$b)
    if (any(ap <= 0)) {
      stop("DPprior_fit object contains non-positive Gamma parameters.",
           call. = FALSE)
    }
    return(ap)
  }

  # Check if numeric vector c(a, b)
  if (is.numeric(alpha_prior) && length(alpha_prior) == 2) {
    names(alpha_prior) <- c("a", "b")
    if (any(alpha_prior <= 0)) {
      stop("Alpha prior parameters must be positive.", call. = FALSE)
    }
    return(alpha_prior)
  }

  stop("Invalid alpha_prior specification. Use NULL, c(a, b), or a ",
       "DPprior_fit object.", call. = FALSE)
}


# ============================================================================
# S3 Methods for dpmirt_spec
# ============================================================================

#' @rdname dpmirt_spec
#' @param x A \code{dpmirt_spec} object.
#' @param ... Additional arguments (currently unused).
#' @export
print.dpmirt_spec <- function(x, ...) {
  cat("DPMirt Model Specification\n")
  cat("==========================\n")
  cat("Model:           ", toupper(x$config$model), "\n")
  cat("Prior:           ", x$config$prior, "\n")
  if (x$config$model != "rasch") {
    cat("Parameterization:", x$config$parameterization, "\n")
  }
  cat("Identification:  ", x$config$identification, "\n")
  cat("Persons (N):     ", x$config$N, "\n")
  cat("Items (I):       ", x$config$I, "\n")
  if (x$config$prior == "dpm") {
    cat("Max clusters (M):", x$config$M, "\n")
    cat("Alpha prior:      Gamma(",
        x$config$alpha_prior[1], ", ",
        x$config$alpha_prior[2], ")\n", sep = "")
  }
  cat("Monitors:        ", paste(x$monitors, collapse = ", "), "\n")
  cat("Monitors2:       ", paste(x$monitors2, collapse = ", "), "\n")
  invisible(x)
}
