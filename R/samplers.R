# ============================================================================
# Module 11: Custom NIMBLE Components
# ============================================================================
# Purpose: Define custom NIMBLE samplers and distributions.
# Blueprint: Section 6 - Module 11, Section 7.5, 7.6
#
# Components adapted from Paganin et al. (2023):
#   - logProb_summer:  Pseudo-sampler for log-probability monitoring
#   - sampler_centered: Joint adaptive MH for (log_lambda, gamma) [SI param]
#   - dBernoulliVector/rBernoulliVector: Vectorized Bernoulli distribution
#
# DESIGN PATTERN:
# We use lazy initialization — NIMBLE components are created on first use
# via .get_*() functions, cached in a package-level environment, and
# passed as objects (not strings) to avoid namespace lookup issues.
# ============================================================================


# Package-level cache environment for NIMBLE components
.nimble_cache <- new.env(parent = emptyenv())


# ============================================================================
# logProb_summer: Log-Probability Monitoring Sampler
# ============================================================================

#' Get or create the logProb_summer sampler
#'
#' Lazily creates the nimbleFunction on first call, then caches it.
#' This avoids executing nimbleFunction() at package load time while
#' ensuring the sampler is available when needed.
#'
#' @return A nimbleFunction class (sampler definition)
#' @noRd
.get_logProb_summer <- function() {
  if (is.null(.nimble_cache$logProb_summer)) {
    .nimble_cache$logProb_summer <- nimbleFunction(
      name = "logProb_summer",
      contains = sampler_BASE,
      setup = function(model, mvSaved, target, control) {
        if (!is.null(control$nodeList)) {
          nodes <- model$expandNodeNames(control$nodeList)
        } else {
          nodes <- model$getNodeNames()
        }
        nodes <- nodes[nodes != target]
      },
      run = function() {
        model[[target]] <<- model$getLogProb(nodes)
        copy(from = model, to = mvSaved,
             row = 1, nodes = target, logProb = TRUE)
      },
      methods = list(
        reset = function() {}
      )
    )
  }
  .nimble_cache$logProb_summer
}


# ============================================================================
# dBernoulliVector: Vectorized Bernoulli Distribution
# ============================================================================
# Used for constrained_item identification where the likelihood for a
# person's full response vector is computed as:
#   y[j, 1:I] ~ dBernoulliVector(pi[j, 1:I])
#
# This allows vectorized item constraints (mean-center beta, etc.) because
# the logit link can be computed on the entire vector at once.
# ============================================================================

#' Register dBernoulliVector custom distribution with NIMBLE
#'
#' Creates and registers the vectorized Bernoulli distribution functions.
#' Must be called before building a NIMBLE model that uses dBernoulliVector.
#'
#' @noRd
.register_dBernoulliVector <- function() {
  if (is.null(.nimble_cache$dBernoulliVector_registered)) {

    # Density function
    .nimble_cache$dBernoulliVector <- nimbleFunction(
      run = function(x    = double(1),
                     prob = double(1),
                     log  = integer(0, default = 0)) {
        returnType(double(0))
        logProb <- sum(dbinom(x, size = 1, prob = prob, log = TRUE))
        if (log) return(logProb) else return(exp(logProb))
      }
    )

    # Random generation function
    .nimble_cache$rBernoulliVector <- nimbleFunction(
      run = function(n    = integer(0),
                     prob = double(1)) {
        returnType(double(1))
        n <- length(prob)
        return(rbinom(n, size = 1, prob = prob))
      }
    )

    # Register with NIMBLE's distribution system
    registerDistributions(list(
      dBernoulliVector = list(
        BUGSdist = "dBernoulliVector(prob)",
        types    = c("value = double(1)", "prob = double(1)"),
        discrete = TRUE
      )
    ))

    .nimble_cache$dBernoulliVector_registered <- TRUE
  }
  invisible(TRUE)
}


# ============================================================================
# sampler_centered: Joint Adaptive MH for SI Parameterization
# ============================================================================
# Reduces posterior correlation between log_lambda[i] and gamma[i] by
# centering gamma proposals around mean(eta). This is Paganin et al.'s
# (2023) centered sampler, critical for efficient mixing in SI parameterization.
#
# Proposal mechanism for each item i:
#   1. Propose: log_lambda'[i] = log_lambda[i] + epsilon,  epsilon ~ N(0, scale^2)
#   2. Center:  gamma'[i] = gamma[i] + mean(eta) * (exp(log_lambda[i]) - exp(log_lambda'[i]))
#   3. MH accept/reject
#
# The centering mean is updated each iteration to mean(eta) across all persons.
# Adaptive scaling targets acceptance rate = 0.44 (optimal for univariate MH).
# ============================================================================

#' Get or create the centered sampler components
#'
#' @return A list with $single (per-item sampler) and $multi (coordinator)
#' @noRd
.get_centered_sampler <- function() {
  if (is.null(.nimble_cache$sampler_centered)) {

    # --- Virtual base for single-item centered sampler ---
    sampler_centered_single_BASE <- nimbleFunctionVirtual(
      methods = list(
        reset    = function() {},
        set_mean = function(m = double()) {}
      )
    )
    .nimble_cache$sampler_centered_single_BASE <- sampler_centered_single_BASE

    # --- Single-item centered sampler ---
    .nimble_cache$sampler_centered_single <- nimbleFunction(
      name     = "centered_single",
      contains = sampler_centered_single_BASE,
      setup = function(model, mvSaved, target, control) {
        centering_mean <- 0
        scale         <- if (!is.null(control$scale)) control$scale else 1
        adaptive      <- if (!is.null(control$adaptive)) control$adaptive else TRUE
        adaptInterval <- if (!is.null(control$adaptInterval)) control$adaptInterval else 200

        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        if (length(targetAsScalar) != 2) {
          stop("centered sampler requires exactly two nodes")
        }

        target1 <- targetAsScalar[1]   # log_lambda[i] or logLambda.tmp[i]
        target2 <- targetAsScalar[2]   # gamma[i] or gamma.tmp[i]
        calcNodes       <- model$getDependencies(target)
        calcNodesNoSelf <- model$getDependencies(target, self = FALSE)

        scaleOriginal <- scale
        timesRan      <- 0
        timesAccepted <- 0
        timesAdapted  <- 0
        scaleHistory       <- c(0, 0)
        acceptanceHistory  <- c(0, 0)
        saveMCMChistory <- nimbleOptions("MCMCsaveHistory")
        optimalAR     <- 0.44
        gamma1        <- 0

        if (any(model$isDiscrete(targetAsScalar))) {
          stop("centered sampler cannot be used on discrete-valued target")
        }
      },
      run = function() {
        currentValue1 <- model[[target1]]   # current log_lambda
        currentValue2 <- model[[target2]]   # current gamma

        # Propose new log_lambda
        propLogScale <- rnorm(1, mean = 0, sd = scale)
        propValue1 <- currentValue1 + propLogScale

        # Center gamma: gamma' = gamma + mean(eta)*(lambda - lambda')
        propValue2 <- currentValue2 +
          centering_mean * (exp(currentValue1) - exp(propValue1))

        model[[target1]] <<- propValue1
        model[[target2]] <<- propValue2

        logMHR <- calculateDiff(model, target)
        if (logMHR == -Inf) {
          nimCopy(from = mvSaved, to = model,
                  row = 1, nodes = target, logProb = TRUE)
          jump <- FALSE
        } else {
          logMHR <- logMHR + calculateDiff(model, calcNodesNoSelf)
          jump <- decide(logMHR)
          if (jump) {
            nimCopy(from = model, to = mvSaved,
                    row = 1, nodes = calcNodes, logProb = TRUE)
          } else {
            nimCopy(from = mvSaved, to = model,
                    row = 1, nodes = calcNodes, logProb = TRUE)
          }
        }
        if (adaptive) adaptiveProcedure(jump)
      },
      methods = list(
        set_mean = function(m = double()) {
          centering_mean <<- m
        },
        adaptiveProcedure = function(jump = logical()) {
          timesRan <<- timesRan + 1
          if (jump) timesAccepted <<- timesAccepted + 1
          if (timesRan %% adaptInterval == 0) {
            acceptanceRate <- timesAccepted / timesRan
            timesAdapted <<- timesAdapted + 1
            gamma1 <<- 1 / ((timesAdapted + 3)^0.8)
            gamma2 <- 10 * gamma1
            adaptFactor <- exp(gamma2 * (acceptanceRate - optimalAR))
            scale <<- scale * adaptFactor
            timesRan <<- 0
            timesAccepted <<- 0
          }
        },
        reset = function() {
          scale         <<- scaleOriginal
          timesRan      <<- 0
          timesAccepted <<- 0
          timesAdapted  <<- 0
          gamma1        <<- 0
        }
      )
    )

    # sampler_centered_single is a local variable in this function scope.
    # The coordinator nimbleFunction defined below can resolve it via
    # lexical scoping (enclosing environment).
    sampler_centered_single <- .nimble_cache$sampler_centered_single

    # --- Multi-item coordinator sampler ---
    .nimble_cache$sampler_centered <- nimbleFunction(
      name     = "centered",
      contains = sampler_BASE,
      setup = function(model, mvSaved, target, control) {
        nodesToCenter <- control$nodesToCenter   # "eta"
        scale <- if (!is.null(control$scale)) control$scale else 1

        targetNodes <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        numPairs <- length(targetNodes) / 2
        samplers <- nimbleFunctionList(sampler_centered_single_BASE)

        for (i in 1:numPairs) {
          samplers[[i]] <- sampler_centered_single(
            model, mvSaved,
            c(targetNodes[i], targetNodes[i + numPairs]),
            control = list(scale = scale)
          )
        }
      },
      run = function() {
        mean_covariate <- mean(model[[nodesToCenter]])
        for (i in 1:numPairs) {
          samplers[[i]]$set_mean(mean_covariate)
          samplers[[i]]$run()
        }
      },
      methods = list(
        reset = function() {
          for (i in 1:numPairs) {
            samplers[[i]]$reset()
          }
        }
      )
    )
  }

  list(
    single = .nimble_cache$sampler_centered_single,
    multi  = .nimble_cache$sampler_centered,
    BASE   = .nimble_cache$sampler_centered_single_BASE
  )
}
