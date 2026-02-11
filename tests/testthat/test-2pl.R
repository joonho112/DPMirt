# ============================================================================
# Unit Tests for Phase 3: 2PL Models
# ============================================================================
# Tests 2PL-specific functionality (model spec, rescaling, monitors, etc.)
# These tests do NOT require NIMBLE (pure R logic testing).
# ============================================================================


# --- Helper: Create test data ---
set.seed(42)
y_test <- matrix(rbinom(50 * 10, 1, 0.5), nrow = 50, ncol = 10)


# ============================================================================
# 2PL IRT Model Specification
# ============================================================================

test_that("dpmirt_spec creates 2PL-IRT-Normal-unconstrained spec", {
  spec <- dpmirt_spec(y_test, model = "2pl", prior = "normal",
                       parameterization = "irt")
  expect_s3_class(spec, "dpmirt_spec")
  expect_equal(spec$config$model, "2pl")
  expect_equal(spec$config$parameterization, "irt")
  expect_equal(spec$config$identification, "unconstrained")  # default
})

test_that("2PL IRT spec has correct monitors (beta, lambda)", {
  spec <- dpmirt_spec(y_test, model = "2pl", prior = "normal",
                       parameterization = "irt")
  expect_true("beta" %in% spec$monitors)
  expect_true("lambda" %in% spec$monitors)
  expect_true("mu" %in% spec$monitors)
  expect_true("s2.eta" %in% spec$monitors)
})

test_that("2PL IRT constrained_item spec works", {
  spec <- dpmirt_spec(y_test, model = "2pl", prior = "normal",
                       parameterization = "irt",
                       identification = "constrained_item")
  expect_equal(spec$config$identification, "constrained_item")
})

test_that("2PL IRT constrained_ability spec works", {
  spec <- dpmirt_spec(y_test, model = "2pl", prior = "normal",
                       parameterization = "irt",
                       identification = "constrained_ability")
  expect_equal(spec$config$identification, "constrained_ability")
  # Should not have mu or s2.eta in monitors
  expect_false("mu" %in% spec$monitors)
  expect_false("s2.eta" %in% spec$monitors)
})


# ============================================================================
# 2PL SI Model Specification
# ============================================================================

test_that("dpmirt_spec creates 2PL-SI-Normal-unconstrained spec", {
  spec <- dpmirt_spec(y_test, model = "2pl", prior = "normal",
                       parameterization = "si")
  expect_equal(spec$config$parameterization, "si")
  expect_true("gamma" %in% spec$monitors)
  expect_true("lambda" %in% spec$monitors)
})

test_that("2PL SI constrained_item spec works", {
  spec <- dpmirt_spec(y_test, model = "2pl", prior = "normal",
                       parameterization = "si",
                       identification = "constrained_item")
  expect_equal(spec$config$identification, "constrained_item")
})


# ============================================================================
# 2PL DPM Models
# ============================================================================

test_that("2PL-IRT-DPM-unconstrained spec works", {
  spec <- dpmirt_spec(y_test, model = "2pl", prior = "dpm",
                       parameterization = "irt")
  expect_true("alpha" %in% spec$monitors)
  expect_true("zi" %in% spec$monitors)
  expect_true("lambda" %in% spec$monitors)
  expect_true("beta" %in% spec$monitors)
})

test_that("2PL-SI-DPM-unconstrained spec works", {
  spec <- dpmirt_spec(y_test, model = "2pl", prior = "dpm",
                       parameterization = "si")
  expect_true("alpha" %in% spec$monitors)
  expect_true("gamma" %in% spec$monitors)
  expect_true("lambda" %in% spec$monitors)
})

test_that("2PL-DPM rejects constrained_ability", {
  expect_error(
    dpmirt_spec(y_test, model = "2pl", prior = "dpm",
                 identification = "constrained_ability"),
    "constrained_ability"
  )
})


# ============================================================================
# 2PL Initialization
# ============================================================================

test_that("2PL IRT unconstrained has lambda in inits", {
  spec <- dpmirt_spec(y_test, model = "2pl", prior = "normal",
                       parameterization = "irt")
  expect_true(!is.null(spec$inits$lambda))
  expect_length(spec$inits$lambda, 10)
  expect_true(all(spec$inits$lambda > 0))  # lambda > 0 (log-normal)
})

test_that("2PL IRT constrained_item has logLambda.tmp and beta.tmp", {
  spec <- dpmirt_spec(y_test, model = "2pl", prior = "normal",
                       parameterization = "irt",
                       identification = "constrained_item")
  expect_true(!is.null(spec$inits[["logLambda.tmp"]]))
  expect_true(!is.null(spec$inits[["beta.tmp"]]))
  expect_length(spec$inits[["logLambda.tmp"]], 10)
})

test_that("2PL SI unconstrained has gamma + lambda in inits", {
  spec <- dpmirt_spec(y_test, model = "2pl", prior = "normal",
                       parameterization = "si")
  expect_true(!is.null(spec$inits$gamma))
  expect_true(!is.null(spec$inits$lambda))
})

test_that("2PL SI constrained_item has gamma.tmp + logLambda.tmp in inits", {
  spec <- dpmirt_spec(y_test, model = "2pl", prior = "normal",
                       parameterization = "si",
                       identification = "constrained_item")
  expect_true(!is.null(spec$inits[["gamma.tmp"]]))
  expect_true(!is.null(spec$inits[["logLambda.tmp"]]))
})


# ============================================================================
# 2PL Constants
# ============================================================================

test_that("2PL IRT constants include sigma2_beta", {
  spec <- dpmirt_spec(y_test, model = "2pl", prior = "normal",
                       parameterization = "irt")
  expect_true(!is.null(spec$constants$sigma2_beta))
})

test_that("2PL SI constants include sigma2_gamma", {
  spec <- dpmirt_spec(y_test, model = "2pl", prior = "normal",
                       parameterization = "si")
  expect_true(!is.null(spec$constants$sigma2_gamma))
})


# ============================================================================
# 2PL Rescaling (Pure R Logic)
# ============================================================================

test_that("IRT rescaling centers beta and normalizes lambda", {
  # Simulate fake MCMC output
  niter <- 100
  I <- 5
  N <- 20

  fake_samples <- matrix(rnorm(niter * (I + I + 2 + 3)), nrow = niter)
  colnames(fake_samples) <- c(
    paste0("beta[", 1:I, "]"),
    paste0("lambda[", 1:I, "]"),
    "mu", "s2.eta",
    "myLogProbAll", "myLogProbSome", "myLogLik"
  )
  # Make lambda positive
  fake_samples[, paste0("lambda[", 1:I, "]")] <- exp(fake_samples[, paste0("lambda[", 1:I, "]")])

  fake_samples2 <- matrix(rnorm(niter * N), nrow = niter)
  colnames(fake_samples2) <- paste0("eta[", 1:N, "]")

  fake_obj <- structure(list(
    samples = fake_samples,
    samples2 = fake_samples2,
    model_config = list(model = "2pl", prior = "normal",
                        parameterization = "irt",
                        identification = "unconstrained",
                        I = I, N = N)
  ), class = "dpmirt_samples")

  result <- .rescale_irt(fake_obj)

  # Check: mean(beta*) should be ~0 per iteration
  beta_means <- rowMeans(result$beta_samp)
  expect_true(all(abs(beta_means) < 1e-10))

  # Check: geom_mean(lambda*) should be ~1
  lambda_gm <- apply(result$lambda_samp, 1, function(x) prod(x)^(1/I))
  expect_true(all(abs(lambda_gm - 1) < 1e-10))
})

test_that("SI rescaling normalizes lambda and adjusts gamma", {
  niter <- 100
  I <- 5
  N <- 20

  fake_samples <- matrix(rnorm(niter * (I + I + 2 + 3)), nrow = niter)
  colnames(fake_samples) <- c(
    paste0("gamma[", 1:I, "]"),
    paste0("lambda[", 1:I, "]"),
    "mu", "s2.eta",
    "myLogProbAll", "myLogProbSome", "myLogLik"
  )
  fake_samples[, paste0("lambda[", 1:I, "]")] <- exp(fake_samples[, paste0("lambda[", 1:I, "]")])

  fake_samples2 <- matrix(rnorm(niter * N), nrow = niter)
  colnames(fake_samples2) <- paste0("eta[", 1:N, "]")

  fake_obj <- structure(list(
    samples = fake_samples,
    samples2 = fake_samples2,
    model_config = list(model = "2pl", prior = "normal",
                        parameterization = "si",
                        identification = "unconstrained",
                        I = I, N = N)
  ), class = "dpmirt_samples")

  result <- .rescale_si(fake_obj)

  # Check: geom_mean(lambda*) should be ~1
  lambda_gm <- apply(result$lambda_samp, 1, function(x) prod(x)^(1/I))
  expect_true(all(abs(lambda_gm - 1) < 1e-10))

  # Check: sum(gamma*)/sum(lambda*) should be ~0
  location_check <- sapply(seq_len(niter), function(s) {
    sum(result$beta_samp[s, ]) / sum(result$lambda_samp[s, ])
  })
  expect_true(all(abs(location_check) < 1e-10))
})


# ============================================================================
# Print Methods for 2PL
# ============================================================================

test_that("print.dpmirt_spec works for 2PL", {
  spec <- dpmirt_spec(y_test, model = "2pl", prior = "normal",
                       parameterization = "irt")
  output <- capture.output(print(spec))
  expect_true(any(grepl("2PL", output)))
  expect_true(any(grepl("irt", output, ignore.case = TRUE)))
})

test_that("print.dpmirt_spec works for 2PL-SI-DPM", {
  spec <- dpmirt_spec(y_test, model = "2pl", prior = "dpm",
                       parameterization = "si")
  output <- capture.output(print(spec))
  expect_true(any(grepl("2PL", output)))
  expect_true(any(grepl("dpm", output, ignore.case = TRUE)))
  expect_true(any(grepl("si", output, ignore.case = TRUE)))
})


# ============================================================================
# Rasch Validation (SI rejected)
# ============================================================================

test_that("Rasch + SI parameterization is rejected", {
  expect_error(
    dpmirt_spec(y_test, model = "rasch", parameterization = "si"),
    "Slope-Intercept"
  )
})
