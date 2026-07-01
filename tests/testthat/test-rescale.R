# ============================================================================
# Tests for R/rescale.R
# Rescaling functions for post-hoc identification
# ============================================================================

# ============================================================================
# Helper: create a mock dpmirt_samples object
# ============================================================================

.mock_samples_rasch <- function(N = 20, I = 5, niter = 100, seed = 42) {
  set.seed(seed)

  # Simulate beta and eta with column names matching NIMBLE output
  beta_samp <- matrix(rnorm(niter * I, mean = 0.5), nrow = niter, ncol = I)
  colnames(beta_samp) <- paste0("beta[", seq_len(I), "]")

  eta_samp <- matrix(rnorm(niter * N), nrow = niter, ncol = N)
  colnames(eta_samp) <- paste0("eta[", seq_len(N), "]")

  # Main samples include beta + monitoring nodes
  log_nodes <- matrix(rnorm(niter * 3), nrow = niter, ncol = 3)
  colnames(log_nodes) <- c("myLogProbAll", "myLogProbSome", "myLogLik")

  samples <- cbind(beta_samp, log_nodes)

  structure(
    list(
      samples  = samples,
      samples2 = eta_samp,
      model_config = list(
        model = "rasch", prior = "normal",
        parameterization = "irt",
        identification = "unconstrained",
        N = N, I = I
      )
    ),
    class = "dpmirt_samples"
  )
}


.mock_samples_2pl <- function(N = 20, I = 5, niter = 100, seed = 42) {
  set.seed(seed)

  beta_samp <- matrix(rnorm(niter * I, mean = 0.5), nrow = niter, ncol = I)
  colnames(beta_samp) <- paste0("beta[", seq_len(I), "]")

  lambda_samp <- matrix(exp(rnorm(niter * I, 0.3, 0.2)),
                         nrow = niter, ncol = I)
  colnames(lambda_samp) <- paste0("lambda[", seq_len(I), "]")

  eta_samp <- matrix(rnorm(niter * N), nrow = niter, ncol = N)
  colnames(eta_samp) <- paste0("eta[", seq_len(N), "]")

  log_nodes <- matrix(rnorm(niter * 3), nrow = niter, ncol = 3)
  colnames(log_nodes) <- c("myLogProbAll", "myLogProbSome", "myLogLik")

  samples <- cbind(beta_samp, lambda_samp, log_nodes)

  structure(
    list(
      samples  = samples,
      samples2 = eta_samp,
      model_config = list(
        model = "2pl", prior = "normal",
        parameterization = "irt",
        identification = "unconstrained",
        N = N, I = I
      )
    ),
    class = "dpmirt_samples"
  )
}


.mock_samples_si <- function(N = 20, I = 5, niter = 100, seed = 42) {
  set.seed(seed)

  gamma_samp <- matrix(rnorm(niter * I), nrow = niter, ncol = I)
  colnames(gamma_samp) <- paste0("gamma[", seq_len(I), "]")

  lambda_samp <- matrix(exp(rnorm(niter * I, 0.3, 0.2)),
                         nrow = niter, ncol = I)
  colnames(lambda_samp) <- paste0("lambda[", seq_len(I), "]")

  eta_samp <- matrix(rnorm(niter * N), nrow = niter, ncol = N)
  colnames(eta_samp) <- paste0("eta[", seq_len(N), "]")

  log_nodes <- matrix(rnorm(niter * 3), nrow = niter, ncol = 3)
  colnames(log_nodes) <- c("myLogProbAll", "myLogProbSome", "myLogLik")

  samples <- cbind(gamma_samp, lambda_samp, log_nodes)

  structure(
    list(
      samples  = samples,
      samples2 = eta_samp,
      model_config = list(
        model = "2pl", prior = "normal",
        parameterization = "si",
        identification = "unconstrained",
        N = N, I = I
      )
    ),
    class = "dpmirt_samples"
  )
}


# ============================================================================
# dpmirt_rescale() input validation
# ============================================================================

test_that("dpmirt_rescale rejects non-samples input", {
  expect_error(
    dpmirt_rescale("not a samples object"),
    "dpmirt_samples"
  )
})


test_that("dpmirt_rescale returns no-rescale for constrained models", {
  samp <- .mock_samples_rasch()
  samp$model_config$identification <- "constrained_item"

  result <- dpmirt_rescale(samp)
  # location_shift should all be 0

  expect_true(all(result$location_shift == 0))
  # scale_shift should all be 1
  expect_true(all(result$scale_shift == 1))
})


test_that("dpmirt_rescale passes through when rescale=FALSE", {
  samp <- .mock_samples_rasch()

  result <- dpmirt_rescale(samp, rescale = FALSE)
  beta_cols <- grep("^beta\\[", colnames(samp$samples))
  eta_cols <- grep("^eta\\[", colnames(samp$samples2))

  expect_true(all(result$location_shift == 0))
  expect_true(all(result$scale_shift == 1))
  expect_equal(result$beta_samp, samp$samples[, beta_cols, drop = FALSE])
  expect_equal(result$theta_samp, samp$samples2[, eta_cols, drop = FALSE])
  expect_false(result$config$rescale)
})


# ============================================================================
# Rasch rescaling (.rescale_rasch)
# ============================================================================

test_that(".rescale_rasch centers beta at zero", {
  samp <- .mock_samples_rasch()
  result <- DPMirt:::.rescale_rasch(samp)

  # After rescaling, the mean of beta should be ~0 for each iteration
  beta_means <- rowMeans(result$beta_samp)
  expect_true(all(abs(beta_means) < 1e-10))
})


test_that(".rescale_rasch shifts theta by location shift", {
  samp <- .mock_samples_rasch()
  result <- DPMirt:::.rescale_rasch(samp)

  # location_shift = rowMeans(original beta)
  beta_cols <- grep("^beta\\[", colnames(samp$samples))
  orig_beta_means <- rowMeans(samp$samples[, beta_cols])

  expect_equal(result$location_shift, orig_beta_means, tolerance = 1e-10)

  # theta should be shifted: theta_rescaled = theta_raw - location_shift
  eta_cols <- grep("^eta\\[", colnames(samp$samples2))
  orig_eta <- samp$samples2[, eta_cols]
  expected_theta <- orig_eta - result$location_shift

  expect_equal(as.vector(result$theta_samp), as.vector(expected_theta),
               tolerance = 1e-10)
})


test_that(".rescale_rasch returns no scale shift", {
  samp <- .mock_samples_rasch()
  result <- DPMirt:::.rescale_rasch(samp)

  # Rasch has no scale ambiguity
  expect_true(all(result$scale_shift == 1))
})


test_that(".rescale_rasch returns NULL for lambda_samp", {
  samp <- .mock_samples_rasch()
  result <- DPMirt:::.rescale_rasch(samp)

  expect_null(result$lambda_samp)
})


test_that(".rescale_rasch preserves raw samples", {
  samp <- .mock_samples_rasch()
  result <- DPMirt:::.rescale_rasch(samp)

  expect_identical(result$samples_raw, samp$samples)
  expect_identical(result$samples2_raw, samp$samples2)
})


# ============================================================================
# 2PL IRT rescaling (.rescale_irt)
# ============================================================================

test_that(".rescale_irt centers beta at zero", {
  samp <- .mock_samples_2pl()
  result <- DPMirt:::.rescale_irt(samp)

  beta_means <- rowMeans(result$beta_samp)
  expect_true(all(abs(beta_means) < 1e-10))
})


test_that(".rescale_irt normalizes geometric mean of lambda to 1", {
  samp <- .mock_samples_2pl()
  result <- DPMirt:::.rescale_irt(samp)

  I <- samp$model_config$I
  geom_means <- apply(result$lambda_samp, 1, function(x) prod(x)^(1/I))

  # Geometric mean should be ~1 after rescaling
  expect_true(all(abs(geom_means - 1) < 1e-10))
})


test_that(".rescale_irt applies correct location and scale shifts", {
  samp <- .mock_samples_2pl()
  result <- DPMirt:::.rescale_irt(samp)

  I <- samp$model_config$I

  # Extract original
  beta_cols <- grep("^beta\\[", colnames(samp$samples))
  lambda_cols <- grep("^lambda\\[", colnames(samp$samples))
  orig_beta <- samp$samples[, beta_cols]
  orig_lambda <- samp$samples[, lambda_cols]

  loc <- rowMeans(orig_beta)
  scl <- apply(orig_lambda, 1, function(x) prod(x)^(-1/I))

  # beta* = (beta - loc) / scl
  expected_beta <- (orig_beta - loc) / scl
  expect_equal(as.vector(result$beta_samp), as.vector(expected_beta),
               tolerance = 1e-10)

  # lambda* = lambda * scl
  expected_lambda <- orig_lambda * scl
  expect_equal(as.vector(result$lambda_samp), as.vector(expected_lambda),
               tolerance = 1e-10)
})


test_that(".rescale_irt rescales theta consistently", {
  samp <- .mock_samples_2pl()
  result <- DPMirt:::.rescale_irt(samp)

  I <- samp$model_config$I
  eta_cols <- grep("^eta\\[", colnames(samp$samples2))
  orig_eta <- samp$samples2[, eta_cols]

  # theta* = (theta - loc) / scl
  expected_theta <- (orig_eta - result$location_shift) / result$scale_shift
  expect_equal(as.vector(result$theta_samp), as.vector(expected_theta),
               tolerance = 1e-10)
})


# ============================================================================
# 2PL SI rescaling (.rescale_si)
# ============================================================================

test_that(".rescale_si applies correct SI location shift", {
  samp <- .mock_samples_si()
  result <- DPMirt:::.rescale_si(samp)

  # SI location: sum(gamma) / sum(lambda) per iteration
  gamma_cols <- grep("^gamma\\[", colnames(samp$samples))
  lambda_cols <- grep("^lambda\\[", colnames(samp$samples))
  orig_gamma <- samp$samples[, gamma_cols]
  orig_lambda <- samp$samples[, lambda_cols]

  niter <- nrow(orig_gamma)
  expected_loc <- sapply(seq_len(niter), function(s) {
    sum(orig_gamma[s, ]) / sum(orig_lambda[s, ])
  })

  expect_equal(result$location_shift, expected_loc, tolerance = 1e-10)
})


test_that(".rescale_si normalizes geometric mean of lambda to 1", {
  samp <- .mock_samples_si()
  result <- DPMirt:::.rescale_si(samp)

  I <- samp$model_config$I
  geom_means <- apply(result$lambda_samp, 1, function(x) prod(x)^(1/I))

  expect_true(all(abs(geom_means - 1) < 1e-10))
})


test_that(".rescale_si applies correct sign convention for theta", {
  # SI: theta* = (theta + location) / scale (note: + not -)
  samp <- .mock_samples_si()
  result <- DPMirt:::.rescale_si(samp)

  eta_cols <- grep("^eta\\[", colnames(samp$samples2))
  orig_eta <- samp$samples2[, eta_cols]

  expected_theta <- (orig_eta + result$location_shift) / result$scale_shift
  expect_equal(as.vector(result$theta_samp), as.vector(expected_theta),
               tolerance = 1e-10)
})


# ============================================================================
# No-rescale passthrough (.no_rescale)
# ============================================================================

test_that(".no_rescale extracts beta for IRT parameterization", {
  samp <- .mock_samples_rasch()
  samp$model_config$identification <- "constrained_item"

  result <- DPMirt:::.no_rescale(samp)

  beta_cols <- grep("^beta\\[", colnames(samp$samples))
  expect_equal(ncol(result$beta_samp), length(beta_cols))
})


test_that(".no_rescale extracts gamma for SI parameterization", {
  samp <- .mock_samples_si()
  samp$model_config$identification <- "constrained_item"

  result <- DPMirt:::.no_rescale(samp)

  gamma_cols <- grep("^gamma\\[", colnames(samp$samples))
  expect_equal(ncol(result$beta_samp), length(gamma_cols))
})


test_that(".no_rescale extracts lambda for 2PL", {
  samp <- .mock_samples_2pl()
  samp$model_config$identification <- "constrained_item"

  result <- DPMirt:::.no_rescale(samp)

  lambda_cols <- grep("^lambda\\[", colnames(samp$samples))
  expect_equal(ncol(result$lambda_samp), length(lambda_cols))
})


test_that(".no_rescale returns NULL lambda for Rasch", {
  samp <- .mock_samples_rasch()
  samp$model_config$identification <- "constrained_item"

  result <- DPMirt:::.no_rescale(samp)

  expect_null(result$lambda_samp)
})


# ============================================================================
# Delta extraction helper
# ============================================================================

test_that(".extract_delta_samp returns NULL when no delta columns", {
  samp <- .mock_samples_rasch()
  result <- DPMirt:::.extract_delta_samp(samp$samples)
  expect_null(result)
})


test_that(".extract_delta_samp extracts delta columns when present", {
  set.seed(42)
  niter <- 50; I <- 3
  delta_cols <- matrix(rbeta(niter * I, 2, 8), nrow = niter, ncol = I)
  colnames(delta_cols) <- paste0("delta[", seq_len(I), "]")

  beta_cols <- matrix(rnorm(niter * I), nrow = niter, ncol = I)
  colnames(beta_cols) <- paste0("beta[", seq_len(I), "]")

  samples <- cbind(beta_cols, delta_cols)
  result <- DPMirt:::.extract_delta_samp(samples)

  expect_equal(ncol(result), I)
  expect_equal(nrow(result), niter)
  expect_equal(as.vector(result), as.vector(delta_cols))
})


# ============================================================================
# 3PL rescaling edge cases
# ============================================================================

test_that(".rescale_irt preserves 3PL delta draws unchanged", {
  set.seed(42)
  samp <- .mock_samples_2pl(niter = 50)
  I <- samp$model_config$I

  delta_samp <- matrix(rbeta(50 * I, 4, 12), nrow = 50, ncol = I)
  colnames(delta_samp) <- paste0("delta[", seq_len(I), "]")

  samp$samples <- cbind(samp$samples, delta_samp)
  samp$model_config$model <- "3pl"

  result <- DPMirt:::.rescale_irt(samp)

  expect_equal(dim(result$delta_samp), dim(delta_samp))
  expect_equal(as.vector(result$delta_samp), as.vector(delta_samp))
})


test_that(".rescale_si preserves 3PL delta draws unchanged", {
  set.seed(42)
  samp <- .mock_samples_si(niter = 50)
  I <- samp$model_config$I

  delta_samp <- matrix(rbeta(50 * I, 4, 12), nrow = 50, ncol = I)
  colnames(delta_samp) <- paste0("delta[", seq_len(I), "]")

  samp$samples <- cbind(samp$samples, delta_samp)
  samp$model_config$model <- "3pl"

  result <- DPMirt:::.rescale_si(samp)

  expect_equal(dim(result$delta_samp), dim(delta_samp))
  expect_equal(as.vector(result$delta_samp), as.vector(delta_samp))
})


# ============================================================================
# Thinning alignment
# ============================================================================

.attach_draw_index <- function(samp, main_iters, theta_iters) {
  samp$draw_index <- list(
    main = data.frame(
      row = seq_along(main_iters),
      chain_id = rep(1L, length(main_iters)),
      within_chain_draw = seq_along(main_iters),
      mcmc_iteration = as.integer(main_iters)
    ),
    theta = data.frame(
      row = seq_along(theta_iters),
      chain_id = rep(1L, length(theta_iters)),
      within_chain_draw = seq_along(theta_iters),
      mcmc_iteration = as.integer(theta_iters)
    )
  )
  samp
}


test_that(".rescale_rasch handles different thinning rates", {
  samp <- .mock_samples_rasch(niter = 200)

  # Simulate different thinning: eta has half the rows
  eta_cols <- grep("^eta\\[", colnames(samp$samples2))
  samp$samples2 <- samp$samples2[seq(1, 200, by = 2), , drop = FALSE]

  result <- DPMirt:::.rescale_rasch(samp)

  # Should still produce valid output
  expect_equal(nrow(result$theta_samp), 100)
  expect_equal(nrow(result$beta_samp), 200)
})


test_that(".rescale_irt aligns integer-ratio theta thinning to main draws", {
  samp <- .mock_samples_2pl(niter = 4, seed = 123)

  eta_cols <- grep("^eta\\[", colnames(samp$samples2))
  theta_main_rows <- c(2L, 4L)
  theta_raw <- samp$samples2[theta_main_rows, eta_cols, drop = FALSE]
  samp$samples2 <- theta_raw
  samp <- .attach_draw_index(samp, main_iters = 1:4, theta_iters = c(2L, 4L))

  result <- DPMirt:::.rescale_irt(samp)
  expected_theta <- (theta_raw - result$location_shift[theta_main_rows]) /
    result$scale_shift[theta_main_rows]

  expect_equal(nrow(result$theta_samp), 2)
  expect_equal(as.vector(result$theta_samp), as.vector(expected_theta),
               tolerance = 1e-10)
})


test_that(".rescale_rasch aligns theta shifts by mcmc iteration metadata", {
  samp <- .mock_samples_rasch(niter = 6, seed = 123)
  main_iters <- c(12, 14, 16, 18, 20, 22)
  theta_rows <- c(2L, 4L, 6L)
  theta_iters <- main_iters[theta_rows]

  eta_cols <- grep("^eta\\[", colnames(samp$samples2))
  theta_raw <- samp$samples2[theta_rows, eta_cols, drop = FALSE]
  samp$samples2 <- theta_raw
  samp <- .attach_draw_index(samp, main_iters, theta_iters)

  result <- DPMirt:::.rescale_rasch(samp)
  expected_theta <- theta_raw - result$location_shift[theta_rows]

  expect_equal(as.vector(result$theta_samp), as.vector(expected_theta),
               tolerance = 1e-10)
})


test_that(".rescale_irt aligns theta shifts by mcmc iteration metadata", {
  samp <- .mock_samples_2pl(niter = 6, seed = 123)
  main_iters <- c(12, 14, 16, 18, 20, 22)
  theta_rows <- c(2L, 4L, 6L)
  theta_iters <- main_iters[theta_rows]

  eta_cols <- grep("^eta\\[", colnames(samp$samples2))
  theta_raw <- samp$samples2[theta_rows, eta_cols, drop = FALSE]
  samp$samples2 <- theta_raw
  samp <- .attach_draw_index(samp, main_iters, theta_iters)

  result <- DPMirt:::.rescale_irt(samp)
  expected_theta <- (theta_raw - result$location_shift[theta_rows]) /
    result$scale_shift[theta_rows]

  expect_equal(as.vector(result$theta_samp), as.vector(expected_theta),
               tolerance = 1e-10)
})


test_that(".rescale_irt aligns theta shifts by chain and mcmc iteration", {
  samp <- .mock_samples_2pl(niter = 4, seed = 123)
  main_iters <- c(1L, 2L, 1L, 2L)
  theta_rows <- c(3L, 4L)

  eta_cols <- grep("^eta\\[", colnames(samp$samples2))
  theta_raw <- samp$samples2[theta_rows, eta_cols, drop = FALSE]
  samp$samples2 <- theta_raw
  samp$draw_index <- list(
    main = data.frame(
      row = 1:4,
      chain_id = c(1L, 1L, 2L, 2L),
      within_chain_draw = c(1L, 2L, 1L, 2L),
      mcmc_iteration = main_iters
    ),
    theta = data.frame(
      row = 1:2,
      chain_id = c(2L, 2L),
      within_chain_draw = 1:2,
      mcmc_iteration = c(1L, 2L)
    )
  )

  result <- DPMirt:::.rescale_irt(samp)
  expected_theta <- (theta_raw - result$location_shift[theta_rows]) /
    result$scale_shift[theta_rows]

  expect_equal(as.vector(result$theta_samp), as.vector(expected_theta),
               tolerance = 1e-10)
})


test_that(".rescale_si aligns theta shifts by mcmc iteration metadata", {
  samp <- .mock_samples_si(niter = 6, seed = 123)
  main_iters <- c(12, 14, 16, 18, 20, 22)
  theta_rows <- c(2L, 4L, 6L)
  theta_iters <- main_iters[theta_rows]

  eta_cols <- grep("^eta\\[", colnames(samp$samples2))
  theta_raw <- samp$samples2[theta_rows, eta_cols, drop = FALSE]
  samp$samples2 <- theta_raw
  samp <- .attach_draw_index(samp, main_iters, theta_iters)

  result <- DPMirt:::.rescale_si(samp)
  expected_theta <- (theta_raw + result$location_shift[theta_rows]) /
    result$scale_shift[theta_rows]

  expect_equal(as.vector(result$theta_samp), as.vector(expected_theta),
               tolerance = 1e-10)
})


test_that("rescaling errors when theta iterations are not in main monitors", {
  samp <- .mock_samples_2pl(niter = 6, seed = 123)
  main_iters <- c(12, 14, 16, 18, 20, 22)
  theta_iters <- c(15, 19, 23)

  samp$samples2 <- samp$samples2[1:3, , drop = FALSE]
  samp <- .attach_draw_index(samp, main_iters, theta_iters)

  expect_error(
    DPMirt:::.rescale_irt(samp),
    "Cannot align theta draws"
  )
})


test_that("rescaling errors on duplicate main chain/iteration keys", {
  samp <- .mock_samples_2pl(niter = 4, seed = 123)
  samp$samples2 <- samp$samples2[1:2, , drop = FALSE]
  samp$draw_index <- list(
    main = data.frame(
      row = 1:4,
      chain_id = c(1L, 1L, 1L, 1L),
      within_chain_draw = 1:4,
      mcmc_iteration = c(1L, 2L, 2L, 4L)
    ),
    theta = data.frame(
      row = 1:2,
      chain_id = c(1L, 1L),
      within_chain_draw = 1:2,
      mcmc_iteration = c(2L, 4L)
    )
  )

  expect_error(
    DPMirt:::.rescale_irt(samp),
    "duplicate main draw"
  )
})


test_that("rescaling errors for unsupported non-integer thinning without metadata", {
  samp <- .mock_samples_2pl(niter = 5, seed = 123)
  samp$samples2 <- samp$samples2[1:2, , drop = FALSE]

  expect_error(
    DPMirt:::.rescale_irt(samp),
    "Cannot align theta draws"
  )
})


# ============================================================================
# dpmirt_rescale() returns dpmirt_fit (pipeline integration)
# ============================================================================

test_that("dpmirt_rescale returns dpmirt_fit class", {
  samp <- .mock_samples_rasch()
  result <- dpmirt_rescale(samp)
  expect_s3_class(result, "dpmirt_fit")
})


test_that("dpmirt_rescale result has all dpmirt_fit fields", {
  samp <- .mock_samples_rasch()
  samp$waic <- 123.45
  result <- dpmirt_rescale(samp)

  # Core samples

  expect_true(!is.null(result$theta_samp))
  expect_true(!is.null(result$beta_samp))
  expect_true(!is.null(result$scale_shift))
  expect_true(!is.null(result$location_shift))
  expect_identical(result$samples_raw, samp$samples)
  expect_identical(result$samples2_raw, samp$samples2)
  expect_equal(result$schema_version, "chain-aware-v1")
  expect_equal(result$chain_info$chain_id, 1L)
  expect_equal(result$run_history$chain_id, 1L)
  expect_equal(nrow(result$draw_index$main), nrow(samp$samples))
  expect_equal(nrow(result$draw_index$theta), nrow(samp$samples2))

  # Diagnostics
  expect_true(!is.null(result$ess))
  expect_identical(result$diagnostics$ess$pooled, result$ess)
  expect_null(result$diagnostics$rhat)
  expect_true(is.list(result$ess))
  expect_equal(result$waic, 123.45)
  expect_equal(result$diagnostics$waic$value, 123.45)
  expect_equal(result$diagnostics$waic$aggregation, "single_chain")
  expect_equal(result$diagnostics$waic$by_chain$waic, 123.45)
  expect_length(result$ess$items, samp$model_config$I)
  expect_length(result$ess$theta, samp$model_config$N)
  expect_named(result$ess$items, colnames(result$beta_samp))
  expect_named(result$ess$theta, colnames(result$theta_samp))
  expect_true(all(is.finite(result$ess$items)))
  expect_true(all(is.finite(result$ess$theta)))
  expect_true(all(result$ess$items >= 0))
  expect_true(all(result$ess$theta >= 0))

  # Configuration
  expect_true(!is.null(result$config))
  expect_equal(result$config$model, "rasch")
  expect_equal(result$config$N, 20)
  expect_equal(result$config$I, 5)
})


test_that("dpmirt_rescale result works with dpmirt_estimates", {
  samp <- .mock_samples_rasch()
  fit <- dpmirt_rescale(samp)

  est <- suppressWarnings(dpmirt_estimates(fit, methods = c("pm", "cb", "gr")))
  expect_s3_class(est, "dpmirt_estimates")
  expect_true("theta_pm" %in% names(est$theta))
  expect_true("theta_cb" %in% names(est$theta))
  expect_true("theta_gr" %in% names(est$theta))
})


test_that("step-by-step pipeline: rescale -> estimates -> loss", {
  samp <- .mock_samples_rasch(N = 50, I = 10, niter = 200)
  fit <- dpmirt_rescale(samp)

  est <- suppressWarnings(dpmirt_estimates(fit, methods = c("pm", "cb", "gr")))
  true_theta <- rnorm(50)

  loss <- dpmirt_loss(est, true_theta = true_theta, metrics = c("msel", "ks"))
  expect_s3_class(loss, "data.frame")
  expect_true("msel" %in% names(loss))
  expect_true("ks" %in% names(loss))
  expect_true(nrow(loss) >= 3)  # pm, cb, gr rows
})


test_that("dpmirt_rescale result for 2PL model", {
  samp <- .mock_samples_2pl()
  result <- dpmirt_rescale(samp)

  expect_s3_class(result, "dpmirt_fit")
  expect_true(!is.null(result$lambda_samp))
  expect_equal(result$config$model, "2pl")
})
