# ============================================================================
# Tests for R/estimates.R, R/losses.R, R/draws.R
# Phase 4: Posterior Summarization (PM, CB, GR)
# ============================================================================

# ============================================================================
# .triple_goal() unit tests
# ============================================================================

test_that(".triple_goal returns correct structure", {
  set.seed(42)
  s <- matrix(rnorm(500 * 10), nrow = 500, ncol = 10)

  result <- DPMirt:::.triple_goal(s)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 10)
  expect_true(all(c("theta_pm", "theta_psd", "theta_cb", "theta_gr",
                     "rbar", "rhat") %in% names(result)))
})


test_that(".triple_goal PM equals colMeans", {
  set.seed(42)
  s <- matrix(rnorm(1000 * 20), nrow = 1000, ncol = 20)

  result <- DPMirt:::.triple_goal(s)
  expected_pm <- colMeans(s)

  expect_equal(result$theta_pm, expected_pm, tolerance = 1e-12)
})


test_that(".triple_goal PSD equals column SDs", {
  set.seed(42)
  s <- matrix(rnorm(1000 * 15), nrow = 1000, ncol = 15)

  result <- DPMirt:::.triple_goal(s)
  expected_psd <- apply(s, 2, sd)

  expect_equal(result$theta_psd, expected_psd, tolerance = 1e-12)
})


test_that(".triple_goal CB satisfies mean constraint", {
  set.seed(42)
  s <- matrix(rnorm(2000 * 30), nrow = 2000, ncol = 30)

  result <- DPMirt:::.triple_goal(s)
  etadot <- mean(result$theta_pm)

  # CB constraint: mean(theta_cb) == mean(theta_pm) = etadot
  expect_equal(mean(result$theta_cb), etadot, tolerance = 1e-10)
})


test_that(".triple_goal CB satisfies variance constraint", {
  set.seed(42)
  s <- matrix(rnorm(2000 * 30), nrow = 2000, ncol = 30)

  result <- DPMirt:::.triple_goal(s)

  # CB constraint: var(theta_cb) == mean(lambda_k) + var(theta_pm)
  lambda_k <- result$theta_psd^2
  expected_var <- mean(lambda_k) + var(result$theta_pm)

  expect_equal(var(result$theta_cb), expected_var, tolerance = 1e-10)
})


test_that(".triple_goal GR rhat is a valid permutation", {
  set.seed(42)
  s <- matrix(rnorm(500 * 10), nrow = 500, ncol = 10)

  result <- DPMirt:::.triple_goal(s)

  # rhat should be a permutation of 1:K
  expect_equal(sort(result$rhat), 1:10)
})


test_that(".triple_goal handles near-zero variance gracefully", {
  # Create data where all posterior means are identical
  s <- matrix(0, nrow = 100, ncol = 5)
  s <- s + rnorm(nrow(s) * ncol(s), sd = 0.001)

  expect_warning(
    result <- DPMirt:::.triple_goal(s),
    "Near-zero variance"
  )

  qf <- attr(result, "quality_flags")
  expect_true(qf$cb_fallback)
})


test_that(".triple_goal stop_if_ties works", {
  # Create data likely to have tied mean ranks
  set.seed(42)
  s <- matrix(rnorm(50 * 3), nrow = 50, ncol = 3)

  # This should work with stop_if_ties = FALSE
  result <- DPMirt:::.triple_goal(s, stop_if_ties = FALSE)
  expect_s3_class(result, "data.frame")
})


test_that(".triple_goal GR preserves rank ordering", {
  # GR should roughly preserve the ordering from posterior means
  set.seed(42)
  s <- matrix(rnorm(2000 * 20, mean = rep(seq(-3, 3, length.out = 20), each = 2000),
                     sd = 0.5),
               nrow = 2000, ncol = 20)

  result <- DPMirt:::.triple_goal(s)

  # GR rank and PM rank should be highly correlated
  pm_rank <- rank(result$theta_pm)
  gr_rank <- rank(result$theta_gr)
  expect_true(cor(pm_rank, gr_rank) > 0.95)
})


test_that(".triple_goal CB expands distribution relative to PM", {
  # CB should have larger variance than PM
  set.seed(42)
  s <- matrix(rnorm(2000 * 30), nrow = 2000, ncol = 30)

  result <- DPMirt:::.triple_goal(s)

  expect_true(var(result$theta_cb) > var(result$theta_pm))
})


# ============================================================================
# dpmirt_estimates() — wrapper tests
# ============================================================================

test_that("dpmirt_estimates rejects non-fit input", {
  expect_error(
    dpmirt_estimates("not a fit"),
    "dpmirt_fit"
  )
})


test_that("dpmirt_estimates works with minimal fake fit (Rasch)", {
  set.seed(42)
  N <- 30; I <- 10; niter <- 500
  fake_fit <- structure(
    list(
      theta_samp = matrix(rnorm(niter * N), nrow = niter, ncol = N),
      beta_samp  = matrix(rnorm(niter * I), nrow = niter, ncol = I),
      lambda_samp = NULL,
      delta_samp  = NULL,
      config = list(model = "rasch", prior = "normal",
                    parameterization = "irt",
                    N = N, I = I)
    ),
    class = "dpmirt_fit"
  )

  est <- dpmirt_estimates(fake_fit)
  expect_s3_class(est, "dpmirt_estimates")
  expect_equal(nrow(est$theta), N)
  expect_equal(nrow(est$beta), I)
  expect_null(est$lambda)
  expect_null(est$delta)
})


test_that("dpmirt_estimates includes lambda for 2PL fit", {
  set.seed(42)
  N <- 30; I <- 10; niter <- 500
  fake_fit <- structure(
    list(
      theta_samp  = matrix(rnorm(niter * N), nrow = niter, ncol = N),
      beta_samp   = matrix(rnorm(niter * I), nrow = niter, ncol = I),
      lambda_samp = matrix(exp(rnorm(niter * I, 0.5, 0.3)),
                            nrow = niter, ncol = I),
      delta_samp  = NULL,
      config = list(model = "2pl", prior = "normal",
                    parameterization = "irt",
                    N = N, I = I)
    ),
    class = "dpmirt_fit"
  )

  est <- dpmirt_estimates(fake_fit)
  expect_s3_class(est, "dpmirt_estimates")
  expect_false(is.null(est$lambda))
  expect_equal(nrow(est$lambda), I)
  expect_true(all(c("lambda_pm", "lambda_psd", "lambda_lower",
                     "lambda_upper") %in% names(est$lambda)))
})


test_that("dpmirt_estimates credible intervals contain posterior mean", {
  set.seed(42)
  N <- 30; niter <- 2000
  fake_fit <- structure(
    list(
      theta_samp = matrix(rnorm(niter * N), nrow = niter, ncol = N),
      beta_samp  = matrix(rnorm(niter * 10), nrow = niter, ncol = 10),
      lambda_samp = NULL,
      delta_samp  = NULL,
      config = list(model = "rasch", prior = "normal",
                    parameterization = "irt",
                    N = N, I = 10)
    ),
    class = "dpmirt_fit"
  )

  est <- dpmirt_estimates(fake_fit, alpha = 0.05)

  # For most persons, PM should fall within 95% CI
  in_ci <- est$theta$theta_pm >= est$theta$theta_lower &
    est$theta$theta_pm <= est$theta$theta_upper
  expect_true(all(in_ci))
})


test_that("dpmirt_estimates methods='pm' skips CB/GR computation", {
  set.seed(42)
  N <- 30; niter <- 500
  fake_fit <- structure(
    list(
      theta_samp = matrix(rnorm(niter * N), nrow = niter, ncol = N),
      beta_samp  = matrix(rnorm(niter * 10), nrow = niter, ncol = 10),
      lambda_samp = NULL,
      delta_samp  = NULL,
      config = list(model = "rasch", prior = "normal",
                    parameterization = "irt",
                    N = N, I = 10)
    ),
    class = "dpmirt_fit"
  )

  est <- dpmirt_estimates(fake_fit, methods = "pm", include_items = FALSE)
  expect_equal(est$methods, "pm")
  # Beta should still have PM (but no CB/GR since include_items = FALSE)
  expect_true("beta_pm" %in% names(est$beta))
  expect_false("beta_cb" %in% names(est$beta))
})


# ============================================================================
# print.dpmirt_estimates tests
# ============================================================================

test_that("print.dpmirt_estimates shows lambda for 2PL", {
  set.seed(42)
  N <- 10; I <- 5; niter <- 100
  fake_fit <- structure(
    list(
      theta_samp  = matrix(rnorm(niter * N), nrow = niter, ncol = N),
      beta_samp   = matrix(rnorm(niter * I), nrow = niter, ncol = I),
      lambda_samp = matrix(exp(rnorm(niter * I)), nrow = niter, ncol = I),
      delta_samp  = NULL,
      config = list(model = "2pl", prior = "normal",
                    parameterization = "irt",
                    N = N, I = I)
    ),
    class = "dpmirt_fit"
  )

  est <- dpmirt_estimates(fake_fit)
  output <- capture.output(print(est))
  expect_true(any(grepl("Lambda", output, ignore.case = TRUE)))
})


# ============================================================================
# Loss function tests
# ============================================================================

test_that("dpmirt_loss computes MSEL correctly", {
  # Create a simple estimates object manually
  theta_df <- data.frame(
    theta_pm = c(0.0, 1.0, -1.0),
    theta_psd = c(0.5, 0.5, 0.5),
    theta_cb = c(0.0, 1.2, -1.2),
    theta_gr = c(-0.1, 0.9, -0.9),
    rbar = c(2, 3, 1),
    rhat = c(2L, 3L, 1L)
  )

  est <- structure(
    list(
      theta = theta_df,
      beta = NULL,
      lambda = NULL,
      delta = NULL,
      methods = c("pm", "cb", "gr"),
      alpha = 0.05,
      quality_flags = list()
    ),
    class = "dpmirt_estimates"
  )

  true_theta <- c(0.1, 0.8, -1.1)

  result <- dpmirt_loss(est, true_theta = true_theta, metrics = "msel")

  expect_true(all(c("parameter", "method", "msel") %in% names(result)))
  expect_equal(nrow(result), 3)  # pm, cb, gr
  expect_true(all(result$msel > 0))

  # Check PM MSEL by hand
  pm_msel <- mean((c(0.0, 1.0, -1.0) - true_theta)^2)
  expect_equal(result$msel[result$method == "pm"], pm_msel)
})


test_that("dpmirt_loss computes KS correctly", {
  theta_df <- data.frame(
    theta_pm = rnorm(50),
    theta_psd = rep(0.5, 50),
    theta_cb = rnorm(50),
    theta_gr = rnorm(50),
    rbar = 1:50,
    rhat = 1:50
  )

  est <- structure(
    list(theta = theta_df, beta = NULL, lambda = NULL,
         delta = NULL, methods = c("pm", "cb", "gr"),
         alpha = 0.05, quality_flags = list()),
    class = "dpmirt_estimates"
  )

  true_theta <- rnorm(50)
  result <- dpmirt_loss(est, true_theta = true_theta, metrics = "ks")

  expect_true(all(result$ks >= 0))
  expect_true(all(result$ks <= 1))
})


test_that("dpmirt_loss computes MSELR correctly", {
  theta_df <- data.frame(
    theta_pm = c(1, 2, 3, 4, 5),
    theta_psd = rep(0.5, 5),
    theta_cb = c(1, 2, 3, 4, 5),
    theta_gr = c(1, 2, 3, 4, 5),
    rbar = 1:5,
    rhat = 1:5
  )

  est <- structure(
    list(theta = theta_df, beta = NULL, lambda = NULL,
         delta = NULL, methods = c("pm"),
         alpha = 0.05, quality_flags = list()),
    class = "dpmirt_estimates"
  )

  # Perfect ranking → MSELR = 0
  true_theta <- c(1, 2, 3, 4, 5)
  result <- dpmirt_loss(est, true_theta = true_theta, metrics = "mselr")
  expect_equal(result$mselr[result$method == "pm"], 0)

  # Reversed ranking → MSELR > 0
  true_theta_rev <- c(5, 4, 3, 2, 1)
  result_rev <- dpmirt_loss(est, true_theta = true_theta_rev, metrics = "mselr")
  expect_true(result_rev$mselr[result_rev$method == "pm"] > 0)
})


test_that("dpmirt_loss computes item beta losses", {
  theta_df <- data.frame(
    theta_pm = rnorm(20),
    theta_psd = rep(0.5, 20),
    theta_cb = rnorm(20),
    theta_gr = rnorm(20),
    rbar = 1:20,
    rhat = 1:20
  )
  beta_df <- data.frame(
    beta_pm = c(-1, 0, 1),
    beta_psd = c(0.3, 0.3, 0.3),
    beta_cb = c(-1.1, 0.1, 1.1),
    beta_gr = c(-0.9, 0, 0.9),
    rbar = c(1, 2, 3),
    rhat = c(1L, 2L, 3L),
    beta_lower = c(-1.5, -0.5, 0.5),
    beta_upper = c(-0.5, 0.5, 1.5)
  )

  est <- structure(
    list(theta = theta_df, beta = beta_df, lambda = NULL,
         delta = NULL, methods = c("pm", "cb", "gr"),
         alpha = 0.05, quality_flags = list()),
    class = "dpmirt_estimates"
  )

  true_theta <- rnorm(20)
  true_beta <- c(-1.2, 0.1, 0.8)

  result <- dpmirt_loss(est, true_theta = true_theta,
                         true_beta = true_beta, metrics = "msel")

  # Should have theta rows + beta rows
  expect_true(any(result$parameter == "theta"))
  expect_true(any(result$parameter == "beta"))
})


test_that("dpmirt_loss computes lambda losses for 2PL", {
  theta_df <- data.frame(
    theta_pm = rnorm(20),
    theta_psd = rep(0.5, 20),
    theta_cb = rnorm(20),
    theta_gr = rnorm(20),
    rbar = 1:20,
    rhat = 1:20
  )

  lambda_df <- data.frame(
    lambda_pm    = c(1.2, 0.8, 1.5),
    lambda_psd   = c(0.1, 0.1, 0.1),
    lambda_lower = c(1.0, 0.6, 1.3),
    lambda_upper = c(1.4, 1.0, 1.7)
  )

  est <- structure(
    list(theta = theta_df, beta = NULL,
         lambda = lambda_df, delta = NULL,
         methods = c("pm"),
         alpha = 0.05, quality_flags = list()),
    class = "dpmirt_estimates"
  )

  true_theta <- rnorm(20)
  true_lambda <- c(1.0, 1.0, 1.5)

  result <- dpmirt_loss(est, true_theta = true_theta,
                         true_lambda = true_lambda, metrics = "msel")

  expect_true(any(result$parameter == "lambda"))
  lambda_row <- result[result$parameter == "lambda", ]
  expect_equal(nrow(lambda_row), 1)
  expect_equal(lambda_row$method, "pm")
  # Check MSEL by hand
  expected_msel <- mean((c(1.2, 0.8, 1.5) - true_lambda)^2)
  expect_equal(lambda_row$msel, expected_msel)
})


test_that("dpmirt_loss custom_loss works", {
  theta_df <- data.frame(
    theta_pm = c(1, 2, 3),
    theta_psd = rep(0.5, 3),
    rbar = 1:3,
    rhat = 1:3
  )

  est <- structure(
    list(theta = theta_df, beta = NULL, lambda = NULL,
         delta = NULL, methods = c("pm"),
         alpha = 0.05, quality_flags = list()),
    class = "dpmirt_estimates"
  )

  # Custom loss: mean absolute error
  mae <- function(est, true) mean(abs(est - true))
  true_theta <- c(1.1, 1.9, 3.2)

  result <- dpmirt_loss(est, true_theta = true_theta,
                         metrics = "msel", custom_loss = mae)

  expect_true("custom" %in% names(result))
  expected_mae <- mean(abs(c(1, 2, 3) - true_theta))
  expect_equal(result$custom[result$method == "pm"], expected_mae)
})


# ============================================================================
# dpmirt_draws() tests
# ============================================================================

test_that("dpmirt_draws rejects non-fit input", {
  expect_error(
    dpmirt_draws("not a fit"),
    "dpmirt_fit"
  )
})


test_that("dpmirt_draws returns matrix format", {
  set.seed(42)
  N <- 20; I <- 5; niter <- 100
  fake_fit <- structure(
    list(
      theta_samp = matrix(rnorm(niter * N), nrow = niter, ncol = N),
      beta_samp  = matrix(rnorm(niter * I), nrow = niter, ncol = I),
      lambda_samp = NULL,
      delta_samp  = NULL,
      config = list(model = "rasch", N = N, I = I)
    ),
    class = "dpmirt_fit"
  )

  draws <- dpmirt_draws(fake_fit, vars = "theta", format = "matrix")
  expect_true(is.matrix(draws))
  expect_equal(dim(draws), c(niter, N))

  draws_beta <- dpmirt_draws(fake_fit, vars = "beta", format = "matrix")
  expect_equal(dim(draws_beta), c(niter, I))
})


test_that("dpmirt_draws returns long format", {
  set.seed(42)
  N <- 20; I <- 5; niter <- 100
  fake_fit <- structure(
    list(
      theta_samp = matrix(rnorm(niter * N), nrow = niter, ncol = N),
      beta_samp  = matrix(rnorm(niter * I), nrow = niter, ncol = I),
      lambda_samp = NULL,
      delta_samp  = NULL,
      config = list(model = "rasch", N = N, I = I)
    ),
    class = "dpmirt_fit"
  )

  draws <- dpmirt_draws(fake_fit, vars = "theta", format = "long")
  expect_s3_class(draws, "data.frame")
  expect_true(all(c("iteration", "index", "value", "variable") %in% names(draws)))
  expect_equal(nrow(draws), niter * N)
})


test_that("dpmirt_draws returns multiple vars as list in matrix format", {
  set.seed(42)
  N <- 20; I <- 5; niter <- 100
  fake_fit <- structure(
    list(
      theta_samp = matrix(rnorm(niter * N), nrow = niter, ncol = N),
      beta_samp  = matrix(rnorm(niter * I), nrow = niter, ncol = I),
      lambda_samp = matrix(exp(rnorm(niter * I)), nrow = niter, ncol = I),
      delta_samp  = NULL,
      config = list(model = "2pl", N = N, I = I)
    ),
    class = "dpmirt_fit"
  )

  draws <- dpmirt_draws(fake_fit, vars = c("theta", "lambda"), format = "matrix")
  expect_true(is.list(draws))
  expect_true("theta" %in% names(draws))
  expect_true("lambda" %in% names(draws))
  expect_equal(dim(draws$theta), c(niter, N))
  expect_equal(dim(draws$lambda), c(niter, I))
})


test_that("dpmirt_draws errors for missing variable", {
  set.seed(42)
  fake_fit <- structure(
    list(
      theta_samp = matrix(rnorm(100), nrow = 50, ncol = 2),
      beta_samp  = NULL,
      lambda_samp = NULL,
      delta_samp  = NULL,
      config = list(model = "rasch", N = 2, I = 0)
    ),
    class = "dpmirt_fit"
  )

  expect_error(
    dpmirt_draws(fake_fit, vars = "beta"),
    "No samples found"
  )
})


test_that("dpmirt_draws long format stacks multiple vars", {
  set.seed(42)
  N <- 10; I <- 5; niter <- 50
  fake_fit <- structure(
    list(
      theta_samp = matrix(rnorm(niter * N), nrow = niter, ncol = N),
      beta_samp  = matrix(rnorm(niter * I), nrow = niter, ncol = I),
      lambda_samp = NULL,
      delta_samp  = NULL,
      config = list(model = "rasch", N = N, I = I)
    ),
    class = "dpmirt_fit"
  )

  draws <- dpmirt_draws(fake_fit, vars = c("theta", "beta"), format = "long")
  expect_s3_class(draws, "data.frame")
  expect_true("theta" %in% draws$variable)
  expect_true("beta" %in% draws$variable)
  expect_equal(nrow(draws), niter * N + niter * I)
})


# ============================================================================
# KS statistic internal tests
# ============================================================================

test_that(".ks_stat returns 0 for identical samples", {
  x <- c(1, 2, 3, 4, 5)
  expect_equal(DPMirt:::.ks_stat(x, x), 0)
})


test_that(".ks_stat returns 1 for non-overlapping samples", {
  x <- c(1, 2, 3)
  y <- c(10, 11, 12)
  ks <- DPMirt:::.ks_stat(x, y)
  expect_equal(ks, 1)
})
