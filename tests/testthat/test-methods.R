# ============================================================================
# Tests for R/methods.R
# S3 methods: print, summary, coef for dpmirt_fit
# ============================================================================

# ============================================================================
# Helper: create a mock dpmirt_fit object
# ============================================================================

.mock_fit_full <- function(model = "rasch", prior = "normal",
                           parameterization = "irt",
                           N = 30, I = 10, niter = 200, seed = 42) {
  set.seed(seed)
  structure(
    list(
      theta_samp = matrix(rnorm(niter * N), nrow = niter, ncol = N),
      beta_samp  = matrix(rnorm(niter * I, sd = 1.5), nrow = niter, ncol = I),
      lambda_samp = if (model %in% c("2pl", "3pl"))
        matrix(exp(rnorm(niter * I, 0, 0.3)), nrow = niter, ncol = I) else NULL,
      delta_samp = if (model == "3pl")
        matrix(rbeta(niter * I, 2, 8), nrow = niter, ncol = I) else NULL,
      config = list(
        model = model, prior = prior,
        parameterization = parameterization,
        identification = "unconstrained",
        N = N, I = I,
        niter = niter, nburnin = 50, thin = 1, nchains = 1,
        rescale = TRUE
      ),
      ess = list(
        items = runif(I, 100, 400),
        theta = runif(N, 80, 350)
      ),
      waic = 2500.3,
      cluster_info = if (prior == "dpm")
        list(
          n_clusters = sample(3:8, niter, replace = TRUE),
          alpha_summary = list(mean = 0.45, median = 0.4, sd = 0.15,
                               q025 = 0.2, q975 = 0.85)
        ) else NULL,
      dp_density = if (prior == "dpm")
        list(grid = seq(-3, 3, length.out = 100),
             density_mean = dnorm(seq(-3, 3, length.out = 100))) else NULL,
      compilation_time = 35.2,
      sampling_time = 8.5,
      total_time = 43.7,
      other_samp = NULL
    ),
    class = "dpmirt_fit"
  )
}


# ============================================================================
# print.dpmirt_fit tests
# ============================================================================

test_that("print.dpmirt_fit produces expected output for Rasch-Normal", {
  fit <- .mock_fit_full()
  output <- capture.output(print(fit))

  expect_true(any(grepl("DPMirt Model Fit", output)))
  expect_true(any(grepl("RASCH", output)))
  expect_true(any(grepl("normal", output)))
  expect_true(any(grepl("WAIC", output)))
  expect_true(any(grepl("Total time", output)))
})


test_that("print.dpmirt_fit shows ESS info", {
  fit <- .mock_fit_full()
  output <- capture.output(print(fit))

  expect_true(any(grepl("Min ESS", output)))
})


test_that("print.dpmirt_fit shows cluster info for DPM", {
  fit <- .mock_fit_full(prior = "dpm")
  output <- capture.output(print(fit))

  expect_true(any(grepl("Cluster", output, ignore.case = TRUE)))
  expect_true(any(grepl("Alpha", output, ignore.case = TRUE)))
})


test_that("print.dpmirt_fit returns invisible(x)", {
  fit <- .mock_fit_full()
  result <- capture.output(ret <- print(fit))
  expect_identical(ret, fit)
})


# ============================================================================
# summary.dpmirt_fit tests
# ============================================================================

test_that("summary.dpmirt_fit produces comprehensive output for Rasch", {
  fit <- .mock_fit_full()
  output <- capture.output(summary(fit))

  expect_true(any(grepl("DPMirt Model Summary", output)))
  expect_true(any(grepl("Model Configuration", output)))
  expect_true(any(grepl("Data", output)))
  expect_true(any(grepl("MCMC Settings", output)))
  expect_true(any(grepl("Timing", output)))
  expect_true(any(grepl("Item Difficulty", output)))
  expect_true(any(grepl("Person Ability", output)))
  expect_true(any(grepl("WAIC", output)))
})


test_that("summary.dpmirt_fit shows lambda for 2PL", {
  fit <- .mock_fit_full(model = "2pl")
  output <- capture.output(summary(fit))

  expect_true(any(grepl("Lambda", output)))
})


test_that("summary.dpmirt_fit shows delta for 3PL", {
  fit <- .mock_fit_full(model = "3pl")
  output <- capture.output(summary(fit))

  expect_true(any(grepl("Delta", output)))
})


test_that("summary.dpmirt_fit shows SI parameterization correctly", {
  fit <- .mock_fit_full(model = "2pl", parameterization = "si")
  output <- capture.output(summary(fit))

  expect_true(any(grepl("Intercept", output)) ||
              any(grepl("gamma", output)))
})


test_that("summary.dpmirt_fit shows DPM diagnostics", {
  fit <- .mock_fit_full(prior = "dpm")
  output <- capture.output(summary(fit))

  expect_true(any(grepl("DPM Diagnostics", output)))
  expect_true(any(grepl("Alpha", output, ignore.case = TRUE)))
  expect_true(any(grepl("cluster", output, ignore.case = TRUE)))
})


test_that("summary.dpmirt_fit shows DP density info", {
  fit <- .mock_fit_full(prior = "dpm")
  output <- capture.output(summary(fit))

  expect_true(any(grepl("DP density", output, ignore.case = TRUE)))
  expect_true(any(grepl("100", output)))  # 100 grid points
})


test_that("summary.dpmirt_fit returns invisible(object)", {
  fit <- .mock_fit_full()
  result <- capture.output(ret <- summary(fit))
  expect_identical(ret, fit)
})


# ============================================================================
# coef.dpmirt_fit tests
# ============================================================================

test_that("coef.dpmirt_fit returns item estimates by default", {
  fit <- .mock_fit_full()
  result <- coef(fit)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 10)  # I = 10
  expect_true("beta" %in% names(result))
  expect_true(all(grepl("^item_", rownames(result))))
})


test_that("coef.dpmirt_fit returns person estimates", {
  fit <- .mock_fit_full()
  result <- coef(fit, type = "persons")

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 30)  # N = 30
  expect_true("theta" %in% names(result))
  expect_true(all(grepl("^person_", rownames(result))))
})


test_that("coef.dpmirt_fit includes lambda for 2PL items", {
  fit <- .mock_fit_full(model = "2pl")
  result <- coef(fit, type = "items")

  expect_true("beta" %in% names(result))
  expect_true("lambda" %in% names(result))
})


test_that("coef.dpmirt_fit includes delta for 3PL items", {
  fit <- .mock_fit_full(model = "3pl")
  result <- coef(fit, type = "items")

  expect_true("delta" %in% names(result))
})


test_that("coef.dpmirt_fit uses gamma for SI parameterization", {
  fit <- .mock_fit_full(model = "2pl", parameterization = "si")
  result <- coef(fit, type = "items")

  expect_true("gamma" %in% names(result))
  expect_false("beta" %in% names(result))
})


test_that("coef.dpmirt_fit values match colMeans", {
  fit <- .mock_fit_full()
  result <- coef(fit, type = "items")

  expected_beta <- colMeans(fit$beta_samp)
  expect_equal(result$beta, expected_beta)

  result_person <- coef(fit, type = "persons")
  expected_theta <- colMeans(fit$theta_samp)
  expect_equal(result_person$theta, expected_theta)
})
