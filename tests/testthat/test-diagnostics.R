# ============================================================================
# Tests for R/diagnostics.R
# dpmirt_diagnostics(), dpmirt_compare(), and internal helpers
# ============================================================================

# ============================================================================
# Helper: create a mock dpmirt_fit object
# ============================================================================

.mock_fit <- function(model = "rasch", prior = "normal",
                      N = 30, I = 10, niter = 500, seed = 42) {
  set.seed(seed)
  structure(
    list(
      theta_samp = matrix(rnorm(niter * N), nrow = niter, ncol = N),
      beta_samp  = matrix(rnorm(niter * I), nrow = niter, ncol = I),
      lambda_samp = if (model %in% c("2pl", "3pl"))
        matrix(exp(rnorm(niter * I)), nrow = niter, ncol = I) else NULL,
      delta_samp  = NULL,
      config = list(
        model = model, prior = prior,
        parameterization = "irt",
        identification = "unconstrained",
        N = N, I = I,
        niter = niter, nburnin = 100, thin = 1, nchains = 1,
        rescale = TRUE
      ),
      ess = list(
        items = runif(I, 200, 500),
        theta = runif(N, 150, 400)
      ),
      waic = runif(1, 2000, 5000),
      loglik_trace = rnorm(niter, -3000, 100),
      cluster_info = if (prior == "dpm")
        list(
          n_clusters = sample(3:10, niter, replace = TRUE),
          alpha_summary = list(
            mean = 0.5, median = 0.45, sd = 0.2,
            q025 = 0.15, q975 = 0.95
          )
        ) else NULL,
      compilation_time = 45.2,
      sampling_time = 12.3,
      total_time = 57.5
    ),
    class = "dpmirt_fit"
  )
}


# ============================================================================
# dpmirt_diagnostics() tests
# ============================================================================

test_that("dpmirt_diagnostics rejects non-fit input", {
  expect_error(
    dpmirt_diagnostics("not a fit"),
    "dpmirt_fit"
  )
})


test_that("dpmirt_diagnostics returns correct structure for Normal prior", {
  fit <- .mock_fit(prior = "normal")
  diag <- dpmirt_diagnostics(fit)

  expect_s3_class(diag, "dpmirt_diagnostics")
  expect_true(!is.null(diag$ess))
  expect_true(!is.null(diag$waic))
  expect_true(!is.null(diag$loglik_trace))
  expect_true(is.numeric(diag$ess_min_items))
  expect_true(is.numeric(diag$ess_min_theta))
  expect_true(is.numeric(diag$compilation_time))
  expect_true(is.numeric(diag$sampling_time))
  expect_true(is.numeric(diag$total_time))
})


test_that("dpmirt_diagnostics includes DPM cluster info for DPM prior", {
  fit <- .mock_fit(prior = "dpm")
  diag <- dpmirt_diagnostics(fit)

  expect_true(!is.null(diag$n_clusters))
  expect_true(!is.null(diag$n_clusters_summary))
  expect_true(!is.null(diag$alpha_summary))
  expect_equal(diag$alpha_summary$mean, 0.5)
})


test_that("dpmirt_diagnostics computes min ESS correctly", {
  fit <- .mock_fit()
  diag <- dpmirt_diagnostics(fit)

  expect_equal(diag$ess_min_items, min(fit$ess$items))
  expect_equal(diag$ess_min_theta, min(fit$ess$theta))
})


test_that("dpmirt_diagnostics handles missing ESS gracefully", {
  fit <- .mock_fit()
  fit$ess$items <- NULL
  fit$ess$theta <- NULL

  diag <- dpmirt_diagnostics(fit)
  expect_true(is.na(diag$ess_min_items))
  expect_true(is.na(diag$ess_min_theta))
})


# ============================================================================
# print.dpmirt_diagnostics tests
# ============================================================================

test_that("print.dpmirt_diagnostics produces output for Normal prior", {
  fit <- .mock_fit(prior = "normal")
  diag <- dpmirt_diagnostics(fit)

  output <- capture.output(print(diag))
  expect_true(any(grepl("DPMirt MCMC Diagnostics", output)))
  expect_true(any(grepl("ESS", output)))
  expect_true(any(grepl("WAIC", output)))
  expect_true(any(grepl("Timing", output)))
})


test_that("print.dpmirt_diagnostics shows cluster info for DPM", {
  fit <- .mock_fit(prior = "dpm")
  diag <- dpmirt_diagnostics(fit)

  output <- capture.output(print(diag))
  expect_true(any(grepl("Cluster", output, ignore.case = TRUE)))
  expect_true(any(grepl("Alpha", output, ignore.case = TRUE)))
})


# ============================================================================
# dpmirt_compare() tests
# ============================================================================

test_that("dpmirt_compare rejects non-fit input", {
  fit <- .mock_fit()
  expect_error(
    dpmirt_compare(fit, "not a fit"),
    "not a dpmirt_fit"
  )
})


test_that("dpmirt_compare returns ordered data.frame", {
  fit1 <- .mock_fit(prior = "normal")
  fit1$waic <- 3000
  fit2 <- .mock_fit(prior = "dpm")
  fit2$waic <- 2800

  result <- dpmirt_compare(fit1, fit2)

  expect_s3_class(result, "data.frame")
  expect_true(all(c("model", "waic", "delta_waic") %in% names(result)))
  expect_equal(nrow(result), 2)

  # Result should be ordered by WAIC (ascending)
  expect_true(result$waic[1] <= result$waic[2])

  # Best model delta_waic should be 0
  expect_equal(result$delta_waic[1], 0)
})


test_that("dpmirt_compare generates correct labels", {
  fit1 <- .mock_fit(model = "rasch", prior = "normal")
  fit2 <- .mock_fit(model = "rasch", prior = "dpm")

  result <- dpmirt_compare(fit1, fit2)

  expect_true("RASCH-normal" %in% result$model)
  expect_true("RASCH-dpm" %in% result$model)
})


test_that("dpmirt_compare handles missing WAIC", {
  fit1 <- .mock_fit()
  fit2 <- .mock_fit()
  fit2$waic <- NULL

  result <- dpmirt_compare(fit1, fit2)
  expect_true(is.na(result$waic[result$model == "RASCH-normal"][1]) ||
              is.na(result$waic[2]) ||
              any(is.na(result$waic)))
})


test_that("dpmirt_compare compares three or more models", {
  fit1 <- .mock_fit(prior = "normal")
  fit1$waic <- 3000
  fit2 <- .mock_fit(prior = "dpm")
  fit2$waic <- 2800
  fit3 <- .mock_fit(model = "2pl", prior = "normal")
  fit3$waic <- 2500

  result <- dpmirt_compare(fit1, fit2, fit3)
  expect_equal(nrow(result), 3)
  expect_true(result$waic[1] <= result$waic[2])
  expect_true(result$waic[2] <= result$waic[3])
})


test_that("dpmirt_compare rejects unsupported criterion", {
  fit1 <- .mock_fit()
  fit2 <- .mock_fit()

  expect_error(
    dpmirt_compare(fit1, fit2, criterion = "loo"),
    "waic"
  )
})


# ============================================================================
# Internal helpers
# ============================================================================

test_that(".summarize_n_clusters returns correct structure", {
  n_clusters <- sample(3:8, 100, replace = TRUE)
  result <- DPMirt:::.summarize_n_clusters(n_clusters)

  expect_true(is.list(result))
  expect_true(all(c("mean", "median", "sd", "q025", "q975",
                     "min", "max", "mode") %in% names(result)))
  expect_equal(result$mean, mean(n_clusters))
  expect_equal(result$median, median(n_clusters))
  expect_equal(result$min, min(n_clusters))
  expect_equal(result$max, max(n_clusters))
})


test_that(".summarize_n_clusters returns NULL for NULL input", {
  result <- DPMirt:::.summarize_n_clusters(NULL)
  expect_null(result)
})


test_that(".numeric_mode returns correct mode", {
  x <- c(3, 3, 3, 5, 5, 7)
  expect_equal(DPMirt:::.numeric_mode(x), 3)

  y <- c(1, 2, 2, 2, 3, 3)
  expect_equal(DPMirt:::.numeric_mode(y), 2)
})


test_that(".extract_cluster_info extracts zi counts", {
  set.seed(42)
  niter <- 50; N <- 10
  zi_samp <- matrix(sample(1:4, niter * N, replace = TRUE),
                     nrow = niter, ncol = N)
  colnames(zi_samp) <- paste0("zi[", seq_len(N), "]")

  alpha_samp <- matrix(rexp(niter, 2), nrow = niter, ncol = 1)
  colnames(alpha_samp) <- "alpha"

  samples <- cbind(zi_samp, alpha_samp)
  result <- DPMirt:::.extract_cluster_info(samples, N)

  expect_equal(length(result$n_clusters), niter)
  expect_true(all(result$n_clusters >= 1))
  expect_true(all(result$n_clusters <= 4))
  expect_true(!is.null(result$alpha_summary))
  expect_equal(result$alpha_summary$mean, mean(alpha_samp))
})


test_that(".extract_cluster_info returns NULL when no zi columns", {
  set.seed(42)
  niter <- 50
  samples <- matrix(rnorm(niter * 3), nrow = niter, ncol = 3)
  colnames(samples) <- c("beta[1]", "beta[2]", "beta[3]")

  result <- DPMirt:::.extract_cluster_info(samples, 10)
  expect_null(result$n_clusters)
  expect_null(result$alpha_summary)
})
