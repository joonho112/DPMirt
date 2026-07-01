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


test_that("dpmirt_diagnostics ignores non-finite values in min ESS", {
  fit <- .mock_fit()
  fit$ess$items <- c(NA_real_, 12, Inf, 8)
  fit$ess$theta <- c(Inf, 5, NA_real_)

  diag <- dpmirt_diagnostics(fit)

  expect_equal(diag$ess_min_items, 8)
  expect_equal(diag$ess_min_theta, 5)
})


test_that("dpmirt_diagnostics handles missing ESS gracefully", {
  fit <- .mock_fit()
  fit$ess$items <- NULL
  fit$ess$theta <- NULL

  diag <- dpmirt_diagnostics(fit)
  expect_true(is.na(diag$ess_min_items))
  expect_true(is.na(diag$ess_min_theta))
})


test_that(".build_chain_diagnostics computes R-hat and by-chain ESS", {
  set.seed(42)
  n_per_chain <- 30
  n_total <- 2 * n_per_chain

  beta_samp <- rbind(
    matrix(rnorm(n_per_chain * 2, mean = 0), nrow = n_per_chain),
    matrix(rnorm(n_per_chain * 2, mean = 0.1), nrow = n_per_chain)
  )
  colnames(beta_samp) <- paste0("beta[", 1:2, "]")

  theta_samp <- rbind(
    matrix(rnorm(n_per_chain * 3, mean = 0), nrow = n_per_chain),
    matrix(rnorm(n_per_chain * 3, mean = 0.1), nrow = n_per_chain)
  )
  colnames(theta_samp) <- paste0("eta[", 1:3, "]")

  make_index <- function() {
    data.frame(
      row = seq_len(n_total),
      chain_id = rep(1:2, each = n_per_chain),
      within_chain_draw = rep(seq_len(n_per_chain), times = 2),
      mcmc_iteration = rep(seq_len(n_per_chain), times = 2)
    )
  }

  draw_index <- list(main = make_index(), theta = make_index())
  chain_info <- data.frame(
    chain_id = 1:2,
    waic = c(100, 102)
  )

  diag <- DPMirt:::.build_chain_diagnostics(
    ess = list(items = rep(NA_real_, 2), theta = rep(NA_real_, 3)),
    waic = 101,
    loglik_trace = rnorm(n_total),
    cluster_info = NULL,
    chain_info = chain_info,
    draw_index = draw_index,
    draws = list(items = beta_samp, theta = theta_samp)
  )

  expected_items <- coda::gelman.diag(
    coda::mcmc.list(
      coda::as.mcmc(beta_samp[1:n_per_chain, , drop = FALSE]),
      coda::as.mcmc(beta_samp[(n_per_chain + 1):n_total, , drop = FALSE])
    ),
    autoburnin = FALSE,
    multivariate = FALSE
  )$psrf[, "Point est."]

  expect_equal(unname(diag$rhat$items), unname(expected_items),
               tolerance = 1e-10)
  expect_true(is.data.frame(diag$ess$by_chain$items))
  expect_equal(unique(diag$ess$by_chain$items$chain_id), 1:2)
  expect_equal(nrow(diag$ess$by_chain$items), 2L * ncol(beta_samp))
  expect_setequal(diag$ess$by_chain$items$parameter, colnames(beta_samp))
  expect_true(all(is.finite(diag$ess$by_chain$items$ess)))
  expect_equal(diag$waic$value, 101)
  expect_equal(diag$waic$by_chain$waic, c(100, 102))
  expect_equal(diag$waic$aggregation, "mean_of_chain_waic")
  expect_equal(nrow(diag$loglik$by_chain), n_total)
})


test_that(".build_chain_diagnostics omits R-hat for single-chain fits", {
  set.seed(42)
  niter <- 20
  beta_samp <- matrix(rnorm(niter * 2), nrow = niter)
  colnames(beta_samp) <- paste0("beta[", 1:2, "]")
  draw_index <- list(
    main = data.frame(
      row = seq_len(niter),
      chain_id = rep(1L, niter),
      within_chain_draw = seq_len(niter),
      mcmc_iteration = seq_len(niter)
    ),
    theta = NULL
  )

  diag <- DPMirt:::.build_chain_diagnostics(
    ess = list(items = rep(NA_real_, 2)),
    waic = NULL,
    loglik_trace = NULL,
    cluster_info = NULL,
    chain_info = data.frame(chain_id = 1L, waic = NA_real_),
    draw_index = draw_index,
    draws = list(items = beta_samp)
  )

  expect_null(diag$rhat)
  expect_true(is.data.frame(diag$ess$by_chain$items))
})


test_that("dpmirt_diagnostics exposes chain-aware diagnostics when present", {
  fit <- .mock_fit(N = 3, I = 2, niter = 20)
  fit$diagnostics <- list(
    ess = list(
      pooled = fit$ess,
      by_chain = list(
        items = data.frame(
          chain_id = c(1L, 2L),
          parameter = c("beta[1]", "beta[1]"),
          ess = c(10, 12)
        )
      )
    ),
    rhat = list(items = c("beta[1]" = 1.03, "beta[2]" = 1.01)),
    waic = list(
      value = fit$waic,
      by_chain = data.frame(chain_id = 1:2, waic = c(100, 102)),
      aggregation = "mean_of_chain_waic"
    ),
    loglik = list(
      trace = fit$loglik_trace,
      by_chain = data.frame(
        chain_id = rep(1:2, each = 10),
        row = 1:20,
        within_chain_draw = rep(1:10, times = 2),
        mcmc_iteration = rep(1:10, times = 2),
        value = fit$loglik_trace
      )
    ),
    clusters = list(trace = NULL, by_chain = NULL)
  )
  fit$chain_info <- data.frame(chain_id = 1:2)

  diag <- dpmirt_diagnostics(fit)

  expect_equal(diag$rhat$items, fit$diagnostics$rhat$items)
  expect_equal(diag$rhat_max, 1.03)
  expect_equal(diag$ess_by_chain$items$ess, c(10, 12))
  expect_equal(diag$waic_by_chain$chain_id, 1:2)
  expect_equal(nrow(diag$loglik_by_chain), 20)
})


test_that("dpmirt_diagnostics and print handle unavailable diagnostics cleanly", {
  fit <- .mock_fit()
  fit$ess$items <- c(NA_real_, Inf)
  fit$ess$theta <- numeric(0)
  fit$waic <- NA_real_
  fit$compilation_time <- NULL
  fit$sampling_time <- NA_real_
  fit$total_time <- NULL

  diag <- dpmirt_diagnostics(fit)

  expect_true(is.na(diag$ess_min_items))
  expect_true(is.na(diag$ess_min_theta))

  output <- capture.output(print(diag))
  expect_false(any(grepl("^WAIC:", output)))
  expect_true(any(grepl("not available", output)))
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


test_that("dpmirt_compare requires at least two fits", {
  fit <- .mock_fit()
  expect_error(
    dpmirt_compare(fit),
    "At least two"
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
  fit2$waic <- NA_real_

  expect_warning(
    result <- dpmirt_compare(fit1, fit2),
    "unavailable WAIC"
  )
  expect_true(any(is.na(result$waic)))
  expect_true(any(is.na(result$delta_waic)))
})


test_that("dpmirt_compare errors clearly when all WAIC values are missing", {
  fit1 <- .mock_fit()
  fit2 <- .mock_fit()
  fit1$waic <- NULL
  fit2$waic <- NA_real_

  expect_error(
    dpmirt_compare(fit1, fit2),
    "No comparable WAIC"
  )
})


test_that("dpmirt_compare warns for mean-of-run WAIC provenance", {
  fit1 <- .mock_fit()
  fit2 <- .mock_fit()
  fit1$diagnostics <- list(
    waic = list(aggregation = "mean_of_chain_waic")
  )

  expect_warning(
    dpmirt_compare(fit1, fit2),
    "mean-of-run WAIC"
  )

  fit1$diagnostics <- NULL
  fit1$chain_info <- data.frame(chain_id = 1:2)
  expect_warning(
    dpmirt_compare(fit1, fit2),
    "mean-of-run WAIC"
  )
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


test_that(".summarize_n_clusters handles length-one traces", {
  result <- DPMirt:::.summarize_n_clusters(4L)

  expect_equal(result$mean, 4)
  expect_equal(result$median, 4)
  expect_true(is.na(result$sd))
  expect_equal(result$min, 4)
  expect_equal(result$max, 4)
  expect_equal(result$mode, 4)
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


test_that(".extract_cluster_info counts unique labels, not label magnitudes", {
  samples <- matrix(
    c(7, 7, 7,
      10, 20, 20,
      30, 30, 30),
    nrow = 3,
    byrow = TRUE
  )
  colnames(samples) <- paste0("zi[", 1:3, "]")

  result <- DPMirt:::.extract_cluster_info(samples, N = 3)

  expect_equal(result$n_clusters, c(1L, 2L, 1L))
  expect_null(result$alpha_summary)
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
