# ============================================================================
# Tests: Visualization (Phase 7B)
# ============================================================================

# --------------------------------------------------------------------------
# Mock object factory
# --------------------------------------------------------------------------

.mock_fit <- function(model = "rasch", prior = "normal",
                      N = 50, I = 10, niter = 100) {
  set.seed(42)
  fit <- list(
    theta_samp  = matrix(rnorm(niter * N), nrow = niter, ncol = N),
    beta_samp   = matrix(rnorm(niter * I, sd = 0.5), nrow = niter, ncol = I),
    lambda_samp = NULL,
    delta_samp  = NULL,
    loglik_trace = cumsum(rnorm(niter, -500, 10)),
    ess = list(items = runif(I, 100, 500), theta = runif(N, 80, 400)),
    dp_density = NULL,
    cluster_info = NULL,
    config = list(model = model, prior = prior, N = N, I = I, niter = niter,
                  parameterization = "irt", identification = "unconstrained",
                  nburnin = 20, thin = 1, rescale = TRUE)
  )

  if (model %in% c("2pl", "3pl")) {
    fit$lambda_samp <- matrix(rlnorm(niter * I, 0, 0.3),
                              nrow = niter, ncol = I)
  }
  if (model == "3pl") {
    fit$delta_samp <- matrix(rbeta(niter * I, 2, 8),
                             nrow = niter, ncol = I)
  }
  if (prior == "dpm") {
    fit$cluster_info <- list(
      n_clusters = sample(3:8, niter, replace = TRUE),
      alpha_summary = list(mean = 1.5, median = 1.3, sd = 0.5)
    )
    grid <- seq(-4, 4, length.out = 100)
    fit$dp_density <- list(
      grid = grid,
      density_mean  = dnorm(grid),
      density_lower = dnorm(grid) * 0.8,
      density_upper = dnorm(grid) * 1.2,
      ci_level = 0.95
    )
    class(fit$dp_density) <- "dpmirt_dp_density"
  }

  class(fit) <- "dpmirt_fit"
  fit
}


.mock_estimates <- function() {
  N <- 50
  set.seed(42)
  theta_pm <- rnorm(N)
  theta_psd <- runif(N, 0.3, 0.6)
  est <- structure(
    list(
      theta = data.frame(
        theta_pm    = theta_pm,
        theta_cb    = theta_pm * 1.1,
        theta_gr    = theta_pm * 0.95,
        theta_psd   = theta_psd,
        theta_lower = theta_pm - 1.96 * theta_psd,
        theta_upper = theta_pm + 1.96 * theta_psd
      ),
      beta = data.frame(
        beta_pm    = rnorm(10, sd = 0.5),
        beta_psd   = runif(10, 0.1, 0.3),
        beta_lower = rnorm(10, sd = 0.5) - 0.5,
        beta_upper = rnorm(10, sd = 0.5) + 0.5
      ),
      methods = c("pm", "cb", "gr"),
      alpha = 0.95,
      quality_flags = list()
    ),
    class = "dpmirt_estimates"
  )
  est
}


.mock_sim <- function(model = "rasch") {
  set.seed(42)
  N <- 50
  I <- 10
  theta <- rnorm(N)
  beta <- rnorm(I, sd = 0.5)
  lambda <- if (model %in% c("2pl", "3pl")) rlnorm(I, 0, 0.3) else rep(1, I)
  delta <- if (model == "3pl") rbeta(I, 2, 8) else rep(0, I)

  P <- outer(theta, seq_len(I), function(th, j) {
    delta[j] + (1 - delta[j]) / (1 + exp(-lambda[j] * (th - beta[j])))
  })
  response <- matrix(rbinom(N * I, 1, P), nrow = N, ncol = I)

  structure(
    list(
      theta = theta, beta = beta, lambda = lambda, delta = delta,
      response = response, n_persons = N, n_items = I,
      model = model, reliability = 0.85, latent_shape = "normal"
    ),
    class = "dpmirt_sim"
  )
}


# ============================================================================
# Tests: Base R plot types (12 types)
# ============================================================================

test_that("base R: density plot works for Rasch/Normal", {
  fit <- .mock_fit("rasch", "normal")
  expect_silent(plot(fit, type = "density", engine = "base"))
})

test_that("base R: items plot works", {
  fit <- .mock_fit("rasch", "normal")
  expect_silent(plot(fit, type = "items", engine = "base"))
})

test_that("base R: trace plot works", {
  fit <- .mock_fit("rasch", "normal")
  expect_silent(plot(fit, type = "trace", engine = "base"))
})

test_that("base R: clusters plot works for DPM", {
  fit <- .mock_fit("rasch", "dpm")
  expect_silent(plot(fit, type = "clusters", engine = "base"))
})

test_that("base R: clusters plot errors for non-DPM", {
  fit <- .mock_fit("rasch", "normal")
  expect_error(plot(fit, type = "clusters", engine = "base"),
               "DPM models")
})

test_that("base R: dp_density plot works", {
  fit <- .mock_fit("rasch", "dpm")
  expect_silent(plot(fit, type = "dp_density", engine = "base"))
})

test_that("base R: dp_density plot errors without dp_density", {
  fit <- .mock_fit("rasch", "normal")
  expect_error(plot(fit, type = "dp_density", engine = "base"),
               "No DP density")
})

test_that("base R: ICC plot works for Rasch", {
  fit <- .mock_fit("rasch", "normal")
  expect_silent(plot(fit, type = "icc", engine = "base"))
})

test_that("base R: ICC plot works for 2PL", {
  fit <- .mock_fit("2pl", "normal")
  expect_silent(plot(fit, type = "icc", engine = "base"))
})

test_that("base R: ICC plot works for 3PL", {
  fit <- .mock_fit("3pl", "normal")
  expect_silent(plot(fit, type = "icc", engine = "base"))
})

test_that("base R: ICC plot with item selection", {
  fit <- .mock_fit("rasch", "normal")
  expect_silent(plot(fit, type = "icc", engine = "base", items = c(1, 3, 5)))
})

test_that("base R: wright_map works", {
  fit <- .mock_fit("rasch", "normal")
  expect_silent(plot(fit, type = "wright_map", engine = "base"))
})

test_that("base R: parameter_trace works for beta", {
  fit <- .mock_fit("rasch", "normal")
  expect_silent(plot(fit, type = "parameter_trace", engine = "base",
                     param = "beta", indices = 1:4))
})

test_that("base R: parameter_trace works for theta", {
  fit <- .mock_fit("rasch", "normal")
  expect_silent(plot(fit, type = "parameter_trace", engine = "base",
                     param = "theta", indices = 1:2))
})

test_that("base R: parameter_trace errors for unavailable param", {
  fit <- .mock_fit("rasch", "normal")
  expect_error(plot(fit, type = "parameter_trace", engine = "base",
                    param = "lambda"),
               "No lambda samples")
})

test_that("base R: caterpillar works for beta", {
  fit <- .mock_fit("rasch", "normal")
  expect_silent(plot(fit, type = "caterpillar", engine = "base",
                     param = "beta"))
})

test_that("base R: caterpillar works for theta with max_show", {
  fit <- .mock_fit("rasch", "normal")
  expect_silent(plot(fit, type = "caterpillar", engine = "base",
                     param = "theta", max_show = 20))
})

test_that("base R: density_compare works (Normal prior)", {
  fit <- .mock_fit("rasch", "normal")
  expect_silent(plot(fit, type = "density_compare", engine = "base"))
})

test_that("base R: density_compare works (DPM prior with dp_density)", {
  fit <- .mock_fit("rasch", "dpm")
  expect_silent(plot(fit, type = "density_compare", engine = "base"))
})

test_that("base R: info plot works for Rasch", {
  fit <- .mock_fit("rasch", "normal")
  expect_silent(plot(fit, type = "info", engine = "base"))
})

test_that("base R: info plot works for 2PL with item info", {
  fit <- .mock_fit("2pl", "normal")
  expect_silent(plot(fit, type = "info", engine = "base",
                     show_items = TRUE))
})

test_that("base R: info plot works for 3PL", {
  fit <- .mock_fit("3pl", "normal")
  expect_silent(plot(fit, type = "info", engine = "base"))
})

test_that("base R: pp_check works (prop_correct)", {
  fit <- .mock_fit("rasch", "normal")
  expect_silent(plot(fit, type = "pp_check", engine = "base",
                     stat = "prop_correct", n_rep = 10))
})

test_that("base R: pp_check works (total_score)", {
  fit <- .mock_fit("rasch", "normal")
  expect_silent(plot(fit, type = "pp_check", engine = "base",
                     stat = "total_score", n_rep = 10))
})


# ============================================================================
# Tests: S3 methods for other classes
# ============================================================================

test_that("plot.dpmirt_estimates: estimates type works", {
  est <- .mock_estimates()
  expect_silent(plot(est, type = "estimates", param = "theta"))
})

test_that("plot.dpmirt_estimates: shrinkage type works", {
  est <- .mock_estimates()
  expect_silent(plot(est, type = "shrinkage", param = "theta"))
})

test_that("plot.dpmirt_estimates: beta param works", {
  est <- .mock_estimates()
  expect_silent(plot(est, type = "estimates", param = "beta"))
})

test_that("plot.dpmirt_sim: parameters works for Rasch", {
  sim <- .mock_sim("rasch")
  expect_silent(plot(sim, type = "parameters"))
})

test_that("plot.dpmirt_sim: parameters works for 2PL", {
  sim <- .mock_sim("2pl")
  expect_silent(plot(sim, type = "parameters"))
})

test_that("plot.dpmirt_sim: response works", {
  sim <- .mock_sim("rasch")
  expect_silent(plot(sim, type = "response"))
})


# ============================================================================
# Tests: ggplot2 engine (skip if not installed)
# ============================================================================

test_that("ggplot2: density plot returns ggplot", {
  skip_if_not_installed("ggplot2")
  fit <- .mock_fit("rasch", "normal")
  p <- dpmirt_plot_density(fit)
  expect_s3_class(p, "ggplot")
})

test_that("ggplot2: items plot returns ggplot", {
  skip_if_not_installed("ggplot2")
  fit <- .mock_fit("rasch", "normal")
  p <- dpmirt_plot_items(fit)
  expect_s3_class(p, "ggplot")
})

test_that("ggplot2: trace plot returns ggplot", {
  skip_if_not_installed("ggplot2")
  fit <- .mock_fit("rasch", "normal")
  p <- dpmirt_plot_trace(fit)
  expect_s3_class(p, "ggplot")
})

test_that("ggplot2: ICC plot returns ggplot for all models", {
  skip_if_not_installed("ggplot2")
  for (m in c("rasch", "2pl", "3pl")) {
    fit <- .mock_fit(m, "normal")
    p <- dpmirt_plot_icc(fit, items = 1:3)
    expect_s3_class(p, "ggplot")
  }
})

test_that("ggplot2: wright_map returns ggplot", {
  skip_if_not_installed("ggplot2")
  fit <- .mock_fit("rasch", "normal")
  p <- dpmirt_plot_wright_map(fit)
  expect_s3_class(p, "ggplot")
})

test_that("ggplot2: parameter_trace returns ggplot", {
  skip_if_not_installed("ggplot2")
  fit <- .mock_fit("rasch", "normal")
  p <- dpmirt_plot_parameter_trace(fit, param = "beta", indices = 1:2)
  expect_s3_class(p, "ggplot")
})

test_that("ggplot2: caterpillar returns ggplot", {
  skip_if_not_installed("ggplot2")
  fit <- .mock_fit("rasch", "normal")
  p <- dpmirt_plot_caterpillar(fit, param = "beta")
  expect_s3_class(p, "ggplot")
})

test_that("ggplot2: density_compare returns ggplot", {
  skip_if_not_installed("ggplot2")
  fit <- .mock_fit("rasch", "dpm")
  p <- dpmirt_plot_density_compare(fit)
  expect_s3_class(p, "ggplot")
})

test_that("ggplot2: info plot returns ggplot", {
  skip_if_not_installed("ggplot2")
  fit <- .mock_fit("2pl", "normal")
  p <- dpmirt_plot_info(fit, show_items = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("ggplot2: pp_check returns ggplot", {
  skip_if_not_installed("ggplot2")
  fit <- .mock_fit("rasch", "normal")
  p <- dpmirt_plot_pp_check(fit, stat = "prop_correct", n_rep = 5)
  expect_s3_class(p, "ggplot")
})

test_that("ggplot2: clusters plot returns ggplot", {
  skip_if_not_installed("ggplot2")
  fit <- .mock_fit("rasch", "dpm")
  p <- dpmirt_plot_clusters(fit)
  expect_s3_class(p, "ggplot")
})

test_that("ggplot2: dp_density plot returns ggplot", {
  skip_if_not_installed("ggplot2")
  fit <- .mock_fit("rasch", "dpm")
  p <- dpmirt_plot_dp_density(fit)
  expect_s3_class(p, "ggplot")
})


# ============================================================================
# Tests: engine dispatch
# ============================================================================

test_that("engine='auto' dispatches correctly", {
  fit <- .mock_fit("rasch", "normal")
  # Should not error regardless of ggplot2 availability
  expect_no_error(plot(fit, type = "density"))
})

test_that("engine='base' forces base R", {
  fit <- .mock_fit("rasch", "normal")
  result <- plot(fit, type = "density", engine = "base")
  expect_null(result)
})

test_that("engine='ggplot2' returns ggplot object", {
  skip_if_not_installed("ggplot2")
  fit <- .mock_fit("rasch", "normal")
  result <- plot(fit, type = "density", engine = "ggplot2")
  expect_s3_class(result, "ggplot")
})


# ============================================================================
# Tests: ICC helper
# ============================================================================

test_that(".icc_prob returns correct values for Rasch", {
  # At theta = beta, P should be 0.5 (for Rasch with delta=0)
  expect_equal(.icc_prob(0, 0), 0.5)
  # Far positive theta -> P ~ 1
  expect_gt(.icc_prob(5, 0), 0.99)
  # Far negative theta -> P ~ 0
  expect_lt(.icc_prob(-5, 0), 0.01)
})

test_that(".icc_prob handles 3PL guessing", {
  # With delta=0.25, minimum P should be 0.25
  expect_gt(.icc_prob(-10, 0, lambda = 1, delta = 0.25), 0.24)
  # At theta = beta with guessing, P > 0.5
  expect_gt(.icc_prob(0, 0, lambda = 1, delta = 0.25), 0.5)
})
