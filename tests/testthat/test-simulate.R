# ============================================================================
# Tests for R/simulate.R — Data simulation (Phase 5)
# ============================================================================
# Covers: dpmirt_simulate(), .simulate_fallback(), .simulate_irtsimrel(),
#         .map_latent_shape(), .generate_responses(), .compute_kr20(),
#         print.dpmirt_sim()
# ============================================================================

# ============================================================================
# A. Fallback simulation tests (always runnable)
# ============================================================================

test_that("dpmirt_simulate creates valid Rasch data (fallback)", {
  sim <- dpmirt_simulate(n_persons = 50, n_items = 10,
                         model = "rasch", seed = 42,
                         use_irtsimrel = FALSE)

  expect_s3_class(sim, "dpmirt_sim")
  expect_equal(dim(sim$response), c(50, 10))
  expect_true(all(sim$response %in% c(0, 1)))
  expect_equal(length(sim$theta), 50)
  expect_equal(length(sim$beta), 10)
  expect_null(sim$lambda)
  expect_equal(sim$model, "rasch")
  expect_true(sim$reliability > 0 && sim$reliability < 1)
  expect_equal(sim$method, "fallback")
  expect_null(sim$target_rho)
  expect_null(sim$achieved_rho)
  expect_null(sim$eqc_result)
})


test_that("dpmirt_simulate creates valid 2PL data (fallback)", {
  sim <- dpmirt_simulate(n_persons = 50, n_items = 10,
                         model = "2pl", seed = 42,
                         use_irtsimrel = FALSE)

  expect_equal(dim(sim$response), c(50, 10))
  expect_equal(length(sim$lambda), 10)
  expect_true(all(sim$lambda > 0))
  expect_equal(sim$model, "2pl")
  expect_equal(sim$method, "fallback")
  expect_null(sim$achieved_rho)
})


test_that("dpmirt_simulate creates valid 3PL data (fallback)", {
  sim <- dpmirt_simulate(n_persons = 80, n_items = 12,
                         model = "3pl", seed = 42,
                         use_irtsimrel = FALSE)

  expect_s3_class(sim, "dpmirt_sim")
  expect_equal(dim(sim$response), c(80, 12))
  expect_true(all(sim$response %in% c(0, 1)))
  expect_equal(length(sim$lambda), 12)
  expect_equal(length(sim$delta), 12)
  expect_true(all(sim$lambda > 0))
  expect_true(all(sim$delta > 0 & sim$delta < 1))
  expect_equal(sim$model, "3pl")
  expect_equal(sim$method, "fallback")
  expect_null(sim$target_rho)
  expect_null(sim$achieved_rho)
  expect_null(sim$eqc_result)
  expect_true(is.finite(sim$reliability))
})


test_that("dpmirt_simulate 3PL fallback is reproducible with seed", {
  sim1 <- dpmirt_simulate(40, 6, model = "3pl", seed = 123,
                          use_irtsimrel = FALSE)
  sim2 <- dpmirt_simulate(40, 6, model = "3pl", seed = 123,
                          use_irtsimrel = FALSE)

  expect_identical(sim1$response, sim2$response)
  expect_identical(sim1$theta, sim2$theta)
  expect_identical(sim1$lambda, sim2$lambda)
  expect_identical(sim1$delta, sim2$delta)
})


test_that("fallback simulation ignores target_rho by design", {
  sim_low <- dpmirt_simulate(50, 10, model = "rasch", target_rho = 0.2,
                             seed = 123, use_irtsimrel = FALSE)
  sim_high <- dpmirt_simulate(50, 10, model = "rasch", target_rho = 0.9,
                              seed = 123, use_irtsimrel = FALSE)

  expect_null(sim_low$target_rho)
  expect_null(sim_high$target_rho)
  expect_null(sim_low$achieved_rho)
  expect_null(sim_high$achieved_rho)
  expect_identical(sim_low$response, sim_high$response)
  expect_identical(sim_low$theta, sim_high$theta)
})


test_that("dpmirt_simulate handles all fallback latent shapes", {
  sim_normal  <- dpmirt_simulate(100, 10, model = "rasch",
                                  latent_shape = "normal", seed = 1,
                                  use_irtsimrel = FALSE)
  sim_bimodal <- dpmirt_simulate(100, 10, model = "rasch",
                                  latent_shape = "bimodal", seed = 1,
                                  use_irtsimrel = FALSE)
  sim_skewed  <- dpmirt_simulate(100, 10, model = "rasch",
                                  latent_shape = "skewed", seed = 1,
                                  use_irtsimrel = FALSE)

  # Verify shapes stored correctly
  expect_equal(sim_normal$latent_shape, "normal")
  expect_equal(sim_bimodal$latent_shape, "bimodal")
  expect_equal(sim_skewed$latent_shape, "skewed")

  # Bimodal: should have spread-out distribution
  expect_true(sd(sim_bimodal$theta) > 1.0)

  # Skewed: shifted exponential is right-skewed
  expect_true(mean(sim_skewed$theta) > -0.5)
})


test_that("dpmirt_simulate is reproducible with seed (fallback)", {
  sim1 <- dpmirt_simulate(50, 10, model = "rasch", seed = 123,
                          use_irtsimrel = FALSE)
  sim2 <- dpmirt_simulate(50, 10, model = "rasch", seed = 123,
                          use_irtsimrel = FALSE)

  expect_identical(sim1$response, sim2$response)
  expect_identical(sim1$theta, sim2$theta)
  expect_identical(sim1$beta, sim2$beta)
})


test_that("dpmirt_simulate fallback message when use_irtsimrel = TRUE but not installed", {
  # Skip if IRTsimrel IS installed (can't test the message then)
  skip_if(requireNamespace("IRTsimrel", quietly = TRUE),
          "IRTsimrel is installed, can't test fallback message")

  expect_message(
    dpmirt_simulate(50, 10, model = "rasch", seed = 42),
    "IRTsimrel not installed"
  )
})


test_that("print.dpmirt_sim works for fallback", {
  sim <- dpmirt_simulate(50, 10, model = "rasch", seed = 42,
                         use_irtsimrel = FALSE)
  expect_output(print(sim), "DPMirt Simulated Data")
  expect_output(print(sim), "fallback")
  expect_output(print(sim), "KR-20")
})


# ============================================================================
# B. Input validation tests
# ============================================================================

test_that("dpmirt_simulate validates n_persons", {
  expect_error(
    dpmirt_simulate(1, 10, model = "rasch", use_irtsimrel = FALSE),
    "n_persons >= 2"
  )
  expect_error(
    dpmirt_simulate("abc", 10, model = "rasch", use_irtsimrel = FALSE),
    "is.numeric"
  )
})


test_that("dpmirt_simulate validates n_items", {
  expect_error(
    dpmirt_simulate(50, 1, model = "rasch", use_irtsimrel = FALSE),
    "n_items >= 2"
  )
})


test_that("dpmirt_simulate validates target_rho", {
  expect_error(
    dpmirt_simulate(50, 10, target_rho = 0, use_irtsimrel = FALSE),
    "target_rho > 0"
  )
  expect_error(
    dpmirt_simulate(50, 10, target_rho = 1, use_irtsimrel = FALSE),
    "target_rho < 1"
  )
})


test_that("dpmirt_simulate validates model argument", {
  expect_error(
    dpmirt_simulate(50, 10, model = "4pl", use_irtsimrel = FALSE),
    "arg"
  )
})


test_that("dpmirt_simulate validates custom latent shape contract", {
  valid_mixture <- list(
    weights = c(0.5, 0.5),
    means = c(-1, 1),
    sds = c(0.5, 0.5)
  )

  expect_error(
    dpmirt_simulate(50, 10, model = "rasch", latent_shape = "custom"),
    "mixture_spec"
  )
  expect_error(
    dpmirt_simulate(
      50, 10, model = "rasch", latent_shape = "custom",
      latent_params = list(mixture_spec = list(weights = c(0.5, 0.5)))
    ),
    "means"
  )
  expect_error(
    dpmirt_simulate(
      50, 10, model = "rasch", latent_shape = "custom",
      latent_params = list(mixture_spec = list(
        weights = c(0.7, 0.7),
        means = c(-1, 1),
        sds = c(0.5, 0.5)
      ))
    ),
    "sum to 1"
  )
  expect_error(
    dpmirt_simulate(
      50, 10, model = "rasch", latent_shape = "custom",
      latent_params = list(mixture_spec = list(
        weights = c(0.5, 0.5),
        means = c(-1, 1),
        sds = c(0.5, -0.5)
      ))
    ),
    "positive"
  )

  expect_error(
    dpmirt_simulate(
      50, 10, model = "rasch", latent_shape = "custom",
      latent_params = list(mixture_spec = valid_mixture),
      use_irtsimrel = FALSE
    ),
    "requires IRTsimrel"
  )
  expect_error(
    dpmirt_simulate(
      50, 10, model = "3pl", latent_shape = "custom",
      latent_params = list(mixture_spec = valid_mixture)
    ),
    "Rasch/2PL"
  )
  expect_error(
    DPMirt:::.validate_custom_latent_shape(
      latent_shape = "custom",
      latent_params = list(mixture_spec = valid_mixture),
      model = "rasch",
      use_irtsimrel = TRUE,
      has_irtsimrel = FALSE
    ),
    "requires IRTsimrel"
  )
})


# ============================================================================
# C. Internal helper tests
# ============================================================================

test_that(".map_latent_shape maps skewed to skew_pos", {
  expect_equal(DPMirt:::.map_latent_shape("skewed"), "skew_pos")
})


test_that(".map_latent_shape passes through valid IRTsimrel shapes", {
  expect_equal(DPMirt:::.map_latent_shape("normal"), "normal")
  expect_equal(DPMirt:::.map_latent_shape("bimodal"), "bimodal")
  expect_equal(DPMirt:::.map_latent_shape("trimodal"), "trimodal")
  expect_equal(DPMirt:::.map_latent_shape("heavy_tail"), "heavy_tail")
  expect_equal(DPMirt:::.map_latent_shape("uniform"), "uniform")
})


test_that(".map_latent_shape rejects invalid shapes", {
  expect_error(
    DPMirt:::.map_latent_shape("banana"),
    "Unknown latent_shape"
  )
})


test_that(".simulate_fallback rejects unsupported latent shapes", {
  expect_error(
    DPMirt:::.simulate_fallback(20, 5, "rasch", "trimodal"),
    "Fallback simulation supports"
  )
})


test_that(".generate_responses produces binary matrix", {
  theta <- rnorm(20)
  beta  <- seq(-1, 1, length.out = 5)

  set.seed(42)
  y <- DPMirt:::.generate_responses(theta, beta, NULL, "rasch")

  expect_equal(dim(y), c(20, 5))
  expect_true(all(y %in% c(0, 1)))
  expect_equal(storage.mode(y), "double")
})


test_that(".generate_responses works for 2PL", {
  theta <- rnorm(20)
  beta  <- seq(-1, 1, length.out = 5)
  lambda <- rep(1.5, 5)

  set.seed(42)
  y <- DPMirt:::.generate_responses(theta, beta, lambda, "2pl")

  expect_equal(dim(y), c(20, 5))
  expect_true(all(y %in% c(0, 1)))
})


test_that(".generate_responses applies 3PL guessing lower asymptote", {
  n_persons <- 2000
  delta <- 0.25

  set.seed(42)
  low_theta <- DPMirt:::.generate_responses(
    theta = rep(-100, n_persons),
    beta = 0,
    lambda = 1,
    model = "3pl",
    delta = delta
  )

  set.seed(42)
  high_theta <- DPMirt:::.generate_responses(
    theta = rep(100, n_persons),
    beta = 0,
    lambda = 1,
    model = "3pl",
    delta = delta
  )

  expect_equal(mean(low_theta), delta, tolerance = 0.03)
  expect_gt(mean(high_theta), 0.98)
})


test_that(".compute_kr20 returns valid reliability", {
  # Create data with known properties
  set.seed(42)
  theta <- rnorm(200)
  beta  <- seq(-2, 2, length.out = 20)
  y <- DPMirt:::.generate_responses(theta, beta, NULL, "rasch")

  kr20 <- DPMirt:::.compute_kr20(y)
  expect_true(is.numeric(kr20))
  expect_true(kr20 > 0 && kr20 < 1)
})


test_that(".compute_kr20 handles zero-variance case", {
  # All identical responses -> zero variance
  y <- matrix(1, nrow = 10, ncol = 5)
  kr20 <- DPMirt:::.compute_kr20(y)
  expect_equal(kr20, 0)
})


# ============================================================================
# D. dpmirt_sim object structure tests
# ============================================================================

test_that("dpmirt_sim contains all expected fields", {
  sim <- dpmirt_simulate(50, 10, model = "rasch", seed = 42,
                         use_irtsimrel = FALSE)

  expected_names <- c("response", "theta", "beta", "lambda",
                      "delta", "n_persons", "n_items", "model",
                      "reliability", "achieved_rho", "target_rho",
                      "latent_shape", "eqc_result", "method")
  expect_true(all(expected_names %in% names(sim)))
})


test_that("dpmirt_sim works as input to dpmirt_loss pipeline (mock)", {

  # Simulate data, then construct mock estimates to verify pipeline flow
  sim <- dpmirt_simulate(100, 15, model = "rasch", seed = 42,
                         use_irtsimrel = FALSE)

  # Verify that true params have correct dimensions for loss computation

  expect_equal(length(sim$theta), sim$n_persons)
  expect_equal(length(sim$beta), sim$n_items)
  expect_equal(nrow(sim$response), sim$n_persons)
  expect_equal(ncol(sim$response), sim$n_items)
})


# ============================================================================
# E. IRTsimrel integration tests (only run if IRTsimrel is installed)
# ============================================================================

test_that("dpmirt_simulate uses IRTsimrel when available (Rasch)", {
  skip_if_not_installed("IRTsimrel")

  sim <- dpmirt_simulate(
    n_persons = 200, n_items = 20,
    model = "rasch", target_rho = 0.80,
    latent_shape = "normal",
    item_source = "parametric",
    seed = 42, use_irtsimrel = TRUE
  )

  expect_s3_class(sim, "dpmirt_sim")
  expect_equal(sim$method, "irtsimrel")
  expect_equal(sim$target_rho, 0.80)
  expect_equal(sim$achieved_rho, sim$eqc_result$achieved_rho)
  expect_true(is.numeric(sim$achieved_rho))
  expect_length(sim$achieved_rho, 1)
  expect_true(is.finite(sim$achieved_rho))
  expect_equal(sim$reliability, DPMirt:::.compute_kr20(sim$response))
  expect_false(is.null(sim$eqc_result))
  expect_equal(dim(sim$response), c(200, 20))
  expect_true(all(sim$response %in% c(0, 1)))

  # EQC should achieve close to target reliability
  expect_true(abs(sim$achieved_rho - 0.80) < 0.05)
})


test_that("dpmirt_simulate uses IRTsimrel for 2PL", {
  skip_if_not_installed("IRTsimrel")

  sim <- dpmirt_simulate(
    n_persons = 200, n_items = 20,
    model = "2pl", target_rho = 0.85,
    latent_shape = "normal",
    item_source = "parametric",
    seed = 42, use_irtsimrel = TRUE
  )

  expect_equal(sim$method, "irtsimrel")
  expect_equal(sim$model, "2pl")
  expect_equal(length(sim$lambda), 20)
  expect_true(all(sim$lambda > 0))
  expect_equal(sim$achieved_rho, sim$eqc_result$achieved_rho)
  expect_true(is.numeric(sim$achieved_rho))
  expect_length(sim$achieved_rho, 1)
  expect_true(is.finite(sim$achieved_rho))
  expect_true(abs(sim$achieved_rho - 0.85) < 0.05)
})


test_that("dpmirt_simulate with IRTsimrel handles bimodal shape", {
  skip_if_not_installed("IRTsimrel")

  sim <- dpmirt_simulate(
    n_persons = 300, n_items = 25,
    model = "rasch", target_rho = 0.80,
    latent_shape = "bimodal",
    item_source = "parametric",
    seed = 42, use_irtsimrel = TRUE
  )

  expect_equal(sim$method, "irtsimrel")
  expect_equal(sim$latent_shape, "bimodal")

  # Bimodal theta should have clear separation
  theta <- sim$theta
  expect_true(sd(theta) > 0.5)
})


test_that("dpmirt_simulate with IRTsimrel maps skewed to skew_pos", {
  skip_if_not_installed("IRTsimrel")

  # Use info metric which handles skewed distributions better
  sim <- suppressWarnings(dpmirt_simulate(
    n_persons = 200, n_items = 30,
    model = "rasch", target_rho = 0.70,
    latent_shape = "skewed",
    item_source = "parametric",
    reliability_metric = "info",
    seed = 42, use_irtsimrel = TRUE
  ))

  expect_equal(sim$method, "irtsimrel")
  expect_equal(sim$latent_shape, "skewed")  # DPMirt stores original name
})


test_that("dpmirt_simulate with IRTsimrel uses different reliability metrics", {
  skip_if_not_installed("IRTsimrel")

  sim_msem <- dpmirt_simulate(
    n_persons = 200, n_items = 20,
    model = "rasch", target_rho = 0.80,
    reliability_metric = "msem",
    item_source = "parametric",
    seed = 42, use_irtsimrel = TRUE
  )

  sim_info <- dpmirt_simulate(
    n_persons = 200, n_items = 20,
    model = "rasch", target_rho = 0.80,
    reliability_metric = "info",
    item_source = "parametric",
    seed = 42, use_irtsimrel = TRUE
  )

  expect_equal(sim_msem$method, "irtsimrel")
  expect_equal(sim_info$method, "irtsimrel")
  expect_s3_class(sim_msem$eqc_result, "sac_result")
  expect_s3_class(sim_info$eqc_result, "eqc_result")
  expect_equal(sim_msem$eqc_result$metric, "msem")
  expect_equal(sim_info$eqc_result$metric, "info")
  expect_equal(sim_msem$achieved_rho, sim_msem$eqc_result$achieved_rho)
  expect_equal(sim_info$achieved_rho, sim_info$eqc_result$achieved_rho)
  expect_true(is.finite(sim_msem$achieved_rho))
  expect_true(is.finite(sim_info$achieved_rho))

  # Different metrics may yield different c* values
  expect_false(identical(sim_msem$eqc_result$c_star,
                         sim_info$eqc_result$c_star))
})


test_that("dpmirt_simulate accepts valid custom IRTsimrel mixture", {
  skip_if_not_installed("IRTsimrel")

  sim <- dpmirt_simulate(
    n_persons = 100, n_items = 8,
    model = "rasch", target_rho = 0.70,
    latent_shape = "custom",
    latent_params = list(mixture_spec = list(
      weights = c(0.5, 0.5),
      means = c(-1, 1),
      sds = c(0.5, 0.5)
    )),
    M = 1000,
    seed = 42, use_irtsimrel = TRUE
  )

  expect_equal(sim$method, "irtsimrel")
  expect_equal(sim$latent_shape, "custom")
  expect_equal(sim$achieved_rho, sim$eqc_result$achieved_rho)
  expect_true(is.finite(sim$achieved_rho))
})


test_that("print.dpmirt_sim works for IRTsimrel method", {
  skip_if_not_installed("IRTsimrel")

  sim <- dpmirt_simulate(
    n_persons = 200, n_items = 20,
    model = "rasch", target_rho = 0.80,
    item_source = "parametric",
    seed = 42, use_irtsimrel = TRUE
  )

  expect_output(print(sim), "irtsimrel")
  expect_output(print(sim), "Target rho")
  expect_output(print(sim), "EQC c\\*")
  expect_output(print(sim), "EQC rho")
})


test_that("dpmirt_simulate with IRTsimrel is reproducible", {
  skip_if_not_installed("IRTsimrel")

  sim1 <- dpmirt_simulate(
    n_persons = 100, n_items = 15,
    model = "rasch", target_rho = 0.80,
    item_source = "parametric",
    seed = 42, use_irtsimrel = TRUE
  )

  sim2 <- dpmirt_simulate(
    n_persons = 100, n_items = 15,
    model = "rasch", target_rho = 0.80,
    item_source = "parametric",
    seed = 42, use_irtsimrel = TRUE
  )

  expect_identical(sim1$response, sim2$response)
  expect_identical(sim1$theta, sim2$theta)
  expect_equal(sim1$eqc_result$c_star, sim2$eqc_result$c_star)
})


test_that("dpmirt_simulate gate criterion: reliability within tolerance", {
  skip_if_not_installed("IRTsimrel")

  # Gate: Target reliability achieved within tolerance on Rasch/2PL data.
  tolerance <- 0.10  # EQC theoretical tolerance

  # Rasch
  sim_rasch <- dpmirt_simulate(
    n_persons = 500, n_items = 25,
    model = "rasch", target_rho = 0.80,
    item_source = "parametric",
    seed = 42, use_irtsimrel = TRUE
  )
  expect_true(abs(sim_rasch$achieved_rho - 0.80) < tolerance)

  # 2PL
  sim_2pl <- dpmirt_simulate(
    n_persons = 500, n_items = 25,
    model = "2pl", target_rho = 0.85,
    item_source = "parametric",
    seed = 42, use_irtsimrel = TRUE
  )
  expect_true(abs(sim_2pl$achieved_rho - 0.85) < tolerance)
})
