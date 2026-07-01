# ============================================================================
# Tests for R/dp_measure.R
# DP density reconstruction helpers
# ============================================================================

.dp_test_samples <- function(n_draws = 3L, N = 2L, M = 2L) {
  cn <- c(
    "alpha",
    paste0("zi[", seq_len(N), "]"),
    paste0("muTilde[", seq_len(M), "]"),
    paste0("s2Tilde[", seq_len(M), "]")
  )
  out <- matrix(
    1,
    nrow = n_draws,
    ncol = length(cn),
    dimnames = list(NULL, cn)
  )
  out[, "alpha"] <- seq_len(n_draws)
  out
}


.dp_test_fit <- function(samples_raw = .dp_test_samples(),
                         other_samp = samples_raw,
                         N = 2L,
                         M = 2L) {
  structure(
    list(
      config = list(
        prior = "dpm",
        model = "rasch",
        parameterization = "irt",
        identification = "constrained_item",
        N = N,
        M = M
      ),
      samples_raw = samples_raw,
      other_samp = other_samp,
      compiled = list(spec = list(constants = list(N = N, M = M)))
    ),
    class = "dpmirt_fit"
  )
}


test_that("dpmirt_dp_density rejects non-fit and non-DPM inputs early", {
  expect_error(
    dpmirt_dp_density("not a fit", verbose = FALSE),
    "dpmirt_fit"
  )

  fake_fit <- structure(
    list(config = list(prior = "normal")),
    class = "dpmirt_fit"
  )

  expect_error(
    dpmirt_dp_density(fake_fit, verbose = FALSE),
    "only available for DPM"
  )

  malformed_fit <- structure(
    list(config = list()),
    class = "dpmirt_fit"
  )
  expect_error(
    dpmirt_dp_density(malformed_fit, verbose = FALSE),
    "only available for DPM"
  )
})


test_that("dpmirt_dp_density validates arguments before NIMBLE work", {
  fake_fit <- .dp_test_fit()

  expect_error(
    dpmirt_dp_density(fake_fit, grid = numeric(0), verbose = FALSE),
    "grid.*length at least 2"
  )
  expect_error(
    dpmirt_dp_density(fake_fit, grid = c(-1, NA_real_), verbose = FALSE),
    "grid.*length at least 2"
  )
  expect_error(
    dpmirt_dp_density(fake_fit, credible_interval = 1, verbose = FALSE),
    "credible_interval.*between 0 and 1"
  )
  expect_error(
    dpmirt_dp_density(fake_fit, apply_rescaling = NA, verbose = FALSE),
    "apply_rescaling"
  )
  expect_error(
    dpmirt_dp_density(fake_fit, n_grid = 0, verbose = FALSE),
    "Unused argument"
  )
})


test_that(".extract_dp_samples requires all DP monitor columns", {
  fake_fit <- structure(
    list(
      config = list(prior = "dpm"),
      other_samp = matrix(
        1,
        nrow = 3,
        ncol = 2,
        dimnames = list(NULL, c("alpha", "zi[1]"))
      )
    ),
    class = "dpmirt_fit"
  )

  expect_error(
    DPMirt:::.extract_dp_samples(fake_fit),
    "muTilde"
  )
})


test_that(".extract_dp_samples rejects missing or malformed DP sample storage", {
  fake_fit <- structure(
    list(config = list(prior = "dpm"), other_samp = NULL),
    class = "dpmirt_fit"
  )

  expect_error(
    DPMirt:::.extract_dp_samples(fake_fit),
    "Cannot find raw posterior samples"
  )

  fake_fit$other_samp <- matrix(numeric(0), nrow = 0, ncol = 4)
  colnames(fake_fit$other_samp) <- c("alpha", "zi[1]", "muTilde[1]", "s2Tilde[1]")
  expect_error(
    DPMirt:::.extract_dp_samples(fake_fit),
    "at least one draw"
  )

  fake_fit$other_samp <- matrix(1, nrow = 2, ncol = 4)
  expect_error(
    DPMirt:::.extract_dp_samples(fake_fit),
    "column names"
  )

  fake_fit$other_samp <- matrix(
    c("a", "1", "0", "1"),
    nrow = 1,
    dimnames = list(NULL, c("alpha", "zi[1]", "muTilde[1]", "s2Tilde[1]"))
  )
  expect_error(
    DPMirt:::.extract_dp_samples(fake_fit),
    "must be numeric"
  )

  fake_fit$other_samp <- matrix(
    c(1, 1, NA_real_, 1),
    nrow = 1,
    dimnames = list(NULL, c("alpha", "zi[1]", "muTilde[1]", "s2Tilde[1]"))
  )
  expect_error(
    DPMirt:::.extract_dp_samples(fake_fit),
    "finite numeric"
  )

  fake_fit$other_samp <- matrix(
    c(1, 1, 0, Inf),
    nrow = 1,
    dimnames = list(NULL, c("alpha", "zi[1]", "muTilde[1]", "s2Tilde[1]"))
  )
  expect_error(
    DPMirt:::.extract_dp_samples(fake_fit),
    "finite numeric"
  )
})


test_that(".extract_dp_samples falls back to samples_raw when other_samp is absent", {
  fake_fit <- .dp_test_fit(other_samp = NULL)

  expect_equal(DPMirt:::.extract_dp_samples(fake_fit), fake_fit$samples_raw)
})


test_that(".extract_dp_samples uses samples_raw fallback when other_samp is stale", {
  fake_fit <- .dp_test_fit()
  fake_fit$other_samp <- matrix(
    1,
    nrow = 3,
    ncol = 1,
    dimnames = list(NULL, "alpha")
  )

  expect_equal(DPMirt:::.extract_dp_samples(fake_fit), fake_fit$samples_raw)
})


test_that(".extract_dp_samples validates N and M cardinality when known", {
  fit_bad_N <- .dp_test_fit(
    samples_raw = .dp_test_samples(N = 1L, M = 2L),
    other_samp = NULL,
    N = 2L,
    M = 2L
  )
  expect_error(
    DPMirt:::.extract_dp_samples(fit_bad_N),
    "Expected 2 zi"
  )

  fit_bad_mu <- .dp_test_fit(
    samples_raw = .dp_test_samples(N = 2L, M = 1L),
    other_samp = NULL,
    N = 2L,
    M = 2L
  )
  expect_error(
    DPMirt:::.extract_dp_samples(fit_bad_mu),
    "Expected 2 muTilde"
  )

  fit_bad_s2 <- .dp_test_fit(
    samples_raw = .dp_test_samples(N = 2L, M = 2L)[, -6, drop = FALSE],
    other_samp = NULL,
    N = 2L,
    M = 2L
  )
  expect_error(
    DPMirt:::.extract_dp_samples(fit_bad_s2),
    "Expected 2 s2Tilde"
  )
})


test_that(".extract_dp_samples reports each missing DP monitor family", {
  make_fit <- function(cols) {
    structure(
      list(
        config = list(prior = "dpm"),
        other_samp = matrix(1, nrow = 2, ncol = length(cols),
                            dimnames = list(NULL, cols))
      ),
      class = "dpmirt_fit"
    )
  }

  expect_error(
    DPMirt:::.extract_dp_samples(
      make_fit(c("zi[1]", "muTilde[1]", "s2Tilde[1]"))
    ),
    "alpha"
  )
  expect_error(
    DPMirt:::.extract_dp_samples(
      make_fit(c("alpha", "muTilde[1]", "s2Tilde[1]"))
    ),
    "zi"
  )
  expect_error(
    DPMirt:::.extract_dp_samples(
      make_fit(c("alpha", "zi[1]", "s2Tilde[1]"))
    ),
    "muTilde"
  )
  expect_error(
    DPMirt:::.extract_dp_samples(
      make_fit(c("alpha", "zi[1]", "muTilde[1]"))
    ),
    "s2Tilde"
  )
})


test_that(".get_dp_measure_samples errors clearly without compiled spec", {
  fake_fit <- .dp_test_fit()
  fake_fit$compiled <- NULL

  expect_error(
    DPMirt:::.get_dp_measure_samples(
      fake_fit,
      DPMirt:::.extract_dp_samples(.dp_test_fit()),
      verbose = FALSE
    ),
    "Cannot find model specification"
  )
})


test_that(".matrix_to_model_values loads DP sample matrices via public accessors", {
  skip_if_not_installed("nimble")

  mv_conf <- nimble::modelValuesConf(
    vars = c("alpha", "zi", "muTilde", "s2Tilde"),
    types = c("double", "int", "double", "double"),
    sizes = list(alpha = 1, zi = 2, muTilde = 2, s2Tilde = 2)
  )
  mv <- nimble::modelValues(mv_conf, m = 1)
  dp_samples <- .dp_test_samples(n_draws = 3, N = 2, M = 2)

  DPMirt:::.matrix_to_model_values(dp_samples, mv)

  mv_mat <- as.matrix(mv)
  expect_equal(nimble::getsize(mv), 3)
  expect_equal(unname(mv_mat[, "alpha[1]"]), unname(dp_samples[, "alpha"]))
  expect_equal(
    unname(mv_mat[, c("zi[1]", "zi[2]")]),
    unname(dp_samples[, c("zi[1]", "zi[2]")])
  )
  expect_equal(
    unname(mv_mat[, c("muTilde[1]", "muTilde[2]")]),
    unname(dp_samples[, c("muTilde[1]", "muTilde[2]")])
  )
  expect_equal(
    unname(mv_mat[, c("s2Tilde[1]", "s2Tilde[2]")]),
    unname(dp_samples[, c("s2Tilde[1]", "s2Tilde[2]")])
  )
})


test_that(".matrix_to_model_values validates modelValues dimensions", {
  skip_if_not_installed("nimble")

  mv_conf <- nimble::modelValuesConf(
    vars = c("alpha", "zi", "muTilde", "s2Tilde"),
    types = c("double", "int", "double", "double"),
    sizes = list(alpha = 1, zi = 2, muTilde = 2, s2Tilde = 2)
  )
  mv <- nimble::modelValues(mv_conf, m = 1)
  dp_samples <- .dp_test_samples(n_draws = 2, N = 2, M = 2)
  bad_samples <- dp_samples[, colnames(dp_samples) != "zi[2]", drop = FALSE]

  expect_error(
    DPMirt:::.matrix_to_model_values(bad_samples, mv),
    "dimensions for variable 'zi'"
  )
})


test_that(".evaluate_dp_density evaluates weighted normal mixtures", {
  dp_measure <- list(
    matrix(
      c(0.25, -1, 1,
        0.75,  1, 4),
      ncol = 3,
      byrow = TRUE
    )
  )
  grid <- c(-1, 0, 1)
  fake_fit <- structure(
    list(config = list(identification = "constrained_item")),
    class = "dpmirt_fit"
  )

  result <- DPMirt:::.evaluate_dp_density(
    dp_measure = dp_measure,
    grid = grid,
    fit = fake_fit,
    apply_rescaling = TRUE
  )

  expected <- 0.25 * dnorm(grid, mean = -1, sd = 1) +
    0.75 * dnorm(grid, mean = 1, sd = 2)
  expect_equal(as.vector(result), expected, tolerance = 1e-12)
})


test_that(".evaluate_dp_density rejects malformed DP measure draws", {
  fake_fit <- structure(
    list(config = list(identification = "constrained_item")),
    class = "dpmirt_fit"
  )
  grid <- c(-1, 0, 1)

  expect_error(
    DPMirt:::.evaluate_dp_density(list(), grid, fake_fit),
    "non-empty list"
  )
  expect_error(
    DPMirt:::.evaluate_dp_density(list(NULL), grid, fake_fit),
    "matrix-like"
  )
  expect_error(
    DPMirt:::.evaluate_dp_density(list(matrix(numeric(0), nrow = 0, ncol = 3)),
                                  grid, fake_fit),
    "at least one mixture component"
  )
  expect_error(
    DPMirt:::.evaluate_dp_density(list(matrix(c(1, 0), ncol = 2)),
                                  grid, fake_fit),
    "at least three columns"
  )
  expect_error(
    DPMirt:::.evaluate_dp_density(list(matrix(c(1, NA_real_, 1), ncol = 3)),
                                  grid, fake_fit),
    "finite numeric values"
  )
  expect_error(
    DPMirt:::.evaluate_dp_density(list(matrix(c(-1, 0, 1), ncol = 3)),
                                  grid, fake_fit),
    "negative mixture weights"
  )
  expect_error(
    DPMirt:::.evaluate_dp_density(list(matrix(c(0, 0, 1), ncol = 3)),
                                  grid, fake_fit),
    "positive total mixture weight"
  )
  expect_error(
    DPMirt:::.evaluate_dp_density(list(matrix(c(1, 0, 0), ncol = 3)),
                                  grid, fake_fit),
    "non-positive variances"
  )
})


test_that(".evaluate_dp_density applies unconstrained location shift", {
  dp_measure <- list(matrix(c(1, 0, 1), ncol = 3))
  fake_fit <- structure(
    list(
      config = list(identification = "unconstrained", model = "rasch"),
      location_shift = 1
    ),
    class = "dpmirt_fit"
  )

  shifted <- DPMirt:::.evaluate_dp_density(
    dp_measure = dp_measure,
    grid = c(0, 1),
    fit = fake_fit,
    apply_rescaling = TRUE
  )
  unshifted <- DPMirt:::.evaluate_dp_density(
    dp_measure = dp_measure,
    grid = c(0, 1),
    fit = fake_fit,
    apply_rescaling = FALSE
  )

  expect_equal(as.numeric(shifted), dnorm(c(1, 2), mean = 0, sd = 1))
  expect_equal(as.numeric(unshifted), dnorm(c(0, 1), mean = 0, sd = 1))
  expect_false(isTRUE(all.equal(shifted, unshifted)))
})


test_that(".evaluate_dp_density requires location shifts for all DP samples", {
  dp_measure <- list(
    matrix(c(1, 0, 1), ncol = 3),
    matrix(c(1, 1, 1), ncol = 3)
  )
  fake_fit <- structure(
    list(
      config = list(identification = "unconstrained", model = "rasch"),
      location_shift = 1
    ),
    class = "dpmirt_fit"
  )

  expect_error(
    DPMirt:::.evaluate_dp_density(
      dp_measure = dp_measure,
      grid = c(-1, 0, 1),
      fit = fake_fit,
      apply_rescaling = TRUE
    ),
    "location_shift"
  )

  fake_fit$location_shift <- NULL
  expect_error(
    DPMirt:::.evaluate_dp_density(
      dp_measure = dp_measure,
      grid = c(-1, 0, 1),
      fit = fake_fit,
      apply_rescaling = TRUE
    ),
    "location_shift"
  )

  fake_fit$location_shift <- c(0, 1, 2)
  expect_error(
    DPMirt:::.evaluate_dp_density(
      dp_measure = dp_measure,
      grid = c(-1, 0, 1),
      fit = fake_fit,
      apply_rescaling = TRUE
    ),
    "location_shift"
  )

  fake_fit$location_shift <- c(0, NA_real_)
  expect_error(
    DPMirt:::.evaluate_dp_density(
      dp_measure = dp_measure,
      grid = c(-1, 0, 1),
      fit = fake_fit,
      apply_rescaling = TRUE
    ),
    "location_shift"
  )

  fake_fit$location_shift <- c(0, Inf)
  expect_error(
    DPMirt:::.evaluate_dp_density(
      dp_measure = dp_measure,
      grid = c(-1, 0, 1),
      fit = fake_fit,
      apply_rescaling = TRUE
    ),
    "location_shift"
  )

  fake_fit$location_shift <- list(0, 1)
  expect_error(
    DPMirt:::.evaluate_dp_density(
      dp_measure = dp_measure,
      grid = c(-1, 0, 1),
      fit = fake_fit,
      apply_rescaling = TRUE
    ),
    "location_shift"
  )
})


test_that(".dp_percentile evaluates raw DP mixture CDF", {
  dp_density <- list(
    dp_samples = list(matrix(c(1, 0, 1), ncol = 3))
  )
  theta_values <- c(-1, 0, 1)

  result <- DPMirt:::.dp_percentile(dp_density, theta_values)

  expect_equal(as.vector(result), pnorm(theta_values), tolerance = 1e-12)
})


test_that(".dp_percentile rejects malformed DP measure draws", {
  theta_values <- c(-1, 0, 1)

  expect_error(
    DPMirt:::.dp_percentile(list(dp_samples = list()), theta_values),
    "non-empty list"
  )
  expect_error(
    DPMirt:::.dp_percentile(
      list(dp_samples = list(matrix(c(1, 0, 1), ncol = 3))),
      c(0, NA_real_)
    ),
    "finite numeric vector"
  )
  expect_error(
    DPMirt:::.dp_percentile(
      list(dp_samples = list(matrix(numeric(0), nrow = 0, ncol = 3))),
      theta_values
    ),
    "at least one mixture component"
  )
  expect_error(
    DPMirt:::.dp_percentile(
      list(dp_samples = list(matrix(c(1, NA_real_, 1), ncol = 3))),
      theta_values
    ),
    "finite numeric values"
  )
})
