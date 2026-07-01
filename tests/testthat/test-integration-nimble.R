# ============================================================================
# Tiny real NIMBLE integration smoke tests
# ============================================================================
# These tests are opt-in because real NIMBLE compilation dominates runtime.
# Run locally with:
#   NOT_CRAN=true DPMIRT_RUN_NIMBLE_INTEGRATION=true devtools::test(filter = "integration-nimble")
# ============================================================================

.skip_if_not_nimble_integration <- function() {
  skip_on_cran()
  skip_if_not_installed("nimble", minimum_version = "1.0.0")

  run_integration <- tolower(Sys.getenv("DPMIRT_RUN_NIMBLE_INTEGRATION")) %in%
    c("1", "true", "yes")
  skip_if_not(
    run_integration,
    "Set DPMIRT_RUN_NIMBLE_INTEGRATION=true to run real NIMBLE integration tests."
  )
}


.tiny_three_item_response <- function() {
  matrix(
    c(
      0, 0, 0,
      0, 0, 1,
      0, 1, 0,
      0, 1, 1,
      1, 0, 0,
      1, 0, 1,
      1, 1, 0,
      1, 1, 1
    ),
    nrow = 8,
    byrow = TRUE
  )
}


.fit_tiny_rasch_dpm <- function() {
  fit <- NULL
  invisible(capture.output({
    fit <- suppressMessages(dpmirt(
      .tiny_three_item_response(),
      model = "rasch",
      prior = "dpm",
      M = 8L,
      alpha_prior = c(1, 3),
      niter = 20L,
      nburnin = 10L,
      thin = 1L,
      thin2 = 1L,
      nchains = 1L,
      seed = 20260201L,
      rescale = TRUE,
      compute_waic = FALSE,
      compute_dp_density = FALSE,
      verbose = FALSE
    ))
  }))
  fit
}


.fit_tiny_rasch_normal_waic <- function() {
  fit <- NULL
  invisible(capture.output({
    fit <- suppressMessages(dpmirt(
      .tiny_three_item_response(),
      model = "rasch",
      prior = "normal",
      identification = "unconstrained",
      niter = 20L,
      nburnin = 10L,
      thin = 1L,
      thin2 = 1L,
      nchains = 1L,
      seed = 20260203L,
      rescale = TRUE,
      compute_waic = TRUE,
      compute_dp_density = FALSE,
      verbose = FALSE
    ))
  }))
  fit
}


.fit_tiny_rasch_normal_two_chain <- function(seed = 20260205L) {
  fit <- NULL
  invisible(capture.output({
    fit <- suppressMessages(dpmirt(
      .tiny_three_item_response(),
      model = "rasch",
      prior = "normal",
      identification = "unconstrained",
      niter = 16L,
      nburnin = 8L,
      thin = 1L,
      thin2 = 1L,
      nchains = 2L,
      seed = seed,
      rescale = TRUE,
      compute_waic = FALSE,
      compute_dp_density = FALSE,
      verbose = FALSE
    ))
  }))
  fit
}


test_that("dpmirt runs a tiny Rasch-DPM NIMBLE fit end to end", {
  .skip_if_not_nimble_integration()

  fit <- .fit_tiny_rasch_dpm()

  expect_s3_class(fit, "dpmirt_fit")
  expect_equal(fit$config$model, "rasch")
  expect_equal(fit$config$prior, "dpm")
  expect_equal(fit$config$identification, "constrained_item")
  expect_equal(fit$config$N, 8L)
  expect_equal(fit$config$I, 3L)
  expect_equal(fit$config$M, 8L)

  expect_equal(dim(fit$theta_samp), c(10L, 8L))
  expect_equal(dim(fit$beta_samp), c(10L, 3L))
  expect_true(all(is.finite(fit$theta_samp)))
  expect_true(all(is.finite(fit$beta_samp)))
  expect_lt(max(abs(rowMeans(fit$beta_samp))), 1e-6)

  expect_true(all(
    c("alpha", "zi[1]", "muTilde[1]", "s2Tilde[1]") %in%
      colnames(fit$samples_raw)
  ))
  expect_true(all(
    c("alpha", "zi[1]", "muTilde[1]", "s2Tilde[1]") %in%
      colnames(fit$other_samp)
  ))
  expect_false(is.null(fit$cluster_info))
  expect_equal(length(fit$cluster_info$n_clusters), 10L)
  expect_true(all(fit$cluster_info$n_clusters >= 1L))
  expect_true(all(fit$cluster_info$n_clusters <= 8L))

  expect_equal(nrow(fit$draw_index$main), 10L)
  expect_equal(fit$draw_index$main$mcmc_iteration, 11:20)

  expect_null(fit$waic)
  expect_null(fit$dp_density)

  dpd <- NULL
  invisible(capture.output({
    dpd <- suppressMessages(dpmirt_dp_density(
      fit,
      grid = c(-1, 0, 1),
      verbose = FALSE
    ))
  }))

  expect_s3_class(dpd, "dpmirt_dp_density")
  expect_equal(length(dpd$dp_samples), nrow(fit$samples_raw))
  expect_equal(dim(dpd$density_samples), c(nrow(fit$samples_raw), 3L))
  expect_true(all(is.finite(dpd$density_samples)))
  expect_true(all(dpd$density_samples >= 0))
  expect_true(all(vapply(dpd$dp_samples, ncol, integer(1)) >= 3L))
  expect_true(all(vapply(dpd$dp_samples, nrow, integer(1)) >= 1L))
})


test_that("tiny Rasch-Normal two-chain fit uses independent reproducible inits", {
  .skip_if_not_nimble_integration()

  fit <- .fit_tiny_rasch_normal_two_chain(seed = 20260205L)
  fit_again <- .fit_tiny_rasch_normal_two_chain(seed = 20260205L)

  expect_s3_class(fit, "dpmirt_fit")
  expect_equal(fit$config$nchains, 2L)
  expect_equal(dim(fit$theta_samp), c(16L, 8L))
  expect_equal(dim(fit$beta_samp), c(16L, 3L))
  expect_true(all(is.finite(fit$theta_samp)))
  expect_true(all(is.finite(fit$beta_samp)))

  expect_equal(nrow(fit$chain_info), 2L)
  expect_equal(fit$chain_info$chain_id, 1:2)
  expect_equal(fit$chain_info$seed, c(20260205L, 20260206L))
  expect_equal(fit$chain_info$init_seed, c(20260205L, 20260206L))
  expect_equal(fit$chain_info$init_strategy, rep("chain_seeded", 2L))
  expect_equal(fit$diagnostics$waic$aggregation, "mean_of_chain_waic")

  expect_equal(nrow(fit$draw_index$main), 16L)
  expect_equal(nrow(fit$draw_index$theta), 16L)
  expect_equal(fit$draw_index$main$chain_id[1:8], rep(1L, 8L))
  expect_equal(fit$draw_index$main$chain_id[9:16], rep(2L, 8L))
  expect_false(identical(
    fit$samples_raw[1:8, , drop = FALSE],
    fit$samples_raw[9:16, , drop = FALSE]
  ))

  expect_equal(fit$samples_raw, fit_again$samples_raw, tolerance = 1e-10)
  expect_equal(fit$theta_samp, fit_again$theta_samp, tolerance = 1e-10)
  expect_equal(fit$beta_samp, fit_again$beta_samp, tolerance = 1e-10)
})


test_that("tiny Rasch-Normal fit exposes WAIC, ESS, and rescale invariants", {
  .skip_if_not_nimble_integration()

  fit <- .fit_tiny_rasch_normal_waic()

  expect_s3_class(fit, "dpmirt_fit")
  expect_equal(fit$config$model, "rasch")
  expect_equal(fit$config$prior, "normal")
  expect_equal(fit$config$identification, "unconstrained")
  expect_equal(fit$config$N, 8L)
  expect_equal(fit$config$I, 3L)

  expect_equal(dim(fit$theta_samp), c(10L, 8L))
  expect_equal(dim(fit$beta_samp), c(10L, 3L))
  expect_true(all(is.finite(fit$theta_samp)))
  expect_true(all(is.finite(fit$beta_samp)))
  expect_lt(max(abs(rowMeans(fit$beta_samp))), 1e-6)
  expect_true(all(fit$scale_shift == 1))
  expect_true(all(is.finite(fit$location_shift)))
  expect_equal(length(fit$location_shift), nrow(fit$beta_samp))

  expect_true(is.numeric(fit$waic))
  expect_length(fit$waic, 1L)
  expect_true(is.finite(fit$waic))
  expect_equal(fit$diagnostics$waic$value, fit$waic)
  expect_equal(fit$diagnostics$waic$aggregation, "single_chain")
  expect_equal(fit$chain_info$waic, fit$waic, tolerance = 1e-8)

  expect_equal(length(fit$ess$items), 3L)
  expect_equal(length(fit$ess$theta), 8L)
  expect_true(all(is.finite(fit$ess$items)))
  expect_true(all(is.finite(fit$ess$theta)))
  expect_true(all(fit$ess$items > 0))
  expect_true(all(fit$ess$theta > 0))
  expect_identical(fit$diagnostics$ess$pooled, fit$ess)
  expect_equal(nrow(fit$diagnostics$ess$by_chain$items), 3L)
  expect_equal(nrow(fit$diagnostics$ess$by_chain$theta), 8L)

  expect_equal(length(fit$loglik_trace), 10L)
  expect_equal(nrow(fit$diagnostics$loglik$by_chain), 10L)
  expect_null(fit$diagnostics$rhat)
  expect_null(fit$cluster_info)
  expect_null(fit$dp_density)
})
