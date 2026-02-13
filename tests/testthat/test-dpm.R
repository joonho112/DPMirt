# ============================================================================
# Unit Tests for Phase 2: DPM Extension
# ============================================================================
# Tests DPM-specific functionality (model spec, cluster diagnostics, etc.)
# These tests do NOT require NIMBLE (pure R logic testing).
# ============================================================================


# --- Helper: Create test data ---
set.seed(42)
y_test <- matrix(rbinom(50 * 10, 1, 0.5), nrow = 50, ncol = 10)


# ============================================================================
# DPM Model Specification
# ============================================================================

test_that("dpmirt_spec creates Rasch-DPM-constrained_item spec", {
  spec <- dpmirt_spec(y_test, model = "rasch", prior = "dpm",
                       identification = "constrained_item")
  expect_s3_class(spec, "dpmirt_spec")
  expect_equal(spec$config$prior, "dpm")
  expect_equal(spec$config$identification, "constrained_item")
})

test_that("dpmirt_spec creates Rasch-DPM-unconstrained spec", {
  spec <- dpmirt_spec(y_test, model = "rasch", prior = "dpm",
                       identification = "unconstrained")
  expect_s3_class(spec, "dpmirt_spec")
  expect_equal(spec$config$identification, "unconstrained")
})

test_that("DPM spec has correct monitors", {
  spec <- dpmirt_spec(y_test, model = "rasch", prior = "dpm")
  expect_true(all(c("alpha", "zi", "muTilde", "s2Tilde") %in% spec$monitors))
  expect_true("eta" %in% spec$monitors2)
})

test_that("DPM spec has correct constants", {
  spec <- dpmirt_spec(y_test, model = "rasch", prior = "dpm")
  expect_true(!is.null(spec$constants$M))
  expect_true(!is.null(spec$constants$a))
  expect_true(!is.null(spec$constants$b))
  expect_true(!is.null(spec$constants$s2_mu))
  expect_true(!is.null(spec$constants$nu1))
  expect_true(!is.null(spec$constants$nu2))
})

test_that("DPM spec has DPM-specific inits", {
  spec <- dpmirt_spec(y_test, model = "rasch", prior = "dpm")
  expect_true(!is.null(spec$inits$zi))
  expect_true(!is.null(spec$inits$alpha))
  expect_true(!is.null(spec$inits$muTilde))
  expect_true(!is.null(spec$inits$s2Tilde))
  expect_length(spec$inits$zi, 50)
  expect_length(spec$inits$muTilde, spec$constants$M)
})

test_that("DPM spec default M = 50", {
  spec <- dpmirt_spec(y_test, model = "rasch", prior = "dpm")
  expect_equal(spec$constants$M, 50)
})

test_that("DPM spec custom M", {
  spec <- dpmirt_spec(y_test, model = "rasch", prior = "dpm", M = 30L)
  expect_equal(spec$constants$M, 30)
})


# ============================================================================
# Alpha Prior Resolution
# ============================================================================

test_that("NULL alpha_prior gives Gamma(1, 3)", {
  spec <- dpmirt_spec(y_test, model = "rasch", prior = "dpm",
                       alpha_prior = NULL)
  expect_equal(unname(spec$config$alpha_prior["a"]), 1)
  expect_equal(unname(spec$config$alpha_prior["b"]), 3)
})

test_that("Numeric c(2, 4) alpha_prior", {
  spec <- dpmirt_spec(y_test, model = "rasch", prior = "dpm",
                       alpha_prior = c(2, 4))
  expect_equal(unname(spec$config$alpha_prior["a"]), 2)
  expect_equal(unname(spec$config$alpha_prior["b"]), 4)
})

test_that("DPprior_fit list alpha_prior", {
  fake_fit <- list(a = 1.5, b = 2.5)
  spec <- dpmirt_spec(y_test, model = "rasch", prior = "dpm",
                       alpha_prior = fake_fit)
  expect_equal(unname(spec$config$alpha_prior["a"]), 1.5)
  expect_equal(unname(spec$config$alpha_prior["b"]), 2.5)
})

test_that("DPprior_fit object with class accepted as alpha_prior", {
  fake_fit <- list(a = 2.0, b = 1.0, J = 200, method = "A2-MN",
                   target = list(mu_K = 5))
  class(fake_fit) <- "DPprior_fit"
  spec <- dpmirt_spec(y_test, model = "rasch", prior = "dpm",
                       alpha_prior = fake_fit)
  expect_equal(unname(spec$config$alpha_prior["a"]), 2.0)
  expect_equal(unname(spec$config$alpha_prior["b"]), 1.0)
})

test_that("Invalid alpha_prior errors", {
  expect_error(
    dpmirt_spec(y_test, model = "rasch", prior = "dpm", alpha_prior = "bad"),
    "Invalid alpha_prior"
  )
})


# ============================================================================
# DPM Identification Constraints
# ============================================================================

test_that("constrained_ability + DPM is rejected", {
  expect_error(
    dpmirt_spec(y_test, model = "rasch", prior = "dpm",
                 identification = "constrained_ability"),
    "constrained_ability"
  )
})

test_that("DPM default identification is constrained_item", {
  spec <- dpmirt_spec(y_test, model = "rasch", prior = "dpm")
  expect_equal(spec$config$identification, "constrained_item")
})


# ============================================================================
# Cluster Diagnostics
# ============================================================================

test_that(".extract_cluster_info counts unique clusters", {
  fake_samples <- matrix(0, nrow = 100, ncol = 55)
  colnames(fake_samples) <- c(
    "alpha",
    paste0("zi[", 1:50, "]"),
    paste0("muTilde[", 1:2, "]"),
    paste0("s2Tilde[", 1:2, "]")
  )
  fake_samples[, "alpha"] <- rgamma(100, 1, 3)
  for (i in 1:100) {
    fake_samples[i, paste0("zi[", 1:50, "]")] <- sample(1:3, 50, replace = TRUE)
  }

  ci <- .extract_cluster_info(fake_samples, N = 50)
  expect_length(ci$n_clusters, 100)
  expect_true(all(ci$n_clusters >= 1))
  expect_true(all(ci$n_clusters <= 50))
  expect_true(ci$alpha_summary$mean > 0)
})

test_that(".summarize_n_clusters returns correct structure", {
  n_cl <- c(3, 4, 3, 5, 4, 3, 4, 4, 3, 5)
  s <- .summarize_n_clusters(n_cl)
  expect_true(!is.null(s$mean))
  expect_true(!is.null(s$median))
  expect_true(!is.null(s$mode))
  expect_true(!is.null(s$q025))
  expect_true(!is.null(s$q975))
  expect_true(!is.null(s$min))
  expect_true(!is.null(s$max))
})

test_that(".numeric_mode returns correct mode", {
  expect_equal(.numeric_mode(c(1, 2, 2, 3, 3, 3, 4)), 3)
  expect_equal(.numeric_mode(c(5, 5, 5)), 5)
})
