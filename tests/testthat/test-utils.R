# ============================================================================
# Tests for R/utils.R — Input validation and utilities
# ============================================================================

test_that(".validate_model accepts valid models", {
  expect_equal(.validate_model("rasch"), "rasch")
  expect_equal(.validate_model("2pl"), "2pl")
  expect_equal(.validate_model("3pl"), "3pl")
  expect_equal(.validate_model("RASCH"), "rasch")
  expect_equal(.validate_model("Rasch"), "rasch")
})

test_that(".validate_model rejects invalid models", {
  expect_error(.validate_model("4pl"))
  expect_error(.validate_model(""))
})

test_that(".validate_prior accepts valid priors", {
  expect_equal(.validate_prior("normal"), "normal")
  expect_equal(.validate_prior("dpm"), "dpm")
  expect_equal(.validate_prior("DPM"), "dpm")
})

test_that(".validate_parameterization validates correctly", {
  expect_equal(.validate_parameterization("irt", "rasch"), "irt")
  expect_equal(.validate_parameterization("irt", "2pl"), "irt")
  expect_equal(.validate_parameterization("si", "2pl"), "si")
  expect_error(.validate_parameterization("si", "rasch"),
               "not meaningful for")
})

test_that(".resolve_identification applies correct defaults", {
  # Rasch default: constrained_item

  expect_equal(.resolve_identification(NULL, "rasch", "normal"),
               "constrained_item")
  expect_equal(.resolve_identification(NULL, "rasch", "dpm"),
               "constrained_item")

  # 2PL/3PL default: unconstrained
  expect_equal(.resolve_identification(NULL, "2pl", "normal"),
               "unconstrained")
  expect_equal(.resolve_identification(NULL, "3pl", "normal"),
               "unconstrained")
})

test_that(".resolve_identification rejects invalid combinations", {
  # constrained_ability + dpm
  expect_error(.resolve_identification("constrained_ability", "rasch", "dpm"),
               "defeats the purpose")

  # constrained_item + 3pl
  expect_error(.resolve_identification("constrained_item", "3pl", "normal"),
               "not implemented")
})

test_that(".validate_data accepts valid matrix data", {
  y <- matrix(c(1,0,1,0, 0,1,0,1, 1,1,0,0), nrow = 3, ncol = 4)
  result <- .validate_data(y, "matrix")

  expect_equal(result$N, 3)
  expect_equal(result$I, 4)
  expect_true(is.matrix(result$y))
  expect_equal(result$data_format, "matrix")
})

test_that(".validate_data rejects non-binary data", {
  y <- matrix(c(1,0,2,0), nrow = 2, ncol = 2)
  expect_error(.validate_data(y, "matrix"), "binary")
})

test_that(".validate_data rejects too-small data", {
  y <- matrix(c(1,0), nrow = 1, ncol = 2)
  expect_error(.validate_data(y, "matrix"), "at least 2 persons")
})

test_that(".detect_data_format works", {
  y_mat <- matrix(0, 10, 5)
  expect_equal(.detect_data_format(y_mat), "matrix")

  y_df <- data.frame(a = 1:10, b = 1:10, c = 1:10, d = 1:10, e = 1:10)
  expect_equal(.detect_data_format(y_df), "matrix")
})

test_that(".format_time produces readable output", {
  expect_match(.format_time(30), "sec")
  expect_match(.format_time(120), "min")
  expect_match(.format_time(7200), "hr")
})

test_that("digest_simple produces consistent hashes", {
  h1 <- digest_simple("test_string")
  h2 <- digest_simple("test_string")
  expect_equal(h1, h2)

  h3 <- digest_simple("different_string")
  expect_false(h1 == h3)
})
