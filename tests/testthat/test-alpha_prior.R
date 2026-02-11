# ============================================================================
# Tests for R/alpha_prior.R
# dpmirt_alpha_prior() and alpha prior resolution
# ============================================================================

# ============================================================================
# dpmirt_alpha_prior() tests
# ============================================================================

test_that("dpmirt_alpha_prior returns named numeric vector", {
  # This tests the graceful fallback when DPprior is NOT installed,

  # or the actual output when it IS installed.
  result <- suppressMessages(dpmirt_alpha_prior(N = 200))

  expect_true(is.numeric(result))
  expect_equal(length(result), 2)
  expect_true("a" %in% names(result))
  expect_true("b" %in% names(result))
  expect_true(result["a"] > 0)
  expect_true(result["b"] > 0)
})


test_that("dpmirt_alpha_prior returns Gamma(1,3) when DPprior not available", {
  # Temporarily hide DPprior by using requireNamespace behavior
  # If DPprior is not installed, we expect the default
  has_dpprior <- requireNamespace("DPprior", quietly = TRUE)

  if (!has_dpprior) {
    result <- suppressMessages(dpmirt_alpha_prior(N = 100))
    expect_equal(unname(result), c(1, 3))
  } else {
    # If DPprior IS installed, just verify it returns valid values
    result <- suppressMessages(dpmirt_alpha_prior(N = 100))
    expect_true(result["a"] > 0)
    expect_true(result["b"] > 0)
  }
})


test_that("dpmirt_alpha_prior produces message about default or result", {
  expect_message(
    dpmirt_alpha_prior(N = 200),
    "."  # At least some message is produced
  )
})


test_that("dpmirt_alpha_prior with custom mu_K produces valid output", {
  result <- suppressMessages(
    dpmirt_alpha_prior(N = 500, mu_K = 5, confidence = "medium")
  )

  expect_true(is.numeric(result))
  expect_equal(length(result), 2)
  expect_true(result["a"] > 0)
  expect_true(result["b"] > 0)
})


test_that("dpmirt_alpha_prior auto-selects mu_K when NULL", {
  # mu_K should default to max(3, ceiling(log(N)))
  # Use N = 200 to stay within DPprior's supported range
  expect_message(
    result <- dpmirt_alpha_prior(N = 200),
    "."
  )

  expect_true(is.numeric(result))
})


# ============================================================================
# .resolve_alpha_prior() tests (from model_spec.R / compile.R)
# ============================================================================

test_that(".resolve_alpha_prior handles NULL (default Gamma(1,3))", {
  result <- DPMirt:::.resolve_alpha_prior(NULL, N = 200)

  expect_true(is.numeric(result))
  expect_equal(length(result), 2)
  expect_equal(unname(result), c(1, 3))
})


test_that(".resolve_alpha_prior handles numeric vector input", {
  result <- DPMirt:::.resolve_alpha_prior(c(2, 5), N = 200)

  expect_equal(unname(result), c(2, 5))
})


test_that(".resolve_alpha_prior handles named vector input", {
  result <- DPMirt:::.resolve_alpha_prior(c(a = 0.5, b = 1.0), N = 200)

  expect_equal(unname(result), c(0.5, 1.0))
})
