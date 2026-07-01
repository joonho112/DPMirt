# ============================================================================
# Tests for R/alpha_prior.R
# dpmirt_alpha_prior() and alpha prior resolution
# ============================================================================

# ============================================================================
# dpmirt_alpha_prior() tests
# ============================================================================

.quiet_alpha_prior <- function(...) {
  # Use only in value-shape tests; condition behavior is asserted separately.
  suppressWarnings(suppressMessages(dpmirt_alpha_prior(...)))
}


.capture_alpha_prior_messages <- function(...) {
  result <- NULL
  messages <- capture.output(
    result <- dpmirt_alpha_prior(...),
    type = "message"
  )
  list(result = result, messages = messages)
}


test_that("dpmirt_alpha_prior returns named numeric vector", {
  # This tests the graceful fallback when DPprior is NOT installed,

  # or the actual output when it IS installed.
  result <- .quiet_alpha_prior(N = 200)

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
    result <- .quiet_alpha_prior(N = 100)
    expect_equal(unname(result), c(1, 3))
  } else {
    # If DPprior IS installed, just verify it returns valid values
    result <- .quiet_alpha_prior(N = 100)
    expect_true(result["a"] > 0)
    expect_true(result["b"] > 0)
  }
})


test_that("dpmirt_alpha_prior produces message about default or result", {
  captured <- .capture_alpha_prior_messages(N = 200)

  expect_true(
    any(grepl("Using default mu_K|Alpha prior|DPprior not installed",
              captured$messages))
  )
})


test_that("dpmirt_alpha_prior suppresses dominance warnings by default", {
  skip_if_not_installed("DPprior")
  skip_if_not("warn_dominance" %in% names(formals(DPprior::DPprior_fit)))

  expect_warning(
    result <- suppressMessages(
      dpmirt_alpha_prior(N = 200, mu_K = 5, confidence = "medium")
    ),
    NA
  )

  expect_true(is.numeric(result))
  expect_equal(length(result), 2)
})


test_that("dpmirt_alpha_prior can opt in to DPprior dominance warnings", {
  skip_if_not_installed("DPprior")
  skip_if_not("warn_dominance" %in% names(formals(DPprior::DPprior_fit)))

  expect_warning(
    suppressMessages(
      dpmirt_alpha_prior(
        N = 200, mu_K = 5, confidence = "medium",
        warn_dominance = TRUE
      )
    ),
    "HIGH DOMINANCE RISK: P\\(w1 > 0\\.5\\)"
  )
})


test_that("dpmirt_alpha_prior with custom mu_K produces valid output", {
  result <- .quiet_alpha_prior(N = 500, mu_K = 5, confidence = "medium")

  expect_true(is.numeric(result))
  expect_equal(length(result), 2)
  expect_true(result["a"] > 0)
  expect_true(result["b"] > 0)
})


test_that("dpmirt_alpha_prior auto-selects mu_K when NULL", {
  # mu_K should default to max(3, ceiling(log(N)))
  # Use N = 200 to stay within DPprior's supported range
  captured <- .capture_alpha_prior_messages(N = 200)
  result <- captured$result

  if (requireNamespace("DPprior", quietly = TRUE)) {
    expect_true(any(grepl("Using default mu_K = 6", captured$messages)))
  } else {
    expect_true(any(grepl("DPprior not installed", captured$messages)))
  }

  expect_true(
    any(grepl("Alpha prior|DPprior not installed", captured$messages))
  )

  expect_true(is.numeric(result))
})


test_that("dpmirt_alpha_prior return_fit returns DPprior_fit object", {
  skip_if_not_installed("DPprior")

  fit <- .quiet_alpha_prior(N = 200, mu_K = 5, confidence = "medium",
                            return_fit = TRUE)

  # Should be a DPprior_fit object
  expect_true(inherits(fit, "DPprior_fit"))
  expect_true(!is.null(fit$a))
  expect_true(!is.null(fit$b))
  expect_true(fit$a > 0)
  expect_true(fit$b > 0)

  # Should have target info
  expect_true(!is.null(fit$target))
  expect_equal(fit$J, 200)
})


test_that("dpmirt_alpha_prior return_fit=FALSE is backward compatible", {
  skip_if_not_installed("DPprior")

  result <- .quiet_alpha_prior(N = 200, mu_K = 5, confidence = "medium",
                               return_fit = FALSE)

  # Default behavior: named numeric vector
  expect_true(is.numeric(result))
  expect_false(inherits(result, "DPprior_fit"))
  expect_equal(length(result), 2)
  expect_equal(names(result), c("a", "b"))
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


test_that(".resolve_alpha_prior handles DPprior_fit object with class", {
  # Simulate a DPprior_fit object with the proper class
  fake_fit <- list(a = 1.5, b = 2.5, J = 200, method = "A2-MN")
  class(fake_fit) <- "DPprior_fit"

  result <- DPMirt:::.resolve_alpha_prior(fake_fit, N = 200)

  expect_equal(unname(result), c(1.5, 2.5))
  expect_equal(names(result), c("a", "b"))
})


test_that(".resolve_alpha_prior handles duck-typed list with a and b", {
  # A plain list with $a and $b (backward compatibility)
  fake_fit <- list(a = 0.8, b = 1.2)

  result <- DPMirt:::.resolve_alpha_prior(fake_fit, N = 200)

  expect_equal(unname(result), c(0.8, 1.2))
})


test_that(".resolve_alpha_prior rejects invalid DPprior_fit with non-positive params", {
  fake_fit <- list(a = -1, b = 2)
  class(fake_fit) <- "DPprior_fit"

  expect_error(
    DPMirt:::.resolve_alpha_prior(fake_fit, N = 200),
    "non-positive"
  )
})


test_that(".resolve_alpha_prior rejects invalid input types", {
  expect_error(
    DPMirt:::.resolve_alpha_prior("invalid", N = 200),
    "Invalid alpha_prior"
  )

  expect_error(
    DPMirt:::.resolve_alpha_prior(c(1, 2, 3), N = 200),
    "Invalid alpha_prior"
  )
})
