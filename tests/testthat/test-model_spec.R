# ============================================================================
# Tests for R/model_spec.R — Model specification and code generation
# ============================================================================

# --- Helper: create small test data ---
make_test_data <- function(N = 20, I = 5, seed = 42) {
  set.seed(seed)
  theta <- rnorm(N)
  beta  <- seq(-1, 1, length.out = I)
  y <- matrix(NA, N, I)
  for (i in seq_len(I)) {
    prob <- 1 / (1 + exp(-(theta - beta[i])))
    y[, i] <- rbinom(N, 1, prob)
  }
  storage.mode(y) <- "double"
  y
}


# ============================================================================
# dpmirt_spec() — Rasch + Normal
# ============================================================================

test_that("dpmirt_spec creates valid spec for Rasch-Normal-constrainedItem", {
  y <- make_test_data()
  spec <- dpmirt_spec(y, model = "rasch", prior = "normal")

  expect_s3_class(spec, "dpmirt_spec")
  expect_equal(spec$config$model, "rasch")
  expect_equal(spec$config$prior, "normal")
  expect_equal(spec$config$identification, "constrained_item")
  expect_equal(spec$config$N, 20)
  expect_equal(spec$config$I, 5)

  # Check code is nimbleCode
  expect_true(inherits(spec$code, "nimbleCode") ||
              is.call(spec$code) || is.expression(spec$code))

  # Check monitors
  expect_true("beta" %in% spec$monitors ||
              "beta.tmp" %in% names(spec$inits))
  expect_true("eta" %in% spec$monitors2)

  # Check inits have required components
  expect_true(!is.null(spec$inits$eta))
  expect_true(!is.null(spec$inits[["beta.tmp"]]) ||
              !is.null(spec$inits$beta))
})


test_that("dpmirt_spec creates valid spec for Rasch-Normal-unconstrained", {
  y <- make_test_data()
  spec <- dpmirt_spec(y, model = "rasch", prior = "normal",
                      identification = "unconstrained")

  expect_equal(spec$config$identification, "unconstrained")
  expect_true("mu" %in% spec$monitors)
  expect_true("s2.eta" %in% spec$monitors)
})


test_that("dpmirt_spec creates valid spec for Rasch-Normal-constrainedAbility", {
  y <- make_test_data()
  spec <- dpmirt_spec(y, model = "rasch", prior = "normal",
                      identification = "constrained_ability")

  expect_equal(spec$config$identification, "constrained_ability")
  # constrained_ability should NOT have mu, s2.eta in monitors
  expect_false("mu" %in% spec$monitors)
  expect_false("s2.eta" %in% spec$monitors)
})


# ============================================================================
# dpmirt_spec() — Rasch + DPM
# ============================================================================

test_that("dpmirt_spec creates valid spec for Rasch-DPM-constrainedItem", {
  y <- make_test_data()
  spec <- dpmirt_spec(y, model = "rasch", prior = "dpm")

  expect_equal(spec$config$prior, "dpm")
  expect_equal(spec$config$identification, "constrained_item")

  # DPM-specific constants
  expect_true(!is.null(spec$constants$M))
  expect_true(!is.null(spec$constants$a))
  expect_true(!is.null(spec$constants$b))

  # DPM monitors
  expect_true("alpha" %in% spec$monitors)
  expect_true("zi" %in% spec$monitors)
  expect_true("muTilde" %in% spec$monitors)
  expect_true("s2Tilde" %in% spec$monitors)

  # DPM inits
  expect_true(!is.null(spec$inits$zi))
  expect_true(!is.null(spec$inits$alpha))
  expect_true(!is.null(spec$inits$muTilde))
  expect_true(!is.null(spec$inits$s2Tilde))
})


test_that("dpmirt_spec creates valid spec for Rasch-DPM-unconstrained", {
  y <- make_test_data()
  spec <- dpmirt_spec(y, model = "rasch", prior = "dpm",
                      identification = "unconstrained")

  expect_equal(spec$config$identification, "unconstrained")
  expect_true("alpha" %in% spec$monitors)
})


# ============================================================================
# Invalid combinations
# ============================================================================

test_that("dpmirt_spec rejects constrained_ability + dpm", {
  y <- make_test_data()
  expect_error(
    dpmirt_spec(y, model = "rasch", prior = "dpm",
                identification = "constrained_ability"),
    "defeats the purpose"
  )
})


test_that("dpmirt_spec rejects SI + Rasch", {
  y <- make_test_data()
  expect_error(
    dpmirt_spec(y, model = "rasch", prior = "normal",
                parameterization = "si"),
    "not meaningful"
  )
})


# ============================================================================
# Alpha prior resolution
# ============================================================================

test_that(".resolve_alpha_prior handles various inputs", {
  # NULL -> default
  result <- .resolve_alpha_prior(NULL, 100)
  expect_equal(result, c(a = 1, b = 3))

  # Numeric vector
  result <- .resolve_alpha_prior(c(2, 4), 100)
  expect_equal(unname(result), c(2, 4))

  # List with a, b (DPprior_fit-like)
  result <- .resolve_alpha_prior(list(a = 1.5, b = 2.5), 100)
  expect_equal(unname(result), c(1.5, 2.5))

  # Invalid
  expect_error(.resolve_alpha_prior("bad", 100))
  expect_error(.resolve_alpha_prior(c(-1, 3), 100))
})


# ============================================================================
# print method
# ============================================================================

test_that("print.dpmirt_spec works", {
  y <- make_test_data()
  spec <- dpmirt_spec(y, model = "rasch", prior = "normal")
  expect_output(print(spec), "DPMirt Model Specification")
})
