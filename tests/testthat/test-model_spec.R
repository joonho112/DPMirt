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


test_that("dpmirt_spec accepts only reserved empty item_priors", {
  y <- make_test_data()

  spec_null <- dpmirt_spec(y, model = "rasch", prior = "normal",
                           item_priors = NULL)
  expect_identical(spec_null$config$item_priors, list())

  spec_empty <- dpmirt_spec(y, model = "rasch", prior = "normal",
                            item_priors = list())
  expect_identical(spec_empty$config$item_priors, list())

  expect_error(
    dpmirt_spec(y, model = "rasch", prior = "normal",
                item_priors = "default"),
    "item_priors.*reserved"
  )

  expect_error(
    dpmirt_spec(
      y, model = "rasch", prior = "normal",
      item_priors = list(beta = list(mean = 0, var = 3))
    ),
    "Non-empty item_priors.*reserved"
  )
})


test_that("chain-specific initial values are reproducible and isolated", {
  y <- make_test_data(N = 16, I = 4, seed = 20260211)
  specs <- list(
    rasch_normal = dpmirt_spec(
      y, model = "rasch", prior = "normal", identification = "unconstrained"
    ),
    rasch_dpm = dpmirt_spec(
      y, model = "rasch", prior = "dpm", M = 8L
    ),
    two_pl = dpmirt_spec(
      y, model = "2pl", prior = "normal", identification = "unconstrained"
    )
  )

  state_before <- .Random.seed
  for (spec in specs) {
    init1 <- DPMirt:::.generate_chain_inits(spec, chain_id = 1L, seed = 700L)
    init1_again <- DPMirt:::.generate_chain_inits(
      spec, chain_id = 1L, seed = 700L
    )
    init2 <- DPMirt:::.generate_chain_inits(spec, chain_id = 2L, seed = 700L)

    expect_equal(init1, init1_again)
    expect_equal(names(init1), names(init2))
    expect_equal(length(init1$eta), spec$config$N)
    expect_false(identical(
      unlist(init1, use.names = FALSE),
      unlist(init2, use.names = FALSE)
    ))
  }
  expect_identical(.Random.seed, state_before)
})


test_that(".chain_initial_values records per-chain init provenance", {
  y <- make_test_data(N = 12, I = 4, seed = 20260212)
  spec <- dpmirt_spec(y, model = "rasch", prior = "normal")

  seeded_1 <- DPMirt:::.chain_initial_values(
    spec, chain_id = 1L, seed = 11L, nchains = 2L
  )
  seeded_2 <- DPMirt:::.chain_initial_values(
    spec, chain_id = 2L, seed = 11L, nchains = 2L
  )
  expect_equal(seeded_1$init_seed, 11L)
  expect_equal(seeded_2$init_seed, 12L)
  expect_equal(seeded_1$init_strategy, "chain_seeded")
  expect_false(identical(seeded_1$inits, seeded_2$inits))

  random_multi <- DPMirt:::.chain_initial_values(
    spec, chain_id = 1L, seed = NULL, nchains = 2L
  )
  expect_null(random_multi$init_seed)
  expect_equal(random_multi$init_strategy, "chain_random")
  expect_true(is.list(random_multi$inits))

  inherited_single <- DPMirt:::.chain_initial_values(
    spec, chain_id = 1L, seed = NULL, nchains = 1L
  )
  expect_null(inherited_single$inits)
  expect_null(inherited_single$init_seed)
  expect_equal(inherited_single$init_strategy, "compiled_spec")
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


test_that("dpmirt_spec rejects constrained_item + 3PL", {
  y <- make_test_data()
  expect_error(
    dpmirt_spec(y, model = "3pl", prior = "normal",
                identification = "constrained_item"),
    "not implemented"
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
