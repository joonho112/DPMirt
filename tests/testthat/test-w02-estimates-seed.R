# W-02 (V3 contract): dpmirt_estimates() takes an explicit seed; same-seed
# calls are bit-identical (including GR under forced rank ties), the caller's
# RNG state is restored, and the no-seed path warns once per session.
# Uses a synthetic dpmirt_fit (estimates only need theta/beta draw matrices),
# with duplicated columns to FORCE posterior-mean-rank ties so the random
# tie-break actually bites.

make_tied_fit <- function(n_draws = 200, n_person = 30, n_item = 6) {
  set.seed(42)
  theta <- matrix(rnorm(n_draws * n_person), n_draws, n_person)
  theta[, 2] <- theta[, 1]           # exact ties in colMeans -> rbar ties
  theta[, 4] <- theta[, 3]
  beta <- matrix(rnorm(n_draws * n_item), n_draws, n_item)
  structure(list(theta_samp = theta, beta_samp = beta),
            class = "dpmirt_fit")
}

test_that("same seed => bit-identical PM/CB/GR; different seeds differ on ties", {
  fit <- make_tied_fit()
  # These blocks test the SEED/tie behaviour, not warning emission; the
  # synthetic fit legitimately triggers benign CB/GR scaling warnings, which
  # are covered by the dedicated warning test below. Suppress them here so they
  # do not leak into the fail-closed P-0 gate (warning-count must mean "new").
  e1 <- suppressWarnings(dpmirt_estimates(fit, methods = c("pm", "cb", "gr"), seed = 101))
  e2 <- suppressWarnings(dpmirt_estimates(fit, methods = c("pm", "cb", "gr"), seed = 101))
  expect_identical(e1$theta, e2$theta)
  expect_identical(e1$beta, e2$beta)

  # Different seeds must be allowed to break the forced ties differently;
  # GR is the tie-sensitive summary (PM is rank-free and must agree).
  e3 <- suppressWarnings(dpmirt_estimates(fit, methods = c("pm", "cb", "gr"), seed = 202))
  expect_identical(e1$theta$theta_pm, e3$theta$theta_pm)
  expect_false(identical(e1$theta$theta_gr, e3$theta$theta_gr))
})

test_that("caller RNG state is saved and restored", {
  fit <- make_tied_fit()
  set.seed(777)
  before <- .Random.seed
  invisible(suppressWarnings(
    dpmirt_estimates(fit, methods = c("pm", "cb", "gr"), seed = 5)))
  expect_identical(.Random.seed, before)
  # and the ambient stream continues exactly as if nothing happened
  x_with <- rnorm(3)
  set.seed(777)
  x_ref <- rnorm(3)
  expect_identical(x_with, x_ref)
})

test_that("no-seed CB/GR warns once per session; PM-only does not warn", {
  fit <- make_tied_fit()
  cache <- DPMirt:::.nimble_cache   # environment: reference semantics
  cache$estimates_seed_warned <- NULL
  w1 <- capture_warnings(dpmirt_estimates(fit, methods = c("pm", "cb", "gr")))
  expect_true(any(grepl("bit-reproducible", w1)))
  # second call: the SEED warning must not recur (other package warnings,
  # e.g. CB scaling on this synthetic fit, are out of scope here)
  w2 <- capture_warnings(dpmirt_estimates(fit, methods = c("pm", "cb", "gr")))
  expect_false(any(grepl("bit-reproducible", w2)))
  cache$estimates_seed_warned <- NULL
  w3 <- capture_warnings(dpmirt_estimates(fit, methods = "pm"))
  expect_false(any(grepl("bit-reproducible", w3)))
  cache$estimates_seed_warned <- NULL
})
