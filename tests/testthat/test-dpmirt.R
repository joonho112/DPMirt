# ============================================================================
# Tests for R/dpmirt.R
# Main entry point validation and fit-object assembly contracts
# ============================================================================

test_that(".validate_draw_storage_args accepts current in-memory mode", {
  expect_true(DPMirt:::.validate_draw_storage_args(TRUE, NULL))
})


test_that(".validate_draw_storage_args rejects unsupported draw-storage modes", {
  expect_error(
    DPMirt:::.validate_draw_storage_args(FALSE, NULL),
    "save_draws.*reserved.*TRUE"
  )

  expect_error(
    DPMirt:::.validate_draw_storage_args(NA, NULL),
    "save_draws.*reserved.*TRUE"
  )

  expect_error(
    DPMirt:::.validate_draw_storage_args(TRUE, tempfile()),
    "save_path.*reserved.*disk-backed"
  )
})


test_that("dpmirt rejects unsupported draw-storage arguments before fitting", {
  response <- matrix(c(0, 1, 1, 0), nrow = 2)

  expect_error(
    dpmirt(response, save_draws = FALSE, verbose = FALSE),
    "save_draws.*reserved.*TRUE"
  )

  draw_path <- tempfile()
  expect_error(
    dpmirt(response, save_path = draw_path, verbose = FALSE),
    "save_path.*reserved.*disk-backed"
  )
  expect_false(file.exists(draw_path))
})


test_that("dpmirt rejects invalid MCMC controls before fitting", {
  response <- matrix(c(0, 1, 1, 0), nrow = 2)

  expect_error(
    dpmirt(response, niter = 0, verbose = FALSE),
    "niter.*positive"
  )
  expect_error(
    dpmirt(response, niter = 10, nburnin = 10, verbose = FALSE),
    "niter.*greater than nburnin"
  )
  expect_error(
    dpmirt(response, niter = 10, nburnin = -1, verbose = FALSE),
    "nburnin.*non-negative"
  )
  expect_error(
    dpmirt(response, thin = 0, verbose = FALSE),
    "thin.*positive"
  )
  expect_error(
    dpmirt(response, niter = 10, nburnin = 0, thin = 20,
           verbose = FALSE),
    "zero main draws"
  )
  expect_error(
    dpmirt(response, thin2 = 0, verbose = FALSE),
    "thin2.*positive"
  )
  expect_error(
    dpmirt(response, niter = 10, nburnin = 0, thin2 = 20,
           verbose = FALSE),
    "zero theta draws"
  )
  expect_error(
    dpmirt(response, nchains = 0, verbose = FALSE),
    "nchains.*positive"
  )
  expect_error(
    dpmirt(response, nchains = 1.5, verbose = FALSE),
    "nchains.*whole"
  )
})


test_that("dpmirt rejects reserved item priors and sampler schemas before fitting", {
  response <- matrix(c(0, 1, 1, 0), nrow = 2)

  expect_error(
    dpmirt(
      response,
      item_priors = list(beta = list(mean = 0, var = 3)),
      verbose = FALSE
    ),
    "item_priors.*reserved"
  )

  expect_error(
    dpmirt(
      response,
      sampler_config = list(type = "slice"),
      verbose = FALSE
    ),
    "sampler_config.*reserved"
  )
})
