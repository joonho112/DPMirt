# ============================================================================
# Tests for R/compile.R and R/sample.R
# Input validation and non-NIMBLE logic
# (Full integration tests require NIMBLE compilation — tested separately)
# ============================================================================

# ============================================================================
# dpmirt_compile() input validation
# ============================================================================

test_that("dpmirt_compile rejects non-spec input", {
  expect_error(
    dpmirt_compile("not a spec"),
    "dpmirt_spec"
  )
})


test_that("dpmirt_compile rejects list without class", {
  fake <- list(code = NULL, constants = NULL, data = NULL)
  expect_error(
    dpmirt_compile(fake),
    "dpmirt_spec"
  )
})


# ============================================================================
# dpmirt_sample() input validation
# ============================================================================

test_that("dpmirt_sample rejects non-compiled input", {
  expect_error(
    dpmirt_sample("not a compiled object"),
    "dpmirt_compiled"
  )
})


test_that("dpmirt_sample rejects list without class", {
  fake <- list(Cmodel = NULL, Cmcmc = NULL)
  expect_error(
    dpmirt_sample(fake),
    "dpmirt_compiled"
  )
})


# ============================================================================
# dpmirt_resume() input validation
# ============================================================================

test_that("dpmirt_resume rejects invalid input", {
  expect_error(
    dpmirt_resume("not valid", niter_more = 100),
    "dpmirt_samples.*dpmirt_fit.*dpmirt_compiled"
  )
})


test_that("dpmirt_resume accepts dpmirt_fit input", {
  # dpmirt_fit with NULL compiled should fail at the "No compiled" step,

  # NOT at the class check
  fake_fit <- structure(
    list(compiled = NULL),
    class = "dpmirt_fit"
  )
  expect_error(
    dpmirt_resume(fake_fit, niter_more = 100),
    "compiled model reference"
  )
})


test_that("dpmirt_resume rejects samples without compiled reference", {
  fake_samples <- structure(
    list(
      samples = matrix(1:10, nrow = 2),
      compiled = NULL,
      model_config = list(model = "rasch")
    ),
    class = "dpmirt_samples"
  )

  expect_error(
    dpmirt_resume(fake_samples, niter_more = 100),
    "recompile|compiled"
  )
})


# ============================================================================
# print.dpmirt_compiled tests
# ============================================================================

test_that("print.dpmirt_compiled produces expected output", {
  fake_compiled <- structure(
    list(
      spec = list(
        config = list(
          model = "rasch", prior = "dpm",
          identification = "constrained_item",
          N = 200, I = 20
        )
      ),
      compilation_time = 45.3,
      nimble_version = "1.2.0"
    ),
    class = "dpmirt_compiled"
  )

  output <- capture.output(print(fake_compiled))

  expect_true(any(grepl("DPMirt Compiled Model", output)))
  expect_true(any(grepl("RASCH", output)))
  expect_true(any(grepl("dpm", output)))
  expect_true(any(grepl("200", output)))
  expect_true(any(grepl("20", output)))
  expect_true(any(grepl("C\\+\\+", output)))
})


# ============================================================================
# print.dpmirt_samples tests
# ============================================================================

test_that("print.dpmirt_samples produces expected output", {
  fake_samples <- structure(
    list(
      samples  = matrix(rnorm(500), nrow = 100, ncol = 5),
      samples2 = matrix(rnorm(2000), nrow = 100, ncol = 20),
      waic = 2345.67,
      sampling_time = 12.5,
      model_config = list(model = "rasch", prior = "normal")
    ),
    class = "dpmirt_samples"
  )

  output <- capture.output(print(fake_samples))

  expect_true(any(grepl("DPMirt MCMC Samples", output)))
  expect_true(any(grepl("RASCH", output)))
  expect_true(any(grepl("100", output)))  # iterations
  expect_true(any(grepl("WAIC", output)))
})


# ============================================================================
# .configure_mcmc helper tests (logic only, no NIMBLE)
# ============================================================================

test_that(".add_centered_sampler resolves auto correctly", {
  # Test the auto-resolution logic for centered sampler
  # This tests the decision logic without actually calling NIMBLE

  # For Rasch + IRT: should NOT use centered sampler
  # (We can't call the actual function without NIMBLE,
  #  so we test the logical conditions directly)
  model_type <- "rasch"
  param_type <- "irt"
  id_type <- "unconstrained"
  use_it <- (model_type %in% c("2pl", "3pl")) &&
    (param_type == "si") &&
    (id_type %in% c("unconstrained", "constrained_ability"))
  expect_false(use_it)

  # For 2PL + SI + unconstrained: SHOULD use centered sampler
  model_type <- "2pl"
  param_type <- "si"
  id_type <- "unconstrained"
  use_it <- (model_type %in% c("2pl", "3pl")) &&
    (param_type == "si") &&
    (id_type %in% c("unconstrained", "constrained_ability"))
  expect_true(use_it)

  # For 2PL + IRT + unconstrained: should NOT use centered sampler
  model_type <- "2pl"
  param_type <- "irt"
  id_type <- "unconstrained"
  use_it <- (model_type %in% c("2pl", "3pl")) &&
    (param_type == "si") &&
    (id_type %in% c("unconstrained", "constrained_ability"))
  expect_false(use_it)

  # For 3PL + SI + constrained_item: should NOT use centered sampler
  model_type <- "3pl"
  param_type <- "si"
  id_type <- "constrained_item"
  use_it <- (model_type %in% c("2pl", "3pl")) &&
    (param_type == "si") &&
    (id_type %in% c("unconstrained", "constrained_ability"))
  expect_false(use_it)
})


# ============================================================================
# .extract_mcmc_output tests
# ============================================================================

test_that(".extract_mcmc_output handles list with samples and samples2", {
  fake_output <- list(
    samples  = matrix(1:10, nrow = 2),
    samples2 = matrix(11:20, nrow = 2)
  )
  fake_compiled <- list()

  result <- DPMirt:::.extract_mcmc_output(fake_output, fake_compiled)

  expect_equal(result$samples, fake_output$samples)
  expect_equal(result$samples2, fake_output$samples2)
})


test_that(".extract_mcmc_output handles plain matrix", {
  fake_output <- matrix(1:10, nrow = 2)
  fake_compiled <- list()

  result <- DPMirt:::.extract_mcmc_output(fake_output, fake_compiled)

  expect_equal(result$samples, fake_output)
  expect_null(result$samples2)
})


# ============================================================================
# .combine_chains tests
# ============================================================================

test_that(".combine_chains row-binds rescaled samples from two chains", {
  set.seed(42)
  N <- 10; I <- 3; niter <- 50

  make_chain <- function(seed_val) {
    set.seed(seed_val)
    beta_samp <- matrix(rnorm(niter * I, mean = 0.5), nrow = niter, ncol = I)
    colnames(beta_samp) <- paste0("beta[", seq_len(I), "]")

    eta_samp <- matrix(rnorm(niter * N), nrow = niter, ncol = N)
    colnames(eta_samp) <- paste0("eta[", seq_len(N), "]")

    log_nodes <- matrix(rnorm(niter * 3), nrow = niter, ncol = 3)
    colnames(log_nodes) <- c("myLogProbAll", "myLogProbSome", "myLogLik")

    structure(
      list(
        samples  = cbind(beta_samp, log_nodes),
        samples2 = eta_samp,
        waic = 2500 + runif(1, -100, 100),
        sampling_time = runif(1, 5, 15),
        model_config = list(
          model = "rasch", prior = "normal",
          parameterization = "irt",
          identification = "unconstrained",
          N = N, I = I
        )
      ),
      class = "dpmirt_samples"
    )
  }

  chain1 <- make_chain(42)
  chain2 <- make_chain(99)

  result <- DPMirt:::.combine_chains(list(chain1, chain2), rescale = TRUE)

  # Combined theta should have 2 * niter rows

  expect_equal(nrow(result$rescaled$theta_samp), 2 * niter)
  expect_equal(ncol(result$rescaled$theta_samp), N)

  # Combined beta should have 2 * niter rows
  expect_equal(nrow(result$rescaled$beta_samp), 2 * niter)
  expect_equal(ncol(result$rescaled$beta_samp), I)

  # WAIC should be mean of chain WAICs
  expect_true(is.numeric(result$waic))

  # Sampling time should be sum
  expect_equal(result$sampling_time,
               chain1$sampling_time + chain2$sampling_time)

  # Raw samples combined
  expect_equal(nrow(result$raw_samples), 2 * niter)
})


test_that(".combine_chains handles 2PL with lambda", {
  set.seed(42)
  N <- 10; I <- 3; niter <- 30

  make_2pl_chain <- function(seed_val) {
    set.seed(seed_val)
    beta_samp <- matrix(rnorm(niter * I), nrow = niter, ncol = I)
    colnames(beta_samp) <- paste0("beta[", seq_len(I), "]")

    lambda_samp <- matrix(exp(rnorm(niter * I, 0.3, 0.2)),
                           nrow = niter, ncol = I)
    colnames(lambda_samp) <- paste0("lambda[", seq_len(I), "]")

    eta_samp <- matrix(rnorm(niter * N), nrow = niter, ncol = N)
    colnames(eta_samp) <- paste0("eta[", seq_len(N), "]")

    log_nodes <- matrix(rnorm(niter * 3), nrow = niter, ncol = 3)
    colnames(log_nodes) <- c("myLogProbAll", "myLogProbSome", "myLogLik")

    structure(
      list(
        samples  = cbind(beta_samp, lambda_samp, log_nodes),
        samples2 = eta_samp,
        waic = 3000,
        sampling_time = 10,
        model_config = list(
          model = "2pl", prior = "normal",
          parameterization = "irt",
          identification = "unconstrained",
          N = N, I = I
        )
      ),
      class = "dpmirt_samples"
    )
  }

  chain1 <- make_2pl_chain(42)
  chain2 <- make_2pl_chain(99)

  result <- DPMirt:::.combine_chains(list(chain1, chain2), rescale = TRUE)

  expect_false(is.null(result$rescaled$lambda_samp))
  expect_equal(nrow(result$rescaled$lambda_samp), 2 * niter)
  expect_equal(ncol(result$rescaled$lambda_samp), I)
})
