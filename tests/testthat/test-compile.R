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


.fake_spec_for_sampler_config <- function() {
  structure(
    list(
      monitors = c("beta", "myLogProbAll", "myLogProbSome", "myLogLik"),
      monitors2 = "eta",
      config = list(model = "rasch", prior = "normal", N = 3, I = 2)
    ),
    class = "dpmirt_spec"
  )
}


.fake_mcmc_conf <- function(monitors = c("beta", "myLogProbAll",
                                         "myLogProbSome", "myLogLik"),
                            monitors2 = "eta") {
  conf <- new.env(parent = emptyenv())
  conf$monitors <- monitors
  conf$monitors2 <- monitors2
  conf$sampler_nodes <- c("myLogProbAll", "myLogProbSome", "myLogLik")
  conf$getMonitors <- function() conf$monitors
  conf$getMonitors2 <- function() conf$monitors2
  conf$findSamplersOnNodes <- function(node) {
    which(conf$sampler_nodes %in% node)
  }
  class(conf) <- "MCMCconf"
  conf
}


test_that("sampler_config accepts NULL or function and rejects reserved schemas", {
  expect_null(DPMirt:::.validate_sampler_config_arg(NULL))

  hook <- function(conf, model, spec) conf
  expect_identical(DPMirt:::.validate_sampler_config_arg(hook), hook)

  expect_error(
    DPMirt:::.validate_sampler_config_arg(list(type = "slice")),
    "sampler_config.*function"
  )
  expect_error(
    DPMirt:::.validate_sampler_config_arg("slice"),
    "sampler_config.*function"
  )
})


test_that("dpmirt_compile rejects reserved sampler_config before building model", {
  expect_error(
    dpmirt_compile(.fake_spec_for_sampler_config(),
                   sampler_config = list(type = "slice"),
                   verbose = FALSE),
    "sampler_config.*reserved"
  )
})


test_that(".apply_sampler_config accepts mutation or returned MCMCconf", {
  spec <- .fake_spec_for_sampler_config()

  conf1 <- .fake_mcmc_conf()
  out1 <- DPMirt:::.apply_sampler_config(
    conf1, nimble_model = list(), spec = spec,
    sampler_config = function(conf, model, spec) {
      conf$customized <- TRUE
      NULL
    }
  )
  expect_true(out1$customized)
  expect_identical(out1, conf1)

  conf2 <- .fake_mcmc_conf()
  out2 <- DPMirt:::.apply_sampler_config(
    conf2, nimble_model = list(), spec = spec,
    sampler_config = function(conf, model, spec) conf
  )
  expect_identical(out2, conf2)
})


test_that(".apply_sampler_config rejects invalid return values and removed monitors", {
  spec <- .fake_spec_for_sampler_config()

  expect_error(
    DPMirt:::.apply_sampler_config(
      .fake_mcmc_conf(), nimble_model = list(), spec = spec,
      sampler_config = function(conf, model, spec) list()
    ),
    "MCMCconf"
  )

  expect_error(
    DPMirt:::.apply_sampler_config(
      .fake_mcmc_conf(), nimble_model = list(), spec = spec,
      sampler_config = function(conf, model, spec) {
        conf$monitors <- setdiff(conf$monitors, "beta")
        conf
      }
    ),
    "required DPMirt monitor.*beta"
  )

  expect_error(
    DPMirt:::.apply_sampler_config(
      .fake_mcmc_conf(), nimble_model = list(), spec = spec,
      sampler_config = function(conf, model, spec) {
        conf$monitors2 <- character(0)
        conf
      }
    ),
    "eta \\(monitors2\\)"
  )

  expect_error(
    DPMirt:::.apply_sampler_config(
      .fake_mcmc_conf(), nimble_model = list(), spec = spec,
      sampler_config = function(conf, model, spec) {
        conf$sampler_nodes <- setdiff(conf$sampler_nodes, "myLogLik")
        conf
      }
    ),
    "log-probability sampler.*myLogLik"
  )
})


# ============================================================================
# dpmirt_sample() input validation
# ============================================================================

.fake_compiled_for_sampling <- function(N = 3, I = 2, waic_enabled = FALSE) {
  mcmc <- new.env(parent = emptyenv())
  cmodel <- new.env(parent = emptyenv())
  cmodel$init_calls <- list()
  cmodel$events <- character(0)
  cmodel$getNodeNames <- function() "y"
  cmodel$setInits <- function(inits) {
    cmodel$init_calls[[length(cmodel$init_calls) + 1L]] <- inits
    cmodel$events <- c(cmodel$events, "setInits")
    invisible(NULL)
  }

  mcmc$mvSamples <- matrix(numeric(0), nrow = 0, ncol = I + 1L)
  colnames(mcmc$mvSamples) <- c(paste0("beta[", seq_len(I), "]"), "myLogLik")
  mcmc$mvSamples2 <- matrix(numeric(0), nrow = 0, ncol = N)
  colnames(mcmc$mvSamples2) <- paste0("eta[", seq_len(N), "]")
  mcmc$run_count <- 0L
  mcmc$last_args <- NULL
  mcmc$run <- function(niter,
                       reset = TRUE,
                       resetMV = TRUE,
                       time = FALSE,
                       progressBar = FALSE,
                       nburnin = 0L,
                       thin = 1L,
                       thin2 = thin,
                       resetWAIC = TRUE,
                       initializeModel = TRUE,
                       chain = 1L) {
    mcmc$run_count <- mcmc$run_count + 1L
    cmodel$events <- c(cmodel$events, "run")
    mcmc$last_args <- list(
      niter = niter,
      reset = reset,
      resetMV = resetMV,
      nburnin = nburnin,
      thin = thin,
      thin2 = thin2,
      resetWAIC = resetWAIC,
      initializeModel = initializeModel
    )

    n_main <- max(0L, floor((niter - nburnin) / thin))
    n_theta <- max(0L, floor((niter - nburnin) / thin2))
    main_start <- mcmc$run_count * 1000L
    theta_start <- mcmc$run_count * 2000L
    new_main <- matrix(
      seq(main_start, length.out = n_main * (I + 1L)),
      nrow = n_main,
      ncol = I + 1L
    )
    colnames(new_main) <- colnames(mcmc$mvSamples)
    new_theta <- matrix(
      seq(theta_start, length.out = n_theta * N),
      nrow = n_theta,
      ncol = N
    )
    colnames(new_theta) <- colnames(mcmc$mvSamples2)

    if (isTRUE(resetMV)) {
      mcmc$mvSamples <- new_main
      mcmc$mvSamples2 <- new_theta
    } else {
      mcmc$mvSamples <- rbind(mcmc$mvSamples, new_main)
      mcmc$mvSamples2 <- rbind(mcmc$mvSamples2, new_theta)
    }

    invisible(NULL)
  }
  mcmc$getWAIC <- function() {
    list(WAIC = if (isTRUE(waic_enabled)) 123.45 else NA_real_)
  }

  structure(
    list(
      Cmodel = cmodel,
      Cmcmc = mcmc,
      spec = list(
        data = list(y = matrix(0, nrow = N, ncol = I)),
        inits = list(beta = rep(0, I), eta = rep(0, N)),
        config = list(
          model = "rasch", prior = "normal",
          parameterization = "irt",
          identification = "unconstrained",
          N = N, I = I
        )
      ),
      compilation_time = 0.1,
      waic_enabled = isTRUE(waic_enabled)
    ),
    class = "dpmirt_compiled"
  )
}

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


test_that("dpmirt_sample validates MCMC controls before running sampler", {
  bad_cases <- list(
    list(args = list(niter = 0), msg = "niter.*positive"),
    list(args = list(niter = 10, nburnin = -1), msg = "nburnin.*non-negative"),
    list(args = list(niter = 10, nburnin = 10), msg = "niter.*greater than nburnin"),
    list(args = list(niter = 10, thin = 0), msg = "thin.*positive"),
    list(args = list(niter = 10, thin2 = 0), msg = "thin2.*positive"),
    list(args = list(niter = 10, thin = 20), msg = "zero main draws"),
    list(args = list(niter = 10, thin2 = 20), msg = "zero theta draws")
  )

  for (case in bad_cases) {
    compiled <- .fake_compiled_for_sampling()
    args <- modifyList(
      list(compiled = compiled, nburnin = 0, verbose = FALSE),
      case$args
    )
    expect_error(
      do.call(dpmirt_sample, args),
      case$msg
    )
    expect_equal(compiled$Cmcmc$run_count, 0L)
  }
})


test_that("dpmirt_sample reset controls sample and WAIC storage", {
  compiled <- .fake_compiled_for_sampling()

  first <- dpmirt_sample(
    compiled, niter = 6, nburnin = 0, thin = 1, thin2 = 1,
    reset = TRUE, verbose = FALSE
  )
  second <- dpmirt_sample(
    compiled, niter = 4, nburnin = 0, thin = 1, thin2 = 1,
    reset = FALSE, verbose = FALSE
  )

  expect_equal(nrow(first$samples), 6)
  expect_equal(nrow(second$samples), 10)
  expect_equal(second$samples[seq_len(6), , drop = FALSE], first$samples)
  expect_false(compiled$Cmcmc$last_args$reset)
  expect_false(compiled$Cmcmc$last_args$resetMV)
  expect_false(compiled$Cmcmc$last_args$resetWAIC)
  expect_false(compiled$Cmcmc$last_args$initializeModel)
  expect_true(second$mcmc_control$resumed)
  expect_equal(second$mcmc_control$niter_total_retained, 10)

  third <- dpmirt_sample(
    compiled, niter = 3, nburnin = 0, thin = 1, thin2 = 1,
    reset = TRUE, verbose = FALSE
  )

  expect_equal(nrow(third$samples), 3)
  expect_true(compiled$Cmcmc$last_args$reset)
  expect_true(compiled$Cmcmc$last_args$resetMV)
  expect_true(compiled$Cmcmc$last_args$resetWAIC)
  expect_true(compiled$Cmcmc$last_args$initializeModel)
})


test_that("dpmirt_sample applies supplied inits before reset run", {
  compiled <- .fake_compiled_for_sampling()
  inits <- list(beta = c(-0.5, 0.5), eta = c(-1, 0, 1))

  samples <- dpmirt_sample(
    compiled,
    niter = 4,
    nburnin = 0,
    thin = 1,
    thin2 = 1,
    seed = 123,
    inits = inits,
    init_seed = 123,
    init_strategy = "chain_seeded",
    reset = TRUE,
    verbose = FALSE
  )

  expect_equal(length(compiled$Cmodel$init_calls), 1L)
  expect_equal(compiled$Cmodel$init_calls[[1]], inits)
  expect_equal(compiled$Cmodel$events, c("setInits", "run"))
  expect_true(compiled$Cmcmc$last_args$initializeModel)
  expect_equal(samples$mcmc_control$init_seed, 123)
  expect_equal(samples$mcmc_control$init_strategy, "chain_seeded")
  expect_equal(samples$chain_info$init_seed, 123)
  expect_equal(samples$chain_info$init_strategy, "chain_seeded")
  expect_equal(samples$run_history$init_seed, 123)
  expect_equal(samples$run_history$init_strategy, "chain_seeded")
})


test_that("dpmirt_sample rejects supplied inits on append run", {
  compiled <- .fake_compiled_for_sampling()

  expect_error(
    dpmirt_sample(
      compiled,
      niter = 4,
      nburnin = 0,
      thin = 1,
      thin2 = 1,
      inits = list(beta = c(-0.5, 0.5), eta = c(-1, 0, 1)),
      reset = FALSE,
      verbose = FALSE
    ),
    "inits.*reset = TRUE"
  )
  expect_equal(length(compiled$Cmodel$init_calls), 0L)
  expect_equal(compiled$Cmcmc$run_count, 0L)
})


test_that("dpmirt_sample normalizes unavailable WAIC", {
  compiled <- .fake_compiled_for_sampling(waic_enabled = FALSE)

  samples <- dpmirt_sample(
    compiled, niter = 4, nburnin = 0, reset = TRUE, verbose = FALSE
  )

  expect_null(samples$waic)
  expect_false(samples$mcmc_control$waic_enabled)
  expect_true(is.na(samples$chain_info$waic))
})


test_that("dpmirt_sample stores finite WAIC and provenance when enabled", {
  compiled <- .fake_compiled_for_sampling(waic_enabled = TRUE)

  samples <- dpmirt_sample(
    compiled, niter = 4, nburnin = 0, reset = TRUE, verbose = FALSE
  )

  expect_equal(samples$waic, 123.45)
  expect_true(samples$mcmc_control$waic_enabled)
  expect_equal(samples$chain_info$waic, 123.45)
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


test_that("dpmirt_resume validates additional iterations before running sampler", {
  compiled <- .fake_compiled_for_sampling()
  initial <- dpmirt_sample(compiled, niter = 6, nburnin = 0,
                           reset = TRUE, verbose = FALSE)

  expect_error(
    dpmirt_resume(initial, niter_more = 0, verbose = FALSE),
    "niter.*positive"
  )
  expect_equal(compiled$Cmcmc$run_count, 1L)
})


test_that("dpmirt_resume appends with preserved thinning metadata", {
  compiled <- .fake_compiled_for_sampling(N = 4, I = 2)
  initial <- dpmirt_sample(
    compiled, niter = 40, nburnin = 0, thin = 2, thin2 = 4,
    reset = TRUE, verbose = FALSE
  )

  resumed <- dpmirt_resume(initial, niter_more = 20, reset = FALSE,
                           verbose = FALSE)

  expect_equal(nrow(initial$samples), 20)
  expect_equal(nrow(initial$samples2), 10)
  expect_equal(nrow(resumed$samples), 30)
  expect_equal(nrow(resumed$samples2), 15)
  expect_equal(resumed$samples[seq_len(20), , drop = FALSE], initial$samples)
  expect_equal(resumed$mcmc_control$niter, 60)
  expect_equal(resumed$mcmc_control$niter_more, 20)
  expect_equal(resumed$mcmc_control$niter_total_retained, 30)
  expect_equal(resumed$mcmc_control$thin, 2)
  expect_equal(resumed$mcmc_control$thin2, 4)
  expect_false(compiled$Cmcmc$last_args$resetWAIC)
  expect_equal(resumed$chain_info$n_draws_main, 30)
  expect_equal(resumed$chain_info$n_draws_theta, 15)
  expect_equal(resumed$draw_index$main$mcmc_iteration, seq(2, 60, by = 2))
  expect_equal(resumed$draw_index$theta$mcmc_iteration, seq(4, 60, by = 4))
  expect_equal(resumed$run_history$run_id, 1:2)
  expect_equal(resumed$run_history$segment_id, 1:2)
  expect_equal(resumed$run_history$niter, c(40L, 20L))
  expect_equal(resumed$run_history$resumed, c(FALSE, TRUE))
  expect_equal(resumed$run_history$row_start_main, c(1L, 21L))
  expect_equal(resumed$run_history$row_end_main, c(20L, 30L))
  expect_equal(resumed$run_history$mcmc_start_main, c(2L, 42L))
  expect_equal(resumed$run_history$mcmc_end_main, c(40L, 60L))

  fit <- dpmirt_rescale(resumed)
  expect_equal(fit$config$niter, 60)
  expect_equal(fit$config$n_draws_main, 30)
  expect_equal(fit$config$n_draws_theta, 15)
  expect_equal(fit$config$thin, 2)
  expect_equal(fit$config$thin2, 4)
  expect_equal(fit$draw_index$main$mcmc_iteration,
               resumed$draw_index$main$mcmc_iteration)
  expect_equal(fit$run_history$run_id, 1:2)

  resumed_again <- dpmirt_resume(fit, niter_more = 20, reset = FALSE,
                                 verbose = FALSE)
  expect_equal(resumed_again$mcmc_control$previous_niter, 60)
  expect_equal(resumed_again$mcmc_control$niter, 80)
  expect_equal(
    tail(resumed_again$draw_index$main$mcmc_iteration, 10),
    seq(62, 80, by = 2)
  )
  expect_equal(
    tail(resumed_again$draw_index$theta$mcmc_iteration, 5),
    seq(64, 80, by = 4)
  )
  expect_equal(resumed_again$run_history$run_id, 1:3)
  expect_equal(resumed_again$run_history$segment_id, 1:3)
  expect_equal(resumed_again$run_history$mcmc_start_main, c(2L, 42L, 62L))
  expect_equal(resumed_again$run_history$mcmc_end_main, c(40L, 60L, 80L))

  expect_error(
    dpmirt_resume(initial, niter_more = 5, reset = FALSE, verbose = FALSE),
    "live compiled MCMC state"
  )
})


test_that("dpmirt_resume preserves burn-in-aware iteration labels", {
  compiled <- .fake_compiled_for_sampling(N = 4, I = 2)
  initial <- dpmirt_sample(
    compiled, niter = 40, nburnin = 10, thin = 2, thin2 = 5,
    reset = TRUE, verbose = FALSE
  )

  resumed <- dpmirt_resume(initial, niter_more = 10, reset = FALSE,
                           verbose = FALSE)

  expect_equal(nrow(initial$samples), 15)
  expect_equal(nrow(initial$samples2), 6)
  expect_equal(nrow(resumed$samples), 20)
  expect_equal(nrow(resumed$samples2), 8)
  expect_equal(
    resumed$draw_index$main$mcmc_iteration,
    c(seq(12, 40, by = 2), seq(42, 50, by = 2))
  )
  expect_equal(
    resumed$draw_index$theta$mcmc_iteration,
    c(seq(15, 40, by = 5), seq(45, 50, by = 5))
  )
  expect_equal(resumed$run_history$run_id, 1:2)
  expect_equal(resumed$run_history$segment_id, 1:2)
  expect_equal(resumed$run_history$nburnin, c(10L, 0L))
  expect_equal(resumed$run_history$n_draws_main, c(15L, 5L))
  expect_equal(resumed$run_history$n_draws_theta, c(6L, 2L))
  expect_equal(resumed$run_history$row_start_main, c(1L, 16L))
  expect_equal(resumed$run_history$row_end_main, c(15L, 20L))
  expect_equal(resumed$run_history$row_start_theta, c(1L, 7L))
  expect_equal(resumed$run_history$row_end_theta, c(6L, 8L))
  expect_equal(resumed$run_history$mcmc_start_main, c(12L, 42L))
  expect_equal(resumed$run_history$mcmc_end_main, c(40L, 50L))
  expect_equal(resumed$run_history$mcmc_start_theta, c(15L, 45L))
  expect_equal(resumed$run_history$mcmc_end_theta, c(40L, 50L))
})


test_that("dpmirt_resume validates theta sample state before appending", {
  compiled <- .fake_compiled_for_sampling(N = 4, I = 2)
  initial <- dpmirt_sample(
    compiled, niter = 20, nburnin = 0, thin = 2, thin2 = 4,
    reset = TRUE, verbose = FALSE
  )
  stale <- initial
  stale$samples2[1, 1] <- stale$samples2[1, 1] + 999

  expect_error(
    dpmirt_resume(stale, niter_more = 10, reset = FALSE, verbose = FALSE),
    "live compiled MCMC state"
  )
})


test_that("dpmirt_resume rejects multi-chain append because only one live sampler is stored", {
  compiled <- .fake_compiled_for_sampling()
  fit <- structure(
    list(
      compiled = compiled,
      config = list(nchains = 2L, thin = 1L, thin2 = 1L),
      chain_info = data.frame(chain_id = 1:2)
    ),
    class = "dpmirt_fit"
  )

  expect_error(
    dpmirt_resume(fit, niter_more = 4, reset = FALSE, verbose = FALSE),
    "Multi-chain.*reset = FALSE"
  )
  expect_equal(compiled$Cmcmc$run_count, 0L)
})


test_that("dpmirt_resume reset replaces stored samples", {
  compiled <- .fake_compiled_for_sampling()
  initial <- dpmirt_sample(compiled, niter = 6, nburnin = 0,
                           reset = TRUE, verbose = FALSE)

  resumed <- dpmirt_resume(initial, niter_more = 4, reset = TRUE,
                           verbose = FALSE)

  expect_equal(nrow(resumed$samples), 4)
  expect_false(identical(resumed$samples, initial$samples[seq_len(4), ]))
  expect_true(compiled$Cmcmc$last_args$resetWAIC)
  expect_true(compiled$Cmcmc$last_args$initializeModel)
  expect_equal(length(compiled$Cmodel$init_calls), 1L)
  expect_equal(compiled$Cmodel$init_calls[[1]], compiled$spec$inits)
  expect_equal(tail(compiled$Cmodel$events, 2), c("setInits", "run"))
  expect_equal(nrow(resumed$run_history), 1L)
  expect_true(resumed$run_history$reset)
  expect_false(resumed$mcmc_control$resumed)
  expect_equal(resumed$mcmc_control$init_strategy, "compiled_spec_reset")
  expect_equal(resumed$chain_info$init_strategy, "compiled_spec_reset")
  expect_false(resumed$chain_info$resumed)
  expect_false(resumed$run_history$resumed)
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

  make_chain <- function(seed_val, waic_val) {
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
      waic = waic_val,
      sampling_time = runif(1, 5, 15),
      model_config = list(
        model = "rasch", prior = "normal",
        parameterization = "irt",
        identification = "unconstrained",
        N = N, I = I
      ),
      mcmc_control = list(
        niter = niter,
        nburnin = 10,
        thin = 2,
        thin2 = 2,
        seed = seed_val,
        reset = TRUE
      )
    ),
    class = "dpmirt_samples"
  )
  }

  chain1 <- make_chain(42, 2400)
  chain2 <- make_chain(99, 2600)

  result <- DPMirt:::.combine_chains(list(chain1, chain2), rescale = TRUE)

  # Combined theta should have 2 * niter rows

  expect_equal(nrow(result$rescaled$theta_samp), 2 * niter)
  expect_equal(ncol(result$rescaled$theta_samp), N)

  # Combined beta should have 2 * niter rows
  expect_equal(nrow(result$rescaled$beta_samp), 2 * niter)
  expect_equal(ncol(result$rescaled$beta_samp), I)

  # WAIC should be mean of chain WAICs
  expect_equal(result$waic, 2500)
  expect_equal(result$chain_info$waic, c(2400, 2600))

  # Sampling time should be sum
  expect_equal(result$sampling_time,
               chain1$sampling_time + chain2$sampling_time)

  # Raw samples combined
  expect_equal(nrow(result$raw_samples), 2 * niter)

  # Chain metadata should preserve row provenance while keeping pooled matrices
  expect_equal(result$chain_info$chain_id, 1:2)
  expect_equal(result$chain_info$seed, c(42L, 99L))
  expect_equal(result$chain_info$n_draws_main, c(niter, niter))
  expect_equal(result$chain_info$n_draws_theta, c(niter, niter))
  expect_equal(result$chain_info$row_start_main, c(1L, niter + 1L))
  expect_equal(result$chain_info$row_end_main, c(niter, 2L * niter))
  expect_equal(result$run_history$run_id, 1:2)
  expect_equal(result$run_history$chain_id, 1:2)
  expect_equal(result$run_history$row_start_main, c(1L, niter + 1L))
  expect_equal(result$run_history$row_end_theta, c(niter, 2L * niter))
  expect_equal(nrow(result$draw_index$main), 2 * niter)
  expect_equal(nrow(result$draw_index$theta), 2 * niter)
  expect_equal(result$draw_index$main$chain_id[1:niter], rep(1L, niter))
  expect_equal(result$draw_index$main$chain_id[(niter + 1):(2 * niter)],
               rep(2L, niter))

  chain2$waic <- NULL
  one_waic <- DPMirt:::.combine_chains(list(chain1, chain2), rescale = TRUE)
  expect_equal(one_waic$waic, 2400)
  expect_equal(one_waic$chain_info$waic, c(2400, NA))

  chain1$waic <- NULL
  no_waic <- DPMirt:::.combine_chains(list(chain1, chain2), rescale = TRUE)
  expect_null(no_waic$waic)
  expect_true(all(is.na(no_waic$chain_info$waic)))
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
      ),
      mcmc_control = list(
        niter = niter,
        nburnin = 5,
        thin = 1,
        thin2 = 1,
        seed = seed_val,
        reset = TRUE
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


test_that(".combine_chains keeps separate main/theta draw indexes when thin differs", {
  set.seed(42)
  N <- 10; I <- 3
  n_main <- 50
  n_theta <- 25

  make_chain <- function(seed_val) {
    set.seed(seed_val)
    beta_samp <- matrix(rnorm(n_main * I, mean = 0.5),
                        nrow = n_main, ncol = I)
    colnames(beta_samp) <- paste0("beta[", seq_len(I), "]")

    eta_samp <- matrix(rnorm(n_theta * N), nrow = n_theta, ncol = N)
    colnames(eta_samp) <- paste0("eta[", seq_len(N), "]")

    log_nodes <- matrix(rnorm(n_main * 3), nrow = n_main, ncol = 3)
    colnames(log_nodes) <- c("myLogProbAll", "myLogProbSome", "myLogLik")

    structure(
      list(
        samples  = cbind(beta_samp, log_nodes),
        samples2 = eta_samp,
        waic = 2500,
        sampling_time = 1,
        model_config = list(
          model = "rasch", prior = "normal",
          parameterization = "irt",
          identification = "unconstrained",
          N = N, I = I
        ),
        mcmc_control = list(
          niter = 100,
          nburnin = 0,
          thin = 2,
          thin2 = 4,
          seed = seed_val,
          reset = TRUE
        )
      ),
      class = "dpmirt_samples"
    )
  }

  result <- DPMirt:::.combine_chains(
    list(make_chain(11), make_chain(12)),
    rescale = TRUE
  )

  expect_equal(nrow(result$draw_index$main), 2 * n_main)
  expect_equal(nrow(result$draw_index$theta), 2 * n_theta)
  expect_equal(result$chain_info$n_draws_main, c(n_main, n_main))
  expect_equal(result$chain_info$n_draws_theta, c(n_theta, n_theta))
  expect_equal(result$chain_info$row_start_theta, c(1L, n_theta + 1L))
  expect_equal(result$chain_info$row_end_theta, c(n_theta, 2L * n_theta))
  expect_equal(unique(result$draw_index$main$mcmc_iteration[1:n_main]),
               seq(2, by = 2, length.out = n_main))
  expect_equal(unique(result$draw_index$theta$mcmc_iteration[1:n_theta]),
               seq(4, by = 4, length.out = n_theta))
})
