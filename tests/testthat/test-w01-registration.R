# W-01 (V3 contract): dBernoulliVector registration must be deterministic —
# namespace bindings resolvable by NIMBLE at registration AND model-build time.
# The fresh-Rscript and N-way parallel checks run in the V3 P-0 fixture suite;
# this in-package test asserts the structural preconditions and an in-session
# model build through the registration path.

test_that("dBernoulliVector/rBernoulliVector are namespace bindings", {
  ns <- asNamespace("DPMirt")
  expect_true(exists("dBernoulliVector", envir = ns, inherits = FALSE))
  expect_true(exists("rBernoulliVector", envir = ns, inherits = FALSE))
  # install-time nimbleFunctions serialize as RC functions (no setup code)
  expect_true(nimble:::is.rcf(get("dBernoulliVector", envir = ns)))
  expect_true(nimble:::is.rcf(get("rBernoulliVector", envir = ns)))
})

test_that("registration uses a dedicated environment and never globalenv", {
  binding_names <- c("dBernoulliVector", "rBernoulliVector")
  had_global <- vapply(binding_names, exists, logical(1),
                       envir = globalenv(), inherits = FALSE)
  old_global <- lapply(binding_names[had_global], get,
                       envir = globalenv(), inherits = FALSE)
  names(old_global) <- binding_names[had_global]
  on.exit({
    for (nm in binding_names) {
      if (exists(nm, envir = globalenv(), inherits = FALSE)) {
        rm(list = nm, envir = globalenv())
      }
    }
    for (nm in names(old_global)) {
      assign(nm, old_global[[nm]], envir = globalenv())
    }
  }, add = TRUE)
  if (any(had_global)) {
    rm(list = binding_names[had_global], envir = globalenv())
  }

  env1 <- DPMirt:::.register_dBernoulliVector()
  env2 <- DPMirt:::.register_dBernoulliVector()
  expect_identical(env1, env2)
  expect_identical(environmentName(env1),
                   DPMirt:::.dpmirt_distribution_env_name())
  expect_identical(search()[[2L]],
                   DPMirt:::.dpmirt_distribution_env_name())
  expect_true(isTRUE(DPMirt:::.nimble_cache$dBernoulliVector_registered))
  expect_false(any(vapply(binding_names, exists, logical(1),
                          envir = globalenv(), inherits = FALSE)))
  ns <- asNamespace("DPMirt")
  for (nm in binding_names) {
    expect_identical(get(nm, envir = env1, inherits = FALSE),
                     get(nm, envir = ns, inherits = FALSE))
  }
})

test_that("a foreign same-name environment is a hard collision", {
  env_name <- DPMirt:::.dpmirt_distribution_env_name()
  DPMirt:::.detach_dpmirt_distribution_env(deregister = TRUE)
  attach(list(fake_owner = TRUE), pos = 2L, name = env_name,
         warn.conflicts = FALSE)
  on.exit({
    if (env_name %in% search()) {
      detach(env_name, character.only = TRUE)
    }
    DPMirt:::.register_dBernoulliVector()
  }, add = TRUE)
  expect_error(
    DPMirt:::.register_dBernoulliVector(),
    "search-path collision"
  )
})

test_that("a non-identical global d/r binding is a hard collision", {
  binding_names <- c("dBernoulliVector", "rBernoulliVector")
  expect_false(any(vapply(binding_names, exists, logical(1),
                          envir = globalenv(), inherits = FALSE)))
  assign("dBernoulliVector", function(...) "foreign", envir = globalenv())
  on.exit({
    if (exists("dBernoulliVector", envir = globalenv(), inherits = FALSE)) {
      rm("dBernoulliVector", envir = globalenv())
    }
    DPMirt:::.register_dBernoulliVector()
  }, add = TRUE)
  expect_error(
    DPMirt:::.register_dBernoulliVector(),
    "global custom distribution binding collision.*dBernoulliVector"
  )
})

test_that("deleted binding and NIMBLE registry entry are repaired", {
  env <- DPMirt:::.register_dBernoulliVector()
  rm("dBernoulliVector", envir = env)
  nimble::deregisterDistributions(
    "dBernoulliVector", userEnv = env, warn = FALSE
  )
  expect_false(exists("dBernoulliVector", envir = env, inherits = FALSE))

  repaired <- DPMirt:::.register_dBernoulliVector()
  expect_identical(repaired, env)
  expect_identical(
    get("dBernoulliVector", envir = repaired, inherits = FALSE),
    get("dBernoulliVector", envir = asNamespace("DPMirt"), inherits = FALSE)
  )

  code <- nimble::nimbleCode({
    y[1:I] ~ dBernoulliVector(p[1:I])
  })
  m <- nimble::nimbleModel(
    code,
    constants = list(I = 4),
    data = list(y = c(1, 0, 1, 0)),
    inits = list(p = rep(0.5, 4)),
    calculate = TRUE,
    userEnv = repaired
  )
  expect_true(is.finite(m$calculate("y")))
})

test_that("an altered owned binding is refused rather than overwritten", {
  env <- DPMirt:::.register_dBernoulliVector()
  assign("rBernoulliVector", function(...) "foreign", envir = env)
  on.exit({
    if (exists("rBernoulliVector", envir = env, inherits = FALSE)) {
      rm("rBernoulliVector", envir = env)
    }
    DPMirt:::.register_dBernoulliVector()
  }, add = TRUE)
  expect_error(
    DPMirt:::.register_dBernoulliVector(),
    "binding collision.*rBernoulliVector"
  )
})

test_that("constrained 2PL normal and DPM configurations compile and fit", {
  skip_on_cran()
  # 2PL + constrained_item is the ONLY configuration that routes through
  # dBernoulliVector. Before the registration exposed the d/r nimbleFunctions
  # through an authenticated search environment, this failed at compileNimble()
  # with
  # "rBernoulliVector is not available or its output type is unknown". This is
  # the production-faithful check: the full spec -> build -> compile -> sample
  # pipeline must work for both V3 prior arms and return a discrimination
  # posterior without creating global bindings.
  set.seed(11)
  N <- 24L; I <- 6L
  y <- matrix(rbinom(N * I, 1,
                     plogis(outer(rnorm(N), seq(-1.2, 1.2, length.out = I),
                                  `-`))), N, I)
  for (prior in c("normal", "dpm")) {
    args <- list(
      data = y, model = "2pl", prior = prior,
      identification = "constrained_item",
      niter = 100L, nburnin = 50L, thin = 1L, nchains = 1L, seed = 7L,
      verbose = FALSE, compute_waic = FALSE, compute_dp_density = FALSE
    )
    if (identical(prior, "dpm")) {
      args$alpha_prior <- c(a = 2, b = 4)
      args$M <- 15L
    }
    fit <- suppressWarnings(do.call(dpmirt, args))
    expect_false(is.null(fit$lambda_samp), label = prior)
    expect_identical(ncol(fit$lambda_samp), I, label = prior)
    expect_true(all(is.finite(fit$lambda_samp)), label = prior)
  }
  expect_false(exists("dBernoulliVector", envir = globalenv(),
                      inherits = FALSE))
  expect_false(exists("rBernoulliVector", envir = globalenv(),
                      inherits = FALSE))
})
