# W-03 (V3 contract, decision D-V3-02): the nchains = 4 path must reproduce
# four matched-seed single-chain calls BIT-EXACTLY (per-chain seeds and inits
# verifiably applied; draws chain-labeled; no cross-chain state leakage from
# the reused compiled object). Passing this fixture is what promotes
# dpmirt(nchains = 4) from nchains_independence_unvalidated to the V3
# production path; failure falls back to four sequential one-chain calls.
#
# C-09 hardening: (a) EVERY posterior family is compared (theta, beta, lambda,
# other) — not just theta/beta, so a divergent 2PL discrimination or DPM
# mixture posterior cannot slip through; (b) chain labels are REQUIRED — the
# old positional-block fallback is gone, an unlabeled fit is an explicit
# failure; (c) chain_info must record the exact per-chain seed consumed; and
# (d) the explicit seed-VECTOR contract is exercised.
#
# Small models keep this affordable; it still compiles 4 models per config
# (1 multi-chain + 3 single-chain re-compiles; chain 1 shares the multi-chain
# compile check via seed identity).

# Chain-row lookup: REQUIRE an explicit chain_id label. No positional fallback —
# if a fit cannot say which rows belong to which chain, that is a failure, not
# something to paper over by assuming contiguous appends. The draw index is a
# per-family list (draw_index$main / draw_index$theta), each carrying its own
# chain_id, because main monitors and theta may be thinned independently.
.w03_chain_rows <- function(fit, ch, family = c("main", "theta")) {
  family <- match.arg(family)
  di <- fit$draw_index[[family]]
  if (is.null(di) || is.null(di$chain_id)) {
    stop("W03: fit$draw_index$", family,
         "$chain_id is required (no positional fallback)")
  }
  which(di$chain_id == ch)
}

# Compare every posterior family present in the single-chain fit against the
# corresponding chain-labeled rows of the multi-chain fit. theta_samp follows
# the theta (monitors2) index; the item/hyperparameter families follow the main
# index.
.w03_expect_families <- function(fit4, fit1, main_rows, theta_rows, tag) {
  fam_rows <- list(theta_samp = theta_rows, beta_samp = main_rows,
                   lambda_samp = main_rows, other_samp = main_rows)
  for (fam in names(fam_rows)) {
    a <- fit1[[fam]]
    if (is.null(a) || (is.matrix(a) && ncol(a) == 0)) next
    expect_identical(
      unname(fit4[[fam]][fam_rows[[fam]], , drop = FALSE]), unname(a),
      label = paste0(tag, " ", fam)
    )
  }
}

test_that("nchains=4 bit-equals four matched-seed single-chain fits (all families)", {
  skip_on_cran()

  set.seed(31)
  N <- 24; I <- 8
  theta0 <- rnorm(N); beta0 <- seq(-1.2, 1.2, length.out = I)
  p <- plogis(outer(theta0, beta0, `-`))
  y <- matrix(rbinom(N * I, 1, as.vector(p)), N, I)

  configs <- list(
    list(model = "rasch", prior = "normal"),
    list(model = "rasch", prior = "dpm"),
    list(model = "2pl",   prior = "normal"),
    list(model = "2pl",   prior = "dpm")
  )

  base_seed <- 777L
  niter <- 300L; nburnin <- 100L; thin <- 1L
  draws_per_chain <- (niter - nburnin) %/% thin

  for (cfg in configs) {
    args <- list(data = y, model = cfg$model, prior = cfg$prior,
                 niter = niter, nburnin = nburnin, thin = thin,
                 seed = base_seed, verbose = FALSE,
                 compute_waic = FALSE, compute_dp_density = FALSE)
    if (cfg$prior == "dpm") {
      args$alpha_prior <- c(a = 2, b = 4)
      args$M <- 15L
    }

    fit4 <- suppressWarnings(do.call(dpmirt, c(args, list(nchains = 4L))))
    expect_identical(nrow(fit4$theta_samp), 4L * draws_per_chain)

    # chain labels are mandatory (C-09)
    expect_false(is.null(fit4$draw_index))
    expect_false(is.null(fit4$draw_index$main$chain_id))
    expect_false(is.null(fit4$draw_index$theta$chain_id))
    expect_identical(sort(unique(fit4$draw_index$main$chain_id)), 1:4)
    expect_identical(sort(unique(fit4$draw_index$theta$chain_id)), 1:4)

    for (ch in 1:4) {
      args1 <- args
      args1$seed <- DPMirt:::.chain_seed(base_seed, ch)
      fit1 <- suppressWarnings(do.call(dpmirt, c(args1, list(nchains = 1L))))
      main_rows  <- .w03_chain_rows(fit4, ch, "main")
      theta_rows <- .w03_chain_rows(fit4, ch, "theta")
      expect_identical(length(main_rows), as.integer(draws_per_chain))
      .w03_expect_families(fit4, fit1, main_rows, theta_rows,
                           tag = paste0(cfg$model, "/", cfg$prior, " ch", ch))
      # the seed the multi-chain run recorded for this chain matches the
      # single-chain seed that reproduced it (seed reconciliation).
      ci <- fit4$chain_info
      expect_identical(as.integer(ci$seed[ch]),
                       as.integer(DPMirt:::.chain_seed(base_seed, ch)))
    }
  }
})

test_that("explicit per-chain seed VECTOR is consumed verbatim (C-09 contract)", {
  skip_on_cran()

  set.seed(7)
  N <- 20; I <- 6
  theta0 <- rnorm(N); beta0 <- seq(-1, 1, length.out = I)
  y <- matrix(rbinom(N * I, 1, plogis(outer(theta0, beta0, `-`))), N, I)

  # four arbitrary, NON-contiguous seeds (so an internal `+ch-1` rule could not
  # coincidentally reproduce them) — the V3 seed-tree case.
  seed_vec <- c(4001L, 8887L, 123L, 55501L)
  niter <- 260L; nburnin <- 100L; dpc <- niter - nburnin

  fit4 <- suppressWarnings(dpmirt(y, model = "2pl", prior = "normal",
                                  niter = niter, nburnin = nburnin, thin = 1L,
                                  nchains = 4L, seed = seed_vec, verbose = FALSE,
                                  compute_waic = FALSE, compute_dp_density = FALSE))
  # chain_info records exactly the vector we passed
  expect_identical(as.integer(fit4$chain_info$seed), seed_vec)

  for (ch in 1:4) {
    fit1 <- suppressWarnings(dpmirt(y, model = "2pl", prior = "normal",
                                    niter = niter, nburnin = nburnin, thin = 1L,
                                    nchains = 1L, seed = seed_vec[ch],
                                    verbose = FALSE, compute_waic = FALSE,
                                    compute_dp_density = FALSE))
    main_rows  <- .w03_chain_rows(fit4, ch, "main")
    theta_rows <- .w03_chain_rows(fit4, ch, "theta")
    expect_identical(length(main_rows), as.integer(dpc))
    .w03_expect_families(fit4, fit1, main_rows, theta_rows,
                         tag = paste0("vec ch", ch))
  }
})

test_that("a malformed seed length is rejected", {
  expect_error(DPMirt:::.resolve_chain_seeds(c(1L, 2L, 3L), 4L),
               "length 1 or nchains")
  # length-1 and length-nchains are both fine
  expect_length(DPMirt:::.resolve_chain_seeds(5L, 4L), 4L)
  expect_identical(DPMirt:::.resolve_chain_seeds(c(9L, 8L), 2L), list(9L, 8L))
})

test_that("duplicate and overflowing effective chain seeds are refused", {
  expect_error(
    DPMirt:::.resolve_chain_seeds(c(4001L, 4001L), 2L),
    "must be unique"
  )
  expect_error(
    DPMirt:::.resolve_chain_seeds(.Machine$integer.max, 2L),
    "exceeds the maximum"
  )
  expect_error(
    DPMirt:::.resolve_chain_seeds(c(1, 2.5), 2L),
    "finite whole"
  )
  expect_error(
    DPMirt:::.resolve_chain_seeds(-1L, 1L),
    "must be in"
  )
})

test_that("seed collisions are refused before model specification or compile", {
  y <- matrix(c(0, 1, 1, 0), nrow = 2L)

  # A valid response matrix would proceed into specification/compilation if
  # seed resolution were still late. The collision-specific refusal proves the
  # launch is stopped at the API boundary instead.
  expect_error(
    dpmirt(y, nchains = 2L, seed = c(77L, 77L),
           niter = 4L, nburnin = 2L, verbose = FALSE),
    "duplicate effective chain seeds"
  )
})
