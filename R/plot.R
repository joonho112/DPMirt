# ============================================================================
# Module 13: Visualization
# ============================================================================
# Blueprint: Section 6 - Module 13
# Phase 2: Basic DPM plots (density overlay, cluster count, alpha trace)
# Phase 7B: Full visualization suite with base R + ggplot2 backends
# ============================================================================


# --------------------------------------------------------------------------
# Package-level color palette (Paganin-consistent)
# --------------------------------------------------------------------------

#' DPMirt color palette (internal)
#' @noRd
.dpmirt_colors <- list(
  parametric     = "#56B4E9",
  semiparametric = "#E69F00",
  reference      = "gray50",
  ci_fill        = "steelblue",
  trace          = "steelblue",
  items          = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a",
                     "#66a61e", "#e6ab02", "#a6761d", "#666666")
)


#' Check if ggplot2 is available
#' @noRd
.check_gg <- function() {
  requireNamespace("ggplot2", quietly = TRUE)
}


# --------------------------------------------------------------------------
# S3 dispatcher: plot.dpmirt_fit
# --------------------------------------------------------------------------

#' Plot a DPMirt fit object
#'
#' Produces visualizations for fitted IRT models. Supports 12 plot types
#' with automatic selection between base R and ggplot2 backends.
#'
#' @param x A \code{dpmirt_fit} object.
#' @param type Character. Plot type. One of:
#'   \describe{
#'     \item{\code{"density"}}{Kernel density of posterior mean theta.}
#'     \item{\code{"items"}}{Item difficulty estimates with error bars.}
#'     \item{\code{"trace"}}{Log-likelihood MCMC trace.}
#'     \item{\code{"clusters"}}{Cluster count trace and histogram (DPM only).}
#'     \item{\code{"dp_density"}}{DP mixture density with credible band (DPM only).}
#'     \item{\code{"icc"}}{Item Characteristic Curves.}
#'     \item{\code{"wright_map"}}{Person-Item map (Wright map).}
#'     \item{\code{"parameter_trace"}}{Individual parameter MCMC traces.}
#'     \item{\code{"caterpillar"}}{Sorted estimates with credible intervals.}
#'     \item{\code{"density_compare"}}{Posterior density vs reference overlay.}
#'     \item{\code{"info"}}{Test Information Function.}
#'     \item{\code{"pp_check"}}{Posterior predictive check.}
#'   }
#' @param engine Character. Plotting backend: \code{"auto"} (default) uses
#'   ggplot2 if available, \code{"base"} forces base R, \code{"ggplot2"}
#'   requires ggplot2.
#' @param ... Additional arguments passed to the specific plotting function.
#'
#' @return Invisibly returns the plot object (ggplot) or NULL (base R).
#'
#' @examples
#' \dontrun{
#' sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
#' fit <- dpmirt(sim$response, model = "rasch", prior = "normal",
#'               niter = 5000, nburnin = 1000, seed = 123)
#'
#' # Theta density
#' plot(fit, type = "density")
#'
#' # Item difficulty estimates
#' plot(fit, type = "items")
#'
#' # MCMC trace
#' plot(fit, type = "trace")
#'
#' # Force base R backend
#' plot(fit, type = "density", engine = "base")
#' }
#'
#' @family visualization
#' @seealso \code{\link{dpmirt_plot_density}}, \code{\link{dpmirt_plot_items}},
#'   \code{\link{dpmirt_plot_trace}}, \code{\link{dpmirt_plot_icc}}
#' @export
plot.dpmirt_fit <- function(x,
                            type = c("density", "items", "trace",
                                     "clusters", "dp_density",
                                     "icc", "wright_map", "parameter_trace",
                                     "caterpillar", "density_compare",
                                     "info", "pp_check"),
                            engine = c("auto", "base", "ggplot2"),
                            ...) {

  type   <- match.arg(type)
  engine <- match.arg(engine)

  # Resolve engine
  use_gg <- switch(engine,
    "auto"    = .check_gg(),
    "ggplot2" = {
      if (!.check_gg())
        stop("ggplot2 is required. Install with: install.packages('ggplot2')",
             call. = FALSE)
      TRUE
    },
    "base"    = FALSE
  )

  if (use_gg) {
    .dispatch_gg(x, type, ...)
  } else {
    .dispatch_base(x, type, ...)
  }
}


#' Dispatch to base R plot functions
#' @noRd
.dispatch_base <- function(x, type, ...) {
  switch(type,
    "density"          = .plot_density(x, ...),
    "items"            = .plot_items(x, ...),
    "trace"            = .plot_trace(x, ...),
    "clusters"         = .plot_clusters(x, ...),
    "dp_density"       = .plot_dp_density(x, ...),
    "icc"              = .plot_icc(x, ...),
    "wright_map"       = .plot_wright_map(x, ...),
    "parameter_trace"  = .plot_parameter_trace(x, ...),
    "caterpillar"      = .plot_caterpillar(x, ...),
    "density_compare"  = .plot_density_compare(x, ...),
    "info"             = .plot_info(x, ...),
    "pp_check"         = .plot_pp_check(x, ...)
  )
  invisible(NULL)
}


#' Dispatch to ggplot2 plot functions
#' @noRd
.dispatch_gg <- function(x, type, ...) {
  fn_name <- paste0("dpmirt_plot_", type)
  fn <- get(fn_name, envir = asNamespace("DPMirt"))
  p <- fn(x, ...)
  print(p)
  invisible(p)
}


# ============================================================================
# Base R: Existing Plot Functions (Phase 2)
# ============================================================================

#' Posterior mean density plot (base R)
#' @noRd
.plot_density <- function(fit, ...) {
  if (is.null(fit$theta_samp)) {
    stop("No theta samples available for density plot.", call. = FALSE)
  }

  theta_pm <- colMeans(fit$theta_samp)

  plot(density(theta_pm),
       main = paste0(toupper(fit$config$model), " (", fit$config$prior,
                     ") -- Posterior Mean Density"),
       xlab = expression(theta),
       ylab = "Density",
       lwd = 2,
       ...)
}


#' Item parameter plot (base R)
#' @noRd
.plot_items <- function(fit, ...) {
  if (is.null(fit$beta_samp)) {
    stop("No beta samples available.", call. = FALSE)
  }

  beta_pm  <- colMeans(fit$beta_samp)
  beta_psd <- apply(fit$beta_samp, 2, sd)
  I <- length(beta_pm)

  plot(seq_len(I), beta_pm,
       ylim = range(c(beta_pm - 2 * beta_psd, beta_pm + 2 * beta_psd)),
       pch = 19, xlab = "Item", ylab = expression(beta),
       main = "Item Difficulty Estimates",
       ...)
  arrows(seq_len(I), beta_pm - 2 * beta_psd,
         seq_len(I), beta_pm + 2 * beta_psd,
         length = 0.05, angle = 90, code = 3, col = "gray50")
}


#' Log-likelihood trace plot (base R)
#' @noRd
.plot_trace <- function(fit, ...) {
  if (is.null(fit$loglik_trace)) {
    stop("No log-likelihood trace available.", call. = FALSE)
  }

  plot(fit$loglik_trace,
       type = "l",
       xlab = "Iteration", ylab = "Log-Likelihood",
       main = "Log-Likelihood Trace",
       col = .dpmirt_colors$trace,
       ...)
}


#' Cluster count trace and histogram (base R)
#' @noRd
.plot_clusters <- function(fit, ...) {
  if (fit$config$prior != "dpm") {
    stop("Cluster plot is only available for DPM models.", call. = FALSE)
  }

  if (is.null(fit$cluster_info) || is.null(fit$cluster_info$n_clusters)) {
    stop("No cluster information available. ",
         "Ensure the model was fit with prior = 'dpm'.", call. = FALSE)
  }

  n_cl <- fit$cluster_info$n_clusters

  old_par <- par(mfrow = c(1, 2))
  on.exit(par(old_par))

  plot(n_cl, type = "l",
       xlab = "Iteration", ylab = "Number of Clusters",
       main = "Cluster Count Trace",
       col = .dpmirt_colors$trace, ...)

  hist(n_cl,
       breaks = seq(min(n_cl) - 0.5, max(n_cl) + 0.5, by = 1),
       xlab = "Number of Clusters",
       main = "Cluster Count Distribution",
       col = "lightblue", border = "white",
       ...)
  abline(v = mean(n_cl), col = "red", lwd = 2, lty = 2)
}


#' DP mixture density with credible band (base R)
#' @noRd
.plot_dp_density <- function(fit, ...) {
  if (is.null(fit$dp_density)) {
    stop("No DP density available. ",
         "Compute it with dpmirt_dp_density(fit).", call. = FALSE)
  }

  dpd <- fit$dp_density

  plot(dpd$grid, dpd$density_mean,
       type = "l", lwd = 2,
       xlab = expression(theta),
       ylab = "Density",
       main = "DP Mixture Density (Posterior)",
       ylim = c(0, max(dpd$density_upper) * 1.05),
       col = .dpmirt_colors$ci_fill,
       ...)

  polygon(c(dpd$grid, rev(dpd$grid)),
          c(dpd$density_lower, rev(dpd$density_upper)),
          col = adjustcolor(.dpmirt_colors$ci_fill, alpha.f = 0.2),
          border = NA)

  curve(dnorm(x), add = TRUE, col = .dpmirt_colors$reference, lty = 2,
        lwd = 1.5)

  legend("topright",
         legend = c("DP Posterior Mean",
                    paste0(dpd$ci_level * 100, "% CI"),
                    "N(0,1) Reference"),
         col = c(.dpmirt_colors$ci_fill,
                 adjustcolor(.dpmirt_colors$ci_fill, alpha.f = 0.2),
                 .dpmirt_colors$reference),
         lty = c(1, NA, 2),
         lwd = c(2, NA, 1.5),
         pch = c(NA, 15, NA),
         pt.cex = 2,
         bty = "n")
}


# ============================================================================
# Base R: New Plot Functions (Phase 7B)
# ============================================================================

# --------------------------------------------------------------------------
# ICC helper: compute P(correct | theta) for one item
# --------------------------------------------------------------------------

#' Compute ICC probability
#' @noRd
.icc_prob <- function(theta, beta, lambda = 1, delta = 0) {
  eta <- lambda * (theta - beta)
  p_star <- 1 / (1 + exp(-eta))
  delta + (1 - delta) * p_star
}


#' Item Characteristic Curves (base R)
#'
#' @param fit A \code{dpmirt_fit} object.
#' @param items Integer vector of item indices to plot. Default: all items
#'   (up to 10).
#' @param theta_range Numeric vector of length 2. Range of theta axis.
#' @param n_points Number of grid points.
#' @param ... Additional graphical parameters.
#'
#' @noRd
.plot_icc <- function(fit, items = NULL, theta_range = c(-4, 4),
                      n_points = 201, ...) {
  if (is.null(fit$beta_samp)) {
    stop("No beta samples available for ICC plot.", call. = FALSE)
  }

  beta_pm <- colMeans(fit$beta_samp)
  I <- length(beta_pm)

  lambda_pm <- if (!is.null(fit$lambda_samp)) colMeans(fit$lambda_samp) else
    rep(1, I)
  delta_pm  <- if (!is.null(fit$delta_samp)) colMeans(fit$delta_samp) else
    rep(0, I)

  if (is.null(items)) items <- seq_len(min(I, 10))
  items <- items[items >= 1 & items <= I]

  theta_grid <- seq(theta_range[1], theta_range[2], length.out = n_points)

  # Compute probabilities
  prob_mat <- sapply(items, function(j) {
    .icc_prob(theta_grid, beta_pm[j], lambda_pm[j], delta_pm[j])
  })

  cols <- rep_len(.dpmirt_colors$items, length(items))

  matplot(theta_grid, prob_mat,
          type = "l", lty = 1, lwd = 2, col = cols,
          xlab = expression(theta), ylab = "P(correct)",
          main = paste0("Item Characteristic Curves (",
                        toupper(fit$config$model), ")"),
          ylim = c(0, 1),
          ...)
  abline(h = 0.5, col = "gray80", lty = 3)
  legend("bottomright",
         legend = paste0("Item ", items),
         col = cols, lty = 1, lwd = 2,
         bty = "n", cex = 0.8,
         ncol = if (length(items) > 5) 2 else 1)
}


#' Person-Item Map / Wright Map (base R)
#'
#' @param fit A \code{dpmirt_fit} object.
#' @param bins Number of histogram bins for person distribution.
#' @param item_labels Optional character vector of item labels.
#' @param ... Additional graphical parameters.
#'
#' @noRd
.plot_wright_map <- function(fit, bins = 30, item_labels = NULL, ...) {
  if (is.null(fit$theta_samp) || is.null(fit$beta_samp)) {
    stop("No theta/beta samples available for Wright map.", call. = FALSE)
  }

  theta_pm <- colMeans(fit$theta_samp)
  beta_pm  <- colMeans(fit$beta_samp)
  I <- length(beta_pm)

  if (is.null(item_labels)) item_labels <- paste0("I", seq_len(I))

  # Shared y-axis range
  y_range <- range(c(theta_pm, beta_pm)) + c(-0.5, 0.5)

  old_par <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 0.5))
  on.exit(par(old_par))

  # Left panel: Person ability histogram (horizontal)
  h <- hist(theta_pm, breaks = bins, plot = FALSE)
  barplot(h$counts, horiz = TRUE, space = 0,
          names.arg = round(h$mids, 1),
          las = 1, col = adjustcolor(.dpmirt_colors$parametric, 0.6),
          border = "white",
          xlab = "Frequency", ylab = expression(theta),
          main = "Persons")

  # Right panel: Item difficulties
  par(mar = c(4, 0.5, 3, 2))
  beta_ord <- order(beta_pm)
  plot(rep(0, I), beta_pm[beta_ord],
       pch = 18, cex = 1.5,
       col = .dpmirt_colors$semiparametric,
       xlim = c(-1, 1), ylim = y_range,
       xlab = "", ylab = "", yaxt = "n",
       main = "Items", xaxt = "n", ...)
  text(rep(0.3, I), beta_pm[beta_ord],
       labels = item_labels[beta_ord],
       cex = 0.7, adj = 0)
  abline(h = mean(theta_pm), col = "gray60", lty = 2)
}


#' Individual parameter MCMC trace (base R)
#'
#' @param fit A \code{dpmirt_fit} object.
#' @param param Character. Which parameter: \code{"beta"}, \code{"theta"},
#'   \code{"lambda"}, \code{"delta"}.
#' @param indices Integer vector of parameter indices to plot. Default: 1:4.
#' @param ... Additional graphical parameters.
#'
#' @noRd
.plot_parameter_trace <- function(fit,
                                  param = c("beta", "theta",
                                            "lambda", "delta"),
                                  indices = 1:4,
                                  ...) {
  param <- match.arg(param)

  samp <- switch(param,
    "beta"   = fit$beta_samp,
    "theta"  = fit$theta_samp,
    "lambda" = fit$lambda_samp,
    "delta"  = fit$delta_samp
  )

  if (is.null(samp)) {
    stop(sprintf("No %s samples available.", param), call. = FALSE)
  }

  n_cols <- ncol(samp)
  indices <- indices[indices >= 1 & indices <= n_cols]
  if (length(indices) == 0) {
    stop("No valid indices specified.", call. = FALSE)
  }

  n_panels <- length(indices)
  n_row <- ceiling(n_panels / 2)
  n_col <- min(n_panels, 2)

  old_par <- par(mfrow = c(n_row, n_col), mar = c(3, 4, 2, 1))
  on.exit(par(old_par))

  for (idx in indices) {
    plot(samp[, idx], type = "l",
         col = adjustcolor(.dpmirt_colors$trace, 0.7),
         xlab = "", ylab = "",
         main = paste0(param, "[", idx, "]"),
         ...)
    abline(h = mean(samp[, idx]), col = "red", lty = 2)
    mtext("Iteration", side = 1, line = 2, cex = 0.7)
  }
}


#' Caterpillar plot with sorted estimates and CI (base R)
#'
#' @param fit A \code{dpmirt_fit} object.
#' @param param Character. Which parameter: \code{"beta"}, \code{"theta"},
#'   \code{"lambda"}, \code{"delta"}.
#' @param sort Logical. Sort by estimate magnitude. Default TRUE.
#' @param ci_level Numeric. Credible interval level. Default 0.95.
#' @param max_show Maximum number of parameters to display. Default 50.
#' @param ... Additional graphical parameters.
#'
#' @noRd
.plot_caterpillar <- function(fit,
                              param = c("beta", "theta",
                                        "lambda", "delta"),
                              sort = TRUE, ci_level = 0.95,
                              max_show = 50, ...) {
  param <- match.arg(param)

  samp <- switch(param,
    "beta"   = fit$beta_samp,
    "theta"  = fit$theta_samp,
    "lambda" = fit$lambda_samp,
    "delta"  = fit$delta_samp
  )

  if (is.null(samp)) {
    stop(sprintf("No %s samples available.", param), call. = FALSE)
  }

  lower_q <- (1 - ci_level) / 2
  upper_q <- 1 - lower_q

  pm    <- colMeans(samp)
  lower <- apply(samp, 2, quantile, probs = lower_q)
  upper <- apply(samp, 2, quantile, probs = upper_q)

  n_total <- length(pm)
  show_idx <- seq_len(n_total)

  if (n_total > max_show) {
    show_idx <- round(seq(1, n_total, length.out = max_show))
  }

  if (sort) {
    ord <- order(pm[show_idx])
    show_idx <- show_idx[ord]
  }

  n <- length(show_idx)
  y_pos <- seq_len(n)

  plot(pm[show_idx], y_pos,
       xlim = range(c(lower[show_idx], upper[show_idx])),
       pch = 19, cex = 0.6,
       yaxt = "n", xlab = param, ylab = "",
       main = paste0("Caterpillar: ", param,
                     " (", round(ci_level * 100), "% CI)"),
       ...)
  segments(lower[show_idx], y_pos, upper[show_idx], y_pos,
           col = adjustcolor(.dpmirt_colors$ci_fill, 0.5))
  axis(2, at = y_pos, labels = show_idx, las = 1, cex.axis = 0.6)
  abline(v = 0, col = "gray70", lty = 2)
}


#' Posterior density comparison (base R)
#'
#' Overlays the posterior mean theta density with a Normal reference,
#' or with the DP mixture density if available.
#'
#' @param fit A \code{dpmirt_fit} object.
#' @param reference Character. Reference distribution: \code{"normal"}
#'   (default) overlays N(0,1).
#' @param ... Additional graphical parameters.
#'
#' @noRd
.plot_density_compare <- function(fit,
                                  reference = c("normal"),
                                  ...) {
  reference <- match.arg(reference)

  if (is.null(fit$theta_samp)) {
    stop("No theta samples available.", call. = FALSE)
  }

  theta_pm <- colMeans(fit$theta_samp)
  d <- density(theta_pm)

  if (!is.null(fit$dp_density)) {
    dpd <- fit$dp_density
    y_max <- max(c(d$y, dpd$density_mean, dnorm(0))) * 1.1

    plot(d, lwd = 2, col = .dpmirt_colors$parametric,
         main = "Posterior Density Comparison",
         xlab = expression(theta), ylab = "Density",
         ylim = c(0, y_max), ...)
    lines(dpd$grid, dpd$density_mean,
          col = .dpmirt_colors$semiparametric, lwd = 2)
    polygon(c(dpd$grid, rev(dpd$grid)),
            c(dpd$density_lower, rev(dpd$density_upper)),
            col = adjustcolor(.dpmirt_colors$semiparametric, 0.15),
            border = NA)
    curve(dnorm(x), add = TRUE, col = .dpmirt_colors$reference,
          lty = 2, lwd = 1.5)
    legend("topright",
           legend = c("Empirical (PM)", "DP Posterior", "N(0,1)"),
           col = c(.dpmirt_colors$parametric,
                   .dpmirt_colors$semiparametric,
                   .dpmirt_colors$reference),
           lty = c(1, 1, 2), lwd = c(2, 2, 1.5), bty = "n")
  } else {
    y_max <- max(c(d$y, dnorm(0))) * 1.1
    plot(d, lwd = 2, col = .dpmirt_colors$ci_fill,
         main = "Posterior Mean Density vs Reference",
         xlab = expression(theta), ylab = "Density",
         ylim = c(0, y_max), ...)
    curve(dnorm(x), add = TRUE, col = .dpmirt_colors$reference,
          lty = 2, lwd = 1.5)
    legend("topright",
           legend = c("Posterior Mean", "N(0,1)"),
           col = c(.dpmirt_colors$ci_fill, .dpmirt_colors$reference),
           lty = c(1, 2), lwd = c(2, 1.5), bty = "n")
  }
}


#' Test Information Function (base R)
#'
#' @param fit A \code{dpmirt_fit} object.
#' @param theta_range Numeric vector of length 2. Range of theta axis.
#' @param n_points Number of grid points.
#' @param show_items Logical. Show individual item information curves.
#' @param show_density Logical. Overlay person ability density (scaled).
#' @param ... Additional graphical parameters.
#'
#' @noRd
.plot_info <- function(fit, theta_range = c(-4, 4), n_points = 201,
                       show_items = FALSE, show_density = TRUE, ...) {
  if (is.null(fit$beta_samp)) {
    stop("No beta samples available for information plot.", call. = FALSE)
  }

  beta_pm   <- colMeans(fit$beta_samp)
  I <- length(beta_pm)
  lambda_pm <- if (!is.null(fit$lambda_samp)) colMeans(fit$lambda_samp) else
    rep(1, I)
  delta_pm  <- if (!is.null(fit$delta_samp)) colMeans(fit$delta_samp) else
    rep(0, I)

  theta_grid <- seq(theta_range[1], theta_range[2], length.out = n_points)

  # Compute item information
  info_mat <- sapply(seq_len(I), function(j) {
    P <- .icc_prob(theta_grid, beta_pm[j], lambda_pm[j], delta_pm[j])
    P <- pmax(P, delta_pm[j] + 1e-10)  # avoid division by zero
    P <- pmin(P, 1 - 1e-10)

    if (delta_pm[j] > 0) {
      # 3PL information
      lambda_pm[j]^2 * (P - delta_pm[j])^2 /
        ((1 - delta_pm[j])^2 * P * (1 - P))
    } else {
      # Rasch or 2PL information
      lambda_pm[j]^2 * P * (1 - P)
    }
  })

  test_info <- rowSums(info_mat)

  y_max <- max(test_info) * 1.1
  plot(theta_grid, test_info,
       type = "l", lwd = 2.5, col = .dpmirt_colors$ci_fill,
       xlab = expression(theta), ylab = "Information",
       main = paste0("Test Information Function (",
                     toupper(fit$config$model), ")"),
       ylim = c(0, y_max), ...)

  if (show_items) {
    cols <- rep_len(.dpmirt_colors$items, I)
    for (j in seq_len(I)) {
      lines(theta_grid, info_mat[, j], col = adjustcolor(cols[j], 0.5),
            lty = 2)
    }
  }

  if (show_density && !is.null(fit$theta_samp)) {
    theta_pm <- colMeans(fit$theta_samp)
    d <- density(theta_pm)
    # Scale density to fit on same axis
    d_scaled <- d$y / max(d$y) * y_max * 0.3
    polygon(c(d$x, rev(d$x)), c(d_scaled, rep(0, length(d$x))),
            col = adjustcolor(.dpmirt_colors$reference, 0.2), border = NA)
    lines(d$x, d_scaled, col = .dpmirt_colors$reference, lty = 3)
  }
}


#' Posterior predictive check (base R)
#'
#' @param fit A \code{dpmirt_fit} object.
#' @param stat Character. Summary statistic for comparison:
#'   \code{"prop_correct"} (default), \code{"total_score"}.
#' @param n_rep Number of replicated datasets. Default 50.
#' @param ... Additional graphical parameters.
#'
#' @noRd
.plot_pp_check <- function(fit,
                           stat = c("prop_correct", "total_score"),
                           n_rep = 50, ...) {
  stat <- match.arg(stat)

  if (is.null(fit$theta_samp) || is.null(fit$beta_samp)) {
    stop("Need theta and beta samples for PP check.", call. = FALSE)
  }

  niter <- nrow(fit$theta_samp)
  N <- ncol(fit$theta_samp)
  I <- ncol(fit$beta_samp)

  # Select thinned iterations
  rep_iters <- sort(sample.int(niter, min(n_rep, niter)))

  # Compute replicated statistics
  rep_stats <- lapply(rep_iters, function(s) {
    theta_s  <- fit$theta_samp[s, ]
    beta_s   <- fit$beta_samp[s, ]
    lambda_s <- if (!is.null(fit$lambda_samp)) fit$lambda_samp[s, ] else
      rep(1, I)
    delta_s  <- if (!is.null(fit$delta_samp)) fit$delta_samp[s, ] else
      rep(0, I)

    # Compute P matrix
    P_mat <- outer(theta_s, seq_len(I), function(th, j) {
      .icc_prob(th, beta_s[j], lambda_s[j], delta_s[j])
    })

    # Simulate replicated data
    y_rep <- matrix(rbinom(N * I, 1, P_mat), nrow = N, ncol = I)

    if (stat == "prop_correct") {
      colMeans(y_rep)  # per-item proportion correct
    } else {
      rowSums(y_rep)   # total scores
    }
  })

  if (stat == "prop_correct") {
    # Overlay replicated item proportions
    rep_mat <- do.call(rbind, rep_stats)
    rep_mean <- colMeans(rep_mat)
    rep_lower <- apply(rep_mat, 2, quantile, 0.025)
    rep_upper <- apply(rep_mat, 2, quantile, 0.975)

    # Observed: we don't have original data, use theta_pm-based approx
    # Use the posterior mean to compute expected proportions
    theta_pm <- colMeans(fit$theta_samp)
    beta_pm  <- colMeans(fit$beta_samp)
    lambda_pm <- if (!is.null(fit$lambda_samp)) colMeans(fit$lambda_samp) else
      rep(1, I)
    delta_pm <- if (!is.null(fit$delta_samp)) colMeans(fit$delta_samp) else
      rep(0, I)

    expected <- sapply(seq_len(I), function(j) {
      mean(.icc_prob(theta_pm, beta_pm[j], lambda_pm[j], delta_pm[j]))
    })

    plot(seq_len(I), rep_mean, pch = 19, cex = 0.8,
         ylim = c(0, 1),
         xlab = "Item", ylab = "Proportion Correct",
         main = "Posterior Predictive Check: Item Proportions",
         col = .dpmirt_colors$ci_fill, ...)
    segments(seq_len(I), rep_lower, seq_len(I), rep_upper,
             col = adjustcolor(.dpmirt_colors$ci_fill, 0.4), lwd = 2)
    points(seq_len(I), expected, pch = 4, cex = 1, col = "red", lwd = 2)
    legend("topright",
           legend = c("Replicated Mean", "95% PI", "Expected"),
           pch = c(19, NA, 4),
           col = c(.dpmirt_colors$ci_fill, .dpmirt_colors$ci_fill, "red"),
           lty = c(NA, 1, NA), lwd = c(NA, 2, 2),
           bty = "n")
  } else {
    # Total score histograms
    rep_scores <- unlist(rep_stats)
    hist(rep_scores, breaks = 30,
         col = adjustcolor(.dpmirt_colors$ci_fill, 0.4),
         border = "white",
         xlab = "Total Score", main = "PP Check: Total Score Distribution",
         freq = FALSE, ...)
  }
}


# ============================================================================
# S3 Methods for Other Classes
# ============================================================================

#' Plot dpmirt_estimates object
#'
#' @param x A \code{dpmirt_estimates} object.
#' @param type Character. Plot type: \code{"estimates"} (caterpillar of
#'   PM/CB/GR) or \code{"shrinkage"} (PM vs CB scatter).
#' @param param Character. Which parameter: \code{"theta"} or \code{"beta"}.
#' @param ... Additional graphical parameters.
#'
#' @examples
#' \dontrun{
#' sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
#' fit <- dpmirt(sim$response, model = "rasch", prior = "normal",
#'               niter = 5000, nburnin = 1000, seed = 123)
#' est <- dpmirt_estimates(fit)
#'
#' # PM/CB/GR caterpillar plot
#' plot(est, type = "estimates")
#'
#' # Shrinkage plot (PM vs CB)
#' plot(est, type = "shrinkage")
#' }
#'
#' @family visualization
#' @seealso \code{\link{dpmirt_estimates}}, \code{\link{dpmirt_plot_caterpillar}}
#' @export
plot.dpmirt_estimates <- function(x,
                                  type = c("estimates", "shrinkage"),
                                  param = c("theta", "beta"),
                                  ...) {
  type  <- match.arg(type)
  param <- match.arg(param)

  est <- x[[param]]
  if (is.null(est)) {
    stop(sprintf("No %s estimates available.", param), call. = FALSE)
  }

  pm_col <- paste0(param, "_pm")
  if (!pm_col %in% names(est)) {
    stop("Posterior mean column not found.", call. = FALSE)
  }

  pm_vals <- est[[pm_col]]

  if (type == "estimates") {
    # Caterpillar of PM (and CB/GR if available)
    n <- length(pm_vals)
    ord <- order(pm_vals)

    lower_col <- paste0(param, "_lower")
    upper_col <- paste0(param, "_upper")

    plot(pm_vals[ord], seq_len(n),
         pch = 19, cex = 0.5,
         yaxt = "n", xlab = param, ylab = "",
         main = paste0(param, " Estimates"),
         ...)

    if (lower_col %in% names(est)) {
      segments(est[[lower_col]][ord], seq_len(n),
               est[[upper_col]][ord], seq_len(n),
               col = adjustcolor(.dpmirt_colors$ci_fill, 0.3))
    }

    cb_col <- paste0(param, "_cb")
    gr_col <- paste0(param, "_gr")
    if (cb_col %in% names(est)) {
      points(est[[cb_col]][ord], seq_len(n),
             pch = 3, cex = 0.5, col = .dpmirt_colors$semiparametric)
    }
    if (gr_col %in% names(est)) {
      points(est[[gr_col]][ord], seq_len(n),
             pch = 4, cex = 0.5, col = .dpmirt_colors$parametric)
    }

    methods_used <- c("PM")
    cols <- c("black")
    pchs <- c(19)
    if (cb_col %in% names(est)) {
      methods_used <- c(methods_used, "CB")
      cols <- c(cols, .dpmirt_colors$semiparametric)
      pchs <- c(pchs, 3)
    }
    if (gr_col %in% names(est)) {
      methods_used <- c(methods_used, "GR")
      cols <- c(cols, .dpmirt_colors$parametric)
      pchs <- c(pchs, 4)
    }
    legend("bottomright", legend = methods_used,
           col = cols, pch = pchs, bty = "n", cex = 0.8)

    abline(v = 0, col = "gray70", lty = 2)

  } else if (type == "shrinkage") {
    cb_col <- paste0(param, "_cb")
    if (!cb_col %in% names(est)) {
      stop("CB estimates not available for shrinkage plot.", call. = FALSE)
    }

    cb_vals <- est[[cb_col]]
    lims <- range(c(pm_vals, cb_vals))

    plot(pm_vals, cb_vals,
         pch = 19, cex = 0.6,
         col = adjustcolor(.dpmirt_colors$ci_fill, 0.6),
         xlab = "Posterior Mean (PM)",
         ylab = "Constrained Bayes (CB)",
         main = paste0(param, " Shrinkage: PM vs CB"),
         xlim = lims, ylim = lims,
         ...)
    abline(0, 1, col = "gray50", lty = 2)
  }

  invisible(x)
}


#' Plot dpmirt_sim object
#'
#' @param x A \code{dpmirt_sim} object.
#' @param type Character. Plot type: \code{"parameters"} (histograms of true
#'   parameters) or \code{"response"} (heatmap of response matrix).
#' @param ... Additional graphical parameters.
#'
#' @examples
#' \dontrun{
#' sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
#'
#' # True parameter histograms
#' plot(sim, type = "parameters")
#'
#' # Response matrix heatmap
#' plot(sim, type = "response")
#' }
#'
#' @family visualization
#' @seealso \code{\link{dpmirt_simulate}}
#' @export
plot.dpmirt_sim <- function(x,
                            type = c("parameters", "response"),
                            ...) {
  type <- match.arg(type)

  if (type == "parameters") {
    # Count panels needed
    has_lambda <- !is.null(x$lambda) && !all(x$lambda == 1)
    has_delta  <- !is.null(x$delta) && !all(x$delta == 0)
    n_panels <- 2 + has_lambda + has_delta

    old_par <- par(mfrow = c(1, n_panels), mar = c(4, 4, 2, 1))
    on.exit(par(old_par))

    hist(x$theta, breaks = 30, main = expression(theta ~ "(true)"),
         xlab = expression(theta),
         col = adjustcolor(.dpmirt_colors$parametric, 0.5),
         border = "white", ...)
    hist(x$beta, breaks = 20, main = expression(beta ~ "(true)"),
         xlab = expression(beta),
         col = adjustcolor(.dpmirt_colors$semiparametric, 0.5),
         border = "white", ...)
    if (has_lambda) {
      hist(x$lambda, breaks = 20, main = expression(lambda ~ "(true)"),
           xlab = expression(lambda),
           col = adjustcolor("#7570b3", 0.5),
           border = "white", ...)
    }
    if (has_delta) {
      hist(x$delta, breaks = 20, main = expression(delta ~ "(true)"),
           xlab = expression(delta),
           col = adjustcolor("#e7298a", 0.5),
           border = "white", ...)
    }

  } else if (type == "response") {
    resp <- x$response
    image(t(resp[nrow(resp):1, ]),
          col = c("white", .dpmirt_colors$ci_fill),
          xlab = "Item", ylab = "Person",
          main = "Response Matrix",
          axes = FALSE, ...)
    axis(1, at = seq(0, 1, length.out = ncol(resp)),
         labels = seq_len(ncol(resp)), cex.axis = 0.7)
    axis(2, at = seq(0, 1, length.out = nrow(resp)),
         labels = seq_len(nrow(resp)), cex.axis = 0.5, las = 1)
  }

  invisible(x)
}
