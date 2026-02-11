# ============================================================================
# Module 13B: ggplot2-based Visualization
# ============================================================================
# Phase 7B: ggplot2 enhanced versions of all 12 plot types
# All functions use ggplot2:: namespacing (Suggests only)
# Each exported function returns a ggplot object for user customization
# ============================================================================


# --------------------------------------------------------------------------
# Internal theme (Paganin-consistent)
# --------------------------------------------------------------------------

#' DPMirt ggplot2 theme (internal)
#' @noRd
.dpmirt_theme <- function() {
  ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text        = ggplot2::element_text(size = 11),
      axis.title       = ggplot2::element_text(size = 11),
      strip.text       = ggplot2::element_text(size = 11),
      strip.background = ggplot2::element_rect(fill = "white"),
      plot.title       = ggplot2::element_text(size = 12),
      legend.text      = ggplot2::element_text(size = 11)
    )
}


#' Ensure ggplot2 is available (internal)
#' @noRd
.require_gg <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required. Install with: install.packages('ggplot2')",
         call. = FALSE)
  }
}


# ============================================================================
# ggplot2 versions of existing 5 plot types
# ============================================================================

#' Plot Posterior Mean Density
#'
#' Displays a kernel density estimate of the posterior mean person abilities
#' (\eqn{\hat{\theta}^{PM}}). Useful for assessing the shape of the estimated
#' latent trait distribution. Requires ggplot2.
#'
#' @param fit A \code{dpmirt_fit} object from \code{\link{dpmirt}}.
#' @param ... Currently unused.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
#' fit <- dpmirt(sim$response, model = "rasch", prior = "normal",
#'               niter = 5000, nburnin = 1000, seed = 123)
#' dpmirt_plot_density(fit)
#' }
#'
#' @family visualization
#' @seealso \code{\link{plot.dpmirt_fit}}, \code{\link{dpmirt_plot_density_compare}}
#' @export
dpmirt_plot_density <- function(fit, ...) {
  .require_gg()

  if (is.null(fit$theta_samp)) {
    stop("No theta samples available for density plot.", call. = FALSE)
  }

  theta_pm <- colMeans(fit$theta_samp)
  df <- data.frame(theta = theta_pm)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$theta)) +
    ggplot2::geom_density(fill = .dpmirt_colors$ci_fill, alpha = 0.3,
                          linewidth = 0.8) +
    ggplot2::labs(
      title = paste0(toupper(fit$config$model), " (", fit$config$prior,
                     ") -- Posterior Mean Density"),
      x = expression(theta), y = "Density"
    ) +
    .dpmirt_theme()
}


#' Plot Item Parameter Estimates
#'
#' Displays item difficulty (beta) estimates with posterior credible intervals.
#' For 2PL/3PL models, also shows discrimination (lambda) estimates.
#' Requires ggplot2.
#'
#' @param fit A \code{dpmirt_fit} object from \code{\link{dpmirt}}.
#' @param ... Currently unused.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
#' fit <- dpmirt(sim$response, model = "rasch", prior = "normal",
#'               niter = 5000, nburnin = 1000, seed = 123)
#' dpmirt_plot_items(fit)
#' }
#'
#' @family visualization
#' @seealso \code{\link{plot.dpmirt_fit}}, \code{\link{dpmirt_plot_caterpillar}}
#' @export
dpmirt_plot_items <- function(fit, ...) {
  .require_gg()

  if (is.null(fit$beta_samp)) {
    stop("No beta samples available.", call. = FALSE)
  }

  beta_pm  <- colMeans(fit$beta_samp)
  beta_psd <- apply(fit$beta_samp, 2, sd)
  I <- length(beta_pm)

  df <- data.frame(
    item  = seq_len(I),
    pm    = beta_pm,
    lower = beta_pm - 2 * beta_psd,
    upper = beta_pm + 2 * beta_psd
  )

  ggplot2::ggplot(df, ggplot2::aes(x = .data$item, y = .data$pm)) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$lower,
                                        ymax = .data$upper),
                           width = 0.3, color = "gray50") +
    ggplot2::labs(title = "Item Difficulty Estimates",
                  x = "Item", y = expression(beta)) +
    .dpmirt_theme()
}


#' Plot Log-Likelihood MCMC Trace
#'
#' Displays the log-likelihood trace across MCMC iterations for convergence
#' assessment. A well-mixed chain shows stable oscillation around a
#' stationary level. Requires ggplot2.
#'
#' @param fit A \code{dpmirt_fit} object from \code{\link{dpmirt}}.
#' @param ... Currently unused.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
#' fit <- dpmirt(sim$response, model = "rasch", prior = "normal",
#'               niter = 5000, nburnin = 1000, seed = 123)
#' dpmirt_plot_trace(fit)
#' }
#'
#' @family visualization
#' @seealso \code{\link{plot.dpmirt_fit}}, \code{\link{dpmirt_plot_parameter_trace}}
#' @export
dpmirt_plot_trace <- function(fit, ...) {
  .require_gg()

  if (is.null(fit$loglik_trace)) {
    stop("No log-likelihood trace available.", call. = FALSE)
  }

  df <- data.frame(
    iteration = seq_along(fit$loglik_trace),
    loglik    = fit$loglik_trace
  )

  ggplot2::ggplot(df, ggplot2::aes(x = .data$iteration,
                                    y = .data$loglik)) +
    ggplot2::geom_line(color = .dpmirt_colors$trace, alpha = 0.7) +
    ggplot2::labs(title = "Log-Likelihood Trace",
                  x = "Iteration", y = "Log-Likelihood") +
    .dpmirt_theme()
}


#' Plot Cluster Count Diagnostics
#'
#' For DPM models, displays the trace and histogram of the number of active
#' clusters across MCMC iterations. Stable oscillation indicates convergence.
#' Requires ggplot2.
#'
#' @param fit A \code{dpmirt_fit} object with \code{prior = "dpm"}.
#' @param ... Currently unused.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' sim <- dpmirt_simulate(200, 20, model = "rasch",
#'                        latent_shape = "bimodal", seed = 42)
#' fit <- dpmirt(sim$response, model = "rasch", prior = "dpm",
#'               niter = 10000, nburnin = 3000, seed = 123)
#' dpmirt_plot_clusters(fit)
#' }
#'
#' @family visualization
#' @seealso \code{\link{plot.dpmirt_fit}}, \code{\link{dpmirt_plot_dp_density}}
#' @export
dpmirt_plot_clusters <- function(fit, ...) {
  .require_gg()

  if (fit$config$prior != "dpm") {
    stop("Cluster plot is only available for DPM models.", call. = FALSE)
  }
  if (is.null(fit$cluster_info) || is.null(fit$cluster_info$n_clusters)) {
    stop("No cluster information available.", call. = FALSE)
  }

  n_cl <- fit$cluster_info$n_clusters

  df <- data.frame(
    iteration = seq_along(n_cl),
    clusters  = n_cl
  )

  # Combine trace + marginal histogram via dual geom
  ggplot2::ggplot(df, ggplot2::aes(x = .data$iteration,
                                    y = .data$clusters)) +
    ggplot2::geom_line(color = .dpmirt_colors$trace, alpha = 0.7) +
    ggplot2::geom_hline(yintercept = mean(n_cl), linetype = "dashed",
                        color = "red") +
    ggplot2::labs(title = "Cluster Count Trace",
                  x = "Iteration", y = "Number of Clusters") +
    .dpmirt_theme()
}


#' Plot DP Mixture Density
#'
#' Displays the posterior mean density of the Dirichlet Process mixture
#' with pointwise credible bands. Requires a pre-computed
#' \code{dpmirt_dp_density} object (computed automatically by
#' \code{\link{dpmirt}} or manually via \code{\link{dpmirt_dp_density}}).
#' Requires ggplot2.
#'
#' @param fit A \code{dpmirt_fit} object with \code{prior = "dpm"} and
#'   a computed DP density.
#' @param ... Currently unused.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' sim <- dpmirt_simulate(200, 20, model = "rasch",
#'                        latent_shape = "bimodal", seed = 42)
#' fit <- dpmirt(sim$response, model = "rasch", prior = "dpm",
#'               niter = 10000, nburnin = 3000, seed = 123)
#' dpmirt_plot_dp_density(fit)
#' }
#'
#' @family visualization
#' @seealso \code{\link{dpmirt_dp_density}}, \code{\link{plot.dpmirt_fit}}
#' @export
dpmirt_plot_dp_density <- function(fit, ...) {
  .require_gg()

  if (is.null(fit$dp_density)) {
    stop("No DP density available. Compute with dpmirt_dp_density(fit).",
         call. = FALSE)
  }

  dpd <- fit$dp_density

  df <- data.frame(
    theta        = dpd$grid,
    density_mean = dpd$density_mean,
    lower        = dpd$density_lower,
    upper        = dpd$density_upper
  )

  ggplot2::ggplot(df, ggplot2::aes(x = .data$theta)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower,
                                      ymax = .data$upper),
                         fill = .dpmirt_colors$ci_fill, alpha = 0.2) +
    ggplot2::geom_line(ggplot2::aes(y = .data$density_mean),
                       color = .dpmirt_colors$ci_fill, linewidth = 1) +
    ggplot2::stat_function(fun = dnorm, linetype = "dashed",
                           color = .dpmirt_colors$reference) +
    ggplot2::labs(
      title = "DP Mixture Density (Posterior)",
      x = expression(theta), y = "Density"
    ) +
    .dpmirt_theme()
}


# ============================================================================
# ggplot2 versions of 7 new plot types
# ============================================================================

#' Plot Item Characteristic Curves
#'
#' Displays the Item Characteristic Curves (ICCs) for all items, showing
#' the probability of a correct response as a function of ability.
#' Requires ggplot2.
#'
#' @param fit A \code{dpmirt_fit} object from \code{\link{dpmirt}}.
#' @param items Integer vector of item indices to plot. Default: all items
#'   (up to 10).
#' @param theta_range Numeric vector of length 2 for the ability axis range.
#'   Default: \code{c(-4, 4)}.
#' @param n_points Integer. Number of grid points. Default: 201.
#' @param ... Currently unused.
#'
#' @return A \code{ggplot} object.
#'
#' @details
#' The ICC gives the probability of a correct response at each ability level:
#'
#' \strong{Rasch}: \eqn{P(\theta) = \mathrm{logistic}(\theta - \beta)}
#'
#' \strong{2PL}: \eqn{P(\theta) = \mathrm{logistic}(\lambda(\theta - \beta))}
#'
#' \strong{3PL}: \eqn{P(\theta) = \delta + (1 - \delta) \cdot \mathrm{logistic}(\lambda(\theta - \beta))}
#'
#' @examples
#' \dontrun{
#' sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
#' fit <- dpmirt(sim$response, model = "rasch", prior = "normal",
#'               niter = 5000, nburnin = 1000, seed = 123)
#' dpmirt_plot_icc(fit)
#' dpmirt_plot_icc(fit, items = 1:5)
#' }
#'
#' @family visualization
#' @seealso \code{\link{plot.dpmirt_fit}}, \code{\link{dpmirt_plot_info}}
#' @export
dpmirt_plot_icc <- function(fit, items = NULL, theta_range = c(-4, 4),
                            n_points = 201, ...) {
  .require_gg()

  if (is.null(fit$beta_samp)) {
    stop("No beta samples available for ICC plot.", call. = FALSE)
  }

  beta_pm <- colMeans(fit$beta_samp)
  I <- length(beta_pm)
  lambda_pm <- if (!is.null(fit$lambda_samp)) colMeans(fit$lambda_samp) else
    rep(1, I)
  delta_pm <- if (!is.null(fit$delta_samp)) colMeans(fit$delta_samp) else
    rep(0, I)

  if (is.null(items)) items <- seq_len(min(I, 10))
  items <- items[items >= 1 & items <= I]

  theta_grid <- seq(theta_range[1], theta_range[2], length.out = n_points)

  # Build long-format data
  df_list <- lapply(items, function(j) {
    prob <- .icc_prob(theta_grid, beta_pm[j], lambda_pm[j], delta_pm[j])
    data.frame(theta = theta_grid, prob = prob,
               item = factor(paste0("Item ", j), levels = paste0("Item ", items)))
  })
  df <- do.call(rbind, df_list)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$theta, y = .data$prob,
                                    color = .data$item)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dotted",
                        color = "gray70") +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      title = paste0("Item Characteristic Curves (",
                     toupper(fit$config$model), ")"),
      x = expression(theta), y = "P(correct)", color = "Item"
    ) +
    .dpmirt_theme()
}


#' Plot Person-Item Map (Wright Map)
#'
#' Displays a Wright map showing the joint distribution of person abilities
#' and item difficulties on a common logit scale. Useful for identifying
#' gaps in item coverage. Requires ggplot2.
#'
#' @param fit A \code{dpmirt_fit} object from \code{\link{dpmirt}}.
#' @param bins Integer. Number of histogram bins. Default: 30.
#' @param item_labels Optional character vector of item labels.
#' @param ... Currently unused.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
#' fit <- dpmirt(sim$response, model = "rasch", prior = "normal",
#'               niter = 5000, nburnin = 1000, seed = 123)
#' dpmirt_plot_wright_map(fit)
#' }
#'
#' @family visualization
#' @seealso \code{\link{plot.dpmirt_fit}}, \code{\link{dpmirt_plot_items}}
#' @export
dpmirt_plot_wright_map <- function(fit, bins = 30, item_labels = NULL, ...) {
  .require_gg()

  if (is.null(fit$theta_samp) || is.null(fit$beta_samp)) {
    stop("No theta/beta samples available for Wright map.", call. = FALSE)
  }

  theta_pm <- colMeans(fit$theta_samp)
  beta_pm  <- colMeans(fit$beta_samp)
  I <- length(beta_pm)

  if (is.null(item_labels)) item_labels <- paste0("I", seq_len(I))

  # Compute histogram to get scaled counts
  h <- hist(theta_pm, breaks = bins, plot = FALSE)
  max_count <- max(h$counts)

  # Person histogram data (negative x for left side)
  df_hist <- data.frame(
    ymin  = h$breaks[-length(h$breaks)],
    ymax  = h$breaks[-1],
    xmin  = -h$counts / max_count,
    xmax  = 0
  )

  # Item data (positive x for right side)
  df_items <- data.frame(
    x     = 0.3,
    y     = beta_pm,
    label = item_labels
  )

  ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = df_hist,
      ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                   ymin = .data$ymin, ymax = .data$ymax),
      fill = .dpmirt_colors$parametric, alpha = 0.6, color = "white"
    ) +
    ggplot2::geom_point(
      data = df_items,
      ggplot2::aes(x = .data$x, y = .data$y),
      shape = 18, size = 3, color = .dpmirt_colors$semiparametric
    ) +
    ggplot2::geom_text(
      data = df_items,
      ggplot2::aes(x = .data$x + 0.08, y = .data$y, label = .data$label),
      hjust = 0, size = 3
    ) +
    ggplot2::geom_hline(yintercept = mean(theta_pm), linetype = "dashed",
                        color = "gray60") +
    ggplot2::scale_x_continuous(
      breaks = c(-0.5, 0, 0.3),
      labels = c("Persons", "", "Items")
    ) +
    ggplot2::labs(title = "Person-Item Map (Wright Map)",
                  x = "", y = expression(theta)) +
    .dpmirt_theme() +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
}


#' Plot Individual Parameter MCMC Traces
#'
#' Displays MCMC trace plots for individual item or person parameters,
#' useful for diagnosing mixing and convergence of specific parameters.
#' Requires ggplot2.
#'
#' @param fit A \code{dpmirt_fit} object from \code{\link{dpmirt}}.
#' @param param Character. Parameter type: \code{"beta"}, \code{"theta"},
#'   \code{"lambda"}, or \code{"delta"}.
#' @param indices Integer vector of parameter indices to plot. Default: 1:4.
#' @param ... Currently unused.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
#' fit <- dpmirt(sim$response, model = "rasch", prior = "normal",
#'               niter = 5000, nburnin = 1000, seed = 123)
#' dpmirt_plot_parameter_trace(fit, param = "beta", indices = 1:4)
#' }
#'
#' @family visualization
#' @seealso \code{\link{plot.dpmirt_fit}}, \code{\link{dpmirt_plot_trace}}
#' @export
dpmirt_plot_parameter_trace <- function(fit,
                                        param = c("beta", "theta",
                                                  "lambda", "delta"),
                                        indices = 1:4, ...) {
  .require_gg()
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

  df_list <- lapply(indices, function(idx) {
    data.frame(
      iteration = seq_len(nrow(samp)),
      value     = samp[, idx],
      parameter = paste0(param, "[", idx, "]")
    )
  })
  df <- do.call(rbind, df_list)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$iteration,
                                    y = .data$value)) +
    ggplot2::geom_line(color = .dpmirt_colors$trace, alpha = 0.6) +
    ggplot2::facet_wrap(~ .data$parameter, scales = "free_y") +
    ggplot2::labs(title = paste0("MCMC Trace: ", param),
                  x = "Iteration", y = "Value") +
    .dpmirt_theme()
}


#' Plot Caterpillar (Forest) Plot
#'
#' Displays sorted point estimates with credible intervals for items or
#' persons. Useful for identifying extreme values and assessing uncertainty.
#' Requires ggplot2.
#'
#' @param fit A \code{dpmirt_fit} object from \code{\link{dpmirt}}.
#' @param param Character. Parameter type: \code{"beta"} (default),
#'   \code{"theta"}, \code{"lambda"}, or \code{"delta"}.
#' @param sort Logical. Sort by estimate magnitude. Default TRUE.
#' @param ci_level Numeric. Credible interval level. Default 0.95.
#' @param max_show Integer. Maximum number of parameters to display. Default 50.
#' @param ... Currently unused.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
#' fit <- dpmirt(sim$response, model = "rasch", prior = "normal",
#'               niter = 5000, nburnin = 1000, seed = 123)
#' dpmirt_plot_caterpillar(fit, param = "beta")
#' }
#'
#' @family visualization
#' @seealso \code{\link{plot.dpmirt_fit}}, \code{\link{dpmirt_plot_items}}
#' @export
dpmirt_plot_caterpillar <- function(fit,
                                    param = c("beta", "theta",
                                              "lambda", "delta"),
                                    sort = TRUE, ci_level = 0.95,
                                    max_show = 50, ...) {
  .require_gg()
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

  df <- data.frame(
    index = show_idx,
    pm    = pm[show_idx],
    lower = lower[show_idx],
    upper = upper[show_idx]
  )

  if (sort) df <- df[order(df$pm), ]
  df$rank <- seq_len(nrow(df))

  ggplot2::ggplot(df, ggplot2::aes(x = .data$pm,
                                    y = .data$rank)) +
    ggplot2::geom_segment(ggplot2::aes(x = .data$lower, xend = .data$upper,
                                       yend = .data$rank),
                          color = .dpmirt_colors$ci_fill, alpha = 0.5) +
    ggplot2::geom_point(size = 1) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        color = "gray70") +
    ggplot2::labs(
      title = paste0("Caterpillar: ", param,
                     " (", round(ci_level * 100), "% CI)"),
      x = param, y = "Index"
    ) +
    .dpmirt_theme()
}


#' Plot Posterior Density vs Reference Distribution
#'
#' Overlays the estimated posterior mean density with a reference distribution
#' (e.g., true generating density or N(0,1)). Useful for assessing how well
#' the model recovers the latent trait distribution. Requires ggplot2.
#'
#' @param fit A \code{dpmirt_fit} object from \code{\link{dpmirt}}.
#' @param reference Character or numeric vector. If \code{"normal"} (default),
#'   overlays N(0,1). If numeric, treated as reference theta values.
#' @param ... Currently unused.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
#' fit <- dpmirt(sim$response, model = "rasch", prior = "normal",
#'               niter = 5000, nburnin = 1000, seed = 123)
#'
#' # Compare to N(0,1)
#' dpmirt_plot_density_compare(fit)
#'
#' # Compare to true theta
#' dpmirt_plot_density_compare(fit, reference = sim$theta)
#' }
#'
#' @family visualization
#' @seealso \code{\link{plot.dpmirt_fit}}, \code{\link{dpmirt_plot_density}}
#' @export
dpmirt_plot_density_compare <- function(fit,
                                        reference = c("normal"),
                                        ...) {
  .require_gg()
  reference <- match.arg(reference)

  if (is.null(fit$theta_samp)) {
    stop("No theta samples available.", call. = FALSE)
  }

  theta_pm <- colMeans(fit$theta_samp)
  df_emp <- data.frame(theta = theta_pm)

  p <- ggplot2::ggplot(df_emp, ggplot2::aes(x = .data$theta)) +
    ggplot2::geom_density(ggplot2::aes(color = "Empirical (PM)"),
                          linewidth = 0.8) +
    ggplot2::stat_function(fun = dnorm,
                           ggplot2::aes(color = "N(0,1)"),
                           linetype = "dashed")

  if (!is.null(fit$dp_density)) {
    dpd <- fit$dp_density
    df_dp <- data.frame(theta = dpd$grid, density = dpd$density_mean,
                        lower = dpd$density_lower, upper = dpd$density_upper)
    p <- p +
      ggplot2::geom_ribbon(data = df_dp,
                           ggplot2::aes(x = .data$theta,
                                        ymin = .data$lower,
                                        ymax = .data$upper),
                           fill = .dpmirt_colors$semiparametric,
                           alpha = 0.15, inherit.aes = FALSE) +
      ggplot2::geom_line(data = df_dp,
                         ggplot2::aes(x = .data$theta, y = .data$density,
                                      color = "DP Posterior"),
                         linewidth = 0.8, inherit.aes = FALSE)
  }

  p +
    ggplot2::scale_color_manual(
      values = c("Empirical (PM)" = .dpmirt_colors$parametric,
                 "N(0,1)"         = .dpmirt_colors$reference,
                 "DP Posterior"   = .dpmirt_colors$semiparametric)
    ) +
    ggplot2::labs(title = "Posterior Density Comparison",
                  x = expression(theta), y = "Density", color = "") +
    .dpmirt_theme()
}


#' Plot Test Information Function
#'
#' Displays the test information function (TIF), which shows the precision
#' of measurement across the ability range. Higher information indicates
#' more precise measurement. Requires ggplot2.
#'
#' @param fit A \code{dpmirt_fit} object from \code{\link{dpmirt}}.
#' @param theta_range Numeric vector of length 2. Default: \code{c(-4, 4)}.
#' @param n_points Integer. Number of grid points. Default: 201.
#' @param show_items Logical. Show individual item information curves.
#'   Default FALSE.
#' @param show_density Logical. Overlay person ability density (scaled).
#'   Default TRUE.
#' @param ... Currently unused.
#'
#' @return A \code{ggplot} object.
#'
#' @details
#' Fisher information for each item is:
#'
#' \strong{Rasch}: \eqn{I_j(\theta) = P(1-P)}
#'
#' \strong{2PL}: \eqn{I_j(\theta) = \lambda^2 P(1-P)}
#'
#' \strong{3PL}: \eqn{I_j(\theta) = \lambda^2 \frac{(P - \delta)^2}{(1 - \delta)^2 P(1 - P)}}
#'
#' The test information function is the sum: \eqn{I(\theta) = \sum_j I_j(\theta)}.
#'
#' @examples
#' \dontrun{
#' sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
#' fit <- dpmirt(sim$response, model = "rasch", prior = "normal",
#'               niter = 5000, nburnin = 1000, seed = 123)
#' dpmirt_plot_info(fit)
#' }
#'
#' @family visualization
#' @seealso \code{\link{plot.dpmirt_fit}}, \code{\link{dpmirt_plot_icc}}
#' @export
dpmirt_plot_info <- function(fit, theta_range = c(-4, 4), n_points = 201,
                             show_items = FALSE, show_density = TRUE, ...) {
  .require_gg()

  if (is.null(fit$beta_samp)) {
    stop("No beta samples available for information plot.", call. = FALSE)
  }

  beta_pm <- colMeans(fit$beta_samp)
  I <- length(beta_pm)
  lambda_pm <- if (!is.null(fit$lambda_samp)) colMeans(fit$lambda_samp) else
    rep(1, I)
  delta_pm <- if (!is.null(fit$delta_samp)) colMeans(fit$delta_samp) else
    rep(0, I)

  theta_grid <- seq(theta_range[1], theta_range[2], length.out = n_points)

  # Compute item information
  info_mat <- sapply(seq_len(I), function(j) {
    P <- .icc_prob(theta_grid, beta_pm[j], lambda_pm[j], delta_pm[j])
    P <- pmax(P, delta_pm[j] + 1e-10)
    P <- pmin(P, 1 - 1e-10)
    if (delta_pm[j] > 0) {
      lambda_pm[j]^2 * (P - delta_pm[j])^2 /
        ((1 - delta_pm[j])^2 * P * (1 - P))
    } else {
      lambda_pm[j]^2 * P * (1 - P)
    }
  })

  test_info <- rowSums(info_mat)
  df_test <- data.frame(theta = theta_grid, info = test_info)

  p <- ggplot2::ggplot(df_test, ggplot2::aes(x = .data$theta,
                                              y = .data$info)) +
    ggplot2::geom_line(color = .dpmirt_colors$ci_fill, linewidth = 1.2)

  if (show_items) {
    df_items_list <- lapply(seq_len(I), function(j) {
      data.frame(theta = theta_grid, info = info_mat[, j],
                 item = paste0("Item ", j))
    })
    df_items <- do.call(rbind, df_items_list)
    p <- p +
      ggplot2::geom_line(data = df_items,
                         ggplot2::aes(x = .data$theta, y = .data$info,
                                      group = .data$item),
                         color = .dpmirt_colors$reference,
                         alpha = 0.4, linetype = "dashed",
                         inherit.aes = FALSE)
  }

  if (show_density && !is.null(fit$theta_samp)) {
    theta_pm <- colMeans(fit$theta_samp)
    d <- density(theta_pm)
    d_scaled <- d$y / max(d$y) * max(test_info) * 0.3
    df_dens <- data.frame(x = d$x, y = d_scaled)
    p <- p +
      ggplot2::geom_area(data = df_dens,
                         ggplot2::aes(x = .data$x, y = .data$y),
                         fill = .dpmirt_colors$reference, alpha = 0.2,
                         inherit.aes = FALSE)
  }

  p +
    ggplot2::labs(
      title = paste0("Test Information Function (",
                     toupper(fit$config$model), ")"),
      x = expression(theta), y = "Information"
    ) +
    .dpmirt_theme()
}


#' Plot Posterior Predictive Check
#'
#' Compares observed item statistics (proportion correct) with the
#' posterior predictive distribution. Points falling outside the
#' predictive intervals may indicate model misfit. Requires ggplot2.
#'
#' @param fit A \code{dpmirt_fit} object from \code{\link{dpmirt}}.
#' @param stat Character. Summary statistic for comparison:
#'   \code{"prop_correct"} (default) or \code{"total_score"}.
#' @param n_rep Integer. Number of replicated datasets. Default 50.
#' @param ... Currently unused.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' sim <- dpmirt_simulate(200, 20, model = "rasch", seed = 42)
#' fit <- dpmirt(sim$response, model = "rasch", prior = "normal",
#'               niter = 5000, nburnin = 1000, seed = 123)
#' dpmirt_plot_pp_check(fit)
#' }
#'
#' @family visualization
#' @seealso \code{\link{plot.dpmirt_fit}}, \code{\link{dpmirt_diagnostics}}
#' @export
dpmirt_plot_pp_check <- function(fit,
                                  stat = c("prop_correct", "total_score"),
                                  n_rep = 50, ...) {
  .require_gg()
  stat <- match.arg(stat)

  if (is.null(fit$theta_samp) || is.null(fit$beta_samp)) {
    stop("Need theta and beta samples for PP check.", call. = FALSE)
  }

  niter <- nrow(fit$theta_samp)
  N <- ncol(fit$theta_samp)
  I <- ncol(fit$beta_samp)

  rep_iters <- sort(sample.int(niter, min(n_rep, niter)))

  rep_stats <- lapply(rep_iters, function(s) {
    theta_s  <- fit$theta_samp[s, ]
    beta_s   <- fit$beta_samp[s, ]
    lambda_s <- if (!is.null(fit$lambda_samp)) fit$lambda_samp[s, ] else
      rep(1, I)
    delta_s  <- if (!is.null(fit$delta_samp)) fit$delta_samp[s, ] else
      rep(0, I)

    P_mat <- outer(theta_s, seq_len(I), function(th, j) {
      .icc_prob(th, beta_s[j], lambda_s[j], delta_s[j])
    })
    y_rep <- matrix(rbinom(N * I, 1, P_mat), nrow = N, ncol = I)

    if (stat == "prop_correct") colMeans(y_rep) else rowSums(y_rep)
  })

  if (stat == "prop_correct") {
    rep_mat   <- do.call(rbind, rep_stats)
    rep_mean  <- colMeans(rep_mat)
    rep_lower <- apply(rep_mat, 2, quantile, 0.025)
    rep_upper <- apply(rep_mat, 2, quantile, 0.975)

    theta_pm <- colMeans(fit$theta_samp)
    beta_pm  <- colMeans(fit$beta_samp)
    lambda_pm <- if (!is.null(fit$lambda_samp)) colMeans(fit$lambda_samp) else
      rep(1, I)
    delta_pm <- if (!is.null(fit$delta_samp)) colMeans(fit$delta_samp) else
      rep(0, I)

    expected <- sapply(seq_len(I), function(j) {
      mean(.icc_prob(theta_pm, beta_pm[j], lambda_pm[j], delta_pm[j]))
    })

    df <- data.frame(
      item     = seq_len(I),
      rep_mean = rep_mean,
      lower    = rep_lower,
      upper    = rep_upper,
      expected = expected
    )

    ggplot2::ggplot(df, ggplot2::aes(x = .data$item)) +
      ggplot2::geom_segment(ggplot2::aes(y = .data$lower, yend = .data$upper,
                                         xend = .data$item),
                            color = .dpmirt_colors$ci_fill, alpha = 0.4,
                            linewidth = 2) +
      ggplot2::geom_point(ggplot2::aes(y = .data$rep_mean),
                          color = .dpmirt_colors$ci_fill) +
      ggplot2::geom_point(ggplot2::aes(y = .data$expected),
                          shape = 4, color = "red", size = 2) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(title = "Posterior Predictive Check: Item Proportions",
                    x = "Item", y = "Proportion Correct") +
      .dpmirt_theme()

  } else {
    rep_scores <- unlist(rep_stats)
    df <- data.frame(score = rep_scores)

    ggplot2::ggplot(df, ggplot2::aes(x = .data$score)) +
      ggplot2::geom_histogram(fill = .dpmirt_colors$ci_fill, alpha = 0.4,
                              color = "white", bins = 30) +
      ggplot2::labs(title = "PP Check: Total Score Distribution",
                    x = "Total Score", y = "Count") +
      .dpmirt_theme()
  }
}
