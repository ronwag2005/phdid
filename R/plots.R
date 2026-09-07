## Plot methods reproducing the paper's figures. Base graphics only, so the
## package carries no plotting dependency.

#' Plot cohort-time effects grouped by the selected partition
#'
#' Figure 4 (left panel) of the paper: each cohort-time effect with its
#' interval, coloured by the group it was assigned to, with horizontal bars
#' marking the grouped effects. Cells that share a colour and a bar were judged
#' to share a common effect.
#'
#' @param x an [l0_ph] or [ph_fit] object.
#' @param sort order cells by their flexible estimate rather than by cell label.
#' @param level interval level for the flexible estimates.
#' @param col a vector of group colours, recycled as needed.
#' @param main,xlab,ylab plot labels.
#' @param ... passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @export
#' @examples
#' d <- ph_data(c(0.10, 0.11, 0.42, 0.40), Sigma = 0.02^2 * (0.3 + 0.7 * diag(4)))
#' plot(l0_ph(d))
plot.ph_fit <- function(x, sort = TRUE, level = 0.95,
                        col = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a",
                                "#66a61e", "#e6ab02", "#a6761d"),
                        main = NULL, xlab = "cohort-time cell",
                        ylab = "effect", ...) {
  tau <- x$data$tau
  se <- sqrt(diag(x$data$Sigma))
  z <- stats::qnorm(1 - (1 - level) / 2)
  K <- length(tau)
  ord <- if (sort) order(tau) else seq_len(K)
  cols <- rep(col, length.out = max(x$partition))

  if (is.null(main)) {
    main <- sprintf("Cohort-time effects grouped into %d level%s",
                    x$m, if (x$m == 1L) "" else "s")
  }

  graphics::plot(seq_len(K), tau[ord], pch = 19, cex = 1.2,
                 col = cols[x$partition[ord]],
                 ylim = range(tau - z * se, tau + z * se),
                 xaxt = "n", xlab = xlab, ylab = ylab, main = main, ...)
  graphics::axis(1, at = seq_len(K), labels = x$data$cells$label[ord],
                 las = 2, cex.axis = 0.75)
  graphics::arrows(seq_len(K), (tau - z * se)[ord],
                   seq_len(K), (tau + z * se)[ord],
                   length = 0, col = "grey70")
  graphics::segments(seq_len(K) - 0.3, x$tau[ord],
                     seq_len(K) + 0.3, x$tau[ord],
                     col = cols[x$partition[ord]], lwd = 3)
  graphics::abline(h = 0, lty = 3)
  graphics::legend("bottomright", paste("group", seq_len(x$m)),
                   col = cols[seq_len(x$m)], pch = 19, bty = "n")
  invisible(x)
}


#' Plot the posterior co-clustering matrix
#'
#' Figure 4 (right panel): darker cells co-cluster more often. Read this rather
#' than a single partition when the grouping is uncertain --- a diffuse matrix
#' is the correct report that the data do not pin the fine structure down.
#'
#' @param x a [bayes_ph] or `ph_enumeration` object.
#' @param sort order cells by their flexible estimate.
#' @param main plot title.
#' @param digits if not `NULL`, overlay the probabilities rounded to this many
#'   digits.
#' @param ... passed to [graphics::image()].
#' @return `x`, invisibly.
#' @export
#' @examples
#' d <- ph_data(c(0.10, 0.11, 0.42, 0.40), Sigma = 0.02^2 * (0.3 + 0.7 * diag(4)))
#' plot_coclustering(bayes_ph(d, iters = 500, burn = 100, seed = 1))
plot_coclustering <- function(x, sort = TRUE,
                              main = "Posterior co-clustering probability",
                              digits = NULL, ...) {
  P <- x$coclust
  labels <- x$data$cells$label
  K <- ncol(P)
  ord <- if (sort) order(x$data$tau) else seq_len(K)
  P <- P[ord, ord]

  op <- graphics::par(mar = c(6, 6, 3, 1))
  on.exit(graphics::par(op), add = TRUE)

  graphics::image(seq_len(K), seq_len(K), P, axes = FALSE, xlab = "", ylab = "",
                  zlim = c(0, 1),
                  col = grDevices::grey.colors(24, start = 1, end = 0),
                  main = main, ...)
  graphics::axis(1, seq_len(K), labels[ord], las = 2, cex.axis = 0.75)
  graphics::axis(2, seq_len(K), labels[ord], las = 1, cex.axis = 0.75)
  graphics::box()
  if (!is.null(digits)) {
    for (i in seq_len(K)) for (j in seq_len(K)) {
      graphics::text(i, j, formatC(P[i, j], format = "f", digits = digits),
                     cex = 0.6, col = if (P[i, j] > 0.5) "white" else "black")
    }
  }
  invisible(x)
}

#' @rdname plot_coclustering
#' @export
plot.bayes_ph <- function(x, ...) plot_coclustering(x, ...)


#' Plot a regularisation path
#'
#' Draws the table returned by [alpha_sensitivity()] or [lambda_sensitivity()].
#' The two views answer different questions.
#'
#' With an aggregate table (`type = "overall"` and friends) the top panel shows
#' the estimate and its band against the fully pooled and fully flexible
#' anchors, and the lower panel the number of groups. Where the path is flat the
#' tuning parameter does not matter; where it slides, report the path rather
#' than a single value.
#'
#' With a per-cell table (`type = "cells"`) each cohort-time effect gets its own
#' line, drawn from its flexible value on the loose-penalty side toward the
#' pooled value as the penalty tightens. Lines that converge are cells the
#' method is merging, and where they meet is the penalty at which it does so.
#' This is the view that shows what a flat aggregate path can hide: the overall
#' ATT is robust to over-pooling, while the individual effects are not.
#'
#' @param x a data frame from [alpha_sensitivity()] or [lambda_sensitivity()].
#' @param term for an aggregate table with several terms, which one to plot.
#' @param band draw the interval band. On a per-cell plot with many cells the
#'   bands overlap badly, so the default omits them there.
#' @param label write the cell name at the right edge of each line.
#' @param col a vector of line colours, recycled across cells.
#' @param main plot title.
#' @param ... passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @export
#' @examples
#' d <- ph_data(c(0.10, 0.11, 0.42, 0.40), Sigma = 0.02^2 * (0.3 + 0.7 * diag(4)))
#' plot_sensitivity(lambda_sensitivity(d, type = "cells"))
plot_sensitivity <- function(x, term = NULL, band = NULL, label = TRUE,
                             col = c("#1b9e77", "#d95f02", "#7570b3",
                                     "#e7298a", "#66a61e", "#e6ab02",
                                     "#a6761d"),
                             main = NULL, ...) {
  param <- attr(x, "param")
  if (is.null(param)) {
    stop("`x` does not look like a sensitivity table; use ",
         "alpha_sensitivity() or lambda_sensitivity().", call. = FALSE)
  }
  if (identical(attr(x, "type"), "cells")) {
    plot_sensitivity_cells(x, param, band, label, col, main, ...)
  } else {
    plot_sensitivity_aggregate(x, param, term, main, ...)
  }
}

#' @noRd
plot_sensitivity_aggregate <- function(x, param, term, main, ...) {
  if (!is.null(term)) x <- x[x$term == term, , drop = FALSE]
  if (length(unique(x$term)) > 1L) {
    stop("The table holds several terms; pick one with `term`.", call. = FALSE)
  }
  anchors <- attr(x, "anchors")
  tune <- x[[param]]
  logx <- if (all(tune > 0)) "x" else ""

  op <- graphics::par(mfrow = c(2, 1), mar = c(4, 4, 2.5, 1))
  on.exit(graphics::par(op), add = TRUE)

  if (is.null(main)) {
    main <- sprintf("Sensitivity to %s", param)
  }
  ylim <- range(x$conf.low, x$conf.high, anchors, na.rm = TRUE)
  graphics::plot(tune, x$estimate, type = "n", log = logx, ylim = ylim,
                 xlab = param, ylab = "estimate", main = main, ...)
  graphics::polygon(c(tune, rev(tune)), c(x$conf.low, rev(x$conf.high)),
                    col = grDevices::adjustcolor("steelblue", 0.2),
                    border = NA)
  graphics::lines(tune, x$estimate, lwd = 2, col = "steelblue")
  if (!is.null(anchors)) {
    graphics::abline(h = anchors[["pooled"]], lty = 3, col = "firebrick")
    graphics::abline(h = anchors[["flexible"]], lty = 2, col = "navy")
    graphics::legend("topright", c("estimate", "pooled", "flexible"),
                     col = c("steelblue", "firebrick", "navy"),
                     lty = c(1, 3, 2), lwd = c(2, 1, 1), bty = "n")
  }

  graphics::plot(tune, x$m, type = "b", pch = 19, log = logx,
                 xlab = param, ylab = "groups", main = "")
  if (!is.null(x$prior_m)) {
    graphics::lines(tune, x$prior_m, lty = 2, col = "grey50")
    graphics::legend("topleft", c("posterior", "prior E[m]"),
                     col = c("black", "grey50"), lty = c(1, 2),
                     pch = c(19, NA), bty = "n")
  }
  invisible(x)
}

#' @noRd
plot_sensitivity_cells <- function(x, param, band, label, col, main, ...) {
  cells <- unique(x$term)
  K <- length(cells)
  if (is.null(band)) band <- K <= 6L
  cols <- rep(col, length.out = K)
  tune <- sort(unique(x[[param]]))
  logx <- if (all(tune > 0)) "x" else ""

  op <- graphics::par(mar = c(4, 4, 2.5, if (label) 6 else 1))
  on.exit(graphics::par(op), add = TRUE)

  if (is.null(main)) {
    main <- sprintf("Cohort-time effects across %s", param)
  }
  ylim <- range(if (band) c(x$conf.low, x$conf.high) else x$estimate,
                attr(x, "flexible"), attr(x, "pooled"), na.rm = TRUE)

  graphics::plot(range(tune), ylim, type = "n", log = logx,
                 xlab = param, ylab = "effect", main = main, ...)
  graphics::abline(h = attr(x, "pooled"), lty = 3, col = "firebrick")
  graphics::abline(h = 0, lty = 3, col = "grey60")

  for (i in seq_len(K)) {
    s <- x[x$term == cells[i], , drop = FALSE]
    s <- s[order(s[[param]]), ]
    if (band) {
      graphics::polygon(c(s[[param]], rev(s[[param]])),
                        c(s$conf.low, rev(s$conf.high)),
                        col = grDevices::adjustcolor(cols[i], 0.13),
                        border = NA)
    }
    graphics::lines(s[[param]], s$estimate, lwd = 2, col = cols[i])
    if (label) {
      graphics::text(max(tune), s$estimate[nrow(s)], paste0(" ", cells[i]),
                     col = cols[i], cex = 0.7, adj = 0, xpd = NA)
    }
  }
  graphics::legend("bottomleft", "fully pooled", col = "firebrick", lty = 3,
                   bty = "n", cex = 0.8)
  invisible(x)
}

#' Plot the l0 solution path
#'
#' Figure 3 (left): how the fit and the criteria move along the agglomeration
#' path, indexed by the number of groups. Small `m` is a strong penalty. The
#' dashed line marks the selected point.
#'
#' @param x an [l0_ph] object.
#' @param ... unused.
#' @return `x`, invisibly.
#' @export
#' @examples
#' d <- ph_data(c(0.10, 0.11, 0.42, 0.40), Sigma = 0.02^2 * (0.3 + 0.7 * diag(4)))
#' plot_l0_path(l0_ph(d))
plot_l0_path <- function(x, ...) {
  stopifnot(inherits(x, "l0_ph"))
  p <- x$path

  op <- graphics::par(mar = c(4, 4, 3, 4))
  on.exit(graphics::par(op), add = TRUE)

  graphics::plot(p$m, p$bic, type = "b", pch = 19, col = "firebrick",
                 xlab = "number of groups m  (small m = strong penalty)",
                 ylab = "BIC", main = "l0 agglomeration path")
  graphics::abline(v = x$selected, lty = 2)
  graphics::par(new = TRUE)
  graphics::plot(p$m, p$deviance, type = "b", pch = 17, col = "steelblue",
                 axes = FALSE, xlab = "", ylab = "", log = "y")
  graphics::axis(4)
  graphics::mtext("GLS deviance (log scale)", side = 4, line = 2.5)
  graphics::legend("topright", c("BIC (left)", "deviance (right)", "selected"),
                   col = c("firebrick", "steelblue", "black"),
                   pch = c(19, 17, NA), lty = c(1, 1, 2), bty = "n")
  invisible(x)
}


#' Plot event effects against a placebo noise band
#'
#' Figure 7: the cohort-time effects sorted, with the pooled effect as a dashed
#' line and a shaded band showing the spread expected under a common effect.
#' When the estimates stay inside the band there is no recoverable
#' heterogeneity, and any partition returned by the estimators is a partition
#' of noise.
#'
#' The band is the central `level` interval of the reference distribution: from
#' the supplied placebo estimates when [homogeneity_test()] was given them, and
#' otherwise from the first-stage covariance under the common-effect null.
#'
#' @param x a [homogeneity_test] object.
#' @param level width of the band.
#' @param main,xlab,ylab plot labels.
#' @param ... passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @export
#' @examples
#' d <- ph_data(c(0.10, 0.11, 0.42, 0.40), Sigma = 0.02^2 * (0.3 + 0.7 * diag(4)))
#' plot_placebo_band(homogeneity_test(d))
plot_placebo_band <- function(x, level = 0.95,
                              main = "Cell effects against the noise band",
                              xlab = "cell (sorted by effect)",
                              ylab = "effect", ...) {
  stopifnot(inherits(x, "homogeneity_test"))
  d <- x$data
  tau <- d$tau
  K <- length(tau)
  ord <- order(tau)
  pooled <- as.numeric(crossprod(d$weights, pooled_twfe(d)$tau))

  # Spread expected around the common effect, per cell, under the null.
  noise_sd <- if (!is.null(x$placebo)) {
    rep(x$placebo$s_pre, K)
  } else {
    sqrt(pmax(diag(d$Sigma), 0))
  }
  z <- stats::qnorm(1 - (1 - level) / 2)

  graphics::plot(seq_len(K), tau[ord], pch = 19, cex = 0.9,
                 ylim = range(tau, pooled + z * noise_sd,
                              pooled - z * noise_sd),
                 xlab = xlab, ylab = ylab, main = main, ...)
  graphics::polygon(
    c(seq_len(K), rev(seq_len(K))),
    c(pooled - z * noise_sd[ord], rev(pooled + z * noise_sd[ord])),
    col = grDevices::adjustcolor("grey60", 0.3), border = NA)
  graphics::points(seq_len(K), tau[ord], pch = 19, cex = 0.9)
  graphics::abline(h = pooled, lty = 2)
  graphics::abline(h = 0, lty = 3, col = "grey50")
  graphics::legend("topleft",
                   c("cell effect", "pooled effect",
                     sprintf("%.0f%% band under a common effect",
                             100 * level)),
                   pch = c(19, NA, 15), lty = c(NA, 2, NA),
                   col = c("black", "black", grDevices::adjustcolor("grey60", 0.5)),
                   bty = "n")
  invisible(x)
}


#' Plot Monte Carlo results across the true number of effects
#'
#' Figure 1: average sampling variance of the cohort-time estimates against
#' \eqn{m^*}. The feasible estimators track the infeasible oracle in the
#' partial-homogeneity region and carry a selection overhead when every cell is
#' genuinely distinct.
#'
#' @param x a data frame of [sim_study()] results, stacked across `m_star`.
#' @param y_var which column to plot.
#' @param ... unused.
#' @return `x`, invisibly.
#' @export
#' @examples
#' \donttest{
#' des <- ph_design(N = 300, T = 6, cohorts = c(0, 3, 5))
#' res <- do.call(rbind, lapply(c(1, 3), function(m) {
#'   sim_study(des, m_star = m, delta = 6, R = 3,
#'             methods = c("flexible", "oracle"), progress = FALSE, seed = m)
#' }))
#' plot_sim_study(res)
#' }
plot_sim_study <- function(x, y_var = "var_ratio", ...) {
  if (!y_var %in% names(x)) {
    stop(sprintf("`%s` is not a column of the results table.", y_var),
         call. = FALSE)
  }
  methods <- unique(x$method)
  cols <- grDevices::hcl.colors(length(methods), "Dark 3")
  names(cols) <- methods

  graphics::plot(range(x$m_star), range(x[[y_var]], na.rm = TRUE), type = "n",
                 xlab = expression(m^"*" ~ "(true number of distinct effects)"),
                 ylab = y_var,
                 main = sprintf("Monte Carlo: %s", y_var))
  if (y_var == "var_ratio") graphics::abline(h = 1, lty = 3)
  for (mth in methods) {
    s <- x[x$method == mth, ]
    s <- s[order(s$m_star), ]
    graphics::lines(s$m_star, s[[y_var]], type = "b", pch = 19, lwd = 2,
                    col = cols[mth])
  }
  graphics::legend("topleft", methods, col = cols, lwd = 2, pch = 19,
                   bty = "n")
  invisible(x)
}
