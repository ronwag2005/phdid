#' Convergence diagnostics for the Gibbs sampler
#'
#' Reports the Gelman--Rubin statistic \eqn{\hat R} and the effective sample
#' size for the two scalar summaries that matter: the overall aggregate effect
#' and the number of groups. This is the check of Appendix E, where four chains
#' from dispersed initialisations (all singletons, one pooled group, and two
#' random partitions) return \eqn{\hat R = 1.00} for both quantities.
#'
#' \eqn{\hat R} needs at least two chains, so run [bayes_ph()] with
#' `chains >= 2`. With a single chain only the effective sample sizes are
#' reported.
#'
#' @param x a [bayes_ph] object.
#' @param split split each chain in half before computing \eqn{\hat R}, which
#'   detects within-chain trends that the classic statistic can miss.
#' @return An object of class `ph_rhat`; a data frame of diagnostics is stored
#'   in its `stats` component.
#'
#' @references
#' Arora, P. and Wagle, R. (2026). Appendix E.
#' See \code{citation("phdid")} for the full reference.
#'
#' @export
#' @examples
#' d <- ph_data(c(0.10, 0.11, 0.42, 0.40), Sigma = 0.02^2 * (0.3 + 0.7 * diag(4)))
#' fit <- bayes_ph(d, alpha = 1, iters = 1000, burn = 200, chains = 4, seed = 1)
#' ph_rhat(fit)
ph_rhat <- function(x, split = TRUE) {
  stopifnot(inherits(x, "bayes_ph"))
  chain <- x$draws$chain
  agg <- as.numeric(x$draws$tau %*% x$data$weights)
  quantities <- list(`overall effect` = agg, `number of groups` = x$draws$m)

  stats_df <- do.call(rbind, lapply(names(quantities), function(nm) {
    v <- quantities[[nm]]
    mats <- split(v, chain)
    n_chains <- length(mats)
    data.frame(
      quantity = nm,
      n_chains = n_chains,
      rhat = if (n_chains >= 2L) gelman_rubin(mats, split = split) else NA_real_,
      ess = ess_sum(mats),
      mean = mean(v),
      sd = stats::sd(v),
      stringsAsFactors = FALSE
    )
  }))

  structure(list(stats = stats_df, chains = x$chains, iters = x$iters,
                 burn = x$burn),
            class = "ph_rhat")
}

#' Gelman-Rubin potential scale reduction factor
#' @noRd
gelman_rubin <- function(chains, split = TRUE) {
  if (split) {
    chains <- unlist(lapply(chains, function(v) {
      h <- floor(length(v) / 2)
      if (h < 2L) return(list(v))
      list(v[seq_len(h)], v[seq(h + 1L, 2 * h)])
    }), recursive = FALSE)
  }
  n <- min(lengths(chains))
  if (n < 2L) return(NA_real_)
  chains <- lapply(chains, function(v) v[seq_len(n)])
  m <- length(chains)

  means <- vapply(chains, mean, numeric(1))
  vars <- vapply(chains, stats::var, numeric(1))
  W <- mean(vars)
  B <- n * stats::var(means)
  if (W <= 0) return(NA_real_)
  var_hat <- ((n - 1) / n) * W + B / n
  sqrt(var_hat / W)
}

#' Effective sample size, summed across chains
#'
#' Uses the initial positive sequence estimator: autocorrelations are summed
#' until the first pair of consecutive lags whose sum is negative.
#'
#' @noRd
ess_sum <- function(chains) {
  sum(vapply(chains, ess_one, numeric(1)))
}

#' @noRd
ess_one <- function(v) {
  n <- length(v)
  if (n < 10L || stats::sd(v) == 0) return(as.numeric(n))
  max_lag <- min(n - 2L, 500L)
  ac <- stats::acf(v, lag.max = max_lag, plot = FALSE, demean = TRUE)$acf[-1L]
  total <- 0
  for (i in seq(1L, length(ac) - 1L, by = 2L)) {
    pair <- ac[i] + ac[i + 1L]
    if (pair < 0) break
    total <- total + pair
  }
  as.numeric(min(n, n / (1 + 2 * total)))
}

#' @export
print.ph_rhat <- function(x, ...) {
  cat("<ph_rhat>  sampler convergence diagnostics\n\n")
  cat(sprintf("  %d chain(s), %d sweeps each, %d burn-in\n\n",
              x$chains, x$iters, x$burn))
  df <- x$stats
  df$rhat <- ifelse(is.na(df$rhat), "  --", sprintf("%.3f", df$rhat))
  df$ess <- sprintf("%.0f", df$ess)
  df$mean <- sprintf("%.4f", df$mean)
  df$sd <- sprintf("%.4f", df$sd)
  print(df[, c("quantity", "rhat", "ess", "mean", "sd")], row.names = FALSE)
  if (x$chains < 2L) {
    cat("\n  R-hat needs at least two chains; refit with `chains = 4` for the\n",
        "  dispersed-initialisation check of Appendix E.\n", sep = "")
  } else {
    worst <- max(x$stats$rhat, na.rm = TRUE)
    note <- if (worst < 1.01) {
      "  (converged)"
    } else {
      "  (above 1.01: the chains have not mixed; run longer)"
    }
    cat(sprintf("\n  Largest R-hat: %.3f%s\n", worst, note))
  }
  invisible(x)
}


#' Regularisation paths: how the estimates move with the tuning parameter
#'
#' Both estimators have a single knob that slides them between full pooling and
#' the fully flexible fit --- the Dirichlet Process concentration `alpha` for
#' [bayes_ph()], and the \eqn{\ell_0} penalty `lambda` for [l0_ph()]. These
#' functions trace what the knob does.
#'
#' Reporting the path rather than a single value is the paper's own
#' recommendation, and for a good reason. In its simulation the path is nearly
#' flat, so the choice is immaterial; but in its seven-cell application the
#' overall effect slides from about -0.017 under heavy pooling to about -0.038
#' as the model approaches the flexible fit. Choosing the tuning parameter from
#' the data would turn it into a data-dependent object, with the usual
#' empirical-Bayes consequences of understated posterior uncertainty and a
#' double use of the data, so the honest report is the whole path.
#'
#' @section Aggregate or per-cell:
#' With `type = "overall"` (the default) each row is one value of the tuning
#' parameter and the reported estimand is the overall ATT. `"dynamic"`,
#' `"group"` and `"calendar"` behave the same way for the corresponding
#' aggregations.
#'
#' With `type = "cells"` the path is reported **for every cohort-time effect
#' separately**, one row per (tuning value, cell). This is usually the more
#' informative view: the aggregate is famously robust to over-pooling, so a
#' flat overall path can hide cells that move a great deal. Watching the
#' individual CATTs collapse toward one another as the penalty tightens shows
#' exactly which cells the grouping is merging, and at what point.
#'
#' @param object a [ph_data] object.
#' @param alpha_grid concentration values to trace.
#' @param lambda_grid penalty values to trace. `NULL` derives a grid from the
#'   agglomeration path so that every partition the penalty can reach is
#'   visited.
#' @param by index the \eqn{\ell_0} path by the penalty (`"lambda"`, the
#'   default) or by the number of groups (`"m"`). See the note below on why
#'   these are not the same path.
#' @param type `"overall"` (default), `"cells"` for one row per cell, or any
#'   aggregation accepted by [aggregate.bayes_ph()].
#' @param level credible / confidence level.
#' @param n_bic passed to [l0_ph()].
#' @param ... further arguments to [bayes_ph()], such as `iters` and `burn`.
#' @param seed optional integer for reproducibility.
#'
#' @return A data frame with one row per tuning value (and per cell when
#'   `type = "cells"`), holding the estimate and its interval bounds. The
#'   `anchors` attribute carries the fully pooled and fully flexible reference
#'   values, and `param` names the tuning parameter. Plot with
#'   [plot_sensitivity()].
#'
#' @seealso [plot_sensitivity()].
#'
#' @export
#' @examples
#' d <- ph_data(c(0.10, 0.11, 0.42, 0.40), Sigma = 0.02^2 * (0.3 + 0.7 * diag(4)))
#'
#' # The overall ATT across the concentration parameter.
#' alpha_sensitivity(d, alpha_grid = c(0.5, 1, 5), iters = 400, burn = 100,
#'                   seed = 1)
#'
#' # The same path for every cohort-time effect separately.
#' alpha_sensitivity(d, alpha_grid = c(0.5, 5), type = "cells",
#'                   iters = 400, burn = 100, seed = 1)
alpha_sensitivity <- function(object, alpha_grid = c(0.1, 0.25, 0.5, 1, 2, 5,
                                                     10, 25, 50, 100),
                              type = "overall", level = 0.95, ...,
                              seed = NULL) {
  stopifnot(inherits(object, "ph_data"))
  if (!is.null(seed)) set.seed(seed)

  rows <- lapply(alpha_grid, function(a) {
    fit <- bayes_ph(object, alpha = a, progress = FALSE, ...)
    body <- if (identical(type, "cells")) {
      ci <- stats::confint(fit, level = level)
      data.frame(term = object$cells$label,
                 estimate = fit$tau,
                 conf.low = ci[, 1L], conf.high = ci[, 2L],
                 stringsAsFactors = FALSE)
    } else {
      agg <- aggregate(fit, type = type, level = level)
      data.frame(term = agg$term, estimate = agg$estimate,
                 conf.low = agg$conf.low, conf.high = agg$conf.high,
                 stringsAsFactors = FALSE)
    }
    data.frame(alpha = a, body,
               m = fit$m_mean,
               prior_m = expected_groups(object$K, a),
               stringsAsFactors = FALSE, row.names = NULL)
  })
  finish_sensitivity(do.call(rbind, rows), object, type, "alpha")
}


#' @rdname alpha_sensitivity
#'
#' @details
#' `lambda_sensitivity()` is the frequentist counterpart. For each penalty it
#' refits [l0_ph()] under the Appendix B stopping rule and records the grouped
#' effects. Because the \eqn{\ell_0} path is a step function --- the selected
#' partition changes only when `lambda` crosses a merge threshold --- the
#' default grid is derived from the agglomeration path itself, so that every
#' distinct partition on it is visited exactly once rather than being missed
#' between grid points.
#'
#' Its intervals are the plug-in intervals of [l0_ph()] and condition on the
#' selected partition being correct, so they narrow as the penalty pools more
#' cells even where that pooling is wrong. Remark 3 of the paper is the caveat;
#' the Bayesian path does not have this problem because it averages over
#' partitions.
#'
#' @section Some group counts are unreachable by any lambda:
#' Indexing by `lambda` and indexing by `m` are genuinely different paths, and
#' the first can skip values the second visits. The Appendix B stopping rule
#' merges while a pair's cost falls below \eqn{\lambda |A| |B|}, and the
#' factor \eqn{|A||B|} grows as groups absorb one another. The effective
#' per-merge threshold \eqn{\mathrm{cost}/(|A||B|)} is therefore not monotone
#' in the merge order, so a later merge can be cheaper in threshold terms than
#' an earlier one --- and when that happens, raising `lambda` past the earlier
#' merge triggers both at once and no value of `lambda` selects the partition
#' in between.
#'
#' On the paper's first application this is not hypothetical: the thresholds
#' run 0.079, 0.686, 0.397, 1.350, 2.765, 3.571, so nothing between the second
#' and third selects five groups. Use `by = "m"` to walk the agglomeration path
#' one merge at a time, which visits every group count from `K` down to 1.
#'
#' @export
#' @examples
#'
#' # The l0 path, per cell, across the penalty.
#' lambda_sensitivity(d, type = "cells")
lambda_sensitivity <- function(object, lambda_grid = NULL, type = "overall",
                               by = c("lambda", "m"), level = 0.95,
                               n_bic = NULL, ...) {
  stopifnot(inherits(object, "ph_data"))
  by <- match.arg(by)

  base <- l0_ph(object, n_bic = n_bic)
  grid <- if (identical(by, "m")) {
    seq_len(object$K)
  } else if (is.null(lambda_grid)) {
    default_lambda_grid(base)
  } else {
    lambda_grid
  }

  rows <- lapply(grid, function(lam) {
    fit <- if (identical(by, "m")) {
      l0_ph(object, select = "m", m = lam, n_bic = n_bic)
    } else {
      l0_ph(object, select = "lambda", lambda = lam, n_bic = n_bic)
    }
    body <- if (identical(type, "cells")) {
      ci <- stats::confint(fit, level = level)
      data.frame(term = object$cells$label,
                 estimate = fit$tau,
                 conf.low = ci[, 1L], conf.high = ci[, 2L],
                 stringsAsFactors = FALSE)
    } else {
      agg <- aggregate(fit, type = type, level = level)
      data.frame(term = agg$term, estimate = agg$estimate,
                 conf.low = agg$conf.low, conf.high = agg$conf.high,
                 stringsAsFactors = FALSE)
    }
    if (identical(by, "m")) {
      data.frame(m = fit$m, body, stringsAsFactors = FALSE, row.names = NULL)
    } else {
      data.frame(lambda = lam, body, m = fit$m,
                 stringsAsFactors = FALSE, row.names = NULL)
    }
  })
  out <- finish_sensitivity(do.call(rbind, rows), object, type, by)
  attr(out, "selected") <- base$m
  out
}

#' Attach the anchors and metadata every sensitivity table carries
#' @noRd
finish_sensitivity <- function(out, object, type, param) {
  anchors <- if (identical(type, "cells")) {
    NULL
  } else {
    c(pooled = as.numeric(crossprod(object$weights, pooled_twfe(object)$tau)),
      flexible = as.numeric(crossprod(object$weights, object$tau)))
  }
  attr(out, "anchors") <- anchors
  attr(out, "param") <- param
  attr(out, "type") <- type
  attr(out, "flexible") <- stats::setNames(object$tau, object$cells$label)
  attr(out, "pooled") <- pooled_twfe(object)$tau[1L]
  out
}

#' A lambda grid that visits every partition on the agglomeration path
#'
#' The path is a step function in `lambda`: the selected partition changes only
#' when the penalty crosses a merge threshold, which for the merge producing
#' `m` groups is `merge_cost / (|A| * |B|)`. Placing a grid point just above
#' each threshold, plus one below the smallest, visits every distinct partition
#' exactly once.
#'
#' @noRd
default_lambda_grid <- function(fit) {
  thresh <- fit$path$merge_cost / fit$path$pair_product
  if (is.null(fit$path$pair_product)) {
    thresh <- fit$path$merge_cost
  }
  thresh <- sort(unique(thresh[is.finite(thresh) & thresh > 0]))
  if (!length(thresh)) return(c(0, 1))
  # Just above each threshold, plus a point below the smallest (no merges).
  eps <- 1 + 1e-6
  c(thresh[1L] / 2, thresh * eps)
}


#' Compare exact and diagonal treatment of the covariance
#'
#' Runs the sampler twice, once with the exact collapsed marginal likelihood in
#' the cluster-assignment moves and once with the diagonal (independence)
#' shortcut, and reports how far apart they land. Appendix E finds the two agree
#' on aggregates but that the shortcut can misstate individual co-clustering
#' probabilities by up to 0.39 in a correlated design.
#'
#' Worth running once on any new design: it tells you whether the shortcut,
#' which is much faster, is safe for the quantities you intend to report.
#'
#' @param object a [ph_data] object.
#' @param ... further arguments to [bayes_ph()].
#' @param seed optional integer for reproducibility.
#' @return A list with the two fits and a summary of the discrepancies.
#' @export
#' @examples
#' d <- ph_data(c(0.10, 0.11, 0.42, 0.40), Sigma = 0.02^2 * (0.3 + 0.7 * diag(4)))
#' cmp <- covariance_check(d, iters = 500, burn = 100, seed = 1)
#' cmp$summary
covariance_check <- function(object, ..., seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  exact <- bayes_ph(object, marginal = "exact", progress = FALSE, ...)
  if (!is.null(seed)) set.seed(seed)
  diag_fit <- bayes_ph(object, marginal = "diagonal", progress = FALSE, ...)

  a_exact <- aggregate(exact, "overall")
  a_diag <- aggregate(diag_fit, "overall")

  summary <- data.frame(
    quantity = c("overall effect", "overall CI width",
                 "E[# groups]", "max co-clustering gap"),
    exact = c(a_exact$estimate, a_exact$conf.high - a_exact$conf.low,
              exact$m_mean, NA_real_),
    diagonal = c(a_diag$estimate, a_diag$conf.high - a_diag$conf.low,
                 diag_fit$m_mean, NA_real_),
    gap = c(abs(a_exact$estimate - a_diag$estimate),
            abs((a_exact$conf.high - a_exact$conf.low) -
                  (a_diag$conf.high - a_diag$conf.low)),
            abs(exact$m_mean - diag_fit$m_mean),
            max(abs(exact$coclust - diag_fit$coclust))),
    stringsAsFactors = FALSE
  )

  list(exact = exact, diagonal = diag_fit, summary = summary)
}
