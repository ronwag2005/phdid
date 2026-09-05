#' Aggregate cohort-time effects into a reported estimand
#'
#' The cohort-time effects are rarely the final object; papers report an
#' overall ATT, an event-study profile, or effects by cohort. All of these are
#' linear aggregators \eqn{\sum_k w_k \tau_k} of the cell effects, so the
#' specification question is orthogonal to the choice of aggregator.
#'
#' For a [bayes_ph] posterior the aggregator is applied to **each draw** and
#' summarised across the chain. This is the point of the Bayesian route: the
#' reported interval propagates uncertainty about the partition itself, not
#' merely sampling uncertainty given a grouping.
#'
#' For an [l0_ph] or [ph_fit] object the aggregate is computed from the fitted
#' effects with the plug-in variance \eqn{w' R S^{-1} R' w}, which conditions on
#' the selected partition being correct. Remark 3 of the paper is the caveat.
#'
#' @param x a [bayes_ph], [l0_ph] or [ph_fit] object.
#' @param type the aggregation scheme. `"overall"` uses the weights carried by
#'   the [ph_data] object (eq. 2). `"dynamic"` averages within event time
#'   \eqn{t - g}, the event-study profile. `"group"` averages within cohort,
#'   `"calendar"` within calendar period. `"dynamic"`, `"group"` and
#'   `"calendar"` need cohort and time labels on the cells.
#' @param weights optional numeric vector of length K overriding the object's
#'   weights, or a K x J matrix whose columns are J separate aggregators.
#' @param level credible / confidence level.
#' @param ... unused.
#'
#' @return A data frame with one row per reported estimand, holding the
#'   estimate, standard error or posterior standard deviation, and interval
#'   bounds.
#'
#' @references
#' Arora, P. and Wagle, R. (2026). Section 3.3 and Remark 3.
#' See \code{citation("phdid")} for the full reference.
#'
#' @export
#' @examples
#' d <- ph_data(c(0.10, 0.11, 0.42, 0.40), Sigma = 0.02^2 * (0.3 + 0.7 * diag(4)),
#'              cells = c("2:2", "2:3", "3:3", "3:4"))
#' fit <- bayes_ph(d, alpha = 1, iters = 400, burn = 100, seed = 1)
#' aggregate(fit, "overall")
#' aggregate(fit, "dynamic")
aggregate.bayes_ph <- function(x, type = c("overall", "dynamic", "group",
                                           "calendar"),
                               weights = NULL, level = 0.95, ...) {
  W <- aggregator_matrix(x$data, type, weights)
  draws <- x$draws$tau %*% W                       # S x J
  a <- (1 - level) / 2

  data.frame(
    term = colnames(W),
    estimate = colMeans(draws),
    std.error = apply(draws, 2L, stats::sd),
    conf.low = apply(draws, 2L, stats::quantile, probs = a, names = FALSE),
    conf.high = apply(draws, 2L, stats::quantile, probs = 1 - a, names = FALSE),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' @rdname aggregate.bayes_ph
#' @export
aggregate.ph_fit <- function(x, type = c("overall", "dynamic", "group",
                                         "calendar"),
                             weights = NULL, level = 0.95, ...) {
  W <- aggregator_matrix(x$data, type, weights)
  est <- as.numeric(crossprod(W, x$tau))
  V <- crossprod(W, x$vcov %*% W)
  se <- sqrt(pmax(diag(V), 0))
  z <- stats::qnorm(1 - (1 - level) / 2)

  data.frame(
    term = colnames(W),
    estimate = est,
    std.error = se,
    conf.low = est - z * se,
    conf.high = est + z * se,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' Build the K x J matrix of aggregation weights
#'
#' Each column is one reported estimand and sums to one over the cells it
#' touches.
#'
#' @noRd
aggregator_matrix <- function(object, type, weights) {
  K <- object$K

  if (!is.null(weights)) {
    W <- as.matrix(weights)
    if (nrow(W) != K) {
      stop(sprintf("`weights` must have %d rows, one per cell, not %d.",
                   K, nrow(W)), call. = FALSE)
    }
    if (is.null(colnames(W))) {
      colnames(W) <- if (ncol(W) == 1L) "custom" else paste0("w", seq_len(ncol(W)))
    }
    return(W)
  }

  type <- match.arg(type, c("overall", "dynamic", "group", "calendar"))
  if (type == "overall") {
    W <- matrix(object$weights, ncol = 1L)
    colnames(W) <- "overall"
    return(W)
  }

  g <- object$cells$g
  t <- object$cells$t
  if (anyNA(g) || anyNA(t)) {
    stop(sprintf(
      paste0("`type = \"%s\"` needs cohort and calendar labels for every ",
             "cell.\n  Supply them via the `cells` argument of ph_data(), or ",
             "pass explicit `weights`."), type), call. = FALSE)
  }

  key <- switch(type,
    dynamic = t - g,
    group = g,
    calendar = t
  )
  levels_key <- sort(unique(key))
  W <- matrix(0, K, length(levels_key))
  for (j in seq_along(levels_key)) {
    sel <- key == levels_key[j]
    w <- object$weights * sel
    total <- sum(w)
    # If the base weights vanish on this slice, fall back to equal weighting
    # within it rather than returning a column of zeros.
    if (total <= 0) {
      w <- as.numeric(sel)
      total <- sum(w)
    }
    W[, j] <- w / total
  }
  colnames(W) <- paste0(switch(type, dynamic = "e", group = "g", calendar = "t"),
                        "=", levels_key)
  W
}


#' Posterior co-clustering probabilities
#'
#' Extracts the \eqn{K \times K} posterior similarity matrix
#' \eqn{\hat\Pi = [\hat\pi_{jk}]} of eq. (39), the probability that cells
#' \eqn{j} and \eqn{k} are placed in the same group.
#'
#' This is the honest summary of the grouping structure. A single partition,
#' such as the one [l0_ph()] returns, states that certain cells are equal; the
#' co-clustering matrix says how firmly the data support each such statement.
#' A diffuse matrix is not a failure of the method: with few, correlated cells
#' it is the correct report that the fine grouping is genuinely uncertain.
#'
#' @param x a [bayes_ph] object.
#' @return a K x K matrix with cell labels as dimnames.
#' @export
#' @examples
#' d <- ph_data(c(0.10, 0.11, 0.42, 0.40), Sigma = 0.02^2 * (0.3 + 0.7 * diag(4)))
#' fit <- bayes_ph(d, alpha = 1, iters = 500, burn = 100, seed = 1)
#' round(coclustering(fit), 2)
coclustering <- function(x) {
  stopifnot(inherits(x, "bayes_ph"))
  P <- x$coclust
  dimnames(P) <- list(x$data$cells$label, x$data$cells$label)
  P
}


#' A single representative partition from the posterior
#'
#' The posterior mean cohort-time effects average over partitions and so match
#' no single grouping. When one grouping must be reported, it should be chosen
#' by minimising a loss against the whole posterior rather than by taking the
#' most-visited partition, which is unstable. This implements the loss-based
#' summaries the paper points to.
#'
#' Two losses are available. Binder's loss counts disagreeing pairs, weighting
#' each pair by its co-clustering probability. The variation of information
#' (VI) loss of Wade and Ghahramani (2018) is information-theoretic and tends
#' to be less prone to returning too many small clusters.
#'
#' Candidate partitions are the ones actually visited by the sampler, together
#' with every cut of a hierarchical clustering of \eqn{1 - \hat\Pi}, which
#' supplies sensible candidates the chain may not have visited.
#'
#' @param x a [bayes_ph] object.
#' @param loss `"VI"` (default) or `"binder"`.
#' @return an integer vector of group labels, with attributes `loss` and
#'   `n_groups`.
#'
#' @references
#' Wade, S. and Ghahramani, Z. (2018). Bayesian Cluster Analysis: Point
#' Estimation and Credible Balls. \emph{Bayesian Analysis} 13(2), 559--626.
#'
#' Lau, J. W. and Green, P. J. (2007). Bayesian Model-Based Clustering
#' Procedures. \emph{JCGS} 16(3), 526--558.
#'
#' @export
#' @examples
#' d <- ph_data(c(0.10, 0.11, 0.42, 0.40), Sigma = 0.02^2 * (0.3 + 0.7 * diag(4)))
#' fit <- bayes_ph(d, alpha = 1, iters = 500, burn = 100, seed = 1)
#' point_partition(fit)
point_partition <- function(x, loss = c("VI", "binder")) {
  stopifnot(inherits(x, "bayes_ph"))
  loss <- match.arg(loss)
  Pi <- x$coclust
  K <- ncol(Pi)

  candidates <- unique(as.data.frame(t(apply(x$draws$z, 1L, relabel))))
  candidates <- lapply(seq_len(nrow(candidates)), function(i) {
    as.integer(unlist(candidates[i, ], use.names = FALSE))
  })

  # Hierarchical clustering of 1 - Pi supplies well-formed candidates the
  # sampler may never have visited exactly.
  d_pi <- stats::as.dist(1 - Pi)
  hc <- stats::hclust(d_pi, method = "average")
  for (m in seq_len(K)) {
    candidates[[length(candidates) + 1L]] <- relabel(stats::cutree(hc, k = m))
  }

  losses <- vapply(candidates, function(z) {
    if (loss == "binder") binder_loss(z, Pi) else vi_loss(z, Pi)
  }, numeric(1))

  best <- candidates[[which.min(losses)]]
  structure(best, loss = min(losses), n_groups = max(best), loss_type = loss)
}

#' Posterior expected Binder loss of a candidate partition
#'
#' \eqn{\sum_{j<k} [1(z_j = z_k)(1 - \pi_{jk}) + 1(z_j \ne z_k)\pi_{jk}]}.
#'
#' @noRd
binder_loss <- function(z, Pi) {
  same <- outer(z, z, "==")
  ut <- upper.tri(Pi)
  sum(same[ut] * (1 - Pi[ut]) + (!same[ut]) * Pi[ut])
}

#' Wade-Ghahramani lower bound on the posterior expected VI loss
#' @noRd
vi_loss <- function(z, Pi) {
  K <- length(z)
  same <- outer(z, z, "==")
  a <- rowSums(same)
  b <- rowSums(Pi)
  cc <- rowSums(same * Pi)
  sum(log2(a) + log2(b) - 2 * log2(cc)) / K
}


#' @export
confint.bayes_ph <- function(object, parm, level = 0.95, ...) {
  a <- (1 - level) / 2
  qs <- apply(object$draws$tau, 2L, stats::quantile, probs = c(a, 1 - a),
              names = FALSE)
  out <- cbind(qs[1L, ], qs[2L, ])
  dimnames(out) <- list(object$data$cells$label,
                        paste0(round(100 * c(a, 1 - a), 1), "%"))
  if (!missing(parm)) out <- out[parm, , drop = FALSE]
  out
}

#' @export
confint.ph_fit <- function(object, parm, level = 0.95, ...) {
  se <- sqrt(pmax(diag(object$vcov), 0))
  z <- stats::qnorm(1 - (1 - level) / 2)
  a <- (1 - level) / 2
  out <- cbind(object$tau - z * se, object$tau + z * se)
  dimnames(out) <- list(object$data$cells$label,
                        paste0(round(100 * c(a, 1 - a), 1), "%"))
  if (!missing(parm)) out <- out[parm, , drop = FALSE]
  out
}
