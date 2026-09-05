#' Exact partition posterior by enumeration
#'
#' When \eqn{K} is small enough, the posterior over partitions can be computed
#' in closed form rather than sampled: evaluate the collapsed marginal
#' likelihood and the Chinese Restaurant Process prior at every one of the
#' \eqn{B_K} partitions and normalise. This is the benchmark of Appendix E,
#' where the paper's first application has \eqn{K = 7} and hence \eqn{B_7 =
#' 877} partitions, and the sampler is shown to agree with the exact answer to
#' within Monte Carlo error.
#'
#' Use it for two things: to validate a sampler run on a small design, and, on
#' such designs, simply to report the exact posterior instead of an approximate
#' one. The Bell numbers grow faster than exponentially, so this is infeasible
#' beyond about eleven cells.
#'
#' @param object a [ph_data] object.
#' @param alpha,mu0,sigma0_sq prior settings, matching [bayes_ph()].
#' @param max_cells refuse to enumerate beyond this many cells. Raise it
#'   deliberately; \eqn{B_{12}} is over four million.
#'
#' @return An object of class `ph_enumeration` with the exact posterior mean
#'   effects, credible intervals, co-clustering matrix, expected number of
#'   groups, and the full table of partitions with their posterior
#'   probabilities.
#'
#' @references
#' Arora, P. and Wagle, R. (2026). Appendix E and eq. (25).
#' See \code{citation("phdid")} for the full reference.
#'
#' @export
#' @examples
#' d <- ph_data(c(0.10, 0.11, 0.42, 0.40), Sigma = 0.02^2 * (0.3 + 0.7 * diag(4)))
#' ex <- enumerate_partitions(d, alpha = 1)
#' ex
#'
#' # The sampler should agree with this.
#' fit <- bayes_ph(d, alpha = 1, iters = 2000, burn = 500, seed = 1)
#' max(abs(coclustering(fit) - ex$coclust))
enumerate_partitions <- function(object, alpha = 1, mu0 = NULL,
                                 sigma0_sq = NULL, max_cells = 11L) {
  stopifnot(inherits(object, "ph_data"))
  K <- object$K
  if (K > max_cells) {
    stop(sprintf(
      paste0("Exact enumeration needs all %s partitions of %d cells, which is ",
             "beyond the\n  `max_cells` limit of %d. Use bayes_ph() and check ",
             "convergence with ph_rhat()."),
      format(bell_number(K), big.mark = ",", scientific = FALSE), K,
      max_cells), call. = FALSE)
  }
  if (is.null(mu0)) mu0 <- stats::median(object$tau)
  if (is.null(sigma0_sq)) sigma0_sq <- (10 * stats::sd(object$tau))^2

  Z <- all_partitions(K)
  n_part <- nrow(Z)
  Otau <- as.numeric(object$Omega %*% object$tau)

  log_post <- numeric(n_part)
  m_vec <- integer(n_part)
  # Posterior of the aggregate is a mixture of normals, one per partition.
  agg_mean <- numeric(n_part)
  agg_sd <- numeric(n_part)
  tau_mean <- matrix(NA_real_, n_part, K)
  w <- object$weights

  for (i in seq_len(n_part)) {
    z <- Z[i, ]
    m_vec[i] <- max(z)
    log_post[i] <- log_marginal(z, object$Omega, Otau, mu0, sigma0_sq) +
      log_crp_prior(z, alpha)

    R <- indicator_matrix(z)
    S <- group_precision(object$Omega, z)
    u <- group_score(Otau, z)
    P <- S + diag(ncol(S)) / sigma0_sq
    Sstar <- safe_inverse(P, "group-effect posterior precision")
    mu_star <- as.numeric(Sstar %*% (u + mu0 / sigma0_sq))

    tau_mean[i, ] <- as.numeric(R %*% mu_star)
    agg_mean[i] <- sum(w * tau_mean[i, ])
    Rw <- crossprod(R, w)
    agg_sd[i] <- sqrt(max(0, as.numeric(crossprod(Rw, Sstar %*% Rw))))
  }

  prob <- exp(log_post - max(log_post))
  prob <- prob / sum(prob)

  coclust <- matrix(0, K, K)
  for (i in seq_len(n_part)) {
    if (prob[i] < 1e-12) next
    z <- Z[i, ]
    coclust <- coclust + prob[i] * outer(z, z, `==`)
  }

  post_tau <- as.numeric(crossprod(prob, tau_mean))
  agg <- mixture_summary(prob, agg_mean, agg_sd)

  structure(
    list(
      tau = post_tau,
      coclust = coclust,
      m_mean = sum(prob * m_vec),
      m_table = tapply(prob, m_vec, sum),
      overall = agg,
      partitions = Z,
      prob = prob,
      n_partitions = n_part,
      alpha = alpha, mu0 = mu0, sigma0_sq = sigma0_sq,
      data = object
    ),
    class = "ph_enumeration"
  )
}

#' Log Chinese Restaurant Process prior on a partition
#'
#' Eq. (25): \eqn{\Pr(\mathcal{P}) = \alpha^m \prod_l (|C_l| - 1)! /
#' \alpha^{(K)}} with \eqn{\alpha^{(K)} = \alpha(\alpha+1)\cdots(\alpha+K-1)}.
#' The normalising constant does not depend on the partition, but it is cheap
#' and keeping it makes the returned log-posterior comparable across `alpha`.
#'
#' @noRd
log_crp_prior <- function(z, alpha) {
  sizes <- tabulate(relabel(z))
  m <- length(sizes)
  K <- length(z)
  m * log(alpha) + sum(lgamma(sizes)) -
    sum(log(alpha + seq_len(K) - 1))
}

#' All set partitions of 1..K as restricted growth strings
#'
#' Row `i` is a canonical label vector: `z[1] == 1` and `z[j] <= max(z[1:(j-1)])
#' + 1`, which enumerates each partition exactly once.
#'
#' @noRd
all_partitions <- function(K) {
  n <- bell_number(K)
  out <- matrix(NA_integer_, n, K)
  z <- integer(K)
  z[1L] <- 1L
  row <- 0L

  recurse <- function(i, mx) {
    if (i > K) {
      row <<- row + 1L
      out[row, ] <<- z
      return(invisible(NULL))
    }
    for (v in seq_len(mx + 1L)) {
      z[i] <<- v
      recurse(i + 1L, max(mx, v))
    }
    invisible(NULL)
  }
  if (K == 1L) return(matrix(1L, 1L, 1L))
  recurse(2L, 1L)
  out
}

#' Bell number B_K, the count of partitions of K elements
#' @noRd
bell_number <- function(K) {
  row <- 1
  for (i in seq_len(K)) {
    new <- numeric(i + 1L)
    new[1L] <- row[length(row)]
    for (j in seq_len(i)) new[j + 1L] <- new[j] + row[j]
    row <- new
  }
  row[1L]
}

#' Mean and quantiles of a finite mixture of normals
#'
#' The exact posterior of a linear aggregate is
#' \eqn{\sum_{\mathcal{P}} \Pr(\mathcal{P} \mid y) \, N(m_\mathcal{P},
#' s_\mathcal{P}^2)}. The mean is the weighted mean of the components;
#' quantiles come from inverting the mixture CDF numerically.
#'
#' @noRd
mixture_summary <- function(prob, mu, sd, probs = c(0.025, 0.975)) {
  keep <- prob > 1e-12
  prob <- prob[keep] / sum(prob[keep])
  mu <- mu[keep]
  sd <- pmax(sd[keep], 1e-12)

  mean_val <- sum(prob * mu)
  var_val <- sum(prob * (sd^2 + mu^2)) - mean_val^2

  cdf <- function(x) sum(prob * stats::pnorm(x, mu, sd))
  lo <- min(mu - 8 * sd)
  hi <- max(mu + 8 * sd)
  qs <- vapply(probs, function(q) {
    stats::uniroot(function(x) cdf(x) - q, lower = lo, upper = hi,
                   tol = 1e-10)$root
  }, numeric(1))

  list(mean = mean_val, sd = sqrt(max(0, var_val)),
       lower = qs[1L], upper = qs[2L])
}

#' @export
print.ph_enumeration <- function(x, ...) {
  cat("<ph_enumeration>  exact partition posterior\n\n")
  cat(sprintf("  Cells           : %d\n", x$data$K))
  cat(sprintf("  Partitions       : %s (enumerated in full)\n",
              format(x$n_partitions, big.mark = ",")))
  cat(sprintf("  Concentration   : alpha = %g\n", x$alpha))
  cat(sprintf("  Base measure    : N(%.4g, %.4g)\n", x$mu0, x$sigma0_sq))
  cat(sprintf("\n  Posterior E[# groups] : %.3f\n", x$m_mean))
  cat(sprintf("  Overall effect        : %.4f  [%.4f, %.4f]\n",
              x$overall$mean, x$overall$lower, x$overall$upper))

  ord <- order(x$prob, decreasing = TRUE)
  cat("\n  Most probable partitions:\n")
  for (i in utils::head(ord, 5L)) {
    cat(sprintf("    %.3f  %s\n", x$prob[i],
                format_partition(x$partitions[i, ], x$data$cells$label)))
  }
  invisible(x)
}
