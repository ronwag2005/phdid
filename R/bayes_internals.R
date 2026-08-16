## Internals for the collapsed Gibbs sampler. Not exported.
##
## Everything here works in the two-stage coordinates of Lemma 1: the data
## enter only through the group-summed precision S = R'Omega R and the
## group-summed score u = R'Omega tau. Both are computed with `rowsum()`, which
## does the summation in C, so a candidate assignment costs O(K^2 m) rather
## than the O(r) implied by working with the full r x m design.

#' Group-summed precision, S = R' Omega R
#'
#' @param Omega the K x K precision matrix.
#' @param z partition labels, canonicalised.
#' @return an m x m symmetric matrix.
#' @noRd
group_precision <- function(Omega, z) {
  A <- rowsum(Omega, z, reorder = TRUE)          # m x K
  S <- rowsum(t(A), z, reorder = TRUE)           # m x m
  (S + t(S)) / 2
}

#' Group-summed score, u = R' Omega tau
#' @noRd
group_score <- function(Otau, z) {
  as.numeric(rowsum(Otau, z, reorder = TRUE))
}

#' Collapsed log marginal likelihood of a partition
#'
#' Evaluates eq. (30), in the general-covariance form of eq. (54),
#' \deqn{\log p(\tilde y \mid \mathcal{P}) \propto
#'   -\tfrac{1}{2}\log|I_m + \sigma_0^2 S|
#'   + \tfrac{\sigma_0^2}{2} u_0' (I_m + \sigma_0^2 S)^{-1} u_0,}
#' with \eqn{u_0 = u - \mu_0 S 1_m}. Terms constant across partitions are
#' dropped; this is legitimate because \eqn{R 1_m = 1_K} for every partition,
#' so the centred outcome does not depend on \eqn{\mathcal{P}}.
#'
#' The Woodbury identity and the matrix determinant lemma are what make this
#' an m x m calculation instead of an r x r one.
#'
#' @noRd
log_marginal <- function(z, Omega, Otau, mu0, s0sq) {
  S <- group_precision(Omega, z)
  u <- group_score(Otau, z)
  m <- length(u)

  u0 <- u - mu0 * rowSums(S)
  A <- diag(m) + s0sq * S

  ch <- tryCatch(chol(A), error = function(e) NULL)
  if (is.null(ch)) return(-Inf)

  quad <- sum(u0 * backsolve(ch, backsolve(ch, u0, transpose = TRUE)))
  -sum(log(diag(ch))) + 0.5 * s0sq * quad
}

#' Diagonal (independence) approximation to the assignment move
#'
#' The fast conjugate normal-normal update used by the paper's replication
#' scripts for the cluster-assignment step: cell \eqn{k} is scored against each
#' cluster using only its own precision and the cluster's pooled precision,
#' ignoring the cross-cell covariance.
#'
#' Appendix E reports that this reproduces the same aggregates but can misstate
#' individual co-clustering probabilities by up to 0.39 in a correlated design,
#' which is why [bayes_ph()] uses the exact form by default.
#'
#' @noRd
log_marginal_diag_move <- function(k, cluster_members, tau, prec, mu0, s0sq) {
  b0 <- 1 / s0sq
  Pc <- b0 + sum(prec[cluster_members])
  mc <- (mu0 * b0 + sum(prec[cluster_members] * tau[cluster_members])) / Pc
  stats::dnorm(tau[k], mc, sqrt(1 / prec[k] + 1 / Pc), log = TRUE)
}

#' @noRd
log_marginal_diag_new <- function(k, tau, prec, mu0, s0sq) {
  stats::dnorm(tau[k], mu0, sqrt(1 / prec[k] + s0sq), log = TRUE)
}

#' Draw the group effects from their conjugate posterior
#'
#' Eq. (31) / (55): \eqn{\phi \mid \mathcal{P}, \Omega \sim N(\mu_*,
#' \Sigma_*)} with \eqn{\Sigma_* = (S + I_m/\sigma_0^2)^{-1}} and
#' \eqn{\mu_* = \Sigma_*(u + \mu_0 1_m/\sigma_0^2)}. This is the
#' ridge-regularised GLS estimator of the group effects; in the diffuse limit
#' it is the restricted GLS fit of [ph_fit()].
#'
#' @noRd
draw_phi <- function(z, Omega, Otau, mu0, s0sq) {
  S <- group_precision(Omega, z)
  u <- group_score(Otau, z)
  m <- length(u)

  P <- S + diag(m) / s0sq
  ch <- chol(P)
  mu <- backsolve(ch, backsolve(ch, u + mu0 / s0sq, transpose = TRUE))
  # x = mu + L^-T eps has covariance P^-1 when P = R'R with R = chol(P).
  as.numeric(mu + backsolve(ch, stats::rnorm(m)))
}

#' Sample one sweep of cluster assignments (Neal 2000, Algorithm 3)
#'
#' Cycles over the K cells. For each, the cell is removed from its cluster
#' (deleting the cluster if it empties), the conditional probabilities of
#' eq. (37) are formed against every surviving cluster and against a fresh
#' cluster weighted by alpha, and a new assignment is drawn.
#'
#' @noRd
sweep_assignments <- function(z, Omega, Otau, tau, prec, alpha, mu0, s0sq,
                              exact) {
  K <- length(z)
  for (k in seq_len(K)) {
    z_minus <- z
    z_minus[k] <- NA_integer_
    others <- which(!is.na(z_minus))
    labs <- sort(unique(z_minus[others]))
    n_cl <- length(labs)

    logp <- numeric(n_cl + 1L)
    for (j in seq_len(n_cl)) {
      members <- others[z_minus[others] == labs[j]]
      logp[j] <- log(length(members)) +
        if (exact) {
          trial <- z_minus
          trial[k] <- labs[j]
          log_marginal(relabel(trial), Omega, Otau, mu0, s0sq)
        } else {
          log_marginal_diag_move(k, members, tau, prec, mu0, s0sq)
        }
    }
    logp[n_cl + 1L] <- log(alpha) +
      if (exact) {
        trial <- z_minus
        trial[k] <- max(labs, 0L) + 1L
        log_marginal(relabel(trial), Omega, Otau, mu0, s0sq)
      } else {
        log_marginal_diag_new(k, tau, prec, mu0, s0sq)
      }

    if (all(!is.finite(logp))) {
      # Degenerate sweep; leave the cell where it was rather than crash.
      next
    }
    p <- exp(logp - max(logp[is.finite(logp)]))
    p[!is.finite(p)] <- 0
    if (sum(p) <= 0) next
    pick <- sample.int(length(p), 1L, prob = p)

    z <- z_minus
    z[k] <- if (pick <= n_cl) labs[pick] else max(labs, 0L) + 1L
    z <- relabel(z)
  }
  z
}

#' Dispersed initial partitions for multi-chain diagnostics
#'
#' Appendix E starts four chains from all singletons, one pooled group, and
#' two random partitions, so that convergence is judged from genuinely
#' different corners of the partition space.
#'
#' @noRd
initial_partition <- function(init, K, chain) {
  if (is.numeric(init)) return(relabel(as_partition(init, K)))
  if (identical(init, "singletons")) return(seq_len(K))
  if (identical(init, "pooled")) return(rep(1L, K))
  if (identical(init, "random")) {
    return(relabel(sample.int(max(2L, K %/% 2L), K, replace = TRUE)))
  }
  # "dispersed": cycle deterministically across chains.
  switch(((chain - 1L) %% 4L) + 1L,
    seq_len(K),
    rep(1L, K),
    relabel(sample.int(max(2L, min(9L, K)), K, replace = TRUE)),
    relabel(sample.int(min(3L, K), K, replace = TRUE))
  )
}
