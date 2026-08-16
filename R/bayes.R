#' The Dirichlet Process partial-homogeneity estimator
#'
#' Places a Dirichlet Process mixture prior on the cohort-time effects and
#' samples the posterior with the collapsed Gibbs sampler of Neal (2000,
#' Algorithm 3). Because draws from a Dirichlet Process are almost surely
#' discrete, realisations put positive probability on ties among the
#' \eqn{\tau_k}, and those ties *are* the partition. The prior assigns positive
#' probability to every partition of the \eqn{K} cells while favouring
#' parsimonious groupings, and never fixes their number.
#'
#' This is the estimator the paper recommends for inference. Unlike
#' [l0_ph()], which commits to one partition and reports intervals conditional
#' on it, the posterior here averages over partitions, so the reported
#' uncertainty includes uncertainty about the grouping itself. In the paper's
#' simulations that distinction is decisive: where the partition is uncertain,
#' the \eqn{\ell_0} plug-in intervals collapse to 0.57--0.62 coverage while
#' these credible intervals degrade gracefully to 0.79--0.81, and in the
#' well-separated regime they attain near-nominal 0.93--0.94.
#'
#' @section Priors:
#' The concentration `alpha` controls the prior expected number of groups
#' through eq. (23), \eqn{E[m] \approx \alpha \log(1 + K/\alpha)}; small
#' `alpha` favours pooling, large `alpha` favours the flexible fit. Use
#' [alpha_for_groups()] to solve that relation for a target. The paper's
#' simulations use `alpha = 7` at \eqn{K = 18} (so \eqn{E[m] \approx K/2}) and
#' its applications use `alpha = 1`.
#'
#' The base measure is \eqn{G_0 = N(\mu_0, \sigma_0^2)}. The defaults are
#' data-scaled, `mu0 = median(tau)` and `sigma0_sq = (10 * sd(tau))^2`,
#' matching the paper's applications and keeping the prior diffuse on whatever
#' scale the outcome happens to have. Note that a data-scaled base measure is
#' mildly empirical-Bayes; supply fixed values if you want a prior that is
#' genuinely independent of the data.
#'
#' @section How much the answer depends on alpha:
#' Not much for aggregates, and potentially a lot for the grouping when
#' \eqn{K} is small. The paper's solution paths are nearly flat in `alpha` for
#' the overall ATT in the simulation, but its first application, with only
#' seven cells, sees the overall effect slide from about -0.017 to -0.038
#' across `alpha`. Rather than fix a value, report the path: see
#' [alpha_sensitivity()].
#'
#' @section On the error variance:
#' When the object came from a micro panel, \eqn{\sigma^2} is a real parameter
#' and `sigma2 = "random"` draws it from its inverse-gamma full conditional
#' each sweep, which is the paper's primary specification. When the object
#' came from a first-stage estimator, the covariance \eqn{\hat\Sigma} is a
#' fixed plug-in and there is no scalar variance left to sample, so the
#' sampler conditions on it (Appendix C). Appendix D shows the two treatments
#' give the same coverage and interval length to two decimals, so nothing
#' hinges on this.
#'
#' @param object a [ph_data] object.
#' @param alpha the Dirichlet Process concentration parameter.
#' @param mu0,sigma0_sq mean and variance of the Gaussian base measure. `NULL`
#'   uses the data-scaled defaults described above.
#' @param iters,burn,thin sweeps per chain, burn-in to discard, and thinning
#'   interval.
#' @param chains number of chains. More than one enables the convergence
#'   diagnostics of [ph_rhat()]; chains start from dispersed partitions.
#' @param sigma2 `"auto"` (default) draws the error variance each sweep when
#'   the object carries a micro panel and conditions on the supplied covariance
#'   otherwise. `"random"` and `"fixed"` force the choice; asking for
#'   `"random"` without a micro panel is not possible and falls back with a
#'   message.
#' @param marginal `"exact"` (default) evaluates the full collapsed marginal
#'   likelihood, eq. (54), at every cluster-assignment move, so the moves see
#'   the cross-cell covariance. `"diagonal"` uses the fast conjugate
#'   normal-normal shortcut that ignores it. Appendix E reports that the
#'   shortcut reproduces the same aggregates but can misstate individual
#'   co-clustering probabilities by up to 0.39, so prefer the default unless
#'   you are deliberately reproducing the shortcut.
#' @param init starting partition: `"dispersed"` (varies by chain),
#'   `"singletons"`, `"pooled"`, `"random"`, or an explicit label vector.
#' @param a0,b0 shape and rate of the inverse-gamma prior on \eqn{\sigma^2},
#'   eq. (28). Used only when `sigma2 = "random"`.
#' @param seed optional integer for reproducibility.
#' @param progress print a progress bar.
#'
#' @return An object of class `bayes_ph` with components `tau` (the posterior
#'   mean cohort-time effects, eq. 38), `lower` and `upper` (2.5th and 97.5th
#'   posterior percentiles), `coclust` (the K x K posterior co-clustering
#'   matrix of eq. 39), `m_mean` (the posterior expected number of groups), and
#'   `draws` (the retained draws, for [aggregate.bayes_ph()] and the
#'   diagnostics).
#'
#' @references
#' Arora, P. and Wagle, R. (2026). Partial Homogeneity in Staggered
#' Difference-in-Differences, Section 3.
#'
#' Neal, R. M. (2000). Markov Chain Sampling Methods for Dirichlet Process
#' Mixture Models. \emph{JCGS} 9(2), 249--265.
#'
#' @seealso [aggregate.bayes_ph()] for partition-aware aggregate estimands,
#'   [coclustering()] for the grouping structure, [alpha_sensitivity()] for the
#'   prior path.
#'
#' @export
#' @examples
#' d <- ph_data(c(0.10, 0.11, 0.42, 0.40), Sigma = 0.02^2 * (0.3 + 0.7 * diag(4)))
#' fit <- bayes_ph(d, alpha = 1, iters = 400, burn = 100, seed = 1)
#' fit
#'
#' # The co-clustering matrix is the honest summary of the grouping: cells 1
#' # and 2 recur together, and so do 3 and 4.
#' round(coclustering(fit), 2)
bayes_ph <- function(object, alpha = 1, mu0 = NULL, sigma0_sq = NULL,
                     iters = 4000L, burn = 1000L, thin = 1L, chains = 1L,
                     sigma2 = c("auto", "random", "fixed"),
                     marginal = c("exact", "diagonal"),
                     init = "dispersed", a0 = 0.001, b0 = 0.001,
                     seed = NULL, progress = interactive()) {
  stopifnot(inherits(object, "ph_data"))
  marginal <- match.arg(marginal)
  sigma2 <- match.arg(sigma2)
  if (!is.null(seed)) set.seed(seed)

  K <- object$K
  if (burn >= iters) {
    stop(sprintf("`burn` (%d) must be smaller than `iters` (%d).", burn, iters),
         call. = FALSE)
  }
  if (alpha <= 0) stop("`alpha` must be positive.", call. = FALSE)

  if (is.null(mu0)) mu0 <- stats::median(object$tau)
  if (is.null(sigma0_sq)) {
    spread <- stats::sd(object$tau)
    if (!is.finite(spread) || spread <= 0) spread <- max(sqrt(diag(object$Sigma)))
    sigma0_sq <- (10 * spread)^2
  }
  if (sigma0_sq <= 0) stop("`sigma0_sq` must be positive.", call. = FALSE)

  # Random sigma^2 needs the outcome-scale quantities, which only the micro
  # panel route carries. On the two-stage route the first-stage covariance is a
  # fixed plug-in and there is no scalar variance left to sample (Appendix C).
  has_panel <- !is.null(object$G)
  random_sigma <- switch(sigma2,
    auto = has_panel,
    random = has_panel,
    fixed = FALSE
  )
  if (identical(sigma2, "random") && !has_panel) {
    message("`sigma2 = \"random\"` needs the micro panel, which this object ",
            "does not carry;\n  conditioning on the supplied first-stage ",
            "covariance instead (Appendix C).\n  Appendix D shows the two ",
            "give the same coverage and interval length.")
  }
  r_dim <- if (random_sigma) object$df_resid + K else NA_real_

  exact <- identical(marginal, "exact")
  prec <- diag(object$Omega)

  n_keep_per_chain <- length(seq(burn + 1L, iters, by = thin))
  n_keep <- n_keep_per_chain * chains

  tau_draws <- matrix(NA_real_, n_keep, K)
  z_draws <- matrix(NA_integer_, n_keep, K)
  m_draws <- numeric(n_keep)
  s2_draws <- rep(NA_real_, n_keep)
  chain_id <- integer(n_keep)

  pb <- NULL
  if (isTRUE(progress)) {
    pb <- utils::txtProgressBar(min = 0, max = chains * iters, style = 3)
  }
  on.exit(if (!is.null(pb)) close(pb), add = TRUE)

  kept <- 0L
  for (ch in seq_len(chains)) {
    z <- initial_partition(init, K, ch)
    s2 <- if (random_sigma) object$sigma2 else NA_real_
    Omega <- if (random_sigma) object$G / s2 else object$Omega
    Otau <- as.numeric(Omega %*% object$tau)
    prec_now <- diag(Omega)

    for (it in seq_len(iters)) {
      z <- sweep_assignments(z, Omega, Otau, object$tau, prec_now, alpha,
                             mu0, sigma0_sq, exact)
      phi <- draw_phi(z, Omega, Otau, mu0, sigma0_sq)
      tau_it <- phi[z]

      if (random_sigma) {
        # RSS(P, phi) = RSS_flex + (tau - R phi)' G (tau - R phi), eq. (56).
        resid <- object$tau - tau_it
        rss <- object$rss_flex + as.numeric(crossprod(resid, object$G %*% resid))
        s2 <- 1 / stats::rgamma(1, shape = a0 + r_dim / 2, rate = b0 + rss / 2)
        Omega <- object$G / s2
        Otau <- as.numeric(Omega %*% object$tau)
        prec_now <- diag(Omega)
      }

      if (it > burn && ((it - burn - 1L) %% thin == 0L)) {
        kept <- kept + 1L
        tau_draws[kept, ] <- tau_it
        z_draws[kept, ] <- z
        m_draws[kept] <- max(z)
        s2_draws[kept] <- s2
        chain_id[kept] <- ch
      }
      if (!is.null(pb)) utils::setTxtProgressBar(pb, (ch - 1L) * iters + it)
    }
  }

  coclust <- coclustering_from_draws(z_draws)
  qs <- apply(tau_draws, 2L, stats::quantile, probs = c(0.025, 0.975),
              names = FALSE)

  structure(
    list(
      tau = colMeans(tau_draws),
      lower = qs[1L, ],
      upper = qs[2L, ],
      sd = apply(tau_draws, 2L, stats::sd),
      coclust = coclust,
      m_mean = mean(m_draws),
      m_table = prop.table(table(m_draws)),
      draws = list(tau = tau_draws, z = z_draws, m = m_draws,
                   sigma2 = s2_draws, chain = chain_id),
      alpha = alpha, mu0 = mu0, sigma0_sq = sigma0_sq,
      iters = iters, burn = burn, thin = thin, chains = chains,
      marginal = marginal, sigma2_treatment = if (random_sigma) "random" else "fixed",
      data = object
    ),
    class = "bayes_ph"
  )
}

#' Posterior co-clustering matrix from partition draws
#'
#' Eq. (39): the posterior probability that cells \eqn{j} and \eqn{k} belong to
#' the same group, averaged over the retained draws.
#'
#' @noRd
coclustering_from_draws <- function(z_draws) {
  K <- ncol(z_draws)
  S <- nrow(z_draws)
  out <- matrix(0, K, K)
  for (s in seq_len(S)) {
    z <- z_draws[s, ]
    out <- out + outer(z, z, `==`)
  }
  out / S
}

#' Concentration parameter giving a target prior number of groups
#'
#' Inverts eq. (23), \eqn{E[m \mid \alpha, K] \approx \alpha \log(1 + K/\alpha)},
#' numerically. Useful for setting a prior that is neutral between pooling and
#' the flexible fit: `alpha_for_groups(K, K / 2)` reproduces the choice made in
#' the paper's simulations.
#'
#' @param K the number of cohort-time cells.
#' @param target the desired prior expected number of groups, between 1 and K.
#' @return the concentration parameter.
#' @export
#' @examples
#' alpha_for_groups(18, 9)   # the paper's simulation setting, alpha ~ 7
#' expected_groups(7, 1)     # E[m] at the applications' alpha = 1
alpha_for_groups <- function(K, target) {
  if (target <= 1 || target >= K) {
    stop(sprintf("`target` must lie strictly between 1 and K = %d.", K),
         call. = FALSE)
  }
  f <- function(a) a * log(1 + K / a) - target
  stats::uniroot(f, interval = c(1e-8, 1e8), tol = .Machine$double.eps^0.75)$root
}

#' @rdname alpha_for_groups
#' @param alpha the concentration parameter.
#' @export
expected_groups <- function(K, alpha) {
  alpha * log(1 + K / alpha)
}

#' @export
print.bayes_ph <- function(x, ...) {
  cat("<bayes_ph>  Dirichlet Process partial-homogeneity posterior\n\n")
  cat(sprintf("  Cells           : %d\n", x$data$K))
  cat(sprintf("  Concentration   : alpha = %g  (prior E[m] = %.2f)\n",
              x$alpha, expected_groups(x$data$K, x$alpha)))
  cat(sprintf("  Base measure    : N(%.4g, %.4g)\n", x$mu0, x$sigma0_sq))
  cat(sprintf("  Sampler         : %d chain(s), %d sweeps, %d burn-in%s\n",
              x$chains, x$iters, x$burn,
              if (x$thin > 1L) sprintf(", thin %d", x$thin) else ""))
  cat(sprintf("  Assignment moves: %s covariance\n", x$marginal))
  cat(sprintf("  Error variance  : %s\n", x$sigma2_treatment))
  cat(sprintf("\n  Posterior E[# groups] : %.2f\n", x$m_mean))

  tab <- data.frame(
    cell = x$data$cells$label,
    flexible = round(x$data$tau, 4),
    posterior = round(x$tau, 4),
    lower = round(x$lower, 4),
    upper = round(x$upper, 4)
  )
  cat("\n")
  print(utils::head(tab, 20L), row.names = FALSE)
  if (x$data$K > 20L) cat(sprintf("  ... %d more cells\n", x$data$K - 20L))
  cat("\n  Intervals are 2.5-97.5% posterior percentiles and already ",
      "marginalize\n  over the partition. Use coclustering() to see the ",
      "grouping structure and\n  aggregate() for the overall ATT or an ",
      "event study.\n", sep = "")
  invisible(x)
}

#' @export
coef.bayes_ph <- function(object, ...) {
  stats::setNames(object$tau, object$data$cells$label)
}
