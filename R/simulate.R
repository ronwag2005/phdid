#' A calibrated partial-homogeneity design
#'
#' Sets up the balanced staggered panel of the paper's Section 4.1: \eqn{N}
#' units split equally across a never-treated group and several treated
#' cohorts, observed over \eqn{T} periods. The default is the paper's headline
#' design, \eqn{N = 2000}, \eqn{T = 10}, cohorts entering at \eqn{t = 3, 5, 7},
#' which yields \eqn{K = 18} post-treatment cohort-time cells.
#'
#' The design is built once and reused across replications, since the unit and
#' time fixed effects are removed by the within transformation and only the
#' errors vary. Everything the estimators need reduces to the \eqn{K \times K}
#' within cross-product \eqn{\tilde D'\tilde D}.
#'
#' @param N number of units; must divide evenly among the cohorts.
#' @param T number of periods.
#' @param cohorts treatment entry periods, with `0` marking the never-treated
#'   group.
#' @param sigma the error standard deviation.
#'
#' @return An object of class `ph_design` carrying the cell list, the within
#'   cross-product, the effective sample sizes \eqn{n_k = \|\tilde D_k\|^2},
#'   and the mean flexible standard error used to calibrate separation.
#'
#' @seealso [ph_truth()] to build a partial-homogeneity parameter vector on
#'   this design, [ph_sample()] to draw one replication, and [sim_study()] to
#'   run a Monte Carlo.
#'
#' @export
#' @examples
#' des <- ph_design()
#' des
ph_design <- function(N = 2000L, T = 10L, cohorts = c(0, 3, 5, 7), sigma = 1) {
  if (N %% length(cohorts) != 0L) {
    stop(sprintf("`N` (%d) must divide evenly among the %d cohorts.",
                 N, length(cohorts)), call. = FALSE)
  }
  treated <- cohorts[cohorts > 0]
  if (!length(treated)) stop("At least one treated cohort is needed.",
                             call. = FALSE)

  cohort_full <- rep(cohorts, each = N / length(cohorts))
  cells <- do.call(rbind, lapply(treated, function(g) {
    data.frame(g = g, t = seq(g, T))
  }))
  K <- nrow(cells)

  Dt <- matrix(0, N * T, K)
  for (k in seq_len(K)) {
    Dk <- matrix(0, N, T)
    Dk[cohort_full == cells$g[k], cells$t[k]] <- 1
    Dt[, k] <- as.vector(double_demean(Dk))
  }
  G <- crossprod(Dt)
  nk <- diag(G)

  structure(
    list(N = N, T = T, cohorts = cohorts, sigma = sigma,
         cohort = cohort_full, cells = cells, K = K,
         Dt = Dt, G = G, Ginv = safe_inverse(G, "within cross-product"),
         nk = nk,
         nobs = N * T,
         df_resid = N * T - N - T + 1L - K,
         sd_flex = mean(sigma / sqrt(nk))),
    class = "ph_design"
  )
}

#' Double-demeaning, the within transformation for a balanced panel
#' @noRd
double_demean <- function(M) {
  M - rowMeans(M) - rep(colMeans(M), each = nrow(M)) + mean(M)
}

#' @export
print.ph_design <- function(x, ...) {
  cat("<ph_design>  simulated staggered panel\n\n")
  cat(sprintf("  N = %d units, T = %d periods, %d observations\n",
              x$N, x$T, x$nobs))
  cat(sprintf("  Cohorts: %s  (0 = never treated)\n",
              paste(x$cohorts, collapse = ", ")))
  cat(sprintf("  K = %d post-treatment cohort-time cells\n", x$K))
  cat(sprintf("  sigma = %g;  mean flexible SE = %.4f\n", x$sigma, x$sd_flex))
  cat("\n  Separation delta is measured in units of that flexible SE,\n")
  cat("  so ph_truth(des, m_star = 6, delta = 6) puts adjacent group means\n")
  cat("  six flexible standard errors apart.\n")
  invisible(x)
}


#' A partial-homogeneity parameter vector
#'
#' Partitions the \eqn{K} cells into `m_star` near-equal groups and assigns
#' each group a mean on an evenly spaced grid, with adjacent gap
#' \eqn{\Delta = \delta \times \mathrm{sd}(\hat\tau_{\mathrm{flex}})}. The
#' unit-free separation \eqn{\delta} is the distance between adjacent group
#' means measured in standard errors of the flexible estimates, and it is what
#' governs whether the groups can be told apart at all.
#'
#' The paper sweeps \eqn{m^* \in \{1, 3, 6, 9, 18\}} and
#' \eqn{\delta \in \{3, 6, 12\}}, with \eqn{\delta = 6} as the headline. At
#' \eqn{\delta = 3} the partition is recovered only about half the time and the
#' feasible estimators are *less* precise than flexible TWFE; at
#' \eqn{\delta = 12} they essentially attain the oracle.
#'
#' @param design a [ph_design] object.
#' @param m_star the true number of distinct effects; `1` is fully homogeneous
#'   and `design$K` fully heterogeneous.
#' @param delta the separation, in flexible standard errors.
#' @return A list with `labels` (the true partition) and `tau` (the true
#'   cohort-time effects).
#' @export
#' @examples
#' des <- ph_design()
#' truth <- ph_truth(des, m_star = 6, delta = 6)
#' table(truth$labels)
ph_truth <- function(design, m_star, delta) {
  stopifnot(inherits(design, "ph_design"))
  K <- design$K
  if (m_star < 1L || m_star > K) {
    stop(sprintf("`m_star` must be between 1 and K = %d.", K), call. = FALSE)
  }
  Delta <- delta * design$sd_flex

  sizes <- rep(K %/% m_star, m_star)
  rem <- K %% m_star
  if (rem > 0L) sizes[seq_len(rem)] <- sizes[seq_len(rem)] + 1L
  labels <- rep(seq_len(m_star), times = sizes)
  mus <- (seq_len(m_star) - (m_star + 1) / 2) * Delta

  list(labels = labels, tau = mus[labels], m_star = m_star, delta = delta)
}


#' Draw one replication from a partial-homogeneity design
#'
#' Generates within-transformed outcomes under the given true effects, fits the
#' fully flexible model, and returns the result as a [ph_data] object ready for
#' any estimator in the package. Unit and time fixed effects are removed by the
#' transformation, so the sampling behaviour is driven entirely by the errors,
#' as in eq. (40).
#'
#' @param design a [ph_design] object.
#' @param tau_true the true cohort-time effects, or a [ph_truth()] list.
#' @return A [ph_data] object.
#' @export
#' @examples
#' des <- ph_design(N = 300, T = 6, cohorts = c(0, 3, 5))
#' truth <- ph_truth(des, m_star = 2, delta = 6)
#' d <- ph_sample(des, truth)
#' l0_ph(d)$m
ph_sample <- function(design, tau_true) {
  stopifnot(inherits(design, "ph_design"))
  if (is.list(tau_true)) tau_true <- tau_true$tau
  if (length(tau_true) != design$K) {
    stop(sprintf("`tau_true` must have length K = %d.", design$K),
         call. = FALSE)
  }

  TAU <- matrix(0, design$N, design$T)
  for (k in seq_len(design$K)) {
    TAU[design$cohort == design$cells$g[k], design$cells$t[k]] <- tau_true[k]
  }
  noise <- matrix(stats::rnorm(design$N * design$T, 0, design$sigma),
                  design$N, design$T)
  ytil <- as.vector(double_demean(TAU + noise))

  cvec <- crossprod(design$Dt, ytil)
  tau <- as.numeric(design$Ginv %*% cvec)
  rss_flex <- sum(ytil^2) - sum(cvec * tau)
  sigma2 <- rss_flex / design$df_resid

  new_ph_data(
    tau = tau,
    Sigma = sigma2 * design$Ginv,
    cells = data.frame(label = paste0(design$cells$g, ":", design$cells$t),
                       g = design$cells$g, t = design$cells$t,
                       stringsAsFactors = FALSE),
    weights = "equal",
    n_units = design$N,
    source = "panel",
    G = design$G, sigma2 = sigma2, rss_flex = rss_flex,
    nobs = design$nobs, df_resid = design$df_resid,
    warn_diagonal = FALSE
  )
}


#' Monte Carlo study of the partial-homogeneity estimators
#'
#' Reproduces the experiments of Section 4: for a given true number of distinct
#' effects and separation, it draws replications and reports, for each
#' estimator, the sampling variance of the cohort-time estimates relative to
#' flexible TWFE, the average absolute bias, how well the true partition is
#' recovered (adjusted Rand index), and the coverage and length of nominal
#' intervals for both the cohort-time effects and the overall ATT.
#'
#' The oracle estimator is handed the true partition and is infeasible; it is
#' the efficiency floor the feasible estimators are chasing.
#'
#' @section Cost:
#' The Bayesian estimator dominates the run time, since every replication runs
#' a full chain. The paper uses 500 replications; the defaults here are much
#' smaller so that an example finishes quickly. Reproducing a published row
#' takes hours, not seconds.
#'
#' @param design a [ph_design] object.
#' @param m_star,delta the truth, passed to [ph_truth()].
#' @param R number of replications.
#' @param methods which estimators to run.
#' @param level nominal interval level.
#' @param bayes_args a list of arguments forwarded to [bayes_ph()].
#' @param l0_args a list of arguments forwarded to [l0_ph()].
#' @param seed optional integer for reproducibility.
#' @param progress print a progress bar.
#'
#' @return A data frame with one row per estimator.
#'
#' @references
#' Arora, P. and Wagle, R. (2026). Partial Homogeneity in Staggered
#' Difference-in-Differences, Tables 2 to 5.
#'
#' @export
#' @examples
#' \donttest{
#' des <- ph_design(N = 300, T = 6, cohorts = c(0, 3, 5))
#' sim_study(des, m_star = 2, delta = 6, R = 10,
#'           bayes_args = list(iters = 200, burn = 50), seed = 1)
#' }
sim_study <- function(design, m_star, delta, R = 100L,
                      methods = c("pooled", "flexible", "oracle", "l0",
                                  "bayes"),
                      level = 0.95,
                      bayes_args = list(iters = 600L, burn = 200L),
                      l0_args = list(), seed = NULL,
                      progress = interactive()) {
  stopifnot(inherits(design, "ph_design"))
  methods <- match.arg(methods, several.ok = TRUE)
  if (!is.null(seed)) set.seed(seed)

  truth <- ph_truth(design, m_star, delta)
  tau_true <- truth$tau
  K <- design$K
  w <- rep(1 / K, K)
  att_true <- sum(w * tau_true)
  z_crit <- stats::qnorm(1 - (1 - level) / 2)

  est <- lapply(methods, function(.) matrix(NA_real_, R, K))
  names(est) <- methods
  cover_cell <- len_cell <- cover_att <- len_att <-
    matrix(NA_real_, R, length(methods), dimnames = list(NULL, methods))
  ari <- mhat <- matrix(NA_real_, R, length(methods),
                        dimnames = list(NULL, methods))

  pb <- if (isTRUE(progress)) utils::txtProgressBar(0, R, style = 3) else NULL
  on.exit(if (!is.null(pb)) close(pb), add = TRUE)

  for (r in seq_len(R)) {
    d <- ph_sample(design, tau_true)

    for (mth in methods) {
      out <- switch(mth,
        pooled = summarise_fit(pooled_twfe(d), w, z_crit),
        flexible = summarise_fit(flex_twfe(d), w, z_crit),
        oracle = summarise_fit(ph_fit(d, truth$labels), w, z_crit),
        l0 = summarise_fit(do.call(l0_ph, c(list(d), l0_args)), w, z_crit),
        bayes = summarise_bayes(
          do.call(bayes_ph, c(list(d, progress = FALSE), bayes_args)),
          w, level)
      )
      est[[mth]][r, ] <- out$tau
      cover_cell[r, mth] <- mean(tau_true >= out$lower & tau_true <= out$upper)
      len_cell[r, mth] <- mean(out$upper - out$lower)
      cover_att[r, mth] <- as.numeric(att_true >= out$att_lower &&
                                        att_true <= out$att_upper)
      len_att[r, mth] <- out$att_upper - out$att_lower
      mhat[r, mth] <- out$m
      ari[r, mth] <- adjusted_rand(out$partition, truth$labels)
    }
    if (!is.null(pb)) utils::setTxtProgressBar(pb, r)
  }

  var_flex <- if ("flexible" %in% methods) {
    mean(apply(est[["flexible"]], 2L, stats::var))
  } else NA_real_

  out <- do.call(rbind, lapply(methods, function(mth) {
    E <- est[[mth]]
    cell_var <- apply(E, 2L, stats::var)
    cell_bias <- colMeans(E) - tau_true
    data.frame(
      method = mth,
      m_star = m_star, delta = delta, R = R,
      variance = mean(cell_var),
      var_ratio = mean(cell_var) / var_flex,
      abs_bias = mean(abs(cell_bias)),
      rmse = mean(sqrt(cell_bias^2 + cell_var)),
      m_hat = mean(mhat[, mth]),
      ari = mean(ari[, mth]),
      cover_CATT = mean(cover_cell[, mth]),
      len_CATT = mean(len_cell[, mth]),
      cover_ATT = mean(cover_att[, mth]),
      len_ATT = mean(len_att[, mth]),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

#' Common summary of a frequentist fit for the Monte Carlo
#' @noRd
summarise_fit <- function(fit, w, z_crit) {
  se <- sqrt(pmax(diag(fit$vcov), 0))
  att <- sum(w * fit$tau)
  att_se <- sqrt(max(0, as.numeric(crossprod(w, fit$vcov %*% w))))
  list(tau = fit$tau,
       lower = fit$tau - z_crit * se, upper = fit$tau + z_crit * se,
       att_lower = att - z_crit * att_se, att_upper = att + z_crit * att_se,
       m = fit$m, partition = fit$partition)
}

#' Common summary of a posterior for the Monte Carlo
#' @noRd
summarise_bayes <- function(fit, w, level) {
  a <- (1 - level) / 2
  att_draws <- as.numeric(fit$draws$tau %*% w)
  qs <- stats::quantile(att_draws, c(a, 1 - a), names = FALSE)
  cell_qs <- apply(fit$draws$tau, 2L, stats::quantile, probs = c(a, 1 - a),
                   names = FALSE)
  list(tau = fit$tau,
       lower = cell_qs[1L, ], upper = cell_qs[2L, ],
       att_lower = qs[1L], att_upper = qs[2L],
       m = fit$m_mean,
       partition = point_partition(fit, loss = "binder"))
}
