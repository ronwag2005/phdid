#' Fit the grouped model under a given partition
#'
#' The workhorse of the package. Given a partition \eqn{\mathcal{P}} of the
#' \eqn{K} cohort-time cells with group-indicator matrix \eqn{R}, this solves
#' the restricted generalised least squares problem of eq. (50) and (55),
#' \deqn{\hat\phi = (R'\hat\Sigma^{-1}R)^{-1} R'\hat\Sigma^{-1}\hat\tau,
#'   \qquad \mathrm{Var}(\hat\phi) = (R'\hat\Sigma^{-1}R)^{-1},}
#' and assigns each cell its group effect, \eqn{\hat\tau_k = \hat\phi_p} for
#' \eqn{k \in C_p}.
#'
#' Two special cases are worth naming. The partition of \eqn{K} singletons
#' returns the fully flexible estimator and its covariance unchanged; the
#' partition with a single group returns the fully pooled estimator of eq. (6).
#' [flex_twfe()] and [pooled_twfe()] are thin wrappers on those two calls.
#'
#' @section Why this is the estimator and not the group means:
#' The obvious shortcut --- averaging the flexible estimates within each group
#' --- coincides with this fit only when the within-transformed cohort-time
#' dummies are orthogonal (eq. 49). In a general panel, cells that share a
#' cohort share unit fixed effects and cells in the same period share time
#' fixed effects, so \eqn{R'\hat\Sigma^{-1}R} has non-zero off-diagonal
#' entries, the group means fail the normal equations, and they are no longer
#' best linear unbiased. Appendix B of the paper makes this argument in full.
#' It is also why [l0_ph()] re-estimates by GLS after its greedy search rather
#' than reporting the merged means.
#'
#' @param object a [ph_data] object.
#' @param partition an integer vector of length K assigning each cell to a
#'   group. Labels need not be consecutive; they are canonicalised internally.
#'   A list of index vectors is also accepted.
#'
#' @return An object of class `ph_fit`: a list with `tau` (the K grouped cell
#'   effects), `vcov` (their K x K covariance \eqn{R S^{-1} R'}), `phi` and
#'   `phi_vcov` (the m distinct group effects), `partition`, `m`, and
#'   `deviance` (the GLS deviance under the partition, used by the information
#'   criteria).
#'
#' @references
#' Arora, P. and Wagle, R. (2026). Partial Homogeneity in Staggered
#' Difference-in-Differences. Eq. (12), (50), (55) and Appendix B.
#'
#' @export
#' @examples
#' tau <- c(0.10, 0.11, 0.42, 0.40)
#' d <- ph_data(tau, Sigma = 0.02^2 * (0.3 + 0.7 * diag(4)))
#'
#' # The oracle fit, handed the true partition: cells 1-2 share an effect and
#' # cells 3-4 share another.
#' fit <- ph_fit(d, c(1, 1, 2, 2))
#' fit
#'
#' # Pooling two cells halves the variance in a balanced orthogonal design.
#' sqrt(diag(fit$vcov)) / sqrt(diag(d$Sigma))
ph_fit <- function(object, partition) {
  stopifnot(inherits(object, "ph_data"))
  z <- as_partition(partition, object$K)

  R <- indicator_matrix(z)
  OR <- object$Omega %*% R
  S <- crossprod(R, OR)
  u <- as.numeric(crossprod(OR, object$tau))

  Sinv <- safe_inverse(S, "grouped precision R'Sigma^-1 R")
  phi <- as.numeric(Sinv %*% u)
  tau_hat <- as.numeric(R %*% phi)
  V <- R %*% Sinv %*% t(R)
  deviance <- object$tOt - sum(u * phi)

  structure(
    list(
      tau = tau_hat,
      vcov = V,
      phi = phi,
      phi_vcov = Sinv,
      partition = z,
      m = max(z),
      deviance = deviance,
      rss = rss_from_deviance(object, deviance),
      data = object
    ),
    class = "ph_fit"
  )
}

#' Coerce a user-supplied partition to a canonical label vector
#' @noRd
as_partition <- function(partition, K) {
  if (is.list(partition)) {
    z <- integer(K)
    for (p in seq_along(partition)) {
      idx <- partition[[p]]
      if (any(idx < 1 | idx > K)) {
        stop("Partition contains cell indices outside 1..K.", call. = FALSE)
      }
      z[idx] <- p
    }
    if (any(z == 0L)) {
      stop("The partition does not cover every cell; cells ",
           paste(which(z == 0L), collapse = ", "), " are unassigned.",
           call. = FALSE)
    }
    return(relabel(z))
  }
  if (length(partition) != K) {
    stop(sprintf("`partition` must have length %d, one label per cell, not %d.",
                 K, length(partition)), call. = FALSE)
  }
  if (anyNA(partition)) {
    stop("`partition` contains missing labels.", call. = FALSE)
  }
  relabel(partition)
}

#' Residual sum of squares implied by a GLS deviance
#'
#' On the micro-panel route the outcome-scale RSS of eq. (32) is available,
#' because \eqn{\mathrm{RSS}(\mathcal{P}) = \mathrm{RSS}_{\mathrm{flex}} +
#' \hat\sigma^2 \, \mathrm{dev}(\mathcal{P})}. On the two-stage route only the
#' deviance is defined, and this returns `NA`.
#'
#' @noRd
rss_from_deviance <- function(object, deviance) {
  if (is.na(object$rss_flex) || is.na(object$sigma2)) return(NA_real_)
  object$rss_flex + object$sigma2 * deviance
}


#' The fully flexible and fully pooled benchmarks
#'
#' The two corners of the partition problem, provided so that every reported
#' comparison runs through the same code path as the estimators being compared.
#'
#' `flex_twfe()` estimates every cohort-time effect as its own parameter, the
#' partition into \eqn{K} singletons. It is unbiased but, when some effects are
#' in fact equal, wastes \eqn{K - m} degrees of freedom and reports
#' unnecessarily wide intervals. It returns the first-stage estimates and
#' covariance unchanged, which makes it a useful identity check.
#'
#' `pooled_twfe()` imposes a single common coefficient on all \eqn{K} cells,
#' eq. (6). Because the grouped regressor is the sum of the clean cohort-time
#' dummies, \eqn{\tilde D_{\mathrm{pool}} = \sum_k \tilde D_k}, this is the
#' pooled TWFE coefficient. It is the most precise estimator available and is
#' biased for every individual effect whenever the effects genuinely differ,
#' by eq. (18); its intervals are short but centred on the wrong estimand,
#' which is why the paper's coverage table reports 0.04 for it.
#'
#' @param object a [ph_data] object.
#' @return A [ph_fit] object.
#' @export
#' @examples
#' d <- ph_data(c(0.10, 0.11, 0.42, 0.40), Sigma = 0.02^2 * (0.3 + 0.7 * diag(4)))
#' flex_twfe(d)$tau     # unchanged first-stage estimates
#' pooled_twfe(d)$tau   # one number, repeated
flex_twfe <- function(object) {
  ph_fit(object, seq_len(object$K))
}

#' @rdname flex_twfe
#' @export
pooled_twfe <- function(object) {
  ph_fit(object, rep(1L, object$K))
}


#' Information criteria for a partition
#'
#' Scores a partition by fit against complexity. Two forms are used, matching
#' the two routes into the package.
#'
#' On the micro-panel route the BIC of eq. (32) is available on the outcome
#' scale,
#' \deqn{\mathrm{BIC}(\mathcal{P}) = NT \log(\mathrm{RSS}(\mathcal{P})/NT)
#'   + m \log(NT),}
#' which is the Schwarz approximation to the log marginal density with
#' \eqn{\sigma^2} profiled out.
#'
#' On the two-stage route the first-stage covariance is taken as known, so the
#' scale is fixed and the deviance form applies,
#' \deqn{\mathrm{BIC}(\mathcal{P}) = \mathrm{dev}(\mathcal{P}) + m \log(n).}
#'
#' @section Choosing `n_bic`:
#' The deviance form needs a sample size for the penalty, and the answer is not
#' forced by the theory: the \eqn{K} first-stage estimates are the
#' observations, while the asymptotics that justify Schwarz's rate run in the
#' number of units. The default is `K`, which is what the paper's replication
#' code uses and what reproduces its reported group counts. Passing
#' `n_bic = object$n_units` gives a markedly stronger penalty and hence coarser
#' partitions. Because the choice moves the answer, it is exposed rather than
#' buried, and it is worth reporting alongside any selected partition.
#'
#' @param object a [ph_data] object.
#' @param deviance the GLS deviance under the partition.
#' @param m the number of groups.
#' @param n_bic the sample size entering the penalty; see above.
#' @return the BIC value; lower is better.
#' @noRd
partition_bic <- function(object, deviance, m, n_bic = NULL) {
  if (identical(object$source, "panel") && !is.na(object$rss_flex)) {
    NT <- object$nobs
    rss <- object$rss_flex + object$sigma2 * deviance
    return(NT * log(rss / NT) + m * log(NT))
  }
  if (is.null(n_bic)) n_bic <- object$K
  deviance + m * log(n_bic)
}


#' @export
print.ph_fit <- function(x, ...) {
  cat("<ph_fit>  grouped cohort-time effects under a fixed partition\n\n")
  cat(sprintf("  Cells      : %d\n", length(x$tau)))
  cat(sprintf("  Groups (m) : %d\n", x$m))
  cat(sprintf("  Partition  : %s\n", format_partition(x$partition,
                                                      x$data$cells$label)))
  cat("\n")

  se_flex <- sqrt(diag(x$data$Sigma))
  se <- sqrt(diag(x$vcov))
  tab <- data.frame(
    cell = x$data$cells$label,
    flexible = round(x$data$tau, 4),
    group = x$partition,
    grouped = round(x$tau, 4),
    se = round(se, 4),
    var_ratio = round((se / se_flex)^2, 2)
  )
  print(utils::head(tab, 20L), row.names = FALSE)
  if (length(x$tau) > 20L) {
    cat(sprintf("  ... %d more cells\n", length(x$tau) - 20L))
  }
  cat("\n  var_ratio is the grouped sampling variance relative to flexible;",
      "\n  below 1 means pooling bought precision for that cell (eq. 16).\n")
  cat("\n  Intervals from this object condition on the partition being ",
      "correct.\n  See Remark 3, and prefer bayes_ph() for honest ",
      "uncertainty.\n", sep = "")
  invisible(x)
}

#' @export
coef.ph_fit <- function(object, ...) {
  stats::setNames(object$tau, object$data$cells$label)
}

#' @export
vcov.ph_fit <- function(object, ...) {
  V <- object$vcov
  dimnames(V) <- list(object$data$cells$label, object$data$cells$label)
  V
}
