#' Assemble the inputs for a partial-homogeneity analysis
#'
#' Every estimator in \pkg{phdid} operates on the same two objects: a vector
#' \eqn{\hat\tau} of first-stage cohort-time effects and their joint sampling
#' covariance \eqn{\hat\Sigma}. This is not a convenience --- it is Lemma 1 of
#' the paper, which shows that the collapsed marginal likelihood and the
#' group-effect posterior depend on the data only through
#' \eqn{R'\hat\Sigma^{-1}\hat\tau} and \eqn{R'\hat\Sigma^{-1}R}. Any consistent
#' first-stage estimator may therefore be used, and nothing downstream depends
#' on the first stage beyond that pair.
#'
#' `ph_data()` is generic, with three entry points:
#'
#' \describe{
#'   \item{a numeric vector}{`ph_data(tau, Sigma)` takes the estimates and
#'     covariance directly, for any first stage you have run yourself
#'     (\pkg{fixest}, a stacked event study, an imputation estimator, ...).}
#'   \item{an object of class `MP`}{the output of [did::att_gt()]. The
#'     post-treatment cells are selected and the exact joint covariance is
#'     taken from the estimator's influence functions.}
#'   \item{a data frame}{a micro panel, from which the fully flexible
#'     (interacted) TWFE model of eq. (3) is fitted internally.}
#' }
#'
#' @section The covariance matters:
#' \eqn{\hat\Sigma} is generally **not** diagonal: cohort-time cells share
#' units and share the unit and time fixed effects, so their estimates are
#' correlated. Remark 1 of the paper is explicit that using the exact
#' cross-cell covariance is what makes the reported intervals honest. Treating
#' the cells as independent understates the posterior variance of linear
#' aggregates such as the overall ATT --- by a factor of about 3.5 in the
#' paper's design --- and misstates individual co-clustering probabilities by
#' up to 0.39. `ph_data()` warns when it is handed a diagonal covariance so
#' that the choice is at least deliberate.
#'
#' @param x a numeric vector of cohort-time effects, an `MP` object from
#'   [did::att_gt()], or a data frame holding a micro panel.
#' @param ... passed to methods.
#'
#' @return An object of class `ph_data`, a list with components `tau`
#'   (the K first-stage effects), `Sigma` (their K x K covariance), `Omega`
#'   (its inverse), `cells` (a data frame of cohort and time labels),
#'   `weights` (the aggregation weights defining the overall ATT), and
#'   book-keeping used by the information criteria.
#'
#' @references
#' Arora, P. and Wagle, R. (2026). Lemma 1 and Remark 1.
#' See \code{citation("phdid")} for the full reference.
#'
#' @seealso [l0_ph()] and [bayes_ph()], which consume this object.
#' @export
ph_data <- function(x, ...) {
  UseMethod("ph_data")
}

#' @rdname ph_data
#'
#' @param Sigma the K x K joint sampling covariance of `x`. A K-vector of
#'   standard errors is accepted as a shorthand for the diagonal covariance,
#'   but see the warning above.
#' @param cells optional data frame or character vector labelling the K cells.
#'   If a data frame, columns `g` (cohort) and `t` (time) are used where
#'   present.
#' @param weights the aggregation weights defining the overall ATT, eq. (2).
#'   Either `"equal"` (the default: each cell gets weight 1/K, as in the
#'   simulations), or a numeric vector of length K, which is normalised to sum
#'   to one.
#' @param n_units the number of units underlying the first stage. Used only to
#'   set the default sample size in the BIC penalty; see [l0_ph()].
#' @param check_diagonal warn when `Sigma` is diagonal. Set to `FALSE` when the
#'   design genuinely delivers uncorrelated cells and you do not want the
#'   reminder.
#'
#' @export
#' @examples
#' # A first stage you have run yourself: four cohort-time effects, two of
#' # which share a common value. The cells share units, so the covariance is
#' # correlated rather than diagonal.
#' tau <- c(0.10, 0.11, 0.42, 0.40)
#' Sigma <- 0.02^2 * (0.3 + 0.7 * diag(4))
#' d <- ph_data(tau, Sigma, cells = c("2004:2004", "2004:2005",
#'                                    "2006:2006", "2006:2007"))
#' d
ph_data.default <- function(x, Sigma, cells = NULL, weights = "equal",
                            n_units = NA_integer_, check_diagonal = TRUE,
                            ...) {
  tau <- as.numeric(x)
  K <- length(tau)
  if (K < 2L) {
    stop("Need at least 2 cohort-time cells; there is no partition problem ",
         "with fewer.", call. = FALSE)
  }
  if (anyNA(tau)) {
    stop("`x` contains missing values. Cohort-time cells that were not ",
         "estimated should be dropped before calling ph_data().", call. = FALSE)
  }

  if (missing(Sigma)) {
    stop("`Sigma`, the joint covariance of the cohort-time effects, is ",
         "required.\n  If you only have standard errors, pass them as a ",
         "vector -- but note that\n  the independence this implies is ",
         "usually wrong and understates aggregate\n  uncertainty (Remark 1).",
         call. = FALSE)
  }
  if (is.null(dim(Sigma))) {
    if (length(Sigma) != K) {
      stop(sprintf("`Sigma` given as a vector of standard errors must have ",
                   "length %d, not %d.", K, length(Sigma)), call. = FALSE)
    }
    Sigma <- diag(as.numeric(Sigma)^2, nrow = K)
  }
  Sigma <- as.matrix(Sigma)
  if (!identical(dim(Sigma), c(K, K))) {
    stop(sprintf("`Sigma` must be %d x %d to match the %d effects in `x`.",
                 K, K, K), call. = FALSE)
  }

  new_ph_data(tau = tau, Sigma = Sigma, cells = cells, weights = weights,
              n_units = n_units, source = "supplied",
              warn_diagonal = isTRUE(check_diagonal))
}

#' Internal constructor
#'
#' Validates and completes the fields shared by all three entry points.
#'
#' @noRd
new_ph_data <- function(tau, Sigma, cells, weights, n_units,
                        source = "supplied", G = NULL, sigma2 = NA_real_,
                        rss_flex = NA_real_, nobs = NA_integer_,
                        df_resid = NA_integer_, warn_diagonal = TRUE) {
  K <- length(tau)

  Sigma <- (Sigma + t(Sigma)) / 2
  Omega <- safe_inverse(Sigma, "first-stage covariance `Sigma`")

  if (warn_diagonal && K > 1L && is_diagonal(Sigma)) {
    warning(
      "`Sigma` is diagonal, so the cohort-time effects are being treated as ",
      "independent.\n  Cells that share units, cohorts or calendar periods ",
      "are generally correlated.\n  Remark 1 of the paper notes that ignoring ",
      "this understates the posterior\n  variance of aggregates such as the ",
      "overall ATT (by ~3.5x in the paper's\n  design) and can misstate ",
      "co-clustering probabilities by up to 0.39.\n  Supply the exact joint ",
      "covariance where you can.",
      call. = FALSE)
  }

  cells <- normalise_cells(cells, K)
  weights <- normalise_weights(weights, K)

  structure(
    list(
      tau = tau,
      Sigma = Sigma,
      Omega = Omega,
      tOt = as.numeric(crossprod(tau, Omega %*% tau)),
      cells = cells,
      weights = weights,
      K = K,
      n_units = n_units,
      source = source,
      # Micro-panel book-keeping; NA for the two-stage entry points. `G` is the
      # within cross-product D~'D~, so that Omega = G / sigma2 and
      # RSS(P) = rss_flex + sigma2 * deviance(P).
      G = G,
      sigma2 = sigma2,
      rss_flex = rss_flex,
      nobs = nobs,
      df_resid = df_resid
    ),
    class = "ph_data"
  )
}

#' @noRd
normalise_cells <- function(cells, K) {
  if (is.null(cells)) {
    return(data.frame(label = paste0("cell", seq_len(K)),
                      g = NA_real_, t = NA_real_,
                      stringsAsFactors = FALSE))
  }
  if (is.data.frame(cells)) {
    if (nrow(cells) != K) {
      stop(sprintf("`cells` has %d rows but there are %d effects.",
                   nrow(cells), K), call. = FALSE)
    }
    g <- if ("g" %in% names(cells)) cells$g else NA_real_
    t <- if ("t" %in% names(cells)) cells$t else NA_real_
    label <- if ("label" %in% names(cells)) {
      as.character(cells$label)
    } else if (!all(is.na(g)) && !all(is.na(t))) {
      paste0(g, ":", t)
    } else {
      paste0("cell", seq_len(K))
    }
    return(data.frame(label = label, g = g, t = t, stringsAsFactors = FALSE))
  }
  cells <- as.character(cells)
  if (length(cells) != K) {
    stop(sprintf("`cells` has length %d but there are %d effects.",
                 length(cells), K), call. = FALSE)
  }
  parts <- strsplit(cells, ":", fixed = TRUE)
  ok <- lengths(parts) == 2L
  g <- t <- rep(NA_real_, K)
  if (all(ok)) {
    g <- suppressWarnings(as.numeric(vapply(parts, `[`, character(1), 1L)))
    t <- suppressWarnings(as.numeric(vapply(parts, `[`, character(1), 2L)))
  }
  data.frame(label = cells, g = g, t = t, stringsAsFactors = FALSE)
}

#' @noRd
normalise_weights <- function(weights, K) {
  if (is.character(weights)) {
    weights <- match.arg(weights, c("equal"))
    return(rep(1 / K, K))
  }
  weights <- as.numeric(weights)
  if (length(weights) != K) {
    stop(sprintf("`weights` must have length %d, not %d.", K, length(weights)),
         call. = FALSE)
  }
  if (anyNA(weights)) stop("`weights` contains missing values.", call. = FALSE)
  total <- sum(weights)
  if (isTRUE(all.equal(total, 0))) {
    stop("`weights` sum to zero and cannot be normalised.", call. = FALSE)
  }
  weights / total
}

#' @export
print.ph_data <- function(x, ...) {
  cat("<ph_data>  first-stage cohort-time effects for a partition analysis\n\n")
  cat(sprintf("  Cells (K)      : %d\n", x$K))
  cat(sprintf("  Source         : %s\n", switch(x$source,
    supplied = "supplied directly",
    did = "did::att_gt() influence functions",
    panel = "flexible TWFE fitted on a micro panel",
    x$source)))
  if (!is.na(x$n_units)) cat(sprintf("  Units          : %d\n", x$n_units))
  se <- sqrt(diag(x$Sigma))
  cat(sprintf("  Covariance     : %s\n",
              if (is_diagonal(x$Sigma)) "diagonal (independence assumed)"
              else "full (cross-cell correlation retained)"))
  cat(sprintf("  sd(tau) / median(se) : %.2f", stats::sd(x$tau) / stats::median(se)))
  cat("   <- descriptive signal-to-noise\n\n")

  tab <- data.frame(
    cell = x$cells$label,
    estimate = round(x$tau, 4),
    se = round(se, 4),
    weight = round(x$weights, 4)
  )
  print(utils::head(tab, 12L), row.names = FALSE)
  if (x$K > 12L) cat(sprintf("  ... %d more cells\n", x$K - 12L))
  cat("\n  Use homogeneity_test() to ask whether there is heterogeneity to ",
      "recover,\n  then l0_ph() or bayes_ph() to recover it.\n", sep = "")
  invisible(x)
}
