#' @rdname ph_data
#'
#' @details
#' ## From a `did::att_gt()` fit
#'
#' The `MP` method keeps the post-treatment cells (\eqn{t \ge g}) and takes the
#' exact joint covariance from the estimator's influence functions, as
#' \eqn{\hat\Sigma = V_{\mathrm{analytical}} / n}. Because the cells share
#' units, this covariance is dense, which is precisely what Remark 1 asks for.
#' Fit with `bstrap = FALSE` so that `V_analytical` is populated.
#'
#' The default weights are cohort sizes, so the overall ATT is the
#' population-weighted average reported in the paper's first application.
#'
#' @export
ph_data.MP <- function(x, weights = c("cohort", "equal"), ...) {
  if (is.null(x$V_analytical)) {
    stop("This `att_gt()` fit has no analytical covariance.\n",
         "  Refit with `bstrap = FALSE` so that `V_analytical` is computed; ",
         "the partition\n  model needs the exact joint covariance of the ",
         "cohort-time effects (Remark 1).",
         call. = FALSE)
  }
  weights <- if (is.character(weights)) match.arg(weights) else weights

  keep <- x$t >= x$group
  if (sum(keep) < 2L) {
    stop("Fewer than two post-treatment cohort-time cells were estimated, ",
         "so there is no partition problem to solve.", call. = FALSE)
  }

  g <- x$group[keep]
  t <- x$t[keep]
  tau <- x$att[keep]
  Sigma <- (as.matrix(x$V_analytical) / x$n)[keep, keep, drop = FALSE]

  if (anyNA(tau) || anyNA(Sigma)) {
    stop("The `att_gt()` fit contains missing cohort-time effects or ",
         "covariances.\n  These cells are typically unidentified for want of ",
         "a clean comparison group\n  (Appendix A); drop them before ",
         "constructing the partition model.", call. = FALSE)
  }

  # `did` reports se = sqrt(diag(V_analytical) / n); confirm we have recovered
  # the same first stage before proceeding.
  se_gap <- max(abs(sqrt(diag(Sigma)) - x$se[keep]))
  if (se_gap > 1e-6) {
    warning(sprintf(
      paste0("Reconstructed standard errors differ from those reported by ",
             "`did` by up to %.2e.\n  The covariance may not correspond to ",
             "the reported effects."), se_gap), call. = FALSE)
  }

  w <- weights
  if (identical(w, "cohort")) {
    w <- cohort_weights_from_mp(x, g)
  }

  new_ph_data(
    tau = tau, Sigma = Sigma,
    cells = data.frame(label = paste0(g, ":", t), g = g, t = t,
                       stringsAsFactors = FALSE),
    weights = w,
    n_units = x$n,
    source = "did"
  )
}

#' Cohort-size weights for the cells of an `MP` fit
#'
#' Weight each cohort-time cell by the size of its cohort, so that the overall
#' ATT is population weighted. Falls back to equal weights, with a warning, if
#' the cohort counts cannot be recovered from the fit.
#'
#' @noRd
cohort_weights_from_mp <- function(x, g) {
  counts <- tryCatch(x$DIDparams$cohort_counts, error = function(e) NULL)
  sizes <- NULL
  if (!is.null(counts) && all(c("cohort", "cohort_size") %in% names(counts))) {
    counts <- as.data.frame(counts)
    sizes <- counts$cohort_size[match(g, counts$cohort)]
  }
  if (is.null(sizes) || anyNA(sizes)) {
    warning("Could not recover cohort sizes from the `att_gt()` fit; ",
            "falling back to equal weights.\n  Pass `weights` explicitly to ",
            "control the overall ATT aggregator.", call. = FALSE)
    return(rep(1 / length(g), length(g)))
  }
  as.numeric(sizes)
}


#' @rdname ph_data
#'
#' @details
#' ## From a micro panel
#'
#' The data frame method fits the fully flexible (interacted) TWFE model of
#' eq. (3),
#' \deqn{Y_{igt} = \alpha_i + \lambda_t + \sum_{g}\sum_{t \ge g}
#'   \tau_{gt} D_{igt} + \varepsilon_{it},}
#' and returns \eqn{\hat\tau_{\mathrm{flex}}} together with
#' \eqn{\hat\sigma^2 (\tilde D'\tilde D)^{-1}}. The unit and time fixed effects
#' (and any covariates) are partialled out by alternating projections, which
#' handles unbalanced panels; by the Frisch--Waugh--Lovell theorem the CATT
#' estimates are the same whether covariates enter the regression directly or
#' are residualised beforehand, which is eq. (8)--(9) of the paper.
#'
#' This route assumes spherical errors, since it estimates a single
#' \eqn{\hat\sigma^2}. For clustered or serially correlated errors, run a
#' first-stage estimator that reports a robust joint covariance and pass the
#' pair `(tau, Sigma)` to the default method instead --- Lemma 1 says the
#' downstream analysis is unchanged.
#'
#' @param yname,idname,tname,gname column names holding the outcome, the unit
#'   identifier, the time period, and the cohort (period of first treatment,
#'   with `0`, `Inf` or `NA` marking the never-treated).
#' @param xformla an optional one-sided formula of covariates to partial out,
#'   e.g. `~ lpop`.
#'
#' @export
ph_data.data.frame <- function(x, yname, idname, tname, gname,
                               xformla = NULL, weights = c("cohort", "equal"),
                               ...) {
  weights <- if (is.character(weights)) match.arg(weights) else weights
  need <- c(yname, idname, tname, gname)
  missing_cols <- setdiff(need, names(x))
  if (length(missing_cols)) {
    stop(sprintf("Column(s) not found in the data: %s.",
                 paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  d <- x[, need, drop = FALSE]
  names(d) <- c("y", "id", "t", "g")

  X <- NULL
  if (!is.null(xformla)) {
    mf <- stats::model.frame(xformla, data = x, na.action = stats::na.pass)
    X <- stats::model.matrix(xformla, mf)
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  }

  ok <- stats::complete.cases(d) & (if (is.null(X)) TRUE else stats::complete.cases(X))
  if (any(!ok)) {
    d <- d[ok, , drop = FALSE]
    if (!is.null(X)) X <- X[ok, , drop = FALSE]
    message(sprintf("Dropped %d rows with missing values.", sum(!ok)))
  }

  # Never-treated units may be coded 0, Inf or NA; normalise to Inf.
  d$g[is.na(d$g) | d$g == 0] <- Inf

  periods <- sort(unique(d$t))
  cohorts <- sort(unique(d$g[is.finite(d$g)]))
  # A cohort first treated in the first period has no pre-period, so it
  # contributes no identified cell.
  cohorts <- cohorts[cohorts > min(periods)]
  if (!length(cohorts)) {
    stop("No treated cohort has a pre-treatment period, so no cohort-time ",
         "effect is identified.", call. = FALSE)
  }

  cells <- do.call(rbind, lapply(cohorts, function(gg) {
    tt <- periods[periods >= gg]
    if (!length(tt)) return(NULL)
    data.frame(g = gg, t = tt)
  }))
  # Keep only cells that actually occur in the data.
  present <- mapply(function(gg, tt) any(d$g == gg & d$t == tt),
                    cells$g, cells$t)
  cells <- cells[present, , drop = FALSE]
  K <- nrow(cells)
  if (K < 2L) {
    stop("Fewer than two cohort-time cells are present in the data.",
         call. = FALSE)
  }

  NT <- nrow(d)
  if (as.numeric(NT) * K > 5e7) {
    stop(sprintf(
      paste0("This panel would need a dense %d x %d design matrix, which is ",
             "too large to\n  build in memory. Run a first-stage estimator ",
             "such as did::att_gt() or a\n  fixest interaction model, then ",
             "pass its estimates and joint covariance to\n  ph_data() ",
             "directly -- Lemma 1 guarantees the same answer."), NT, K),
      call. = FALSE)
  }

  # Build the raw cohort-time treatment dummies.
  D <- matrix(0, NT, K)
  for (k in seq_len(K)) {
    D[, k] <- as.numeric(d$g == cells$g[k] & d$t == cells$t[k])
  }

  # Partial out unit and time fixed effects (and covariates) from both the
  # outcome and the dummies -- eq. (9).
  M <- cbind(d$y, D, X)
  M <- demean_factors(M, list(d$id, d$t))
  ytil <- M[, 1L]
  Dtil <- M[, 1L + seq_len(K), drop = FALSE]
  p_cov <- 0L
  if (!is.null(X)) {
    Xtil <- M[, 1L + K + seq_len(ncol(X)), drop = FALSE]
    keep_x <- apply(Xtil, 2L, function(cc) stats::sd(cc) > 1e-12)
    Xtil <- Xtil[, keep_x, drop = FALSE]
    if (ncol(Xtil)) {
      qx <- qr(Xtil)
      p_cov <- qx$rank
      ytil <- qr.resid(qx, ytil)
      Dtil <- qr.resid(qx, Dtil)
    }
  }

  G <- crossprod(Dtil)
  cvec <- crossprod(Dtil, ytil)
  qrG <- qr(G)
  if (qrG$rank < K) {
    stop(sprintf(
      paste0("The flexible design is rank deficient (rank %d of %d cells), so ",
             "at least one\n  cohort-time effect is not identified -- ",
             "typically a cell with no clean\n  comparison group (Appendix A ",
             "on limited overlap). Drop the affected cells,\n  or supply a ",
             "partition that lets them inherit a group effect."),
      qrG$rank, K), call. = FALSE)
  }

  tau <- as.numeric(solve(qrG, cvec))
  yty <- sum(ytil^2)
  rss_flex <- yty - sum(cvec * tau)

  n_units <- length(unique(d$id))
  n_periods <- length(periods)
  df_resid <- NT - n_units - n_periods + 1L - K - p_cov
  if (df_resid <= 0L) {
    stop("No residual degrees of freedom remain after the fixed effects and ",
         "cohort-time\n  effects; the flexible model is saturated.",
         call. = FALSE)
  }
  sigma2 <- rss_flex / df_resid
  Sigma <- sigma2 * safe_inverse(G, "within cross-product D~'D~")

  w <- weights
  if (identical(w, "cohort")) {
    sizes <- vapply(cells$g, function(gg) {
      length(unique(d$id[d$g == gg]))
    }, numeric(1))
    w <- sizes
  }

  new_ph_data(
    tau = tau, Sigma = Sigma,
    cells = data.frame(label = paste0(cells$g, ":", cells$t),
                       g = cells$g, t = cells$t, stringsAsFactors = FALSE),
    weights = w,
    n_units = n_units,
    source = "panel",
    G = G, sigma2 = sigma2, rss_flex = rss_flex,
    nobs = NT, df_resid = df_resid,
    # A balanced orthogonal design legitimately produces a near-diagonal
    # covariance here, so the Remark 1 warning would be noise.
    warn_diagonal = FALSE
  )
}


#' Residualise columns on a set of factors by alternating projections
#'
#' Sweeps out each factor's group means in turn until the columns stop
#' changing. For a balanced two-way panel this converges to the familiar
#' double-demeaning \eqn{Y - \bar Y_i - \bar Y_t + \bar Y} in a couple of
#' passes; for unbalanced panels the iteration is what makes the projection
#' exact.
#'
#' @param M numeric matrix whose columns are to be residualised.
#' @param factors list of grouping vectors, each of length `nrow(M)`.
#' @param tol,maxit convergence tolerance and iteration cap.
#' @return `M` with the factor projections removed.
#' @noRd
demean_factors <- function(M, factors, tol = 1e-10, maxit = 500L) {
  M <- as.matrix(M)
  codes <- lapply(factors, function(f) as.integer(factor(f)))
  counts <- lapply(codes, function(cc) tabulate(cc))

  scale <- max(abs(M))
  if (scale == 0) return(M)

  for (it in seq_len(maxit)) {
    delta <- 0
    for (j in seq_along(codes)) {
      cc <- codes[[j]]
      means <- rowsum(M, cc, reorder = TRUE) / counts[[j]]
      step <- means[cc, , drop = FALSE]
      M <- M - step
      delta <- max(delta, max(abs(step)))
    }
    if (delta / scale < tol) break
  }
  if (it == maxit) {
    warning("Alternating projections did not fully converge when removing the ",
            "fixed effects.\n  This can happen when the panel is not ",
            "connected across units and periods.", call. = FALSE)
  }
  M
}
