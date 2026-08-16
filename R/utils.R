## Internal helpers. Not exported.
##
## Throughout, a *partition* of the K cohort-time cells is represented as an
## integer vector `z` of length K taking values in 1..m, canonicalised so that
## group labels appear in order of first occurrence (see `relabel()`). This is
## the representation used by the Gibbs sampler and by the agglomerative search.

#' Canonicalise a partition label vector
#'
#' Relabels so that groups are numbered 1..m in order of first appearance,
#' which makes two encodings of the same partition compare equal.
#'
#' @param z integer vector of group labels.
#' @return integer vector of the same length with labels in 1..m.
#' @noRd
relabel <- function(z) {
  match(z, unique(z))
}

#' Group-indicator (membership) matrix
#'
#' Builds the \eqn{K \times m} matrix \eqn{R} of the paper's Section 5, with
#' \eqn{R_{kp} = 1} when cell \eqn{k} belongs to group \eqn{p}. Partial
#' homogeneity is the restriction \eqn{\tau = R\phi}.
#'
#' @param z integer vector of group labels, length K.
#' @return a K x m indicator matrix.
#' @noRd
indicator_matrix <- function(z) {
  z <- relabel(z)
  m <- max(z)
  R <- matrix(0, length(z), m)
  R[cbind(seq_along(z), z)] <- 1
  R
}

#' Number of cross-group pairs, c(P)
#'
#' The \eqn{\ell_0} penalty's complexity term of eq. (33):
#' \eqn{c(\mathcal{P}) = \binom{K}{2} - \sum_p \binom{|C_p|}{2}}.
#'
#' @param z integer vector of group labels.
#' @return the number of pairs of cells assigned to different groups.
#' @noRd
cross_group_pairs <- function(z) {
  sizes <- tabulate(relabel(z))
  choose2 <- function(x) x * (x - 1) / 2
  choose2(length(z)) - sum(choose2(sizes))
}

#' Symmetric solve with a informative error
#'
#' Wraps `chol2inv(chol(.))` so that a non-positive-definite covariance or
#' precision produces a message naming the likely cause rather than a bare
#' LAPACK error.
#'
#' @param A a symmetric positive-definite matrix.
#' @param what a short name used in the error message.
#' @return the inverse of `A`.
#' @noRd
safe_inverse <- function(A, what = "matrix") {
  A <- (A + t(A)) / 2
  ch <- tryCatch(chol(A), error = function(e) NULL)
  if (is.null(ch)) {
    stop(sprintf(
      paste0("The %s is not positive definite, so it cannot be inverted.\n",
             "  This usually means the first-stage covariance is singular, ",
             "which happens when\n  a cohort-time cell has no clean comparison ",
             "and so carries no independent\n  identifying variation (see the ",
             "limited-overlap discussion in Appendix A).\n",
             "  Drop the affected cell, or supply a regularised covariance."),
      what), call. = FALSE)
  }
  chol2inv(ch)
}

#' Log-determinant of a symmetric positive-definite matrix
#' @noRd
log_det_spd <- function(A) {
  A <- (A + t(A)) / 2
  2 * sum(log(diag(chol(A))))
}

#' Adjusted Rand index between two partitions
#'
#' Used to score recovery of the true partition in the simulations
#' (paper Tables 2 and 3).
#'
#' @param a,b integer vectors of group labels of equal length.
#' @return the adjusted Rand index; 1 for identical partitions, 0 in
#'   expectation for independent ones.
#' @export
#' @examples
#' adjusted_rand(c(1, 1, 2, 2), c(2, 2, 1, 1))  # same partition, relabelled
#' adjusted_rand(c(1, 1, 2, 2), c(1, 2, 1, 2))
adjusted_rand <- function(a, b) {
  stopifnot(length(a) == length(b))
  tab <- table(a, b)
  n <- sum(tab)
  choose2 <- function(x) x * (x - 1) / 2
  si <- sum(choose2(tab))
  sa <- sum(choose2(rowSums(tab)))
  sb <- sum(choose2(colSums(tab)))
  expected <- sa * sb / choose2(n)
  maximum <- (sa + sb) / 2
  if (isTRUE(all.equal(maximum, expected))) return(1)
  (si - expected) / (maximum - expected)
}

#' Is a matrix diagonal (up to tolerance)?
#'
#' Used to raise the Remark 1 warning: treating the cohort-time estimates as
#' independent understates the posterior variance of linear aggregates and
#' misstates co-clustering probabilities.
#'
#' @noRd
is_diagonal <- function(A, tol = 1e-10) {
  off <- A - diag(diag(A), nrow = nrow(A))
  scale <- max(abs(diag(A)))
  if (scale == 0) return(TRUE)
  max(abs(off)) / scale < tol
}

#' Format a partition compactly for printing, e.g. "{1,3} {2} {4,5}"
#' @noRd
format_partition <- function(z, labels = NULL) {
  z <- relabel(z)
  if (is.null(labels)) labels <- as.character(seq_along(z))
  groups <- split(labels, z)
  paste(vapply(groups, function(g) paste0("{", paste(g, collapse = ","), "}"),
               character(1)), collapse = " ")
}
