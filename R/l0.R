#' The l0-penalised partial-homogeneity estimator
#'
#' Recovers a partition of the cohort-time cells by minimising
#' \deqn{Q(\mathcal{P}) = \mathrm{RSS}(\mathcal{P}) + \lambda \, c(\mathcal{P}),}
#' where \eqn{c(\mathcal{P})} counts the pairs of cells placed in different
#' groups (eq. 34--35). Proposition 2 of the paper shows this is not an
#' unrelated penalised estimator but the fixed-variance maximum a posteriori
#' partition of the same Bayesian model that [bayes_ph()] samples from,
#' differing only in the partition prior --- a pairwise prior here, the Chinese
#' Restaurant Process there. The \eqn{\ell_0} form is used rather than an
#' \eqn{\ell_1} fusion penalty because it produces exact equality of
#' coefficients rather than continuous shrinkage.
#'
#' @section The two steps:
#' Exact minimisation over partitions would require searching all \eqn{B_K}
#' of them, so the estimator is computed in two steps (Appendix B):
#'
#' 1. **Search.** A greedy agglomerative pass starting from \eqn{K} singletons.
#'    At each stage the pair of groups minimising
#'    \deqn{\Delta\mathrm{Obj}(A,B) = \frac{n_A n_B}{n_A + n_B}
#'      (\hat\tau_A - \hat\tau_B)^2 - \lambda |A| |B|}
#'    is merged. This runs in \eqn{O(K^3)} time.
#' 2. **Re-estimation.** The grouped model is re-fitted by GLS under the
#'    selected partition, [ph_fit()].
#'
#' Step 2 is not a refinement. The merge criterion in step 1 treats the design
#' as orthogonal, and under orthogonality the greedy group means coincide with
#' the GLS fit (eq. 49). In a general panel they do not: cells sharing a cohort
#' or a period leave non-zero off-diagonal entries, the greedy means fail the
#' normal equations, and only the re-estimated fit is best linear unbiased.
#'
#' @section Choosing the penalty:
#' With `select = "bic"` (the default) the estimator sweeps the whole
#' agglomeration path and picks the number of groups minimising the BIC, which
#' needs no tuning grid. Proposition 3 reads the threshold as a rate: a correct
#' merge costs \eqn{O_p(\sigma^2)} while a wrong merge costs
#' \eqn{O(N \Delta\tau^2)}, so any \eqn{\lambda} between the two separates them,
#' and BIC's effective threshold \eqn{\lambda \asymp \sigma^2 \log(NT)} sits in
#' that window. Passing a numeric `lambda` instead applies the stopping rule of
#' Appendix B directly.
#'
#' @section What it cannot do:
#' The recovery window narrows as the distinct effects get closer together.
#' The paper's simulations show that at a separation of three standard errors
#' the partition is recovered only about half the time and the selection noise
#' erases the precision gain entirely. Run [homogeneity_test()] first: if the
#' effects carry no recoverable heterogeneity, the groups this function returns
#' are quantiles of noise, not effect clusters.
#'
#' @param object a [ph_data] object.
#' @param select how to choose the point on the agglomeration path. `"bic"`
#'   (default) minimises the BIC; `"lambda"` applies the Appendix B stopping
#'   rule at the supplied `lambda`; `"m"` takes a fixed number of groups.
#' @param lambda the \eqn{\ell_0} penalty, on the scale of eq. (36): a pair of
#'   singleton cells is merged when their squared difference, weighted by the
#'   harmonic mean of their effective sample sizes, falls below `lambda`.
#'   Required when `select = "lambda"`.
#' @param m the number of groups, when `select = "m"`.
#' @param search `"greedy"` for the paper's Appendix B algorithm, whose merge
#'   costs use the diagonal (orthogonal) approximation. `"exact-cost"` is an
#'   extension, not in the paper: it scores each candidate merge by its true
#'   GLS deviance increment under the full covariance. It is still a greedy
#'   search and so still a heuristic, but it removes the orthogonality
#'   approximation from step 1. Costs \eqn{O(K^5)}; use only for small `K`.
#' @param n_bic sample size in the BIC penalty on the two-stage route. See
#'   details under `partition_bic`; the default is `K`, which reproduces the
#'   paper's reported group counts.
#'
#' @return An object of class `l0_ph`, extending [ph_fit], with the additional
#'   components `path` (a data frame of deviance, RSS and BIC at every point on
#'   the agglomeration path), `partitions` (the partition at each `m`),
#'   `selected` (the chosen `m`) and `lambda`.
#'
#' @references
#' Arora, P. and Wagle, R. (2026). Partial Homogeneity in Staggered
#' Difference-in-Differences. Propositions 2 and 3, and Appendix B.
#'
#' @seealso [bayes_ph()], which averages over partitions instead of committing
#'   to one, and [homogeneity_test()], which asks whether to group at all.
#'
#' @export
#' @examples
#' # Four cells, two well-separated effect levels.
#' d <- ph_data(c(0.10, 0.11, 0.42, 0.40), Sigma = 0.02^2 * (0.3 + 0.7 * diag(4)))
#' fit <- l0_ph(d)
#' fit
#'
#' # The whole agglomeration path, which is what BIC is selecting over.
#' fit$path
l0_ph <- function(object, select = c("bic", "lambda", "m"), lambda = NULL,
                  m = NULL, search = c("greedy", "exact-cost"),
                  n_bic = NULL) {
  stopifnot(inherits(object, "ph_data"))
  select <- match.arg(select)
  search <- match.arg(search)
  K <- object$K

  if (select == "lambda" && is.null(lambda)) {
    stop("`select = \"lambda\"` needs a numeric `lambda`.\n",
         "  Use `select = \"bic\"` to choose the penalty automatically.",
         call. = FALSE)
  }
  if (select == "m") {
    if (is.null(m)) stop("`select = \"m\"` needs a number of groups `m`.",
                         call. = FALSE)
    if (m < 1L || m > K) {
      stop(sprintf("`m` must be between 1 and %d.", K), call. = FALSE)
    }
  }
  if (search == "exact-cost" && K > 40L) {
    warning(sprintf(
      "`search = \"exact-cost\"` is O(K^5) and K is %d; this may take a very ",
      "long time.\n  Consider the default `search = \"greedy\"`.", K),
      call. = FALSE)
  }

  path <- agglomeration_path(object, search = search)

  # Step 2: exact GLS re-estimation at every point on the path, so that the
  # criteria are computed from the true fit rather than the greedy means.
  fits <- lapply(path$partitions, function(z) ph_fit(object, z))
  deviance <- vapply(fits, `[[`, numeric(1), "deviance")
  rss <- vapply(fits, `[[`, numeric(1), "rss")
  bic <- vapply(seq_len(K), function(i) {
    partition_bic(object, deviance[i], i, n_bic = n_bic)
  }, numeric(1))

  path_df <- data.frame(
    m = seq_len(K),
    deviance = deviance,
    rss = rss,
    bic = bic,
    cross_group_pairs = vapply(path$partitions, cross_group_pairs, numeric(1)),
    merge_cost = path$merge_cost
  )

  selected <- switch(
    select,
    bic = which.min(bic),
    m = as.integer(m),
    lambda = select_by_lambda(path, lambda, K)
  )

  fit <- fits[[selected]]
  fit$path <- path_df
  fit$partitions <- path$partitions
  fit$selected <- selected
  fit$lambda <- lambda
  fit$select <- select
  fit$search <- search
  fit$n_bic <- if (is.null(n_bic)) object$K else n_bic
  fit$effective_size <- path$effective_size
  class(fit) <- c("l0_ph", "ph_fit")
  fit
}

#' Apply the Appendix B stopping rule at a fixed lambda
#'
#' Merging continues while some pair has \eqn{\Delta\mathrm{Obj} < 0}, i.e.
#' while the recorded merge cost falls below \eqn{\lambda |A||B|}. The path
#' records the cost and the pair sizes at each step, so the rule is read off
#' the path rather than re-run.
#'
#' @noRd
select_by_lambda <- function(path, lambda, K) {
  # path$merge_cost[i] is the cost of the merge that produced m = i groups,
  # and path$pair_product[i] the |A||B| of that merge. Entry K is NA (no merge).
  m <- K
  for (i in seq(K - 1L, 1L)) {
    if (is.na(path$merge_cost[i])) break
    if (path$merge_cost[i] < lambda * path$pair_product[i]) {
      m <- i
    } else {
      break
    }
  }
  m
}

#' Greedy agglomerative merge path
#'
#' Builds the full sequence of partitions from K singletons down to one group,
#' recording the cost of each merge. Returns partitions indexed by the number
#' of groups, so `partitions[[m]]` is the partition with `m` groups.
#'
#' The cell weights \eqn{n_k} are the diagonal of the precision matrix
#' \eqn{\hat\Sigma^{-1}}, rescaled by \eqn{\hat\sigma^2} when it is known so
#' that they are the effective sample sizes \eqn{\|\tilde D_k\|^2} of the
#' paper. On the two-stage route no such scale exists and the precision itself
#' is used, which is what the paper's application code does.
#'
#' @noRd
agglomeration_path <- function(object, search = "greedy") {
  K <- object$K
  prec <- diag(object$Omega)
  if (!is.na(object$sigma2)) prec <- prec * object$sigma2

  groups <- as.list(seq_len(K))
  g_tau <- object$tau
  g_n <- prec
  g_size <- rep(1L, K)

  partitions <- vector("list", K)
  partitions[[K]] <- seq_len(K)
  merge_cost <- rep(NA_real_, K)
  pair_product <- rep(NA_real_, K)

  for (step in seq_len(K - 1L)) {
    ng <- length(groups)
    best <- Inf
    bi <- bj <- NA_integer_

    if (search == "greedy") {
      for (a in seq_len(ng - 1L)) {
        for (b in seq(a + 1L, ng)) {
          cost <- (g_n[a] * g_n[b] / (g_n[a] + g_n[b])) *
            (g_tau[a] - g_tau[b])^2
          if (cost < best) { best <- cost; bi <- a; bj <- b }
        }
      }
    } else {
      # Exact GLS deviance increment for each candidate merge. Extension to
      # the paper's algorithm; removes the orthogonality approximation from
      # the search, at O(K^2) refits per step.
      current <- ph_fit(object, labels_from_groups(groups, K))$deviance
      for (a in seq_len(ng - 1L)) {
        for (b in seq(a + 1L, ng)) {
          trial <- groups
          trial[[a]] <- c(trial[[a]], trial[[b]])
          trial[[b]] <- NULL
          cost <- ph_fit(object, labels_from_groups(trial, K))$deviance - current
          if (cost < best) { best <- cost; bi <- a; bj <- b }
        }
      }
      if (!is.na(object$sigma2)) best <- best * object$sigma2
    }

    merge_cost[ng - 1L] <- best
    pair_product[ng - 1L] <- g_size[bi] * g_size[bj]

    new_n <- g_n[bi] + g_n[bj]
    g_tau[bi] <- (g_n[bi] * g_tau[bi] + g_n[bj] * g_tau[bj]) / new_n
    g_n[bi] <- new_n
    g_size[bi] <- g_size[bi] + g_size[bj]
    groups[[bi]] <- c(groups[[bi]], groups[[bj]])

    groups[[bj]] <- NULL
    g_tau <- g_tau[-bj]
    g_n <- g_n[-bj]
    g_size <- g_size[-bj]

    partitions[[ng - 1L]] <- labels_from_groups(groups, K)
  }

  list(partitions = partitions, merge_cost = merge_cost,
       pair_product = pair_product, effective_size = prec)
}

#' @noRd
labels_from_groups <- function(groups, K) {
  z <- integer(K)
  for (p in seq_along(groups)) z[groups[[p]]] <- p
  relabel(z)
}


#' @export
print.l0_ph <- function(x, ...) {
  cat("<l0_ph>  l0-penalised partial-homogeneity estimator\n\n")
  cat(sprintf("  Cells      : %d\n", length(x$tau)))
  cat(sprintf("  Groups (m) : %d", x$m))
  cat(sprintf("   (selected by %s)\n", switch(x$select,
    bic = sprintf("BIC, n = %g in the penalty", x$n_bic),
    lambda = sprintf("the Appendix B rule at lambda = %g", x$lambda),
    m = "hand")))
  cat(sprintf("  Search     : %s\n", switch(x$search,
    greedy = "greedy agglomerative (Appendix B)",
    `exact-cost` = "greedy with exact GLS merge costs (extension)")))
  cat(sprintf("  Partition  : %s\n\n", format_partition(x$partition,
                                                        x$data$cells$label)))

  se_flex <- sqrt(diag(x$data$Sigma))
  se <- sqrt(diag(x$vcov))
  ord <- order(x$partition, x$tau)
  tab <- data.frame(
    cell = x$data$cells$label,
    flexible = round(x$data$tau, 4),
    se = round(se_flex, 4),
    group = x$partition,
    grouped = round(x$tau, 4),
    var_ratio = round((se / se_flex)^2, 2)
  )[ord, ]
  print(utils::head(tab, 20L), row.names = FALSE)
  if (length(x$tau) > 20L) {
    cat(sprintf("  ... %d more cells\n", length(x$tau) - 20L))
  }

  pooled_cells <- tabulate(x$partition)[x$partition] > 1L
  if (any(pooled_cells)) {
    cat(sprintf("\n  Mean variance ratio among pooled cells: %.2f\n",
                mean(((se / se_flex)^2)[pooled_cells])))
  }
  cat("\n  These intervals condition on the selected partition being ",
      "correct.\n  Remark 3: they can under-cover, badly so when the ",
      "partition is uncertain\n  (0.57-0.62 in the paper's hard regimes). ",
      "Use bayes_ph() for inference.\n", sep = "")
  invisible(x)
}
