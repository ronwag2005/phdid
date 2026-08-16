#' Is there heterogeneity to recover?
#'
#' Ask this before grouping anything. Both [l0_ph()] and [bayes_ph()] will
#' happily return a partition of pure noise: the paper's second application
#' finds the \eqn{\ell_0} estimator splitting 138 event effects into five
#' "bands" and the DP estimator reporting between 2.7 and 7.6 groups, when a
#' randomization test cannot reject a single common effect. At a signal-to-noise
#' ratio below one, such groups are quantiles of noise, not effect clusters.
#' This function is the guard against reading them as structure.
#'
#' Three assessments are reported, the first two always and the third when
#' placebo estimates are supplied.
#'
#' @section 1. Model-implied homogeneity test:
#' Under the null that all \eqn{K} cohort-time effects share a common value,
#' \eqn{\hat\tau \sim N(\theta 1_K, \hat\Sigma)}, so the GLS deviance of the
#' fully pooled fit is \eqn{\chi^2_{K-1}}. This is a genuine test with an exact
#' reference distribution and it uses the full cross-cell covariance, which is
#' what makes it trustworthy where a test built on independent standard errors
#' would not be (Remark 1). It is the natural test implied by the paper's own
#' model, and is reported here for every design.
#'
#' @section 2. Excess-dispersion heterogeneity share:
#' The cross-cell dispersion of the estimates mixes real heterogeneity with
#' sampling noise. Writing \eqn{M = I - 11'/K} for the centring projection, the
#' dispersion expected under a common effect is \eqn{\mathrm{tr}(M\hat\Sigma M)
#' / (K-1)}, so
#' \deqn{\hat\sigma^2_{\mathrm{het}} = \max\left\{0, \;
#'   \frac{\sum_k (\hat\tau_k - \bar\tau)^2}{K-1}
#'   - \frac{\mathrm{tr}(M \hat\Sigma M)}{K-1}\right\}}
#' estimates the genuine component, and its ratio to the observed dispersion is
#' the heterogeneity share. This is the covariance-based counterpart of the
#' paper's \eqn{\hat\sigma^2_{\mathrm{het}} = \max\{0, s^2_{\mathrm{post}} -
#' s^2_{\mathrm{pre}}\}}, using the estimator's own covariance in place of a
#' placebo as the noise gauge.
#'
#' @section 3. Placebo gauge and randomization test (Section 5.2.3):
#' When the design supplies a pre-treatment summary for each cell, constructed
#' symmetrically with the post-treatment one, no-anticipation makes it a direct
#' internal estimate of sampling noise. The paper's second application finds
#' \eqn{s_{\mathrm{pre}} = 0.0065} exceeding \eqn{s_{\mathrm{post}} = 0.0050},
#' a signal-to-noise ratio of 0.77 and a zero heterogeneity share.
#'
#' The randomization test then exchanges each cell's pre- and post-treatment
#' summaries. Because the exchange is internal to a cell, it holds the
#' cross-cell dependence structure fixed by construction, which is what makes
#' the test valid when the cells share controls. The observed post-treatment
#' dispersion is compared to the resulting reference distribution.
#'
#' The summaries are centred within each series before exchanging. Under the
#' null the effects share a *common* value, not a zero one, so the raw
#' post-treatment summaries are shifted by that common effect while the
#' pre-treatment ones are not; exchanging them unshifted would manufacture
#' dispersion and make the test conservative. Centring removes the shift and
#' leaves the deviations exchangeable, which is the null being tested. Set
#' `center = FALSE` for the literal unshifted exchange.
#'
#' @param object a [ph_data] object.
#' @param pre optional numeric vector of length K holding a pre-treatment
#'   (placebo) summary for each cell, constructed symmetrically with the
#'   post-treatment estimates.
#' @param nsim number of randomization draws.
#' @param center centre each series before exchanging; see above.
#' @param seed optional integer for reproducibility.
#'
#' @return An object of class `homogeneity_test`. Printing it gives the three
#'   assessments and a plain reading of which regime the design is in.
#'
#' @references
#' Arora, P. and Wagle, R. (2026). Partial Homogeneity in Staggered
#' Difference-in-Differences, Sections 5.1 and 5.2.3.
#'
#' @export
#' @examples
#' # Genuinely heterogeneous effects: the test rejects a common value.
#' d <- ph_data(c(0.10, 0.11, 0.42, 0.40), Sigma = 0.02^2 * (0.3 + 0.7 * diag(4)))
#' homogeneity_test(d)
#'
#' # Noise around a common value: it does not.
#' set.seed(1)
#' d0 <- ph_data(rnorm(8, 0.05, 0.02), Sigma = 0.02^2 * (0.3 + 0.7 * diag(8)))
#' homogeneity_test(d0)
homogeneity_test <- function(object, pre = NULL, nsim = 10000L, center = TRUE,
                             seed = NULL) {
  stopifnot(inherits(object, "ph_data"))
  if (!is.null(seed)) set.seed(seed)
  K <- object$K
  tau <- object$tau
  se <- sqrt(diag(object$Sigma))

  ## 1. Model-implied chi-square test on the pooled GLS deviance.
  pooled <- pooled_twfe(object)
  stat <- pooled$deviance
  df <- K - 1L
  p_chisq <- stats::pchisq(stat, df = df, lower.tail = FALSE)

  ## 2. Excess dispersion relative to what the covariance predicts.
  M <- diag(K) - matrix(1 / K, K, K)
  observed_disp <- sum((tau - mean(tau))^2) / (K - 1)
  expected_disp <- sum(diag(M %*% object$Sigma %*% M)) / (K - 1)
  het_var <- max(0, observed_disp - expected_disp)
  het_share <- if (observed_disp > 0) het_var / observed_disp else 0
  snr_model <- sqrt(observed_disp) / sqrt(expected_disp)

  ## Descriptive signal-to-noise of Section 5.1.
  snr_descriptive <- stats::sd(tau) / stats::median(se)

  ## 3. Placebo gauge and within-cell randomization test.
  placebo <- NULL
  if (!is.null(pre)) {
    pre <- as.numeric(pre)
    if (length(pre) != K) {
      stop(sprintf("`pre` must have length %d, one placebo summary per cell, ",
                   "not %d.", K, length(pre)), call. = FALSE)
    }
    if (anyNA(pre)) stop("`pre` contains missing values.", call. = FALSE)

    s_post <- stats::sd(tau)
    s_pre <- stats::sd(pre)
    het_var_pl <- max(0, s_post^2 - s_pre^2)
    het_share_pl <- if (s_post > 0) het_var_pl / s_post^2 else 0

    post_c <- if (center) tau - mean(tau) else tau
    pre_c <- if (center) pre - mean(pre) else pre

    ref <- numeric(nsim)
    for (b in seq_len(nsim)) {
      swap <- stats::runif(K) < 0.5
      drawn <- ifelse(swap, pre_c, post_c)
      ref[b] <- stats::sd(drawn)
    }
    observed <- stats::sd(post_c)
    # One-sided: how often does the reference produce dispersion at least as
    # large as observed? Small p means unusually dispersed, i.e. heterogeneity.
    p_rand <- (1 + sum(ref >= observed)) / (nsim + 1)

    placebo <- list(
      s_pre = s_pre, s_post = s_post,
      snr = s_post / s_pre,
      het_var = het_var_pl, het_share = het_share_pl,
      observed = observed, reference = ref,
      reference_median = stats::median(ref),
      p_value = p_rand, center = center, nsim = nsim
    )
  }

  structure(
    list(
      K = K,
      chisq = list(statistic = stat, df = df, p_value = p_chisq),
      dispersion = list(observed = observed_disp, expected = expected_disp,
                        het_var = het_var, het_share = het_share,
                        snr = snr_model),
      snr_descriptive = snr_descriptive,
      placebo = placebo,
      data = object
    ),
    class = "homogeneity_test"
  )
}

#' @export
print.homogeneity_test <- function(x, ...) {
  cat("<homogeneity_test>  is there heterogeneity to recover?\n\n")

  cat("1. Common-effect test (pooled GLS deviance, exact covariance)\n")
  cat(sprintf("     chi-squared = %.2f on %d df,  p = %s\n",
              x$chisq$statistic, x$chisq$df, format.pval(x$chisq$p_value,
                                                         digits = 3)))

  cat("\n2. Dispersion decomposition\n")
  cat(sprintf("     observed cross-cell sd    : %.5f\n",
              sqrt(x$dispersion$observed)))
  cat(sprintf("     expected under a common effect : %.5f\n",
              sqrt(x$dispersion$expected)))
  cat(sprintf("     heterogeneity sd (excess) : %.5f\n",
              sqrt(x$dispersion$het_var)))
  cat(sprintf("     heterogeneity share       : %.0f%%\n",
              100 * x$dispersion$het_share))
  cat(sprintf("     signal-to-noise           : %.2f\n", x$dispersion$snr))
  cat(sprintf("     sd(tau)/median(se)        : %.2f   (descriptive)\n",
              x$snr_descriptive))

  if (!is.null(x$placebo)) {
    p <- x$placebo
    cat("\n3. Placebo gauge and within-cell randomization test\n")
    cat(sprintf("     s_pre  (noise)            : %.5f\n", p$s_pre))
    cat(sprintf("     s_post (observed spread)  : %.5f\n", p$s_post))
    cat(sprintf("     signal-to-noise           : %.2f\n", p$snr))
    cat(sprintf("     heterogeneity share       : %.0f%%\n",
                100 * p$het_share))
    cat(sprintf("     randomization p (one-sided): %.3f   (%d draws%s)\n",
                p$p_value, p$nsim, if (p$center) ", centred" else ""))
    cat(sprintf("     observed sits %s the reference median\n",
                if (p$observed < p$reference_median) "below" else "above"))
  }

  cat("\n", strrep("-", 68), "\n", sep = "")
  cat(verdict_text(x))
  invisible(x)
}

#' Plain reading of which regime the design is in
#'
#' Deliberately conservative: the paper is emphatic that a failure to reject
#' does not establish a literally common effect, and that recovered groups at
#' low separation are noise quantiles.
#'
#' @noRd
verdict_text <- function(x) {
  p <- x$chisq$p_value
  share <- x$dispersion$het_share
  snr <- x$dispersion$snr
  p_rand <- if (!is.null(x$placebo)) x$placebo$p_value else NA_real_

  reject <- p < 0.05 && (is.na(p_rand) || p_rand < 0.10)

  if (!reject && share < 0.10) {
    return(paste0(
      "Reading: no recoverable heterogeneity.\n",
      "  A common effect cannot be rejected and the estimated heterogeneity ",
      "share is\n  near zero, so pooling is an adequate specification for the ",
      "aggregate. Groups\n  returned by l0_ph() or bayes_ph() here would be ",
      "quantiles of noise; do not\n  read them as effect clusters. Note this ",
      "is a failure to reject, not proof of\n  a literally common effect.\n"))
  }

  if (reject && snr >= 2) {
    return(paste0(
      "Reading: recoverable heterogeneity.\n",
      "  The spread across cells clearly exceeds what sampling noise would ",
      "produce, so\n  there is structure for the partition model to find. ",
      "Proceed to l0_ph() for a\n  point partition and bayes_ph() for ",
      "inference, and check the co-clustering\n  matrix to see which groupings ",
      "are firm.\n"))
  }

  paste0(
    "Reading: borderline, the low-separation regime.\n",
    "  There is some evidence of heterogeneity, but the effects are not ",
    "clearly\n  separated. This is the regime the paper identifies as the ",
    "boundary of reliable\n  recovery, where selection noise can erase or ",
    "even reverse the precision gain\n  and the two estimators may disagree. ",
    "Treat any single partition with caution;\n  report the co-clustering ",
    "matrix and the alpha path rather than one grouping.\n")
}
