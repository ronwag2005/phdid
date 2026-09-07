#' @keywords internal
#' @aliases phdid-package
#'
#' @importFrom stats aggregate coef confint vcov
#'
#' @description
#' In a staggered difference-in-differences design the treatment effect is not
#' a single number but a vector of cohort-time effects (CATTs), one per
#' cohort-time cell. Estimating every one separately is unbiased but
#' inefficient when some are in fact equal; pooling them all into a single TWFE
#' coefficient is efficient but biased whenever the heterogeneity is genuine.
#' \pkg{phdid} treats the choice between these extremes as a *partition
#' selection* problem on the cells and provides two answers, along with the
#' tests and diagnostics needed to know whether either should be used at all.
#'
#' @section Where to start:
#' \describe{
#'   \item{[ph_data()]}{Assemble the inputs. Accepts a `did::att_gt()` fit, a
#'     micro panel, or any vector of first-stage effects with their joint
#'     covariance.}
#'   \item{[homogeneity_test()]}{Ask whether there is heterogeneity to recover
#'     before recovering any. Both estimators will partition pure noise if
#'     asked to.}
#'   \item{[bayes_ph()]}{The Dirichlet Process estimator. Averages over
#'     partitions, so its intervals include uncertainty about the grouping.
#'     This is the recommended route for inference.}
#'   \item{[l0_ph()]}{The \eqn{\ell_0}-penalised estimator, which returns a
#'     single interpretable partition. Its intervals condition on that
#'     partition being right.}
#'   \item{[aggregate.bayes_ph()]}{Overall ATT, event study, or effects by
#'     cohort, computed per posterior draw.}
#' }
#'
#' @section Two things the package will keep reminding you of:
#' First, the joint covariance of the first-stage effects is not diagonal, and
#' using the exact one is what makes the reported intervals honest. Second, a
#' partition is only worth reporting when the effects are separated enough to
#' be recovered; below that threshold the estimators return groupings that are
#' quantiles of noise, and the package says so rather than letting the output
#' speak for itself.
#'
#' @section Authors:
#' The methods implemented here are joint work by Parush Arora and Rohan Wagle.
#' See `citation("phdid")`.
#'
#' @references
#' Arora, P. and Wagle, R. (2026). A Bayesian Approach to Partial Homogeneity
#' in Staggered Difference-in-Difference. Ashoka University Economics
#' Discussion Paper 166.
#' \url{https://www.ashoka.edu.in/research/a-bayesian-approach-to-partial-homogeneity-in-staggered-difference-in-difference/}
#'
#' Arora, P. and Wagle, R. (2026). A Bayesian Approach to Partial Homogeneity
#' in Staggered Difference-in-Difference. SSRN Working Paper 7207083
#' (1 March 2026). \doi{10.2139/ssrn.7207083}
#'
#' Callaway, B. and Sant'Anna, P. H. C. (2021). Difference-in-Differences with
#' Multiple Time Periods. \emph{Journal of Econometrics} 225(2), 200--230.
#'
#' Neal, R. M. (2000). Markov Chain Sampling Methods for Dirichlet Process
#' Mixture Models. \emph{JCGS} 9(2), 249--265.
#'
#' Wooldridge, J. M. (2025). Two-Way Fixed Effects, the Two-Way Mundlak
#' Regression, and Difference-in-Differences Estimators. \emph{Empirical
#' Economics} 69, 2545--2587.
"_PACKAGE"
