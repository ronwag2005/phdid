# phdid 0.1.0

Released on CRAN 2026-09-15. First release, implementing Arora and Wagle
(2026), "A Bayesian Approach to Partial Homogeneity in Staggered
Difference-in-Difference" (doi:10.2139/ssrn.7207083).

## Estimation

* `ph_data()` assembles the inputs from a `did::att_gt()` fit, a micro panel,
  or any vector of first-stage effects with their joint covariance. Lemma 1 of
  the paper guarantees the three routes give identical answers, and the test
  suite asserts it.
* `bayes_ph()` fits the Dirichlet Process mixture by collapsed Gibbs sampling
  (Neal 2000, Algorithm 3), with the exact collapsed marginal likelihood in the
  cluster-assignment moves by default.
* `l0_ph()` computes the two-step l0 estimator of Appendix B, supporting both
  BIC selection along the agglomeration path and the fixed-`lambda` stopping
  rule.
* `ph_fit()`, `flex_twfe()` and `pooled_twfe()` provide the restricted GLS fit
  under a known partition and the two benchmark corners.

## Inference

* `aggregate()` applies overall, event-study, by-cohort and by-calendar
  aggregators to each posterior draw, so reported intervals propagate
  uncertainty about the partition itself.
* `coclustering()` and `point_partition()` summarise the grouping structure,
  the latter by Wade-Ghahramani VI or Binder loss.

## Specification testing

* `homogeneity_test()` reports a model-implied common-effect chi^2 test on the
  pooled GLS deviance, an excess-dispersion heterogeneity share computed
  against the exact covariance, and -- when placebo estimates are supplied -- the
  pre/post noise gauge and within-cell randomization test of Section 5.2.3.

## Diagnostics

* `enumerate_partitions()` computes the exact partition posterior for small
  designs, the Appendix E benchmark.
* `ph_rhat()` reports Gelman-Rubin statistics and effective sample sizes across
  dispersed chains.
* `covariance_check()` compares exact against diagonal handling of the
  covariance in the cluster-assignment moves.

## Regularisation paths

* `alpha_sensitivity()` traces the DP estimator across the concentration
  parameter, and `lambda_sensitivity()` traces the `l0` estimator across its
  penalty. Both take `type = "cells"`, reporting the path for every cohort-time
  effect separately rather than only for the aggregate. This is usually the
  more informative view: the overall ATT is robust to over-pooling, so a flat
  aggregate path can hide cells that move a great deal.
* `plot_sensitivity()` draws either table; a per-cell table gets one labelled
  line per cohort-time effect.
* `lambda_sensitivity()` can index the `l0` path by the penalty (`by =
  "lambda"`) or by the number of groups (`by = "m"`). These are different
  paths: the Appendix B stopping rule compares each merge cost against
  `lambda * |A| * |B|`, and that factor grows as groups absorb one another, so
  the effective per-merge threshold is not monotone in the merge order and some
  group counts are unreachable by any single `lambda`. On the paper's first
  application the thresholds run 0.079, 0.686, 0.397, 1.350, 2.765, 3.571, so
  nothing selects five groups. `by = "m"` walks the path one merge at a time.
* `l0_ph()` returns `pair_product` in its `path`, the `|A| * |B|` of each
  merge.

## Simulation

* `ph_design()`, `ph_truth()`, `ph_sample()` and `sim_study()` reproduce the
  calibrated Monte Carlo of Section 4.

## Notes on choices that move results

* The BIC penalty on the two-stage route needs a sample size that the theory
  does not fix. The default is `K`, which reproduces the paper's reported group
  counts; `n_bic` exposes the choice because it changes the selected partition.
* The default base measure is data-scaled (`mu0 = median(tau)`,
  `sigma0_sq = (10 * sd(tau))^2`), matching the paper's applications. This is
  mildly empirical-Bayes; supply fixed values for a prior independent of the
  data.
* `bayes_ph(marginal = "diagonal")` reproduces the faster shortcut used in the
  original replication scripts for the cluster-assignment moves. The default is
  `"exact"`, matching what the paper reports.
