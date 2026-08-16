# phdid 0.1.0

First release, implementing Arora and Wagle (2026), "Partial Homogeneity in
Staggered Difference-in-Differences".

## Estimation

* `ph_data()` assembles the inputs from a `did::att_gt()` fit, a micro panel,
  or any vector of first-stage effects with their joint covariance. Lemma 1 of
  the paper guarantees the three routes give identical answers, and the test
  suite asserts it.
* `bayes_ph()` fits the Dirichlet Process mixture by collapsed Gibbs sampling
  (Neal 2000, Algorithm 3), with the exact collapsed marginal likelihood in the
  cluster-assignment moves by default.
* `l0_ph()` computes the two-step ℓ₀ estimator of Appendix B, supporting both
  BIC selection along the agglomeration path and the fixed-`lambda` stopping
  rule.
* `ph_fit()`, `flex_twfe()` and `pooled_twfe()` provide the restricted GLS fit
  under a known partition and the two benchmark corners.

## Inference

* `aggregate()` applies overall, event-study, by-cohort and by-calendar
  aggregators to each posterior draw, so reported intervals propagate
  uncertainty about the partition itself.
* `coclustering()` and `point_partition()` summarise the grouping structure,
  the latter by Wade–Ghahramani VI or Binder loss.

## Specification testing

* `homogeneity_test()` reports a model-implied common-effect χ² test on the
  pooled GLS deviance, an excess-dispersion heterogeneity share computed
  against the exact covariance, and — when placebo estimates are supplied — the
  pre/post noise gauge and within-cell randomization test of Section 5.2.3.

## Diagnostics

* `enumerate_partitions()` computes the exact partition posterior for small
  designs, the Appendix E benchmark.
* `ph_rhat()` reports Gelman–Rubin statistics and effective sample sizes across
  dispersed chains.
* `alpha_sensitivity()` traces the posterior across the concentration
  parameter, and `covariance_check()` compares exact against diagonal handling
  of the covariance.

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
