test_that("the default design reproduces the paper's Section 4.1 geometry", {
  des <- ph_design()
  expect_equal(des$K, 18L)
  expect_equal(des$N, 2000L)
  expect_equal(des$T, 10L)
  # Three treated cohorts entering at 3, 5, 7 in a 10-period panel:
  # (10-3+1) + (10-5+1) + (10-7+1) = 8 + 6 + 4 = 18.
  expect_equal(nrow(des$cells), 18L)
  expect_true(all(des$cells$t >= des$cells$g))
})

test_that("ph_truth partitions the cells into near-equal groups at the requested separation", {
  des <- ph_design(N = 400, T = 8, cohorts = c(0, 3, 5, 7))
  for (m in c(1, 2, 3, des$K)) {
    truth <- ph_truth(des, m_star = m, delta = 6)
    expect_equal(max(truth$labels), m)
    expect_equal(length(truth$tau), des$K)
    sizes <- tabulate(truth$labels)
    expect_lte(max(sizes) - min(sizes), 1L)
  }

  # Adjacent group means sit delta flexible standard errors apart.
  truth <- ph_truth(des, m_star = 3, delta = 6)
  means <- sort(unique(truth$tau))
  expect_equal(diff(means), rep(6 * des$sd_flex, 2), tolerance = 1e-10)
})

test_that("the flexible estimator is unbiased for the truth", {
  set.seed(41)
  des <- ph_design(N = 600, T = 6, cohorts = c(0, 3, 5))
  truth <- ph_truth(des, m_star = 2, delta = 6)

  R <- 400
  draws <- replicate(R, ph_sample(des, truth)$tau)
  bias <- rowMeans(draws) - truth$tau
  mc_se <- apply(draws, 1L, stats::sd) / sqrt(R)

  # Judge against Monte Carlo error rather than a fixed tolerance: with K
  # cells a fixed threshold is either vacuous or flaky depending on the design.
  expect_true(all(abs(bias / mc_se) < 4))
})

test_that("the oracle estimator beats flexible on variance under partial homogeneity", {
  set.seed(42)
  des <- ph_design(N = 600, T = 6, cohorts = c(0, 3, 5))
  truth <- ph_truth(des, m_star = 2, delta = 6)

  res <- sim_study(des, m_star = 2, delta = 6, R = 40,
                   methods = c("flexible", "oracle"), progress = FALSE,
                   seed = 42)
  oracle <- res$var_ratio[res$method == "oracle"]
  expect_lt(oracle, 1)
  expect_equal(res$var_ratio[res$method == "flexible"], 1)
  # The oracle knows the truth, so it is unbiased.
  expect_lt(res$abs_bias[res$method == "oracle"], 0.02)
})

test_that("the pooled estimator is precise and biased under heterogeneity", {
  set.seed(43)
  des <- ph_design(N = 600, T = 6, cohorts = c(0, 3, 5))
  res <- sim_study(des, m_star = 3, delta = 6, R = 30,
                   methods = c("pooled", "flexible"), progress = FALSE,
                   seed = 43)
  pooled <- res[res$method == "pooled", ]
  expect_lt(pooled$var_ratio, 0.5)
  expect_gt(pooled$abs_bias, 0.05)
  # Short intervals centred on the wrong estimand: coverage collapses.
  expect_lt(pooled$cover_CATT, 0.5)
})

test_that("flexible intervals are correctly calibrated", {
  set.seed(44)
  des <- ph_design(N = 600, T = 6, cohorts = c(0, 3, 5))
  res <- sim_study(des, m_star = 3, delta = 6, R = 60,
                   methods = "flexible", progress = FALSE, seed = 44)
  expect_equal(res$cover_CATT, 0.95, tolerance = 0.05)
})

test_that("the adjusted Rand index behaves at its extremes", {
  expect_equal(adjusted_rand(c(1, 1, 2, 2), c(1, 1, 2, 2)), 1)
  expect_equal(adjusted_rand(c(1, 1, 2, 2), c(2, 2, 1, 1)), 1)
  expect_equal(adjusted_rand(rep(1, 4), rep(1, 4)), 1)
  expect_lt(adjusted_rand(c(1, 1, 2, 2), c(1, 2, 1, 2)), 0.1)
})
