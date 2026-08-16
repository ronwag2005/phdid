## The sampler is validated against an exact computation wherever one exists.

test_that("the Gibbs sampler agrees with exact enumeration (Appendix E)", {
  skip_on_cran()
  set.seed(21)
  tau <- c(0.10, 0.12, 0.40, 0.42, 0.41)
  Sigma <- 0.02^2 * (0.3 + 0.7 * diag(5))
  d <- ph_data(tau, Sigma)

  ex <- enumerate_partitions(d, alpha = 1)
  fit <- bayes_ph(d, alpha = 1, iters = 20000, burn = 2000, seed = 21,
                  progress = FALSE)

  expect_equal(fit$m_mean, ex$m_mean, tolerance = 0.05)
  expect_lt(max(abs(fit$coclust - ex$coclust)), 0.03)

  agg <- aggregate(fit, "overall")
  expect_equal(agg$estimate, ex$overall$mean, tolerance = 0.002)
})

test_that("enumeration counts exactly the Bell number of partitions", {
  # B_1..B_7
  bells <- c(1, 2, 5, 15, 52, 203, 877)
  for (K in 2:7) {
    d <- ph_data(seq_len(K) / 10, Sigma = 0.01 * diag(K),
                 check_diagonal = FALSE)
    ex <- enumerate_partitions(d, alpha = 1)
    expect_equal(ex$n_partitions, bells[K])
    expect_equal(sum(ex$prob), 1, tolerance = 1e-12)
  }
})

test_that("enumerated partitions are unique and canonical", {
  d <- ph_data(seq_len(5) / 10, Sigma = 0.01 * (0.3 + 0.7 * diag(5)))
  ex <- suppressWarnings(enumerate_partitions(d, alpha = 1))
  Z <- ex$partitions
  expect_equal(nrow(unique(Z)), nrow(Z))
  # Restricted growth: z[1] == 1 and each label is at most one more than the
  # running maximum.
  expect_true(all(Z[, 1] == 1L))
  for (j in 2:ncol(Z)) {
    running_max <- apply(Z[, seq_len(j - 1), drop = FALSE], 1, max)
    expect_true(all(Z[, j] <= running_max + 1L))
  }
})

test_that("the CRP prior integrates to one over all partitions", {
  for (K in c(3, 5, 6)) {
    Z <- phdid:::all_partitions(K)
    for (alpha in c(0.5, 1, 4)) {
      total <- sum(exp(apply(Z, 1, phdid:::log_crp_prior, alpha = alpha)))
      expect_equal(total, 1, tolerance = 1e-10)
    }
  }
})

test_that("expected_groups and alpha_for_groups invert each other", {
  for (K in c(7, 18, 50)) {
    for (target in c(2, K / 2, K - 1)) {
      a <- alpha_for_groups(K, target)
      expect_equal(expected_groups(K, a), target, tolerance = 1e-5)
    }
  }
})

test_that("alpha moves the posterior between pooling and the flexible fit", {
  set.seed(22)
  # Effects only mildly separated, so the prior has room to matter. With
  # overwhelming separation the likelihood dominates and both priors return K
  # groups, which is correct but tests nothing.
  d <- ph_data(c(0.00, 0.02, 0.04, 0.06),
               Sigma = 0.02^2 * (0.3 + 0.7 * diag(4)))
  hi <- bayes_ph(d, alpha = 500, iters = 3000, burn = 500, seed = 1,
                 progress = FALSE)
  lo <- bayes_ph(d, alpha = 0.001, iters = 3000, burn = 500, seed = 1,
                 progress = FALSE)
  expect_gt(hi$m_mean, lo$m_mean)
  expect_gt(hi$m_mean, 3)   # near-flexible
  expect_lt(lo$m_mean, 2)   # near-pooled

  # With overwhelming separation the data win regardless of the prior.
  sharp <- ph_data(c(0, 0.3, 0.6, 0.9), Sigma = 0.01^2 * (0.3 + 0.7 * diag(4)))
  sharp_fit <- bayes_ph(sharp, alpha = 0.001, iters = 2000, burn = 500,
                        seed = 2, progress = FALSE)
  expect_equal(sharp_fit$m_mean, 4, tolerance = 0.05)
})

test_that("the co-clustering matrix is a valid similarity matrix", {
  d <- ph_data(c(0.1, 0.11, 0.4, 0.42), Sigma = 0.02^2 * (0.3 + 0.7 * diag(4)))
  fit <- bayes_ph(d, iters = 800, burn = 200, seed = 3, progress = FALSE)
  P <- coclustering(fit)
  expect_equal(diag(P), rep(1, 4), ignore_attr = TRUE)
  expect_equal(P, t(P))
  expect_true(all(P >= 0 & P <= 1))
})

test_that("aggregation over draws propagates partition uncertainty", {
  d <- ph_data(c(0.1, 0.11, 0.4, 0.42), Sigma = 0.02^2 * (0.3 + 0.7 * diag(4)),
               cells = c("2:2", "2:3", "3:3", "3:4"))
  fit <- bayes_ph(d, iters = 1000, burn = 200, seed = 4, progress = FALSE)

  overall <- aggregate(fit, "overall")
  expect_equal(nrow(overall), 1L)
  # Aggregating each draw then averaging equals averaging then aggregating,
  # for the point estimate.
  expect_equal(overall$estimate, sum(d$weights * fit$tau), tolerance = 1e-12)

  dyn <- aggregate(fit, "dynamic")
  expect_setequal(dyn$term, c("e=0", "e=1"))
})

test_that("point_partition returns a valid partition of the cells", {
  d <- ph_data(c(0.1, 0.11, 0.4, 0.42), Sigma = 0.02^2 * (0.3 + 0.7 * diag(4)))
  fit <- bayes_ph(d, iters = 800, burn = 200, seed = 5, progress = FALSE)
  for (loss in c("VI", "binder")) {
    z <- point_partition(fit, loss = loss)
    expect_length(z, 4L)
    expect_true(all(z >= 1 & z <= 4))
    expect_equal(sort(unique(z)), seq_len(max(z)))
  }
})
