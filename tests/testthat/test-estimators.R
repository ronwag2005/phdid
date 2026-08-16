## These tests check identities the theory forces, not stored numbers. If an
## identity breaks, the estimator is wrong, not merely different.

test_that("the singleton partition returns the first stage unchanged", {
  set.seed(1)
  K <- 6
  tau <- rnorm(K)
  A <- matrix(rnorm(K * K), K)
  Sigma <- crossprod(A) + diag(K)
  d <- ph_data(tau, Sigma)

  fit <- flex_twfe(d)
  # tau = R phi with R = I gives back the GLS solution to an exactly
  # identified system, which is the first-stage estimate itself.
  expect_equal(fit$tau, tau)
  expect_equal(fit$vcov, Sigma, ignore_attr = TRUE)
  expect_equal(fit$deviance, 0)
  expect_equal(fit$m, K)
})

test_that("pooling K cells in a balanced orthogonal design divides the variance by |C_p| (eq. 16)", {
  K <- 8
  s2 <- 0.04
  d <- ph_data(rep(0.1, K), Sigma = s2 * diag(K), check_diagonal = FALSE)

  # One group of 3, one of 5.
  z <- c(1, 1, 1, 2, 2, 2, 2, 2)
  fit <- ph_fit(d, z)
  ratio <- diag(fit$vcov) / diag(d$Sigma)
  sizes <- tabulate(z)[z]

  expect_equal(ratio, 1 / sizes, tolerance = 1e-12)
})

test_that("the grouped fit solves the grouped normal equations", {
  set.seed(2)
  K <- 7
  A <- matrix(rnorm(K * K), K)
  Sigma <- crossprod(A) + diag(K)
  d <- ph_data(rnorm(K), Sigma)
  z <- c(1, 1, 2, 2, 3, 3, 3)

  fit <- ph_fit(d, z)
  R <- matrix(0, K, 3)
  R[cbind(seq_len(K), z)] <- 1
  # R' Omega (tau - R phi) = 0
  resid <- d$tau - fit$tau
  expect_equal(as.numeric(crossprod(R, d$Omega %*% resid)), rep(0, 3),
               tolerance = 1e-10)
})

test_that("in an orthogonal design the greedy group means equal the GLS refit (eq. 49)", {
  set.seed(3)
  K <- 6
  # Orthogonal design: diagonal covariance, so Sigma^-1 is diagonal and the
  # GLS fit reduces to the precision-weighted group mean.
  se <- runif(K, 0.01, 0.05)
  tau <- rnorm(K)
  d <- ph_data(tau, Sigma = diag(se^2), check_diagonal = FALSE)
  z <- c(1, 1, 1, 2, 2, 2)

  fit <- ph_fit(d, z)
  prec <- 1 / se^2
  manual <- tapply(prec * tau, z, sum) / tapply(prec, z, sum)
  expect_equal(fit$tau, as.numeric(manual)[z], tolerance = 1e-12)
})

test_that("the pooled fit is the GLS-weighted average of all cells", {
  set.seed(4)
  K <- 5
  A <- matrix(rnorm(K * K), K)
  Sigma <- crossprod(A) + diag(K)
  d <- ph_data(rnorm(K), Sigma)

  fit <- pooled_twfe(d)
  one <- rep(1, K)
  manual <- sum(d$Omega %*% d$tau) / sum(d$Omega)
  expect_equal(unique(round(fit$tau, 12)), round(manual, 12))
  expect_equal(fit$m, 1L)
})

test_that("l0_ph re-estimates rather than reporting greedy means in a correlated design", {
  set.seed(5)
  K <- 6
  tau <- c(0.10, 0.11, 0.12, 0.50, 0.51, 0.52)
  rho <- 0.6
  Sigma <- 0.02^2 * (rho + (1 - rho) * diag(K))
  d <- ph_data(tau, Sigma)

  fit <- l0_ph(d)
  # Whatever partition is chosen, the reported effects must satisfy the exact
  # grouped normal equations under the full covariance.
  R <- matrix(0, K, fit$m)
  R[cbind(seq_len(K), fit$partition)] <- 1
  resid <- d$tau - fit$tau
  expect_equal(as.numeric(crossprod(R, d$Omega %*% resid)), rep(0, fit$m),
               tolerance = 1e-10)
})

test_that("the l0 path runs from K groups down to 1 and is nested", {
  set.seed(6)
  d <- ph_data(rnorm(7), Sigma = 0.1 * (0.3 + 0.7 * diag(7)))
  fit <- l0_ph(d)

  expect_equal(nrow(fit$path), 7L)
  expect_equal(fit$path$m, 1:7)
  expect_equal(max(fit$partitions[[7]]), 7L)
  expect_equal(max(fit$partitions[[1]]), 1L)

  # Agglomerative merges are nested: cells together at m must stay together
  # at m - 1.
  for (m in 7:2) {
    fine <- fit$partitions[[m]]
    coarse <- fit$partitions[[m - 1]]
    same_fine <- outer(fine, fine, "==")
    same_coarse <- outer(coarse, coarse, "==")
    expect_true(all(same_coarse[same_fine]))
  }
})

test_that("deviance decreases monotonically as groups are added", {
  set.seed(7)
  d <- ph_data(rnorm(8), Sigma = 0.05 * (0.3 + 0.7 * diag(8)))
  fit <- l0_ph(d)
  expect_true(all(diff(fit$path$deviance) <= 1e-9))
  expect_equal(fit$path$deviance[8], 0, tolerance = 1e-8)
})
