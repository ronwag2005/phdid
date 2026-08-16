## Draw estimates that actually have the covariance we hand to ph_data().
## Declaring one covariance while simulating another would make the tests
## measure the mismatch rather than the estimator.
rmvn <- function(n, mu, Sigma) {
  L <- t(chol(Sigma))
  replicate(n, as.numeric(mu + L %*% stats::rnorm(length(mu))))
}

test_that("the common-effect test rejects real heterogeneity and not noise", {
  # Four cells, two well-separated levels: should reject.
  d1 <- ph_data(c(0.10, 0.11, 0.42, 0.40), Sigma = 0.01^2 * (0.3 + 0.7 * diag(4)))
  h1 <- homogeneity_test(d1)
  expect_lt(h1$chisq$p_value, 0.001)
  expect_gt(h1$dispersion$het_share, 0.9)

  # Cells drawn around a single common value: the test should hold its size.
  set.seed(31)
  K <- 8
  Sigma <- 0.02^2 * (0.3 + 0.7 * diag(K))
  draws <- rmvn(500, rep(0.05, K), Sigma)
  rejections <- apply(draws, 2L, function(tau) {
    homogeneity_test(ph_data(tau, Sigma))$chisq$p_value < 0.05
  })
  # Nominal 5%; with 500 replications the standard error is about 1 point.
  expect_lt(mean(rejections), 0.09)
  expect_gt(mean(rejections), 0.02)
})

test_that("the chi-square statistic is the pooled GLS deviance on K-1 df", {
  set.seed(32)
  K <- 6
  A <- matrix(rnorm(K * K), K)
  Sigma <- crossprod(A) + diag(K)
  d <- ph_data(rnorm(K), Sigma)

  h <- homogeneity_test(d)
  expect_equal(h$chisq$statistic, pooled_twfe(d)$deviance)
  expect_equal(h$chisq$df, K - 1L)
})

test_that("the heterogeneity share is zero when dispersion matches the noise", {
  # Under a literally common effect the excess dispersion should be near zero
  # on average.
  set.seed(33)
  K <- 10
  Sigma <- 0.05^2 * (0.3 + 0.7 * diag(K))
  draws <- rmvn(200, rep(0, K), Sigma)
  shares <- apply(draws, 2L, function(tau) {
    homogeneity_test(ph_data(tau, Sigma))$dispersion$het_share
  })
  expect_lt(mean(shares), 0.35)
  expect_true(all(shares >= 0))
})

test_that("the randomization test is calibrated under the null", {
  # Pre and post summaries drawn from the same distribution: p should be
  # roughly uniform, so rejections at 5% should be rare.
  set.seed(34)
  p_values <- replicate(60, {
    tau <- rnorm(12, 0.02, 0.01)
    pre <- rnorm(12, 0, 0.01)
    d <- suppressWarnings(ph_data(tau, Sigma = 0.01^2 * (0.3 + 0.7 * diag(12))))
    homogeneity_test(d, pre = pre, nsim = 500)$placebo$p_value
  })
  expect_lt(mean(p_values < 0.05), 0.2)
})

test_that("the randomization test detects real heterogeneity", {
  set.seed(35)
  # Post-treatment effects far more dispersed than the placebo.
  tau <- c(rep(-0.10, 6), rep(0.10, 6))
  pre <- rnorm(12, 0, 0.005)
  d <- suppressWarnings(ph_data(tau, Sigma = 0.005^2 * (0.3 + 0.7 * diag(12))))
  h <- homogeneity_test(d, pre = pre, nsim = 2000)
  expect_lt(h$placebo$p_value, 0.05)
  expect_gt(h$placebo$snr, 2)
})

test_that("the placebo gauge reproduces the paper's arithmetic", {
  # s_post below s_pre gives a zero heterogeneity share and SNR below one,
  # the configuration of the paper's second application.
  set.seed(36)
  pre <- rnorm(50); pre <- pre / sd(pre) * 0.0065
  tau <- rnorm(50); tau <- tau / sd(tau) * 0.0050
  d <- suppressWarnings(ph_data(tau, Sigma = 0.0065^2 * (0.3 + 0.7 * diag(50))))
  h <- homogeneity_test(d, pre = pre, nsim = 1000)

  expect_equal(h$placebo$s_pre, 0.0065, tolerance = 1e-9)
  expect_equal(h$placebo$s_post, 0.0050, tolerance = 1e-9)
  expect_equal(h$placebo$snr, 0.0050 / 0.0065, tolerance = 1e-6)
  expect_equal(h$placebo$het_share, 0)
})
