## Lemma 1 says the analysis depends on the data only through R'Sigma^-1 tau
## and R'Sigma^-1 R. These tests hold the package to that claim: the two entry
## points must give not merely similar but identical answers.

test_that("the micro-panel and two-stage routes give the same fit (Lemma 1)", {
  set.seed(11)
  N <- 120
  T <- 6
  cohorts <- c(0, 3, 5)
  cohort <- rep(cohorts, each = N / length(cohorts))

  panel <- expand.grid(id = seq_len(N), t = seq_len(T))
  panel$g <- cohort[panel$id]
  panel$y <- rnorm(N)[panel$id] + 0.1 * panel$t +
    0.5 * (panel$g > 0 & panel$t >= panel$g) + rnorm(nrow(panel))

  from_panel <- ph_data(panel, yname = "y", idname = "id", tname = "t",
                        gname = "g")

  # Feed the fitted estimates and covariance back in as if they had come from
  # an external first stage.
  from_pair <- ph_data(from_panel$tau, Sigma = from_panel$Sigma,
                       cells = from_panel$cells)

  z <- rep(1:2, length.out = from_panel$K)
  a <- ph_fit(from_panel, z)
  b <- ph_fit(from_pair, z)

  expect_equal(a$tau, b$tau)
  expect_equal(a$vcov, b$vcov)
  expect_equal(a$deviance, b$deviance)
})

test_that("scaling the covariance leaves the partition search unchanged", {
  set.seed(12)
  K <- 7
  tau <- c(0.1, 0.11, 0.12, 0.4, 0.41, 0.8, 0.81)
  Sigma <- 0.02^2 * diag(K)
  mk <- function(S) ph_data(tau, S, check_diagonal = FALSE)

  a <- l0_ph(mk(Sigma))
  # The BIC penalty is on the deviance scale, so a common rescaling of the
  # covariance genuinely changes the trade-off; the *path* however is scale
  # free, since every merge cost scales identically.
  b <- l0_ph(mk(4 * Sigma))
  expect_equal(a$partitions, b$partitions)
})

test_that("standard errors supplied as a vector build the diagonal covariance", {
  se <- c(0.01, 0.02, 0.03)
  expect_warning(d <- ph_data(c(1, 2, 3), Sigma = se), "diagonal")
  expect_equal(diag(d$Sigma), se^2)
})

test_that("a diagonal covariance triggers the Remark 1 warning", {
  expect_warning(ph_data(c(1, 2, 3), Sigma = diag(3)), "independent")
  # A correlated covariance should not.
  S <- matrix(0.3, 3, 3) + diag(3) * 0.7
  expect_silent(ph_data(c(1, 2, 3), Sigma = S))
})

test_that("weights are normalised and used by the overall aggregator", {
  d <- ph_data(c(0, 1), Sigma = 0.01 * (0.3 + 0.7 * diag(2)), weights = c(3, 1))
  expect_equal(d$weights, c(0.75, 0.25))
  agg <- aggregate(flex_twfe(d), "overall")
  expect_equal(agg$estimate, 0.25)
})
