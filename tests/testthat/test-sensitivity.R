## Regularisation paths, aggregate and per-cell.

make_d <- function() {
  ph_data(c(0.10, 0.11, 0.42, 0.40, 0.41),
          Sigma = 0.02^2 * (0.3 + 0.7 * diag(5)),
          cells = c("2:2", "2:3", "3:3", "3:4", "4:4"))
}

test_that("lambda_sensitivity spans flexible to pooled", {
  d <- make_d()
  out <- lambda_sensitivity(d, type = "cells", by = "m")

  expect_setequal(unique(out$m), seq_len(d$K))
  expect_equal(nrow(out), d$K^2)

  # At m = K the grouped effects are the flexible estimates.
  at_K <- out[out$m == d$K, ]
  at_K <- at_K[match(d$cells$label, at_K$term), ]
  expect_equal(at_K$estimate, d$tau, tolerance = 1e-10)

  # At m = 1 every cell shares the pooled value.
  at_1 <- out[out$m == 1L, ]
  expect_equal(length(unique(round(at_1$estimate, 12))), 1L)
  expect_equal(at_1$estimate[1], pooled_twfe(d)$tau[1], tolerance = 1e-10)
})

test_that("by = 'm' is gap-free where by = 'lambda' need not be", {
  d <- make_d()
  by_m <- lambda_sensitivity(d, type = "overall", by = "m")
  expect_setequal(unique(by_m$m), seq_len(d$K))

  # The lambda path visits a subset: the stopping rule's effective threshold
  # cost/(|A||B|) is not monotone in the merge order, so some group counts can
  # be unreachable by any single lambda.
  by_lam <- lambda_sensitivity(d, type = "overall", by = "lambda")
  expect_true(all(unique(by_lam$m) %in% seq_len(d$K)))
  expect_true(1L %in% by_lam$m)
  expect_true(d$K %in% by_lam$m)
})

test_that("a supplied lambda grid is honoured and is monotone in pooling", {
  d <- make_d()
  out <- lambda_sensitivity(d, lambda_grid = c(1e-8, 1e8), type = "overall")
  expect_equal(nrow(out), 2L)
  expect_equal(out$m[1], d$K)   # no merges
  expect_equal(out$m[2], 1L)    # everything merged
})

test_that("alpha_sensitivity returns one row per alpha and cell", {
  d <- make_d()
  grid <- c(0.1, 5)
  out <- alpha_sensitivity(d, alpha_grid = grid, type = "cells",
                           iters = 400, burn = 100, seed = 1)

  expect_equal(nrow(out), length(grid) * d$K)
  expect_setequal(unique(out$term), d$cells$label)
  expect_setequal(unique(out$alpha), grid)
  expect_true(all(out$conf.low <= out$estimate + 1e-12))
  expect_true(all(out$estimate <= out$conf.high + 1e-12))
})

test_that("larger alpha spreads the cells apart", {
  d <- make_d()
  out <- alpha_sensitivity(d, alpha_grid = c(0.05, 200), type = "cells",
                           iters = 1500, burn = 300, seed = 2)
  spread <- tapply(out$estimate, out$alpha, function(v) diff(range(v)))
  expect_lt(spread[["0.05"]], spread[["200"]])
})

test_that("the aggregate path still works and carries its anchors", {
  d <- make_d()
  out <- alpha_sensitivity(d, alpha_grid = c(0.5, 2), iters = 400, burn = 100,
                           seed = 3)
  expect_equal(unique(out$term), "overall")
  anchors <- attr(out, "anchors")
  expect_named(anchors, c("pooled", "flexible"))
  expect_equal(anchors[["flexible"]], sum(d$weights * d$tau), tolerance = 1e-12)
})

test_that("sensitivity tables carry the metadata the plotter needs", {
  d <- make_d()
  a <- alpha_sensitivity(d, alpha_grid = c(1, 4), type = "cells",
                         iters = 300, burn = 100, seed = 4)
  l <- lambda_sensitivity(d, type = "cells", by = "m")

  expect_equal(attr(a, "param"), "alpha")
  expect_equal(attr(l, "param"), "m")
  expect_equal(attr(a, "type"), "cells")
  expect_equal(attr(l, "flexible"), stats::setNames(d$tau, d$cells$label))
  expect_equal(attr(l, "pooled"), pooled_twfe(d)$tau[1], tolerance = 1e-12)
})

test_that("plot_sensitivity rejects a table it cannot read", {
  expect_error(plot_sensitivity(data.frame(a = 1)), "sensitivity table")
})
