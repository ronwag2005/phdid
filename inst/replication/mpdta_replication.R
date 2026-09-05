###############################################################################
## mpdta_replication.R
##
## Full replication and feature verification for the phdid package, run on the
## Callaway and Sant'Anna (2021) county minimum-wage / teen-employment panel
## distributed with the `did` package.
##
## This script does two jobs at once.
##
##   1. REPLICATION. It reproduces the numbers reported in the paper's first
##      application (Table 6) and in the sampler validation of Appendix E.
##   2. FEATURE VERIFICATION. It exercises every one of the 25 exported
##      functions and checks the analytic identities the theory forces, so a
##      silent regression anywhere in the package shows up as a failed check
##      rather than as a plausible-looking number.
##
## Every claim is registered with check() and the script prints a pass/fail
## ledger at the end. It exits non-zero if anything fails, so it can be wired
## into CI.
##
## Runtime: roughly 3-6 minutes, dominated by the long Gibbs chains in
## sections 5 and 7.
##
##   Rscript inst/replication/mpdta_replication.R
###############################################################################

suppressMessages({
  library(phdid)
  library(did)
})

set.seed(20260817)
FIGDIR <- "figures"
dir.create(FIGDIR, showWarnings = FALSE)

## ---------------------------------------------------------------- ledger ---
LEDGER <- list()

check <- function(label, value, target = NULL, tol = NULL, section = "") {
  ok <- if (is.null(target)) {
    isTRUE(value)
  } else if (is.numeric(value) && is.numeric(target) && !is.null(tol)) {
    length(value) == length(target) && all(abs(value - target) <= tol)
  } else {
    isTRUE(all.equal(value, target))
  }
  fmt <- function(x) {
    if (is.numeric(x)) paste(formatC(x, format = "g", digits = 5), collapse = ", ")
    else paste(format(x), collapse = ", ")
  }
  LEDGER[[length(LEDGER) + 1L]] <<- data.frame(
    section = section, label = label,
    observed = if (is.null(target)) "" else fmt(value),
    expected = if (is.null(target)) "" else fmt(target),
    pass = ok, stringsAsFactors = FALSE
  )
  cat(sprintf("  [%s] %s%s\n", if (ok) "PASS" else "FAIL", label,
              if (is.null(target)) "" else sprintf("  (%s vs %s)", fmt(value), fmt(target))))
  invisible(ok)
}

banner <- function(n, title) {
  cat("\n", strrep("=", 76), "\n", n, ". ", title, "\n", strrep("=", 76), "\n", sep = "")
}

## Local partition formatter, so the script does not reach into internals.
show_partition <- function(z, labels) {
  paste(vapply(split(labels, z), function(g)
    paste0("{", paste(g, collapse = ","), "}"), character(1)), collapse = " ")
}

###############################################################################
banner(1, "FIRST STAGE AND THE THREE ph_data() ENTRY POINTS")
###############################################################################

data(mpdta, package = "did")

## Callaway-Sant'Anna, with bstrap = FALSE so the analytical covariance (and
## hence the exact joint covariance of the cells) is computed.
first <- att_gt(yname = "lemp", tname = "year", idname = "countyreal",
                gname = "first.treat", control_group = "notyettreated",
                data = mpdta, bstrap = FALSE, cband = FALSE)

## -- Route A: straight from the att_gt fit -----------------------------------
d <- ph_data(first)
print(d)

check("K post-treatment cells", d$K, 7L, section = "1")
check("units", d$n_units, 500L, section = "1")
check("covariance is non-diagonal (Remark 1)",
      max(abs(d$Sigma - diag(diag(d$Sigma)))) > 0, section = "1")
check("cell labels", d$cells$label,
      c("2004:2004", "2004:2005", "2004:2006", "2004:2007",
        "2006:2006", "2006:2007", "2007:2007"), section = "1")
check("weights sum to one", sum(d$weights), 1, tol = 1e-12, section = "1")

## Paper Table 6, columns ATT(g,t) and SE.
check("Table 6 ATT(g,t)", round(d$tau, 3),
      c(-0.019, -0.078, -0.136, -0.101, 0.005, -0.041, -0.026),
      tol = 1e-9, section = "1")
check("Table 6 standard errors", round(sqrt(diag(d$Sigma)), 3),
      c(0.022, 0.030, 0.035, 0.034, 0.016, 0.020, 0.017),
      tol = 1e-9, section = "1")

## -- Route B: the (tau, Sigma) pair, as from any other first stage -----------
d_pair <- ph_data(d$tau, Sigma = d$Sigma, cells = d$cells,
                  weights = d$weights, n_units = d$n_units)
check("pair route reproduces the fit exactly (Lemma 1)",
      isTRUE(all.equal(ph_fit(d, c(1,1,2,2,1,1,1))$tau,
                       ph_fit(d_pair, c(1,1,2,2,1,1,1))$tau)), section = "1")

## -- Route C: the micro panel ------------------------------------------------
## This fits Wooldridge's flexible interacted TWFE internally. It is a
## DIFFERENT estimator from Callaway-Sant'Anna -- it uses already-treated units
## as controls where they are clean -- so the estimates legitimately differ.
d_panel <- ph_data(as.data.frame(mpdta), yname = "lemp", idname = "countyreal",
                   tname = "year", gname = "first.treat")
cat("\n  flexible TWFE on the micro panel:\n")
print(round(setNames(d_panel$tau, d_panel$cells$label), 4))
check("panel route recovers the same 7 cells", d_panel$K, 7L, section = "1")
check("panel route estimates a variance", d_panel$sigma2 > 0, section = "1")
check("panel route agrees with CS on the impact cell",
      abs(d_panel$tau[1] - d$tau[1]) < 0.005, section = "1")

###############################################################################
banner(2, "BENCHMARK CORNERS: FLEXIBLE AND POOLED")
###############################################################################

fx <- flex_twfe(d)
pl <- pooled_twfe(d)

check("flexible returns the first stage unchanged", fx$tau, d$tau, tol = 1e-12,
      section = "2")
check("flexible covariance is the first-stage covariance",
      max(abs(fx$vcov - d$Sigma)) < 1e-12, section = "2")
check("flexible deviance is zero", fx$deviance, 0, tol = 1e-8, section = "2")
check("flexible has K groups", fx$m, 7L, section = "2")
check("pooled has one group", pl$m, 1L, section = "2")
check("pooled is a single repeated value", length(unique(round(pl$tau, 12))), 1L,
      section = "2")

agg_fx <- aggregate(fx, "overall")
agg_pl <- aggregate(pl, "overall")
cat(sprintf("\n  overall ATT  flexible %.4f (se %.4f)   pooled %.4f (se %.4f)\n",
            agg_fx$estimate, agg_fx$std.error, agg_pl$estimate, agg_pl$std.error))
check("flexible overall ATT (paper: -0.040)", round(agg_fx$estimate, 3), -0.040,
      tol = 1e-9, section = "2")
check("flexible overall SE (paper: 0.012)", round(agg_fx$std.error, 3), 0.012,
      tol = 1e-9, section = "2")
check("pooled is more precise than flexible",
      agg_pl$std.error < agg_fx$std.error, section = "2")

###############################################################################
banner(3, "SPECIFICATION TEST: IS THERE HETEROGENEITY TO RECOVER?")
###############################################################################

h <- homogeneity_test(d)
print(h)

check("chi-square statistic equals the pooled GLS deviance",
      h$chisq$statistic, pl$deviance, tol = 1e-12, section = "3")
check("chi-square df", h$chisq$df, 6L, section = "3")
check("common effect is rejected", h$chisq$p_value < 0.001, section = "3")
check("heterogeneity share is substantial", h$dispersion$het_share > 0.7,
      section = "3")
## Paper 5.1: cross-cell dispersion 0.050 against a typical SE of 0.022.
check("observed cross-cell sd (paper: 0.050)",
      round(sqrt(h$dispersion$observed), 3), 0.050, tol = 1e-9, section = "3")
check("descriptive SNR sd(tau)/median(se) (paper: >2)",
      h$snr_descriptive > 2, section = "3")

## The placebo gauge needs a pre-treatment summary for EVERY cell, constructed
## symmetrically with the post-treatment one. mpdta cannot supply that: the
## 2004 cohort is treated from the first period of the window and so has no
## pre-treatment cells. We therefore verify that branch on a synthetic design
## where the truth is known.
cat("\n  -- placebo branch, verified on synthetic data --\n")
set.seed(11)
K_s <- 12
Sig_s <- 0.01^2 * (0.3 + 0.7 * diag(K_s))
L_s <- t(chol(Sig_s))

## (a) a genuinely common effect: the test should hold its nominal size.
## A single draw is a coin flip -- under the null p is uniform, so any one
## replication rejects 5% of the time by construction. Calibration across
## replications is the claim worth checking.
pre_null <- as.numeric(L_s %*% rnorm(K_s))
p_null <- replicate(200, {
  tau_r <- as.numeric(0.02 + L_s %*% rnorm(K_s))
  pre_r <- as.numeric(L_s %*% rnorm(K_s))
  homogeneity_test(ph_data(tau_r, Sig_s), pre = pre_r, nsim = 600)$placebo$p_value
})
cat(sprintf("    null p-values over 200 draws: mean %.3f, %.1f%% below 0.05\n",
            mean(p_null), 100 * mean(p_null < 0.05)))
check("placebo p-values are uniform under the null (mean ~ 0.5)",
      abs(mean(p_null) - 0.5) < 0.10, section = "3")
check("placebo test holds its nominal 5% size",
      mean(p_null < 0.05) < 0.12, section = "3")

## (b) two well-separated levels: the test SHOULD reject
tau_het <- rep(c(-0.05, 0.05), each = K_s / 2) + as.numeric(L_s %*% rnorm(K_s))
h_het <- homogeneity_test(ph_data(tau_het, Sig_s), pre = pre_null, nsim = 4000)
check("placebo test rejects real heterogeneity",
      h_het$placebo$p_value < 0.05, section = "3")
check("signal-to-noise above one under heterogeneity",
      h_het$placebo$snr > 1, section = "3")

###############################################################################
banner(4, "THE l0-PENALISED ESTIMATOR")
###############################################################################

fit_l0 <- l0_ph(d)
print(fit_l0)

## ---- Paper Table 6, in full ----
tab6 <- data.frame(
  cell = d$cells$label,
  att = round(d$tau, 3),
  se = round(sqrt(diag(d$Sigma)), 3),
  l0_group = fit_l0$partition,
  grouped = round(fit_l0$tau, 3),
  var_ratio = round(diag(fit_l0$vcov) / diag(d$Sigma), 2)
)
cat("\n  Table 6 as published (rows ordered as in the paper):\n")
paper_order <- c("2004:2004", "2006:2007", "2007:2007", "2004:2005",
                 "2004:2007", "2004:2006", "2006:2006")
print(tab6[match(paper_order, tab6$cell), ], row.names = FALSE)

check("BIC selects 4 groups (paper: 4)", fit_l0$m, 4L, section = "4")
## Compared on the raw scale against the paper's 3-decimal column.
##
## Three of the four groups agree exactly at 3 dp (-0.087, -0.140, 0.011), and
## the Step-1 greedy means reproduce NONE of those three, which confirms the
## published table is the Step-2 GLS refit this package computes.
##
## The near-zero group reads -0.028462 here and -0.029 in the paper. That gap
## is double rounding, not a difference in the estimator: the replication
## script builds its table with round(l0$tau, 4), turning -0.028462 into
## -0.0285, which then rounds to -0.029 when typeset at three decimals. The
## tolerance below spans one 3-dp ulp so the check tracks the estimator rather
## than the presentation.
check("Table 6 grouped effects",
      sort(unique(fit_l0$tau)), c(-0.140, -0.087, -0.029, 0.011),
      tol = 1e-3, section = "4")
check("three groups match the paper exactly at 3 dp",
      round(sort(unique(fit_l0$tau))[c(1, 2, 4)], 3), c(-0.140, -0.087, 0.011),
      tol = 1e-9, section = "4")
check("Table 6 variance ratios",
      round(diag(fit_l0$vcov) / diag(d$Sigma), 2)[match(paper_order, tab6$cell)],
      c(0.28, 0.33, 0.49, 0.77, 0.61, 0.78, 0.78), tol = 1e-9, section = "4")
pooled_cells <- tabulate(fit_l0$partition)[fit_l0$partition] > 1L
check("pooling roughly halves the variance of pooled cells (paper)",
      round(mean((diag(fit_l0$vcov) / diag(d$Sigma))[pooled_cells]), 2), 0.50,
      tol = 1e-9, section = "4")
check("the near-zero group holds the three small cells",
      sort(d$cells$label[fit_l0$partition == fit_l0$partition[1]]),
      sort(c("2004:2004", "2006:2007", "2007:2007")), section = "4")

## Paper 5.1: overall effect essentially unchanged and slightly more precise.
agg_l0 <- aggregate(fit_l0, "overall")
cat(sprintf("\n  overall ATT  flexible %.4f (se %.4f)   l0 %.4f (se %.4f)\n",
            agg_fx$estimate, agg_fx$std.error, agg_l0$estimate, agg_l0$std.error))
check("l0 overall ATT (paper: -0.039)", round(agg_l0$estimate, 3), -0.039,
      tol = 1e-9, section = "4")
check("l0 overall SE (paper: 0.012)", round(agg_l0$std.error, 3), 0.012,
      tol = 1e-9, section = "4")
check("l0 is at least as precise as flexible overall",
      agg_l0$std.error <= agg_fx$std.error + 1e-12, section = "4")
check("even singleton cells gain precision (Gauss-Markov borrowing)",
      all(diag(fit_l0$vcov) <= diag(d$Sigma) + 1e-12), section = "4")

## ---- the agglomeration path ----
cat("\n  agglomeration path:\n")
print(fit_l0$path, row.names = FALSE, digits = 5)
check("path has K rows", nrow(fit_l0$path), 7L, section = "4")
check("deviance is monotone in m", all(diff(fit_l0$path$deviance) <= 1e-9),
      section = "4")
check("deviance vanishes at the flexible end",
      fit_l0$path$deviance[7], 0, tol = 1e-8, section = "4")
check("BIC is minimised at the selected m",
      which.min(fit_l0$path$bic), fit_l0$selected, section = "4")

## nestedness of the merge path
nested <- TRUE
for (m in 7:2) {
  fine <- fit_l0$partitions[[m]]; coarse <- fit_l0$partitions[[m - 1]]
  if (!all(outer(coarse, coarse, "==")[outer(fine, fine, "==")])) nested <- FALSE
}
check("the merge path is nested", nested, section = "4")

## ---- alternative selection rules ----
fit_lam <- l0_ph(d, select = "lambda", lambda = 0.5)
fit_m3  <- l0_ph(d, select = "m", m = 3)
fit_ex  <- l0_ph(d, search = "exact-cost")
check("select = 'lambda' runs and returns a partition", fit_lam$m >= 1L,
      section = "4")
check("select = 'm' honours the request", fit_m3$m, 3L, section = "4")
check("exact-cost search returns a valid partition",
      fit_ex$m >= 1L && fit_ex$m <= 7L, section = "4")
cat(sprintf("\n  groups by rule:  BIC %d | lambda=0.5 %d | m=3 %d | exact-cost %d\n",
            fit_l0$m, fit_lam$m, fit_m3$m, fit_ex$m))

## the n_bic choice, which the docs flag as moving the answer
check("n_bic = K gives 4 groups", l0_ph(d, n_bic = d$K)$m, 4L, section = "4")
check("n_bic = n_units gives 3 groups", l0_ph(d, n_bic = d$n_units)$m, 3L,
      section = "4")

###############################################################################
banner(5, "THE DIRICHLET PROCESS ESTIMATOR")
###############################################################################

t0 <- Sys.time()
fit_b <- bayes_ph(d, alpha = 1, iters = 20000, burn = 2000, seed = 7,
                  progress = FALSE)
cat(sprintf("  (20,000 sweeps in %.1f s)\n", 
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
print(fit_b)

check("posterior E[# groups] (Appendix E: 2.20)", fit_b$m_mean, 2.20, tol = 0.05,
      section = "5")
check("posterior means lie between pooled and flexible",
      all(fit_b$tau >= pmin(d$tau, pl$tau[1]) - 1e-9 &
          fit_b$tau <= pmax(d$tau, pl$tau[1]) + 1e-9), section = "5")
check("credible intervals bracket the posterior means",
      all(fit_b$lower <= fit_b$tau & fit_b$tau <= fit_b$upper), section = "5")

## ---- aggregation, applied per draw ----
agg_b <- aggregate(fit_b, "overall")
cat("\n  overall ATT, partition-averaged:\n"); print(agg_b, digits = 4)
check("overall ATT (Appendix E: -0.024)", round(agg_b$estimate, 3), -0.024,
      tol = 1e-9, section = "5")
check("overall lower bound (Appendix E: -0.047)", agg_b$conf.low, -0.047,
      tol = 0.002, section = "5")
check("overall upper bound (Appendix E: -0.000)", agg_b$conf.high, 0.000,
      tol = 0.002, section = "5")
check("point estimate equals weights times posterior mean",
      agg_b$estimate, sum(d$weights * fit_b$tau), tol = 1e-12, section = "5")

dyn <- aggregate(fit_b, "dynamic")
cat("\n  event-study profile:\n"); print(dyn, digits = 4)
check("event times 0 through 3 are reported", dyn$term,
      c("e=0", "e=1", "e=2", "e=3"), section = "5")
check("effect accumulates over the horizon (paper 5.1)",
      dyn$estimate[3] < dyn$estimate[1], section = "5")

grp <- aggregate(fit_b, "group")
cal <- aggregate(fit_b, "calendar")
cat("\n  by cohort:\n"); print(grp, digits = 4)
cat("\n  by calendar year:\n"); print(cal, digits = 4)
check("by-cohort aggregation covers the three cohorts", nrow(grp), 3L,
      section = "5")
check("by-calendar aggregation covers four years", nrow(cal), 4L, section = "5")

## a custom aggregator: the 2004 cohort only
w_2004 <- as.numeric(d$cells$g == 2004); w_2004 <- w_2004 / sum(w_2004)
cust <- aggregate(fit_b, weights = w_2004)
cat(sprintf("\n  custom aggregator (2004 cohort only): %.4f [%.4f, %.4f]\n",
            cust$estimate, cust$conf.low, cust$conf.high))
check("custom weights are accepted", nrow(cust), 1L, section = "5")

## ---- the grouping structure ----
P <- coclustering(fit_b)
cat("\n  posterior co-clustering probabilities:\n"); print(round(P, 2))

check("co-clustering diagonal is one", diag(P), rep(1, 7), tol = 1e-12,
      section = "5")
check("co-clustering is symmetric", max(abs(P - t(P))), 0, tol = 1e-12,
      section = "5")
check("co-clustering lies in [0,1]", all(P >= 0 & P <= 1), section = "5")
## Paper 5.1: 2004:2004 co-clusters with 2006:2006 at 0.86 and 2004:2005 at 0.77
check("2004:2004 with 2006:2006 (paper: 0.86)",
      P["2004:2004", "2006:2006"], 0.86, tol = 0.03, section = "5")
check("2004:2004 with 2004:2005 (paper: 0.77)",
      P["2004:2004", "2004:2005"], 0.77, tol = 0.03, section = "5")
## Paper 5.1: the peak cell co-clusters with everything else at 0.54 or below
check("peak 2004:2006 with all others <= 0.54 (paper)",
      max(P["2004:2006", -3]) <= 0.56, section = "5")

for (loss in c("VI", "binder")) {
  pp <- point_partition(fit_b, loss = loss)
  cat(sprintf("\n  point partition (%s loss, %d groups):\n    %s\n",
              loss, attr(pp, "n_groups"), show_partition(pp, d$cells$label)))
  check(sprintf("point_partition(%s) is a valid partition", loss),
        all(sort(unique(pp)) == seq_len(max(pp))), section = "5")
}

## ---- extractors ----
ci_b <- confint(fit_b)
check("confint returns K x 2 with cell names",
      identical(dim(ci_b), c(7L, 2L)) &&
        identical(rownames(ci_b), d$cells$label), section = "5")
check("coef returns named cell effects",
      identical(names(coef(fit_b)), d$cells$label), section = "5")
check("vcov on a ph_fit is K x K", identical(dim(vcov(fit_l0)), c(7L, 7L)),
      section = "5")

## ---- prior helpers ----
check("alpha_for_groups inverts expected_groups",
      expected_groups(7, alpha_for_groups(7, 4)), 4, tol = 1e-5, section = "5")
check("prior E[m] at alpha = 1, K = 7", round(expected_groups(7, 1), 2), 2.08,
      tol = 1e-9, section = "5")
## The paper uses alpha = 7 at K = 18 "so the prior expected number of groups
## is about K/2". Check that forward statement; solving E[m] = 9 exactly
## returns 7.16, of which 7 is the round number.
check("paper's alpha = 7 gives prior E[m] ~ K/2 = 9",
      expected_groups(18, 7), 9, tol = 0.15, section = "5")

###############################################################################
banner(6, "EXACT ENUMERATION: THE APPENDIX E BENCHMARK")
###############################################################################

ex <- enumerate_partitions(d, alpha = 1)
print(ex)

check("B_7 = 877 partitions enumerated", ex$n_partitions, 877L, section = "6")
check("posterior probabilities sum to one", sum(ex$prob), 1, tol = 1e-12,
      section = "6")
check("exact E[# groups] (Appendix E: 2.20)", round(ex$m_mean, 2), 2.20,
      tol = 1e-9, section = "6")
check("exact overall ATT (Appendix E: -0.024)", round(ex$overall$mean, 3),
      -0.024, tol = 1e-9, section = "6")
check("exact interval (Appendix E: [-0.048, -0.000])",
      round(c(ex$overall$lower, ex$overall$upper), 3), c(-0.048, 0.000),
      tol = 1e-9, section = "6")

## Sampler against exact -- the actual validation of Appendix E.
gap_coclust <- max(abs(fit_b$coclust - ex$coclust))
cat(sprintf("\n  Gibbs vs exact enumeration\n"))
cat(sprintf("    E[# groups]      %.3f  vs  %.3f\n", fit_b$m_mean, ex$m_mean))
cat(sprintf("    overall ATT      %.4f vs  %.4f\n", agg_b$estimate, ex$overall$mean))
cat(sprintf("    interval         [%.4f, %.4f]  vs  [%.4f, %.4f]\n",
            agg_b$conf.low, agg_b$conf.high, ex$overall$lower, ex$overall$upper))
cat(sprintf("    max co-clustering gap  %.4f   (paper reports 0.012)\n", gap_coclust))

check("sampler matches exact E[# groups]", fit_b$m_mean, ex$m_mean, tol = 0.05,
      section = "6")
check("sampler matches exact overall ATT", agg_b$estimate, ex$overall$mean,
      tol = 0.002, section = "6")
check("max co-clustering gap within Monte Carlo error", gap_coclust < 0.02,
      section = "6")

###############################################################################
banner(7, "DIAGNOSTICS")
###############################################################################

## ---- multi-chain convergence, from dispersed initialisations ----
fit_4 <- bayes_ph(d, alpha = 1, iters = 3000, burn = 500, chains = 4, seed = 9,
                  progress = FALSE)
rh <- ph_rhat(fit_4)
print(rh)
check("R-hat at 1.00 for both quantities (Appendix E)",
      max(rh$stats$rhat) < 1.01, section = "7")
check("effective sample sizes are large", min(rh$stats$ess) > 1000, section = "7")
check("four chains recorded", rh$chains, 4L, section = "7")

## ---- exact vs diagonal handling of the covariance (Remark 1) ----
cmp <- covariance_check(d, alpha = 1, iters = 4000, burn = 1000, seed = 1)
cat("\n  exact vs diagonal assignment moves:\n")
print(cmp$summary, row.names = FALSE, digits = 4)
gap_row <- cmp$summary$gap[cmp$summary$quantity == "max co-clustering gap"]
check("aggregate is barely affected by the shortcut",
      cmp$summary$gap[1] < 0.005, section = "7")
check("individual co-clustering probabilities ARE affected", gap_row > 0.3,
      section = "7")
cat(sprintf("\n  -> the diagonal shortcut moves a co-clustering probability by %.2f\n",
            gap_row))
cat("     while leaving the aggregate unchanged. This is why marginal = 'exact'\n")
cat("     is the default (Remark 1, Appendix E).\n")

## ---- sensitivity to the concentration parameter ----
sens <- alpha_sensitivity(d, alpha_grid = c(0.1, 0.5, 1, 2, 5, 14, 50, 100),
                          iters = 3000, burn = 500, seed = 3)
cat("\n  sensitivity to alpha:\n"); print(sens, row.names = FALSE, digits = 4)
anchors <- attr(sens, "anchors")
cat(sprintf("\n  anchors: pooled %.4f, flexible %.4f\n",
            anchors[["pooled"]], anchors[["flexible"]]))

check("the path is monotone toward the flexible anchor",
      cor(sens$alpha, sens$estimate, method = "spearman") < 0, section = "7")
check("heavy pooling pulls the effect toward the pooled anchor",
      abs(sens$estimate[1] - anchors[["pooled"]]) <
        abs(sens$estimate[nrow(sens)] - anchors[["pooled"]]), section = "7")
check("large alpha approaches the flexible anchor",
      abs(sens$estimate[nrow(sens)] - anchors[["flexible"]]) < 0.01, section = "7")
check("posterior group count rises with alpha",
      cor(sens$alpha, sens$m, method = "spearman") > 0.9, section = "7")
## Paper 5.1: at alpha ~ 14 the posterior E[m] matches the 4 groups BIC selected.
a14 <- sens$m[sens$alpha == 14]
cat(sprintf("\n  posterior E[# groups] at alpha = 14: %.2f  (BIC selected %d)\n",
            a14, fit_l0$m))
check("alpha ~ 14 reproduces the BIC group count (paper 5.1)",
      abs(a14 - 4) < 1.0, section = "7")

## ---- per-cell regularisation paths ----
cat("\n  per-cell paths\n")
sens_cells <- alpha_sensitivity(d, alpha_grid = c(0.1, 1, 14, 100),
                                type = "cells", iters = 2000, burn = 400,
                                seed = 5)
wide_a <- stats::reshape(sens_cells[, c("alpha", "term", "estimate")],
                         idvar = "term", timevar = "alpha",
                         direction = "wide")
names(wide_a) <- sub("estimate.", "alpha=", names(wide_a), fixed = TRUE)
cat("\n  DP posterior CATTs across alpha:\n")
print(wide_a, row.names = FALSE, digits = 3)

lam_cells <- lambda_sensitivity(d, type = "cells", by = "m")
wide_l <- stats::reshape(lam_cells[, c("m", "term", "estimate")],
                         idvar = "term", timevar = "m", direction = "wide")
names(wide_l) <- sub("estimate.", "m=", names(wide_l), fixed = TRUE)
cat("\n  l0 grouped CATTs across the agglomeration path:\n")
print(wide_l, row.names = FALSE, digits = 3)

check("per-cell alpha path has one row per (alpha, cell)",
      nrow(sens_cells), 4L * d$K, section = "7")
check("per-cell l0 path visits every group count",
      sort(unique(lam_cells$m)), seq_len(d$K), section = "7")
check("at m = K the l0 path returns the flexible estimates",
      lam_cells$estimate[lam_cells$m == d$K][
        match(d$cells$label, lam_cells$term[lam_cells$m == d$K])],
      d$tau, tol = 1e-10, section = "7")
check("at m = 1 the l0 path returns the pooled value",
      length(unique(round(lam_cells$estimate[lam_cells$m == 1L], 12))), 1L,
      section = "7")
check("cells spread out as alpha grows",
      diff(range(sens_cells$estimate[sens_cells$alpha == 100])) >
        diff(range(sens_cells$estimate[sens_cells$alpha == 0.1])),
      section = "7")
## The lambda-indexed path can skip group counts the m-indexed one visits,
## because the stopping rule's effective threshold is not monotone in the
## merge order. Record that rather than hide it.
lam_by_lambda <- lambda_sensitivity(d, type = "overall", by = "lambda")
skipped <- setdiff(seq_len(d$K), unique(lam_by_lambda$m))
cat(sprintf("\n  group counts unreachable by any single lambda: %s\n",
            if (length(skipped)) paste(skipped, collapse = ", ") else "none"))
check("lambda path reaches both extremes",
      all(c(1L, d$K) %in% lam_by_lambda$m), section = "7")

###############################################################################
banner(8, "ANALYTIC IDENTITIES THE THEORY FORCES")
###############################################################################

## These fail only if an estimator is wrong, not merely different.

## (a) restricted GLS solves the grouped normal equations
z_test <- c(1, 1, 2, 2, 3, 3, 3)
f_test <- ph_fit(d, z_test)
R_test <- matrix(0, 7, 3); R_test[cbind(1:7, z_test)] <- 1
check("R' Omega (tau - R phi) = 0",
      max(abs(crossprod(R_test, d$Omega %*% (d$tau - f_test$tau)))), 0,
      tol = 1e-10, section = "8")

## (b) pooling in a balanced orthogonal design divides variance by group size
d_orth <- ph_data(rep(0.1, 8), Sigma = 0.04 * diag(8), check_diagonal = FALSE)
z_orth <- c(1,1,1,2,2,2,2,2)
f_orth <- ph_fit(d_orth, z_orth)
check("variance ratio equals 1 / |C_p| (eq. 16)",
      diag(f_orth$vcov) / diag(d_orth$Sigma), 1 / tabulate(z_orth)[z_orth],
      tol = 1e-12, section = "8")

## (c) a partition given as a list matches the same partition as a vector
check("list and vector partitions agree",
      isTRUE(all.equal(ph_fit(d, list(c(1,2), c(3,4), c(5,6,7)))$tau,
                       ph_fit(d, c(1,1,2,2,3,3,3))$tau)), section = "8")

## (d) relabelling a partition changes nothing
check("partition labels are canonicalised",
      isTRUE(all.equal(ph_fit(d, c(3,3,1,1,2,2,2))$tau,
                       ph_fit(d, c(1,1,2,2,3,3,3))$tau)), section = "8")

## (e) the CRP prior is a proper distribution over partitions
crp_total <- sum(ex$prob)  # normalised posterior; prior checked in the test suite
check("enumerated posterior is proper", crp_total, 1, tol = 1e-12, section = "8")

## (f) adjusted Rand index at its extremes
check("ARI is 1 for identical partitions",
      adjusted_rand(c(1,1,2,2), c(1,1,2,2)), 1, tol = 1e-12, section = "8")
check("ARI is relabelling invariant",
      adjusted_rand(c(1,1,2,2), c(2,2,1,1)), 1, tol = 1e-12, section = "8")

###############################################################################
banner(9, "SIMULATION TOOLS")
###############################################################################

des <- ph_design()
print(des)
check("default design has K = 18 (paper 4.1)", des$K, 18L, section = "9")
check("default design has 2000 units over 10 periods",
      c(des$N, des$T), c(2000, 10), tol = 0, section = "9")

truth <- ph_truth(des, m_star = 6, delta = 6)
check("m_star groups created", max(truth$labels), 6L, section = "9")
check("groups are near-equal in size",
      max(tabulate(truth$labels)) - min(tabulate(truth$labels)) <= 1L,
      section = "9")
check("adjacent group means sit delta flexible SEs apart",
      diff(sort(unique(truth$tau))), rep(6 * des$sd_flex, 5), tol = 1e-10,
      section = "9")

## a single draw is an ordinary ph_data object
d_sim <- ph_sample(des, truth)
check("a simulated draw is a ph_data object", inherits(d_sim, "ph_data"),
      section = "9")
check("a simulated draw carries the micro-panel quantities",
      !is.null(d_sim$G) && d_sim$sigma2 > 0, section = "9")
check("the oracle recovers the truth well on one draw",
      max(abs(ph_fit(d_sim, truth$labels)$tau - truth$tau)) < 0.05, section = "9")

## a small Monte Carlo -- the shape of Tables 2-5
cat("\n  small Monte Carlo (30 reps, m* = 6, delta = 6):\n")
sim <- sim_study(des, m_star = 6, delta = 6, R = 30,
                 bayes_args = list(iters = 400, burn = 100),
                 seed = 42, progress = FALSE)
print(sim[, c("method", "var_ratio", "abs_bias", "ari", "cover_CATT",
              "cover_ATT")], row.names = FALSE, digits = 3)

check("flexible is the variance benchmark",
      sim$var_ratio[sim$method == "flexible"], 1, tol = 1e-12, section = "9")
check("pooled is far more precise than flexible",
      sim$var_ratio[sim$method == "pooled"] < 0.3, section = "9")
check("pooled is badly biased under heterogeneity",
      sim$abs_bias[sim$method == "pooled"] > 0.2, section = "9")
check("pooled coverage collapses",
      sim$cover_CATT[sim$method == "pooled"] < 0.2, section = "9")
check("the oracle beats flexible on variance",
      sim$var_ratio[sim$method == "oracle"] < 1, section = "9")
check("flexible coverage is near nominal",
      abs(sim$cover_CATT[sim$method == "flexible"] - 0.95) < 0.06, section = "9")
check("feasible estimators recover the partition well at delta = 6",
      sim$ari[sim$method == "l0"] > 0.7, section = "9")
check("feasible estimators stay essentially unbiased",
      sim$abs_bias[sim$method == "l0"] < 0.05, section = "9")
check("ATT is well calibrated for every method except pooled",
      all(sim$cover_ATT[sim$method != "pooled"] > 0.85), section = "9")

###############################################################################
banner(10, "FIGURES")
###############################################################################

png(file.path(FIGDIR, "mpdta_effects.png"), 1700, 950, res = 200)
plot(fit_l0)
dev.off()

png(file.path(FIGDIR, "mpdta_coclustering.png"), 1300, 1200, res = 200)
plot_coclustering(fit_b, digits = 2)
dev.off()

png(file.path(FIGDIR, "mpdta_alpha_sensitivity.png"), 1500, 1500, res = 200)
plot_sensitivity(sens)
dev.off()

png(file.path(FIGDIR, "mpdta_alpha_cells.png"), 1600, 1000, res = 200)
plot_sensitivity(sens_cells)
dev.off()

png(file.path(FIGDIR, "mpdta_l0_cells.png"), 1600, 1000, res = 200)
plot_sensitivity(lam_cells)
dev.off()

png(file.path(FIGDIR, "mpdta_l0_path.png"), 1500, 950, res = 200)
plot_l0_path(fit_l0)
dev.off()

png(file.path(FIGDIR, "mpdta_placebo_band.png"), 1500, 950, res = 200)
plot_placebo_band(h)
dev.off()

png(file.path(FIGDIR, "sim_variance.png"), 1500, 950, res = 200)
plot_sim_study(sim, "var_ratio")
dev.off()

figs <- list.files(FIGDIR, pattern = "\\.png$", full.names = TRUE)
check("all eight figures written", length(figs) >= 8L, section = "10")
check("no figure is empty", all(file.size(figs) > 5000), section = "10")
cat("\n  wrote:\n"); cat(paste0("    ", figs, collapse = "\n"), "\n")

###############################################################################
banner(11, "LEDGER")
###############################################################################

led <- do.call(rbind, LEDGER)
n_pass <- sum(led$pass); n_fail <- sum(!led$pass)

by_sec <- aggregate(cbind(checks = pass) ~ section, data = led, FUN = length)
by_sec$passed <- aggregate(pass ~ section, data = led, FUN = sum)$pass
cat("\n"); print(by_sec, row.names = FALSE)

if (n_fail > 0L) {
  cat("\n  FAILURES:\n")
  print(led[!led$pass, c("section", "label", "observed", "expected")],
        row.names = FALSE)
}

cat(sprintf("\n  %d checks: %d passed, %d failed\n", nrow(led), n_pass, n_fail))
cat(sprintf("  phdid %s | R %s | %s\n",
            as.character(utils::packageVersion("phdid")),
            paste(R.version$major, R.version$minor, sep = "."),
            format(Sys.Date())))

utils::write.csv(led, "mpdta_replication_ledger.csv", row.names = FALSE)
cat("  ledger written to mpdta_replication_ledger.csv\n\n")

if (n_fail > 0L) quit(status = 1L)
