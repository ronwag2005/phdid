# phdid

**Partial homogeneity in staggered difference-in-differences.**

<!-- badges: start -->
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE.md)
[![R >= 4.1](https://img.shields.io/badge/R-%3E%3D%204.1-blue.svg)](https://cran.r-project.org/)
<!-- badges: end -->

In a staggered DiD design the treatment effect is not a single number but a
vector of cohort-time effects (CATTs), one per cohort-time cell. Estimating
every one separately is unbiased but inefficient when some are in fact equal;
pooling them all into a single TWFE coefficient is efficient but biased
whenever the heterogeneity is genuine.

`phdid` treats the choice between these extremes as a **partition-selection
problem** on the cells, and implements the two estimators of Arora and Wagle
(2026) that solve it, along with the specification tests and diagnostics needed
to know whether either should be used at all.

---

## Installation

```r
# install.packages("remotes")
remotes::install_github("ronwag2005/phdid")
```

## The idea in sixty seconds

With `G` treated cohorts and `T` periods there are `K` cohort-time cells. The
true effect vector exhibits *partial homogeneity* if there is a partition of
those cells into `m < K` groups within which the effects are equal. Under the
true partition, the grouped estimator is both unbiased and strictly more
efficient than the flexible one -- in a balanced design, exactly `|C_p|` times
more efficient for every cell in a group of size `|C_p|`.

The partition is unknown, so the package recovers it two ways:

| | `l0_ph()` | `bayes_ph()` |
|---|---|---|
| **Returns** | one partition | a posterior over partitions |
| **Selects by** | BIC on the agglomeration path, or a fixed `lambda` | Dirichlet Process prior, collapsed Gibbs |
| **Intervals** | condition on the selected partition | marginalize over the partition |
| **Coverage when the partition is uncertain** | 0.57-0.62 | 0.79-0.81 |
| **Coverage when effects are separated** | 0.95 | 0.94 |

The two are not rivals: `l0_ph()` is the fixed-variance MAP of the same
Bayesian model, differing only in the partition prior. Use `l0_ph()` when you
want one interpretable grouping to report and `bayes_ph()` for inference.

## Quick start

```r
library(phdid)
library(did)
data(mpdta, package = "did")

# Any heterogeneity-robust first stage will do.
first <- att_gt(yname = "lemp", tname = "year", idname = "countyreal",
                gname = "first.treat", control_group = "notyettreated",
                data = mpdta, bstrap = FALSE)

d <- ph_data(first)   # keeps post-treatment cells + the exact joint covariance

# 1. Is there heterogeneity to recover at all?
homogeneity_test(d)
#> chi-squared = 25.29 on 6 df, p = 0.000302
#> heterogeneity share: 80%
#> Reading: recoverable heterogeneity.

# 2. One partition to report.
l0_ph(d)
#> Groups (m) : 4  (selected by BIC)
#> Mean variance ratio among pooled cells: 0.50

# 3. Inference that accounts for not knowing the partition.
fit <- bayes_ph(d, alpha = 1, iters = 20000, burn = 2000)
aggregate(fit, "overall")
#>     term  estimate  conf.low  conf.high
#>  overall   -0.0239   -0.0476    -0.0004

coclustering(fit)      # which groupings are firm
plot(fit)              # the co-clustering heatmap
```

## Two things the package will keep telling you

**The covariance is not diagonal.** Cohort-time cells share units and share the
unit and time fixed effects, so their estimates are correlated. Using the exact
joint covariance is what makes the reported intervals honest: treating the
cells as independent understates the posterior variance of aggregates by about
3.5x in the paper's design and can misstate individual co-clustering
probabilities by up to 0.39. `ph_data()` warns when handed a diagonal
covariance, and `covariance_check()` quantifies the difference on your design.

**A partition is only worth reporting when the effects are separated enough to
be recovered.** Below that threshold both estimators return groupings that are
quantiles of noise. In the paper's second application the l0 estimator splits
138 event effects into five neat "bands" while a randomization test cannot
reject a single common effect. Run `homogeneity_test()` first; it says so
plainly rather than letting the output speak for itself.

## What's included

**Estimation** -- `ph_data()` (three entry points: a `did::att_gt()` fit, a
micro panel, or any `(tau, Sigma)` pair), `l0_ph()`, `bayes_ph()`, `ph_fit()`,
`flex_twfe()`, `pooled_twfe()`.

**Inference** -- `aggregate()` for overall ATT, event-study, by-cohort and
by-calendar estimands computed per posterior draw; `coclustering()`;
`point_partition()` (Wade-Ghahramani VI and Binder losses); `confint()`.

**Specification testing** -- `homogeneity_test()`, bundling a model-implied
common-effect chi^2 test, an excess-dispersion heterogeneity share, and the
placebo gauge and within-cell randomization test.

**Diagnostics** -- `enumerate_partitions()` for the exact posterior on small
designs, `ph_rhat()` for multi-chain convergence, `alpha_sensitivity()` for the
prior path, `covariance_check()` for exact vs. diagonal.

**Simulation** -- `ph_design()`, `ph_truth()`, `ph_sample()`, `sim_study()`,
reproducing the paper's Monte Carlo tables.

## Validation

The package is checked against the paper rather than against stored output. On
the Callaway-Sant'Anna minimum-wage data it reproduces Table 6 exactly (all
seven cells, the four selected groups, the grouped effects, and every variance
ratio), and the sampler matches the Appendix E exact-enumeration benchmark
(posterior E[groups] 2.20, overall ATT -0.024 with interval [-0.048, -0.000]).
The test suite asserts analytic identities -- that the singleton partition
returns the first stage unchanged, that pooling in a balanced orthogonal design
divides the variance by the group size, that the micro-panel and two-stage
routes agree exactly (Lemma 1), and that the Gibbs sampler agrees with exact
enumeration -- rather than golden numbers.

## Authorship and citation

The methods implemented here are joint work by **Parush Arora** and **Rohan
Wagle** (Department of Economics, Ashoka University). The simulation and
application code that this package generalizes was written by Parush Arora; the
package itself is written and maintained by Rohan Wagle. Both are copyright
holders under the MIT license.

If you use `phdid`, please cite the paper:

> Arora, Parush and Rohan Wagle (2026). "A Bayesian Approach to Partial
> Homogeneity in Staggered Difference-in-Difference." Ashoka University
> Economics Discussion Paper 166.
> [Ashoka](https://www.ashoka.edu.in/research/a-bayesian-approach-to-partial-homogeneity-in-staggered-difference-in-difference/) | [SSRN](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=7207083)

```r
citation("phdid")
```

## References

Arora, P. and Wagle, R. (2026). A Bayesian Approach to Partial Homogeneity in
Staggered Difference-in-Difference. Ashoka University Economics Discussion
Paper 166.
<https://www.ashoka.edu.in/research/a-bayesian-approach-to-partial-homogeneity-in-staggered-difference-in-difference/>
<https://papers.ssrn.com/sol3/papers.cfm?abstract_id=7207083>

Callaway, B. and Sant'Anna, P. H. C. (2021). Difference-in-Differences with
Multiple Time Periods. *Journal of Econometrics* 225(2), 200-230.

Neal, R. M. (2000). Markov Chain Sampling Methods for Dirichlet Process Mixture
Models. *JCGS* 9(2), 249-265.

Wade, S. and Ghahramani, Z. (2018). Bayesian Cluster Analysis: Point Estimation
and Credible Balls. *Bayesian Analysis* 13(2), 559-626.

Wooldridge, J. M. (2025). Two-Way Fixed Effects, the Two-Way Mundlak
Regression, and Difference-in-Differences Estimators. *Empirical Economics* 69,
2545-2587.

## License

MIT (c) 2026 Rohan Wagle and Parush Arora. See [LICENSE.md](LICENSE.md).
