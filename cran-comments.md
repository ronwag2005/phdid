## Test environments

All of the following were run on the submitted sources and returned OK.

* local: macOS 26.6 (aarch64), R 4.6.1 -- `R CMD check --as-cran`
* GitHub Actions, `--as-cran`, erroring on warnings:
  * windows-latest, R release
  * ubuntu-latest, R devel
  * ubuntu-latest, R release
  * ubuntu-latest, R oldrel-1
  * macos-latest, R release
* R-hub v2: windows (R-devel), linux (R-devel), macos (R-devel), ubuntu-clang,
  noSuggests
* win-builder: R-devel and R-release

## R CMD check results

0 errors | 0 warnings | 1 NOTE

The NOTE on win-builder is the CRAN incoming-feasibility check, and has three
parts:

* "New submission" -- expected; this is the package's first release.
* "Possibly misspelled words in DESCRIPTION: Arora, Wagle" -- these are the
  surnames of the paper's authors, spelled correctly.
* "Found the following (possibly) invalid file URI: LICENSE.md from README.md"
  -- fixed. README.md now links to the file by its full URL rather than by a
  relative path, since LICENSE.md is not shipped in the tarball.

Locally the check reports a different NOTE, "Skipping checking HTML validation: 'tidy'
doesn't look like recent enough HTML Tidy" together with "Skipping checking
math rendering: package 'V8' unavailable". Both are properties of the local
machine rather than the package, and neither arises on a machine with a current
HTML Tidy and V8 installed.

## Notes for the reviewer

This is a new submission.

Earlier drafts of this submission were checked and revised in response to three
findings, all now fixed and re-verified:

* R-hub `noSuggests` failed at "re-building of vignette outputs" with "there is
  no package called 'rmarkdown'". The vignettes use the `knitr::rmarkdown`
  engine, so per Writing R Extensions both packages belong in
  `VignetteBuilder`; the field now reads `knitr, rmarkdown`, and `noSuggests`
  passes.
* win-builder R-devel reported `ph_rhat` taking 14.4s in examples. The example
  ran four MCMC chains, which is inexpensive locally but not on Windows. All
  MCMC-bearing examples were trimmed; the slowest now runs in well under a
  second locally and the whole examples block takes about 3s.
* win-builder reported an invalid file URI for `LICENSE.md` referenced from
  `README.md`. That file is in `.Rbuildignore`, so the relative link could not
  resolve in the tarball; the README now uses the full URL.

Vignette build time was also reduced from about 20 minutes to under a minute by
shortening the MCMC chains and Monte Carlo replication counts used for
illustration. The settings that reproduce the paper's published tables remain
documented inside the vignette.

One R-hub container, `gcc15`, fails while installing dependencies, before this
package is checked at all: the suggested package `did` cannot be built there
because `vctrs` fails to load on that toolchain ("symbol bindings not supported
yet"). This is an upstream incompatibility with the experimental gcc15
configuration and is unrelated to phdid.

The `Description` field cites the working paper the methods come from. It has
no DOI yet, so the field carries its two stable landing pages in angle brackets
as `<https://...>` per the Writing R Extensions guidance on references without
a DOI. Both URLs point to the same paper, one at the publishing institution and
one on SSRN.

Long-running examples are wrapped in `\donttest{}`. These are the Monte Carlo
driver `sim_study()` and the plot method that consumes its output; a single
replication runs a full Gibbs chain, so an honest example exceeds five seconds.
All other examples run in well under a second, using short chains
(`iters = 400-500`) that are adequate to demonstrate the interface.

The package writes no files. Examples, tests and vignettes create nothing
outside `tempdir()`, and no function writes to the user's home or working
directory. `inst/replication/mpdta_replication.R` does write figures and a CSV
of its results, but it is a script shipped for reproducibility and is never
executed by `R CMD check`.

Functions that set graphical parameters restore them with
`on.exit(par(op))`. No function calls `set.seed()` unless the user supplies a
`seed` argument.

`did` is used only in one vignette and one optional entry point, and is
declared in `Suggests` with the vignette guarded by `requireNamespace()`.
