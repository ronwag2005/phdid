## Test environments

* local: macOS 26.6 (aarch64), R 4.6.1 -- `R CMD check --as-cran`
* GitHub Actions, all `--as-cran`, all OK:
  * windows-latest, R release
  * ubuntu-latest, R devel / release / oldrel-1
  * macos-latest, R release
* R-hub v2, all OK: linux (R-devel), windows (R-devel), macos (R-devel),
  ubuntu-clang
* win-builder, R-release (R 4.6.1 ucrt, Windows Server 2022): 1 NOTE, described
  below

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

Two R-hub containers did not produce a clean result, and neither reflects a
problem in this package.

`noSuggests` initially failed at "re-building of vignette outputs" with "there
is no package called 'rmarkdown'". That was a real defect in DESCRIPTION and is
fixed here: the vignettes use the `knitr::rmarkdown` engine, so per Writing R
Extensions both packages belong in `VignetteBuilder`, and the field now reads
`knitr, rmarkdown`.

`gcc15` fails while installing dependencies, before this package is checked at
all. The suggested package `did` cannot be built there because `vctrs` fails to
load on that toolchain ("symbol bindings not supported yet"). This is an
upstream incompatibility between `vctrs` and the experimental gcc15
configuration, unrelated to phdid.

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
