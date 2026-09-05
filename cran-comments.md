## Test environments

* local: macOS 26.6 (aarch64), R 4.6.1, `R CMD check --as-cran`

<!-- BEFORE SUBMITTING: run the checks below and list the results here.
     CRAN expects at least one Windows and one Linux result for a new
     submission; the lines are commented out until they have actually been run.

* win-builder: devel and release
     devtools::check_win_devel(); devtools::check_win_release()
* R-hub: windows-x86_64-devel, ubuntu-gcc-release, fedora-clang-devel
     rhub::rhub_check(platforms = c("windows", "ubuntu-release", "fedora"))
-->

## R CMD check results

0 errors | 0 warnings | 0 notes

Locally the check reports one NOTE, "Skipping checking HTML validation: 'tidy'
doesn't look like recent enough HTML Tidy" together with "Skipping checking
math rendering: package 'V8' unavailable". Both are properties of the local
machine rather than the package, and neither arises on a machine with a current
HTML Tidy and V8 installed.

## Notes for the reviewer

This is a new submission.

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
