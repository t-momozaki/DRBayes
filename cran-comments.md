# DRBayes 0.1.1

This is a new submission: the first release of the package.

## Test environments

The package contains compiled code whose results are meant to match, bit for
bit, the R implementations it replaces. That equality depends on the width R
accumulates sums in and on which BLAS it reaches, neither of which varies on
the development machine, so the test suite was run on four combinations rather
than one.

| environment | `capabilities("long.double")` | BLAS | result |
| --- | --- | --- | --- |
| macOS 26, aarch64, R 4.5.1 | FALSE (`long double` is `double` here) | Accelerate | 2409 pass, 0 fail |
| Debian, aarch64, R 4.6.1 | TRUE, 128-bit | reference | 2371 pass, 0 fail |
| r-hub `nold`, x86_64, R-devel | FALSE, though the platform has an 80-bit `long double` | reference | 2371 pass, 0 fail |
| r-hub `clang-asan`, x86_64, R-devel | TRUE, 80-bit | reference | 2371 pass, 0 fail |

The `clang-asan` run is R built with `-fsanitize=address,undefined`, and it
reported nothing: no AddressSanitizer finding, no undefined behaviour, over the
whole suite. UndefinedBehaviorSanitizer was additionally run on macOS, where
AddressSanitizer cannot be used because its runtime must be loaded before the
process starts and injecting it with `DYLD_INSERT_LIBRARIES` hangs in dyld's
own initialiser; it reported nothing there either.

`tools/platform-checks/` in the sources holds the three container definitions
and the script, so these runs can be repeated.

Submitted to the two builders CRAN maintainers are asked to use before
submission, and checked on the five platforms GitHub Actions offers:

* win-builder, R-devel: 1 NOTE, the new submission one.
* win-builder, R 4.6.1: submitted.
* macOS builder: builds.
* GitHub Actions, R-CMD-check on macOS, Windows and Ubuntu at release, devel
  and oldrel-1: all pass.

The first win-builder run reported the `drbayes_pc` examples at 18.6 seconds
against the 10 second guidance. They run in 5.5 seconds on the development
machine, so the margin was not visible locally; they have been reduced, and
the one that supplies draws from elsewhere, the most expensive on the page,
moved to `\donttest`. What is timed now takes 1.7 seconds here, and the
subsequent R-devel run reports no examples NOTE.

## R CMD check results

    Status: 2 NOTEs
    0 errors | 0 warnings | 2 notes

Run with `_R_CHECK_FORCE_SUGGESTS_=false` and `_R_CHECK_LIMIT_CORES_=TRUE` on
macOS 26 (aarch64), R 4.5.1. Examples take 11s in total and the vignettes
rebuild in 68s, so the whole check finishes well inside the CRAN budget. The
PDF manual builds. Every example runs: there is no `\dontrun` in the package,
and the two `\donttest` blocks are the ones that need CmdStan or spawn
parallel workers.

The second of the two NOTEs is local: this machine's HTML Tidy is older than
the one the HTML validation step wants, so that step is skipped rather than
failing. It is not a property of the package.

## Notes we expect, and why

The NOTEs below follow from the package rather than from anything
fixable, and are listed so that any NOTE not among them can be treated as a
real finding.

* **New submission.** This is the first release of DRBayes, so the note
  "New submission" is expected.

* **Suggested package not in a mainstream repository: `cmdstanr`.**
  `bayes_stan()` fits the package's six models with CmdStan as an alternative
  to the built-in Gibbs samplers. `cmdstanr` is distributed from
  <https://stan-dev.r-universe.dev/>, which `Additional_repositories` in
  DESCRIPTION declares. It is used conditionally throughout:

  - every use is guarded by `requireNamespace()`, and a missing CmdStan
    produces an error naming both installation steps rather than a silent
    substitution;
  - the tests that need it call `skip_if_not_installed()` and skip on a machine
    without CmdStan;
  - the vignette chunks that compile or fit a Stan model are guarded by
    `instantiate::stan_cmdstan_exists()`, so the vignette knits on a machine
    with no CmdStan and shows the code without executing it.

  The package installs, checks and runs its whole documented workflow with
  neither `cmdstanr` nor CmdStan present. It is never the default backend, and
  will not become one, because CRAN's build machines have no CmdStan and the
  binaries CRAN distributes would therefore carry no compiled Stan models.

* **Possibly misspelled words in DESCRIPTION.** If the spell check flags
  "Bayesian", "Gibbs", "horseshoe", "probit", "Orihara", "Momozaki" or
  "Sugasawa", these are respectively standard statistical terminology and the
  surnames of the authors of the paper the package implements, all spelled as
  intended.

Anything else the check reports is a real finding and is not accounted for
here.

## Other submission notes

* The DOI in the Description field, `<doi:10.48550/arXiv.2506.04868>`, points
  to the arXiv preprint that the package implements. It resolves.

* The examples that fit a model use small sample sizes and short chains. The
  slowest are those of `summary.DRBayes`, `drbayes_sensitivity` and
  `print.DRBayes_sensitivity`, each about 7 to 8 seconds of elapsed time on the
  machine above; the rest are under 3 seconds. Longer demonstrations are inside
  `\dontrun{}`, and the paper reproduction scripts live in `scripts/`, which
  `.Rbuildignore` excludes from the build.

* No functions write to the user's file space, home directory or working
  directory. `bayes_stan()` caches its compiled Stan model under
  `tempdir()` unless the user names a directory with `cache_dir`.

* Functions that change graphics parameters restore them with `on.exit()`, and
  functions that take a `seed` restore the caller's random number stream.

* Both vignettes fit models while knitting. Over repeated builds on the machine
  above, which was running other jobs at the same time, "getting-started" took
  49 to 94 seconds and "external-software-integration" 83 to 134 seconds; the
  spread is contention, not variation in the work. Compiling and fitting the
  three Stan models of the second vignette adds about 20 seconds where CmdStan
  is present, and is skipped where it is not. The dominant cost in both is the
  sequential Monte Carlo sweep, whose number of steps depends on how far the
  tilting parameter has to travel for the data at hand.

* This is a new package, so there are no reverse dependencies to check.
