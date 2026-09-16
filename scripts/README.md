# Reproduction scripts

Scripts that reproduce the results of Orihara, Momozaki and Sugasawa (2025),
*Bayesian Doubly Robust Causal Inference via Posterior Coupling*
([arXiv:2506.04868](https://doi.org/10.48550/arXiv.2506.04868)).

They are excluded from the package build by `.Rbuildignore` and are meant to
be read as well as run. Run them from the top of the package:

```sh
Rscript scripts/simulation_table1.R
```

Each is set up to finish in a few minutes on one core. The settings the paper
used sit next to the reduced ones, commented out, so what was run for
publication is visible without digging. At the reduced settings the numbers
carry far too much Monte Carlo error to be read as a reproduction, and each
script says so where it prints its results.

| Script | Reproduces | Minutes |
| --- | --- | --- |
| `simulation_table1.R` | Sections 5.1 to 5.3, Table 1, Figure 1 | 1.8-1.9 |
| `simulation_table2.R` | Section 5.4, Table 2, BART outcome model | pending |
| `simulation_appendix_d1.R` | Appendix D.1, high dimensions | 0.4 |
| `simulation_appendix_d2.R` | Appendix D.2, Table 3 | 0.6 |
| `application_rhc.R` | Section 7.3, right heart catheterization | 1.8-1.9 |

Appendix D.2 is Kang and Schafer's (2007) misspecification of the outcome
model. The minutes are wall clock over two runs of each script on one core of
one machine; repeated runs there have varied by up to a third, so read them as
approximate. The download `application_rhc.R` makes the first time is not
counted.

`simulation_common.R` is not a script. It holds what the simulation scripts
share: the two estimators of Table 1 that are not part of the package, the
five performance measures of section 5.2, and the seeding.

Every fit is run under a seed derived from the one `study_seed` at the top of
each script, so a replication is reproducible on its own and the studies
return the same numbers however many workers run them. Running a script twice
gives byte-identical tables.

## What is not reproduced

`simulation_table2.R` does not run yet. `drbayes_pc()` takes posterior draws of
the outcome model's coefficients, and a BART posterior has none; coupling it
needs an interface that takes draws of the fitted mean instead. The script is
written against that interface, stops with an explanation until it exists, and
has never been run.

The `Proposed (pruning)` rows of Table 1 and of Table 3 are not reproduced.
Both scripts do run `drbayes_control(pruning = TRUE)` and print the row, but
at the reduced settings the sweep discards little or nothing, so the row comes
out on top of the unpruned one rather than separating from it as the published
pair does. Section 5.3.1 gives no discard threshold, so what counts as a small
sampling weight is `drbayes_control()`'s reading rather than the paper's, and
tuning it until the published number appeared would be fitting the setting to
the answer. Each script says so where it prints the row.

`simulation_appendix_d1.R` reproduces the appendix's comparison but not its
standard deviations: both shrunk estimators come out well below the published
0.46 and 0.45, and the reduced draw count is not the reason. The script states
that rather than offering a cause it has not shown.

Section 6, the Health and Retirement Study application, has no script here
because the data cannot be downloaded without registration. `data-raw/hrs.R`
writes the derivation out in full for a reader who has the files.

Two of the six rows of Table 1 are not package code. `DR` is
`aipw_estimate()` and `Luo` is `luo_estimate()`, both in
`simulation_common.R`, and `Luo` is a reimplementation of the method described
in section 4.3 of the paper rather than the authors' own code. The other four
rows, `G-formula`, `Saarela`, `Proposed` and `Proposed (pruning)`, come from
the package.

## Environment variables

| Variable | Default | Used by |
| --- | --- | --- |
| `DRBAYES_OUTPUT_DIR` | `tempdir()` | where tables and figures go |
| `DRBAYES_DATA_DIR` | `tempdir()/drbayes-data` | where downloads are cached |

Both default to temporary directories, so nothing lands in the repository
unless you ask for it.

## Packages

The simulation scripts need furrr and future whatever `n_workers` is set to,
because `run_parallel_simulation()` generates the datasets through them. Both
are dependencies of DRBayes, so nothing further is needed.

`simulation_table2.R` additionally needs dbarts, which is not a dependency of
DRBayes: the package fits parametric outcome models only.

## Data

`data-raw/rhc.R` downloads the right heart catheterization data of Connors et
al. (1996) from <https://hbiostat.org/data/> and prepares it;
`scripts/application_rhc.R` calls it the first time it needs the data.

`data-raw/hrs.R` downloads nothing. Every Health and Retirement Study product
needs a registered account and the genetic file needs a further agreement, so
the script documents which products to obtain and derives the section 6 cohort
from files you already have.

Neither directory ships in the package tarball.
