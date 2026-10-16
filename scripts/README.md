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
| `reproduce_table1.R` | Table 1, the package's rows, with the authors' settings | 20 at 200 replications on 12 cores |
| `simulation_table2.R` | Section 5.4, Table 2, BART outcome model | hours; see below |
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

## Reproducing the published numbers

The numbers in the paper come from the authors' own simulation code, whose
settings the paper does not state: one chain of 20,000 draws, an outcome prior
of N(0, 20^2) on the coefficients, a propensity score prior close to flat, the
grid 0, 0.05, ..., 5 for Algorithm 2 with a stop only when the mean of B_n
changes sign, and the pruning rule that `drbayes_control()` now uses by
default. `reproduce_table1.R` gives the package those settings. At n = 500 and
200 replications it reproduces every one of the package's rows of Table 1
within Monte Carlo error, the pruning row included:

| Block | Method | ABias | ESE | CP | AvL | published |
| --- | --- | --- | --- | --- | --- | --- |
| Correct / Correct | Proposed | 0.007 | 0.115 | 90.0 | 0.396 | 0.001, 0.110, 92.2, 0.396 |
| Incorrect / Correct | Proposed | 0.006 | 0.113 | 91.0 | 0.394 | 0.000, 0.108, 92.8, 0.396 |
| Correct / Incorrect | G-formula | 2.705 | 2.263 | 77.5 | 8.998 | 2.535, 2.286, 80.8, 8.969 |
| Correct / Incorrect | Saarela | 0.006 | 0.993 | 93.5 | 3.757 | 0.008, 1.056, 92.2, 3.807 |
| Correct / Incorrect | Proposed | 1.343 | 1.265 | 99.0 | 8.979 | 1.267, 1.271, 98.9, 8.960 |
| Correct / Incorrect | Proposed (pruning) | 0.103 | 1.439 | 99.0 | 7.555 | 0.086, 1.428, 98.0, 7.539 |

With 200 replications the Monte Carlo standard error of ABias in the
Correct / Incorrect block is 0.09 to 0.16 (ESE over the square root of 200),
and of CP up to about 2 points, so these differences are within it. The authors did not fix their seed, so exact agreement is not
expected. The package's defaults differ from these settings on purpose: one
chain cannot be checked for convergence and the outcome prior is not scaled to
the data, so `simulation_table1.R` keeps the defaults and `reproduce_table1.R`
holds the paper's settings.

## What is not reproduced

`simulation_table2.R` runs, through the `outcome.mu`, `outcome.mu1` and
`outcome.mu0` arguments of `drbayes_pc()`, which take draws of the fitted
means of a BART outcome model. It does not reproduce Table 2. Section 5.4
names no BART software, number of trees, burn-in, prior or number of chains.
The authors' own code, obtained since, used the BART package's wbart with its
default 200 trees, one chain of 20,000 draws after 5,000 burn-in, and the
treatment interacted with the covariates. This script instead uses dbarts,
runs four chains and sweeps the number of trees, because a single chain leaves
much of the reported spread to the sampler; its header sets out each choice.
Unlike the other scripts it takes hours rather than minutes at the full
settings; environment variables set the replications, sample sizes, tree
counts and workers, so a short check runs the same code.

With the package's defaults, the `Proposed (pruning)` rows of Table 1 and of
Table 3 are not reproduced; with the authors' settings Table 1's is, as above.
`simulation_table1.R` and `simulation_appendix_d2.R` run
`drbayes_control(pruning = TRUE)`, whose default is the
pruning rule of the authors' own simulation code, and print the row. At the
reduced settings the rule acts in most fits but moves the estimate by less than
its Monte Carlo error, so the row comes out on top of the unpruned one rather
than separating from it as the published pair does. Section 5.3.1 gives no
discard threshold; the rule's settings are left at the authors' values rather
than tuned until the published number appeared, which would be fitting the
setting to the answer. Each script says so where it prints the row.

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
