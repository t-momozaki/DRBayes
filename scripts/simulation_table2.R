# Reproduces Table 2 of section 5.4 of Orihara, Momozaki and Sugasawa (2025):
# the same data generating mechanism as Table 1, but with the outcome model
# replaced by BART. Section 4.2 argues that coupling suppresses the initial
# variability of the outcome posterior and so improves the convergence rate of
# a nonparametric outcome model; Table 2 is where that is measured.
#
# ---------------------------------------------------------------------------
# PENDING. This script does not run yet, and has never been run.
#
# drbayes_pc() takes posterior draws of the outcome COEFFICIENTS. A BART
# posterior has no coefficients, only draws of the fitted mean, so coupling it
# needs an interface that takes those directly: outcome.mu, outcome.mu1 and
# outcome.mu0, each a draws by n matrix of posterior draws of the fitted mean
# at the observed covariates and at the two counterfactual treatment values.
# That interface is being added. The script is written against it and stops
# with an explanation until it exists, rather than being held back until then,
# so that the shape of what section 5.4 needs is on record.
# ---------------------------------------------------------------------------
#
# Run it, once the interface exists, with
#   Rscript scripts/simulation_table2.R
# from the top of the package, which is also where scripts/simulation_common.R
# is read from. It needs the dbarts package, which is not a dependency of
# DRBayes: the package fits parametric outcome models only.
# Set DRBAYES_OUTPUT_DIR to choose where the table is written.
#
# The published numbers, over 2000 replications, are printed again when the
# script finishes. Table 2 reports no interval, so no coverage probability and
# no average length: the comparison section 5.4 makes is about bias and
# variance.
#
# References
#   Chipman, H. A., George, E. I., & McCulloch, R. E. (2010). BART: Bayesian
#     additive regression trees. Annals of Applied Statistics, 4(1), 266-298.
#   Orihara, S., Momozaki, T., & Sugasawa, S. (2025). Bayesian Doubly Robust
#     Causal Inference via Posterior Coupling. arXiv:2506.04868.

library(DRBayes)

common <- "scripts/simulation_common.R"
if (!file.exists(common)) {
  stop("Run this script from the top of the package, as in ",
       "Rscript scripts/simulation_table2.R: it reads ", common, ".",
       call. = FALSE)
}
source(common)

if (!all(c("outcome.mu", "outcome.mu1", "outcome.mu0") %in%
           names(formals(drbayes_pc)))) {
  stop("drbayes_pc() does not yet accept outcome.mu, outcome.mu1 and ",
       "outcome.mu0, the posterior draws of a nonparametric outcome model's ",
       "fitted mean, so section 5.4 cannot be run. Nothing below has ever ",
       "been run either; treat it as a specification of what Table 2 needs, ",
       "not as a reproduction.", call. = FALSE)
}

if (!requireNamespace("dbarts", quietly = TRUE)) {
  stop("This script needs the dbarts package to fit the BART outcome model. ",
       "Install it with install.packages(\"dbarts\"). It is not a dependency ",
       "of DRBayes, which fits parametric outcome models only.", call. = FALSE)
}


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

# Section 5: J = 2000 replications.
n_replications <- 8L
# n_replications <- 2000L

# Table 2 reports three sample sizes side by side.
sample_sizes <- 200L
# sample_sizes <- c(200L, 500L, 1500L)

# The footnote to Table 2: "For all methods, 10,000 samples out of 20,000
# posterior draws are used". The propensity score model is sampled to the same
# length, so that the coupling has one particle per outcome draw.
n_bart_draws <- 2000L
# n_bart_draws <- 20000L
n_bart_burn <- 1000L

# Every second draw, which spreads the subsample over the whole chain.
used   <- seq(2L, n_bart_draws, by = 2L)
n_used <- length(used)

true_ate <- 110

study_seed <- 20250604L

n_workers <- 1L

output_dir <- Sys.getenv("DRBAYES_OUTPUT_DIR", unset = tempdir())

outcome.formula <- YY ~ AA + WW1 + WW2 + WW3 + WW4

# Table 2 reports the proposed method under a correct and under a misspecified
# propensity score model. The outcome model is nonparametric in both, so there
# is no outcome misspecification to vary.
ps_formulas <- list(
  "Proposed (PS correct)"   = AA ~ WW1 + WW2 + WW3 + WW4,
  "Proposed (PS incorrect)" = AA ~ WW1
)

method_order <- c("G-formula", names(ps_formulas))


# ---------------------------------------------------------------------------
# The BART outcome posterior
# ---------------------------------------------------------------------------

#' Posterior draws of the BART fitted mean, observed and counterfactual
#'
#' The three matrices are what a nonparametric outcome model can hand posterior
#' coupling in place of coefficients: the moment condition of equation (3.4)
#' needs the fitted mean at the observed treatment, and the g-computation
#' estimand of equation (3.7) needs it at both counterfactual values.
#'
#' The counterfactual predictions are made by setting the treatment in the DATA
#' and rebuilding the design matrix, so that any transformed or interacted term
#' is recomputed rather than left at its observed value.
#'
#' @param outcome.formula Formula for the outcome model.
#' @param data Data frame holding every variable in it.
#' @param treatment Name of the treatment column.
#' @param n_draws Number of posterior draws to keep.
#' @param n_burn Number of leading draws to discard.
#'
#' @return A list of three \code{n_draws} by \code{nrow(data)} matrices,
#'   \code{mu}, \code{mu1} and \code{mu0}.
#'
#' @references
#' Chipman, H. A., George, E. I., & McCulloch, R. E. (2010). BART: Bayesian
#' additive regression trees. Annals of Applied Statistics, 4(1), 266-298.
bart_outcome_draws <- function(outcome.formula, data, treatment,
                               n_draws, n_burn) {
  response <- stats::model.response(stats::model.frame(outcome.formula, data))

  design <- function(value) {
    modified <- data
    if (!is.null(value)) modified[[treatment]] <- value
    stats::model.matrix(outcome.formula, modified)[, -1L, drop = FALSE]
  }

  n <- nrow(data)
  fit <- dbarts::bart(x.train = design(NULL), y.train = response,
                      x.test  = rbind(design(1), design(0)),
                      ndpost = n_draws, nskip = n_burn,
                      keeptrees = FALSE, verbose = FALSE)

  list(mu  = fit$yhat.train,
       mu1 = fit$yhat.test[, seq_len(n), drop = FALSE],
       mu0 = fit$yhat.test[, n + seq_len(n), drop = FALSE])
}


#' Both rows of Table 2, on one dataset
#'
#' @param data One dataset from \code{generate_dataset()}.
#' @param seed Seed for this replication.
#'
#' @return A data frame with one row per method. \code{lower} and \code{upper}
#'   are NA because Table 2 reports no interval.
run_replication <- function(data, seed) {
  with_seed(seed, {
    draws <- bart_outcome_draws(outcome.formula, data, "AA",
                                n_bart_draws, n_bart_burn)

    # "10,000 samples out of 20,000 posterior draws are used".
    draws <- lapply(draws, function(m) m[used, , drop = FALSE])

    couple <- function(ps.formula) {
      # Algorithm 1, the importance sampling of section 3.3, which the
      # footnote to Table 2 says was used here in place of the sweep.
      drbayes_pc(outcome.formula = outcome.formula, ps.formula = ps.formula,
                 data = data, method = "is",
                 outcome.mu  = draws$mu,
                 outcome.mu1 = draws$mu1,
                 outcome.mu0 = draws$mu0,
                 mc = n_used + n_bart_burn, bn = n_bart_burn,
                 thin = 1L, chains = 1L,
                 verbose = FALSE)
    }

    fits <- lapply(ps_formulas, couple)

    # The untilted posterior does not depend on the propensity score model, so
    # the g-formula row is read off either fit.
    estimates <- rbind(
      summarise_draws(fits[[1L]]$g.comp, "G-formula"),
      do.call(rbind, Map(function(fit, method) {
        summarise_draws(fit$pc, method)
      }, fits, names(ps_formulas))))

    # Table 2 reports no interval, so the bounds are dropped rather than
    # reported as though the table had them.
    estimates$lower <- NA_real_
    estimates$upper <- NA_real_
    estimates$block <- "BART outcome model"
    estimates
  })
}


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

started <- Sys.time()

#' Run every replication at one sample size
#'
#' @param n Sample size.
#'
#' @return A data frame with one row per replication and method.
run_sample_size <- function(n) {
  message("Sample size ", n, ": generating ", n_replications, " datasets.")
  datasets <- run_parallel_simulation(num_simulations = n_replications,
                                      nn = n, pp = 0, seed = study_seed,
                                      verbose = FALSE)
  seeds <- replication_seeds(study_seed, n_replications)

  results <- map_replications(seq_len(n_replications), function(replication) {
    out <- run_replication(datasets[[replication]], seeds[replication, 1L])
    out$replication <- replication
    out
  }, n_workers = n_workers, seed = study_seed)

  results <- do.call(rbind, results)
  results$n <- n
  results
}

results <- do.call(rbind, lapply(sample_sizes, run_sample_size))

summary_table <- performance_table(results, true_ate, "BART outcome model",
                                   method_order)

cat("\nTable 2: BART outcome model, section 5.4.\n")
print_performance_table(summary_table, measures = c("ABias", "ESE", "RMSE"))

writeLines(c(
  "",
  "G-formula: the untilted posterior, the Bayesian g-formula over the BART",
  "  fitted means. Proposed: the same draws after coupling, under a correct",
  "  and under a misspecified propensity score model.",
  "",
  "Published, over 2000 replications:",
  "                          ABias    ESE   RMSE",
  "  n = 200  G-formula      0.428  0.949  1.042",
  "           Proposed (c.)  0.424  0.921  1.014",
  "           Proposed (i.)  0.422  0.921  1.013",
  "  n = 500  G-formula      0.189  0.365  0.410",
  "           Proposed (c.)  0.183  0.355  0.400",
  "           Proposed (i.)  0.181  0.355  0.398",
  "  n = 1500 G-formula      0.076  0.139  0.159",
  "           Proposed (c.)  0.067  0.138  0.153",
  "           Proposed (i.)  0.066  0.138  0.153"))

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

table_file <- file.path(output_dir, "table2.csv")
utils::write.csv(summary_table, table_file, row.names = FALSE)

draws_file <- file.path(output_dir, "table2-replications.csv")
utils::write.csv(results, draws_file, row.names = FALSE)

cat("\nWrote", table_file, "\n     ", draws_file, "\n")
cat("Elapsed:", format(round(Sys.time() - started, 1)), "\n")
