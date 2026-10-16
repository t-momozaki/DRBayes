# Reproduces the package's rows of Table 1 of Orihara, Momozaki and Sugasawa
# (2025), "Bayesian Doubly Robust Causal Inference via Posterior Coupling", by
# giving the package the settings the authors' own simulation code used.
#
# scripts/simulation_table1.R runs the same design with the package's defaults,
# which the paper does not use. This script runs it as the paper was run, so
# that the published numbers can be checked against the package:
#
#   * one chain of 21,000 iterations, the first 1000 discarded, so 20,000
#     posterior draws, as the Table 1 footnote states;
#   * outcome coefficients N(0, 20^2), error variance InvGamma(2.5, 2.5), and a
#     propensity score prior close to flat (precision 1e-8);
#   * Algorithm 2 on the grid 0, 0.05, ..., 5, stopping only when the posterior
#     mean of B_n changes sign (a tolerance of 1e-12 is never met first);
#   * the pruning rule of the authors' code, which is drbayes_control()'s
#     default rule;
#   * Saarela's estimator with 20,000 Dirichlet draws.
#
# None of these is stated in the paper, and several are poor defaults for
# general use, which is why the package does not adopt them: one chain cannot
# be checked for convergence, and the outcome prior is not scaled to the data.
#
# The DR and Luo rows are not package code; scripts/simulation_table1.R runs
# them. Section 5 used J = 2000 replications; the default here is 200, enough
# to compare with the published values within Monte Carlo error, and the
# authors did not fix their seed, so exact agreement is not expected at any J.
#
# Run it from the top of the package:
#
#   Rscript scripts/reproduce_table1.R
#
# Environment variables: DRBAYES_REPS (default 200), DRBAYES_N (default 500;
# Table 1 also reports 1500), DRBAYES_WORKERS (default 1), DRBAYES_OUTPUT_DIR
# (default tempdir()). One replication at n = 500 takes about two minutes on
# one core of an Apple M4 Max.
#
# References
#   Orihara, S., Momozaki, T., & Sugasawa, S. (2025). Bayesian Doubly Robust
#     Causal Inference via Posterior Coupling. arXiv:2506.04868.
#   Saarela, O., Belzile, L. R., & Stephens, D. A. (2016). A Bayesian view of
#     doubly robust causal inference. Biometrika, 103(3), 667-681.

library(DRBayes)

common <- "scripts/simulation_common.R"
if (!file.exists(common)) {
  stop("Run this script from the top of the package, as in ",
       "Rscript scripts/reproduce_table1.R: it reads ", common, ".",
       call. = FALSE)
}
source(common)

env_integer <- function(name, default) {
  value <- Sys.getenv(name, "")
  if (identical(value, "")) default else as.integer(strsplit(value, ",")[[1]])
}

n_replications <- env_integer("DRBAYES_REPS", 200L)
sample_sizes   <- env_integer("DRBAYES_N", 500L)
n_workers      <- env_integer("DRBAYES_WORKERS", 1L)
output_dir     <- Sys.getenv("DRBAYES_OUTPUT_DIR", unset = tempdir())
study_seed     <- 20261016L
true_ate       <- 110

# The authors' settings, given to drbayes_pc() as arguments.
paper_settings <- list(
  mc = 21000L, bn = 1000L, thin = 1L, chains = 1L,
  outcome.priors = list(theta_prior = 1 / 400, sigma_prior = c(2.5, 2.5)),
  ps.priors      = list(theta_prior = 1e-8),
  diagnostics    = "none"
)
paper_control <- function(pruning) {
  drbayes_control(lambda_max = 5, n_steps = 100, tol = 1e-12,
                  smoothing = 0.99, pruning = pruning)
}

outcome_formulas <- list(correct   = YY ~ AA + WW1 + WW2 + WW3 + WW4,
                         incorrect = YY ~ AA + WW1)
ps_formulas      <- list(correct   = AA ~ WW1 + WW2 + WW3 + WW4,
                         incorrect = AA ~ WW1)
specifications <- list(
  list(ps = "correct",   outcome = "correct",   pruning = FALSE),
  list(ps = "incorrect", outcome = "correct",   pruning = FALSE),
  list(ps = "correct",   outcome = "incorrect", pruning = TRUE)
)
block_label <- function(specification) {
  capitalise <- function(word) {
    paste0(toupper(substring(word, 1L, 1L)), substring(word, 2L))
  }
  paste(capitalise(specification$ps), "/", capitalise(specification$outcome))
}
block_labels <- vapply(specifications, block_label, character(1))
method_order <- c("G-formula", "Saarela", "Proposed", "Proposed (pruning)")

#' The package's estimators of Table 1, on one dataset under one specification
run_replication <- function(data, specification, seed) {
  outcome.formula <- outcome_formulas[[specification$outcome]]
  ps.formula      <- ps_formulas[[specification$ps]]
  couple <- function(pruning) {
    do.call(drbayes_pc, c(
      list(outcome.formula = outcome.formula, ps.formula = ps.formula,
           data = data, control = paper_control(pruning), seed = seed,
           verbose = FALSE),
      paper_settings))
  }

  with_seed(seed, {
    fit <- suppressWarnings(couple(FALSE))
    saarela <- drbayes_bb(outcome.formula, ps.formula, data = data,
                          num_iterations = 20000L, verbose = FALSE)
    estimates <- rbind(summarise_draws(fit$g.comp, "G-formula"),
                       summarise_draws(saarela, "Saarela"),
                       summarise_draws(fit$pc, "Proposed"))
    if (isTRUE(specification$pruning)) {
      pruned <- suppressWarnings(couple(TRUE))
      estimates <- rbind(estimates,
                         summarise_draws(pruned$pc, "Proposed (pruning)"))
    }
    estimates$block <- block_label(specification)
    estimates
  })
}

run_sample_size <- function(n) {
  message("Sample size ", n, ": ", n_replications, " replications.")
  datasets <- run_parallel_simulation(num_simulations = n_replications,
                                      nn = n, pp = 0, seed = study_seed,
                                      verbose = FALSE)
  seeds <- replication_seeds(study_seed, n_replications,
                             length(specifications))
  results <- map_replications(seq_len(n_replications), function(replication) {
    do.call(rbind, lapply(seq_along(specifications), function(k) {
      run_replication(datasets[[replication]], specifications[[k]],
                      seeds[replication, k])
    }))
  }, n_workers = n_workers, seed = study_seed)
  results <- do.call(rbind, results)
  results$n <- n
  results
}

started <- proc.time()
results <- do.call(rbind, lapply(sample_sizes, run_sample_size))
table1  <- performance_table(results, true_ate, block_labels, method_order)

utils::write.csv(results, file.path(output_dir, "reproduce-table1-replications.csv"),
                 row.names = FALSE)
utils::write.csv(table1, file.path(output_dir, "reproduce-table1.csv"),
                 row.names = FALSE)

cat("\nTable 1 with the authors' settings:", n_replications,
    "replications, the paper used 2000.\n\n")
print_performance_table(table1)
cat("\nWrote", file.path(output_dir, "reproduce-table1.csv"), "\n")
cat("Elapsed:", round((proc.time() - started)[["elapsed"]] / 60, 1), "minutes\n")
