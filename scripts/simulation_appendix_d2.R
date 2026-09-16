# Reproduces Table 3 of Appendix D.2 of Orihara, Momozaki and Sugasawa (2025):
# the same study as Table 1, but with the outcome model misspecified in the
# severe way of Kang and Schafer (2007) rather than by dropping three
# covariates. The propensity score model stays correct, so the block is the one
# where posterior coupling has to do the work, and the sample pruning variant
# of section 5.3.1 is reported beside it.
#
# Run it with
#   Rscript scripts/simulation_appendix_d2.R
# from the top of the package, which is also where scripts/simulation_common.R
# is read from. As written it takes about half a minute on one core, 37 seconds
# on each of two runs of the machine it was written on.
# Set DRBAYES_OUTPUT_DIR to choose where the table is written.
#
# The published numbers, at n = 500 over 2000 replications:
#
#   Method               ABias    ESE   RMSE     CP    AvL
#   DR                   0.019  1.266  1.266  100.0  7.103
#   G-formula            3.439  1.494  3.750   37.7  5.927
#   Proposed             1.565  1.092  1.909   88.8  5.879
#   Proposed (pruning)   0.327  1.333  1.372   91.3  4.640
#
# The pruning row is run and printed but does not reproduce the published one:
# at the reduced settings below pruning discards little or nothing, so it comes
# out on top of the unpruned row instead of separating from it. The output says
# so where it prints the table.
#
# References
#   Kang, J. D. Y., & Schafer, J. L. (2007). Demystifying double robustness.
#     Statistical Science, 22(4), 523-539.
#   Orihara, S., Momozaki, T., & Sugasawa, S. (2025). Bayesian Doubly Robust
#     Causal Inference via Posterior Coupling. arXiv:2506.04868.

library(DRBayes)

common <- "scripts/simulation_common.R"
if (!file.exists(common)) {
  stop("Run this script from the top of the package, as in ",
       "Rscript scripts/simulation_appendix_d2.R: it reads ", common, ".",
       call. = FALSE)
}
source(common)


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

# Appendix D.2: "The other settings are the same as those in the simulation
# experiments in the main manuscript", so J = 2000 replications and 20,000
# posterior draws. Table 3 reports one sample size.
n_replications <- 8L
# n_replications <- 2000L

n_observations <- 500L

mcmc <- list(mc = 1200L, bn = 400L, thin = 2L, chains = 2L)
# mcmc <- list(mc = 6000L, bn = 1000L, thin = 1L, chains = 4L)

true_ate <- 110

study_seed <- 20250604L

n_workers <- 1L

output_dir <- Sys.getenv("DRBAYES_OUTPUT_DIR", unset = tempdir())


# ---------------------------------------------------------------------------
# Kang and Schafer's misspecification
# ---------------------------------------------------------------------------

#' Replace the four confounders by the transformed covariates of Appendix D.2
#'
#' Appendix D.2 hands the outcome model
#' \deqn{\left(\exp(X_1 / 2),\; 10 + X_2 / (1 + \exp(X_1)),\;
#'   (0.6 + X_1 X_3 / 25)^3,\; (20 + X_1 + X_4)^2\right)}
#' "after standardization" in place of the four covariates the outcome was
#' generated from. The propensity score model still sees the untransformed
#' ones, so this is a misspecification of the outcome model alone, and a far
#' harsher one than dropping three covariates: none of the four is a linear
#' function of what generated the outcome.
#'
#' Note for anyone comparing this with the wider literature: the fourth
#' transform of Kang and Schafer (2007) is \eqn{(20 + X_2 + X_4)^2}. The
#' expression above is the one Appendix D.2 prints, and is what this script
#' reproduces.
#'
#' @param data One dataset from \code{generate_dataset()}, with covariates
#'   named WW1 to WW4.
#'
#' @return \code{data} with four further columns ZZ1 to ZZ4, each standardised
#'   to mean zero and unit standard deviation.
#'
#' @references
#' Kang, J. D. Y., & Schafer, J. L. (2007). Demystifying double robustness.
#' Statistical Science, 22(4), 523-539.
add_misspecified_covariates <- function(data) {
  transformed <- data.frame(
    ZZ1 = exp(data$WW1 / 2),
    ZZ2 = 10 + data$WW2 / (1 + exp(data$WW1)),
    ZZ3 = (0.6 + data$WW1 * data$WW3 / 25)^3,
    ZZ4 = (20 + data$WW1 + data$WW4)^2
  )
  cbind(data, scale(transformed))
}

outcome.formula <- YY ~ AA + ZZ1 + ZZ2 + ZZ3 + ZZ4
ps.formula      <- AA ~ WW1 + WW2 + WW3 + WW4

method_order <- c("DR", "G-formula", "Proposed", "Proposed (pruning)")
block_label  <- "Correct / Incorrect"


# ---------------------------------------------------------------------------
# One replication
# ---------------------------------------------------------------------------

#' Every estimator of Table 3, on one dataset
#'
#' @param data One dataset from \code{generate_dataset()}, already carrying the
#'   transformed covariates.
#' @param seed Seed for this replication.
#'
#' @return A data frame with one row per estimator, carrying whether the sweep
#'   met the moment condition and, for the pruning row, whether pruning
#'   discarded anything.
run_replication <- function(data, seed) {
  couple <- function(control) {
    do.call(drbayes_pc, c(
      list(outcome.formula = outcome.formula, ps.formula = ps.formula,
           data = data, control = control, seed = seed, verbose = FALSE),
      mcmc))
  }

  with_seed(seed, {
    # Everything below draws from this one seed. drbayes_pc() puts the caller's
    # stream back on exit, so giving it the same seed twice hands both fits
    # identical Step 1 draws and the pruning row differs from the Proposed row
    # through the sweep alone.
    fit    <- couple(drbayes_control())
    pruned <- couple(drbayes_control(pruning = TRUE))

    estimates <- rbind(
      aipw_estimate(outcome.formula, ps.formula, data),
      summarise_draws(fit$g.comp, "G-formula"),
      summarise_draws(fit$pc, "Proposed"),
      summarise_draws(pruned$pc, "Proposed (pruning)"))

    estimates$block <- block_label
    estimates$smc_converged <- c(NA, NA, isTRUE(fit$smc$converged),
                                 isTRUE(pruned$smc$converged))
    # The two sweeps draw the same random numbers until the first step that
    # actually discards a particle, so identical draws mean pruning never
    # fired. That belongs beside the row it produced.
    estimates$pruning_fired <- c(NA, NA, NA, !identical(pruned$pc, fit$pc))
    estimates
  })
}


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

started <- Sys.time()

message("Generating ", n_replications, " datasets of ", n_observations,
        " observations.")
datasets <- run_parallel_simulation(num_simulations = n_replications,
                                    nn = n_observations, pp = 0,
                                    seed = study_seed, verbose = FALSE)
datasets <- lapply(datasets, add_misspecified_covariates)

seeds <- replication_seeds(study_seed, n_replications)

results <- map_replications(seq_len(n_replications), function(replication) {
  out <- run_replication(datasets[[replication]], seeds[replication, 1L])
  out$replication <- replication
  out
}, n_workers = n_workers, seed = study_seed)

results <- do.call(rbind, results)
results$n <- n_observations

summary_table <- performance_table(results, true_ate, block_label,
                                   method_order)

writeLines(c(
  "",
  "Table 3, Appendix D.2: outcome model misspecified in the manner of Kang",
  "and Schafer (2007), propensity score model correct."))
print_performance_table(summary_table)

writeLines(c(
  "",
  "Published, at n = 500 over 2000 replications:",
  "  DR                   0.019  1.266  1.266  100.0  7.103",
  "  G-formula            3.439  1.494  3.750   37.7  5.927",
  "  Proposed             1.565  1.092  1.909   88.8  5.879",
  "  Proposed (pruning)   0.327  1.333  1.372   91.3  4.640",
  "",
  "The reading Appendix D.2 draws from the table is that coupling cuts the",
  "  bias of the g-formula, that pruning cuts it further, and that neither",
  "  repairs the coverage: under this misspecification the credible interval",
  "  is in the wrong place, not merely too short."))

coupled <- results[!is.na(results$smc_converged), ]
for (method in unique(coupled$method)) {
  part <- coupled[coupled$method == method, ]
  cat("\nThe sequential Monte Carlo sweep met the moment condition in ",
      sum(part$smc_converged), " of ", nrow(part), " ", method, " fits.",
      sep = "")
}
cat("\n")

# How far the pruning row is from the row it is meant to improve on, measured
# against its own Monte Carlo error, is the whole of what that row is for.
pruned_rows <- results[results$method == "Proposed (pruning)", ]
plain_rows  <- results[results$method == "Proposed", ]
order_pruned <- order(pruned_rows$replication)
order_plain  <- order(plain_rows$replication)
gap <- pruned_rows$estimate[order_pruned] - plain_rows$estimate[order_plain]

writeLines(c(
  "",
  sprintf(paste0("Proposed (pruning) minus Proposed, paired by replication: ",
                 "%.3f, with a"), mean(gap)),
  sprintf(paste0("standard error of %.3f. Pruning left the draws untouched ",
                 "altogether in %d"),
          stats::sd(gap) / sqrt(length(gap)),
          sum(!pruned_rows$pruning_fired)),
  sprintf(paste0("of %d fits. Table 3 puts the absolute bias of the two rows ",
                 "at 1.565 and 0.327."), nrow(pruned_rows))))

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

table_file <- file.path(output_dir, "table3.csv")
utils::write.csv(summary_table, table_file, row.names = FALSE)

draws_file <- file.path(output_dir, "table3-replications.csv")
utils::write.csv(results, draws_file, row.names = FALSE)

cat("\nWrote", table_file, "\n     ", draws_file, "\n")
cat("Elapsed:", format(round(Sys.time() - started, 1)), "\n")
writeLines(c(
  "",
  sprintf("These are %d replications, not the 2000 of the paper, and %d",
          n_replications, mcmc$chains * length(seq(mcmc$bn + 1L, mcmc$mc,
                                                   by = mcmc$thin))),
  "posterior draws rather than 20,000. ESE, RMSE and CP carry too much Monte",
  "Carlo error to be read as a reproduction, and the spread of the tilted",
  "posterior depends on the number of draws as well; see the note at the end",
  "of scripts/simulation_table1.R.",
  "",
  "The Proposed (pruning) row does not reproduce the last row of Table 3, and",
  "the ABias column is not an exception to that: the published pair separates",
  "by 1.238 in ABias where the two rows above are nearly on top of each other.",
  "Pruning discards particles whose weight against the untilted posterior",
  "falls below prune.w / S, and at these settings that is little or nothing.",
  "Section 5.3.1 gives no threshold, so what counts as small is",
  "drbayes_control()'s reading rather than the paper's; ?drbayes_control says",
  "so under the pruning threshold. Choosing prune.w or prune.q to make the",
  "published number appear would be fitting the setting to the answer, so",
  "they are left at their defaults and the row is reported as not reproduced."))
