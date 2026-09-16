# Measures this package against every block of Table 1 of Orihara, Momozaki and
# Sugasawa (2025), "Bayesian Doubly Robust Causal Inference via Posterior
# Coupling", and prints the published value beside each number it computes.
#
# scripts/simulation_table1.R reproduces the whole of Table 1, including the two
# estimators that are not package code. This script is narrower and answers a
# different question: how far the four estimators the package does provide are
# from the published ones, measured against their own Monte Carlo error. It
# therefore reports three things Table 1 does not.
#
#   * The signed bias and its Monte Carlo standard error, so that a reader can
#     tell a real gap from noise. The standard error is the empirical standard
#     error over the square root of the number of replications, and it belongs
#     to the signed bias, not to the absolute one Table 1 reports. The two are
#     not interchangeable: where the true bias is zero, taking the absolute
#     value leaves ABias centred on about 0.8 MCSE rather than on zero, so an
#     ABias of that size is what no bias at all looks like.
#   * The ratio of average interval length to empirical standard error. An
#     equal-tailed 95 percent interval whose width is right for the sampling
#     variability of the estimator has AvL / ESE = 2 * qnorm(0.975) = 3.92, so
#     the ratio says whether the reported uncertainty is trustworthy in a way
#     that the bias and the standard error on their own do not. Coverage
#     conflates a wrong width with a wrong centre; this separates them.
#   * A direct check of Lemma 1 in the two blocks where the outcome model is
#     correctly specified. Lemma 1 says the tilted posterior converges to the
#     untilted one there, so fit$pc and fit$g.comp should agree in both location
#     and spread. That is a statement about the two distributions and not about
#     the tilting parameter, which can be large without the tilt mattering: with
#     a correct outcome model the moment condition has almost no posterior
#     spread, so a large lambda is what it takes to move its posterior mean by
#     one Monte Carlo standard error, and a large lambda applied to a nearly
#     degenerate quantity reweights almost nothing. lambda and the number of
#     steps are reported beside the two distances anyway, because a sweep that
#     reaches the end of its grid has stopped for want of grid and not because
#     the constraint was met. This is a regression test on the package rather
#     than a comparison with the paper, and it is reported separately.
#
# Run it with
#   Rscript scripts/validate_paper.R
# from the top of the package, which is also where scripts/simulation_common.R
# is read from. Set DRBAYES_OUTPUT_DIR to choose where the tables are written.
#
# The settings below are reduced from the ones the paper used, which are given
# beside each one and are commented out. It is a long run even so: two sample
# sizes by 200 replications by three specifications, each fitting two coupled
# models and one Bayesian bootstrap, is 1200 of each of those three fits.
# Nearly all of the cost is in the sequential Monte Carlo sweep, and the sweep
# is longest where the outcome model is correctly specified, for the reason the
# Lemma 1 note above gives. No wall clock figure is quoted, since it depends on
# the machine and on how many workers it has; the script prints its own elapsed
# time when it finishes, and the Lemma 1 table reports the median number of
# grid steps per fit for the blocks where the sweep is longest.
#
# The reduction is not free, and the output says where it bites: the sequential
# Monte Carlo sweep resamples its particles at every step, so the number of
# posterior draws feeds into the spread of the tilted posterior itself and not
# only into the Monte Carlo error of the summary.
#
# References
#   Orihara, S., Momozaki, T., & Sugasawa, S. (2025). Bayesian Doubly Robust
#     Causal Inference via Posterior Coupling. arXiv:2506.04868.
#   Saarela, O., Belzile, L. R., & Stephens, D. A. (2016). A Bayesian view of
#     doubly robust causal inference. Biometrika, 103(3), 667-681.
#   Daniels, M. J., Linero, A., & Roy, J. (2023). Bayesian Nonparametrics for
#     Causal Inference and Missing Data. Chapman and Hall/CRC.

library(DRBayes)

common <- "scripts/simulation_common.R"
if (!file.exists(common)) {
  stop("Run this script from the top of the package, as in ",
       "Rscript scripts/validate_paper.R: it reads ", common, ".",
       call. = FALSE)
}
source(common)


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

# Section 5.2: J = 2000 replications. The reduced setting leaves the Monte Carlo
# standard error of a bias at a fourteenth of the empirical standard error,
# which is small enough for the comparisons below to be read.
n_replications <- 200L
# n_replications <- 2000L

# Table 1 reports these two sample sizes side by side.
sample_sizes <- c(500L, 1500L)

# The footnote to Table 1: "For all Bayesian methods, 20,000 posterior draws
# are used". Four chains of 6000 iterations, the first 1000 discarded and none
# thinned away, leave exactly that many; the reduced setting leaves 4000.
mcmc <- list(mc = 1500L, bn = 500L, thin = 1L, chains = 4L)
# mcmc <- list(mc = 6000L, bn = 1000L, thin = 1L, chains = 4L)

# Section 5.1: Y_a = 100 + 110a + 13.7(2 X1 + X2 + X3 + X4) + eps, so the
# average treatment effect is 110 exactly.
true_ate <- 110

# Table 1 reports 95 percent credible intervals. That level is fixed inside
# summarise_draws() in scripts/simulation_common.R, which is where every
# interval below is formed, so it is not a setting this script can vary. A
# calibrated interval of that level is this many standard errors wide.
calibrated_ratio <- 2 * stats::qnorm(0.975)

# One seed controls the whole study, so a replication is reproducible on its
# own and the study returns the same numbers however many workers run it.
study_seed <- 20250604L

# Replications are independent, so they parallelise perfectly, and a study this
# size needs them to. One core is left free so that the machine stays usable.
# parallel ships with R, so reading the core count adds no dependency.
physical_cores <- parallel::detectCores(logical = FALSE)
n_workers <- if (is.na(physical_cores)) 2L else max(1L, physical_cores - 1L)

output_dir <- Sys.getenv("DRBAYES_OUTPUT_DIR", unset = tempdir())


# ---------------------------------------------------------------------------
# The three specifications of section 5.2
# ---------------------------------------------------------------------------
#
# "To evaluate the four methods, we consider three situations: 1) both the
# propensity score and outcome model is correctly specified, 2) only the
# propensity score model is correctly specified, and 3) only the outcome model
# is correctly specified. For misspecified model, only covariate X1 is used for
# each model."
#
# generate_dataset() names the covariates WW1 to WW4 and the treatment AA, so a
# misspecified model is one whose only covariate is WW1.

outcome_formulas <- list(
  correct   = YY ~ AA + WW1 + WW2 + WW3 + WW4,
  incorrect = YY ~ AA + WW1
)

ps_formulas <- list(
  correct   = AA ~ WW1 + WW2 + WW3 + WW4,
  incorrect = AA ~ WW1
)

# In the order the blocks appear in Table 1. Table 1 reports sample pruning in
# its last row only, for the block where the residual bias that pruning
# addresses appears; it is run in all three here, because a validation has to
# know what pruning does where the paper does not report it as well as where it
# does.
specifications <- list(
  list(ps = "correct",   outcome = "correct",   label = "C/C"),
  list(ps = "incorrect", outcome = "correct",   label = "I/C"),
  list(ps = "correct",   outcome = "incorrect", label = "C/I")
)

block_labels <- vapply(specifications, function(s) s$label, character(1))

# The blocks whose outcome model is correctly specified, which is the condition
# of Lemma 1.
lemma1_blocks <- vapply(
  specifications[vapply(specifications,
                        function(s) s$outcome == "correct", logical(1))],
  function(s) s$label, character(1))

method_order <- c("G-formula", "Saarela", "Proposed", "Proposed (pruning)")


# ---------------------------------------------------------------------------
# Table 1 as published
# ---------------------------------------------------------------------------
#
# Read off the table around line 486 of the arXiv v3 manuscript. Only the four
# rows this package provides are carried; the DR and Luo rows of the paper are
# not package code and scripts/simulation_table1.R is where they are run.
#
# Table 1 omits the G-formula row from the Incorrect / Correct block, since the
# g-formula never looks at the propensity score model and the row would repeat
# the Correct / Correct one, and it reports sample pruning for Correct /
# Incorrect only. Those cells are NA here and print as a dash.

published <- rbind(
  data.frame(n = 500L, block = "C/C", method = "G-formula",
             ABias = 0.001, ESE = 0.105, RMSE = 0.105, CP = 93.6, AvL = 0.397),
  data.frame(n = 500L, block = "C/C", method = "Saarela",
             ABias = 0.002, ESE = 0.113, RMSE = 0.113, CP = 92.5, AvL = 0.413),
  data.frame(n = 500L, block = "C/C", method = "Proposed",
             ABias = 0.001, ESE = 0.110, RMSE = 0.110, CP = 92.2, AvL = 0.396),
  data.frame(n = 500L, block = "I/C", method = "Saarela",
             ABias = 0.001, ESE = 0.108, RMSE = 0.108, CP = 93.3, AvL = 0.401),
  data.frame(n = 500L, block = "I/C", method = "Proposed",
             ABias = 0.000, ESE = 0.108, RMSE = 0.108, CP = 92.8, AvL = 0.396),
  data.frame(n = 500L, block = "C/I", method = "G-formula",
             ABias = 2.535, ESE = 2.286, RMSE = 3.414, CP = 80.8, AvL = 8.969),
  data.frame(n = 500L, block = "C/I", method = "Saarela",
             ABias = 0.008, ESE = 1.056, RMSE = 1.056, CP = 92.2, AvL = 3.807),
  data.frame(n = 500L, block = "C/I", method = "Proposed",
             ABias = 1.267, ESE = 1.271, RMSE = 1.795, CP = 98.9, AvL = 8.960),
  data.frame(n = 500L, block = "C/I", method = "Proposed (pruning)",
             ABias = 0.086, ESE = 1.428, RMSE = 1.430, CP = 98.0, AvL = 7.539),
  data.frame(n = 1500L, block = "C/C", method = "G-formula",
             ABias = 0.001, ESE = 0.060, RMSE = 0.060, CP = 94.7, AvL = 0.228),
  data.frame(n = 1500L, block = "C/C", method = "Saarela",
             ABias = 0.000, ESE = 0.064, RMSE = 0.064, CP = 94.8, AvL = 0.243),
  data.frame(n = 1500L, block = "C/C", method = "Proposed",
             ABias = 0.001, ESE = 0.061, RMSE = 0.061, CP = 94.2, AvL = 0.227),
  data.frame(n = 1500L, block = "I/C", method = "Saarela",
             ABias = 0.001, ESE = 0.062, RMSE = 0.062, CP = 94.8, AvL = 0.234),
  data.frame(n = 1500L, block = "I/C", method = "Proposed",
             ABias = 0.001, ESE = 0.061, RMSE = 0.061, CP = 94.2, AvL = 0.227),
  data.frame(n = 1500L, block = "C/I", method = "G-formula",
             ABias = 2.108, ESE = 1.315, RMSE = 2.485, CP = 63.4, AvL = 5.208),
  data.frame(n = 1500L, block = "C/I", method = "Saarela",
             ABias = 0.012, ESE = 0.572, RMSE = 0.572, CP = 92.4, AvL = 2.135),
  data.frame(n = 1500L, block = "C/I", method = "Proposed",
             ABias = 1.068, ESE = 0.735, RMSE = 1.296, CP = 97.0, AvL = 5.186),
  data.frame(n = 1500L, block = "C/I", method = "Proposed (pruning)",
             ABias = 0.154, ESE = 0.777, RMSE = 0.792, CP = 98.4, AvL = 4.340))


# ---------------------------------------------------------------------------
# One replication
# ---------------------------------------------------------------------------

#' The four package estimators of Table 1, on one dataset
#'
#' @param data One dataset from \code{generate_dataset()}.
#' @param specification One element of \code{specifications}.
#' @param seed Seed for this replication. Everything drawn below comes from it,
#'   so the rows this returns depend on nothing but its three arguments.
#'
#' @return A data frame with one row per estimator, carrying the posterior mean
#'   and interval of section 5.2 and, for the two coupled rows, the state of the
#'   sequential Monte Carlo sweep that produced them.
run_replication <- function(data, specification, seed) {
  outcome.formula <- outcome_formulas[[specification$outcome]]
  ps.formula      <- ps_formulas[[specification$ps]]

  couple <- function(control) {
    do.call(drbayes_pc, c(
      list(outcome.formula = outcome.formula, ps.formula = ps.formula,
           data = data, control = control, seed = seed, verbose = FALSE),
      mcmc))
  }

  with_seed(seed, {
    # Both coupled fits are given the same seed, so they share their Step 1
    # draws and the pruning row differs from the Proposed row through the sweep
    # alone. drbayes_pc() puts the caller's stream back on exit, so the
    # estimators after them are unaffected by either.
    fit    <- couple(drbayes_control())
    pruned <- couple(drbayes_control(pruning = TRUE))

    # Saarela's estimator is given the same number of draws as the posteriors
    # above, so that a comparison of interval lengths is not a comparison of
    # Monte Carlo error.
    saarela <- drbayes_bb(outcome.formula, ps.formula, data = data,
                          num_iterations = length(fit$pc), verbose = FALSE)

    estimates <- rbind(
      summarise_draws(fit$g.comp, "G-formula"),
      summarise_draws(saarela, "Saarela"),
      summarise_draws(fit$pc, "Proposed"),
      summarise_draws(pruned$pc, "Proposed (pruning)"))

    estimates$block     <- specification$label
    estimates$lambda    <- c(NA, NA, fit$smc$lambda, pruned$smc$lambda)
    estimates$n_steps   <- c(NA, NA, fit$smc$n_steps, pruned$smc$n_steps)
    estimates$converged <- c(NA, NA, isTRUE(fit$smc$converged),
                             isTRUE(pruned$smc$converged))

    # Lemma 1 is a statement about the whole tilted posterior, not only its
    # mean, so both the shift of the centre and the change of spread are
    # recorded. Each is in units of the untilted posterior's own standard
    # deviation, which is the scale on which "the two posteriors agree" means
    # anything.
    untilted_sd <- stats::sd(fit$g.comp)
    tilt_shift <- function(draws) {
      (mean(draws) - mean(fit$g.comp)) / untilted_sd
    }
    estimates$tilt_shift <- c(NA, NA, tilt_shift(fit$pc),
                              tilt_shift(pruned$pc))
    estimates$tilt_sd_ratio <- c(NA, NA, stats::sd(fit$pc) / untilted_sd,
                                 stats::sd(pruned$pc) / untilted_sd)

    # The two sweeps draw the same random numbers until the first step that
    # actually discards a particle, so identical draws mean pruning never fired.
    estimates$pruning_fired <- c(NA, NA, NA, !identical(pruned$pc, fit$pc))
    estimates
  })
}


#' Run every replication of every specification at one sample size
#'
#' @param n Sample size.
#'
#' @return A data frame with one row per replication, specification and
#'   estimator.
run_sample_size <- function(n) {
  message("Sample size ", n, ": ", n_replications, " replications of ",
          length(specifications), " specifications.")

  # pp = 0 gives the four covariates of section 5.1 and nothing else. The
  # irrelevant covariates belong to the high-dimensional setting of Appendix
  # D.1.
  datasets <- run_parallel_simulation(num_simulations = n_replications,
                                      nn = n, pp = 0, seed = study_seed,
                                      verbose = FALSE)

  seeds <- replication_seeds(study_seed, n_replications,
                             length(specifications))

  results <- map_replications(seq_len(n_replications), function(replication) {
    per_specification <- lapply(seq_along(specifications), function(k) {
      run_replication(datasets[[replication]], specifications[[k]],
                      seeds[replication, k])
    })
    out <- do.call(rbind, per_specification)
    out$replication <- replication
    out
  }, n_workers = n_workers, seed = study_seed)

  results <- do.call(rbind, results)
  results$n <- n
  results
}


# ---------------------------------------------------------------------------
# Summaries
# ---------------------------------------------------------------------------

#' Table 1's five measures, with the Monte Carlo error and the calibration
#'
#' @param results Output of \code{run_sample_size()}.
#'
#' @return A data frame with one row per sample size, block and method, holding
#'   the five measures of section 5.2, the signed bias, the Monte Carlo
#'   standard error of that bias, and the ratio of average interval length to
#'   empirical standard error.
validation_table <- function(results) {
  measures <- performance_table(results, true_ate, block_labels, method_order)

  keys <- list(results$n, results$block, results$method)
  extra <- do.call(rbind, lapply(split(results, keys, drop = TRUE),
                                 function(part) {
    data.frame(n      = part$n[1],
               block  = part$block[1],
               method = part$method[1],
               J      = nrow(part),
               Bias   = mean(part$estimate) - true_ate,
               MCSE   = stats::sd(part$estimate) / sqrt(nrow(part)),
               stringsAsFactors = FALSE)
  }))

  table <- merge(measures, extra, by = c("n", "block", "method"),
                 sort = FALSE)
  table$AvL_ESE <- table$AvL / table$ESE
  table[order(table$n,
              match(table$block, block_labels),
              match(table$method, method_order)), ]
}


#' Print each measured row above the published row it is meant to reproduce
#'
#' @param table Output of \code{validation_table()}.
#'
#' @return \code{table}, invisibly.
print_against_published <- function(table) {
  # MCSE is the standard error of Bias, the signed one, so the two are printed
  # side by side. ABias is |Bias| and is the measure Table 1 reports; MCSE is
  # not its standard error, and Bias / MCSE rather than ABias / MCSE is the z a
  # reader can act on. Table 1 gives no sign, so the paper rows leave Bias
  # blank.
  digits <- c(ABias = 3L, Bias = 3L, MCSE = 3L, ESE = 3L, RMSE = 3L, CP = 1L,
              AvL = 3L, AvL_ESE = 2L)

  format_rows <- function(part, source) {
    shown <- data.frame(Method = part$method, Source = source,
                        stringsAsFactors = FALSE)
    for (measure in names(digits)) {
      value <- part[[measure]]
      if (is.null(value)) value <- rep(NA_real_, nrow(part))
      shown[[measure]] <- ifelse(is.na(value), "-",
                                 formatC(value, format = "f",
                                         digits = digits[[measure]]))
    }
    shown
  }

  for (n in unique(table$n)) {
    for (block in block_labels) {
      part <- table[table$n == n & table$block == block, ]
      if (nrow(part) == 0) next

      # A published row exists for some methods and not others, so it is looked
      # up per method rather than merged, and stays a dash where Table 1 is
      # silent.
      paper <- published[published$n == n & published$block == block, ]
      paper <- paper[match(part$method, paper$method), ]
      paper$method <- part$method
      paper$AvL_ESE <- paper$AvL / paper$ESE

      cat("\nn = ", n, ", PS / outcome ", block, "\n", sep = "")
      rows <- rbind(format_rows(part, "package"), format_rows(paper, "paper"))
      # Two rows per method, the measured one above the published one.
      rows <- rows[order(match(rows$Method, method_order),
                         rows$Source == "paper"), ]
      print(rows, row.names = FALSE)
    }
  }

  invisible(table)
}


#' What the sweep did in the blocks Lemma 1 covers
#'
#' Lemma 1 says the tilted posterior converges to the untilted one when the
#' outcome model is correctly specified, so in those blocks the coupled draws
#' should sit on top of the g-formula draws in both location and spread. This
#' reports how far from that the sweep actually got, and what it cost in
#' tilting parameter and steps to get there.
#'
#' @param results Output of \code{run_sample_size()}.
#'
#' @return A data frame with one row per sample size, block and coupled method.
lemma1_table <- function(results) {
  covered <- results$block %in% lemma1_blocks & !is.na(results$lambda)
  coupled <- results[covered, ]
  keys <- list(coupled$n, coupled$block, coupled$method)

  rows <- lapply(split(coupled, keys, drop = TRUE), function(part) {
    data.frame(
      n            = part$n[1],
      block        = part$block[1],
      method       = part$method[1],
      shift        = mean(abs(part$tilt_shift)),
      shift_max    = max(abs(part$tilt_shift)),
      sd_ratio     = mean(part$tilt_sd_ratio),
      sd_ratio_min = min(part$tilt_sd_ratio),
      abs_lambda   = stats::median(abs(part$lambda)),
      steps        = stats::median(part$n_steps),
      converged    = 100 * mean(part$converged),
      stringsAsFactors = FALSE)
  })

  table <- do.call(rbind, rows)
  table[order(table$n, match(table$block, block_labels),
              match(table$method, method_order)), ]
}


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

started <- Sys.time()

results <- do.call(rbind, lapply(sample_sizes, run_sample_size))

summary_table <- validation_table(results)

writeLines(c(
  "",
  "Table 1 against this package. Each measured row is printed above the row",
  "of Table 1 it is meant to reproduce; a dash is a cell Table 1 does not",
  "report. C/C is a correct propensity score model with a correct outcome",
  "model, I/C an incorrect propensity score model with a correct outcome",
  "model, and C/I a correct propensity score model with an incorrect outcome",
  "model. Bias is the signed bias and MCSE its Monte Carlo standard error.",
  "ABias, the measure Table 1 reports, is the absolute value of that bias, and",
  "MCSE is not its standard error: where the true bias is zero, ABias averages",
  "about 0.8 MCSE rather than zero, so read Bias / MCSE and not ABias / MCSE.",
  "AvL_ESE is the ratio of average interval length to empirical standard",
  sprintf("error, %.2f for an interval whose width matches the sampling",
          calibrated_ratio),
  "variability of its own estimator."))

print_against_published(summary_table)

writeLines(c(
  "",
  "G-formula is the untilted posterior of drbayes_pc(), the Bayesian",
  "  g-formula of Daniels et al. (2023). Saarela is drbayes_bb(), the Bayesian",
  "  bootstrap estimator of Saarela et al. (2016). Proposed is the tilted",
  "  posterior of drbayes_pc(); (pruning) adds the sample pruning of section",
  "  5.3.1, drbayes_control(pruning = TRUE).",
  "",
  "Table 1 omits the G-formula row from the I/C block, since the g-formula",
  "  never looks at the propensity score model and the row would repeat the",
  "  C/C one. It is kept here as a check that the two blocks do agree.",
  "",
  "Table 1 reports sample pruning for the C/I block only. It is run in all",
  "  three blocks here, because what pruning does where the paper does not",
  "  report it is part of what a validation has to establish."))

writeLines(c(
  "",
  "Lemma 1: with the outcome model correctly specified the tilted posterior",
  "converges to the untilted one, so the coupled draws should sit on the",
  "g-formula draws. shift is the distance between the two posterior means in",
  "units of the untilted standard deviation, averaged over replications, and",
  "shift_max is its largest value; sd_ratio is the tilted standard deviation",
  "over the untilted one. Lemma 1 asks for a shift of 0 and an sd_ratio of 1.",
  "It does not ask for a small abs_lambda, which is reported here, as a median",
  "with the median number of steps, only because a sweep that reaches",
  "lambda_max has run out of grid rather than met the constraint. converged is",
  "the percentage of sweeps that met it."))

lemma1 <- lemma1_table(results)
print(data.frame(
  n = lemma1$n, Block = lemma1$block, Method = lemma1$method,
  shift        = formatC(lemma1$shift, format = "f", digits = 3),
  shift_max    = formatC(lemma1$shift_max, format = "f", digits = 3),
  sd_ratio     = formatC(lemma1$sd_ratio, format = "f", digits = 3),
  sd_ratio_min = formatC(lemma1$sd_ratio_min, format = "f", digits = 3),
  abs_lambda   = formatC(lemma1$abs_lambda, format = "f", digits = 2),
  steps        = lemma1$steps,
  converged    = formatC(lemma1$converged, format = "f", digits = 1),
  check.names = FALSE, stringsAsFactors = FALSE), row.names = FALSE)

# The doubly robust guarantee holds only where the sweep met the moment
# condition, so how often it did is part of the result.
coupled <- results[!is.na(results$converged), ]
for (method in intersect(method_order, unique(coupled$method))) {
  part <- coupled[coupled$method == method, ]
  cat("\nThe sweep met the moment condition in ", sum(part$converged), " of ",
      nrow(part), " ", method, " fits.", sep = "")
}

pruned_rows <- results[results$method == "Proposed (pruning)", ]
cat("\nPruning left the draws untouched altogether in ",
    sum(!pruned_rows$pruning_fired), " of ", nrow(pruned_rows), " fits.\n",
    sep = "")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

table_file <- file.path(output_dir, "validate-paper.csv")
utils::write.csv(summary_table, table_file, row.names = FALSE)

draws_file <- file.path(output_dir, "validate-paper-replications.csv")
utils::write.csv(results, draws_file, row.names = FALSE)

cat("\nWrote", table_file, "\n     ", draws_file, "\n")
cat("Elapsed:", format(round(Sys.time() - started, 1)), "\n")

draws_kept <- length(seq(mcmc$bn + 1L, mcmc$mc, by = mcmc$thin)) * mcmc$chains

writeLines(c(
  "",
  sprintf("These are %d replications, not the 2000 of the paper, and %d",
          n_replications, draws_kept),
  "posterior draws, not the 20,000 of its footnote. The MCSE column says what",
  "the replication count costs; the draw count is the more serious of the two",
  "and its cost does not appear in any column. The sweep resamples its",
  "particles at every step and rejuvenates them from a kernel fitted to the",
  "particles it has, so with fewer particles the tilted posterior is narrower",
  "than it should be however many replications are run. Read the AvL, CP and",
  "AvL_ESE of the two Proposed rows at the paper's draw count before",
  "concluding anything about them. The G-formula and Saarela rows are not",
  "affected: neither is touched by the sweep."))
