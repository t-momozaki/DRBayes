# Reproduces Table 2 of section 5.4 of Orihara, Momozaki and Sugasawa (2025):
# the data generating mechanism of section 5.1 with the outcome model replaced
# by BART. Section 4.2 argues that coupling suppresses the initial variability
# of a nonparametric outcome posterior and so improves its convergence rate;
# Table 2 is where that is measured.
#
# Run it with
#   Rscript scripts/simulation_table2.R
# from the top of the package, which is also where scripts/simulation_common.R
# is read from. It needs the dbarts package, which is not a dependency of
# DRBayes: the package fits parametric outcome models only.
#
# The settings below follow the paper wherever it states one; what it leaves
# open, and what was chosen instead, is set out next. Nothing here should be
# edited to run a smaller study: the replications, sample sizes, numbers of
# trees and workers can each be overridden by an environment variable, listed
# with the setting, so that a quick check and the real run take the same code
# path.
#
# What the paper leaves open, and what this script therefore has to choose
# -----------------------------------------------------------------------
# Section 5.4 is five sentences and the Table 2 footnote two lines. Between
# them they name no software, no number of trees, no draws, no burn-in and no
# prior hyperparameter, and never say whether BART was fit once with the
# treatment as a covariate or separately within each arm. Three consequences:
#
#   * Every BART hyperparameter below is written out even where it equals the
#     dbarts default, so that the record is complete.
#   * The number of trees is swept rather than fixed. With four chains, neither
#     count tried reproduces the published spread: at n = 200 the published
#     empirical standard error 0.949 compares with 0.641 at 20 trees and 0.414
#     at the dbarts default of 200 (2,000 replications at 20 trees and 500 at
#     200, every other setting as below).
#   * "Importance sampling without replacement" does not name either algorithm
#     of section 3.3: Algorithm 1 has no resampling step and Algorithm 2
#     resamples with replacement. All three readings are computed from the same
#     draws and reported side by side.
#
# The reported spread is a property of the sampler unless BART is run with more
# than one chain. On this design at n = 200, the posterior mean of the ATE
# moves by 0.473 between chains on the same data while the sampling variability
# across data sets is 0.395, and the effective sample size of the ATE along one
# chain is 11 draws out of 1000. Four chains are therefore not optional here.
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

if (!requireNamespace("dbarts", quietly = TRUE)) {
  stop("This script needs the dbarts package to fit the BART outcome model. ",
       "Install it with install.packages(\"dbarts\"). It is not a dependency ",
       "of DRBayes, which fits parametric outcome models only.", call. = FALSE)
}


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

#' One setting, overridable by an environment variable
#'
#' @param name Environment variable to read.
#' @param default Value to use when it is unset.
#' @param as Coercion to apply to the variable's value.
#'
#' @return The setting.
setting <- function(name, default, as = as.integer) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) default else as(strsplit(value, ",")[[1]])
}

# Section 5: J = 2000 replications, three sample sizes.
n_replications <- setting("DRBAYES_T2_REPS", 2000L)
sample_sizes   <- setting("DRBAYES_T2_N", c(200L, 500L, 1500L))

# The number of trees, which the paper does not state. See the header.
tree_counts <- setting("DRBAYES_T2_NTREE", c(200L, 20L))

# The Table 2 footnote: "For all methods, 10,000 samples out of 20,000
# posterior draws are used". Four chains, so that the spread reported is the
# estimator's and not the sampler's, and 2,500 draws from each after 1,000
# discarded.
#
# That is 10,000 draws, not 20,000. Generating 20,000 and keeping every second
# one costs twice as much and buys nothing: thinning a chain by two is not a
# better 10,000 draws than sampling 10,000. What the footnote's ratio does
# describe is the resampling step, and the half-without-replacement reading
# below keeps it, taking 5,000 of these 10,000. The Monte Carlo error of the
# draw budget is not what limits this study in any case: the effective sample
# size of the ATE along one chain is about 11 draws, so the chains are the
# binding constraint and four of them is the setting that matters.
n_bart_chains <- 4L
n_bart_burn   <- 1000L
n_bart_draws  <- 2500L                                  # per chain
n_total       <- n_bart_chains * n_bart_draws           # 10,000
used          <- seq_len(n_total)
n_used        <- length(used)

# Chipman, George and McCulloch's own defaults, written out rather than left
# implicit: the paper states none, so the record has to come from here.
bart_k      <- 2.0
bart_power  <- 2.0
bart_base   <- 0.95
bart_numcut <- 100L

# The propensity score posterior, matched to the outcome posterior: four
# chains, 1,000 discarded from each, 2,500 kept from each.
ps_chains <- 4L
ps_burn   <- 1000L
ps_mc     <- ps_burn + n_used / ps_chains

true_ate   <- 110
study_seed <- 20250604L

# Twelve of this machine's sixteen cores are performance cores, and a worker
# holds about two gigabytes at n = 1500. Ten leaves room for both.
n_workers  <- setting("DRBAYES_T2_WORKERS", 10L)

# Data sets and per-replication seeds are always drawn for the paper's J,
# whatever J this run uses, and the first n_replications of them are taken. A
# shorter run is then a prefix of the longer one rather than a different study,
# so it can be extended later by running only the replications it lacks.
n_pool <- 2000L

output_dir <- Sys.getenv("DRBAYES_OUTPUT_DIR", unset = tempdir())

outcome.formula <- YY ~ AA + WW1 + WW2 + WW3 + WW4

# Table 2 reports the proposed method under a correct and under a misspecified
# propensity score model. The outcome model is nonparametric in both, so there
# is no outcome misspecification to vary.
ps_formulas <- list(
  "PS correct"   = AA ~ WW1 + WW2 + WW3 + WW4,
  "PS incorrect" = AA ~ WW1
)

# The three readings of "importance sampling without replacement", in the order
# they are reported.
resamplings <- c("with replacement", "weighted mean", "half without replacement")

method_order <- c("G-formula",
                  as.vector(t(outer(names(ps_formulas), resamplings,
                                    function(p, r) paste0("Proposed (", p, ", ", r, ")")))))


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
#' @param ntree Number of trees.
#'
#' @return A list of three \code{n_used} by \code{nrow(data)} matrices,
#'   \code{mu}, \code{mu1} and \code{mu0}.
bart_outcome_draws <- function(outcome.formula, data, treatment, ntree) {
  response <- stats::model.response(stats::model.frame(outcome.formula, data))

  design <- function(value) {
    modified <- data
    if (!is.null(value)) modified[[treatment]] <- value
    stats::model.matrix(outcome.formula, modified)[, -1L, drop = FALSE]
  }

  n <- nrow(data)
  fit <- dbarts::bart(x.train = design(NULL), y.train = response,
                      x.test  = rbind(design(1), design(0)),
                      ntree = ntree, k = bart_k, power = bart_power,
                      base = bart_base, numcut = bart_numcut,
                      nchain = n_bart_chains, ndpost = n_bart_draws,
                      nskip = n_bart_burn,
                      keeptrees = FALSE, verbose = FALSE)

  # dbarts stacks the chains for us, returning draws by observations with the
  # four chains one after another; the guard is for the array layout older
  # versions return. Either way the 10,000 draws are one pool, 2,500 from each
  # of the four chains, and all of them are used.
  flatten <- function(x) if (length(dim(x)) == 3L) {
    matrix(aperm(x, c(2L, 1L, 3L)), nrow = dim(x)[1] * dim(x)[2])
  } else x

  train <- flatten(fit$yhat.train)[used, , drop = FALSE]
  test  <- flatten(fit$yhat.test)[used, , drop = FALSE]

  list(mu  = train,
       mu1 = test[, seq_len(n), drop = FALSE],
       mu0 = test[, n + seq_len(n), drop = FALSE])
}


# ---------------------------------------------------------------------------
# The moment condition, recomputed so that the weights are available
# ---------------------------------------------------------------------------

#' Equation (3.4) at every posterior draw
#'
#' \code{drbayes_pc()} reports the tilting parameter and the posterior mean of
#' the moment condition, but not its value at each draw, which is what the
#' other two readings of the footnote need. It is recomputed here from the
#' propensity score draws the fit returns, and checked against the fit: at the
#' tilting parameter the fit solved for, these weights must average the moment
#' condition to zero.
#'
#' @param alpha Draws by parameters matrix of propensity score coefficients.
#' @param Zps Design matrix of the propensity score model, intercept included.
#'   Both come from the fit, not from a second construction of them.
#' @param A,Y Treatment and outcome vectors.
#' @param mu Draws by observations matrix of fitted outcome means.
#' @param block Number of draws to hold at once.
#'
#' @return A vector with one value per draw.
moment_draws <- function(alpha, Zps, A, Y, mu, block = 1000L) {
  stopifnot(ncol(alpha) == ncol(Zps), nrow(alpha) == nrow(mu),
            ncol(mu) == length(Y))
  # (A - e)/(e(1 - e)) with e = plogis(eta), written so that no division by a
  # saturated probability can arise: s(1 + exp(-s * eta)) with s = 2A - 1.
  s <- 2 * A - 1
  out <- numeric(nrow(alpha))
  for (from in seq(1L, nrow(alpha), by = block)) {
    rows <- from:min(from + block - 1L, nrow(alpha))
    eta  <- tcrossprod(alpha[rows, , drop = FALSE], Zps)
    w    <- rep(s, each = length(rows)) * (1 + exp(-rep(s, each = length(rows)) * eta))
    out[rows] <- rowMeans(w * (rep(Y, each = length(rows)) - mu[rows, , drop = FALSE]))
  }
  out
}


# ---------------------------------------------------------------------------
# One replication
# ---------------------------------------------------------------------------

#' Every row of Table 2, on one dataset, at one number of trees
#'
#' @param data One dataset from \code{generate_dataset()}.
#' @param ntree Number of trees.
#' @param seed Seed for this replication.
#'
#' @return A data frame with one row per method. \code{lower} and \code{upper}
#'   are NA because Table 2 reports no interval.
run_replication <- function(data, ntree, seed) {
  with_seed(seed, {
    draws <- bart_outcome_draws(outcome.formula, data, "AA", ntree)
    ate   <- rowMeans(draws$mu1) - rowMeans(draws$mu0)   # equation (3.7)

    one_ps <- function(ps.formula, label) {
      # Algorithm 1, the importance sampling of section 3.3, which the footnote
      # to Table 2 says was used here in place of the sweep. The package
      # resamples systematically and with replacement.
      fit <- drbayes_pc(outcome.formula = outcome.formula,
                        ps.formula = ps.formula, data = data, method = "is",
                        outcome.mu  = draws$mu,
                        outcome.mu1 = draws$mu1,
                        outcome.mu0 = draws$mu0,
                        mc = ps_mc, bn = ps_burn, thin = 1L,
                        chains = ps_chains, diagnostics = "none",
                        control = drbayes_control(keep_particles = TRUE),
                        verbose = FALSE)

      lambda <- fit$smc$lambda
      # Read from the fit rather than rebuilt, so that no difference in how the
      # design matrix is formed can separate these weights from the fit's.
      td <- fit$particles$model_data
      BB <- moment_draws(fit$particles$ps_untilted, td$Z.ps, td$A, td$Y,
                         td$outcome.mu)

      # The fit reports the moment condition averaged over the draws it kept.
      # Averaging the recomputed values over the same draws has to return it;
      # if it does not, these are not the fit's weights.
      if (!isTRUE(all.equal(mean(BB[fit$particles$index]), fit$smc$B_mean))) {
        stop("The recomputed moment condition does not agree with the fit.",
             call. = FALSE)
      }
      ww <- exp(lambda * BB - max(lambda * BB))

      S <- length(ate)
      estimates <- c(
        "with replacement"         = mean(fit$pc),
        "weighted mean"            = sum(ww * ate) / sum(ww),
        "half without replacement" =
          mean(ate[sample.int(S, S %/% 2L, replace = FALSE, prob = ww)]))

      data.frame(method   = paste0("Proposed (", label, ", ", names(estimates), ")"),
                 estimate = as.vector(estimates),
                 lambda   = lambda,
                 no_op    = identical(fit$pc, fit$g.comp),
                 stringsAsFactors = FALSE)
    }

    rbind(
      data.frame(method = "G-formula", estimate = mean(ate),
                 lambda = NA_real_, no_op = NA, stringsAsFactors = FALSE),
      do.call(rbind, Map(one_ps, ps_formulas, names(ps_formulas))))
  })
}


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

started <- Sys.time()

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#' Run every replication at one sample size and number of trees
#'
#' @param n Sample size.
#' @param ntree Number of trees.
#'
#' @return A data frame with one row per replication and method.
run_cell <- function(n, ntree) {
  message("n = ", n, ", ", ntree, " trees: ", n_replications, " replications.")
  datasets <- run_parallel_simulation(num_simulations = n_pool, nn = n, pp = 0,
                                      seed = study_seed, verbose = FALSE)
  datasets <- datasets[seq_len(n_replications)]
  seeds <- replication_seeds(study_seed, n_pool)

  results <- map_replications(seq_len(n_replications), function(replication) {
    warnings_seen <- 0L
    out <- withCallingHandlers(
      run_replication(datasets[[replication]], ntree, seeds[replication, 1L]),
      warning = function(w) {
        warnings_seen <<- warnings_seen + 1L
        invokeRestart("muffleWarning")
      })
    out$replication <- replication
    out$n_warnings  <- warnings_seen
    out
  }, n_workers = n_workers, seed = study_seed)

  results <- do.call(rbind, results)
  results$n     <- n
  results$ntree <- ntree
  results$block <- paste0("BART outcome model, ", ntree, " trees")
  results$lower <- NA_real_
  results$upper <- NA_real_

  # Written as soon as the cell is done rather than at the end of the study: a
  # run that is interrupted then still leaves every cell that finished.
  utils::write.csv(results,
                   file.path(output_dir,
                             sprintf("cell-n%d-trees%d.csv", n, ntree)),
                   row.names = FALSE)
  message("  done in ", format(round(Sys.time() - started, 1)), ".")
  results
}

cells <- expand.grid(n = sample_sizes, ntree = tree_counts)
results <- do.call(rbind, Map(run_cell, cells$n, cells$ntree))

summary_table <- do.call(rbind, lapply(tree_counts, function(ntree) {
  performance_table(results[results$ntree == ntree, ], true_ate,
                    paste0("BART outcome model, ", ntree, " trees"),
                    method_order)
}))

cat("\nTable 2: BART outcome model, section 5.4.\n")
print_performance_table(summary_table, measures = c("ABias", "ESE", "RMSE"))


# ---------------------------------------------------------------------------
# The diagnostics that explain the table
# ---------------------------------------------------------------------------

# Coupling reduces the spread only if it moves the estimate against its error.
# Reporting the correlation between the two says whether it did, which the
# three summary measures on their own cannot.
wide <- stats::reshape(results[, c("n", "ntree", "replication", "method", "estimate")],
                       idvar = c("n", "ntree", "replication"),
                       timevar = "method", direction = "wide")
names(wide) <- sub("estimate.", "", names(wide), fixed = TRUE)

cat("\nWhat the coupling did, against the untilted g-formula:\n\n")
diagnostics <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
  part <- wide[wide$n == cells$n[i] & wide$ntree == cells$ntree[i], ]
  g <- part[["G-formula"]]
  do.call(rbind, lapply(method_order[-1L], function(m) {
    tilt <- part[[m]] - g
    flags <- results[results$n == cells$n[i] & results$ntree == cells$ntree[i] &
                       results$method == m, ]
    data.frame(n = cells$n[i], ntree = cells$ntree[i], method = m,
               signed_bias = mean(part[[m]]) - true_ate,
               sd_tilt = stats::sd(tilt), corr = stats::cor(g, tilt),
               ESE_change_pct = 100 * (stats::sd(part[[m]]) / stats::sd(g) - 1),
               no_op_pct = 100 * mean(flags$no_op),
               mean_lambda = mean(flags$lambda),
               stringsAsFactors = FALSE)
  }))
}))
print(diagnostics, row.names = FALSE, digits = 3)

cat("\nSigned bias of the untilted g-formula:\n")
print(do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
  part <- wide[wide$n == cells$n[i] & wide$ntree == cells$ntree[i], ]
  data.frame(n = cells$n[i], ntree = cells$ntree[i],
             signed_bias = mean(part[["G-formula"]]) - true_ate,
             mcse = stats::sd(part[["G-formula"]]) / sqrt(nrow(part)))
})), row.names = FALSE, digits = 3)

writeLines(c(
  "",
  "G-formula: the untilted posterior, the Bayesian g-formula over the BART",
  "  fitted means. Proposed: the same draws after coupling, under a correct",
  "  and under a misspecified propensity score model, under each of the three",
  "  readings of the Table 2 footnote.",
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
  "           Proposed (i.)  0.066  0.138  0.153",
  "",
  "The published run's BART configuration is not recorded in the paper, so",
  "these are not the same fit. See the header."))


# ---------------------------------------------------------------------------
# The record
# ---------------------------------------------------------------------------

utils::write.csv(summary_table, file.path(output_dir, "table2.csv"), row.names = FALSE)
utils::write.csv(diagnostics, file.path(output_dir, "table2-diagnostics.csv"), row.names = FALSE)
utils::write.csv(results, file.path(output_dir, "table2-replications.csv"), row.names = FALSE)

# Everything needed to say what produced the numbers, next to the numbers.
invisible(file.copy("scripts/simulation_table2.R",
                    file.path(output_dir, "table2-script.R"), overwrite = TRUE))
writeLines(c(
  paste("replications      ", n_replications),
  paste("sample sizes      ", paste(sample_sizes, collapse = ", ")),
  paste("trees             ", paste(tree_counts, collapse = ", ")),
  paste("BART chains       ", n_bart_chains),
  paste("BART burn-in      ", n_bart_burn, "per chain"),
  paste("BART draws        ", n_bart_draws, "per chain,", n_total, "in all,",
        n_used, "used"),
  paste("data set pool     ", n_pool),
  paste("BART k, power, base, numcut", bart_k, bart_power, bart_base, bart_numcut),
  paste("PS chains, mc, burn-in", ps_chains, ps_mc, ps_burn),
  paste("study seed        ", study_seed),
  paste("workers           ", n_workers),
  paste("elapsed           ", format(round(Sys.time() - started, 1))),
  "",
  utils::capture.output(utils::sessionInfo())),
  file.path(output_dir, "table2-settings.txt"))

cat("\nWrote", output_dir, "\n")
cat("Elapsed:", format(round(Sys.time() - started, 1)), "\n")
