# Reproduces section 7.3 of Orihara, Momozaki and Sugasawa (2025), "Impact of
# Right Heart Catheterization on Mortality with Confounder Selection": the
# effect of a right heart catheter on death within 30 days of admission,
# estimated three ways, with confounders chosen from fifty candidates.
#
# The three estimates the section reports, all risk differences with 95 percent
# credible intervals:
#
#   G-formula, no selection           0.229 ( 0.032, 0.395)  length 0.363
#   Posterior coupling, no selection  0.043 (-0.075, 0.205)  length 0.280
#   Posterior coupling, Algorithm 4   0.188 ( 0.122, 0.312)  length 0.190
#
# The reading the section draws from them is about the interval, not the point:
# selecting confounders and then coupling gives the shortest of the three,
# 0.68 of the width of coupling alone and 0.52 of the g-formula's.
#
# Run it with
#   Rscript scripts/application_rhc.R
# from the top of the package. It calls data-raw/rhc.R the first time, which
# downloads about 2 MB from https://hbiostat.org/data/ into a cache directory
# and never into the repository. Set DRBAYES_DATA_DIR to keep that cache; the
# default is the session's temporary directory. As written the script takes
# about two minutes on one core, 1.8 and 1.9 minutes on two runs of the machine
# it was written on, plus the download the first time.
#
# The settings below are reduced so that it finishes in that time. They are not
# the paper's, and the script says so rather than reading its own numbers
# against the published ones: the settings the paper used sit beside each one,
# commented out, and the comparison is printed only when they are restored and
# the fits meet their convergence criteria.
#
# References
#   Connors, A. F., et al. (1996). The effectiveness of right heart
#     catheterization in the initial care of critically ill patients. JAMA,
#     276(11), 889-897.
#   Ning, Y., Sida, P., & Imai, K. (2020). Robust estimation of causal effects
#     via a high-dimensional covariate balancing propensity score. Biometrika,
#     107(3), 533-554.
#   Orihara, S., Momozaki, T., & Sugasawa, S. (2025). Bayesian Doubly Robust
#     Causal Inference via Posterior Coupling. arXiv:2506.04868.

library(DRBayes)


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

# Both the Gibbs sweeps and the sequential Monte Carlo cost time linear in the
# number of patients, so analysing a random part of the cohort is the setting
# that buys the most. NULL analyses all 5735, which is what section 7.3 does.
n_patients <- 1500L
# n_patients <- NULL

# The footnote to Table 1 asks for 20,000 posterior draws, and section 7 gives
# no reason to depart from it. Sixty-five horseshoe coefficients on each side
# need them; the reduced setting below will warn that the chains have not
# converged, and it is right to.
mcmc <- list(mc = 1000L, bn = 200L, thin = 1L, chains = 2L)
paper_mcmc <- list(mc = 6000L, bn = 1000L, thin = 1L, chains = 4L)
# mcmc <- paper_mcmc

# Horseshoe posteriors are eroded by repeated kernel smoothing, so the sweep
# takes fewer and larger steps in the tilting parameter. See ?drbayes_control.
control <- drbayes_control(n_steps = 200)

# Section 7.3: "we define the subset of selected covariates as
# S = {j : |alpha_bar_j| >= 0.01}".
threshold <- 0.01

analysis_seed <- 20250604L

data_dir <- Sys.getenv("DRBAYES_DATA_DIR", unset = "")
if (!nzchar(data_dir)) {
  data_dir <- file.path(tempdir(), "drbayes-data")
}


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

rhc_file    <- file.path(data_dir, "rhc.rds")
prepare_rhc <- "data-raw/rhc.R"

if (!file.exists(rhc_file)) {
  if (!file.exists(prepare_rhc)) {
    stop(prepare_rhc, " is not here, so the data cannot be prepared. Run ",
         "this script from the top of the package, as in ",
         "Rscript scripts/application_rhc.R.", call. = FALSE)
  }
  message("Preparing the data with ", prepare_rhc, ".")
  # In its own environment, so that the preparation script's helpers do not
  # end up beside the analysis and shadow something later.
  source(prepare_rhc, local = new.env())
}
if (!file.exists(rhc_file)) {
  stop(prepare_rhc, " did not produce ", rhc_file, ". Run it on its own to ",
       "see why, then rerun this script.", call. = FALSE)
}

rhc <- readRDS(rhc_file)
cohort_size <- nrow(rhc)

if (!is.null(n_patients) && n_patients < nrow(rhc)) {
  set.seed(analysis_seed)
  rhc <- rhc[sample.int(nrow(rhc), n_patients), , drop = FALSE]
  # Rare categories such as the seven colon cancer admissions can vanish from
  # a subsample. A factor level with no observations left becomes a column of
  # zeros in the design matrix, which is not estimable and stops the sampler.
  rhc <- droplevels(rhc)
  message("Reduced to a random ", n_patients, " of the ", cohort_size,
          " patients section 7.3 analyses.")
}

# rhc_prepare() puts the treatment and the outcome first and nothing but
# candidate confounders after them, so the candidates need no separate list.
candidates <- setdiff(names(rhc), c("rhc", "death30"))

outcome.formula <- stats::as.formula(
  paste("death30 ~ rhc +", paste(candidates, collapse = " + ")))
ps.formula <- stats::as.formula(
  paste("rhc ~", paste(candidates, collapse = " + ")))

message("Analysing ", nrow(rhc), " patients and ", length(candidates),
        " candidate confounders.")


# ---------------------------------------------------------------------------
# Estimation
# ---------------------------------------------------------------------------

#' Posterior mean and 95 percent credible interval, as section 7.3 reports them
#'
#' @param draws Numeric vector of posterior draws of the risk difference.
#' @param method Label for the row.
#' @param moment One of "met", "NOT MET" or "not applicable", carried into the
#'   printed table so that a row resting on a sweep that failed is marked where
#'   it is read rather than further down the output.
#'
#' @return A one row data frame.
summarise_draws <- function(draws, method, moment) {
  bounds <- stats::quantile(draws, c(0.025, 0.975), names = FALSE)
  data.frame(method = method,
             mean   = mean(draws),
             lower  = bounds[1],
             upper  = bounds[2],
             length = bounds[2] - bounds[1],
             moment = moment,
             stringsAsFactors = FALSE)
}


#' Whether the Step 1 draws met the thresholds their control object sets
#'
#' Posterior coupling tilts draws from both models, so a fit is only as good as
#' the worse of the two posteriors it starts from.
#'
#' @param fit A fit from \code{drbayes_pc()} or \code{drbayes_select()}.
#'
#' @return TRUE when every parameter of both models met the R-hat and
#'   effective sample size thresholds of \code{fit$control}.
step1_converged <- function(fit) {
  diagnostics <- fit$diagnostics
  if (is.null(diagnostics)) {
    return(FALSE)
  }
  all(diagnostics$rhat <= fit$control$rhat_max,
      diagnostics$ess_bulk >= fit$control$ess_min,
      diagnostics$ess_tail >= fit$control$ess_min)
}


#' How a fit describes its own trustworthiness, in one line
#'
#' @param fit A fit from \code{drbayes_pc()} or \code{drbayes_select()}.
#' @param label Name of the fit.
#'
#' @return A string.
describe_fit <- function(fit, label) {
  diagnostics <- fit$diagnostics
  sprintf("  %-13s R-hat %5.3f  ESS %5.0f  moment %s at lambda %s",
          label,
          max(diagnostics$rhat),
          min(diagnostics$ess_bulk, diagnostics$ess_tail),
          if (isTRUE(fit$smc$converged)) "met    " else "NOT MET",
          format(fit$smc$lambda, digits = 3))
}

started <- Sys.time()

# Without selection. Horseshoe priors on both sides, but nothing is dropped and
# the moment condition ranges over every covariate, which is the "simple
# procedure" section 7.2 warns about: a weak confounder is shrunk away in both
# models at once and its confounding is absorbed by the treatment effect.
message("Fitting without confounder selection.")
plain <- do.call(drbayes_pc, c(
  list(outcome.formula = outcome.formula, ps.formula = ps.formula,
       data = rhc, outcome.prior = "horseshoe", ps.prior = "horseshoe",
       control = control, seed = analysis_seed, verbose = FALSE),
  mcmc))

# Algorithm 4. Step 1 fits both models under horseshoe priors and selects the
# covariates whose propensity score coefficient survives; step 2 tilts only
# that block, which gives a covariate discarded from the outcome model a second
# chance through its propensity score coefficient.
message("Fitting Algorithm 4, with confounder selection.")
selected <- do.call(drbayes_select, c(
  list(outcome.formula = outcome.formula, ps.formula = ps.formula,
       data = rhc, threshold = threshold, control = control,
       seed = analysis_seed, verbose = FALSE),
  mcmc))

#' The moment condition cell of the printed table
#'
#' @param fit A fit from \code{drbayes_pc()} or \code{drbayes_select()}.
#'
#' @return "met" when the sweep satisfied the moment condition, "NOT MET"
#'   otherwise. A fit whose sweep never met it is not doubly robust, which is
#'   why this is reported in the table rather than below it.
moment_state <- function(fit) {
  if (isTRUE(fit$smc$converged)) "met" else "NOT MET"
}

results <- rbind(
  summarise_draws(plain$g.comp, "G-formula, no selection", "not applicable"),
  summarise_draws(plain$pc, "Posterior coupling, no selection",
                  moment_state(plain)),
  summarise_draws(selected$pc, "Posterior coupling, Algorithm 4",
                  moment_state(selected)))


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

cat("\nSection 7.3: effect of right heart catheterization on 30-day",
    "mortality\n")
cat("Risk difference, posterior mean and 95 percent credible interval\n\n")

shown <- data.frame(
  Method = results$method,
  Mean   = formatC(results$mean, format = "f", digits = 3),
  `95% CI` = sprintf("(%.3f, %.3f)", results$lower, results$upper),
  Length = formatC(results$length, format = "f", digits = 3),
  `Moment condition` = results$moment,
  check.names = FALSE, stringsAsFactors = FALSE)
print(shown, row.names = FALSE)

writeLines(c(
  "",
  "A tilted posterior whose sweep never met the moment condition is not",
  "doubly robust, so the column is part of the answer rather than a",
  "diagnostic to look up afterwards. The g-formula row is untilted and has no",
  "moment condition to meet."))

fits <- list(`no selection` = plain, `Algorithm 4` = selected)
cat("\nStep 1 convergence and the sweep:\n")
for (label in names(fits)) {
  cat(describe_fit(fits[[label]], label), "\n", sep = "")
}
cat("  thresholds: R-hat <= ", control$rhat_max, ", ESS >= ", control$ess_min,
    "\n", sep = "")

cat("\nSelected ", length(selected$selected), " of ", selected$n_candidates,
    " design matrix columns at |alpha_bar| >= ", threshold, ",",
    "\n  where the paper selects 44 of 58.\n", sep = "")

cat("\nElapsed:", format(round(Sys.time() - started, 1)), "\n")

writeLines(c(
  "",
  "The candidate list in data-raw/rhc.R is the standard Connors et al.",
  "(1996) set and gives 50 candidates in 65 columns, where section 7.3",
  "reports 48 in 58 following Harada and Taguri (2025). The paper does not",
  "restate their list, so neither the selected counts nor the estimates can",
  "be expected to match exactly."))

# Reading these numbers against the published ones takes a run that could
# support the comparison. Say why not, rather than drawing a conclusion from a
# configuration the same output has just reported as untrustworthy.
reservations <- character(0)
if (!is.null(n_patients) && n_patients < cohort_size) {
  reservations <- c(reservations, sprintf(
    "  - it analyses %d of the %d patients section 7.3 analyses",
    n_patients, cohort_size))
}
if (!identical(mcmc, paper_mcmc)) {
  reservations <- c(reservations, sprintf(
    "  - it draws %d posterior samples rather than the 20,000 of the paper",
    mcmc$chains * length(seq(mcmc$bn + 1L, mcmc$mc, by = mcmc$thin))))
}
for (label in names(fits)) {
  if (!step1_converged(fits[[label]])) {
    reservations <- c(reservations, paste0(
      "  - the Step 1 chains of the ", label, " fit did not converge"))
  }
  if (!isTRUE(fits[[label]]$smc$converged)) {
    reservations <- c(reservations, paste0(
      "  - the sweep of the ", label,
      " fit did not meet the moment condition"))
  }
}

if (length(reservations) > 0) {
  writeLines(c(
    "",
    "The estimates above are NOT compared with the published ones, because",
    "this run cannot support the comparison:",
    reservations,
    "",
    "Restore the settings marked in the block at the top and rerun. At the",
    "paper's settings this takes hours rather than minutes, which is what the",
    "reduced ones are for."))
} else {
  crude <- mean(rhc$death30[rhc$rhc == 1]) - mean(rhc$death30[rhc$rhc == 0])
  writeLines(c(
    "",
    "This run is at the paper's settings and both fits met their convergence",
    "criteria, so the numbers above can be read against the published ones.",
    sprintf("For reference, the crude risk difference in this cohort is %.3f.",
            crude),
    "The reading of section 7.3 is about the interval: selection followed by",
    "coupling should give the shortest of the three. Published, it is 0.68 of",
    "the width of coupling alone and 0.52 of the g-formula's.",
    sprintf("Here it is %.3f against %.3f for coupling alone and %.3f for the",
            results$length[3], results$length[2], results$length[1]),
    "g-formula."))
}
