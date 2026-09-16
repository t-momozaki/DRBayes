# The printed output is the whole of what most users ever see of a fit, so its
# formatting carries results rather than decoration. Three ways for it to hide
# one are easy to fall into and are pinned below: a credible interval rounded
# away on the paper's scale, an R-hat displayed as 1, and a sequential Monte
# Carlo sweep that failed without saying so.


# ---------------------------------------------------------------------------
# Fixtures. Building the object by hand rather than fitting one keeps these
# tests fast and lets the awkward states, a failed sweep or an R-hat that could
# not be computed, be constructed on demand. The last test in the file checks
# the fixtures against a real fit.
# ---------------------------------------------------------------------------

example_smc <- function(converged = TRUE) {
  if (converged) {
    list(lambda = 1.24, B_mean = 2.5e-04, B_mean_weighted = 2.5e-04,
         tol = 0.0047, n_steps = 14, max_steps = 200, converged = TRUE,
         ess = 1780, n_ancestors_min = 910, n_distinct = 1502)
  } else {
    list(lambda = 5, B_mean = 0.42, B_mean_weighted = 0.39,
         tol = 0.0047, n_steps = 200, max_steps = 200, converged = FALSE,
         ess = 12, n_ancestors_min = 3, n_distinct = 9)
  }
}

example_diagnostics <- function(rhat = c(1.003, 1.000),
                                ess_bulk = c(1200, 900),
                                ess_tail = c(1100, 950)) {
  data.frame(model = c("outcome", "propensity score"),
             parameter = c("AA", "(Intercept)"),
             mean = c(110.2, 0.11), sd = c(0.15, 0.2), rhat = rhat,
             ess_bulk = ess_bulk, ess_tail = ess_tail,
             mcse_mean = c(0.004, 0.006), stringsAsFactors = FALSE)
}

# The scale of Section 6 of the paper: an average treatment effect near 110
# carrying a posterior standard deviation near 0.15.
example_fit <- function(...) {
  set.seed(20)
  defaults <- list(
    pc      = stats::rnorm(4000, 110, 0.15),
    g.comp  = stats::rnorm(4000, 109.7, 0.15),
    family  = "gaussian",
    link    = "identity",
    ps.link = "logit",
    method  = "smc",
    smc     = example_smc(),
    diagnostics = example_diagnostics(),
    data_info = list(n_observations = 200, n_treated = 104, n_control = 96,
                     missing_observations = 0,
                     outcome_formula = YY ~ AA + WW1,
                     ps_formula = AA ~ WW1),
    call = quote(drbayes_pc(YY ~ AA + WW1, AA ~ WW1, data = trial)))
  structure(utils::modifyList(defaults, list(...)), class = "DRBayes")
}

example_draws <- function(n_iter = 300L, n_chain = 3L, shift = 0) {
  set.seed(31)
  chains <- paste0("chain", seq_len(n_chain))
  make <- function(parameters) {
    array(stats::rnorm(n_iter * n_chain * 2L), c(n_iter, n_chain, 2L),
          dimnames = list(NULL, chains, parameters))
  }
  outcome <- make(c("(Intercept)", "AA"))
  outcome[, n_chain, 2L] <- outcome[, n_chain, 2L] + shift
  list(outcome = outcome, ps = make(c("(Intercept)", "WW1")))
}

printed <- function(x, ...) {
  paste(utils::capture.output(print(x, ...)), collapse = "\n")
}

# The explanatory paragraphs are wrapped to the console width, so a phrase in
# one of them can be broken across a line. Matching against the squeezed text
# tests the wording without also testing where the wrap happened to fall.
squished <- function(x, ...) {
  gsub("[[:space:]]+", " ", printed(x, ...))
}

# The numbers of one row of the estimand table, as they were printed. Anchoring
# on a digit keeps the prose that mentions pc and g.comp out of the match.
estimand_strings <- function(text, name) {
  lines   <- strsplit(text, "\n", fixed = TRUE)[[1L]]
  pattern <- paste0("^", gsub(".", "[.]", name, fixed = TRUE), " +-?[0-9]")
  line    <- grep(pattern, lines, value = TRUE)
  expect_length(line, 1L)
  strsplit(trimws(substring(line, nchar(name) + 1L)), " +")[[1L]]
}

estimand_row <- function(text, name) {
  as.numeric(estimand_strings(text, name))
}

decimals_shown <- function(strings) {
  ifelse(grepl(".", strings, fixed = TRUE),
         nchar(sub("^[^.]*[.]", "", strings)), 0L)
}

# Drawing with no device open writes an Rplots.pdf into the working directory
# under Rscript, which R CMD check reports as a stray file. Everything here
# that draws goes to a device that discards its output.
with_null_device <- function(code) {
  grDevices::pdf(file = NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(code)
}


# ---------------------------------------------------------------------------
# print.DRBayes: the estimand table
# ---------------------------------------------------------------------------

test_that("the credible interval survives an effect on the paper's scale", {
  # Rounding to significant digits, which is what print.data.frame(digits = 3)
  # does, prints an effect of 110 with a standard deviation of 0.15 as three
  # copies of 110, and the interval vanishes.
  fit <- example_fit()
  out <- printed(fit)

  values <- estimand_row(out, "pc")
  expect_length(values, 3L)
  expect_lt(values[2L], values[1L])
  expect_lt(values[1L], values[3L])
  # Four posterior standard deviations of separation is the width to expect;
  # anything much narrower means the digits were spent on the mean.
  expect_gt(values[3L] - values[2L], 0.4)

  # The ends have to be the posterior quantiles rather than a rounded copy of
  # the mean, so they must agree with a direct computation from the draws to
  # within the rounding of the last digit printed.
  target <- c(mean(fit$pc),
              unname(stats::quantile(fit$pc, c(0.025, 0.975))))
  expect_lt(max(abs(values - target)), 1e-03)
})

test_that("both estimands are printed on one shared scale", {
  # A per-row choice of decimals would leave two estimands that cannot be read
  # against each other column by column, which is the only reason to print them
  # together.
  fit <- example_fit(g.comp = stats::rnorm(4000, 110, 15))
  out <- printed(fit)

  pc_shown <- estimand_strings(out, "pc")
  g_shown  <- estimand_strings(out, "g.comp")
  expect_identical(decimals_shown(pc_shown), decimals_shown(g_shown))

  # The wider of the two standard deviations sets the scale: three significant
  # digits of an sd near 15 is one decimal place.
  expect_identical(unique(decimals_shown(pc_shown)), 1L)
})

test_that("digits chooses how finely the standard deviation is resolved", {
  fit <- example_fit()
  decimals_at <- function(...) {
    unique(decimals_shown(estimand_strings(printed(fit, ...), "pc")))
  }
  expect_identical(decimals_at(), 3L)
  # Five significant digits of an sd near 0.15 is five decimal places.
  expect_identical(decimals_at(digits = 5), 5L)
  expect_error(print(fit, digits = 0),
               "digits must be a single number of at least 1")
  expect_error(print(fit, digits = c(3, 4)),
               "digits must be a single number of at least 1")
})

test_that("a degenerate posterior falls back to the requested digits", {
  # Every draw identical gives a standard deviation of zero, and the number of
  # decimal places cannot be read off its magnitude.
  fit <- example_fit(pc = rep(110, 100), g.comp = rep(110, 100))
  out <- printed(fit)
  expect_identical(unique(decimals_shown(estimand_strings(out, "pc"))), 3L)
  expect_identical(estimand_row(out, "pc"), rep(110, 3L))
})


# ---------------------------------------------------------------------------
# print.DRBayes: the sweep and the diagnostics
# ---------------------------------------------------------------------------

test_that("R-hat is printed at fixed decimals", {
  # Three significant digits of 1.003 is 1.00, which prints as 1 once the
  # trailing zeros are dropped, and a breach of the threshold becomes
  # invisible.
  fit <- example_fit(diagnostics = example_diagnostics(rhat = c(1.003, 1.000)))
  expect_match(printed(fit), "Largest R-hat: 1.003 (outcome: AA)",
               fixed = TRUE)
})

test_that("a sweep that missed the moment condition says so prominently", {
  fit <- example_fit(smc = example_smc(converged = FALSE))
  out <- squished(fit)

  expect_match(out, "*** THE MOMENT CONDITION WAS NOT MET ***", fixed = TRUE)
  expect_match(out, "NOT doubly robust", fixed = TRUE)
  # The numbers that make the failure actionable: where the sweep stopped, and
  # how far the constraint still is from the tolerance it needed to reach.
  expect_match(out, "lambda = 5", fixed = TRUE)
  expect_match(out, "200 of 200 steps", fixed = TRUE)
  expect_match(out, "|mean B_n| = 0.42", fixed = TRUE)
  expect_match(out, "above the tolerance 0.0047", fixed = TRUE)

  # summary() reports the same failure, and states it as an answer to a
  # question rather than only in prose.
  summary_out <- printed(summary(fit))
  expect_match(summary_out, "Moment condition met:     NO", fixed = TRUE)
  expect_match(summary_out, "*** THE MOMENT CONDITION WAS NOT MET ***",
               fixed = TRUE)

  # A sweep that did converge must not carry the warning, or it says nothing.
  converged <- squished(example_fit())
  expect_false(grepl("MOMENT CONDITION WAS NOT MET", converged, fixed = TRUE))
  expect_match(converged, "moment condition met", fixed = TRUE)
})

test_that("a fit with no recorded sweep is reported rather than skipped", {
  fit <- example_fit(smc = NULL)
  expect_match(printed(fit), "Sequential Monte Carlo: no information recorded",
               fixed = TRUE)
})

test_that("a diagnostic that could not be computed outranks a bad one", {
  # NA means the parameter was constant or held a non-finite draw, so it could
  # not be assessed at all. Ranking it below a finite 1.5 would let the
  # unassessable parameter pass unmentioned.
  fit <- example_fit(
    diagnostics = example_diagnostics(rhat = c(1.5, NA),
                                      ess_bulk = c(120, 900)))
  out <- printed(fit)

  expect_match(out, "Largest R-hat: NA (propensity score: (Intercept))",
               fixed = TRUE)
  expect_match(out, "could not be computed, treat as a failure", fixed = TRUE)
  # The effective sample size is the worst over bulk and tail of either model.
  expect_match(out, "Smallest ESS:  120 (outcome: AA)", fixed = TRUE)
  expect_match(out, "below 400", fixed = TRUE)
  expect_match(squished(fit), "do not meet the convergence criteria",
               fixed = TRUE)
})

test_that("a single model's diagnostics are labelled by parameter alone", {
  # convergence_diagnostics() on one model returns a table with no model
  # column. Naming the offending parameter is the whole use of the line, so it
  # has to survive the missing column rather than becoming "NA: AA".
  one_model <- example_diagnostics(rhat = c(1.5, 1.0))
  one_model$model <- NULL
  worst <- drbayes_worst_diagnostics(one_model)
  expect_identical(worst$rhat$where, "AA")
  expect_identical(worst$ess$where, "(Intercept)")
})

test_that("the convergence flags are made against the thresholds given", {
  # Whether a number is flagged is a comparison with a threshold, and writing
  # the threshold into the comparison instead of taking the one supplied would
  # make the printed verdict disagree with the warning the fit itself issued.
  # The same diagnostics are judged twice, against two different pairs.
  borderline <- example_diagnostics(rhat = c(1.03, 1.00),
                                    ess_bulk = c(150, 900),
                                    ess_tail = c(150, 950))
  worst <- drbayes_worst_diagnostics(borderline)
  judged <- function(rhat_max, ess_min) {
    thresholds <- list(rhat_max = rhat_max, ess_min = ess_min)
    paste(c(utils::capture.output(print_worst_diagnostics(worst, thresholds)),
            utils::capture.output(print_threshold_breaches(borderline,
                                                           thresholds))),
          collapse = " ")
  }

  strict <- gsub("[[:space:]]+", " ", judged(1.01, 400))
  expect_match(strict, "above 1.01", fixed = TRUE)
  expect_match(strict, "below 400", fixed = TRUE)
  expect_match(strict, "R-hat too high: AA", fixed = TRUE)
  expect_match(strict, "ESS too low: AA", fixed = TRUE)
  # A low effective sample size makes R-hat itself unreliable, and a reader who
  # is not told that reads the small R-hat of the other parameter as evidence.
  expect_match(strict, "makes R-hat itself unreliable", fixed = TRUE)

  relaxed <- gsub("[[:space:]]+", " ", judged(1.05, 100))
  expect_false(grepl("do not meet the convergence criteria", relaxed,
                     fixed = TRUE))
  expect_match(relaxed, "All parameters meet R-hat < 1.05 and ESS > 100",
               fixed = TRUE)
})

test_that("an ordinary fit is judged against the published thresholds", {
  # Vehtari et al. (2021) recommend R-hat below 1.01 and an effective sample
  # size above 400, and those are what drbayes_pc() itself warns against.
  expect_identical(drbayes_thresholds(example_fit()),
                   list(rhat_max = 1.01, ess_min = 400))
  strict <- example_fit(
    diagnostics = example_diagnostics(rhat = c(1.03, 1.00),
                                      ess_bulk = c(150, 900),
                                      ess_tail = c(150, 950)))
  expect_match(printed(strict), "above 1.01", fixed = TRUE)
  expect_match(squished(summary(strict)),
               "Thresholds: R-hat < 1.01, ESS > 400", fixed = TRUE)
})

test_that("diagnostics = none is reported as absent, not as passing", {
  fit <- example_fit(diagnostics = NULL)
  expect_match(printed(fit), "not computed (diagnostics = \"none\")",
               fixed = TRUE)
  expect_match(printed(summary(fit)), "Not computed (diagnostics = \"none\")",
               fixed = TRUE)
  expect_null(summary(fit)$diagnostics)
})

test_that("print describes the models, the data and the draws", {
  out <- printed(example_fit())
  expect_match(out, "gaussian family, identity link", fixed = TRUE)
  expect_match(out, "Propensity score: logit link", fixed = TRUE)
  expect_match(out, "Observations:     200 (104 treated, 96 control)",
               fixed = TRUE)
  expect_match(out, "Posterior draws:  4000", fixed = TRUE)
  expect_match(out, "drbayes_pc(YY ~ AA + WW1, AA ~ WW1, data = trial)",
               fixed = TRUE)

  # Rows lost to missing data change what the estimate is about, so they are
  # named rather than left to be inferred from the total.
  dropped <- example_fit(
    data_info = utils::modifyList(example_fit()$data_info,
                                  list(missing_observations = 7L)))
  expect_match(printed(dropped), "7 dropped for missing values", fixed = TRUE)
})

test_that("a fit that predates the recorded metadata still prints", {
  fit <- example_fit(family = NULL, link = NULL, ps.link = NULL,
                     data_info = NULL, call = NULL)
  out <- printed(fit)
  expect_match(out, "unknown family, unknown link", fixed = TRUE)
  expect_match(out, "Propensity score: unknown link", fixed = TRUE)
  expect_match(out, "Observations:     unknown", fixed = TRUE)
  expect_false(grepl("Call:", out, fixed = TRUE))
})

test_that("print returns its argument invisibly", {
  fit <- example_fit()
  utils::capture.output(returned <- expect_invisible(print(fit)))
  expect_identical(returned, fit)
  utils::capture.output(
    returned <- expect_invisible(print(summary(fit))))
  expect_identical(returned, summary(fit))
})


# ---------------------------------------------------------------------------
# summary.DRBayes
# ---------------------------------------------------------------------------

test_that("summary reports the posterior summaries the draws actually have", {
  # pc and g.comp are both average treatment effects of similar magnitude, so a
  # swapped row order would look entirely plausible. Separating their means
  # makes the row labels testable.
  fit <- example_fit(pc = stats::rnorm(3000, 110, 0.15),
                     g.comp = stats::rnorm(3000, 90, 0.20))
  s <- summary(fit)

  expect_s3_class(s, "summary.DRBayes")
  expect_identical(rownames(s$estimands), c("pc", "g.comp"))
  expect_named(s$estimands, c("mean", "sd", "2.5%", "50%", "97.5%"))
  expect_identical(s$n_draws, length(fit$pc))

  for (name in c("pc", "g.comp")) {
    draws <- fit[[name]]
    expect_equal(s$estimands[name, "mean"], mean(draws))
    expect_equal(s$estimands[name, "sd"], stats::sd(draws))
    expect_equal(unlist(s$estimands[name, c("2.5%", "50%", "97.5%")],
                        use.names = FALSE),
                 unname(stats::quantile(draws, c(0.025, 0.5, 0.975))))
  }
})

test_that("summary carries the fit's own description forward unchanged", {
  fit <- example_fit()
  s <- summary(fit)
  expect_identical(s$call, fit$call)
  expect_identical(s$family, fit$family)
  expect_identical(s$link, fit$link)
  expect_identical(s$ps.link, fit$ps.link)
  expect_identical(s$smc, fit$smc)
  expect_identical(s$diagnostics, fit$diagnostics)
  expect_identical(s$thresholds, list(rhat_max = 1.01, ess_min = 400))

  out <- printed(s)
  expect_match(out, "Outcome formula:  YY ~ AA + WW1", fixed = TRUE)
  expect_match(out, "PS formula:       AA ~ WW1", fixed = TRUE)
  # The full table, which is what summary() adds over print().
  expect_match(out, "ess_bulk", fixed = TRUE)
  expect_match(out, "Tilting parameter lambda: 1.24", fixed = TRUE)
  expect_match(out, "Steps taken:              14 of 200", fixed = TRUE)
})

test_that("draws that are missing are left out of the summaries", {
  # Otherwise a single non-finite draw turns every number in the table into NA,
  # which reads as a broken fit rather than as a damaged draw.
  fit <- example_fit()
  fit$pc[c(1L, 2L)] <- NA
  s <- summary(fit)
  expect_true(is.finite(s$estimands["pc", "mean"]))
  expect_equal(s$estimands["pc", "mean"], mean(fit$pc, na.rm = TRUE))
})


# ---------------------------------------------------------------------------
# Objects the methods must refuse
# ---------------------------------------------------------------------------

test_that("an object that is not a fit is refused by name", {
  # A message naming the wrong cause sends the reader to the wrong place, so
  # each failure has to say which component is at fault and why.
  expect_error(print(structure(list(g.comp = 1:3), class = "DRBayes")),
               "no 'pc' component")
  expect_error(print(structure(list(pc = 1:3), class = "DRBayes")),
               "no 'g.comp' component")
  expect_error(
    print(structure(list(pc = "a", g.comp = 1:3), class = "DRBayes")),
    "must be a non-empty numeric vector of posterior draws, but is character")
  expect_error(
    print(structure(list(pc = numeric(0), g.comp = 1:3), class = "DRBayes")),
    "of length 0")
  expect_error(
    print(structure(list(pc = 1:3, g.comp = 1:4), class = "DRBayes")),
    "hold 3 and 4 draws")
  expect_error(summary(structure(list(g.comp = 1:3), class = "DRBayes")),
               "no 'pc' component")
})


# ---------------------------------------------------------------------------
# plot.DRBayes, density branch
# ---------------------------------------------------------------------------

test_that("the density plot leaves the graphics parameters as it found them", {
  # CRAN policy, and a plot that silently keeps its own margins would reshape
  # every figure the user draws afterwards.
  fit <- example_fit()
  with_null_device({
    before <- graphics::par(no.readonly = TRUE)
    returned <- expect_invisible(plot(fit))
    after <- graphics::par(no.readonly = TRUE)
    expect_identical(after, before)
    expect_identical(returned, fit)
  })
})

test_that("the density plot refuses arguments it cannot honour", {
  fit <- example_fit()
  with_null_device({
    expect_error(plot(fit, col = "red"),
                 "col must give two colours, one for pc and one for g.comp")
    expect_error(plot(fit, legend = "yes"), "legend must be TRUE or FALSE")
    expect_error(plot(fit, legend = NA), "legend must be TRUE or FALSE")
    expect_error(plot(fit, type = "trace"), "'arg' should be one of")
  })
})

test_that("non-finite draws are reported rather than passed to density()", {
  # stats::density() fails on a missing value with a message about an argument
  # the caller never supplied, which sends the reader hunting in the wrong
  # place. How many draws were dropped is the fact that matters.
  fit <- example_fit()
  fit$pc[1:5] <- NA
  fit$pc[6L] <- Inf
  with_null_device(
    expect_warning(plot(fit),
                   "6 of 4000 draws of pc are missing or infinite"))

  all_missing <- example_fit(pc = rep(NA_real_, 100),
                             g.comp = stats::rnorm(100))
  with_null_device(
    expect_error(plot(all_missing),
                 "All 100 draws of pc are missing or infinite"))
})


# ---------------------------------------------------------------------------
# plot.DRBayes, rank branch, and the rank_plot(, plot = FALSE) computation behind it
# ---------------------------------------------------------------------------

test_that("plot(type = 'rank') reaches rank_plot and pools the ranks", {
  # Every bin holds the same number of pooled draws by construction. That is
  # what makes a departure from flatness in one panel a statement about that
  # chain rather than about the binning, and so it is the property that makes
  # the plot mean anything at all.
  fit <- example_fit(draws = example_draws())
  result <- plot(fit, type = "rank", parameter = "AA", plot = FALSE)

  expect_named(result, c("counts", "breaks", "expected", "parameter"))
  expect_identical(result$parameter, "AA")
  expect_identical(dim(result$counts), c(20L, 3L))
  expect_identical(colnames(result$counts), paste0("chain", 1:3))
  expect_true(all(rowSums(result$counts) == 300L * 3L / 20L))
  expect_true(all(colSums(result$counts) == 300L))
  expect_identical(result$expected, 300 / 20)
  expect_identical(result$breaks, seq(0, 900, length.out = 21L))
})

test_that("a chain in the wrong place leans in its rank panel", {
  # The failure the plot exists to show. Pooled uniformity is untouched by the
  # shift, so only the split between chains can reveal it.
  fit <- example_fit(draws = example_draws(shift = 3))
  counts <- plot(fit, type = "rank", parameter = "AA",
                 plot = FALSE)$counts

  expect_true(all(rowSums(counts) == 45L))
  # The shifted chain holds the high ranks and almost none of the low ones.
  expect_gt(counts[20L, 3L], counts[1L, 3L])
  expect_lt(counts[1L, 3L], 45 / 3)
  # The chains left behind take the low ranks the shifted one vacated.
  expect_gt(counts[1L, 1L], counts[20L, 1L])
  expect_gt(counts[1L, 2L], counts[20L, 2L])
})

test_that("plot = FALSE computes the counts without opening a device", {
  # Both because a rank plot is worth testing where no device exists, and
  # because opening one under Rscript leaves an Rplots.pdf behind.
  fit <- example_fit(draws = example_draws())
  before <- grDevices::dev.cur()
  result <- plot(fit, type = "rank", plot = FALSE)
  expect_identical(grDevices::dev.cur(), before)
  expect_true(is.matrix(result$counts))
})

test_that("the model and parameter asked for are the ones plotted", {
  # An ignored model argument would plot the outcome coefficients under the
  # propensity score's name, which is exactly the kind of silently wrong answer
  # a diagnostic must not give. The two models carry disjoint parameter names,
  # so the name that comes back identifies which array was read.
  fit <- example_fit(draws = example_draws())

  outcome <- plot(fit, type = "rank", model = "outcome", parameter = "AA",
                  plot = FALSE)
  ps <- plot(fit, type = "rank", model = "ps", parameter = "WW1",
             plot = FALSE)
  expect_identical(outcome$parameter, "AA")
  expect_identical(ps$parameter, "WW1")
  expect_false(identical(outcome$counts, ps$counts))

  # A name and the index of that name select the same column.
  by_index <- plot(fit, type = "rank", model = "ps", parameter = 2L,
                   plot = FALSE)
  expect_identical(by_index, ps)

  expect_error(plot(fit, type = "rank", parameter = "WW1", plot = FALSE),
               "There is no parameter called 'WW1' in draws")
  expect_error(plot(fit, type = "rank", parameter = 5L, plot = FALSE),
               "parameter index 5 is outside draws")
})

test_that("a fit without per-chain draws says how to get a rank plot", {
  # The draws of the treatment effect that a fit does carry are not chains, and
  # ranking them would produce a plot that looks like a diagnostic and is not
  # one. Refusing is only useful if the message names the way forward.
  fit <- example_fit()
  expect_error(plot(fit, type = "rank", plot = FALSE),
               "Rank plots need the per-chain draws of the model coefficients")
  expect_error(plot(fit, type = "rank", plot = FALSE),
               "rank_plot(", fixed = TRUE)

  # A fit carrying the other model's draws only must not fall back to them.
  half <- example_fit(draws = example_draws()["ps"])
  expect_error(plot(half, type = "rank", model = "outcome", plot = FALSE),
               "does not carry them")
  expect_silent(plot(half, type = "rank", model = "ps", plot = FALSE))
})

test_that("drawing the rank panels restores the graphics parameters", {
  # The panels are laid out with mfrow, which would otherwise leave every later
  # plot of the session in a grid.
  fit <- example_fit(draws = example_draws())
  with_null_device({
    before <- graphics::par(no.readonly = TRUE)
    result <- plot(fit, type = "rank", plot = FALSE)
    after <- graphics::par(no.readonly = TRUE)
    expect_identical(after$mfrow, before$mfrow)
    expect_identical(after$mar, before$mar)
    expect_identical(after$oma, before$oma)
    # Drawing must not change what is computed.
    expect_identical(result$counts,
                     plot(fit, type = "rank", plot = FALSE)$counts)
  })
})


# ---------------------------------------------------------------------------
# rank_plot(, plot = FALSE) called directly, which is the route the error above points at
# ---------------------------------------------------------------------------

test_that("a matrix of draws is read as iterations by chains", {
  # It is the shape reached for when there is one quantity of interest. Reading
  # it the other way round, as one chain of several parameters, would draw a
  # single panel that cannot disagree with anything and so always looks well
  # mixed.
  set.seed(41)
  theta <- matrix(stats::rnorm(200 * 4), 200, 4)
  theta[, 4L] <- theta[, 4L] + 3
  result <- rank_plot(theta, bins = 10L, plot = FALSE)

  expect_identical(dim(result$counts), c(10L, 4L))
  expect_identical(result$parameter, "parameter")
  expect_true(all(rowSums(result$counts) == 200L * 4L / 10L))
  # The shifted column is a chain, so the shift shows.
  expect_gt(result$counts[10L, 4L], result$counts[1L, 4L])
})

test_that("tied draws are ranked by average, not by position in the array", {
  # Ranking ties in array order would order chains by their position: chain 1
  # would take the lowest ranks of every tie group and its panel would lean,
  # reporting a disagreement between chains that the draws do not contain.
  # Draws that are all equal make the difference exact.
  tied <- array(1, dim = c(200L, 4L, 1L),
                dimnames = list(NULL, paste0("chain", 1:4), "theta"))
  result <- rank_plot(tied, bins = 20L, plot = FALSE)

  expect_true(all(result$counts[, 1L] == result$counts[, 2L]))
  expect_true(all(result$counts[, 1L] == result$counts[, 4L]))
  # Every draw shares the middle rank, so one bin holds them all.
  expect_identical(sum(result$counts > 0L), 4L)
})

test_that("rank_plot refuses draws and settings it cannot rank", {
  set.seed(42)
  draws <- array(stats::rnorm(100 * 2), c(100L, 2L, 1L),
                 dimnames = list(NULL, NULL, "theta"))

  # A non-finite draw sorts to one end of the rank axis, which looks exactly
  # like a chain stuck in the tail. Refusing names how many there are.
  spoiled <- draws
  spoiled[1:3, 1L, 1L] <- c(NA, NaN, Inf)
  expect_error(rank_plot(spoiled, plot = FALSE),
               "Parameter 'theta' contains 3 missing or infinite draw")

  # More bins than iterations leaves most bins empty in every panel, so
  # flatness stops meaning anything.
  expect_error(rank_plot(draws, bins = 200L, plot = FALSE),
               "but each chain has only 100 iterations")

  expect_error(rank_plot(draws, bins = 1L, plot = FALSE),
               "bins must be a single whole number of at least 2")
  expect_error(rank_plot(draws, bins = 10.5, plot = FALSE),
               "bins must be a single whole number of at least 2")
  expect_error(rank_plot(draws, plot = NA), "plot must be TRUE or FALSE")
  expect_error(rank_plot("not draws", plot = FALSE), "draws must be a numeric array")
  expect_error(rank_plot(stats::rnorm(100), plot = FALSE), "draws must have 3 dimensions")
  expect_error(rank_plot(array(1, c(1L, 2L, 1L))),
               "at least 2 iterations per chain")
  expect_error(rank_plot(draws, parameter = c(1L, 2L), plot = FALSE),
               "parameter must be a single name or index")
  expect_error(rank_plot(draws, parameter = 1.5, plot = FALSE),
               "parameter must be a parameter name or a whole number index")
  expect_error(rank_plot(array(numeric(0), c(2L, 0L, 1L))),
               "at least one chain and one parameter")
})

test_that("unnamed parameters are given a name to be plotted under", {
  # The name becomes the title of the figure, so leaving it NULL would produce
  # a rank plot of something the reader cannot identify.
  set.seed(43)
  draws <- array(stats::rnorm(100 * 2 * 2), c(100L, 2L, 2L))
  result <- rank_plot(draws, parameter = 2L, bins = 10L, plot = FALSE)
  expect_identical(result$parameter, "parameter 2")
  expect_identical(colnames(result$counts), c("chain 1", "chain 2"))
})


# ---------------------------------------------------------------------------
# The fixtures above stand in for a fit. This is the check that they still do.
# ---------------------------------------------------------------------------

test_that("the methods read an object drbayes_pc actually produced", {
  skip_on_cran()

  set.seed(1)
  data <- generate_dataset(nn = 150, pp = 0, seed = 1)
  utils::capture.output(
    fit <- suppressWarnings(
      drbayes_pc(YY ~ AA + WW1 + WW2 + WW3 + WW4,
                 AA ~ WW1 + WW2 + WW3 + WW4, data = data,
                 outcome.model = bayes_lm, ps.model = bayes_logit,
                 mc = 1500, bn = 500, chains = 2L, verbose = FALSE,
                 control = drbayes_control(n_steps = 50))),
    file = nullfile())

  out <- printed(fit)
  s <- summary(fit)

  # The table print() shows is the same one summary() returns, so a reader who
  # takes the numbers from either gets the same answer.
  expect_lt(max(abs(estimand_row(out, "pc") -
                      unlist(s$estimands["pc", c("mean", "2.5%", "97.5%")],
                             use.names = FALSE))),
            1e-03)

  # The headline diagnostics are the worst of the full table, not a separate
  # computation that could drift from it.
  worst_rhat <- max(fit$diagnostics$rhat)
  expect_match(out,
               paste0("Largest R-hat: ",
                      formatC(worst_rhat, format = "f", digits = 3)),
               fixed = TRUE)
  worst_ess <- min(fit$diagnostics$ess_bulk, fit$diagnostics$ess_tail)
  expect_match(out, paste0("Smallest ESS:  ", round(worst_ess)), fixed = TRUE)

  expect_identical(s$n_draws, length(fit$pc))
  expect_match(out, "Posterior draws:  ", fixed = TRUE)
  expect_match(out, paste0("Observations:     ",
                           fit$data_info$n_observations), fixed = TRUE)
})


test_that("the rank plot can actually be reached from a fit", {
  # A rank histogram is a statement about how chains interleave, so it needs
  # the per-chain draws rather than the tilted particles, which are resampled
  # and have no chains left. Three things have to line up: the control flag
  # that asks for them, the fit that stores them, and a refusal naming an
  # argument the reader can actually type when they were not asked for.
  skip_on_cran()
  set.seed(11)
  n <- 150L
  d <- data.frame(X1 = stats::rnorm(n), X2 = stats::rnorm(n))
  d$A <- stats::rbinom(n, 1, 0.5)
  d$Y <- stats::rnorm(n)
  fit_of <- function(keep) {
    suppressWarnings(suppressMessages(drbayes_pc(
      Y ~ A + X1, A ~ X1 + X2, d,
      outcome.model = bayes_lm, ps.model = bayes_logit,
      control = drbayes_control(mc = 800, bn = 200, chains = 4L,
                                keep_draws = keep),
      diagnostics = "none")))
  }

  kept <- fit_of(TRUE)
  expect_named(kept$draws, c("outcome", "ps"))
  expect_length(dim(kept$draws$outcome), 3L)

  with_null_device({
    result <- plot(kept, type = "rank", parameter = "A")
    expect_named(result, c("counts", "breaks", "expected", "parameter"))
    expect_equal(sum(result$counts), prod(dim(kept$draws$outcome)[1:2]))
    expect_false(is.null(plot(kept, type = "rank", model = "ps",
                             parameter = 1L)$counts))
  })

  # And when they were not kept, the refusal has to name something the reader
  # can actually type.
  bare <- fit_of(FALSE)
  expect_null(bare$draws)
  with_null_device(
    expect_error(plot(bare, type = "rank"), "keep_draws = TRUE"))
})


test_that("starting values reach the samplers", {
  # The same control object has to mean the same thing to drbayes_pc() and to
  # drbayes_select(), both of which build sampler arguments of their own. A
  # starting value that reaches one and not the other is the failure this
  # guards against.
  skip_on_cran()
  set.seed(12)
  n <- 150L
  d <- data.frame(X1 = stats::rnorm(n), X2 = stats::rnorm(n))
  d$A <- stats::rbinom(n, 1, 0.5)
  d$Y <- stats::rnorm(n)
  run <- function(init) {
    suppressWarnings(suppressMessages(drbayes_pc(
      Y ~ A + X1, A ~ X1 + X2, d,
      outcome.model = bayes_lm, ps.model = bayes_logit,
      control = drbayes_control(mc = 600, bn = 200, chains = 2L, init = init),
      diagnostics = "none")))
  }
  expect_false(identical(run(NULL)$pc,
                         run(list(c(0, 0, 0), c(50, 50, 50)))$pc))

  # The length check the documentation promises happens when the sampler runs,
  # which it never did.
  expect_error(run(list(c(0, 0), c(1, 1))))
})


test_that("a fit says which average its interval is over", {
  # Two posteriors, three estimands in the package, and nothing named any of
  # them. A reader seeing mean(fit$pc) and a credible interval had no way to
  # tell whether it was for the covariates observed or for the population they
  # came from, and the two differ by 11 to 21 percent in width wherever the
  # fitted contrast varies from unit to unit.
  # The text is wrapped to the console width, so a phrase to match has to be
  # one that cannot be split across a line.
  fit <- example_fit()
  expect_output(print(fit), "average causal effect")
  expect_output(print(fit), "drbayes_pc")
  expect_output(print(summary(fit)), "observed")
})
