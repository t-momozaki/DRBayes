#' Methods for posterior coupling fits
#'
#' @description
#' Printing, summarising and plotting for the objects of class `"DRBayes"`
#' returned by [drbayes_pc()].
#'
#' @details
#' Every method reports two posteriors for the average treatment effect, and it
#' matters which one is which.
#'
#' `g.comp` is the untilted Bayesian g-formula posterior. It is what the outcome
#' model alone implies about the average treatment effect, obtained by averaging
#' the fitted contrast over the observed covariates. It uses the propensity
#' score model not at all, so it is consistent only if the outcome model is
#' correctly specified. In the notation of the paper it is the tilted posterior
#' at `lambda = 0`.
#'
#' `pc` is the posterior coupling estimand of equation (3.8), and the quantity
#' to report. It is the same posterior tilted until the doubly robust moment
#' condition holds in posterior mean, which links the outcome model draws to the
#' propensity score model draws. It is consistent if either model is correctly
#' specified, which is the double robustness property, and it is only doubly
#' robust when the sweep actually met the moment condition. Both posteriors
#' average over the same empirical distribution of the covariates, so their
#' locations and their spreads are directly comparable.
#'
#' That empirical distribution is the third thing worth being clear about. Both
#' are posteriors for the average causal effect over the covariate vectors that
#' were observed, which is the average equation (3.7) takes, and not for the
#' average over the population those vectors came from. The two differ whenever
#' the fitted contrast varies from unit to unit; [drbayes_pc()] says how much,
#' and when.
#'
#' @name DRBayes-methods
#' @seealso [drbayes_pc()] for the fitting function, [rank_plot()] and
#'   [convergence_diagnostics()] for the underlying convergence tools.
NULL


#' Print a posterior coupling fit
#'
#' @description
#' A one screen summary of a [drbayes_pc()] fit: the call, the two models, the
#' data, the posterior mean and a 95 percent credible interval for both
#' estimands, the state of the sequential Monte Carlo sweep, and the worst
#' convergence diagnostics.
#'
#' @details
#' See [DRBayes-methods] for what separates `pc` from `g.comp`. Two things in
#' the output are warnings rather than description. If the sweep did not meet
#' the moment condition, the `pc` draws are not doubly robust and the message
#' says so; treat them as no better than `g.comp`. If a convergence threshold is
#' breached, the offending R-hat or effective sample size is marked, and the
#' full table is available through [summary.DRBayes()].
#'
#' @param x An object of class `"DRBayes"` from [drbayes_pc()].
#' @param digits Number of significant digits. The table of estimands is shown
#'   to the number of decimal places that gives the posterior standard deviation
#'   this many significant digits, so that a wide credible interval around a
#'   large effect does not round away. Defaults to 3.
#' @param ... Ignored, present for consistency with the generic.
#'
#' @return `x`, invisibly.
#'
#' @examples
#' # generate_dataset() leaves the caller's random number stream alone, so
#' # seed the coupling itself to make the fit below reproducible.
#' set.seed(1)
#' data <- generate_dataset(nn = 200, pp = 0, seed = 1)
#' fit  <- drbayes_pc(YY ~ AA + WW1 + WW2 + WW3 + WW4,
#'                    AA ~ WW1 + WW2 + WW3 + WW4,
#'                    data = data,
#'                    outcome.model = bayes_lm, ps.model = bayes_logit,
#'                    mc = 2000, bn = 500, verbose = FALSE,
#'                    control = drbayes_control(n_steps = 200))
#' fit
#'
#' @seealso [summary.DRBayes()], [plot.DRBayes()]
#' @export
print.DRBayes <- function(x, digits = 3, ...) {

  validate_drbayes_object(x)
  digits <- validate_digits(digits)

  cat("Bayesian doubly robust causal inference via posterior coupling\n")

  if (!is.null(x$call)) {
    cat("\nCall:\n")
    cat(paste0(deparse(x$call), collapse = "\n"), "\n", sep = "")
  }

  cat("\n")
  cat("Outcome model:    ", describe_outcome_model(x), "\n", sep = "")
  cat("Propensity score: ", describe_ps_model(x), "\n", sep = "")
  cat(describe_data(x), "\n", sep = "")
  cat("Posterior draws:  ", length(x$pc), "\n", sep = "")

  cat("\nAverage treatment effect\n")
  estimands <- drbayes_estimand_table(x)
  print(format_estimand_table(estimands, digits,
                              c("mean", "2.5%", "97.5%")))
  cat("\n")
  cat_wrapped(paste("Report pc, the doubly robust estimand. g.comp is the",
                    "untilted g-formula posterior, which relies on the outcome",
                    "model alone."))
  cat("\n")
  cat_wrapped(paste("Both are posteriors for the average causal effect over",
                    "the", x$data_info$n_observations, "covariate vectors",
                    "observed, which is what equation (3.7) averages. See",
                    "?drbayes_pc on when that differs from the average over",
                    "the population they were drawn from."))

  cat("\n")
  print_smc(x$smc, digits = digits)

  worst <- drbayes_worst_diagnostics(x$diagnostics)
  if (is.null(worst)) {
    cat("\nConvergence diagnostics: not computed (diagnostics = \"none\").\n")
  } else {
    cat("\n")
    print_worst_diagnostics(worst, drbayes_thresholds(x))
  }

  invisible(x)
}


#' Summarise a posterior coupling fit
#'
#' @description
#' Collects the posterior summaries of both estimands into a small data frame,
#' together with the sequential Monte Carlo diagnostics and the full convergence
#' diagnostics table, ready to be printed.
#'
#' @details
#' See [DRBayes-methods] for what separates `pc` from `g.comp`. The quantiles
#' are the usual equal-tailed ones, so `2.5\%` and `97.5\%` bound a 95 percent
#' credible interval and `50\%` is the posterior median.
#'
#' @param object An object of class `"DRBayes"` from [drbayes_pc()].
#' @param ... Ignored, present for consistency with the generic.
#'
#' @return An object of class `"summary.DRBayes"`, a list with components
#'   \describe{
#'     \item{estimands}{Data frame with one row per estimand, `pc` first, and
#'       columns `mean`, `sd`, `2.5\%`, `50\%` and `97.5\%`.}
#'     \item{call, family, link, ps.link, data_info, smc, diagnostics}{Taken
#'       unchanged from `object`.}
#'     \item{thresholds}{The R-hat and effective sample size thresholds the
#'       diagnostics are judged against.}
#'     \item{n_draws}{Number of posterior draws behind each estimand.}
#'   }
#'
#' @examples
#' # generate_dataset() leaves the caller's random number stream alone, so
#' # seed the coupling itself to make the fit below reproducible.
#' set.seed(1)
#' data <- generate_dataset(nn = 200, pp = 0, seed = 1)
#' fit  <- drbayes_pc(YY ~ AA + WW1 + WW2 + WW3 + WW4,
#'                    AA ~ WW1 + WW2 + WW3 + WW4,
#'                    data = data,
#'                    outcome.model = bayes_lm, ps.model = bayes_logit,
#'                    mc = 2000, bn = 500, verbose = FALSE,
#'                    control = drbayes_control(n_steps = 200))
#' summary(fit)
#'
#' # The table on its own, for example to put in a report.
#' summary(fit)$estimands
#'
#' @seealso [print.DRBayes()], [plot.DRBayes()]
#' @export
summary.DRBayes <- function(object, ...) {

  validate_drbayes_object(object)

  result <- list(
    estimands   = drbayes_estimand_table(object),
    call        = object$call,
    family      = object$family,
    link        = object$link,
    ps.link     = object$ps.link,
    data_info   = object$data_info,
    smc         = object$smc,
    diagnostics = object$diagnostics,
    thresholds  = drbayes_thresholds(object),
    n_draws     = length(object$pc)
  )

  class(result) <- c("summary.DRBayes", "list")
  result
}


#' @rdname summary.DRBayes
#'
#' @param x An object of class `"summary.DRBayes"`.
#' @param digits Number of significant digits. The table of estimands is shown
#'   to the number of decimal places that gives the posterior standard deviation
#'   this many significant digits, so that a wide credible interval around a
#'   large effect does not round away. Defaults to 3.
#'
#' @return `print.summary.DRBayes` returns `x` invisibly.
#'
#' @export
print.summary.DRBayes <- function(x, digits = 3, ...) {

  digits <- validate_digits(digits)

  cat("Bayesian doubly robust causal inference via posterior coupling\n")

  if (!is.null(x$call)) {
    cat("\nCall:\n")
    cat(paste0(deparse(x$call), collapse = "\n"), "\n", sep = "")
  }

  cat("\n")
  cat("Outcome model:    ", describe_outcome_model(x), "\n", sep = "")
  cat("Propensity score: ", describe_ps_model(x), "\n", sep = "")
  cat(describe_data(x), "\n", sep = "")
  cat("Posterior draws:  ", x$n_draws, "\n", sep = "")

  if (!is.null(x$data_info$outcome_formula)) {
    cat("Outcome formula:  ",
        paste0(deparse(x$data_info$outcome_formula), collapse = " "), "\n",
        sep = "")
  }
  if (!is.null(x$data_info$ps_formula)) {
    cat("PS formula:       ",
        paste0(deparse(x$data_info$ps_formula), collapse = " "), "\n", sep = "")
  }

  cat("\nPosterior summary of the average treatment effect\n")
  print(format_estimand_table(x$estimands, digits, names(x$estimands)))
  cat("\n")
  cat_wrapped(paste("pc is the doubly robust estimand of equation (3.8) and the",
                    "one to report. g.comp is the untilted g-formula posterior,",
                    "lambda = 0, which relies on the outcome model alone. Both",
                    "average over the same covariates, so their spreads are",
                    "comparable."))
  cat("\n")
  cat_wrapped(paste("The covariates they average over are the observed ones.",
                    "Where the fitted contrast varies from unit to unit, which",
                    "a treatment interaction or a non-identity link both cause,",
                    "this is not the same quantity as the average over the",
                    "population, and these intervals do not carry the",
                    "uncertainty of not knowing the covariate distribution.",
                    "?drbayes_pc gives the size of the difference."))

  cat("\nSequential Monte Carlo\n")
  print_smc(x$smc, digits = digits, verbose = TRUE)

  cat("\nConvergence diagnostics\n")
  if (is.null(x$diagnostics)) {
    cat("  Not computed (diagnostics = \"none\").\n")
  } else {
    print(x$diagnostics, digits = digits, row.names = FALSE)
    cat("\n")
    print_threshold_breaches(x$diagnostics, x$thresholds)
  }

  invisible(x)
}


#' Plot a posterior coupling fit
#'
#' @description
#' Posterior densities of the two average treatment effect estimands, or rank
#' plots of the underlying MCMC draws.
#'
#' @details
#' With `type = "density"` the posterior of `pc` is drawn over the posterior of
#' `g.comp` on a common axis. The gap between the two curves is the tilt: it is
#' how far the doubly robust moment condition has had to move the g-formula
#' posterior, and its size is summarised by the tilting parameter `lambda`
#' reported by [print.DRBayes()]. A gap much wider than the posterior spread
#' means the outcome model and the propensity score model disagree strongly, and
#' both should be inspected before either estimate is reported. Curves that lie
#' on top of one another mean the outcome model already satisfied the moment
#' condition, and `lambda` will be at or near zero.
#'
#' With `type = "rank"` the rank histograms of [rank_plot()] are drawn for one
#' coefficient of one of the two models. These need the per-chain draws of the
#' model coefficients, which a fit only carries when [drbayes_pc()] was asked to
#' keep them.
#'
#' Base graphics are used throughout, and the graphics parameters are restored
#' on exit.
#'
#' @param x An object of class `"DRBayes"` from [drbayes_pc()].
#' @param type `"density"` (default) for the overlaid posterior densities, or
#'   `"rank"` for rank plots of the MCMC draws.
#' @param model For `type = "rank"`, which model to plot, `"outcome"` (default)
#'   or `"ps"`.
#' @param parameter For `type = "rank"`, the coefficient to plot, given as a
#'   name or as an index. Defaults to the first.
#' @param col Two colours, for `pc` and for `g.comp` in that order. Used by
#'   `type = "density"` only.
#' @param lwd Line width for the density curves.
#' @param main Plot title. Defaults to a description of the estimand.
#' @param xlab Label of the horizontal axis.
#' @param legend Set to FALSE to leave the legend off the density plot.
#' @param ... Further graphical parameters passed to the underlying plot.
#'
#' @return Invisibly, `x` for `type = "density"`, or the list returned by
#'   [rank_plot()] for `type = "rank"`.
#'
#' @examples
#' # generate_dataset() leaves the caller's random number stream alone, so
#' # seed the coupling itself to make the fit below reproducible.
#' set.seed(1)
#' data <- generate_dataset(nn = 200, pp = 0, seed = 1)
#' fit  <- drbayes_pc(YY ~ AA + WW1 + WW2 + WW3 + WW4,
#'                    AA ~ WW1 + WW2 + WW3 + WW4,
#'                    data = data,
#'                    outcome.model = bayes_lm, ps.model = bayes_logit,
#'                    mc = 2000, bn = 500, verbose = FALSE,
#'                    control = drbayes_control(n_steps = 200))
#' plot(fit)
#'
#' @seealso [print.DRBayes()], [summary.DRBayes()], [rank_plot()]
#' @export
plot.DRBayes <- function(x, type = c("density", "rank"),
                         model = c("outcome", "ps"), parameter = 1L,
                         col = c("#2166AC", "#B2182B"), lwd = 2, main = NULL,
                         xlab = "Average treatment effect", legend = TRUE,
                         ...) {

  validate_drbayes_object(x)
  type <- match.arg(type)

  if (type == "rank") {
    return(invisible(plot_drbayes_rank(x, model = match.arg(model),
                                       parameter = parameter, ...)))
  }

  if (length(col) < 2L) {
    stop("col must give two colours, one for pc and one for g.comp")
  }
  if (!is.logical(legend) || length(legend) != 1L || is.na(legend)) {
    stop("legend must be TRUE or FALSE")
  }

  pc_density <- finite_density(x$pc, "pc")
  g_density  <- finite_density(x$g.comp, "g.comp")

  # CRAN policy: leave the user's graphics state as we found it, whatever
  # happens in between.
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mar = c(4.1, 4.1, 3.1, 1.1), mgp = c(2.4, 0.8, 0))

  # A shared axis in both directions, since the point of the plot is the size
  # of the shift from one curve to the other.
  xlim <- range(pc_density$x, g_density$x)
  ylim <- c(0, max(pc_density$y, g_density$y) * 1.05)

  if (is.null(main)) {
    main <- "Posterior of the average treatment effect"
  }

  graphics::plot.default(NULL, xlim = xlim, ylim = ylim, xlab = xlab,
                         ylab = "Posterior density", main = main,
                         yaxs = "i", bty = "n", ...)
  graphics::lines(g_density, col = col[2L], lwd = lwd, lty = 2)
  graphics::lines(pc_density, col = col[1L], lwd = lwd)
  graphics::box(bty = "l")

  if (legend) {
    graphics::legend("topright",
                     legend = c("pc (doubly robust)", "g.comp (untilted)"),
                     col = col[1:2], lwd = lwd, lty = c(1, 2), bty = "n")
  }

  invisible(x)
}


#' Rank plot branch of plot.DRBayes
#'
#' Split out so that the density branch, which is the common case, stays
#' readable. The per-chain coefficient draws are needed here and are not part of
#' a fit unless they were explicitly kept, so the failure has to name the way to
#' get them rather than just reporting that they are absent.
#'
#' @return The list returned by [rank_plot()].
#'
#' @keywords internal
#' @noRd
plot_drbayes_rank <- function(x, model, parameter, ...) {

  stored <- x$draws
  draws  <- if (is.list(stored) && model %in% names(stored)) {
    stored[[model]]
  } else {
    NULL
  }

  if (is.null(draws)) {
    stop("Rank plots need the per-chain draws of the model coefficients, and ",
         "this fit does not carry them: drbayes_pc() returns only the draws of ",
         "the treatment effect unless it is asked to keep them. Refit with ",
         "drbayes_pc(..., control = drbayes_control(keep_draws = TRUE)), or ",
         "call rank_plot() directly on ",
         "the iterations by chains by parameters array that the sampler, such ",
         "as bayes_lm(), returns.", call. = FALSE)
  }

  rank_plot(draws, parameter = parameter, ...)
}


#' Posterior summaries of both estimands
#'
#' `pc` comes first because it is the estimand to report.
#'
#' @return A data frame with rows `pc` and `g.comp`.
#'
#' @keywords internal
#' @noRd
drbayes_estimand_table <- function(x) {
  probs <- c(0.025, 0.5, 0.975)

  summarise_one <- function(draws) {
    quantiles <- stats::quantile(draws, probs = probs, names = FALSE,
                                 na.rm = TRUE)
    c(mean = mean(draws, na.rm = TRUE), sd = stats::sd(draws, na.rm = TRUE),
      quantiles)
  }

  rows <- rbind(pc = summarise_one(x$pc), g.comp = summarise_one(x$g.comp))
  out  <- as.data.frame(rows)
  names(out) <- c("mean", "sd", "2.5%", "50%", "97.5%")
  out
}


#' Round an estimand table to a readable number of decimal places
#'
#' Rounding to significant digits, which is what `print.data.frame()` does, is
#' the wrong choice here: an average treatment effect of 110 with a posterior
#' standard deviation of 0.15 prints as three copies of 110 and the credible
#' interval disappears. The number of decimal places is instead set by the
#' posterior standard deviation, and shared by both rows so that the two
#' estimands stay comparable column by column.
#'
#' @param columns Which columns of `estimands` to keep, in order.
#'
#' @return A data frame of formatted strings, printed right aligned.
#'
#' @keywords internal
#' @noRd
format_estimand_table <- function(estimands, digits, columns) {
  scale <- max(estimands$sd, na.rm = TRUE)
  decimals <- if (!is.finite(scale) || scale <= 0) {
    digits
  } else {
    min(10L, max(0L, digits - 1L - as.integer(floor(log10(scale)))))
  }

  out <- lapply(estimands[, columns, drop = FALSE], formatC,
                format = "f", digits = decimals)
  out <- as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
  names(out)    <- columns
  rownames(out) <- rownames(estimands)
  out
}


#' Convergence thresholds this fit was judged against
#'
#' The thresholds live in the recorded call when the user set them, and default
#' to the values of Vehtari et al. (2021) that [drbayes_pc()] uses otherwise.
#' Reading them back rather than hard coding them here keeps the printed flags
#' in step with the warning the fit itself issued.
#'
#' @return A list with `rhat_max` and `ess_min`.
#'
#' @keywords internal
#' @noRd
drbayes_thresholds <- function(x) {
  from_call <- function(name, default) {
    value <- x$call[[name]]
    if (is.numeric(value) && length(value) == 1L && is.finite(value)) {
      return(value)
    }
    default
  }
  list(rhat_max = from_call("rhat_max", 1.01),
       ess_min  = from_call("ess_min", 400))
}


#' Worst R-hat and smallest effective sample size over both models
#'
#' A missing diagnostic is a failure rather than an absence, so it outranks any
#' finite value here: dropping it would let a parameter that could not be
#' assessed pass silently.
#'
#' @return NULL when there are no diagnostics, otherwise a list with components
#'   `rhat` and `ess`, each a list of `value` and `where`.
#'
#' @keywords internal
#' @noRd
drbayes_worst_diagnostics <- function(diagnostics) {
  if (!is.data.frame(diagnostics) || nrow(diagnostics) == 0L) {
    return(NULL)
  }

  label <- if (is.null(diagnostics$model)) {
    as.character(diagnostics$parameter)
  } else {
    paste0(diagnostics$model, ": ", diagnostics$parameter)
  }

  worst <- function(values, largest_is_worst) {
    if (anyNA(values)) {
      return(list(value = NA_real_, where = label[which(is.na(values))[1L]]))
    }
    index <- if (largest_is_worst) which.max(values) else which.min(values)
    list(value = values[index], where = label[index])
  }

  list(rhat = worst(diagnostics$rhat, TRUE),
       ess  = worst(pmin(diagnostics$ess_bulk, diagnostics$ess_tail), FALSE))
}


#' One line descriptions of the fitted models and the data
#'
#' @keywords internal
#' @noRd
describe_outcome_model <- function(x) {
  family <- if (is.null(x$family)) "unknown" else x$family
  link   <- if (is.null(x$link)) "unknown" else x$link
  paste0(family, " family, ", link, " link")
}

#' @keywords internal
#' @noRd
describe_ps_model <- function(x) {
  paste0(if (is.null(x$ps.link)) "unknown" else x$ps.link, " link")
}

#' @keywords internal
#' @noRd
describe_data <- function(x) {
  info <- x$data_info
  if (is.null(info)) {
    return("Observations:     unknown")
  }
  line <- paste0("Observations:     ", info$n_observations, " (",
                 info$n_treated, " treated, ", info$n_control, " control)")
  if (!is.null(info$missing_observations) && info$missing_observations > 0) {
    line <- paste0(line, ", ", info$missing_observations,
                   " dropped for missing values")
  }
  line
}


#' Report the state of the sequential Monte Carlo sweep
#'
#' A sweep that never met the moment condition is the one result in the object
#' that invalidates the headline estimand, so it is stated as a warning rather
#' than as another diagnostic number.
#'
#' @param verbose Add the individual quantities on their own lines, for
#'   `summary()`.
#'
#' @keywords internal
#' @noRd
print_smc <- function(smc, digits, verbose = FALSE) {
  if (is.null(smc)) {
    cat("Sequential Monte Carlo: no information recorded.\n")
    return(invisible(NULL))
  }

  fmt <- function(value) format(value, digits = digits)

  if (verbose) {
    cat("  Tilting parameter lambda: ", fmt(smc$lambda), "\n", sep = "")
    cat("  Posterior mean of B_n:    ", fmt(smc$B_mean), "\n", sep = "")
    cat("  Tolerance:                ", fmt(smc$tol), "\n", sep = "")
    cat("  Steps taken:              ", smc$n_steps, " of ", smc$max_steps,
        "\n", sep = "")
    cat("  Moment condition met:     ",
        if (isTRUE(smc$converged)) "yes" else "NO", "\n", sep = "")
  }

  if (isTRUE(smc$converged)) {
    if (!verbose) {
      cat_wrapped(paste0("Posterior coupling: lambda = ", fmt(smc$lambda),
                         ", moment condition met, |mean B_n| = ",
                         fmt(abs(smc$B_mean)), " below the tolerance ",
                         fmt(smc$tol), ", after ", smc$n_steps, " of ",
                         smc$max_steps, " steps."))
    }
    return(invisible(NULL))
  }

  cat("\n")
  cat("*** THE MOMENT CONDITION WAS NOT MET ***\n")
  cat_wrapped(paste0("The sweep stopped at lambda = ", fmt(smc$lambda),
                     " after ", smc$n_steps, " of ", smc$max_steps,
                     " steps with |mean B_n| = ", fmt(abs(smc$B_mean)),
                     ", still above the tolerance ", fmt(smc$tol),
                     ". The pc draws are therefore NOT doubly robust and carry ",
                     "no guarantee beyond what the outcome model alone gives. ",
                     "Check the propensity score model for near zero or near ",
                     "one fitted values, and see ?drbayes_pc."),
              prefix = "  ")
  invisible(NULL)
}


#' Report the worst convergence diagnostics, flagged against the thresholds
#'
#' @keywords internal
#' @noRd
print_worst_diagnostics <- function(worst, thresholds) {
  rhat_bad <- is.na(worst$rhat$value) || worst$rhat$value > thresholds$rhat_max
  ess_bad  <- is.na(worst$ess$value)  || worst$ess$value  < thresholds$ess_min

  # R-hat lives within a percent or so of 1, so significant digits would print
  # the interesting part away; fixed decimals keep it readable.
  rhat_text <- trimws(formatC(worst$rhat$value, format = "f", digits = 3))
  ess_text  <- trimws(formatC(round(worst$ess$value), format = "d"))

  flag <- function(bad, value, side, threshold) {
    if (!bad) return("")
    if (is.na(value)) return("  <- could not be computed, treat as a failure")
    paste0("  <- ", side, " ", threshold)
  }

  cat("Convergence diagnostics (worst over both models)\n")
  cat("  Largest R-hat: ", rhat_text, " (", worst$rhat$where, ")",
      flag(rhat_bad, worst$rhat$value, "above", thresholds$rhat_max),
      "\n", sep = "")
  cat("  Smallest ESS:  ", ess_text, " (", worst$ess$where, ")",
      flag(ess_bad, worst$ess$value, "below", thresholds$ess_min),
      "\n", sep = "")

  if (rhat_bad || ess_bad) {
    cat_wrapped(paste("The draws do not meet the convergence criteria, so the",
                      "estimates above are not to be trusted. Raise mc, or use",
                      "more chains from more dispersed starting values. See",
                      "summary() for the full table."),
                prefix = "  ")
  }
  invisible(NULL)
}


#' Name the parameters that breach a convergence threshold
#'
#' @keywords internal
#' @noRd
print_threshold_breaches <- function(diagnostics, thresholds) {
  failures <- flag_convergence_failures(diagnostics,
                                        rhat_max = thresholds$rhat_max,
                                        ess_min  = thresholds$ess_min)

  if (isTRUE(failures$ok)) {
    cat_wrapped(paste0("All parameters meet R-hat < ", thresholds$rhat_max,
                       " and ESS > ", thresholds$ess_min, "."))
    return(invisible(NULL))
  }

  cat_wrapped(paste0("Thresholds: R-hat < ", thresholds$rhat_max, ", ESS > ",
                     thresholds$ess_min, "."))
  if (length(failures$rhat) > 0) {
    cat_wrapped(paste0("R-hat too high: ", toString(failures$rhat), "."),
                prefix = "  ")
  }
  if (length(failures$ess) > 0) {
    cat_wrapped(paste0("ESS too low: ", toString(failures$ess), "."),
                prefix = "  ")
  }
  if (isTRUE(failures$rhat_unreliable)) {
    cat_wrapped(paste("An effective sample size below the threshold also makes",
                      "R-hat itself unreliable, so a small R-hat here is not",
                      "evidence of convergence."),
                prefix = "  ")
  }
  invisible(NULL)
}


#' Write a paragraph wrapped to the console width
#'
#' @keywords internal
#' @noRd
cat_wrapped <- function(text, prefix = "") {
  width <- max(40L, getOption("width", 80L) - nchar(prefix))
  cat(paste0(strwrap(text, width = width, prefix = prefix), collapse = "\n"),
      "\n", sep = "")
}


#' Kernel density of a set of draws, with the non-finite ones removed
#'
#' [stats::density()] fails outright on a missing value, which would turn a
#' plotting call into an error message about an argument the caller never
#' supplied. Dropping them and saying how many were dropped is more useful.
#'
#' @keywords internal
#' @noRd
finite_density <- function(draws, name) {
  keep <- is.finite(draws)
  if (!any(keep)) {
    stop("All ", length(draws), " draws of ", name, " are missing or infinite, ",
         "so there is nothing to plot.")
  }
  if (!all(keep)) {
    warning(sum(!keep), " of ", length(draws), " draws of ", name,
            " are missing or infinite and are left out of the density.")
  }
  stats::density(draws[keep])
}


#' Check that an object carries what the methods need
#'
#' @keywords internal
#' @noRd
validate_drbayes_object <- function(x) {
  for (name in c("pc", "g.comp")) {
    draws <- x[[name]]
    if (is.null(draws)) {
      stop("This DRBayes object has no '", name, "' component. It was not ",
           "produced by drbayes_pc(), or it predates the current version of ",
           "the package.", call. = FALSE)
    }
    if (!is.numeric(draws) || length(draws) == 0L) {
      stop("The '", name, "' component of this DRBayes object must be a ",
           "non-empty numeric vector of posterior draws, but is ",
           class(draws)[1L], " of length ", length(draws), ".", call. = FALSE)
    }
  }
  if (length(x$pc) != length(x$g.comp)) {
    stop("The 'pc' and 'g.comp' components hold ", length(x$pc), " and ",
         length(x$g.comp), " draws. They come from the same sweep and must be ",
         "the same length.", call. = FALSE)
  }
  invisible(TRUE)
}


#' @keywords internal
#' @noRd
validate_digits <- function(digits) {
  if (!is.numeric(digits) || length(digits) != 1L || !is.finite(digits) ||
      digits < 1) {
    stop("digits must be a single number of at least 1")
  }
  as.integer(digits)
}
