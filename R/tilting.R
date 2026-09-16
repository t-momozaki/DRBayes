# The seven points is_logistic() has always probed at. Held in one place so the
# propensity score and outcome classifications cannot drift apart.
link_probe <- c(-31, -4.25, -0.5, 0, 1.75, 6, 22)

#' Which inverse link a closure is, if it is one the kernels know
#'
#' The test is identical(), as it has always been, so a closure that merely
#' agrees with plogis to within rounding is not treated as plogis: the fused
#' weight is exact only for the logistic link itself.
#'
#' @param value The closure's value at [link_probe].
#'
#' @return 0 for the identity, 1 for the logit, 2 for the probit, -1 for
#'   anything else, which takes the reference implementation instead.
#'
#' @keywords internal
#' @noRd
classify_link <- function(value) {
  if (identical(value, stats::plogis(link_probe))) return(1L)
  if (identical(value, stats::pnorm(link_probe)))  return(2L)
  if (identical(value, link_probe))                return(0L)
  -1L
}

#' @keywords internal
#' @noRd
ps_link_code <- function(f) classify_link(f(link_probe))

#' The same for the outcome link, which has never been probed before
#'
#' Wrapped, because an inverse link that only accepts the matrix the tilting
#' hands it is legal today and must not start failing now that it is asked a
#' question it was never asked.
#'
#' @keywords internal
#' @noRd
outcome_link_code <- function(f) {
  tryCatch(classify_link(f(link_probe)), error = function(e) -1L)
}

# The cell budget draw_blocks() spends, shared with the compiled kernels so
# that both group the draws the same way.
moment_max_cells <- 2^18

#' Whether R sums in long double on this build
#'
#' `colMeans()`, `rowMeans()` and `mean()` accumulate in `LDOUBLE`, which is
#' `long double` where R was built with it and `double` where it was not. The
#' compiled kernels have to follow that choice or they stop agreeing with the R
#' they replace, and they cannot read it for themselves: the macro that records
#' it lives in R's own `config.h` and is not among the fifteen that the
#' installed `Rconfig.h` exports. Asking R also settles the case a compile time
#' test could not have got right in any form, CRAN's no-long-double flavour,
#' where the platform has a wide `long double` and R has been built not to use
#' it.
#'
#' @return `TRUE` or `FALSE`.
#'
#' @keywords internal
#' @noRd
long_double_sums <- function() isTRUE(capabilities("long.double"))


#' The moment condition of equation (3.4)
#'
#' The mean over observations of
#' \deqn{(A - e(X; \alpha)) / (e(X; \alpha) (1 - e(X; \alpha))) (Y - m_A(X; \beta))}
#' evaluated at every posterior draw, giving one value of \eqn{B_n} per draw.
#'
#' This is the whole cost of the sweep of [tilt_smc()], which evaluates it once
#' a step, so it is written for the two things that cost: the weight is formed
#' by [ipw_weight()], and the draws are taken in blocks by [draw_blocks()] so
#' that nothing the size of the sample by the draws is ever held.
#'
#' @param betas.ps,betas.otc Draws by parameters matrices of coefficients.
#' @param d The data and link functions assembled by [tilting_data()].
#'
#' @return A numeric vector with one entry per draw.
#'
#' @keywords internal
#' @noRd
moment_ipw_r <- function(betas.ps, betas.otc, d) {
  n         <- length(d$Y)
  weight_at <- ipw_weight(d)
  BB        <- numeric(nrow(betas.ps))

  for (draws in draw_blocks(nrow(betas.ps), n)) {
    ww <- weight_at(tcrossprod(d$Z.ps, betas.ps[draws, , drop = FALSE]))
    mu <- d$inverse_link(tcrossprod(d$Z.lm, betas.otc[draws, , drop = FALSE]))
    BB[draws] <- .colMeans(ww * (d$Y - mu), n, length(draws))
  }
  BB
}


#' The inverse probability weight of equation (3.4)
#'
#' The factor \eqn{(A - e(X; \alpha)) / (e(X; \alpha) (1 - e(X; \alpha)))} that
#' equation (3.4) applies to the outcome residual, as a function of the
#' propensity score linear predictor. Both [moment_ipw()] and [moment_np()]
#' weight by this, so the two agree on the moment condition by construction.
#'
#' Under the logistic link the factor collapses. Writing \eqn{s = 2A - 1},
#' which is 1 for a treated observation and -1 for a control one since the
#' treatment is coded 0/1,
#' \deqn{(A - e) / (e (1 - e)) = s (1 + \exp(-s \eta)),}
#' because at \eqn{A = 1} the left hand side is
#' \eqn{(1 - e) / (e (1 - e)) = 1 / e = 1 + \exp(-\eta)} and at \eqn{A = 0} it
#' is \eqn{-e / (e (1 - e)) = -1 / (1 - e) = -(1 + \exp(\eta))}. One
#' exponential then does the work of a [stats::plogis()], two subtractions, a
#' multiplication and a division, and of the matrices each of those would
#' allocate.
#'
#' It is also the more accurate of the two where accuracy is scarce. At
#' \eqn{A = 1} the divided form loses nothing, since the \eqn{1 - e} above and
#' below the line are the same rounded number and cancel. At \eqn{A = 0} it is
#' \eqn{-1 / (1 - e)}, and \eqn{1 - e} is a difference of nearly equal numbers
#' once the score approaches one: at \eqn{\eta = 10} that form is already wrong
#' in the thirteenth digit and at \eqn{\eta = 30} in the third, where
#' \eqn{1 + \exp(\eta)} is right to the last bit throughout. Those are the
#' observations carrying the largest weights, so they are the ones the moment
#' condition is mostly made of.
#'
#' @param d The object returned by [tilting_data()].
#'
#' @return A function of the linear predictors, of any shape, returning the
#'   weights in that shape. It refuses a positivity violation before returning
#'   anything.
#'
#' @keywords internal
#' @noRd
ipw_weight <- function(d) {
  if (!is_logistic(d$ps_inverse_link)) {
    return(function(eta) {
      check_positivity(eta, d)
      ps <- d$ps_inverse_link(eta)
      (d$A - ps) / (ps * (1 - ps))
    })
  }

  s <- 2 * d$A - 1
  function(eta) {
    check_positivity(eta, d)
    s * (1 + exp(-s * eta))
  }
}


#' Whether an inverse link is the logistic one
#'
#' [ipw_weight()] has a cheaper equal form under the logistic link, so the
#' tilting asks the inverse link it was handed, at a spread of linear
#' predictors, whether it is [stats::plogis()]. Every inverse link this package
#' builds is made from a link name and so answers exactly; anything else is
#' taken at its word and weighted by the general form.
#'
#' @param inverse_link The inverse link to identify.
#'
#' @return `TRUE` when the function returns exactly what [stats::plogis()]
#'   returns at each probe.
#'
#' @keywords internal
#' @noRd
is_logistic <- function(inverse_link) {
  eta <- c(-31, -4.25, -0.5, 0, 1.75, 6, 22)
  identical(inverse_link(eta), stats::plogis(eta))
}


#' The draws, in blocks small enough to keep the working matrices small
#'
#' Equation (3.4) is a mean over observations at each draw, so the draws do not
#' interact and the moment condition can be evaluated a block of them at a
#' time. Each block holds a handful of matrices of observations by draws, so
#' blocking bounds those at a few megabytes where the whole sample of draws
#' would need a few hundred.
#'
#' Where the blocks fall is not quite immaterial. Column means are formed one
#' column at a time whatever the block, so the reduction does not care; the
#' matrix product does, because a BLAS may reach a different kernel for a
#' differently shaped product. Two different partitions therefore agree to
#' rounding rather than exactly. The partition is a function of the data alone,
#' so a fit is reproducible, and the compiled kernels use this same one so that
#' they reproduce this implementation exactly.
#'
#' @param S Number of draws.
#' @param n Number of observations.
#' @param max_cells Most entries a working matrix may hold.
#'
#' @return A list of index vectors partitioning `seq_len(S)` in order.
#'
#' @keywords internal
#' @noRd
draw_blocks <- function(S, n, max_cells = 2^18) {
  size <- max(1L, as.integer(max_cells %/% max(n, 1L)))
  lapply(seq_len(ceiling(S / size)),
         function(b) seq.int((b - 1L) * size + 1L, min(b * size, S)))
}


#' The moment condition of equation (3.4) from supplied fitted means
#'
#' The same quantity as [moment_ipw()], for an outcome model that reports the
#' fitted mean \eqn{m_{A_i}(X_i; \beta)} itself rather than the coefficients it
#' was built from, so that a model with no linear predictor to evaluate, such
#' as the BART fit of Section 5.4, can be coupled as well.
#'
#' The fitted means arrive as draws by observations, the transpose of the
#' layout the rest of the tilting works in, and the condition is evaluated a
#' draw at a time rather than for all draws at once. At the scale of Section
#' 5.4 the fitted means run to tens of megabytes each, and transposing them, or
#' forming a matrix of propensity scores to match, would cost as much again;
#' taking a draw at a time leaves the working memory proportional to the number
#' of observations alone.
#'
#' @param betas.ps Draws by parameters matrix of propensity score coefficients.
#' @param d The object returned by [tilting_data()], carrying `outcome.mu` and,
#'   for equation (F.2), the stratum weights added by [subclass_weights()].
#'
#' @return A numeric vector with one entry per draw.
#'
#' @keywords internal
#' @noRd
moment_np_r <- function(betas.ps, d) {
  mu <- d$outcome.mu

  weight_at <- if (is.null(d$subclass_weight)) {
    ipw <- ipw_weight(d)
    function(s) ipw(drop(d$Z.ps %*% betas.ps[s, ]))
  } else {
    # The stratum treated fraction of equation (F.2) is fixed across draws.
    function(s) d$subclass_weight
  }

  vapply(seq_len(nrow(mu)),
         function(s) mean(weight_at(s) * (d$Y - mu[s, ])),
         numeric(1))
}


#' The moment condition of equation (3.4)
#'
#' Each of these is a thin wrapper over the reference implementation in R. The
#' compiled kernels are introduced behind them one at a time, so that the R that
#' preceded them stays available as the thing the compiled path is compared
#' against, and as the path a user-supplied inverse link still takes.
#'
#' @inheritParams moment_ipw_r
#'
#' @return A numeric vector with one entry per draw, unnamed.
#'
#' @keywords internal
#' @noRd
moment_ipw <- function(betas.ps, betas.otc, d) {
  ps_code  <- ps_link_code(d$ps_inverse_link)
  otc_code <- outcome_link_code(d$inverse_link)
  if (ps_code < 0L || otc_code < 0L) {
    return(moment_ipw_r(betas.ps, betas.otc, d))
  }

  out <- drb_moment_ipw(d$Z.ps, d$Z.lm, d$A, d$Y, betas.ps, betas.otc,
                        ps_code, otc_code, moment_max_cells,
                        long_double_sums())

  # min() and max() of a pair are the pair, so the guard sees exactly what it
  # would have seen had it been handed the whole matrix, and its wording stays
  # in the one place it has always been.
  if (length(out$eta_range)) check_positivity(out$eta_range, d)
  out$BB
}


#' The subclassification moment condition of equation (F.2)
#'
#' @inheritParams moment_subclass_r
#' @return A numeric vector with one entry per draw, unnamed.
#'
#' @keywords internal
#' @noRd
moment_subclass <- function(betas.ps, betas.otc, d) {
  if (is.null(d$subclass_weight)) {
    return(moment_subclass_r(betas.ps, betas.otc, d))
  }
  otc_code <- outcome_link_code(d$inverse_link)
  if (otc_code < 0L) {
    return(moment_subclass_r(betas.ps, betas.otc, d))
  }
  drb_moment_subclass(d$Z.lm, d$Y, d$subclass_weight, betas.otc, otc_code,
                      moment_max_cells, long_double_sums())
}


#' The moment condition from supplied draws of the fitted outcome mean
#'
#' @inheritParams moment_np_r
#' @return A numeric vector with one entry per draw, unnamed.
#'
#' @keywords internal
#' @noRd
moment_np <- function(betas.ps, d) {
  if (!is.null(d$subclass_weight)) {
    return(drb_moment_np_subclass(d$Y, d$outcome.mu, d$subclass_weight,
                                  long_double_sums()))
  }
  ps_code <- ps_link_code(d$ps_inverse_link)
  if (ps_code < 0L) {
    return(moment_np_r(betas.ps, d))
  }
  out <- drb_moment_np_ipw(d$Z.ps, d$A, d$Y, d$outcome.mu, betas.ps, ps_code,
                           long_double_sums())
  if (length(out$eta_range)) check_positivity(out$eta_range, d)
  out$BB
}



#' The subclassification moment condition of equation (F.2)
#'
#' Appendix F replaces the inverse propensity score weights of equation (3.4)
#' by the treated fraction of a propensity score stratum, \eqn{n_{k1} / n_{k+}},
#' which is the more stable of the two because an extreme individual score is
#' smoothed by the stratum it falls in.
#'
#' The strata are built once by [subclass_weights()] and carried in `d`, so the
#' propensity score draws enter this condition only through them and
#' `betas.ps` itself is not used here.
#'
#' @param betas.ps,betas.otc Draws by parameters matrices of coefficients.
#' @param d The data assembled by [tilting_data()], carrying the stratum
#'   weights added by [subclass_weights()].
#'
#' @return A numeric vector with one entry per draw.
#'
#' @keywords internal
#' @noRd
moment_subclass_r <- function(betas.ps, betas.otc, d) {
  if (is.null(d$subclass_weight)) {
    stop("The subclassification moment condition needs propensity score ",
         "strata, which are built once before the tilting starts. Reach it ",
         "through run_tilting() with ",
         "control = drbayes_control(moment = \"subclass\") rather than by ",
         "calling this function directly.")
  }
  n  <- length(d$Y)
  BB <- numeric(nrow(betas.otc))

  for (draws in draw_blocks(nrow(betas.otc), n)) {
    mu <- d$inverse_link(tcrossprod(d$Z.lm, betas.otc[draws, , drop = FALSE]))
    BB[draws] <- .colMeans(d$subclass_weight * (d$Y - mu), n, length(draws))
  }
  BB
}


#' Fixed propensity score strata for the subclassification moment condition
#'
#' Equal frequency strata, the quantiles of the estimated propensity score, as
#' in equation (F.1). Within stratum \eqn{k} the treated fraction
#' \eqn{n_{k1} / n_{k+}} plays the part of the propensity score, so an
#' observation in that stratum contributes
#' \eqn{A_i / (n_{k1} / n_{k+}) - (1 - A_i) / (1 - n_{k1} / n_{k+})} to
#' equation (F.2). The stratum sizes \eqn{n_{k+}} are counted rather than taken
#' to be \eqn{n / K}, since quantile strata are exactly equal in size only when
#' the scores have no ties.
#'
#' The boundaries are computed once, from the posterior mean propensity score
#' of the untilted draws, and then held fixed. Restratifying at each draw would
#' make \eqn{B_n} a different function of the parameters at every step, so the
#' sweep over the tilting parameter would no longer be solving one fixed
#' equation.
#'
#' @param betas.ps Draws by parameters matrix of propensity score coefficients.
#' @param d The data assembled by [tilting_data()].
#' @param n_subclass Number of strata, \eqn{K} in equation (F.1).
#'
#' @return A numeric vector of length \eqn{n}, the stratum weight of each
#'   observation.
#'
#' @keywords internal
#' @noRd
subclass_weights <- function(betas.ps, d, n_subclass) {
  ps_code <- ps_link_code(d$ps_inverse_link)
  if (ps_code < 0L) {
    eta <- tcrossprod(d$Z.ps, betas.ps)
    check_positivity(eta, d)
    ps      <- d$ps_inverse_link(eta)
    ps_mean <- .rowMeans(ps, nrow(ps), ncol(ps))
  } else {
    out <- drb_ps_rowmeans(d$Z.ps, betas.ps, ps_code, moment_max_cells,
                           long_double_sums())
    if (length(out$eta_range)) check_positivity(out$eta_range, d)
    ps_mean <- out$ps_mean
  }

  # c_0 = 0 and c_K = 1 as in (F.1); the interior boundaries are the quantiles
  # that make the strata equal in frequency.
  interior <- stats::quantile(ps_mean, names = FALSE,
                              probs = seq_len(n_subclass - 1L) / n_subclass)
  breaks   <- c(0, interior, 1)
  if (any(diff(breaks) <= 0)) {
    stop("The estimated propensity score takes too few distinct values to ",
         "form ", n_subclass, " equal frequency strata: two stratum ",
         "boundaries coincide. Use a smaller n_subclass, or moment = ",
         "\"ipw\", which does not stratify.")
  }

  stratum   <- cut(ps_mean, breaks = breaks, labels = FALSE, right = FALSE,
                   include.lowest = TRUE)
  n_total   <- tabulate(stratum, nbins = n_subclass)
  n_treated <- tabulate(stratum[d$A == 1], nbins = n_subclass)

  unusable <- which(n_total == 0L | n_treated == 0L | n_treated == n_total)
  if (length(unusable) > 0L) {
    k      <- unusable[1L]
    reason <- if (n_total[k] == 0L) {
      "no units at all"
    } else if (n_treated[k] == 0L) {
      "no treated units"
    } else {
      "no control units"
    }
    stop("Propensity score stratum ", k, " of ", n_subclass,
         ", covering scores in [", format(breaks[k], digits = 3), ", ",
         format(breaks[k + 1L], digits = 3), "), contains ", reason, ", so ",
         "equation (F.2) divides by zero there. This is a positivity ",
         "violation within the stratum: use a smaller n_subclass, or ",
         "moment = \"ipw\", which does not stratify.")
  }

  ps_stratum <- (n_treated / n_total)[stratum]
  d$A / ps_stratum - (1 - d$A) / (1 - ps_stratum)
}




#' Refuse propensity scores that equation (3.4) cannot be evaluated at
#'
#' The refusal is decided from the extremes of the linear predictor rather than
#' from every propensity score. The inverse link is increasing, so the smallest
#' and the largest score in a whole matrix of draws are the ones its smallest
#' and largest linear predictors carry, and testing those two settles the same
#' question as testing all \eqn{n S} of them: one pass over the linear
#' predictors, with no logical matrix the size of the draws to hold the answer.
#' `min()` and `max()` carry NA and NaN through, so a linear predictor that is
#' not a number still arrives here as a score that is not finite.
#'
#' That the weight of [ipw_weight()] has no division left to overflow is what
#' makes this necessary rather than redundant. A propensity score of 0 or 1
#' means the covariates separate the treatment groups, and the moment condition
#' there is not finite. Written as the division of equation (3.4) it announces
#' that itself, as an infinity or a 0/0; the fused weight returns an
#' unremarkable looking number instead, so this guard is the only thing between
#' a separated fit and a plausible wrong answer.
#'
#' @param eta Propensity score linear predictors, of any shape.
#' @param d The object returned by [tilting_data()], for the inverse link and
#'   for the formula the error tells the caller to look at.
#'
#' @return `NULL`, invisibly. Called for the error it raises.
#'
#' @keywords internal
#' @noRd
check_positivity <- function(eta, d) {
  ps <- d$ps_inverse_link(c(min(eta), max(eta)))
  if (any(!is.finite(ps)) || any(ps <= 0 | ps >= 1)) {
    stop("The propensity score reached 0 or 1 for some draws, so the moment ",
         "condition is not finite. This is a positivity violation: check ",
         "range(fitted(glm(", d$ps_formula_text, ", binomial, data))) and ",
         "consider dropping covariates that separate the treatment groups.")
  }
  invisible(NULL)
}


#' Everything the tilting needs from the model frame
#'
#' Collecting these into one object keeps the tilting engine independent of how
#' the design matrices were built, so it can be driven by [drbayes_pc()], by the
#' sensitivity analysis, or by the confounder selection without change.
#'
#' `outcome.mu` is the draws by observations matrix of fitted means of the
#' nonparametric path, and is carried only when there is one; the parametric
#' path evaluates `inverse_link(Z.lm %*% beta)` instead.
#'
#' `ps_inverse_link` is an increasing function, as the inverse of a link for a
#' binary response is; [check_positivity()] reads the extremes of a set of
#' propensity scores off the extremes of their linear predictors on that basis.
#'
#' @keywords internal
#' @noRd
tilting_data <- function(Z.lm, Z.ps, A, Y, inverse_link, ps_inverse_link,
                         ps_formula_text, outcome.mu = NULL) {
  d <- list(Z.lm = Z.lm, Z.ps = Z.ps, A = A, Y = Y,
            inverse_link = inverse_link, ps_inverse_link = ps_inverse_link,
            ps_formula_text = ps_formula_text)
  if (!is.null(outcome.mu)) d$outcome.mu <- outcome.mu
  d
}


#' Newton step for the tilting parameter, from lambda = 0
#'
#' One step of the Newton iteration in Algorithm 1, taken where the draws still
#' carry equal weight. Only its sign is used, to fix which half of the lambda
#' grid the sweep travels along. The ridge keeps it defined when the draws put
#' all the mass at \eqn{B_n = 0}.
#'
#' Numerator and denominator are averages over the draws rather than the sums
#' of Algorithm 1, which leaves the step itself unchanged in sign and puts the
#' denominator on the same scale as the weighted second moment [tilt_is()]
#' divides by, so that one value of `ridge` means the same thing in both.
#' Written as sums it would be \eqn{S} times less influential here than there.
#'
#' @param BB The moment condition at the draws.
#' @param ridge Small value added to the denominator.
#'
#' @return The Newton step, a single number.
#'
#' @keywords internal
#' @noRd
newton_lambda <- function(BB, ridge) {
  -mean(BB) / (mean(BB^2) + ridge)
}


#' Tolerance on the posterior mean of the moment condition
#'
#' The posterior mean of the moment condition is itself a Monte Carlo average
#' over the draws, so it cannot be resolved below its own standard error;
#' asking for less than that is asking for noise.
#'
#' @keywords internal
#' @noRd
moment_tolerance <- function(BB) {
  max(1e-8, stats::sd(BB) / sqrt(length(BB)))
}


#' Couple the outcome and propensity score posteriors
#'
#' Tilts the joint posterior of the two models until the moment condition (3.4)
#' holds in posterior mean, which is what makes the resulting g-computation
#' estimand doubly robust (Theorem 1). The moment condition is the inverse
#' probability weighted one of equation (3.4), or the subclassification one of
#' equation (F.2) when `control$moment` is `"subclass"`.
#'
#' @param betas.ps,betas.otc Draws by parameters matrices of coefficients.
#' @param d The object returned by [tilting_data()].
#' @param control A [drbayes_control()] list.
#' @param tol Tolerance on `|mean B_n|`, or `NULL` to use [moment_tolerance()].
#' @param method `"smc"` for the sequential Monte Carlo sweep of Algorithm 2,
#'   `"is"` for the importance sampling of Algorithm 1.
#'
#' @return A list with the tilted draws and the diagnostics recorded in
#'   `fit$smc`:
#'   \describe{
#'     \item{`lambda`}{The tilting parameter reached.}
#'     \item{`B_mean`}{The posterior mean of the moment condition at the draws
#'       that are being returned, so that the number reported describes the
#'       object the caller receives.}
#'     \item{`B_mean_weighted`}{The importance weighted mean of equation (3.9)
#'       that Algorithm 1 solves for, which the resampled draws approximate
#'       with a further Monte Carlo error of their own. `NA` for the sweep,
#'       which resamples to equal weights at every step and so has no
#'       importance weights at the end.}
#'     \item{`tol`}{The tolerance `B_mean` is judged against.}
#'     \item{`n_steps`, `max_steps`}{Steps taken, and the cap they were taken
#'       against: `control$n_steps` for the sweep, `control$newton_steps` for
#'       Algorithm 1.}
#'     \item{`converged`}{Whether the constraint was met on draws that are
#'       still numerous enough to describe a posterior, which needs
#'       `abs(B_mean)` below `tol` and no weight degeneracy.}
#'     \item{`ess`}{Effective sample size of the importance weights of
#'       Algorithm 1. `NA` for the sweep.}
#'     \item{`n_ancestors_min`}{Smallest number of distinct particles to
#'       survive a resampling step of the sweep. `NA` for Algorithm 1, which
#'       does not sweep.}
#'     \item{`n_distinct`}{Number of distinct draws among those Algorithm 1
#'       returns. `NA` for the sweep, whose kernel moves every particle, so
#'       that its returned particles are all distinct by construction.}
#'   }
#'   When the constraint already held and nothing was reweighted, `ess` and
#'   `n_distinct` are the number of draws and there is no sweep to report on.
#'
#' @keywords internal
#' @noRd
run_tilting <- function(betas.ps, betas.otc, d, control, tol = NULL,
                        method = c("smc", "is")) {
  method <- match.arg(method)
  S      <- nrow(betas.ps)

  moment_fn <- moment_ipw
  if (identical(control$moment, "subclass")) {
    d$subclass_weight <- subclass_weights(betas.ps, d, control$n_subclass)
    moment_fn         <- moment_subclass
  }

  BB <- moment_fn(betas.ps, betas.otc, d)
  if (is.null(tol)) tol <- moment_tolerance(BB)

  # Only the sign of the Newton step is used, to decide which way along the
  # lambda grid the constraint is approached.
  lam.new <- newton_lambda(BB, control$ridge)

  # Lemma 1: with a correctly specified outcome model the constraint already
  # holds at lambda = 0 and the tilted posterior is the original one, so the
  # sweep is skipped and the untilted draws are returned unchanged. Nothing is
  # reweighted or resampled, so all S draws remain in play.
  out <- list(betas.ps = betas.ps, betas.otc = betas.otc,
              lambda = 0, n_steps = 0L, BB = BB,
              ess = as.numeric(S), n_distinct = S)

  if (abs(mean(BB)) >= tol && lam.new != 0 && is.finite(lam.new)) {
    out <- if (method == "is") {
      tilt_is(betas.ps, betas.otc, BB, control, tol)
    } else {
      tilt_smc(betas.ps, betas.otc, BB, d, control, tol, sign(lam.new),
               moment_fn)
    }
  }

  list(betas.ps  = out$betas.ps,
       betas.otc = out$betas.otc,
       smc = tilting_report(out, S, control, tol,
                            max_steps = if (method == "is") {
                              control$newton_steps
                            } else {
                              control$n_steps
                            }))
}


#' Judge a finished tilting and say what is wrong with it
#'
#' Shared by [run_tilting()] and [run_tilting_np()] so that the two report the
#' same quantities and refuse the same failures.
#'
#' @param out What [tilt_is()], [tilt_is_index()] or [tilt_smc()] returned, or
#'   the untilted draws when the constraint already held.
#' @param S Number of draws the tilting started from.
#' @param control A [drbayes_control()] list.
#' @param tol Tolerance the posterior mean of the moment condition is judged
#'   against.
#' @param max_steps The cap the steps were taken against.
#'
#' @return The `smc` list documented in [run_tilting()].
#'
#' @keywords internal
#' @noRd
tilting_report <- function(out, S, control, tol, max_steps,
                           nonparametric = FALSE) {
  # Judged at the draws that are handed back, rather than at the weights that
  # produced them, so that the number reported is the number describing the
  # object the caller receives.
  B_mean          <- mean(out$BB)
  ess             <- out$ess %||% NA_real_
  n_ancestors_min <- out$n_ancestors_min %||% NA_integer_
  n_distinct      <- out$n_distinct %||% NA_integer_

  if (isTRUE(out$at_lambda_max)) {
    warning("Algorithm 1 reached the bound lambda_max = ",
            format(control$lambda_max, digits = 3), " without solving ",
            "equation (3.9), so the tilting parameter returned is that bound ",
            "and not a solution. A tilting parameter this large is the ",
            "violation of condition (C.3) that Section 5.3.1 is about, so ",
            "raising lambda_max mostly buys a more degenerate answer. ",
            if (nonparametric) {
              paste0("Algorithm 2 is not available for draws of fitted means, ",
                     "so what is left is to bring the constraint nearer to ",
                     "zero before tilting: fit the propensity score model ",
                     "better, or supply outcome draws that fit the data more ",
                     "closely.")
            } else {
              paste0("Prefer method = \"smc\", with ",
                     "drbayes_control(pruning = TRUE) if the constraint is ",
                     "still not met.")
            })
  }

  degeneracy <- tilting_degeneracy(ess, n_ancestors_min, S, control$ess_frac,
                                   nonparametric = nonparametric)
  if (degeneracy$degenerate) {
    warning(degeneracy$message)
  }
  if (abs(B_mean) >= tol) {
    warning("Posterior coupling did not satisfy the moment condition: ",
            "|mean B_n| = ", format(abs(B_mean), digits = 3), " after ",
            out$n_steps, " steps up to lambda = ",
            format(out$lambda, digits = 3),
            ". The returned 'pc' draws are not doubly robust. Inspect ",
            "$smc, and consider a wider lambda range or a better specified ",
            "propensity score model.")
  }

  list(lambda    = out$lambda,
       B_mean    = B_mean,
       B_mean_weighted = out$B_mean_weighted %||% NA_real_,
       tol       = tol,
       n_steps   = out$n_steps,
       max_steps = max_steps,
       converged = abs(B_mean) < tol && !degeneracy$degenerate,
       ess       = ess,
       n_ancestors_min = n_ancestors_min,
       n_distinct      = n_distinct)
}


#' Couple a nonparametric outcome posterior to the propensity score posterior
#'
#' The counterpart of [run_tilting()] for an outcome model that reports fitted
#' means rather than coefficients, which is what Section 5.4 needs in order to
#' put BART on the outcome side. The moment condition is [moment_np()] and the
#' tilted posterior is described by an index into the draws: entry \eqn{s} of
#' that index names a propensity score draw and, with it, the row of fitted
#' means the same draw carries.
#'
#' Only Algorithm 1 is available here. Algorithm 2 rejuvenates its particles
#' with the Gaussian kernel of Liu and West (2001), which smooths a parameter
#' vector, and a draw of fitted values has no parameter vector to smooth;
#' jittering the fitted values themselves would sample from a posterior nobody
#' specified. Resampling an index instead leaves every draw exactly as its
#' outcome model produced it.
#'
#' The price of dropping that kernel is that nothing replenishes the cloud, so
#' the returned draws are duplicates of draws that were already there. The
#' effective sample size of the importance weights is what measures how few
#' distinct ones are left, and the `control$ess_frac` floor applied by
#' [tilting_report()] is what stops a spent sample being reported as converged.
#'
#' @param betas.ps Draws by parameters matrix of propensity score coefficients.
#' @param d The object returned by [tilting_data()], carrying `outcome.mu`.
#' @param control A [drbayes_control()] list.
#' @param tol Tolerance on `|mean B_n|`, or `NULL` to use [moment_tolerance()].
#'
#' @return A list with the resampled draw index `idx`, the propensity score
#'   draws `betas.ps` it selects, and the diagnostics `smc` of [run_tilting()].
#'
#' @references
#' Liu, J. and West, M. (2001). Combined parameter and state estimation in
#' simulation-based filtering. In Sequential Monte Carlo Methods in Practice,
#' 197-223. Springer.
#'
#' @keywords internal
#' @noRd
run_tilting_np <- function(betas.ps, d, control, tol = NULL) {
  S <- nrow(betas.ps)

  if (identical(control$moment, "subclass")) {
    d$subclass_weight <- subclass_weights(betas.ps, d, control$n_subclass)
  }

  BB <- moment_np(betas.ps, d)
  if (is.null(tol)) tol <- moment_tolerance(BB)
  lam.new <- newton_lambda(BB, control$ridge)

  # Lemma 1, as in run_tilting(): when the constraint already holds the tilted
  # posterior is the original one, and every draw stays in play.
  out <- list(idx = seq_len(S), lambda = 0, n_steps = 0L, BB = BB,
              ess = as.numeric(S), n_distinct = S)

  if (abs(mean(BB)) >= tol && lam.new != 0 && is.finite(lam.new)) {
    # The same bracket [tilt_is_index()] insists on, refused here first because
    # the way out of it that Section 3.3 recommends, Algorithm 2, is the one
    # thing this path cannot offer.
    if (min(BB) >= 0 || max(BB) <= 0) {
      stop("The moment condition has the same sign at all ", S, " posterior ",
           "draws (from ", format(min(BB), digits = 3), " to ",
           format(max(BB), digits = 3), "), so no weighting of them averages ",
           "to zero and equation (3.9) has no solution at any tilting ",
           "parameter. Section 3.3 answers that with Algorithm 2, which moves ",
           "the particles rather than reweighting them, and that needs a ",
           "parameter vector the supplied fitted means do not have. Either ",
           "fit the outcome or the propensity score model so that the two ",
           "sit closer together, or express the outcome model in ",
           "coefficients and pass them as outcome.samples with ",
           "method = \"smc\".")
    }
    out <- tilt_is_index(BB, control, tol)
  }

  list(idx      = out$idx,
       betas.ps = betas.ps[out$idx, , drop = FALSE],
       smc      = tilting_report(out, S, control, tol,
                                 control$newton_steps,
                                 nonparametric = TRUE))
}


#' Whether the tilted posterior is still resting on enough draws to be one
#'
#' Meeting the moment condition is not on its own evidence that the tilting
#' learned anything. The importance weighted mean of \eqn{B_n} tends to
#' \eqn{\min_s B_s} as the tilting parameter runs off to minus infinity, so
#' weights collapsed onto the single smallest draw can sit well inside the
#' tolerance while describing a posterior of one point, whose standard
#' deviation is zero and whose credible interval has zero width. The two
#' methods leave different traces of that collapse, so whichever of the two
#' quantities the method measured is the one the floor applies to.
#'
#' @param ess Effective sample size of the importance weights, or `NA`.
#' @param n_ancestors_min Smallest number of distinct particles to survive a
#'   resampling step of the sweep, or `NA`.
#' @param S Number of draws the tilting started from.
#' @param ess_frac The floor, as a fraction of `S`.
#'
#' @return A list with `degenerate` and, when it is `TRUE`, the `message`
#'   saying which quantity fell below the floor and what to do about it.
#'
#' @keywords internal
#' @noRd
tilting_degeneracy <- function(ess, n_ancestors_min, S, ess_frac,
                               nonparametric = FALSE) {
  from_weights <- !is.na(ess)
  carried      <- if (from_weights) ess else as.numeric(n_ancestors_min)
  if (is.na(carried) || carried >= ess_frac * S) {
    return(list(degenerate = FALSE))
  }

  what <- if (from_weights) {
    paste0("an effective sample size of ", format(carried, digits = 3),
           " out of ", S, " draws")
  } else {
    paste0("only ", carried, " of the ", S, " particles distinct at the ",
           "worst resampling step of the sweep")
  }
  advice <- if (from_weights && nonparametric) {
    paste0("This is the weight degeneracy described in Section 3.3, the ",
           "drawback of Algorithm 1. Algorithm 2 would rejuvenate the ",
           "particles between steps, but it needs a parameter vector to ",
           "smooth and draws of fitted means do not have one, so it is not ",
           "available here. Supply more draws, or bring the moment condition ",
           "nearer to zero before tilting by fitting either model better.")
  } else if (from_weights) {
    paste0("This is the weight degeneracy described in Section 3.3, which is ",
           "the drawback of Algorithm 1: use method = \"smc\", which moves ",
           "the particles along a sequence of tilting parameters instead of ",
           "reweighting the original draws once.")
  } else {
    paste0("The sweep is taking steps in the tilting parameter that the ",
           "particle cloud cannot follow: raise n_steps so that each step ",
           "reweights less.")
  }

  list(degenerate = TRUE,
       message = paste0(
         "The tilting left ", what, ", below the floor of ess_frac = ",
         format(ess_frac), " that drbayes_control() sets, so the tilted ",
         "posterior rests on a handful of draws and its spread, and every ",
         "credible interval built from it, is not trustworthy. It is ",
         "reported as not converged for that reason whatever the moment ",
         "condition says. ", advice))
}


#' Algorithm 1: importance sampling for the tilting parameter
#'
#' Equation (3.9) asks for the tilting parameter at which the draws, weighted
#' by \eqn{\exp(\lambda B_n)}, have weighted mean \eqn{B_n} equal to zero.
#' Algorithm 1 solves it by Newton iteration, dividing that weighted mean by
#' the weighted second moment of \eqn{B_n}. The draws are then reweighted at
#' the solution and resampled once to equal weight.
#'
#' Two things bound the search. The tilting parameter is held inside
#' `control$lambda_max`, because a large \eqn{\lambda} is the violation of
#' condition (C.3) that Section 5.3.1 is about rather than an answer to it, and
#' the caller is told when the bound is what stopped the iteration. And the
#' iteration is capped at `control$newton_steps`.
#'
#' Section 3.3 warns that these weights degenerate when the untilted posteriors
#' have little mass where the constraint holds, and the stopping rule cannot
#' see that happening: as \eqn{\lambda} runs off to minus infinity the weighted
#' mean of \eqn{B_n} tends to \eqn{\min_s B_s}, so weights collapsed onto the
#' single smallest draw satisfy the tolerance whenever that smallest draw is
#' itself small. The effective sample size \eqn{1 / \sum_s w_s^2} is the
#' quantity that does see it, so it is returned and [run_tilting()] holds it to
#' a floor before reporting convergence.
#'
#' Unlike [tilt_smc()] this takes no data: nothing it does asks for the moment
#' condition to be recomputed.
#'
#' @param betas.ps,betas.otc Draws by parameters matrices of coefficients.
#' @param BB The moment condition at those draws.
#' @param control A [drbayes_control()] list.
#' @param tol Tolerance on the weighted mean of the moment condition.
#'
#' @return What [tilt_is_index()] returns, with the coefficient draws it
#'   selects added as `betas.ps` and `betas.otc`.
#'
#' @keywords internal
#' @noRd
tilt_is <- function(betas.ps, betas.otc, BB, control, tol) {
  out <- tilt_is_index(BB, control, tol)
  c(list(betas.ps  = betas.ps[out$idx, , drop = FALSE],
         betas.otc = betas.otc[out$idx, , drop = FALSE]),
    out)
}


#' Algorithm 1 as a resampling of the draw index
#'
#' Everything Algorithm 1 does is a function of the moment condition alone, so
#' the draws it acts on need not be coefficients: [run_tilting_np()] tilts an
#' index into draws whose outcome side is a row of fitted means. The index is
#' returned rather than applied, which is what lets that path resample draws
#' held in matrices too large to copy.
#'
#' @param BB The moment condition at the draws.
#' @param control A [drbayes_control()] list.
#' @param tol Tolerance on the weighted mean of the moment condition.
#'
#' @return As [tilt_smc()], with `n_distinct` in place of `n_ancestors_min`,
#'   plus `idx`, the resampled draw index; `ess`; `B_mean_weighted`, the
#'   weighted mean of equation (3.9) at the chosen lambda; and `at_lambda_max`,
#'   TRUE when the bound is what stopped an iteration that had not reached a
#'   solution. The index comes back in ascending order of `BB`, which is the
#'   order [systematic_resample()] visits the draws in; they are equally
#'   weighted draws from the tilted posterior, so their order carries no
#'   information.
#'
#' @keywords internal
#' @noRd
tilt_is_index <- function(BB, control, tol) {
  S <- length(BB)

  # The weighted mean of B_n rises from min(B_n) to max(B_n) as lambda sweeps
  # the real line, so equation (3.9) has a root only when the draws straddle
  # zero. Without one the Newton iteration walks off towards infinity.
  if (min(BB) >= 0 || max(BB) <= 0) {
    stop("The moment condition has the same sign at all ", S, " posterior ",
         "draws (from ", format(min(BB), digits = 3), " to ",
         format(max(BB), digits = 3), "), so no weighting of them averages to ",
         "zero and equation (3.9) has no solution at any tilting parameter. ",
         "Algorithm 1 cannot proceed here: reweighting draws cannot move them ",
         "to where the constraint holds, and Section 3.3 gives exactly this ",
         "as its drawback and as the reason Algorithm 2 exists. Use ",
         "method = \"smc\", which moves the particles along a sequence of ",
         "tilting parameters, with drbayes_control(pruning = TRUE) if the ",
         "constraint is still not met.")
  }

  lambda  <- 0
  n.steps <- 0L

  for (step in seq_len(control$newton_steps)) {
    ww     <- normalised_weights(lambda * BB)
    B_mean <- sum(ww * BB)
    if (abs(B_mean) < tol) break

    # The step of Algorithm 1: the weighted mean of B_n over its weighted
    # second moment. Ridged on the same scale as newton_lambda(), both
    # denominators being an average of B_n^2 over the draws, so that a
    # weighted second moment of zero is not a division by zero.
    newton   <- -B_mean / (sum(ww * BB^2) + control$ridge)
    proposal <- damped_newton(lambda, newton, BB, abs(B_mean),
                              control$lambda_max)
    if (proposal == lambda) break
    lambda  <- proposal
    n.steps <- step
  }

  # Unless the loop converged it leaves lambda one update ahead of the weights,
  # so the weights are formed once more at the value actually returned.
  ww       <- normalised_weights(lambda * BB)
  weighted <- sum(ww * BB)

  # Resampling duplicates draws rather than moving them, so the moment
  # condition at a resampled draw is the one at its source and follows the
  # index rather than having to be evaluated again.
  idx <- systematic_resample(ww, BB)

  # The bound is worth reporting only when the iteration was still short of a
  # solution when it met it, rather than merely finishing there.
  at_bound <- abs(lambda) >= control$lambda_max && abs(weighted) >= tol

  list(idx       = idx,
       lambda    = lambda,
       n_steps   = n.steps,
       BB        = BB[idx],
       B_mean_weighted = weighted,
       ess       = 1 / sum(ww^2),
       n_distinct = length(unique(idx)),
       at_lambda_max = at_bound)
}


#' Resample to equal weight without spending the tolerance on the resampling
#'
#' Drawing \eqn{S} indices from a multinomial distribution, as step 2 of
#' Algorithm 2 does, leaves the mean of \eqn{B_n} over the resampled draws a
#' Monte Carlo standard error away from the weighted mean the tilting solved
#' for. That standard error is the same size as [moment_tolerance()], so draws
#' resampled from weights that satisfy equation (3.9) exactly would still fail
#' the convergence test about a third of the time, for no reason but the
#' resampling.
#'
#' Systematic resampling takes one evenly spaced sweep through the cumulative
#' weights from a single uniform offset instead, which pins the number of
#' copies of draw \eqn{s} to within one of \eqn{S w_s} while keeping
#' \eqn{S w_s} as its expectation, so the resampled cloud still targets the
#' tilted posterior. Visiting the draws in ascending order of \eqn{B_n} is what
#' turns that into a small error for \eqn{B_n} in particular: neighbouring
#' strata then hold draws with similar \eqn{B_n}, and the rounding they
#' contribute cancels rather than accumulates.
#'
#' @param ww Normalised weights, one per draw.
#' @param BB The moment condition at those draws, used for the ordering only.
#'
#' @return An integer vector of `length(ww)` indices into the draws, in
#'   ascending order of `BB`.
#'
#' @references
#' Kitagawa, G. (1996). Monte Carlo filter and smoother for non-Gaussian
#' nonlinear state space models. Journal of Computational and Graphical
#' Statistics, 5(1), 1-25.
#'
#' @keywords internal
#' @noRd
systematic_resample <- function(ww, BB) {
  S     <- length(ww)
  ranks <- order(BB)
  edges <- cumsum(ww[ranks])

  # Rounding can leave the last cumulative weight a hair under 1, which would
  # put an offset near 1 past the final edge, so the index is capped.
  offsets <- (stats::runif(1) + seq_len(S) - 1) / S
  ranks[pmin(findInterval(offsets, edges) + 1L, S)]
}


#' A Newton step shortened until it improves the weighted mean
#'
#' The weighted mean of \eqn{B_n} flattens out far from the root, where nearly
#' all the weight sits on one draw, so a full Newton step taken from there can
#' overshoot to a point worse than the one it started at and the iteration can
#' then oscillate away. Halving the step until the weighted mean shrinks keeps
#' the iteration inside the bracket the draws provide; the weighted mean is
#' increasing in lambda, so a short enough step in the Newton direction always
#' improves it.
#'
#' What the halving cannot do is notice weight collapse, because collapse is
#' one of the ways the weighted mean shrinks: a step that puts all the weight
#' on the smallest draw passes this test at the first attempt. `bound` is what
#' stops such a step from running away, and the effective sample size floor of
#' [run_tilting()] is what stops the result being called an answer.
#'
#' @param lambda The current tilting parameter.
#' @param step The full Newton step.
#' @param BB The moment condition at the draws.
#' @param best The absolute weighted mean the step has to beat.
#' @param bound Largest tilting parameter in absolute value that the step may
#'   reach, `control$lambda_max`.
#' @param max_halvings Most times the step may be halved.
#'
#' @return The updated tilting parameter, inside `bound`. It is returned
#'   unchanged when the iteration is already at the bound and still pushing
#'   outwards, which is how the caller learns to stop.
#'
#' @keywords internal
#' @noRd
damped_newton <- function(lambda, step, BB, best, bound, max_halvings = 10L) {
  clamp <- function(value) max(-bound, min(bound, value))

  for (halving in 0:max_halvings) {
    candidate <- clamp(lambda + step / 2^halving)
    ww        <- normalised_weights(candidate * BB)
    if (abs(sum(ww * BB)) < best) return(candidate)
  }
  # No improvement even from the shortest step: the iteration is as close to
  # the root as these draws can resolve, and the smallest step is the safest.
  clamp(lambda + step / 2^max_halvings)
}


#' Algorithm 2: sequential Monte Carlo along a grid of tilting parameters
#'
#' @param betas.ps,betas.otc Draws by parameters matrices of coefficients.
#' @param BB The moment condition at those draws.
#' @param d The object returned by [tilting_data()].
#' @param control A [drbayes_control()] list.
#' @param tol Tolerance on the posterior mean of the moment condition.
#' @param direction The sign of the Newton step, fixing which half of the
#'   lambda grid is swept.
#' @param moment_fn The moment condition to hold, [moment_ipw()] or
#'   [moment_subclass()].
#'
#' @return A list with the tilted particles, the tilting parameter reached, the
#'   number of steps taken, the moment condition at the particles, and
#'   `n_ancestors_min`, the smallest number of distinct particles to survive a
#'   resampling step. Every particle is moved by the smoothing kernel after
#'   that resampling, so the particles returned are all distinct however few
#'   ancestors they have between them, and the count of ancestors is the only
#'   record of how much of the cloud the sweep spent.
#'
#' @keywords internal
#' @noRd
tilt_smc <- function(betas.ps, betas.otc, BB, d, control, tol, direction,
                     moment_fn = moment_ipw) {
  S     <- nrow(betas.ps)
  n.ps  <- ncol(betas.ps)
  lam.s <- direction * seq(0, control$lambda_max,
                           length.out = control$n_steps + 1L)

  particles <- cbind(betas.ps, betas.otc)
  aa        <- control$smoothing  # Liu and West (2001) kernel
  BB_mean   <- mean(BB)
  lam.diff  <- diff(lam.s)
  lambda    <- 0
  n.steps   <- 0L
  n.eff.min <- S

  for (ll in seq_along(lam.diff)) {
    # Incremental importance weights exp{(lambda_t - lambda_{t-1}) B_n}
    ww <- normalised_weights(lam.diff[ll] * BB)

    # Section 5.3.1 discards samples with small values of exp(lambda B_n), the
    # weight against the untilted posterior at the tilting parameter reached so
    # far. That is not the incremental weight: an increment of lambda_max /
    # n_steps leaves every incremental weight within a fraction of a percent of
    # 1 / S, so an absolute threshold on those would discard nothing.
    if (isTRUE(control$pruning)) {
      cumulative <- normalised_weights(lam.s[ll + 1L] * BB)
      ww <- prune_weights(ww, cumulative, control)
    }

    resample_idx <- sample.int(S, size = S, replace = TRUE, prob = ww)
    n.eff.min    <- min(n.eff.min, length(unique(resample_idx)))

    ep <- mvtnorm::rmvnorm(S, sigma = (1 - aa^2) * stats::cov(particles))
    particles <- sweep(aa * particles[resample_idx, , drop = FALSE] + ep,
                       2, (1 - aa) * colMeans(particles), "+")

    betas.ps  <- particles[, seq_len(n.ps), drop = FALSE]
    betas.otc <- particles[, -seq_len(n.ps), drop = FALSE]

    BB          <- moment_fn(betas.ps, betas.otc, d)
    BB_mean_new <- mean(BB)
    lambda      <- lam.s[ll + 1L]
    n.steps     <- ll

    # Algorithm 2 step 4 stops once the constraint is met. A sign change also
    # stops the sweep, because it brackets the root and going further would
    # walk away from it, but on its own it does not mean the constraint holds.
    if (abs(BB_mean_new) < tol || BB_mean_new * BB_mean < 0) break
  }

  list(betas.ps = betas.ps, betas.otc = betas.otc,
       lambda = lambda, n_steps = n.steps, BB = BB,
       n_ancestors_min = n.eff.min)
}


#' Discard the least important particles before resampling
#'
#' Section 5.3.1 discards the samples carrying small sampling weights, so that
#' the particles that survive sit closer to \eqn{B_n = 0} and the tilting
#' parameter does not have to grow large to hold the constraint. The paper does
#' not say how small a weight has to be to be discarded. Zeroing the lowest
#' The paper does not state the discard threshold, so both rules below are this
#' package's reading of section 5.3.1 rather than values taken from the paper.
#'
#' @param ww Normalised incremental weights, one per particle. These are what
#'   the resampling draws on, and what is returned.
#' @param cumulative Normalised weights at the tilting parameter reached so far,
#'   which is the quantity section 5.3.1 thresholds. They decide WHICH particles
#'   go; `ww` decides what the sweep then does with the survivors.
#' @param control A [drbayes_control()] list, supplying `prune.rule` and either
#'   `prune.w` or `prune.q`.
#'
#' @return `ww` with the discarded entries set to zero and the rest
#'   renormalised, or `ww` untouched when nothing meets the threshold.
#'
#' @keywords internal
#' @noRd
prune_weights <- function(ww, cumulative, control) {
  S <- length(ww)
  if (identical(control$prune.rule, "quantile")) {
    n_drop <- floor(control$prune.q * S)
    if (n_drop < 1) return(ww)
    drop <- order(cumulative)[seq_len(n_drop)]
    setting <- paste0("prune.q, which is ", format(control$prune.q))
  } else {
    # The threshold is a fraction of the weight a particle would carry if all
    # of them were equal, so it keeps its meaning when the number of draws
    # changes. A bare absolute number would tighten as S grows.
    drop <- which(cumulative < control$prune.w / S)
    setting <- paste0("prune.w, which is ", format(control$prune.w))
  }
  if (length(drop) == 0L) return(ww)

  ww[drop] <- 0
  n_left <- sum(ww > 0)
  if (n_left < 2L) {
    stop("Sample pruning left ", n_left, " particle(s) with a positive ",
         "weight, which is not enough to resample from. Lower ", setting,
         ", or set pruning = FALSE. Weights this concentrated also mean the ",
         "sweep is taking large steps in the tilting parameter, so a larger ",
         "n_steps will help.")
  }
  ww / sum(ww)
}


#' Self-normalised weights from log weights
#'
#' Subtracting the maximum before exponentiating is what keeps `exp(lambda B_n)`
#' representable: the weights are a product over n observations, so the log
#' weights grow with the sample size and overflow long before the normalised
#' weights would.
#'
#' @keywords internal
#' @noRd
normalised_weights <- function(log_w) {
  w <- exp(log_w - max(log_w))
  w / sum(w)
}
