#' Tuning Parameters for Posterior Coupling
#'
#' @description
#' Collects the settings that control the Markov chain Monte Carlo sampling, the
#' sequential Monte Carlo sweep over the tilting parameter, and the convergence
#' check performed by \code{\link{drbayes_pc}}. Pass the result as the
#' \code{control} argument, in the same way \code{\link[stats]{glm.control}} is
#' passed to \code{\link[stats]{glm}}. Only the settings you name have to be
#' given, everything else keeps its default.
#'
#' @param mc Integer. Number of MCMC iterations per chain, before burn-in and
#'   thinning. Ignored when posterior draws are supplied directly. Default 5000.
#' @param bn Integer. Number of burn-in iterations to discard. Applied to
#'   internally generated and externally supplied draws alike, so that feeding a
#'   sampler's own output back in reproduces the internal path exactly.
#'   Default 1000.
#' @param thin Integer. Thinning interval applied after the burn-in. Default 2.
#' @param chains Integer. Number of Markov chains. Default 4. R-hat compares
#'   chains against each other, so at least four are recommended.
#' @param init Starting values for the samplers, or NULL (default) to place them
#'   automatically at a quick frequentist fit and disperse the remaining chains
#'   around it. A list of \code{chains} numeric vectors gives both models the
#'   same starting values, which is only meaningful when they have the same
#'   number of coefficients. A list with components \code{outcome} and \code{ps},
#'   each itself a list of \code{chains} numeric vectors, gives each model its
#'   own. The length of each vector is checked when the sampler runs, since the
#'   number of coefficients is not known until the design matrices have been
#'   built.
#' @param lambda_max Numeric. Largest value of the tilting parameter the sweep
#'   will consider, written \eqn{\bar{\lambda}}{lambda-bar} in the paper.
#'   Default 10.
#' @param n_steps Integer. Number of steps the sweep takes to travel from zero to
#'   \code{lambda_max}, written \eqn{T} in the paper. The increment between
#'   successive values of the tilting parameter is \code{lambda_max / n_steps}.
#'   Default 1000, giving the increment of 0.01 used in the paper. This is the
#'   resolution of the grid of Algorithm 2 only; the iteration cap of
#'   Algorithm 1 is \code{newton_steps}.
#' @param newton_steps Integer. Largest number of Newton iterations Algorithm 1
#'   may take when solving equation (3.9) for the tilting parameter, which is
#'   the work \code{drbayes_pc(method = "is")} does in place of the sweep.
#'   Default 100. The iteration is one dimensional and reaches the solution in
#'   a handful of steps whenever one exists, so the cap is there to end a
#'   search that is not converging rather than to be reached.
#' @param smoothing Numeric. Kernel smoothing coefficient of Liu and West (2001),
#'   written \eqn{a} in the paper. A number in (0, 1]. Default 0.99.
#' @param tol Numeric. Tolerance on the absolute posterior mean of the moment
#'   condition, which is the stopping rule for the sweep. NULL, the default, sets
#'   it to the Monte Carlo standard error of that mean over the retained draws:
#'   the mean is itself an average over a finite number of draws, so it cannot be
#'   resolved below its own sampling error, and asking for less than that is
#'   asking the sweep to chase noise. Give a positive number to fix the tolerance
#'   instead.
#' @param ridge Numeric. Small value added to the denominator of the Newton step
#'   for the tilting parameter, both the step Algorithm 1 iterates and the one
#'   that chooses the initial direction of the sweep. It keeps that step defined
#'   when every draw puts the moment condition at zero, which would otherwise be
#'   a division of zero by zero. Default 1e-5. Both denominators are an average
#'   of \eqn{B_n^2} over the draws, the same average at \eqn{\lambda = 0} and
#'   weighted by \eqn{\exp(\lambda B_n)} thereafter, so one value of
#'   \code{ridge} means the same thing in both. That average shrinks as the
#'   sample grows, since the moment condition itself does, and at a few hundred
#'   thousand observations the default is a noticeable share of it: that
#'   shortens the Newton steps and costs a few more iterations rather than
#'   moving the tilting parameter they converge on, and a smaller \code{ridge}
#'   removes the effect.
#' @param ess_frac Numeric in (0, 1]. Smallest share of the posterior draws that
#'   has to be carrying the tilted posterior for the coupling to be reported as
#'   converged. Default 0.1. Meeting the moment condition is not on its own
#'   evidence that anything was learned: the importance weighted mean of
#'   \eqn{B_n} tends to its smallest value as the tilting parameter runs off to
#'   minus infinity, so weights that have collapsed onto one draw can sit well
#'   inside \code{tol} while describing a posterior of a single point, whose
#'   standard deviation is zero and whose credible interval has zero width. A
#'   fit resting on fewer than \code{ess_frac} of the draws is therefore
#'   reported as not converged whatever the moment condition says. The quantity
#'   compared against it is the effective sample size of the importance weights
#'   under \code{method = "is"}, and the smallest number of distinct particles
#'   to survive a resampling step under the sweep, which resamples to equal
#'   weights and so has no importance weights of its own. A tenth is the usual
#'   rule of thumb for an importance sample that has been spent. Note that this
#'   is a different quantity from \code{ess_min}, which is about the Markov
#'   chains that produced the draws in the first place.
#' @param moment One of \code{"ipw"} (default) or \code{"subclass"}. Which form
#'   of the moment condition the tilting has to satisfy: the inverse probability
#'   weighted one of equation (3.4), or the propensity score subclassification
#'   one of equation (F.2). Subclassification replaces each individual weight by
#'   the treated fraction of the stratum the observation falls in, which is more
#'   stable when a few estimated propensity scores are close to 0 or 1.
#' @param n_subclass Integer. Number of propensity score strata used when
#'   \code{moment = "subclass"}, written \eqn{K} in Appendix F. Default 5,
#'   following Cochran (1968). The strata hold equal numbers of observations, so
#'   a large \code{n_subclass} leaves few in each and makes a stratum with no
#'   treated or no control units, which cannot be used, more likely.
#' @param pruning Logical. Whether to apply the sample pruning of Section 5.3.1,
#'   which discards the particles carrying the smallest weights at each step of
#'   the sweep so that the survivors concentrate where the moment condition
#'   holds. Default FALSE. It reduces the residual bias left when only the
#'   propensity score model is correctly specified, at the cost of discarding
#'   part of the particle cloud at every step, and it also ends the sweep much
#'   sooner: see the note on the pruning threshold below.
#' @param prune.rule Character. How pruning decides which particles to discard.
#'   \code{"weight"}, the default, discards every particle whose weight against
#'   the untilted posterior falls below \code{prune.w}, which is the rule the
#'   authors used. \code{"quantile"} instead discards a fixed fraction
#'   \code{prune.q} at every step, whatever the weights look like.
#' @param prune.w Numeric, positive. Weight below which a particle is discarded,
#'   as a multiple of \eqn{1/S}, the weight a particle would carry if all
#'   \eqn{S} of them were equal. Default 0.3, so a particle is discarded once it
#'   carries less than three tenths of an equal share. Below about 0.12 the
#'   threshold is never reached before the sweep stops, which makes
#'   \code{pruning = TRUE} a no-op; above about 0.5 it discards so much that the
#'   sweep ends after one or two steps. Expressing the threshold
#'   relative to \eqn{1/S} keeps its meaning when the number of draws changes.
#'   The weight tested is \eqn{\exp(\lambda B_n)} at the tilting parameter
#'   reached so far, as in section 5.3.1, not the incremental weight of a single
#'   step: the increments are all within a fraction of a percent of \eqn{1/S},
#'   so an absolute threshold on those would never discard anything. The paper
#'   does not state the value of the threshold, so the default is this package's
#'   choice; treat it as something to vary.
#' @param prune.q Numeric in [0, 1). Fraction of the particles discarded at each
#'   step, smallest weight first, when \code{prune.rule = "quantile"}. Default
#'   0.1. \code{prune.q = 0} leaves the sweep exactly as \code{pruning = FALSE}
#'   does. The particles it removes are those at one tail of the moment
#'   condition, so discarding a fixed fraction of them moves its posterior mean
#'   much further than an increment of the tilting parameter does, which ends
#'   the sweep after very few steps; the note below says how much further.
#' @param rhat_max Numeric. Largest acceptable rank-normalised split R-hat.
#'   Default 1.01, the threshold recommended by Vehtari et al. (2021).
#' @param ess_min Numeric. Smallest acceptable bulk and tail effective sample
#'   size. Default 400, again following Vehtari et al. (2021), who also note that
#'   below this the R-hat estimate is itself unreliable.
#' @param keep_draws Logical. Whether to return the per-chain posterior draws
#'   of the two models' coefficients, as the iterations by chains by parameters
#'   arrays the samplers produce. They are what [rank_plot()] and
#'   `plot(fit, type = "rank")` need, since a rank histogram is a statement
#'   about chains and the tilted particles have none. Default `FALSE`.
#' @param keep_particles Logical. Whether to return the tilted particles
#'   themselves. Default FALSE. See the note on memory below.
#'
#' @return An object of class \code{"drbayes_control"}, a named list with one
#'   element per argument above.
#'
#' @details
#' \strong{Why these are settings rather than constants:}
#'
#' \code{lambda_max}, \code{n_steps} and \code{smoothing} are the tuning
#' constants of Algorithm 2 of Orihara, Momozaki and Sugasawa (2025). The paper
#' leaves the choice of \eqn{\bar{\lambda}}{lambda-bar}, \eqn{T} and the smoothing
#' coefficient \eqn{a} to the user, so hard-coding them would hide exactly the
#' knobs the methodology asks you to tune. The defaults reproduce the settings
#' used in the paper and are a reasonable starting point, not a recommendation to
#' leave them alone.
#'
#' \code{newton_steps} and \code{ess_frac} have no counterpart in the paper,
#' which asks only that Algorithm 1 be iterated "until convergence" and does not
#' say what to do when it does not converge. They are this package's answer to
#' that: a cap on the search, and a floor below which the sample the search
#' arrived at is too thin to be called a posterior. \code{lambda_max} bounds
#' Algorithm 1 as well as the sweep, for the same reason.
#'
#' \strong{Choosing the smoothing coefficient:}
#'
#' The sweep resamples the particles at every step and then applies the Liu and
#' West (2001) kernel, which shrinks each particle towards the cloud mean and
#' adds normal noise carrying the variance that shrinking removed. One
#' application is close to harmless. Applying it \code{n_steps} times is not:
#' repeated shrink-and-jitter pulls the particle cloud towards a multivariate
#' normal with the right first two moments and the wrong shape. That is a poor
#' fit for the heavy-tailed posterior with a spike at zero that the horseshoe
#' priors produce, so a user fitting \code{\link{bayes_lm_hs}},
#' \code{\link{bayes_logit_hs}} or \code{\link{bayes_probit_hs}} may want fewer
#' and larger steps in the tilting parameter, for example \code{n_steps = 200},
#' or a smoothing coefficient closer to 1, which shrinks less at each
#' application.
#'
#' \strong{The pruning threshold:}
#'
#' Section 5.3.1 of the paper describes sample pruning as discarding the samples
#' with small sampling weights, and stops there: it gives no threshold, and no
#' rule for deriving one. \code{prune.q} is this package's reading of that
#' sentence, discarding a fixed fraction of the particles at each step, and its
#' default of 0.1 is a choice made here rather than a value reported in the
#' paper.
#'
#' Pruning ends the sweep far sooner than its name suggests, and it is worth
#' being clear about why. The particles it discards are those at one tail of the
#' moment condition, the tail the tilting is pushing away from, so deleting
#' \code{prune.q} of the cloud moves the posterior mean of \eqn{B_n} by a
#' noticeable fraction of that mean's distance from zero, while one increment of
#' the tilting parameter moves it by far less: at the default settings one
#' pruning step is worth some tens of ordinary steps. Step 4 of Algorithm 2
#' evaluates the constraint on the particles it currently holds, which are the
#' pruned ones, so the sweep typically stops after one or two steps at a tilting
#' parameter close to zero. How many steps it takes is then set by
#' \code{prune.q} and hardly at all by \code{n_steps}.
#'
#' The direction of that is what Section 5.3.1 asks for, since pruning exists
#' precisely so that \eqn{\lambda} need not grow large enough to violate
#' condition (C.3) of Appendix A. What it costs is that the constraint has been
#' met largely by deleting particles rather than by tilting them: the
#' \code{lambda} reported in \code{$smc} is a poor measure of how much tilting
#' was done, and the returned cloud is a truncated one, so its credible
#' intervals are narrower than those of the tilted posterior by an amount that
#' reflects the deletion and not information gained. Fit both ways. A result
#' that changes materially between \code{pruning = FALSE} and
#' \code{pruning = TRUE}, or between \code{prune.q = 0.05} and
#' \code{prune.q = 0.2}, is telling you that the tilting is being carried by a
#' small part of the particle cloud, which is worth knowing.
#'
#' \strong{Keeping the particles:}
#'
#' \code{keep_particles = TRUE} returns the tilted particles, one draws by
#' parameters matrix for the outcome model and one for the propensity score
#' model. They are useful for inspecting what the tilting did to the
#' coefficients, and expensive in memory: with the default settings and a
#' moderate number of covariates they are several times the size of everything
#' else the fit returns.
#'
#' @references
#' Orihara, S., Momozaki, T., & Sugasawa, S. (2025). Bayesian Doubly Robust
#' Causal Inference via Posterior Coupling. arXiv preprint arXiv:2506.04868.
#'
#' Cochran, W. G. (1968). The effectiveness of adjustment by subclassification
#' in removing bias in observational studies. Biometrics, 24(2), 295-313.
#'
#' Liu, J., & West, M. (2001). Combined parameter and state estimation in
#' simulation-based filtering. In A. Doucet, N. de Freitas, & N. Gordon (Eds.),
#' Sequential Monte Carlo Methods in Practice (pp. 197-223). Springer.
#'
#' Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Buerkner, P.-C. (2021).
#' Rank-normalization, folding, and localization: an improved R-hat for assessing
#' convergence of MCMC. Bayesian Analysis, 16(2), 667-718.
#'
#' @seealso \code{\link{drbayes_pc}}, which consumes the result.
#'
#' @examples
#' # The defaults, ready to be modified
#' drbayes_control()
#'
#' # Horseshoe priors give a heavy-tailed posterior that repeated kernel
#' # smoothing distorts, so take fewer and larger steps in lambda.
#' ctrl <- drbayes_control(n_steps = 200)
#' ctrl$n_steps
#'
#' # A plain list of overrides is accepted wherever a control object is, so
#' # control = list(n_steps = 200) means the same thing.
#'
#' set.seed(1)
#' n <- 200
#' dat <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
#' dat$A <- rbinom(n, 1, plogis(0.2 + 0.5 * dat$X1 - 0.3 * dat$X2))
#' dat$Y <- 1 + 1.5 * dat$A + 0.8 * dat$X1 + 0.6 * dat$X2 + rnorm(n)
#'
#' fit <- drbayes_pc(
#'   outcome.formula = Y ~ A + X1 + X2,
#'   ps.formula      = A ~ X1 + X2,
#'   data            = dat,
#'   outcome.model   = bayes_lm_hs,
#'   ps.model        = bayes_logit,
#'   control         = drbayes_control(mc = 1000, bn = 200, n_steps = 60)
#' )
#' fit$smc$n_steps
#'
#' @export
drbayes_control <- function(mc = 5000, bn = 1000, thin = 2, chains = 4L,
                            init = NULL,
                            lambda_max = 10, n_steps = 1000,
                            newton_steps = 100L, smoothing = 0.99,
                            tol = NULL, ridge = 1e-5, ess_frac = 0.1,
                            moment = c("ipw", "subclass"), n_subclass = 5L,
                            pruning = FALSE,
                            prune.rule = c("weight", "quantile"),
                            prune.w = 0.3, prune.q = 0.1,
                            rhat_max = 1.01, ess_min = 400,
                            keep_draws = FALSE,
                            keep_particles = FALSE) {

  mc   <- check_count(mc, "mc", min = 1L,
                      what = "the number of MCMC iterations per chain")
  bn   <- check_count(bn, "bn", min = 0L,
                      what = "the number of burn-in iterations to discard")
  thin <- check_count(thin, "thin", min = 1L,
                      what = "the thinning interval")

  if (mc <= bn + 1L) {
    stop("mc must exceed bn + 1, otherwise no posterior draw survives the ",
         "burn-in (mc = ", mc, ", bn = ", bn, ")")
  }
  n_kept <- length(seq(bn + 1L, mc, by = thin))
  if (n_kept < 2L) {
    stop("mc = ", mc, ", bn = ", bn, " and thin = ", thin, " leave only ",
         n_kept, " draw per chain. Raise mc, or lower bn or thin.")
  }

  chains <- validate_chains(chains)
  init   <- check_control_init(init, chains)

  # lambda_max and n_steps are the range and the resolution of one grid, so a
  # user who changes either has changed the step size the sweep actually takes.
  if (!is.numeric(lambda_max) || length(lambda_max) != 1L ||
      !is.finite(lambda_max) || lambda_max <= 0) {
    stop("lambda_max must be a single positive number, the largest value of ",
         "the tilting parameter the sweep may reach. Received ",
         describe_value(lambda_max))
  }
  n_steps <- check_count(
    n_steps, "n_steps", min = 1L,
    what = "the number of steps the sweep takes to reach lambda_max")
  newton_steps <- check_count(
    newton_steps, "newton_steps", min = 1L,
    what = "the largest number of Newton iterations Algorithm 1 may take")

  if (!is.numeric(smoothing) || length(smoothing) != 1L ||
      !is.finite(smoothing) || smoothing <= 0 || smoothing > 1) {
    stop("smoothing must be a single number in (0, 1], the Liu and West (2001) ",
         "kernel coefficient. Values close to 1 disturb the particle cloud ",
         "least; 1 itself switches the kernel off and leaves plain resampling, ",
         "which lets the particles collapse onto a few distinct values.")
  }

  if (!is.null(tol)) {
    if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
      stop("tol must be NULL or a single positive number, the tolerance on the ",
           "absolute posterior mean of the moment condition. NULL uses that ",
           "mean's own Monte Carlo standard error.")
    }
  }

  if (!is.numeric(ridge) || length(ridge) != 1L || !is.finite(ridge) ||
      ridge <= 0) {
    stop("ridge must be a single positive number. It is added to the ",
         "denominator of the Newton step for the tilting parameter, so that ",
         "the step stays defined when every draw puts the moment condition at ",
         "zero. Use a smaller value such as 1e-10 rather than 0.")
  }

  if (!is.numeric(ess_frac) || length(ess_frac) != 1L ||
      !is.finite(ess_frac) || ess_frac <= 0 || ess_frac > 1) {
    stop("ess_frac must be a single number in (0, 1], the smallest share of ",
         "the posterior draws that has to be carrying the tilted posterior ",
         "for it to be reported as converged. Use a small value such as ",
         "0.01 rather than 0 to accept a nearly collapsed sample. Received ",
         describe_value(ess_frac))
  }

  if (!is.character(moment) || length(moment) == 0L ||
      !all(moment %in% c("ipw", "subclass"))) {
    stop("moment must be \"ipw\", for the inverse probability weighted moment ",
         "condition of equation (3.4), or \"subclass\", for the propensity ",
         "score subclassification one of Appendix F. Received ",
         describe_value(moment))
  }
  moment     <- moment[1L]
  n_subclass <- check_count(
    n_subclass, "n_subclass", min = 2L,
    what = "the number of propensity score strata, which cannot be below two")

  if (!is.logical(pruning) || length(pruning) != 1L || is.na(pruning)) {
    stop("pruning must be TRUE or FALSE, switching the sample pruning of ",
         "Section 5.3.1 on or off")
  }
  prune.rule <- match.arg(prune.rule)
  if (!is.numeric(prune.w) || length(prune.w) != 1L || !is.finite(prune.w) ||
      prune.w <= 0) {
    stop("prune.w must be a single positive number, the weight below which a ",
         "particle is discarded, as a multiple of 1/S. Received ",
         describe_value(prune.w))
  }
  if (!is.numeric(prune.q) || length(prune.q) != 1L || !is.finite(prune.q) ||
      prune.q < 0 || prune.q >= 1) {
    stop("prune.q must be a single number in [0, 1), the fraction of the ",
         "particles that pruning discards at each step. A value of 1 would ",
         "discard every particle; use 0 to leave the sweep untouched. ",
         "Received ", describe_value(prune.q))
  }

  if (!is.numeric(rhat_max) || length(rhat_max) != 1L || !is.finite(rhat_max) ||
      rhat_max <= 1) {
    stop("rhat_max must be a single number greater than 1. R-hat is at least 1 ",
         "by construction, so a threshold of 1 or less can never be met. ",
         "Vehtari et al. (2021) recommend 1.01.")
  }
  if (!is.numeric(ess_min) || length(ess_min) != 1L || !is.finite(ess_min) ||
      ess_min <= 0) {
    stop("ess_min must be a single positive number, the smallest acceptable ",
         "bulk and tail effective sample size. Vehtari et al. (2021) ",
         "recommend 400.")
  }

  if (!is.logical(keep_draws) || length(keep_draws) != 1L ||
      is.na(keep_draws)) {
    stop("keep_draws must be TRUE or FALSE")
  }
  if (!is.logical(keep_particles) || length(keep_particles) != 1L ||
      is.na(keep_particles)) {
    stop("keep_particles must be TRUE or FALSE")
  }

  structure(
    list(
      mc             = mc,
      bn             = bn,
      thin           = thin,
      chains         = chains,
      init           = init,
      lambda_max     = as.numeric(lambda_max),
      n_steps        = n_steps,
      newton_steps   = newton_steps,
      smoothing      = as.numeric(smoothing),
      tol            = if (is.null(tol)) NULL else as.numeric(tol),
      ridge          = as.numeric(ridge),
      ess_frac       = as.numeric(ess_frac),
      moment         = moment,
      n_subclass     = n_subclass,
      pruning        = pruning,
      prune.rule     = prune.rule,
      prune.w        = as.numeric(prune.w),
      prune.q        = as.numeric(prune.q),
      rhat_max       = as.numeric(rhat_max),
      ess_min        = as.numeric(ess_min),
      keep_draws = keep_draws,
      keep_particles = keep_particles
    ),
    class = "drbayes_control")
}


#' Print Posterior Coupling Control Settings
#'
#' @param x An object of class \code{"drbayes_control"}, as returned by
#'   \code{\link{drbayes_control}}.
#' @param ... Ignored, present for compatibility with the generic.
#'
#' @return \code{x}, invisibly. Called for the printed output.
#'
#' @examples
#' print(drbayes_control(n_steps = 200))
#'
#' @export
print.drbayes_control <- function(x, ...) {
  num <- function(value) format(value, digits = 4, trim = TRUE)

  per_model <- intersect(names(x$init), c("outcome", "ps"))
  init_text <- if (is.null(x$init)) {
    "automatic"
  } else if (length(per_model) == 2L) {
    "supplied per model"
  } else if (length(per_model) == 1L) {
    paste0("supplied for the ", per_model, " model")
  } else {
    "supplied"
  }
  tol_text <- if (is.null(x$tol)) "Monte Carlo standard error" else num(x$tol)
  moment_text <- if (identical(x$moment, "subclass")) {
    paste0("subclass (n_subclass = ", num(x$n_subclass), ")")
  } else {
    "ipw"
  }
  prune_text <- if (!isTRUE(x$pruning)) {
    "off"
  } else if (identical(x$prune.rule, "quantile")) {
    paste0("on (quantile, prune.q = ", num(x$prune.q), ")")
  } else {
    paste0("on (weight, prune.w = ", num(x$prune.w), "/S)")
  }

  cat("Control settings for drbayes_pc()\n\n")
  cat("  MCMC         mc = ", num(x$mc), ", bn = ", num(x$bn),
      ", thin = ", num(x$thin), ", chains = ", num(x$chains),
      ", init = ", init_text, "\n", sep = "")
  cat("  Tilting      lambda_max = ", num(x$lambda_max),
      ", n_steps = ", num(x$n_steps),
      " (step ", num(x$lambda_max / x$n_steps), ")",
      ", smoothing = ", num(x$smoothing), "\n", sep = "")
  cat("               tol = ", tol_text, ", ridge = ", num(x$ridge),
      ", ess_frac = ", num(x$ess_frac), "\n", sep = "")
  cat("               moment = ", moment_text, ", pruning = ", prune_text,
      "\n", sep = "")
  cat("               newton_steps = ", num(x$newton_steps),
      " (importance sampling only)\n", sep = "")
  cat("  Diagnostics  rhat_max = ", num(x$rhat_max),
      ", ess_min = ", num(x$ess_min), "\n", sep = "")
  cat("  Output       keep_draws = ", x$keep_draws,
      ", keep_particles = ", x$keep_particles, "\n", sep = "")

  invisible(x)
}


#' Resolve whatever a caller passed as `control`
#'
#' Accepts a `drbayes_control` object, a plain list of overrides, or NULL, and
#' returns a complete validated control object in every case. A supplied object
#' is passed back through [drbayes_control()] rather than trusted, so that an
#' object whose elements were edited by hand after construction is checked
#' before its settings reach the sampler.
#'
#' @param control A `drbayes_control` object, a named list of overrides, or NULL.
#' @param arg Argument name, used in error messages.
#'
#' @return An object of class `"drbayes_control"`.
#'
#' @keywords internal
#' @noRd
resolve_control <- function(control = NULL, arg = "control") {
  if (is.null(control)) {
    return(drbayes_control())
  }

  settings <- unclass(control)
  if (!is.list(settings)) {
    stop(arg, " must be NULL, a list of settings such as list(mc = 2000), or ",
         "an object from drbayes_control(), but is of class ",
         class(control)[1L])
  }
  if (length(settings) == 0L) {
    return(drbayes_control())
  }

  names_given <- names(settings)
  unnamed <- if (is.null(names_given)) {
    length(settings)
  } else {
    sum(is.na(names_given) | !nzchar(names_given))
  }
  if (unnamed > 0L) {
    stop(arg, " must be a NAMED list, as in ", arg, " = list(mc = 2000), so ",
         "that each value can be matched to a setting. ", unnamed, " of its ",
         length(settings), " element(s) have no name.")
  }

  known   <- names(formals(drbayes_control))
  unknown <- setdiff(names_given, known)
  if (length(unknown) > 0) {
    stop(arg, " contains unknown setting(s): ", toString(unknown),
         ". The available settings are: ", toString(known), ".")
  }
  repeated <- unique(names_given[duplicated(names_given)])
  if (length(repeated) > 0) {
    stop(arg, " sets the same setting more than once: ", toString(repeated))
  }

  do.call(drbayes_control, settings)
}


#' Starting values a control object supplies for one of the two models
#'
#' Returns the element of `init` belonging to `model`, or the whole of `init`
#' when one set of starting values was given for both models, or NULL. Isolating
#' this keeps the two accepted shapes of `init` from being reinterpreted at each
#' call site.
#'
#' @param control A validated `drbayes_control` object.
#' @param model Either "outcome" or "ps".
#'
#' @return A list of starting values, one per chain, or NULL.
#'
#' @keywords internal
#' @noRd
control_init_for <- function(control, model = c("outcome", "ps")) {
  model <- match.arg(model)
  init  <- control$init
  if (is.null(init)) {
    return(NULL)
  }
  if (any(c("outcome", "ps") %in% names(init))) {
    return(init[[model]])
  }
  init
}


#' Validate a count-valued setting
#'
#' @param value The value supplied.
#' @param arg Argument name, used in the error message.
#' @param min Smallest permitted value.
#' @param what Plain English description of what the number counts.
#'
#' @return The value, as an integer.
#'
#' @keywords internal
#' @noRd
check_count <- function(value, arg, min = 1L, what = NULL) {
  # The bound is checked here rather than left to as.integer(), which turns
  # anything past it into NA with only a warning, so a setting that passed
  # every guard would come back missing.
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value != round(value) || value < min ||
      value > .Machine$integer.max) {
    stop(arg, " must be a single ",
         if (min >= 1L) "positive" else "non-negative", " integer",
         if (is.null(what)) "" else paste0(", ", what),
         ", and no larger than ", .Machine$integer.max,
         ". Received ", describe_value(value))
  }
  as.integer(value)
}


#' Describe a rejected value for an error message
#'
#' Shows the value itself when it is a single number, since that is what the
#' user needs to see, and falls back on the class and length when printing the
#' value would be unhelpful or long.
#'
#' @param value The value that failed validation.
#'
#' @keywords internal
#' @noRd
describe_value <- function(value) {
  if (is.atomic(value) && length(value) == 1L) {
    return(format(value, digits = 6))
  }
  paste0("an object of class ", class(value)[1L], " and length ", length(value))
}


#' Validate the starting values held in a control object
#'
#' How many coefficients each model has is not known here, since the design
#' matrices are built from the formulas later, so only the shape of `init` is
#' checked. The length of each vector is checked by the sampler.
#'
#' @param init The value supplied for `init`.
#' @param chains Number of chains, one starting value being needed per chain.
#'
#' @return `init`, unchanged, or NULL.
#'
#' @keywords internal
#' @noRd
check_control_init <- function(init, chains) {
  if (is.null(init)) {
    return(NULL)
  }
  if (!is.list(init)) {
    stop("init must be NULL or a list, but is of class ", class(init)[1L],
         ". Give one starting value per chain, as in ",
         "list(rep(0, p), rep(0, p)) for two chains, or name the models, as ",
         "in list(outcome = ..., ps = ...).")
  }

  models <- c("outcome", "ps")
  if (any(names(init) %in% models)) {
    if (!all(names(init) %in% models)) {
      stop("init names the models it applies to, so every element must be ",
           "called \"outcome\" or \"ps\". Found: ",
           toString(setdiff(names(init), models)), ".")
    }
    for (model in names(init)) {
      check_init_chains(init[[model]], chains, paste0("init$", model))
    }
    return(init)
  }

  check_init_chains(init, chains, "init")
  init
}


#' Validate one list of per-chain starting values
#'
#' @param init A list holding one numeric vector per chain.
#' @param chains Number of chains.
#' @param arg Argument name, used in error messages.
#'
#' @keywords internal
#' @noRd
check_init_chains <- function(init, chains, arg) {
  if (!is.list(init)) {
    stop(arg, " must be a list holding one starting value per chain, but is ",
         "of class ", class(init)[1L], ". A single starting value still has ",
         "to be wrapped in a list.")
  }
  if (length(init) != chains) {
    stop(arg, " has ", length(init), " element(s) but chains is ", chains,
         ", so ", chains, " starting values are needed, one per chain")
  }
  for (i in seq_len(chains)) {
    value <- init[[i]]
    if (!is.numeric(value) || length(value) == 0L) {
      stop(arg, "[[", i, "]] must be a non-empty numeric vector of starting ",
           "values, the intercept followed by one value per column of the ",
           "design matrix, but is of class ", class(value)[1L], " and length ",
           length(value))
    }
    if (!all(is.finite(value))) {
      stop(arg, "[[", i, "]] must be finite; found ", sum(!is.finite(value)),
           " missing or infinite value(s)")
    }
  }
  invisible(TRUE)
}
