#' Sensitivity analysis for unmeasured confounding
#'
#' @description
#' Turns the posterior draws of a [drbayes_pc()] fit into the posterior they
#' would have been had the outcome model carried an unmeasured confounding bias
#' \eqn{\xi}. The original draws are reweighted by the importance sampling of
#' Algorithm 3 of Orihara, Momozaki and Sugasawa (2025), or of Algorithm A.1 in
#' their Appendix E, and the reweighted draws are then coupled again by
#' Algorithm 2, which is what keeps the answer doubly robust. No model is
#' refitted.
#'
#' @details
#' \strong{What xi is.}
#'
#' Without strong ignorability the contrast the data identify is not the average
#' treatment effect but
#' \deqn{E[E[Y | A = 1, X] - E[Y | A = 0, X]] = \tau + \xi,}
#' so \eqn{\xi} is the part of the observed contrast that unmeasured confounding
#' contributes, and \eqn{\xi = 0} is strong ignorability (Daniels, Linero and
#' Roy, 2023, and Section 6.2 of the paper). It enters through the shifted
#' outcome model \eqn{m_A(X; \beta, \xi) = m_A(X; \beta) + A\xi}: only treated
#' observations are shifted, so only they contribute to the weights.
#'
#' \strong{The shift is on the conditional mean.}
#'
#' \eqn{m_a(X)} is \eqn{E[Y | A = a, X]}, so the shift is added to the fitted
#' mean and not to the linear predictor. For a gaussian outcome those are the
#' same thing. For a binary one they are not, and it is the mean that is meant:
#' \eqn{\xi} is a shift of the probability, which is what puts it in the units
#' of \eqn{\tau} in the display above and on the risk difference scale that
#' Section 6.3 reports. Added to the log odds instead it would move the average
#' treatment effect by several times less. On the binary example of this
#' package's tests, whose treated observations have a mean fitted risk of 0.17,
#' \eqn{E[\xi] = 0.013} moves the g-formula effect by -0.013 as a shift of the
#' mean and by -0.002 as a shift of the log odds, a factor of 8.
#'
#' A probability cannot leave \eqn{[0, 1]}, so where \eqn{m_A(X; \beta) + \xi}
#' would, the shifted mean is held at the boundary and that observation carries
#' a bias smaller than \eqn{\xi}. How many fitted values that happens to is
#' warned about and reported as `n_saturated`. It is a statement about the
#' outcome and not a numerical convenience: a treated group whose fitted risk
#' is 0.1 has no room for a risk difference of -0.5 of unmeasured confounding,
#' so the \eqn{Tri(-0.5, 0)} of Section 6.3 is a bias a rare binary outcome
#' cannot carry, and this reports how far short of it the analysis fell.
#'
#' A positive \eqn{\xi} means part of what was read as a treatment effect was
#' confounding, so the reweighted posterior of the average treatment effect is
#' smaller than the original one, by about \eqn{E[\xi]} while the reweighting
#' is still well resolved: on the gaussian example of this package's tests
#' \eqn{E[\xi] = 0.033} moves it by -0.033, and on the binary one
#' \eqn{E[\xi] = 0.013} moves it by -0.013. Beyond that the shift falls short,
#' because importance sampling cannot move a posterior further than the draws
#' it started from reach. And the coupling that follows the reweighting gives
#' most of even that back, for the reason set out two sections below.
#'
#' \strong{Which algorithm.}
#'
#' `method = "common"` is Algorithm A.1: one \eqn{\xi} is shared by every
#' observation, the sum over observations happens inside the integral, and the
#' whole uncertainty in \eqn{g()} is carried into the reweighted posterior,
#' which therefore comes out wider than the original by about the variance of
#' \eqn{\xi}. This is the default, and it is what the single sensitivity
#' parameter \eqn{\xi = \Delta} of Section 6.2 means.
#'
#' `method = "per-observation"` is Algorithm 3: each observation has its own
#' \eqn{\xi_i} drawn from \eqn{g()}, and each integral is done separately before
#' the product over observations. Independent biases average out over the
#' sample, so this shifts the posterior by \eqn{E[\xi]} while adding almost none
#' of the spread of \eqn{g()}; use it when the bias really does vary from one
#' observation to the next.
#'
#' \strong{Reweight first, then couple.}
#'
#' The weights of both algorithms are ratios of outcome model likelihoods, so
#' they belong to the untilted posterior \eqn{p_n(\alpha, \beta | D)} of
#' equation (3.3) and that is what they are applied to here. The reweighted
#' draws are then resampled to equal weight and coupled again, which is the
#' order Section 6.2 of the paper prescribes: the posterior samples under
#' unmeasured confounding are obtained first, and the doubly robust posterior of
#' the average treatment effect follows from them through Algorithm 1 or 2.
#'
#' Reweighting the coupled draws instead would leave the moment condition (3.4)
#' violated, because importance weights move the particle cloud and nothing
#' moves it back. On the gaussian example of this package's tests that costs a
#' factor of 14 to 55 on \eqn{|mean B_n|} relative to the tolerance the fit
#' itself was held to, over sensitivity distributions from \eqn{Tri(0, 0.05)}
#' to \eqn{Tri(0, 0.5)}; coupling the reweighted draws instead brings every one
#' of them back inside that tolerance. The recoupling is what the `$smc`
#' component of the result reports, and its `converged` entry is the verdict on
#' the analysis as a whole rather than on the coupling alone; see "Reading the
#' result". It is printed with the result.
#'
#' \strong{What the coupling does to the shift.}
#'
#' The moment condition (3.4) is written with the unshifted mean function
#' \eqn{m_A(X; \beta)} and the propensity score of the same fit, and both of
#' those describe the data as observed, unmeasured confounding included.
#' Imposing it on the reweighted draws therefore pulls the treated fitted
#' values back to what the propensity score weighted data say, and most of the
#' shift with them. Measured on the gaussian example, average treatment effect
#' 2.052:
#'
#' \tabular{lrrrr}{
#'   \strong{g()} \tab \strong{\eqn{E[\xi]}} \tab \strong{after reweighting} \tab
#'     \strong{after coupling} \tab \strong{effective draws} \cr
#'   Tri(0, 0.05) \tab 0.033 \tab -0.033 \tab -0.001 \tab 2355 of 3000 \cr
#'   Tri(0, 0.20) \tab 0.134 \tab -0.103 \tab -0.006 \tab 191 of 3000 \cr
#'   Tri(0, 0.50) \tab 0.334 \tab -0.120 \tab -0.019 \tab 50 of 3000
#' }
#'
#' Each entry is the mean of five runs, because none of these is a fixed
#' number. The sweep moves the fourth column by about 0.002 from one run to the
#' next in the first two rows and by 0.02 in the third, which is the column's
#' own message: what survives the coupling is of the size of the noise of the
#' sweep that produced it.
#'
#' Only the first row is an answer. The other two rest on far fewer effective
#' draws than the floor of [drbayes_control()], so both are reported as not
#' doubly robust and have their credible intervals withheld, whether or not
#' their coupling happens to meet the moment condition. It did in two of the
#' five runs each, which is the coin flip the next section is about. They are
#' tabulated to show what the coupling does to a shift, not as results.
#'
#' The spread of \eqn{g()} survives the coupling but the location shift largely
#' does not, which is why the result reports both rows: `ate` after the
#' coupling and `ate_gcomp` before it. This is not an artefact of this
#' implementation. Section 6.3 of the paper reports the same thing: with
#' \eqn{E[\xi] = 1/3} and with \eqn{E[\xi] = -1/3} the average treatment effect
#' moves from 0.010 to 0.002 and to -0.004, about 0.01 in both cases and in the
#' same direction in both, which is the size of a shift that has been given
#' back rather than one of \eqn{\pm 1/3}. Read `ate` as the doubly robust
#' estimand under the assumed bias and `ate_gcomp` as the g-formula one, and
#' expect the sensitivity of the first to \eqn{\xi} to be the milder of the two.
#'
#' The recoupling is a fresh sequential Monte Carlo sweep, so `xi = 0` returns
#' the original analysis only up to the Monte Carlo error of that sweep, exactly
#' as two calls to [drbayes_pc()] on the same data differ. It is run with the
#' `control` of the original fit and against the tolerance that fit was judged
#' against, so that the same bar is applied to both couplings, and it costs
#' about what the coupling inside [drbayes_pc()] cost. What is not inherited is
#' the number of draws the mean being compared with that bar is taken over; the
#' next section is about what follows from that.
#'
#' \strong{Reading the result.}
#'
#' Importance reweighting can only report what the original draws already
#' covered. The reweighted posterior of the average treatment effect is the
#' original one shifted by about \eqn{-E[\xi]}, so the two overlap only while
#' \eqn{|E[\xi]|} is within a few posterior standard deviations of the effect
#' being estimated. That ratio, and not the size of \eqn{\xi} on its own, is
#' what decides whether the answer means anything: on a fit whose posterior
#' standard deviation is 0.11, the paper's own \eqn{Tri(0, 0.5)} has
#' \eqn{E[\xi] = 1/3}, three posterior standard deviations away, and leaves an
#' effective sample size of under 2 percent of the draws.
#'
#' Two things have to hold before the analysis is still doubly robust, and
#' `smc$converged` is both of them. The coupling of the reweighted draws has to
#' meet the moment condition (3.4), which on its own is
#' `smc$converged_coupling`. And the reweighting that came before it has to
#' have left enough of the original draws in play, which is `ess` against the
#' floor `ess_frac` of [drbayes_control()], the same floor the coupling holds
#' its own weights to. The second is not something the coupling can see: it
#' resamples the reweighted draws and its kernel moves every particle, so the
#' cloud looks fully populated to every diagnostic it computes even where 28 of
#' the 3000 draws are behind it. Nor is the moment condition evidence at that
#' point. The tolerance was set from the fit's own 3000 draws, while the mean
#' being compared with it is now an average over an effective few dozen and
#' carries several times the tolerance as its own Monte Carlo error, so it is
#' met about as often as not: six runs of the identical call at
#' \eqn{Tri(0, 0.3)} on the gaussian example met it four times.
#'
#' The credible interval is reported only where the effective sample size can
#' support it. A 2.5 percent point needs ten effective draws beyond it to be a
#' quantile rather than the most extreme handful that survived the reweighting,
#' so it needs 400 effective draws, and a median needs 20. Below that the
#' quantile is `NA` and the mean and the standard deviation are what there is
#' to read.
#'
#' When the weights degenerate there are three things to do, in order of how
#' much they help: run [drbayes_pc()] with more draws, since the effective
#' sample size grows in proportion to them; assume a sensitivity distribution
#' with less mass far from zero, which is a statement about the science and not
#' a convenience; or trace the analysis over a sequence of smaller \eqn{\xi} and
#' read off where the conclusion changes, which is what a sensitivity analysis
#' is for and needs no single distribution to be believed.
#'
#' \strong{The residual scale of a gaussian outcome.}
#'
#' The gaussian log likelihood ratio is a difference of squared residuals over
#' \eqn{2\sigma^2} summed over the treated observations, so \eqn{\sigma} is not
#' a detail of the weight but the thing that sets its size. A point estimate of
#' it conditions the answer on a number the posterior is uncertain about: on
#' the gaussian example of this package's tests, holding \eqn{\sigma} one
#' posterior standard deviation of its own either side of 1, which is
#' \eqn{\sigma / \sqrt{2n} = 0.035}, moves the shift the reweighting produces
#' by 13 percent and the effective sample size of the weights by a factor of
#' 2. By default \eqn{\sigma} is therefore not fixed: see the `sigma`
#' argument.
#'
#' @param object An object of class `"DRBayes"` from [drbayes_pc()], fitted with
#'   `control = drbayes_control(keep_particles = TRUE)`. The untilted draws it
#'   keeps are what the weights are applied to, and the design matrices and
#'   responses it keeps with them are what the weights are computed from, so a
#'   fit without them cannot be used.
#' @param xi The sensitivity distribution \eqn{g()}, as either a function of `n`
#'   returning `n` draws, such as the closure [xi_triangular()] builds, or a
#'   numeric vector of draws to resample. A single number is a point mass, and
#'   `xi = 0` leaves the weights uniform.
#' @param M Integer. Number of draws from `xi` used for the Monte Carlo integral
#'   in the weight. Default 100. The integrand concentrates far more sharply
#'   than \eqn{g()} does, so far fewer than `M` of the draws carry the integral;
#'   how many is reported as `xi_ess` and warned about when it is small. Raising
#'   `M` raises that count in proportion, and the running time of the weights
#'   with it, while the share of the draws that counts stays where it is: on the
#'   gaussian example of this package's tests, `M` of 100, 1000 and 5000 left a
#'   median of 5.8, 65 and 347 effective draws, 5.8 to 6.9 percent of `M`
#'   throughout. What can be bought is therefore the usual square root: ten
#'   times the work for a third of the Monte Carlo error. A \eqn{g()} with less
#'   mass far from zero costs nothing instead.
#'
#'   Algorithm 3 reads as `M` fresh draws for each of the `S` posterior draws.
#'   They are taken once and reused across the posterior draws here, which makes
#'   the analysis exact for the `M` atom distribution that was actually drawn,
#'   rather than an unbiased estimate of the analysis under \eqn{g()} itself,
#'   and lets the result report which distribution that was. What it costs is
#'   that those atoms are themselves random: at `M = 200` on the gaussian
#'   example, repeating the call moves the g-formula shift by 3 percent of
#'   itself where fresh draws move it by 0.3 percent. In `ate` it is not
#'   detectable, because the coupling gives most of the shift back and most of
#'   that scatter with it, leaving it well under the scatter of the sweep
#'   itself.
#' @param method `"common"` (default) for one \eqn{\xi} shared by all
#'   observations, Algorithm A.1, or `"per-observation"` for one \eqn{\xi_i} per
#'   observation, Algorithm 3. See Details.
#' @param sigma Residual standard deviation of a gaussian outcome model, which
#'   the likelihood ratio needs and the coefficient draws do not carry. Left
#'   NULL, one value is drawn for each coefficient draw from the conditional
#'   posterior of \eqn{\sigma^2} given those coefficients, an inverse gamma with
#'   shape \eqn{n/2} and rate half the residual sum of squares at them, so that
#'   the posterior uncertainty about the residual scale is carried into the
#'   weights rather than assumed away. That conditional is the outcome
#'   likelihood under a flat prior on \eqn{\log \sigma^2}, and it is used
#'   whatever sampler produced the draws: a fit records which prior its sampler
#'   put on \eqn{\sigma^2} but not the prior itself, so that prior is not
#'   reproduced, which at the default of [bayes_lm()] is worth 0.2 percent of a
#'   posterior standard deviation of \eqn{\sigma}. A single number fixes the
#'   scale instead, which understates the spread of the reweighted posterior by
#'   the amount quoted in Details; a vector of one value per posterior draw,
#'   from a sampler that returns its \eqn{\sigma^2} chain, is used as it stands
#'   and is the exact thing the default approximates. Ignored for a binary
#'   outcome.
#' @param ... Not used, and refused rather than ignored so that a mistyped
#'   argument is not swallowed silently.
#'
#' @return An object of class `"DRBayes_sensitivity"`, a list containing:
#' \describe{
#'   \item{ate}{Numeric vector of draws of the average treatment effect under
#'     the assumed unmeasured confounding, reweighted, resampled and coupled
#'     again. They carry equal weight.}
#'   \item{ate_gcomp}{The same draws before the coupling: the g-formula average
#'     treatment effect under the assumed unmeasured confounding, which is
#'     where the whole of the shift is. See Details.}
#'   \item{ate_original}{The `pc` draws of the fit that was reweighted, for
#'     comparison.}
#'   \item{weights}{Numeric vector of self-normalised sensitivity weights on the
#'     untilted draws, one per draw, summing to one. The resampling has already
#'     applied them; they are kept because they are what the diagnostics below
#'     are computed from.}
#'   \item{log_weights}{The weights of Algorithm 3 or A.1 on the log scale,
#'     before normalising. They run to large negative numbers, which is why they
#'     are never exponentiated on their own.}
#'   \item{ess}{Effective sample size of the weights, \eqn{1/\sum w^2}. The
#'     number of the original draws the answer really rests on.}
#'   \item{smc}{The diagnostics of the coupling of the reweighted draws, as
#'     [drbayes_pc()] reports them in its own `smc`: the tilting parameter
#'     `lambda`, the posterior mean `B_mean` of the moment condition (3.4) and
#'     the `tol` it is judged against, the number of steps taken, and the
#'     particle diagnostics. `converged` is whether the sensitivity analysis is
#'     still doubly robust, which takes the reweighting into account as well as
#'     the coupling: `converged_coupling` is the coupling's own verdict, and it
#'     is held to `ess_sensitivity` reaching `ess_floor`, which is `ess_frac`
#'     of the draws. See "Reading the result".}
#'   \item{xi}{The draws from \eqn{g()} that were used: a vector of length `M`
#'     for `method = "common"`, or a treated observations by `M` matrix for
#'     `method = "per-observation"`.}
#'   \item{xi_ess}{The effective number of those draws that carries the Monte
#'     Carlo integral in the weight, as the `median` and the `min` over the
#'     integrals computed. Usually far smaller than `M`, and equal to it only
#'     when every draw of \eqn{\xi} counts equally, as under a point mass.}
#'   \item{n_saturated}{How many (treated observation, posterior draw) pairs
#'     the shift took outside the range of the outcome, `saturated` `of` the
#'     pairs there are, and so zero for a gaussian outcome. See "The shift is
#'     on the conditional mean".}
#'   \item{method, M, n_draws}{The settings the weights were computed with.}
#'   \item{sigma}{The residual standard deviations the gaussian likelihood ratio
#'     used, one per draw unless a single number was supplied, or NA for a
#'     binary outcome, whose likelihood has no scale parameter.}
#'   \item{call}{The matched call.}
#' }
#'
#' @references
#' Orihara, S., Momozaki, T., & Sugasawa, S. (2025). Bayesian Doubly Robust
#' Causal Inference via Posterior Coupling. arXiv preprint arXiv:2506.04868.
#'
#' Daniels, M. J., Linero, A., & Roy, J. (2023). \emph{Bayesian Nonparametrics
#' for Causal Inference and Missing Data}. New York: Chapman and Hall/CRC.
#'
#' @examples
#' # generate_dataset() leaves the caller's random number stream alone, so
#' # seed the coupling itself to make the fit below reproducible.
#' set.seed(1)
#' data <- generate_dataset(nn = 120, pp = 0, seed = 1)
#' fit  <- drbayes_pc(YY ~ AA + WW1 + WW2 + WW3 + WW4,
#'                    AA ~ WW1 + WW2 + WW3 + WW4,
#'                    data = data,
#'                    outcome.model = bayes_lm, ps.model = bayes_logit,
#'                    mc = 900, bn = 300, chains = 2L, verbose = FALSE,
#'                    control = drbayes_control(keep_particles = TRUE,
#'                                              n_steps = 200))
#'
#' # A bias somewhere in (0, 0.1) and most likely at its largest, which is the
#' # shape Section 6.3 of the paper assumes. What a sensitivity parameter has
#' # to be read against is not the size of the effect, 110 here, but the
#' # posterior standard deviation of it, which is the sd column of the summary
#' # below and about 0.16. The mean of xi, 0.067, is a little under half of
#' # that: a large bias to entertain, and one the original draws can still
#' # represent, since two thirds of them stay effective.
#' sens <- drbayes_sensitivity(fit, xi = xi_triangular(0, 0.1), M = 200)
#' sens
#' summary(sens)
#'
#' # Whether the answer is still doubly robust: the coupling met the moment
#' # condition, on a reweighting that had draws left to spare.
#' sens$smc$converged
#'
#' @seealso [xi_triangular()] for the triangular sensitivity distributions of
#'   Section 6.3, [drbayes_pc()] for the fit being reweighted, and
#'   [drbayes_control()] for `keep_particles`.
#'
#' @export
drbayes_sensitivity <- function(object, xi, M = 100,
                                method = c("common", "per-observation"),
                                sigma = NULL, ...) {

  call_info <- match.call()
  method    <- match.arg(method)
  refuse_unused_arguments(list(...))

  if (!inherits(object, "DRBayes")) {
    stop("object must be a fit of class \"DRBayes\" returned by drbayes_pc(), ",
         "but is of class ", class(object)[1L], ".")
  }
  validate_drbayes_object(object)
  M <- validate_draw_count(M)

  parts   <- sensitivity_particles(object)
  d       <- parts$model_data
  treated <- which(d$A == 1L)

  if (length(treated) == 0L) {
    stop("No observation is treated, so the shifted outcome model ",
         "m_A(X; beta) + A xi is the fitted model and there is nothing for a ",
         "sensitivity analysis to vary.")
  }

  sigma <- resolve_sensitivity_sigma(sigma, object$link, d, parts$betas.otc)

  # Only treated observations are shifted, by A xi, so every control
  # observation contributes a likelihood ratio of exactly one and can be left
  # out of the product entirely. Their linear predictors are what the shifted
  # and the unshifted fitted means are both computed from.
  eta      <- tcrossprod(d$Z.lm[treated, , drop = FALSE], parts$betas.otc)
  xi_draws <- draw_xi(xi, method, M, nrow(eta), rownames(eta))

  integral <- sensitivity_log_weights(eta, d$Y[treated], xi_draws, object$link,
                                      sigma, method)
  log_w    <- integral$log_weights

  if (!all(is.finite(log_w))) {
    stop("The shifted outcome model gives ", sum(!is.finite(log_w)), " of ",
         length(log_w), " draws a likelihood ratio of zero, so the reweighted ",
         "posterior is undefined. This happens when xi is far larger than the ",
         "outcome model can absorb: compare it with sigma for an identity ",
         "link, and with the fitted probabilities of the treated for a binary ",
         "one.")
  }

  weights <- normalised_weights(log_w)
  ess     <- 1 / sum(weights^2)
  control <- object$control %||% drbayes_control()

  for (note in c(sensitivity_ess_message(ess, length(weights),
                                         control$ess_frac),
                 sensitivity_xi_message(integral$xi_ess, M),
                 sensitivity_saturation_message(integral$n_saturated,
                                                integral$n_cells))) {
    warning(note, call. = FALSE)
  }

  # Steps (2) and (3) of Section 6.2: the reweighted draws are the posterior
  # under unmeasured confounding, and the doubly robust posterior of the
  # average treatment effect comes from coupling them. The resampling is
  # stratified by the moment condition the coupling is about to drive to zero,
  # which is where its rounding error matters; an equally weighted set comes
  # back untouched, so a sensitivity parameter of zero costs nothing.
  idx      <- systematic_resample(
    weights, moment_ipw(parts$betas.ps, parts$betas.otc, d))
  coupling <- run_tilting(parts$betas.ps[idx, , drop = FALSE],
                          parts$betas.otc[idx, , drop = FALSE], d,
                          control = control,
                          tol = sensitivity_coupling_tol(object),
                          method = "smc")

  result <- list(
    ate          = sensitivity_ate(coupling$betas.otc, parts),
    ate_gcomp    = sensitivity_ate(parts$betas.otc[idx, , drop = FALSE], parts),
    ate_original = object$pc,
    weights      = weights,
    log_weights  = log_w,
    ess          = ess,
    smc          = sensitivity_verdict(coupling$smc, ess, length(weights),
                                       control$ess_frac),
    xi           = xi_draws,
    xi_ess       = integral$xi_ess,
    n_saturated  = c(saturated = integral$n_saturated,
                     of = integral$n_cells),
    method       = method,
    M            = M,
    n_draws      = length(weights),
    sigma        = sigma,
    call         = call_info
  )

  class(result) <- c("DRBayes_sensitivity", "list")
  result
}


#' Triangular sensitivity distributions
#'
#' @description
#' Builds the sensitivity distribution \eqn{g()} that Section 6.3 of Orihara,
#' Momozaki and Sugasawa (2025) uses: a triangular density on `(lower, upper)`
#' peaking at `mode`. The result is a function of `n` returning `n` draws, ready
#' to hand to the `xi` argument of [drbayes_sensitivity()].
#'
#' @details
#' The paper's two scenarios are the triangular distribution on \eqn{(0, 0.5)}
#' peaking at \eqn{0.5} and the one on \eqn{(-0.5, 0)} peaking at \eqn{-0.5},
#' representing unmeasured confounding of a positive and of a negative sign. In
#' both the peak is at the endpoint further from zero, which is the pessimistic
#' reading of a bounded bias, and that endpoint is what `mode` defaults to.
#'
#' Those endpoints are on the scale of the outcome itself, because the shift is
#' on the conditional mean: for the binary outcome of Section 6.3 they are risk
#' differences, and 0.5 is not a moderate bias there but a very large one. That
#' section reports an average treatment effect of 0.010 with 95 percent
#' interval (-0.111, 0.130), a posterior standard deviation of about 0.06, so
#' \eqn{E[\xi] = 1/3} is five of them. What a sensitivity parameter has to be
#' read against is that standard deviation, and a bias several of them wide is
#' one the original draws cannot represent, whatever the units of the outcome;
#' see the "Reading the result" section of [drbayes_sensitivity()].
#'
#' @param lower,upper The support of the sensitivity parameter. Any bias outside
#'   `(lower, upper)` is being ruled out, so these are the assumption the whole
#'   analysis rests on.
#' @param mode Where the density peaks. Defaults to whichever of `lower` and
#'   `upper` is further from zero, the largest bias the support allows.
#'
#' @return A function of `n` returning `n` draws, suitable as the `xi` argument
#'   of [drbayes_sensitivity()].
#'
#' @references
#' Orihara, S., Momozaki, T., & Sugasawa, S. (2025). Bayesian Doubly Robust
#' Causal Inference via Posterior Coupling. arXiv preprint arXiv:2506.04868.
#'
#' @examples
#' positive <- xi_triangular(0, 0.5)
#' negative <- xi_triangular(-0.5, 0)
#' summary(positive(1000))
#'
#' # A bias of unknown sign, symmetric around zero and most likely absent.
#' either <- xi_triangular(-0.5, 0.5, mode = 0)
#'
#' @seealso [drbayes_sensitivity()]
#'
#' @export
xi_triangular <- function(lower, upper, mode = NULL) {

  endpoints <- list(lower = lower, upper = upper)
  for (name in names(endpoints)) {
    value <- endpoints[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
      stop(name, " must be a single finite number, an endpoint of the range ",
           "the sensitivity parameter is allowed to take.")
    }
  }
  if (lower >= upper) {
    stop("lower must be smaller than upper, but they are ", lower, " and ",
         upper, ".")
  }

  if (is.null(mode)) {
    mode <- if (abs(lower) > abs(upper)) lower else upper
  }
  if (!is.numeric(mode) || length(mode) != 1L || !is.finite(mode) ||
      mode < lower || mode > upper) {
    stop("mode must be a single number inside [", lower, ", ", upper,
         "], the value of the sensitivity parameter the density peaks at.")
  }

  function(n) extraDistr::rtriang(n, a = lower, b = upper, c = mode)
}


#' Print a sensitivity analysis
#'
#' @description
#' A one screen summary of a [drbayes_sensitivity()] analysis: the call, the
#' sensitivity distribution, the effective sample size of the weights, whether
#' the analysis is still doubly robust, and the average treatment effect before
#' and after.
#'
#' @details
#' The diagnostics come before the estimates on purpose. A reweighted posterior
#' is only as good as the number of draws that carry weight, and it is doubly
#' robust only if the coupling that follows the reweighting converged on a
#' reweighting that left it draws to work with. Those are two lines of the
#' output: the moment condition, and the verdict on the analysis as a whole.
#'
#' The credible interval is left out where the effective sample size cannot
#' support it. A 2.5 percent point with fewer than ten effective draws beyond
#' it is the most extreme handful of them, not a quantile, which takes 400
#' effective draws.
#'
#' @param x An object of class `"DRBayes_sensitivity"`.
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
#' data <- generate_dataset(nn = 120, pp = 0, seed = 1)
#' fit  <- drbayes_pc(YY ~ AA + WW1 + WW2 + WW3 + WW4,
#'                    AA ~ WW1 + WW2 + WW3 + WW4,
#'                    data = data,
#'                    outcome.model = bayes_lm, ps.model = bayes_logit,
#'                    mc = 900, bn = 300, chains = 2L, verbose = FALSE,
#'                    control = drbayes_control(keep_particles = TRUE,
#'                                              n_steps = 200))
#' drbayes_sensitivity(fit, xi = xi_triangular(0, 0.1))
#'
#' @seealso [summary.DRBayes_sensitivity()]
#' @export
print.DRBayes_sensitivity <- function(x, digits = 3, ...) {

  digits <- validate_digits(digits)

  cat("Sensitivity analysis for unmeasured confounding\n")

  if (!is.null(x$call)) {
    cat("\nCall:\n")
    cat(paste0(deparse(x$call), collapse = "\n"), "\n", sep = "")
  }

  cat("\n")
  cat("Sensitivity parameter: ", describe_sensitivity_method(x$method), "\n",
      sep = "")
  cat("                       ", describe_sensitivity_xi(x, digits), "\n",
      sep = "")
  cat("                       ", describe_sensitivity_atoms(x, digits), "\n",
      sep = "")
  print_sensitivity_sigma(x$sigma, digits)
  print_sensitivity_saturation(x$n_saturated)
  cat("Posterior draws:       ", x$n_draws, "\n", sep = "")
  print_sensitivity_ess(x$ess, x$n_draws, x$smc$ess_frac %||% 0.1)
  print_sensitivity_coupling(x$smc, digits)

  cat("\nAverage treatment effect\n")
  print(format_estimand_table(sensitivity_estimand_table(x), digits,
                              c("mean", "2.5%", "97.5%")))
  cat("\n")
  cat_wrapped(sensitivity_table_note())
  print_sensitivity_interval_note(x$ess)

  invisible(x)
}


#' Summarise a sensitivity analysis
#'
#' @description
#' Collects the posterior summaries of the reweighted and the original average
#' treatment effect into a small data frame, together with the sensitivity
#' distribution the reweighting assumed, the effective sample size of the
#' weights it produced, and the coupling that followed it.
#'
#' @param object An object of class `"DRBayes_sensitivity"`.
#' @param ... Ignored, present for consistency with the generic.
#'
#' @return An object of class `"summary.DRBayes_sensitivity"`, a list with
#'   components
#'   \describe{
#'     \item{estimands}{Data frame with rows `sensitivity`, `reweighted` and
#'       `original` and columns `mean`, `sd`, `2.5%`, `50%` and `97.5%`. A
#'       quantile the effective sample size cannot support is `NA`.}
#'     \item{xi}{Data frame of the mean, standard deviation and range of the
#'       draws from the sensitivity distribution.}
#'     \item{call, method, M, ess, xi_ess, n_draws, sigma, smc}{Taken unchanged
#'       from `object`.}
#'   }
#'
#' @examples
#' # generate_dataset() leaves the caller's random number stream alone, so
#' # seed the coupling itself to make the fit below reproducible.
#' set.seed(1)
#' data <- generate_dataset(nn = 120, pp = 0, seed = 1)
#' fit  <- drbayes_pc(YY ~ AA + WW1 + WW2 + WW3 + WW4,
#'                    AA ~ WW1 + WW2 + WW3 + WW4,
#'                    data = data,
#'                    outcome.model = bayes_lm, ps.model = bayes_logit,
#'                    mc = 900, bn = 300, chains = 2L, verbose = FALSE,
#'                    control = drbayes_control(keep_particles = TRUE,
#'                                              n_steps = 200))
#' summary(drbayes_sensitivity(fit, xi = xi_triangular(0, 0.1)))
#'
#' @seealso [print.DRBayes_sensitivity()]
#' @export
summary.DRBayes_sensitivity <- function(object, ...) {

  values <- as.numeric(object$xi)

  result <- list(
    estimands = sensitivity_estimand_table(object),
    xi        = data.frame(mean = mean(values), sd = stats::sd(values),
                           min = min(values), max = max(values)),
    call        = object$call,
    method      = object$method,
    M           = object$M,
    ess         = object$ess,
    xi_ess      = object$xi_ess,
    n_saturated = object$n_saturated,
    n_draws     = object$n_draws,
    sigma       = object$sigma,
    smc         = object$smc
  )

  class(result) <- c("summary.DRBayes_sensitivity", "list")
  result
}


#' @rdname summary.DRBayes_sensitivity
#'
#' @param x An object of class `"summary.DRBayes_sensitivity"`.
#' @param digits Number of significant digits. The table of estimands is shown
#'   to the number of decimal places that gives the posterior standard deviation
#'   this many significant digits, so that a wide credible interval around a
#'   large effect does not round away. Defaults to 3.
#'
#' @return `print.summary.DRBayes_sensitivity` returns `x` invisibly.
#'
#' @export
print.summary.DRBayes_sensitivity <- function(x, digits = 3, ...) {

  digits <- validate_digits(digits)

  cat("Sensitivity analysis for unmeasured confounding\n")

  if (!is.null(x$call)) {
    cat("\nCall:\n")
    cat(paste0(deparse(x$call), collapse = "\n"), "\n", sep = "")
  }

  cat("\n")
  cat("Sensitivity parameter: ", describe_sensitivity_method(x$method), "\n",
      sep = "")
  cat("Draws from g():        ", x$M, "\n", sep = "")
  cat("                       ", describe_sensitivity_atoms(x, digits), "\n",
      sep = "")
  print_sensitivity_sigma(x$sigma, digits)
  print_sensitivity_saturation(x$n_saturated)
  cat("Posterior draws:       ", x$n_draws, "\n", sep = "")
  print_sensitivity_ess(x$ess, x$n_draws, x$smc$ess_frac %||% 0.1)
  print_sensitivity_coupling(x$smc, digits)

  cat("\nSensitivity parameter xi\n")
  print(x$xi, digits = digits, row.names = FALSE)

  cat("\nPosterior summary of the average treatment effect\n")
  print(format_estimand_table(x$estimands, digits, names(x$estimands)))
  cat("\n")
  cat_wrapped(sensitivity_table_note())
  print_sensitivity_interval_note(x$ess)

  invisible(x)
}


#' What the three rows of the estimand table are
#'
#' @keywords internal
#' @noRd
sensitivity_table_note <- function() {
  paste("The sensitivity row is the posterior under the assumed unmeasured",
        "confounding; the reweighted row is the same draws before the moment",
        "condition was imposed on them again, and holds the whole of the",
        "shift; the original row is the fit that was reweighted. See the",
        "\"Reweight first, then couple\" section of ?drbayes_sensitivity for",
        "why the first two differ.")
}


#' Refuse arguments that a function does not use
#'
#' A dots argument that quietly absorbs everything turns a mistyped argument
#' name into a silently ignored setting, which here would mean reporting a
#' sensitivity analysis the caller did not ask for.
#'
#' @keywords internal
#' @noRd
refuse_unused_arguments <- function(dots) {
  if (length(dots) == 0L) {
    return(invisible(NULL))
  }
  named <- names(dots)
  if (is.null(named)) named <- rep("", length(dots))
  named[!nzchar(named)] <- "<unnamed>"

  if ("data" %in% named) {
    stop("drbayes_sensitivity() takes no data argument. A fit kept with ",
         "keep_particles = TRUE carries the design matrices and the responses ",
         "it was made on, and the weights are computed from those, so that a ",
         "data frame that has since changed cannot be reweighted against by ",
         "mistake. Drop the argument.", call. = FALSE)
  }
  stop("drbayes_sensitivity() does not use the argument(s) ", toString(named),
       ". The sensitivity distribution goes in xi, the number of draws from ",
       "it in M; see ?drbayes_sensitivity.", call. = FALSE)
}


#' Number of Monte Carlo draws from the sensitivity distribution
#'
#' @keywords internal
#' @noRd
validate_draw_count <- function(M) {
  if (!is.numeric(M) || length(M) != 1L || !is.finite(M) || M != round(M) ||
      M < 1L) {
    stop("M must be a single positive integer, the number of draws from the ",
         "sensitivity distribution used for the integral in the weight.",
         call. = FALSE)
  }
  as.integer(M)
}


#' Everything a sensitivity analysis needs from a fit
#'
#' The weights are ratios of outcome model likelihoods at the draws of the
#' original posterior of equation (3.3), so what they act on is the untilted
#' draws and what they are computed from is the model frame the fit was built
#' on. A fit kept with `keep_particles = TRUE` carries both, which is what lets
#' this function refuse to guess at either.
#'
#' @return A list with the untilted draws as `betas.otc` and `betas.ps`, the
#'   [tilting_data()] object as `model_data`, and the design matrices `Z.lm1`
#'   and `Z.lm0` the average treatment effect is averaged over. The untilted
#'   draws are deliberately not called `outcome` and `ps`, which in a fit are
#'   the coupled ones.
#'
#' @keywords internal
#' @noRd
sensitivity_particles <- function(object) {

  parts  <- object$particles
  needed <- c("outcome_untilted", "ps_untilted", "model_data", "Z.lm1", "Z.lm0")

  if (!is.list(parts) || !all(needed %in% names(parts))) {
    stop("A sensitivity analysis reweights the draws of the original ",
         "posterior and needs the data the fit was made on, and this fit ",
         "carries neither. Refit keeping them: drbayes_pc(..., control = ",
         "drbayes_control(keep_particles = TRUE)).", call. = FALSE)
  }

  out <- list(betas.otc  = as.matrix(parts$outcome_untilted),
              betas.ps   = as.matrix(parts$ps_untilted),
              model_data = parts$model_data,
              Z.lm1      = parts$Z.lm1,
              Z.lm0      = parts$Z.lm0)

  if (nrow(out$betas.otc) != length(object$pc)) {
    stop("This fit holds ", length(object$pc), " draws of the average ",
         "treatment effect but ", nrow(out$betas.otc), " outcome model ",
         "draws. They come from the same sweep and must match, so this ",
         "object cannot be reweighted.", call. = FALSE)
  }
  if (!all(is.finite(out$betas.otc)) || !all(is.finite(out$betas.ps))) {
    stop("Some of the particles of this fit are missing or infinite, so no ",
         "likelihood ratio can be computed from them. Refit drbayes_pc() and ",
         "check the convergence diagnostics it reports.", call. = FALSE)
  }

  out
}


#' Residual standard deviations for a gaussian likelihood ratio
#'
#' The gaussian log likelihood ratio divides by the error variance, and the
#' samplers return coefficient draws without the error variance that went with
#' them. Rather than plug in one value for all of them, one is drawn for each
#' coefficient draw from an inverse gamma with shape \eqn{n / 2} and rate
#' \eqn{RSS(\beta) / 2}, the conditional posterior of \eqn{\sigma^2} given
#' those coefficients under the reference prior \eqn{p(\sigma^2) \propto
#' 1/\sigma^2}. Composed with a draw of the coefficients that is a draw from
#' the joint posterior of \eqn{(\beta, \sigma^2)}, so the weights carry the
#' posterior uncertainty about the residual scale instead of conditioning on an
#' estimate of it.
#'
#' It is the likelihood alone. A fit records which prior its outcome sampler
#' used but not the prior itself, so the prior that sampler put on
#' \eqn{\sigma^2} cannot be reproduced here and is not: for the inverse gamma
#' \eqn{(a, b)} that [bayes_lm()] documents, what is dropped is \eqn{a} added
#' to the shape and \eqn{b} to the rate. The two pull against each other at its
#' default \eqn{(1, 1)} and very nearly cancel, moving the mean of \eqn{\sigma}
#' by 0.2 percent of one of its own posterior standard deviations at the 400
#' observations of this package's tests. A sampler that returns its
#' \eqn{\sigma^2} chain should pass it through the `sigma` argument, which is
#' exact and approximates nothing.
#'
#' @param betas Draws by coefficients matrix of the untilted outcome model.
#'
#' @return One standard deviation per draw, or the single number supplied, or
#'   `NA` for a link with no scale parameter.
#'
#' @keywords internal
#' @noRd
resolve_sensitivity_sigma <- function(sigma, link, d, betas) {

  n <- length(d$Y)
  S <- nrow(betas)

  if (link != "identity") {
    if (!is.null(sigma)) {
      warning("sigma is the residual standard deviation of a gaussian outcome ",
              "model and this fit has a ", link, " link, so it is ignored.",
              call. = FALSE)
    }
    return(NA_real_)
  }

  if (!is.null(sigma)) {
    if (!is.numeric(sigma) || !all(is.finite(sigma)) || any(sigma <= 0) ||
        !(length(sigma) %in% c(1L, S))) {
      stop("sigma must be positive and finite, and either a single number or ",
           "one value for each of the ", S, " posterior draws, the residual ",
           "standard deviation of the outcome model.", call. = FALSE)
    }
    return(as.numeric(sigma))
  }

  # Expanded as Y'Y - 2 b'Z'Y + b'Z'Zb so that no draws by observations matrix
  # of fitted values is ever formed.
  ZtZ <- crossprod(d$Z.lm)
  ZtY <- drop(crossprod(d$Z.lm, d$Y))
  rss <- sum(d$Y^2) - 2 * drop(betas %*% ZtY) +
    rowSums((betas %*% ZtZ) * betas)
  rss <- pmax(rss, .Machine$double.eps)

  sigma <- sqrt(1 / stats::rgamma(S, shape = n / 2, rate = rss / 2))
  if (!all(is.finite(sigma)) || any(sigma <= 0)) {
    stop("The residual standard deviation of the outcome model came out as ",
         "zero or infinite for ", sum(!is.finite(sigma) | sigma <= 0), " of ",
         S, " draws, so the gaussian likelihood ratio is undefined. Supply it ",
         "directly with sigma = <a positive number>.", call. = FALSE)
  }
  sigma
}


#' Draws from the sensitivity distribution
#'
#' Algorithm 3 needs one draw per observation per Monte Carlo sample and
#' Algorithm A.1 one draw per Monte Carlo sample, so the two differ only in how
#' many values are drawn and how they are laid out.
#'
#' @param n Number of treated observations.
#' @param observations Their names, used to label the rows of the
#'   per-observation matrix.
#'
#' @return A vector of length `M`, or an `n` by `M` matrix.
#'
#' @keywords internal
#' @noRd
draw_xi <- function(xi, method, M, n, observations) {

  size <- if (method == "common") M else n * M

  if (is.function(xi)) {
    values <- xi(size)
    if (!is.numeric(values) || length(values) != size) {
      stop("xi must be a function of n returning n draws, but xi(", size,
           ") returned ", length(values), " value(s) of class ",
           class(values)[1L], ".", call. = FALSE)
    }
  } else if (is.numeric(xi) && length(xi) > 0L) {
    # A numeric xi is a set of draws from g() to resample, so that draws
    # produced elsewhere can be reused. Indexing rather than sample() itself,
    # because sample(x, ...) on a single number samples from seq_len(x).
    values <- xi[sample.int(length(xi), size, replace = TRUE)]
  } else {
    stop("xi must be either a function of n returning n draws from the ",
         "sensitivity distribution, as xi_triangular() builds, or a non-empty ",
         "numeric vector of draws to resample from.", call. = FALSE)
  }

  if (!all(is.finite(values))) {
    stop("The draws of the sensitivity parameter must all be finite, but ",
         sum(!is.finite(values)), " of ", size, " are missing or infinite.",
         call. = FALSE)
  }

  values <- as.numeric(values)
  if (method == "common") {
    return(values)
  }
  matrix(values, nrow = n, ncol = M, dimnames = list(observations, NULL))
}


#' Log importance weights of Algorithm 3 and of Algorithm A.1
#'
#' Everything here is done on the log scale. The weight is a product over
#' observations of ratios below one, so at any realistic sample size the product
#' itself is not a representable number: at n = 1000 a thousand ratios of 0.9
#' multiply to 1e-46, and it falls away from there.
#'
#' @param eta Treated observations by draws matrix of linear predictors.
#' @param y Outcomes of the treated observations.
#' @param xi The draws from [draw_xi()].
#' @param sigma Residual standard deviation, one per draw or one in all, for the
#'   identity link.
#'
#' @return A list with `log_weights`, one per posterior draw; `xi_ess`, the
#'   median and the smallest effective number of draws from \eqn{g()} behind the
#'   Monte Carlo integrals; and `n_saturated` of `n_cells` fitted means the
#'   shift took outside the range of the outcome.
#'
#' @keywords internal
#' @noRd
sensitivity_log_weights <- function(eta, y, xi, link, sigma, method) {

  # The unshifted outcome model does not depend on xi, so the pieces of it the
  # ratio needs are computed once here rather than M times in the loops below.
  scale <- outcome_scale(eta, link, sigma)

  if (method == "common") {
    # Algorithm A.1: one xi is shared, so the sum over observations happens
    # inside the integral and the average over the M draws is taken once.
    M <- length(xi)
    log_terms <- matrix(NA_real_, nrow = ncol(eta), ncol = M)
    for (m in seq_len(M)) {
      log_terms[, m] <- colSums(outcome_log_ratio(scale, y, xi[m]))
    }
    averaged <- log_mean_exp_rows(log_terms)
    return(list(log_weights = averaged$log_mean,
                xi_ess      = summarise_atoms(averaged$n_effective),
                n_saturated = count_saturated(scale, xi),
                n_cells     = length(eta)))
  }

  # Algorithm 3: every observation has its own xi_i, so every observation gets
  # its own average over the M draws and the sum over observations comes last.
  # The M matrices of log ratios are far too large to hold at once, so their
  # log-sum-exp is accumulated as they are produced: `reference` is the running
  # maximum that keeps each exponential at or below one, `total` the sum
  # measured against it, rescaled whenever the maximum moves, and `squared` the
  # same sum of squares, which is what the effective number of draws needs.
  # The running maximum starts at the first log ratio rather than at zero,
  # because log ratios of a few hundred below zero are ordinary here and every
  # exponential taken against a maximum of zero would underflow.
  M         <- ncol(xi)
  reference <- outcome_log_ratio(scale, y, xi[, 1L])
  total     <- matrix(1, nrow = nrow(eta), ncol = ncol(eta))
  squared   <- total

  for (m in seq_len(M)[-1L]) {
    log_ratio <- outcome_log_ratio(scale, y, xi[, m])
    updated   <- pmax(reference, log_ratio)
    rescale   <- exp(reference - updated)
    term      <- exp(log_ratio - updated)
    total     <- total * rescale + term
    squared   <- squared * rescale^2 + term^2
    reference <- updated
  }

  list(log_weights = colSums(reference + log(total / M)),
       xi_ess      = summarise_atoms(total^2 / squared),
       n_saturated = count_saturated(scale, xi),
       n_cells     = length(eta))
}


#' The unshifted outcome model, in the pieces a shift is applied to
#'
#' The shifted model of Section 6.2 shifts the conditional mean, so what the
#' likelihood ratio needs is the fitted mean and, for a binary outcome, its
#' complement and the logarithms of both. None of them depends on the
#' sensitivity parameter, so they are computed once for all `M` draws of it.
#'
#' The complement of a fitted probability is computed from the linear predictor
#' rather than as `1 - mu`, which for a fitted value near one would be all
#' rounding error.
#'
#' @param eta Treated observations by draws matrix of linear predictors.
#' @param sigma One residual standard deviation, or one per draw.
#'
#' @return A list carrying `link`, either `"identity"` or `"binary"`.
#'
#' @keywords internal
#' @noRd
outcome_scale <- function(eta, link, sigma) {
  binary <- function(mu, mu_c) {
    list(link = "binary", mu = mu, mu_c = mu_c,
         log_mu = log_bounded(mu), log_mu_c = log_bounded(mu_c))
  }
  switch(link,
         "identity" = list(link = "identity", mu = eta, sigma = sigma),
         "logit"    = binary(stats::plogis(eta), stats::plogis(-eta)),
         "probit"   = binary(stats::pnorm(eta), stats::pnorm(-eta)),
         stop("The outcome link \"", link, "\" has no likelihood the ",
              "sensitivity analysis knows how to shift.", call. = FALSE))
}


#' Log likelihood ratio of the shifted outcome model
#'
#' The log of \eqn{\exp\{-f(Y | A, X; \beta, \xi)\} / \exp\{-f(Y | A, X;
#' \beta)\}} for the shifted outcome model \eqn{m_A(X; \beta, \xi) = m_A(X;
#' \beta) + A\xi} of Section 6.2, for treated observations, where the shift
#' applies. \eqn{m_A(X; \beta)} is the conditional mean, so the shift is added
#' to the mean: to the linear predictor for an identity link, and to the fitted
#' probability for a binary one.
#'
#' A shifted probability outside \eqn{[0, 1]} is not a Bernoulli mean, so it is
#' held at the boundary; see [log_bounded()].
#'
#' @param scale The pieces [outcome_scale()] prepared.
#' @param y Outcomes of the treated observations.
#' @param shift One value of the sensitivity parameter, or one per observation.
#'
#' @return A matrix the shape of the linear predictors.
#'
#' @keywords internal
#' @noRd
outcome_log_ratio <- function(scale, y, shift) {
  if (scale$link == "identity") {
    # A shift is one value per observation and recycles down the rows of the
    # matrix, whereas sigma is one value per draw and belongs to a column, so
    # the two cannot be combined into one expression.
    residual <- shift * (2 * (y - scale$mu) - shift)
    return(if (length(scale$sigma) == 1L) {
      residual / (2 * scale$sigma^2)
    } else {
      sweep(residual, 2L, 2 * scale$sigma^2, "/")
    })
  }

  # The log density of a Bernoulli observation with mean mu is log(mu) when y
  # is 1 and log(1 - mu) when y is 0, and the shift moves mu.
  y * (log_bounded(scale$mu + shift) - scale$log_mu) +
    (1 - y) * (log_bounded(scale$mu_c - shift) - scale$log_mu_c)
}


#' Logarithm of a probability, held inside the representable unit interval
#'
#' Two things put a probability on the boundary here: a fitted value so extreme
#' that double precision cannot tell it from zero or one, and a shifted mean
#' \eqn{m_A(X; \beta) + \xi} that leaves \eqn{[0, 1]}, where the shifted model
#' is not a distribution at all. Both are held one machine epsilon inside the
#' boundary, which bounds what a single observation can contribute to a log
#' weight at \eqn{\log(1/\epsilon) = 36} and keeps the weights finite. Applying
#' the same bound to the shifted and the unshifted mean is what makes a
#' sensitivity parameter of zero give a log ratio of exactly zero however
#' extreme the fitted values are.
#'
#' @keywords internal
#' @noRd
log_bounded <- function(p) {
  smallest <- .Machine$double.eps
  log(pmin(pmax(p, smallest), 1 - smallest))
}


#' Fitted values a shift takes outside the range of the outcome model
#'
#' The count of (treated observation, posterior draw) pairs for which some draw
#' of \eqn{\xi} puts the shifted mean outside \eqn{[0, 1]}, where
#' [log_bounded()] holds it at the boundary and the bias the pair carries is
#' smaller than \eqn{\xi}. Always zero for an identity link, whose mean is
#' unbounded.
#'
#' @param xi The draws from [draw_xi()]: a vector shared by every observation,
#'   or one row of draws per observation.
#'
#' @keywords internal
#' @noRd
count_saturated <- function(scale, xi) {
  if (scale$link == "identity") {
    return(0L)
  }
  largest  <- if (is.matrix(xi)) apply(xi, 1L, max) else max(xi)
  smallest <- if (is.matrix(xi)) apply(xi, 1L, min) else min(xi)
  sum(scale$mu + largest > 1 | scale$mu + smallest < 0)
}


#' Log of the mean of exponentials, along the rows of a matrix
#'
#' Taking out the largest value of each row before exponentiating is what makes
#' the average representable: the entries are log likelihood ratios summed over
#' observations, so they run to several hundred in absolute value.
#'
#' @return A list with `log_mean` and, per row, the `n_effective` number of
#'   terms the average rests on.
#'
#' @keywords internal
#' @noRd
log_mean_exp_rows <- function(x) {
  largest <- apply(x, 1L, max)
  scaled  <- exp(x - largest)
  total   <- rowSums(scaled)
  list(log_mean    = largest + log(total / ncol(x)),
       n_effective = total^2 / rowSums(scaled^2))
}


#' How many draws from g() an integral really rests on
#'
#' @param n_effective One effective count per Monte Carlo integral computed.
#'
#' @keywords internal
#' @noRd
summarise_atoms <- function(n_effective) {
  c(median = stats::median(n_effective), min = min(n_effective))
}


#' The tolerance the coupling of the reweighted draws is judged against
#'
#' The tolerance the original fit was held to, so that the same bar is applied
#' to both couplings. Recomputing it from the resampled particles would be
#' optimistic: duplicated particles make the spread of the moment condition
#' look better resolved than the effective sample size behind it really is.
#'
#' What is inherited is the bar and not the noise of the thing being measured
#' against it. [moment_tolerance()] set this to the standard error of the
#' posterior mean of \eqn{B_n} over the fit's own draws, and the mean it now
#' judges is taken over a cloud whose effective size is that of the sensitivity
#' weights, which can be a fiftieth of the draws and carry several times the
#' tolerance as its own Monte Carlo error. `|B_mean| < tol` is then a coin
#' flip rather than evidence, which is why the verdict of
#' [sensitivity_verdict()] holds the effective sample size to a floor as well.
#'
#' @return A positive number, or NULL to let [run_tilting()] set it.
#'
#' @keywords internal
#' @noRd
sensitivity_coupling_tol <- function(object) {
  tol <- object$smc$tol
  if (is.null(tol) || !is.finite(tol) || tol <= 0) {
    return(NULL)
  }
  tol
}


#' The g-computation average treatment effect of a set of draws
#'
#' Equation (3.7) at the design matrices the fit kept, so that the sensitivity
#' analysis and [drbayes_pc()] average the same fitted contrast over the same
#' empirical distribution of the covariates.
#'
#' @keywords internal
#' @noRd
sensitivity_ate <- function(betas, parts) {
  inverse_link <- parts$model_data$inverse_link
  rowMeans(inverse_link(tcrossprod(betas, parts$Z.lm1))) -
    rowMeans(inverse_link(tcrossprod(betas, parts$Z.lm0)))
}


#' Posterior summaries of the reweighted and the original effect
#'
#' The sensitivity row comes first because it is the answer to the question the
#' analysis asked. The middle row is the same draws before the coupling, which
#' is where the whole of the shift is; the distance between the two rows is
#' what imposing the moment condition again gave back. All three are ordinary
#' unweighted summaries: the sensitivity draws were resampled to equal weight
#' before they were coupled, so their weights have already been spent.
#'
#' @return A data frame with rows `sensitivity`, `reweighted` and `original`.
#'
#' @keywords internal
#' @noRd
sensitivity_estimand_table <- function(x) {
  probs <- c(0.025, 0.5, 0.975)

  summarise_one <- function(draws, n_effective) {
    quantiles <- stats::quantile(draws, probs, names = FALSE)
    quantiles[!quantile_is_resolved(n_effective, probs)] <- NA_real_
    c(mean = mean(draws), sd = stats::sd(draws), quantiles)
  }

  out <- as.data.frame(rbind(
    sensitivity = summarise_one(x$ate, x$ess),
    reweighted  = summarise_one(x$ate_gcomp, x$ess),
    original    = summarise_one(x$ate_original, length(x$ate_original))))
  names(out) <- c("mean", "sd", "2.5%", "50%", "97.5%")
  out
}


#' Whether a quantile rests on enough effective draws to be a quantile
#'
#' A probability `p` puts an expected `n_effective * min(p, 1 - p)` draws in the
#' shorter of the two tails it separates, and the reported quantile is that
#' tail's own boundary, so it carries about `1 / sqrt()` of that count as
#' relative Monte Carlo error. Ten of them is the same third of the value that
#' [sensitivity_xi_message()] treats as the point worth saying out loud, and
#' fewer than ten leaves an order statistic of a handful of surviving draws
#' rather than a quantile, so it is not reported at all. A 95 percent credible
#' interval therefore needs 400 effective draws and a median 20.
#'
#' @keywords internal
#' @noRd
quantile_is_resolved <- function(n_effective, probs) {
  n_effective * pmin(probs, 1 - probs) >= 10
}


#' Whether the sensitivity analysis is still doubly robust
#'
#' Two separate things can leave it not being so, and only one of them is what
#' [run_tilting()] reports. The coupling of the reweighted draws can fail to
#' meet the moment condition (3.4), which is its own `converged`. Or the
#' reweighting that came first can have spent the posterior, in which case the
#' coupling is being run on a handful of distinct draws and meets the condition
#' on them; the resampling and the jitter of the sweep hide that from every
#' diagnostic the coupling computes, because they leave the particles distinct
#' again. The effective sample size of the sensitivity weights is the only
#' quantity that sees it, so the verdict is the two together, held to the floor
#' `ess_frac` of [drbayes_control()] that the tilting is held to.
#'
#' @param smc The diagnostics of [run_tilting()].
#' @param ess Effective sample size of the sensitivity weights.
#' @param n_draws Number of posterior draws they were computed on.
#' @param ess_frac The floor, as a fraction of `n_draws`.
#'
#' @return `smc`, with `converged` replaced by the joint verdict and the
#'   coupling's own verdict, the effective sample size and the floor added.
#'
#' @keywords internal
#' @noRd
sensitivity_verdict <- function(smc, ess, n_draws, ess_frac) {
  smc$converged_coupling <- isTRUE(smc$converged)
  smc$ess_sensitivity    <- ess
  smc$ess_frac           <- ess_frac
  smc$ess_floor          <- ess_frac * n_draws
  smc$converged          <- smc$converged_coupling && ess >= smc$ess_floor
  smc
}


#' Whether the weights have degenerated, and what to do about it
#'
#' Shared by the warning the analysis issues, the note its methods print and
#' the verdict of [sensitivity_verdict()], so that the threshold and the advice
#' cannot drift apart.
#'
#' @param ess_frac The floor, as a fraction of `n_draws`.
#'
#' @return NULL when the weights are usable, otherwise the message.
#'
#' @keywords internal
#' @noRd
sensitivity_ess_message <- function(ess, n_draws, ess_frac = 0.1) {
  if (ess >= ess_frac * n_draws) {
    return(NULL)
  }
  paste0("The sensitivity weights have degenerated: their effective sample ",
         "size is ", round(ess), " of ", n_draws, " draws, below the floor of ",
         "ess_frac = ", format(ess_frac), " that drbayes_control() sets, at ",
         "which an importance reweighted posterior stops carrying useful ",
         "information. The reweighted summaries rest on those few draws ",
         "alone, and the analysis is reported as not doubly robust for that ",
         "reason whatever the coupling that follows makes of the moment ",
         "condition. The reweighted posterior sits about E[xi] away from the ",
         "original one, so this happens once that distance is more than a few ",
         "posterior standard deviations of the effect. Refit drbayes_pc() ",
         "with more draws, through a larger mc or thin = 1, which buys ",
         "effective sample size in proportion; assume a sensitivity ",
         "distribution with less mass far from zero; or trace the analysis ",
         "over a sequence of smaller xi and read off where the conclusion ",
         "changes.")
}


#' Whether the shift left the range of the outcome, and what to do about it
#'
#' @param n_saturated Pairs whose shifted mean was held at the boundary.
#' @param n_cells Pairs in all.
#'
#' @return NULL when nothing was held there, otherwise the message.
#'
#' @keywords internal
#' @noRd
sensitivity_saturation_message <- function(n_saturated, n_cells) {
  if (n_saturated == 0L) {
    return(NULL)
  }
  paste0("The shifted mean m_A(X; beta) + xi leaves [0, 1] for ", n_saturated,
         " of the ", n_cells, " (treated observation, posterior draw) pairs, ",
         "where a binary outcome model does not exist. It is held at the ",
         "boundary there, so those pairs carry a bias smaller than xi and the ",
         "analysis understates the sensitivity it was asked for. A fitted ",
         "risk of p has room for a risk difference between -p and 1 - p and ",
         "no more: give g() a support the fitted probabilities of the treated ",
         "leave room for.")
}


#' Whether the integral over g() is resolved, and what to do about it
#'
#' Ten effective draws leave the integral with about a third of its own value
#' as Monte Carlo error, which is the point at which the number of draws is
#' worth saying out loud. A `g()` concentrated enough that every draw of it
#' counts, a point mass above all, is resolved by as many draws as were asked
#' for and is not a case of draws going to waste, so the threshold never
#' exceeds `M`.
#'
#' What raising `M` buys is set out under the `M` argument: the effective count
#' grows in proportion to it and the error as its square root, while the
#' running time grows in proportion too, so the message quotes the factor
#' rather than saying "raise M" and leaving the reader to find out.
#'
#' @param xi_ess The `median` and `min` effective number of draws from g().
#'
#' @return NULL when the integral is resolved, otherwise the message.
#'
#' @keywords internal
#' @noRd
sensitivity_xi_message <- function(xi_ess, M) {
  if (xi_ess[["median"]] >= min(10, M)) {
    return(NULL)
  }
  error  <- 100 / sqrt(xi_ess[["median"]])
  needed <- ceiling(M * 10 / xi_ess[["median"]])
  paste0("The Monte Carlo integral over the sensitivity distribution rests on ",
         "a median of ", format(xi_ess[["median"]], digits = 3), " of the M = ",
         M, " draws taken, so each weight carries roughly ", round(error),
         "% Monte Carlo error. The integrand is far narrower than g() is, and ",
         "narrows further as the sample grows, so the share of the draws that ",
         "counts stays about where it is and only their number can be raised: ",
         "the effective count and the running time both grow in proportion to ",
         "M, and the error as its square root. M = ", needed, " would bring ",
         "this integral to ten effective draws and about 32% error, at ",
         format(round(needed / M, 1)), " times the running time of the ",
         "weights. A g() with less mass far from zero costs nothing instead.")
}


#' @param ess_frac The floor the effective sample size is held to.
#'
#' @keywords internal
#' @noRd
print_sensitivity_ess <- function(ess, n_draws, ess_frac = 0.1) {
  percentage <- formatC(100 * ess / n_draws, format = "f", digits = 1)
  cat("Effective sample size: ", round(ess), " of ", n_draws, " draws (",
      percentage, "%)\n", sep = "")

  degenerate <- sensitivity_ess_message(ess, n_draws, ess_frac)
  if (!is.null(degenerate)) {
    cat("\n*** THE SENSITIVITY WEIGHTS HAVE DEGENERATED ***\n")
    cat_wrapped(degenerate, prefix = "  ")
  }
  invisible(NULL)
}


#' The coupling that followed the reweighting, and the verdict on both halves
#'
#' The moment condition and the verdict are separate lines because they are
#' separate questions: the first is about the coupling alone, the second about
#' the reweighting that came first as well. See [sensitivity_verdict()].
#'
#' @keywords internal
#' @noRd
print_sensitivity_coupling <- function(smc, digits) {
  if (is.null(smc)) {
    return(invisible(NULL))
  }
  met <- if (isTRUE(smc$converged_coupling %||% smc$converged)) {
    "satisfied"
  } else {
    "NOT SATISFIED"
  }
  cat("Coupled again at:      lambda = ", format(smc$lambda, digits = digits),
      "\n", sep = "")
  cat("Moment condition:      |mean B_n| = ",
      format(abs(smc$B_mean), digits = digits), " against a tolerance of ",
      format(smc$tol, digits = digits), ", ", met, "\n", sep = "")
  cat("Still doubly robust:   ", if (isTRUE(smc$converged)) "yes" else "NO",
      "\n", sep = "")

  if (!isTRUE(smc$converged)) {
    cat_wrapped(sensitivity_verdict_note(smc), prefix = "  ")
  }
  invisible(NULL)
}


#' Why a sensitivity analysis is not doubly robust
#'
#' @keywords internal
#' @noRd
sensitivity_verdict_note <- function(smc) {
  if (!isTRUE(smc$converged_coupling %||% smc$converged)) {
    return(paste("The coupling of the reweighted draws did not reach the",
                 "moment condition (3.4), so this sensitivity analysis is not",
                 "doubly robust. Inspect $smc."))
  }
  paste0("The coupling of the reweighted draws met the moment condition ",
         "(3.4), but it was run on a cloud carrying ",
         round(smc$ess_sensitivity), " effective draws, below the floor of ",
         round(smc$ess_floor), ". The resampling makes the particles distinct ",
         "again, so neither the moment condition nor any diagnostic of the ",
         "coupling can see how few of the original draws are behind them: ",
         "the condition holds on a posterior that has been spent, which is ",
         "not the same thing as the analysis being doubly robust.")
}


#' @keywords internal
#' @noRd
print_sensitivity_saturation <- function(n_saturated) {
  if (is.null(n_saturated) || n_saturated[["saturated"]] == 0L) {
    return(invisible(NULL))
  }
  cat("Shift out of range:    ", n_saturated[["saturated"]], " of ",
      n_saturated[["of"]], " fitted means held at 0 or 1\n", sep = "")
  cat_wrapped(sensitivity_saturation_message(n_saturated[["saturated"]],
                                             n_saturated[["of"]]),
              prefix = "  ")
  invisible(NULL)
}


#' @keywords internal
#' @noRd
print_sensitivity_sigma <- function(sigma, digits) {
  if (length(sigma) == 1L && is.na(sigma)) {
    return(invisible(NULL))
  }
  described <- if (length(sigma) == 1L) {
    paste0("fixed at ", format(sigma, digits = digits))
  } else {
    paste0("drawn per posterior draw: mean ",
           format(mean(sigma), digits = digits), ", sd ",
           format(stats::sd(sigma), digits = digits))
  }
  cat("Residual scale sigma:  ", described, "\n", sep = "")
  invisible(NULL)
}


#' @keywords internal
#' @noRd
print_sensitivity_interval_note <- function(ess) {
  if (all(quantile_is_resolved(ess, c(0.025, 0.975)))) {
    return(invisible(NULL))
  }
  cat_wrapped(paste0("The credible interval of the sensitivity row is not ",
                     "reported: an effective sample size of ", round(ess),
                     " puts ", format(round(ess * 0.025, 1)), " effective ",
                     "draws beyond each of its ends, fewer than the ten a ",
                     "quantile of them needs, so what is there is the most ",
                     "extreme handful of draws that survived the reweighting ",
                     "and not a 2.5 percent point. Report the mean and the ",
                     "standard deviation, or reweight more draws."))
  invisible(NULL)
}


#' @keywords internal
#' @noRd
describe_sensitivity_method <- function(method) {
  switch(method,
         "common" = "one xi shared by every observation (Algorithm A.1)",
         "per-observation" = "one xi per observation (Algorithm 3)",
         method)
}


#' @keywords internal
#' @noRd
describe_sensitivity_xi <- function(x, digits) {
  values <- as.numeric(x$xi)
  fmt    <- function(value) format(value, digits = digits)
  drawn  <- if (is.matrix(x$xi)) {
    paste0(x$M, " draws for each of ", nrow(x$xi), " treated observations")
  } else {
    paste0(x$M, " draws")
  }
  paste0(drawn, " from g(): mean ", fmt(mean(values)), ", range ",
         fmt(min(values)), " to ", fmt(max(values)))
}


#' @keywords internal
#' @noRd
describe_sensitivity_atoms <- function(x, digits) {
  paste0("integral resolved by a median of ",
         format(x$xi_ess[["median"]], digits = digits), " of them, ",
         format(x$xi_ess[["min"]], digits = digits), " at worst")
}
