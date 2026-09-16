#' Confounder Selection with Posterior Coupling
#'
#' Implements Algorithm 4 of Orihara, Momozaki and Sugasawa (2025). Both models
#' are fitted under shrinkage priors, the covariates whose propensity score
#' coefficient survives that shrinkage are selected, and the two posteriors are
#' then coupled through a moment condition restricted to the selected set.
#'
#' Applying shrinkage separately to the outcome and the propensity score model
#' is the obvious thing to do and it is what causes the problem this function
#' solves: a weak confounder is shrunk away in both models at once, and the
#' resulting regularization-induced confounding biases the treatment effect
#' (Hahn et al., 2018). Selecting on the treatment model and then coupling
#' gives such a covariate a second chance, because the tilting updates its
#' outcome coefficient jointly with the propensity score.
#'
#' @param outcome.formula Formula for the outcome model, for example
#'   \code{Y ~ A + X1 + X2 + A:X1}, with the response on the left and the
#'   treatment among the predictors on the right.
#' @param ps.formula Formula for the propensity score model, for example
#'   \code{A ~ X1 + X2 + X3}, with the treatment on the left and the candidate
#'   confounders on the right. These candidates are what the selection chooses
#'   from.
#' @param data Data frame holding every variable used by either formula.
#' @param threshold Numeric. A covariate is selected when the posterior mean of
#'   its propensity score coefficient is at least \code{threshold} in absolute
#'   value, that is \eqn{S = \{j : |\bar{\alpha}_j| \ge}{S = {j : |alpha-bar_j|
#'   >=}} \code{threshold}\eqn{\}}{}. Default 0.01, the value used in the right
#'   heart catheterization analysis of the paper. Zero selects every candidate.
#'   The scale is that of the propensity score linear predictor, so it depends
#'   on how the covariates are scaled.
#' @param family Character string for the outcome model family, either
#'   \code{"gaussian"} or \code{"binomial"}. NULL, the default, detects a
#'   binary response as binomial and anything else as gaussian.
#' @param link Character string for the link of the OUTCOME model:
#'   \code{"identity"}, \code{"logit"} or \code{"probit"}. NULL, the default,
#'   takes the canonical link of \code{family}.
#' @param ps.link Character string for the link of the PROPENSITY SCORE model,
#'   either \code{"logit"} (default) or \code{"probit"}. The moment condition
#'   turns the propensity score draws into probabilities with this link.
#' @param outcome.prior,ps.prior Character. Prior for the coefficients of each
#'   model in Step 1, either \code{"horseshoe"} (default) or \code{"normal"}.
#'   Step 1 of Algorithm 4 asks for a shrinkage prior, "such as the horseshoe
#'   prior", which is why the horseshoe is the default rather than the only
#'   option. The selection reads the propensity score coefficients, so
#'   \code{ps.prior = "normal"} shrinks nothing and leaves the threshold
#'   selecting on unregularised posterior means.
#' @param mc,bn,thin,chains Integers. Number of MCMC iterations per chain, the
#'   burn-in to discard, the thinning interval, and the number of chains. They
#'   also live on the \code{control} object; naming one here overrides it.
#' @param method Character. How Step 2 finds the tilting parameter.
#'   \code{"smc"} is Algorithm 2, which walks a grid of tilting parameters and
#'   rejuvenates the particles at each step. \code{"is"} is Algorithm 1, a
#'   single importance sampling step, which is faster but degenerates when the
#'   two posteriors put little mass where the moment condition holds. Step 2 of
#'   Algorithm 4 admits either. Default \code{"smc"}.
#' @param rejuvenate Character, \code{"all"} (default) or \code{"selected"},
#'   choosing which coefficients the Liu and West (2001) kernel refreshes at
#'   each step of the sweep. See "Rejuvenating the carried block" below; the
#'   default departs deliberately from the literal wording of Algorithm 4, and
#'   \code{"selected"} reproduces it. Ignored when \code{method = "is"}, which
#'   moves no coefficient at all and records this argument as NA.
#' @param control List of tuning parameters from \code{\link{drbayes_control}}.
#'   Its \code{n_steps} is the one to think about here: it sets how many times
#'   the sweep resamples, and under \code{rejuvenate = "selected"} the carried
#'   block loses ancestry at every one of those resampling steps. The default
#'   of 1000 is a reasonable choice under the default \code{rejuvenate}, which
#'   refreshes the whole particle and so has no ancestry to lose.
#' @param diagnostics Character. What to do when the Step 1 draws fail the
#'   convergence criteria: \code{"warn"} (default), \code{"error"}, or
#'   \code{"none"} to skip the check.
#' @param outcome.model,ps.model Functions fitting the Step 1 models, for
#'   callers who want a sampler the package does not provide, such as
#'   \code{\link{bayes_stan}} or a wrapper around draws obtained elsewhere.
#'   Each should take \code{Y}, \code{X}, \code{mc} and \code{chains} and
#'   return an iterations by chains by parameters array, the layout the
#'   built-in samplers return. Given, they override \code{outcome.prior} and
#'   \code{ps.prior}, and the fit records the prior it did not use as
#'   \code{"supplied"}. Default NULL, which selects the built-in sampler for
#'   the link and the prior. There is no \code{outcome.samples} argument here
#'   as there is in \code{\link{drbayes_pc}}: draws already in hand are passed
#'   as \code{function(Y, X, mc, chains) draws}.
#' @param outcome.priors,ps.priors Lists of further arguments for the Step 1
#'   samplers, for example \code{list(p0 = 10)} to change the horseshoe's guess
#'   at how many coefficients are non-zero. The treatment columns of the
#'   outcome model are exempted from shrinkage automatically; passing
#'   \code{unshrunk} yourself replaces that choice.
#' @param seed Integer or NULL. When given, the fit is reproducible and the
#'   caller's random number stream is restored afterwards.
#' @param verbose Logical. Whether to report progress, including how many
#'   covariates were selected, through \code{\link[base]{message}}. Default
#'   TRUE.
#' @param na.action One of \code{"na.omit"} (default), \code{"na.fail"} or
#'   \code{"na.exclude"}.
#'
#' @return An object of class \code{"DRBayes_select"}, which inherits from
#'   \code{"DRBayes"} and therefore works with \code{\link{summary.DRBayes}}
#'   and \code{\link{plot.DRBayes}}. It holds everything
#'   \code{\link{drbayes_pc}} returns, and in addition:
#' \describe{
#'   \item{selected}{Character vector of the selected covariates, named as the
#'     columns of the propensity score design matrix, so a factor appears once
#'     per dummy column.}
#'   \item{threshold}{The threshold used.}
#'   \item{n_candidates}{Number of candidate covariates the selection chose
#'     from, again counting dummy columns separately.}
#'   \item{posterior_means}{Named numeric vector of the Step 1 posterior means
#'     \eqn{\bar{\alpha}_j}{alpha-bar_j} the selection was made on, one per
#'     candidate. The intercept is not a candidate and is not included.}
#'   \item{method, rejuvenate}{Which algorithm ran Step 2, and which
#'     coefficients its kernel refreshed. \code{rejuvenate} is NA under
#'     \code{method = "is"}, which applies no kernel and so refreshes none.}
#' }
#'   Its \code{smc} component carries one diagnostic beyond the ones
#'   \code{\link{drbayes_pc}} records: \code{n_ancestors_cumulative}, the
#'   number of distinct Step 1 draws the returned particles still descend from.
#'   Unlike \code{n_ancestors_min}, which counts the survivors of the worst
#'   single resampling step, this accumulates over the whole sweep, and it is
#'   the count that describes a carried block. It is reported whenever some
#'   coordinate of a particle is a verbatim copy of a Step 1 draw: under
#'   \code{method = "is"}, under \code{rejuvenate = "selected"} when the
#'   selection left a coefficient outside the updated block, and when the
#'   constraint already held so that no sweep ran. It is NA otherwise.
#'
#'   \code{smc$converged} means here what it means in
#'   \code{\link{drbayes_pc}}, that the moment condition was met on draws
#'   numerous enough to describe a posterior, and it is FALSE when the carried
#'   block has collapsed however well the constraint is satisfied. The
#'   \code{print} method separates the two cases; \code{summary()} reports only
#'   the flag.
#'
#' @details
#' \strong{The algorithm:}
#'
#' \emph{Step 1} generates posterior samples for both models under shrinkage
#' priors (Carvalho et al., 2010) and forms the posterior means
#' \eqn{\bar{\alpha}}{alpha-bar} of the propensity score coefficients. The
#' selected set is
#' \eqn{S = \{j : |\bar{\alpha}_j| \ge \tau\}}{S = {j : |alpha-bar_j| >= tau}}
#' for the threshold \eqn{\tau}{tau}.
#'
#' \emph{Step 2} replaces the moment condition (3.4) with
#' \deqn{B_n^S(\alpha_S, \beta_S) = \frac{1}{n} \sum_{i=1}^{n}
#'   \frac{A_i - e(X_{Si}; \alpha_S)}{e(X_{Si}; \alpha_S)(1 -
#'   e(X_{Si}; \alpha_S))} (Y_i - m_{A_i}(X_i; \beta_S, \beta_{S^c}))}
#' and runs Algorithm 1 or Algorithm 2 against it. The propensity score is
#' evaluated on the selected covariates alone, while the outcome mean keeps
#' every covariate.
#'
#' \strong{Three points the algorithm turns on:}
#'
#' The selection is made on the coefficients of the PROPENSITY SCORE model, not
#' of the outcome model. That is the whole point: a covariate that predicts
#' treatment is a candidate confounder however weakly it predicts the outcome,
#' and it is exactly such a covariate that separate shrinkage discards.
#'
#' The outcome coefficients outside the selected set,
#' \eqn{\beta_{S^c}}{beta_Sc}, are not dropped, and the g-computation is done
#' with the full \eqn{\beta}{beta}. Using the selected block alone would
#' silently change the estimand into the average treatment effect of a
#' different, smaller outcome model.
#'
#' The treatment is never a candidate for selection and is never shrunk. Its
#' main effect and every interaction with it are exempted from the shrinkage in
#' Step 1 and are always part of the block the sweep updates, since shrinking
#' or freezing them would attenuate the very effect being estimated.
#'
#' \strong{Which outcome coefficients a selected covariate reaches:}
#'
#' A covariate is selected as a column of the propensity score design matrix,
#' and its outcome coefficients are found through the terms structure of the
#' two formulas rather than by matching column names. The propensity score
#' column is traced back to the term it came from, that term to the variables
#' it is built from, and those variables forward to every outcome column they
#' contribute to. So selecting \code{X1} in \code{A ~ X1 + X2} reaches
#' \code{log(X1)}, \code{poly(X1, 2)} and \code{A:X1} in the outcome model,
#' which matching names would not. A selected covariate that reaches no outcome
#' column at all is reported: the outcome model does not adjust for it, so the
#' coupling has no coefficient of it to revive.
#'
#' \strong{Rejuvenating the carried block:}
#'
#' Algorithm 4 says that in Step 2 "only \eqn{(\alpha_S, \beta_S)}{(alpha_S,
#' beta_S)} are updated with \eqn{\beta_{S^c}}{beta_Sc} being unchanged". Taken
#' literally as an instruction to the sequential Monte Carlo sweep, that means
#' applying the Liu and West (2001) kernel to the selected block only. The
#' sweep still resamples whole particles, so the carried block is never
#' refreshed and its ancestry collapses geometrically. Over the thousand steps
#' \code{\link{drbayes_control}} takes by default, a cloud of several thousand
#' particles can come back tracing to a dozen distinct Step 1 draws, and on
#' some data to one, leaving \eqn{\beta_{S^c}}{beta_Sc} a point mass reported
#' as a posterior.
#'
#' The default \code{rejuvenate = "all"} therefore applies the kernel to every
#' coefficient, and departs from the literal wording on purpose. The kernel is
#' a rejuvenation step: it shrinks each particle towards the cloud mean and
#' adds back the variance that shrinking removed, so it preserves the first two
#' moments of the cloud and moves no parameter towards the constraint. What
#' does move a parameter towards the constraint is the importance weight, and
#' \eqn{\beta_{S^c}}{beta_Sc} enters those weights whatever this argument is
#' set to, because it enters \eqn{B_n^S}{B_n^S} through
#' \eqn{m_{A_i}(X_i; \beta_S, \beta_{S^c})}{m}. Rejuvenating the carried block
#' therefore does not update it in the sense Step 2 excludes, while freezing it
#' destroys the posterior it is there to represent. The propensity score
#' coefficients outside \eqn{S}{S} are refreshed on the same grounds, and with
#' more reason: they do not enter \eqn{B_n^S}{B_n^S} at all, so the tilting
#' leaves their conditional distribution alone and a frozen copy of it is the
#' one thing that would not.
#'
#' \code{rejuvenate = "selected"} reproduces the literal reading for anyone who
#' wants it. The kernel then refreshes the selected block alone, so both
#' \eqn{\beta_{S^c}}{beta_Sc} and \eqn{\alpha_{S^c}}{alpha_Sc} are carried
#' verbatim and both collapse together: the propensity score draws in
#' \code{fit$particles$ps} are as much a point mass as the outcome ones, and
#' the warning names how many coefficients of each are affected. It reports
#' \code{fit$smc$n_ancestors_cumulative}, and reports the fit as not converged
#' once the carried block has collapsed onto few enough draws that its
#' posterior spread is not to be believed. That floor is \code{ess_frac} of
#' \code{\link{drbayes_control}}, the same one the tilting holds its own counts
#' to. Fewer steps in the tilting parameter, through
#' \code{drbayes_control(n_steps = )}, mean fewer resampling rounds and so a
#' slower collapse, but they do not remove it.
#'
#' @references
#' Orihara, S., Momozaki, T., & Sugasawa, S. (2025). Bayesian Doubly Robust
#' Causal Inference via Posterior Coupling. arXiv preprint arXiv:2506.04868.
#'
#' Carvalho, C. M., Polson, N. G., & Scott, J. G. (2010). The horseshoe
#' estimator for sparse signals. Biometrika, 97(2), 465-480.
#'
#' Hahn, P. R., Carvalho, C. M., Puelz, D., & He, J. (2018). Regularization and
#' confounding in linear regression for treatment effect estimation. Bayesian
#' Analysis, 13(1), 163-182.
#'
#' Liu, J., & West, M. (2001). Combined parameter and state estimation in
#' simulation-based filtering. In A. Doucet, N. de Freitas, & N. Gordon (Eds.),
#' Sequential Monte Carlo Methods in Practice (pp. 197-223). Springer.
#'
#' Ning, Y., Sida, P., & Imai, K. (2020). Robust estimation of causal effects
#' via a high-dimensional covariate balancing propensity score. Biometrika,
#' 107(3), 533-554.
#'
#' @examples
#' set.seed(1)
#' n <- 150
#' dat <- data.frame(matrix(rnorm(n * 6), n, 6))
#' names(dat) <- paste0("X", 1:6)
#'
#' # X1 predicts treatment strongly and the outcome only weakly, so separate
#' # shrinkage discards it and the treatment effect absorbs its confounding.
#' dat$A <- rbinom(n, 1, plogis(1.5 * dat$X1 + 0.8 * dat$X2))
#' dat$Y <- 2 * dat$A + 0.5 * dat$X1 + 1.2 * dat$X2 + rnorm(n, sd = 3)
#'
#' covariates <- paste(paste0("X", 1:6), collapse = " + ")
#' outcome <- as.formula(paste("Y ~ A +", covariates))
#' ps      <- as.formula(paste("A ~", covariates))
#'
#' fit <- drbayes_select(
#'   outcome.formula = outcome,
#'   ps.formula      = ps,
#'   data            = dat,
#'   threshold       = 0.1,
#'   mc = 800, bn = 200, chains = 2L, seed = 2, verbose = FALSE,
#'   control = drbayes_control(n_steps = 150)
#' )
#' fit
#' fit$selected  # the two covariates that drive treatment
#'
#' # The literal reading of Step 2, which carries the coefficients outside the
#' # selected set at their Step 1 draws. Its ancestry count reports how many
#' # of those draws that block still holds, and fewer steps in the tilting
#' # parameter leave it more of them. The warning it raises is the point of
#' # the comparison, not a problem with the data.
#' literal <- drbayes_select(
#'   outcome.formula = outcome,
#'   ps.formula      = ps,
#'   data            = dat,
#'   threshold       = 0.1,
#'   rejuvenate      = "selected",
#'   mc = 800, bn = 200, chains = 2L, seed = 2, verbose = FALSE,
#'   control = drbayes_control(n_steps = 150)
#' )
#' literal$smc$n_ancestors_cumulative
#'
#' @seealso \code{\link{drbayes_pc}} for posterior coupling without selection,
#'   \code{\link{drbayes_control}} for the tuning parameters, and
#'   \code{\link{bayes_lm_hs}} for the horseshoe samplers Step 1 uses.
#'
#' @export
drbayes_select <- function(outcome.formula, ps.formula, data,
                           threshold = 0.01,
                           family = NULL, link = NULL, ps.link = "logit",
                           outcome.prior = c("horseshoe", "normal"),
                           ps.prior      = c("horseshoe", "normal"),
                           mc = 5000, bn = 1000, thin = 2, chains = 4L,
                           method = c("smc", "is"),
                           rejuvenate = c("all", "selected"),
                           control = drbayes_control(),
                           diagnostics = c("warn", "error", "none"),
                           outcome.model = NULL, ps.model = NULL,
                           outcome.priors = list(), ps.priors = list(),
                           seed = NULL, verbose = TRUE,
                           na.action = "na.omit") {

  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE")
  }
  say <- function(...) if (verbose) message(...)

  # A seed here is local, so a script that seeded itself at the top keeps
  # control of its own stream.
  if (!is.null(seed)) {
    cl <- match.call()
    cl$seed <- NULL
    return(with_preserved_rng(seed, eval(cl, parent.frame())))
  }

  call_info     <- match.call()
  threshold     <- validate_threshold(threshold)
  outcome.prior <- match.arg(outcome.prior)
  ps.prior      <- match.arg(ps.prior)
  method        <- match.arg(method)
  rejuvenate    <- match.arg(rejuvenate)
  diagnostics   <- match.arg(diagnostics)

  # The four MCMC settings are also on the control object; naming one directly
  # is the common case and reads better than wrapping it.
  control  <- resolve_control(control)
  supplied <- names(as.list(match.call())[-1L])
  for (nm in intersect(c("mc", "bn", "thin", "chains"), supplied)) {
    control[[nm]] <- get(nm)
  }
  control <- resolve_control(control)

  md    <- prepare_model_data(outcome.formula, ps.formula, data, na.action)
  Y     <- md$Y
  A     <- md$A
  Z.lm  <- md$Z.lm
  Z.ps  <- md$Z.ps

  # Selection needs candidates on the treatment model and something outside the
  # treatment for the outcome shrinkage to work on. Both are cheaper to refuse
  # here than to let the sampler discover thousands of iterations later.
  if (ncol(md$X.ps) == 0L) {
    stop("ps.formula has no covariates, so there is nothing to select from. ",
         "Put the candidate confounders on the right hand side of ",
         "ps.formula, as in ", md$treatment, " ~ X1 + X2.")
  }
  if (length(md$treatment_columns) >= ncol(md$X.lm)) {
    stop("Every column of the outcome model comes from the treatment, so the ",
         "shrinkage prior of Step 1 has nothing to shrink. Add the candidate ",
         "confounders to outcome.formula, or use drbayes_pc() if no ",
         "selection is wanted.")
  }

  say("Formula processing complete:")
  say("  Observations used: ", length(Y), " out of ", nrow(data))
  say("  Treatment group: ", sum(A), " | Control group: ", sum(1 - A))
  say("  Candidate covariates: ", ncol(md$X.ps))

  links   <- resolve_links(family, link, ps.link, Y)
  family  <- links$family
  link    <- links$link
  ps.link <- links$ps.link

  # ==========================================================================
  # Step 1: shrinkage fits of both models
  # ==========================================================================

  check_sampler_link(outcome.model, link, "outcome.model", "link")
  check_sampler_link(ps.model, ps.link, "ps.model", "ps.link")
  outcome_sampler <- outcome.model %||% select_sampler(link, outcome.prior)
  ps_sampler      <- ps.model %||% select_sampler(ps.link, ps.prior)

  say("Step 1: sampling the outcome model with ",
      step1_description(outcome.model, link, outcome.prior),
      " and the propensity score model with ",
      step1_description(ps.model, ps.link, ps.prior), ".")

  outcome.args <- list(Y = Y, X = md$X.lm, mc = control$mc,
                       chains = control$chains)
  if (has_formal(outcome_sampler, "init")) {
    outcome.args$init <- control_init_for(control, "outcome")
  }
  # The treatment main effect and every interaction with it, found through
  # columns_involving() when the design matrices were built. Shrinking them
  # would attenuate the treatment effect.
  if (has_formal(outcome_sampler, "unshrunk") &&
      !("unshrunk" %in% names(outcome.priors))) {
    outcome.args$unshrunk <- md$treatment_columns
  }
  post.otc <- do.call(outcome_sampler, c(outcome.args, outcome.priors))

  ps.args <- list(Y = A, X = md$X.ps, mc = control$mc,
                  chains = control$chains)
  if (has_formal(ps_sampler, "init")) {
    ps.args$init <- control_init_for(control, "ps")
  }
  post.ps <- do.call(ps_sampler, c(ps.args, ps.priors))

  keep      <- seq(control$bn + 1L, control$mc, by = control$thin)
  draws.otc <- post.otc[keep, , , drop = FALSE]
  draws.ps  <- post.ps[keep, , , drop = FALSE]

  diag_table <- check_step1_convergence(draws.otc, draws.ps, control,
                                        diagnostics)

  # The sweep works on pooled draws: the chain structure has served its purpose
  # once convergence has been assessed.
  betas.otc <- matrix(draws.otc, ncol = dim(draws.otc)[3L],
                      dimnames = list(NULL, dimnames(draws.otc)[[3L]]))
  betas.ps  <- matrix(draws.ps, ncol = dim(draws.ps)[3L],
                      dimnames = list(NULL, dimnames(draws.ps)[[3L]]))

  # ==========================================================================
  # Step 2: select, then couple on the selected block alone
  # ==========================================================================

  # The intercept is not a candidate: it is not a covariate, and the propensity
  # score model needs it whatever the selection decides.
  alpha_bar <- colMeans(betas.ps)[-1L]
  selected  <- select_confounders(alpha_bar, threshold)

  say("Step 1 selected ", length(selected), " of ", length(alpha_bar),
      " candidate covariates at |alpha_bar| >= ", format(threshold), ".")

  ps_columns <- 1L + match(selected, names(alpha_bar))
  ps_update  <- c(1L, ps_columns)
  reach      <- selection_columns(ps_columns, md)
  # Intercept and treatment columns first: neither is ever selected away, and
  # freezing the treatment coefficients would leave the tilting no way to move
  # the treatment effect at all.
  otc_update <- sort(unique(c(1L, md$treatment_columns + 1L, reach$columns)))

  warn_unreached(selected[reach$unmatched], length(reach$columns) == 0L,
                 rejuvenate)

  say("Step 2: the moment condition uses ", length(ps_update),
      " propensity score coefficient(s) and every outcome coefficient.")
  if (method == "smc" && rejuvenate == "selected") {
    say("  The kernel refreshes ", length(otc_update), " outcome ",
        "coefficient(s); the remaining ", ncol(betas.otc) - length(otc_update),
        " are frozen at their Step 1 draws.")
  }

  # B_n^S evaluates the propensity score on the selected covariates alone and
  # the outcome mean on every covariate, exactly as the displayed moment
  # condition of Step 2 does. The positivity check inside the sweep names the
  # restricted model, since that is the one whose fitted scores it saw.
  d <- tilting_data(Z.lm = Z.lm,
                    Z.ps = Z.ps[, ps_update, drop = FALSE],
                    A = A, Y = Y,
                    inverse_link    = links$inverse_link,
                    ps_inverse_link = links$ps_inverse_link,
                    ps_formula_text = selected_formula_text(ps_columns, md))

  tilting <- select_tilting(betas.ps, betas.otc, ps_update, otc_update, d,
                            control, method = method, rejuvenate = rejuvenate,
                            tol = control$tol)

  # ==========================================================================
  # Posterior of the average treatment effect
  # ==========================================================================

  # The FULL beta, the updated block together with the carried beta_{S^c}.
  # Averaging the fitted contrast over only the selected coefficients would be
  # the estimand of a different outcome model.
  gcomp_ate <- function(betas) {
    rowMeans(links$inverse_link(tcrossprod(betas, md$Z.lm1))) -
      rowMeans(links$inverse_link(tcrossprod(betas, md$Z.lm0)))
  }

  result <- list(
    g.comp          = gcomp_ate(betas.otc),
    pc              = gcomp_ate(tilting$betas.otc),
    selected        = selected,
    threshold       = threshold,
    n_candidates    = length(alpha_bar),
    posterior_means = alpha_bar,
    family          = family,
    link            = link,
    ps.link         = ps.link,
    # A supplied sampler carries its own prior, so recording the argument that
    # was then never used would misdescribe the fit.
    outcome.prior   = if (is.null(outcome.model)) outcome.prior
                      else "supplied",
    ps.prior        = if (is.null(ps.model)) ps.prior else "supplied",
    method          = method,
    # Algorithm 1 applies no kernel, so there is no block for this argument
    # to name and the value it was called with would misdescribe the fit.
    rejuvenate      = if (method == "is") NA_character_ else rejuvenate,
    control         = control,
    smc             = tilting$smc,
    diagnostics     = diag_table,
    particles       = if (isTRUE(control$keep_particles)) {
      list(outcome = tilting$betas.otc, ps = tilting$betas.ps)
    } else {
      NULL
    },
    call      = call_info,
    data_info = md$data_info
  )

  class(result) <- c("DRBayes_select", "DRBayes", "list")
  result
}


#' Print a confounder selection fit
#'
#' @description
#' Everything \code{\link{print.DRBayes}} shows, followed by which covariates
#' Step 1 of Algorithm 4 selected, which it left out, and what the sweep was
#' allowed to move.
#'
#' @param x An object of class \code{"DRBayes_select"} from
#'   \code{\link{drbayes_select}}.
#' @param digits Number of significant digits. Default 3.
#' @param ... Ignored, present for consistency with the generic.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' set.seed(1)
#' n <- 250
#' dat <- data.frame(X1 = rnorm(n), X2 = rnorm(n), X3 = rnorm(n))
#' dat$A <- rbinom(n, 1, plogis(dat$X1 + 0.5 * dat$X2))
#' dat$Y <- 2 * dat$A + dat$X1 + dat$X2 + rnorm(n)
#' fit <- drbayes_select(Y ~ A + X1 + X2 + X3, A ~ X1 + X2 + X3, dat,
#'                       mc = 1600, bn = 400, seed = 2, verbose = FALSE,
#'                       control = drbayes_control(n_steps = 150))
#' print(fit)
#'
#' @seealso \code{\link{drbayes_select}}, \code{\link{print.DRBayes}}
#' @export
print.DRBayes_select <- function(x, digits = 3, ...) {

  NextMethod()

  n_selected <- length(x$selected)
  dropped    <- setdiff(names(x$posterior_means), x$selected)

  cat("\nConfounder selection (Algorithm 4)\n")
  cat("  Selected ", n_selected, " of ", x$n_candidates,
      " candidate covariate(s) at |alpha_bar| >= ",
      format(x$threshold, digits = digits), "\n", sep = "")
  cat_wrapped(paste0("Selected: ", toString(x$selected), "."), prefix = "    ")
  if (length(dropped) > 0L) {
    cat_wrapped(paste0("Not selected: ", toString(dropped), "."),
                prefix = "    ")
  }

  if (identical(x$method, "is")) {
    cat_wrapped(paste("Step 2 reweighted the draws against the restricted",
                      "moment condition and resampled them once, so every",
                      "coefficient of a particle is a Step 1 draw and",
                      x$smc$n_distinct, "of the", length(x$pc),
                      "particles are distinct."),
                prefix = "  ")
  } else if (identical(x$rejuvenate, "selected") &&
             is.na(x$smc$n_ancestors_cumulative)) {
    cat_wrapped(paste("The tilting updated the selected block only, the",
                      "literal reading of Step 2, but the selection reached",
                      "every coefficient, so there was nothing outside that",
                      "block to carry and the sweep refreshed the whole",
                      "particle."),
                prefix = "  ")
  } else if (identical(x$rejuvenate, "selected")) {
    cat_wrapped(paste("The tilting updated the selected block only, the",
                      "literal reading of Step 2. The coefficients outside",
                      "that block, in both models, were carried with their",
                      "particle. The outcome ones are still in the",
                      "g-computation, so the estimand is unchanged, but every",
                      "carried coefficient traces back to only",
                      x$smc$n_ancestors_cumulative,
                      "distinct Step 1 draw(s) and its",
                      "posterior spread is understated accordingly."),
                prefix = "  ")
    # The convergence line above reads a fit as having failed the moment
    # condition whenever it is reported as not converged, so a fit stopped by
    # the collapse alone needs its own numbers stating.
    if (carried_collapse(x$smc$n_ancestors_cumulative, length(x$pc),
                         x$control$ess_frac) &&
        abs(x$smc$B_mean) < x$smc$tol) {
      cat_wrapped(paste0("That collapse, and not the moment condition, is why",
                         " this fit is reported as not converged: |mean B_n|",
                         " = ", format(abs(x$smc$B_mean), digits = digits),
                         " is within the tolerance ",
                         format(x$smc$tol, digits = digits),
                         ". Use rejuvenate = \"all\", or take fewer steps",
                         " with drbayes_control(n_steps = )."),
                  prefix = "  ")
    }
  } else {
    cat_wrapped(paste("The moment condition was restricted to the selected",
                      "covariates and the particles were rejuvenated as a",
                      "whole. Every outcome coefficient is in the",
                      "g-computation, so the estimand is unchanged."),
                prefix = "  ")
  }

  invisible(x)
}


#' Covariates whose propensity score coefficient survives the threshold
#'
#' @param alpha_bar Named posterior means of the propensity score coefficients,
#'   without the intercept.
#' @param threshold Smallest absolute posterior mean that is selected.
#'
#' @return A character vector of covariate names, in design matrix order.
#'
#' @keywords internal
#' @noRd
select_confounders <- function(alpha_bar, threshold) {
  keep <- abs(alpha_bar) >= threshold

  if (!any(keep)) {
    largest <- which.max(abs(alpha_bar))
    stop("No covariate reaches the selection threshold of ", format(threshold),
         ", so there is nothing left to couple on. The largest posterior ",
         "mean is ", format(abs(alpha_bar[largest]), digits = 3),
         " in absolute value, for ", names(alpha_bar)[largest],
         ". Lower threshold below that value, or use drbayes_pc() if no ",
         "selection is wanted.")
  }

  names(alpha_bar)[keep]
}


#' Outcome model columns the selected covariates reach
#'
#' A selected covariate is a column of the propensity score design matrix, and
#' the outcome coefficients it owns are found through the terms structure of
#' the two formulas: the column is traced back to its term, the term to the
#' variables it is built from, and those variables forward to every outcome
#' column they contribute to. Matching design matrix column names instead would
#' miss `log(X1)`, `poly(X1, 2)`, `scale(X1)` and any factor coded differently
#' between the two formulas, and freezing those is exactly what Algorithm 4
#' exists to avoid. This is the mapping [columns_involving()] makes for the
#' treatment, generalised to a set of covariates.
#'
#' An outcome column counts as reached when it shares any variable with the
#' selected term, so that a selected interaction reaches the columns of both
#' its parts. Reaching too widely only means the sweep updates a coefficient it
#' could have left alone; reaching too narrowly freezes a confounder.
#'
#' @param ps_columns Integer column indices into `md$Z.ps` of the selected
#'   covariates, the intercept excluded.
#' @param md The object returned by [prepare_model_data()].
#'
#' @return A list with `columns`, the integer column indices into `md$Z.lm`
#'   that the selection reaches, and `unmatched`, a logical vector along
#'   `ps_columns` marking the selected covariates that reach none of them.
#'
#' @keywords internal
#' @noRd
selection_columns <- function(ps_columns, md) {
  # The formulas stored by prepare_model_data() have already had any "." in
  # them expanded, so their terms have the same terms and the same order as the
  # ones the design matrices were built from.
  ps_terms  <- stats::terms(md$data_info$ps_formula)
  otc_terms <- stats::terms(md$data_info$outcome_formula)

  ps_variables  <- term_variables(ps_terms)
  otc_variables <- term_variables(otc_terms)
  ps_assign     <- attr(md$Z.ps, "assign")
  otc_assign    <- attr(md$Z.lm, "assign")

  reached <- lapply(ps_columns, function(column) {
    wanted <- ps_variables[[ps_assign[column]]]
    shares <- vapply(otc_variables, function(v) any(v %in% wanted),
                     logical(1))
    which(otc_assign %in% which(shares))
  })

  list(columns   = sort(unique(as.integer(unlist(reached)))),
       unmatched = lengths(reached) == 0L)
}


#' The variables each term of a model formula is built from
#'
#' A term is named by the expressions it was written with, such as
#' `log(X1):X2`, but what matters here is the underlying data variables, `X1`
#' and `X2`. Those are what let a covariate written one way in one formula be
#' recognised in the other.
#'
#' @param terms_obj A terms object.
#'
#' @return A list with one character vector of variable names per term of
#'   `terms_obj`, in the order the term labels are in.
#'
#' @keywords internal
#' @noRd
term_variables <- function(terms_obj) {
  factors <- attr(terms_obj, "factors")
  if (is.null(factors) || length(factors) == 0L) {
    return(list())
  }
  # The rows of the factors matrix are the variables of the terms object, in
  # the order they appear there, so they can be indexed positionally.
  variables <- as.list(attr(terms_obj, "variables"))[-1L]

  out <- lapply(seq_len(ncol(factors)), function(k) {
    unique(unlist(lapply(variables[factors[, k] != 0], all.vars)))
  })
  names(out) <- colnames(factors)
  out
}


#' The propensity score model the moment condition actually evaluates
#'
#' \eqn{B_n^S} scores the treatment on the selected covariates alone, so a
#' positivity failure inside the sweep belongs to that restricted model and not
#' to the propensity score formula the user wrote.
#'
#' @param ps_columns Integer column indices into `md$Z.ps` of the selected
#'   covariates, the intercept excluded.
#' @param md The object returned by [prepare_model_data()].
#'
#' @return A one line character description of the restricted formula.
#'
#' @keywords internal
#' @noRd
selected_formula_text <- function(ps_columns, md) {
  ps_terms <- stats::terms(md$data_info$ps_formula)
  labels   <- attr(ps_terms, "term.labels")
  used     <- unique(labels[attr(md$Z.ps, "assign")[ps_columns]])
  paste(md$treatment, "~", paste(used, collapse = " + "))
}


#' Report selected covariates the outcome model knows nothing about
#'
#' @param unreached Names of the selected covariates that contribute no column
#'   to the outcome design matrix.
#' @param none_reached TRUE when no selected covariate reaches any outcome
#'   column at all.
#' @param rejuvenate The `rejuvenate` argument of [drbayes_select()].
#'
#' @return NULL, invisibly. Called for the warning.
#'
#' @keywords internal
#' @noRd
warn_unreached <- function(unreached, none_reached, rejuvenate) {
  if (length(unreached) == 0L) {
    return(invisible(NULL))
  }

  text <- paste0("The selected covariate(s) ", toString(unreached),
                 " contribute no column to the outcome model, so the ",
                 "g-computation does not adjust for them and posterior ",
                 "coupling has no outcome coefficient of theirs to revive. ",
                 "Add them to outcome.formula, or drop them from ps.formula ",
                 "if they are not confounders.")
  if (none_reached && rejuvenate == "selected") {
    text <- paste0(text, " No selected covariate reaches an outcome column, ",
                   "so with rejuvenate = \"selected\" the sweep can move ",
                   "nothing but the intercept and the treatment ",
                   "coefficient(s).")
  }
  warning(text)
  invisible(NULL)
}


#' Couple the two posteriors on the selected parameter block
#'
#' The counterpart of [run_tilting()] for Step 2 of Algorithm 4. It restricts
#' the moment condition to the selected propensity score coefficients, decides
#' which way along the tilting parameter grid to sweep, runs Algorithm 1 or
#' Algorithm 2 against the restricted condition, and reports whether it was
#' met.
#'
#' @param betas.ps,betas.otc Draws by parameters matrices of the full
#'   coefficient vectors.
#' @param ps_update,otc_update Column indices of the selected block, into
#'   `betas.ps` and `betas.otc` respectively.
#' @param d The object returned by [tilting_data()], whose `Z.ps` holds the
#'   selected columns only and whose `Z.lm` holds all of them.
#' @param control A [drbayes_control()] list.
#' @param method `"smc"` for Algorithm 2 or `"is"` for Algorithm 1.
#' @param rejuvenate `"all"` or `"selected"`, which coefficients the Liu and
#'   West kernel refreshes.
#' @param tol Tolerance on `|mean B_n^S|`, or NULL for [moment_tolerance()].
#'
#' @return A list with the tilted draws of both models and the diagnostics
#'   recorded in `fit$smc`, in the same shape [run_tilting()] records them
#'   plus `n_ancestors_cumulative`, the number of distinct Step 1 draws the
#'   returned particles descend from. That is `n_ancestors_min` accumulated
#'   over the whole sweep rather than taken at its worst single step, which is
#'   what the frozen block of `rejuvenate = "selected"` holds and what the per
#'   step count cannot see.
#'
#' @keywords internal
#' @noRd
select_tilting <- function(betas.ps, betas.otc, ps_update, otc_update, d,
                           control, method = "smc", rejuvenate = "all",
                           tol = NULL) {
  S    <- nrow(betas.ps)
  n.ps <- ncol(betas.ps)

  # The propensity score strata of equation (F.1) are built from the selected
  # covariates too, since that is the propensity score B_n^S is stated for.
  moment_fn <- moment_ipw
  if (identical(control$moment, "subclass")) {
    d$subclass_weight <- subclass_weights(betas.ps[, ps_update, drop = FALSE],
                                          d, control$n_subclass)
    moment_fn         <- moment_subclass
  }

  # Restricting the propensity score is all that B_n^S changes, so the sweep is
  # handed a moment condition that restricts its own argument and can then be
  # driven exactly as [run_tilting()] drives Algorithm 2.
  restricted <- function(betas.ps, betas.otc, d) {
    moment_fn(betas.ps[, ps_update, drop = FALSE], betas.otc, d)
  }

  BB <- restricted(betas.ps, betas.otc, d)
  if (is.null(tol)) tol <- moment_tolerance(BB)

  # Only the sign of the Newton step is used, to decide which way along the
  # lambda grid the constraint is approached.
  lam.new <- newton_lambda(BB, control$ridge)

  # Lemma 1: with a correctly specified outcome model the constraint already
  # holds at lambda = 0, so the sweep is skipped and the draws are unchanged.
  # Nothing is resampled, so all S particles remain in play and each is still
  # its own Step 1 draw.
  out <- list(betas.ps = betas.ps, betas.otc = betas.otc,
              lambda = 0, n_steps = 0L, BB = BB, ess = as.numeric(S),
              n_distinct = S, n_ancestors_cumulative = S)

  if (abs(mean(BB)) >= tol && lam.new != 0 && is.finite(lam.new)) {
    out <- if (method == "is") {
      tilt_is(betas.ps, betas.otc, BB, control, tol)
    } else if (rejuvenate == "all") {
      tilt_smc(betas.ps, betas.otc, BB, d, control, tol, sign(lam.new),
               restricted)
    } else {
      select_smc(betas.ps, betas.otc, c(ps_update, n.ps + otc_update), BB, d,
                 control, tol, sign(lam.new), restricted)
    }
    if (method == "is") {
      # Algorithm 1 resamples once, so the draws that survive that single step
      # are exactly the ones the returned particles descend from.
      out$n_ancestors_cumulative <- out$n_distinct
    }
  }

  # Judged at the draws that are handed back, as [run_tilting()] judges them.
  B_mean          <- mean(out$BB)
  ess             <- out$ess %||% NA_real_
  n_ancestors_min <- out$n_ancestors_min %||% NA_integer_
  n_distinct      <- out$n_distinct %||% NA_integer_
  n_cumulative    <- out$n_ancestors_cumulative %||% NA_integer_

  if (isTRUE(out$at_lambda_max)) {
    warning("Algorithm 1 reached the bound lambda_max = ",
            format(control$lambda_max, digits = 3), " without solving ",
            "equation (3.9) on the selected covariates, so the tilting ",
            "parameter returned is that bound and not a solution. Prefer ",
            "method = \"smc\", with drbayes_control(pruning = TRUE) if the ",
            "constraint is still not met.")
  }

  degeneracy <- tilting_degeneracy(ess, n_ancestors_min, S, control$ess_frac)
  if (degeneracy$degenerate) {
    warning(degeneracy$message)
  }
  if (abs(B_mean) >= tol) {
    warning("Posterior coupling did not satisfy the moment condition on the ",
            "selected covariates: |mean B_n| = ",
            format(abs(B_mean), digits = 3), " after ", out$n_steps,
            " steps up to lambda = ", format(out$lambda, digits = 3),
            ". The returned 'pc' draws are not doubly robust. Inspect $smc, ",
            "and consider a wider lambda range, a lower selection threshold, ",
            "or a better specified propensity score model.")
  }

  # The carried block is the one degeneracy the counts above cannot see: the
  # particles are all distinct, and a resampling step can look healthy, while
  # the coefficients no kernel refreshes have collapsed onto a handful of
  # values. Held to the same floor the tilting holds its own counts to.
  # Algorithm 1 freezes every coefficient alike, and its own effective sample
  # size has already answered for that above.
  frozen_collapse <- method == "smc" && rejuvenate == "selected" &&
    carried_collapse(n_cumulative, S, control$ess_frac)
  if (frozen_collapse) {
    warning("The coefficients outside the selected set are a ", n_cumulative,
            " point mass: resampling has taken all ", S, " particles back to ",
            "that many Step 1 draws, and rejuvenate = \"selected\" refreshes ",
            "only the selected block. That leaves ",
            ncol(betas.otc) - length(otc_update), " outcome coefficient(s) ",
            "and ", ncol(betas.ps) - length(ps_update), " propensity score ",
            "coefficient(s) carried verbatim, so their posterior spread, and ",
            "every credible interval that depends on them, is understated. ",
            "Use the default rejuvenate = \"all\", which refreshes the whole ",
            "particle, or take fewer steps with drbayes_control(n_steps = ).")
  }

  list(betas.ps  = out$betas.ps,
       betas.otc = out$betas.otc,
       smc = list(lambda    = out$lambda,
                  B_mean    = B_mean,
                  B_mean_weighted = out$B_mean_weighted %||% NA_real_,
                  tol       = tol,
                  n_steps   = out$n_steps,
                  max_steps = if (method == "is") control$newton_steps
                              else control$n_steps,
                  converged = abs(B_mean) < tol && !degeneracy$degenerate &&
                    !frozen_collapse,
                  ess       = ess,
                  n_ancestors_min = n_ancestors_min,
                  n_distinct      = n_distinct,
                  n_ancestors_cumulative = n_cumulative))
}


#' Whether the block the kernel never touched still describes a posterior
#'
#' Held to `ess_frac`, the floor [run_tilting()] holds its own particle counts
#' to, so that a carried block and a resampled one are judged the same way.
#'
#' @param n_cumulative Number of distinct Step 1 draws the returned particles
#'   descend from, or `NA` when nothing was carried.
#' @param S Number of particles the sweep started from.
#' @param ess_frac The floor, as a fraction of `S`.
#'
#' @return TRUE when the carried block has collapsed below the floor.
#'
#' @keywords internal
#' @noRd
carried_collapse <- function(n_cumulative, S, ess_frac) {
  !is.na(n_cumulative) && n_cumulative < ess_frac * S
}


#' Algorithm 2 with the kernel restricted to one block of the particle
#'
#' The sweep of [tilt_smc()], with two differences that only the literal
#' reading of Step 2 of Algorithm 4 needs. The Liu and West (2001) kernel is
#' applied to `update` alone, leaving the rest of each particle at the Step 1
#' draw it was resampled from, and the ancestry of the particles is followed so
#' that the caller can report how many distinct draws that frozen block still
#' holds. [tilt_smc()] does neither, which is why the loop is written out again
#' here rather than reached through [run_tilting()].
#'
#' Resampling still moves the whole particle. Dropping the frozen coefficients,
#' or resampling them separately, would leave each particle a mixture of
#' coefficients from different draws, and the g-computation that follows needs
#' a coherent full \eqn{\beta}.
#'
#' @param betas.ps,betas.otc Draws by parameters matrices of the full
#'   coefficient vectors.
#' @param update Column indices of the block the kernel refreshes, into the
#'   propensity score and outcome coefficients laid side by side.
#' @param BB The moment condition at the initial particles.
#' @param d The object returned by [tilting_data()].
#' @param control A [drbayes_control()] list.
#' @param tol Tolerance on `|mean B_n^S|`.
#' @param direction The sign of the Newton step, fixing which half of the
#'   lambda grid is swept.
#' @param moment_fn The moment condition to hold, already restricted to the
#'   selected covariates by the caller.
#'
#' @return As [tilt_smc()]: the tilted draws, `lambda`, `n_steps`, the final
#'   `BB` and `n_ancestors_min`, the smallest number of distinct particles to
#'   survive a resampling step. Plus `n_ancestors_cumulative`, the number of
#'   distinct initial particles the whole sweep has left the final cloud
#'   descending from, which is what the frozen block holds. That count is `NA`
#'   when `update` covers the whole particle, since then the kernel has
#'   refreshed every coefficient and nothing is frozen for it to describe.
#'
#' @keywords internal
#' @noRd
select_smc <- function(betas.ps, betas.otc, update, BB, d, control, tol,
                       direction, moment_fn = moment_ipw) {
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
  ancestor  <- seq_len(S)

  for (ll in seq_along(lam.diff)) {
    # Incremental importance weights exp{(lambda_t - lambda_{t-1}) B_n^S}
    ww <- normalised_weights(lam.diff[ll] * BB)

    # Pruning judges a particle by its weight against the untilted posterior at
    # the tilting parameter reached so far, not by the incremental weight, so
    # the sweep prunes here by exactly the rule [tilt_smc()] uses.
    if (isTRUE(control$pruning)) {
      cumulative <- normalised_weights(lam.s[ll + 1L] * BB)
      ww <- prune_weights(ww, cumulative, control)
    }

    # The kernel is centred on the cloud at t - 1, so its mean and covariance
    # are taken before the resampling replaces it.
    block <- particles[, update, drop = FALSE]
    mu    <- colMeans(block)
    Sigma <- stats::cov(block)

    resample_idx <- sample.int(S, size = S, replace = TRUE, prob = ww)
    n.eff.min    <- min(n.eff.min, length(unique(resample_idx)))
    ancestor     <- ancestor[resample_idx]
    particles    <- particles[resample_idx, , drop = FALSE]

    ep <- mvtnorm::rmvnorm(S, sigma = (1 - aa^2) * Sigma)
    particles[, update] <-
      sweep(aa * particles[, update, drop = FALSE] + ep, 2, (1 - aa) * mu, "+")

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

  # A selection that reached every coefficient leaves nothing frozen, and then
  # the ancestry is a property of the resampling rather than of the particles
  # handed back, every one of which the kernel has moved.
  carried <- length(update) < ncol(particles)

  list(betas.ps = betas.ps, betas.otc = betas.otc,
       lambda = lambda, n_steps = n.steps, BB = BB,
       n_ancestors_min = n.eff.min,
       n_ancestors_cumulative = if (carried) length(unique(ancestor))
                                else NA_integer_)
}


#' Resolve the family, the two links and their inverses
#'
#' [drbayes_pc()] applies the same rules inline rather than calling this. The
#' two agree on every family and link they accept and differ only in what they
#' say when they refuse one, so having [drbayes_pc()] call this would change
#' the text of errors its users may be matching on.
#'
#' @param family,link,ps.link As supplied by the user, possibly NULL.
#' @param Y The response, used to detect a binary outcome.
#'
#' @return A list with `family`, `link`, `ps.link`, `inverse_link` and
#'   `ps_inverse_link`.
#'
#' @keywords internal
#' @noRd
resolve_links <- function(family, link, ps.link, Y) {
  binary <- all(Y %in% c(0, 1))

  family <- if (is.null(family)) {
    if (binary) "binomial" else "gaussian"
  } else {
    match.arg(family, c("gaussian", "binomial"))
  }
  if (family == "binomial" && !binary) {
    stop("family = \"binomial\" needs a response that is 0 or 1. Use ",
         "family = \"gaussian\", or recode the outcome.")
  }

  valid <- switch(family,
                  gaussian = "identity",
                  binomial = c("logit", "probit"))
  if (is.null(link)) {
    link <- valid[1L]
  } else if (!is.character(link) || length(link) != 1L ||
             !(link %in% valid)) {
    stop("link = ", describe_value(link), " is not one of the links ",
         "available for family = \"", family, "\": ", toString(valid), ".")
  }

  # The propensity score has its own link. Without this the moment condition
  # would evaluate probit draws with plogis(), so the constraint solved would
  # not be the one the double robustness result is stated for.
  ps.link <- match.arg(ps.link, c("logit", "probit"))

  make_inverse <- function(name) {
    switch(name,
           identity = function(eta) eta,
           logit    = function(eta) stats::plogis(eta),
           probit   = function(eta) stats::pnorm(eta))
  }

  list(family = family, link = link, ps.link = ps.link,
       inverse_link = make_inverse(link),
       ps_inverse_link = make_inverse(ps.link))
}


#' Refuse a supplied sampler whose coefficients are on the wrong scale
#'
#' @param model The sampler the user supplied, or NULL.
#' @param expected The link its draws will be evaluated with.
#' @param arg Argument name, used in the error message.
#' @param chosen Name of the argument that set `expected`.
#'
#' @return NULL, invisibly.
#'
#' @keywords internal
#' @noRd
check_sampler_link <- function(model, expected, arg, chosen) {
  actual <- sampler_link(model)
  if (!is.null(actual) && actual != expected) {
    stop(arg, " draws coefficients on the ", actual, " scale, but ", chosen,
         " = \"", expected, "\" evaluates them with the ", expected,
         " inverse link. Set ", chosen, " = \"", actual, "\", or pass the ",
         expected, " sampler instead.")
  }
  invisible(NULL)
}


#' How Step 1 fitted one of the two models, for reporting back to the user
#'
#' @param model The sampler the user supplied, or NULL.
#' @param link The link of the model.
#' @param prior The prior requested, used only when `model` is NULL.
#'
#' @return A one line character description.
#'
#' @keywords internal
#' @noRd
step1_description <- function(model, link, prior) {
  if (is.null(model)) {
    paste0(sampler_name(link, prior), " (", link, " link, ", prior, " prior)")
  } else {
    paste0("the supplied sampler (", link, " link)")
  }
}


#' Whether a function takes an argument of a given name
#'
#' @param fun A function.
#' @param name Name of the argument.
#'
#' @return TRUE or FALSE.
#'
#' @keywords internal
#' @noRd
has_formal <- function(fun, name) {
  name %in% names(formals(args(fun)))
}


#' Assess the Step 1 draws before spending a sweep on them
#'
#' Posterior coupling tilts particles from both posteriors, so the doubly
#' robust estimate is only as trustworthy as the worse of the two. Vehtari et
#' al. (2021).
#'
#' [drbayes_pc()] runs the same two checks inline. This version adds what is
#' specific to Algorithm 4, that the selection itself is read off these draws,
#' so the two messages are not interchangeable.
#'
#' @param draws.otc,draws.ps Iterations by chains by parameters arrays.
#' @param control A validated [drbayes_control()] object.
#' @param diagnostics One of "warn", "error" or "none".
#'
#' @return The combined diagnostics table, or NULL when the check was skipped.
#'
#' @keywords internal
#' @noRd
check_step1_convergence <- function(draws.otc, draws.ps, control,
                                    diagnostics) {
  if (diagnostics == "none") {
    return(NULL)
  }

  tables <- list(outcome = convergence_diagnostics(draws.otc),
                 `propensity score` = convergence_diagnostics(draws.ps))

  messages   <- character(0)
  unreliable <- FALSE
  for (nm in names(tables)) {
    table  <- tables[[nm]]
    failed <- flag_convergence_failures(table,
                                        rhat_max = control$rhat_max,
                                        ess_min  = control$ess_min)
    if (failed$ok) next
    unreliable <- unreliable || isTRUE(failed$rhat_unreliable)

    parts <- character(0)
    if (length(failed$rhat) > 0L) {
      parts <- c(parts, paste0("max R-hat ",
                               format(max(table$rhat, na.rm = TRUE),
                                      digits = 4),
                               " (", toString(failed$rhat), ")"))
    }
    if (length(failed$ess) > 0L) {
      parts <- c(parts, paste0("min ESS ",
                               format(min(c(table$ess_bulk, table$ess_tail),
                                          na.rm = TRUE), digits = 4),
                               " (", toString(failed$ess), ")"))
    }
    messages <- c(messages, paste0(nm, " model: ",
                                   paste(parts, collapse = "; ")))
  }

  combined <- do.call(rbind, Map(function(table, nm) {
    table$model <- nm
    table[, c("model", setdiff(names(table), "model"))]
  }, tables, names(tables)))
  rownames(combined) <- NULL

  if (length(messages) > 0L) {
    advice <- if (unreliable) {
      paste("An effective sample size below", control$ess_min,
            "also makes R-hat itself unreliable, so a small R-hat here is not",
            "evidence of convergence. Raise mc, or thin less.")
    } else {
      "Raise mc, or use more chains from more dispersed starting values."
    }
    text <- paste0("The Step 1 posterior draws do not meet the convergence ",
                   "criteria (R-hat < ", control$rhat_max, ", ESS > ",
                   control$ess_min, "): ", paste(messages, collapse = " | "),
                   ". The selection is made on these draws, so it is not ",
                   "trustworthy either. ", advice, " Inspect $diagnostics.")
    if (diagnostics == "error") stop(text) else warning(text)
  }

  combined
}


#' Validate the selection threshold
#'
#' @param threshold The value supplied for `threshold`.
#'
#' @return The threshold, as a single number.
#'
#' @keywords internal
#' @noRd
validate_threshold <- function(threshold) {
  if (!is.numeric(threshold) || length(threshold) != 1L ||
      !is.finite(threshold) || threshold < 0) {
    stop("threshold must be a single non-negative number, the smallest ",
         "absolute posterior mean of a propensity score coefficient that is ",
         "selected. Received ", describe_value(threshold), ". Use 0 to ",
         "select every candidate covariate.")
  }
  as.numeric(threshold)
}
