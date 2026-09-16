#' Bayesian Doubly Robust Causal Inference via Posterior Coupling
#'
#' This function implements Bayesian doubly robust estimation for average treatment effects (ATE)
#' using posterior coupling. The method incorporates propensity score information via moment
#' conditions directly into the posterior distribution, avoiding the feedback problem and
#' enabling a fully Bayesian interpretation of DR estimation without requiring two-step estimation.
#'
#' The function supports both continuous outcomes (linear regression) and binary outcomes
#' (logistic and probit regression) through the family argument, using an intuitive formula
#' interface that integrates naturally with the R ecosystem.
#'
#' The function supports two modes of operation: internal posterior sampling using built-in
#' model functions, or external posterior samples from other Bayesian software packages
#' (e.g., Stan, JAGS, brms).
#'
#' @param outcome.formula Formula specifying the outcome model.
#'   For example: Y ~ A + X1 + X2 + A:X1. The response variable should be on the left side
#'   and predictors including treatment variable on the right side. Interactions and
#'   transformations are fully supported.
#' @param ps.formula Formula specifying the propensity score model.
#'   For example: A ~ X1 + X2 + X3. The treatment variable should be on the left side
#'   and confounders on the right side.
#' @param data Data.frame containing all variables specified in the formulas.
#' @param family Character string specifying the outcome model family. If NULL (default),
#'   automatically detects based on Y values (binary -> "binomial", otherwise -> "gaussian").
#'   Options are:
#'   \itemize{
#'     \item NULL (default): Automatically detect based on Y values
#'     \item "gaussian": Linear regression for continuous outcomes
#'     \item "binomial": Logistic or probit regression for binary outcomes
#'   }
#' @param link Character string specifying the link function of the OUTCOME model. If NULL
#'   (default), automatically determined from family: "identity" for gaussian, "logit" for
#'   binomial. Options are: "identity", "logit", "probit".
#' @param ps.link Character string specifying the link function of the PROPENSITY SCORE model,
#'   either "logit" (default) or "probit". This is the link the moment condition uses to turn
#'   the propensity score draws into probabilities, so it must match the scale on which
#'   \code{ps.model} or \code{ps.samples} produced them. Using \code{bayes_probit} or
#'   \code{bayes_probit_hs} without setting \code{ps.link = "probit"} evaluates probit coefficients
#'   with the logistic inverse link and silently solves a different constraint.
#' @param outcome.prior Character. Prior for the outcome model coefficients, either
#'   \code{"normal"} (default) or \code{"horseshoe"}. Given the family and the link, the
#'   sampling algorithm follows, so the prior is the only choice left to make and there is
#'   no need to hand in a sampler. Use \code{"horseshoe"} when there are many candidate
#'   confounders; the treatment effect and its interactions are exempted from shrinkage
#'   automatically.
#' @param ps.prior Character. As \code{outcome.prior}, for the propensity score model.
#'   The default \code{"normal"} is the right choice unless the model is high dimensional.
#' @param chains Integer. Number of MCMC chains, started from dispersed values. Default 4,
#'   the minimum Vehtari et al. (2021) recommend for the convergence diagnostics.
#' @param method Character. How the tilting parameter is found. \code{"smc"} is
#'   Algorithm 2 of the paper, which walks a grid of tilting parameters and
#'   rejuvenates the particles at each step. \code{"is"} is Algorithm 1, a single
#'   importance sampling step, which is faster but degenerates when the two
#'   posteriors put little mass where the moment condition holds. Default
#'   \code{"smc"}.
#' @param control List of tuning parameters from \code{\link{drbayes_control}}, covering
#'   the sequential Monte Carlo sweep and the diagnostic thresholds. \code{mc}, \code{bn},
#'   \code{thin} and \code{chains} also live there; naming one of those directly overrides
#'   the control object, which reads better for the settings people change most often.
#' @param seed Integer or NULL. When given, the fit is reproducible and the caller's random
#'   number stream is restored afterwards.
#' @param verbose Logical. Whether to report progress through
#'   \code{\link[base]{message}}. Default TRUE.
#' @param outcome.model Function for fitting the outcome model, for callers who want a
#'   sampler the package does not provide. Overrides \code{outcome.prior}. Should take arguments Y, X, and mc.
#'   Available functions include \code{\link{bayes_lm}} for continuous outcomes, \code{\link{bayes_logit}}
#'   and \code{\link{bayes_probit}} for binary outcomes, and \code{\link{bayes_lm_hs}}, \code{\link{bayes_logit_hs}},
#'   and \code{\link{bayes_probit_hs}} for horseshoe priors. Set to NULL when using external samples. Default is NULL.
#' @param ps.model Function for fitting the propensity score model. Should take arguments Y, X, and mc.
#'   Available functions include \code{\link{bayes_logit}}, \code{\link{bayes_probit}}, \code{\link{bayes_logit_hs}},
#'   and \code{\link{bayes_probit_hs}} for binary treatments. Set to NULL when using external samples. Default is NULL.
#' @param outcome.samples Matrix of posterior samples for outcome model coefficients. Each row
#'   represents one MCMC iteration and each column represents a coefficient. The first column
#'   should be the intercept, followed by coefficients for covariates in the outcome model.
#'   Set to NULL when using internal model functions. Default is NULL.
#' @param ps.samples Matrix of posterior samples for propensity score model coefficients. Each row
#'   represents one MCMC iteration and each column represents a coefficient. The first column
#'   should be the intercept, followed by coefficients for covariates in the propensity score model.
#'   Set to NULL when using internal model functions. Default is NULL.
#' @param outcome.mu Matrix of posterior draws of the fitted outcome mean at the
#'   observed treatment and covariates, with one row per draw and one column per
#'   observation used. Supplying it, together with \code{outcome.mu1} and
#'   \code{outcome.mu0}, is the nonparametric outcome path of Section 5.4, in
#'   which the outcome model is whatever produced the draws and no coefficients
#'   are needed. It is mutually exclusive with \code{outcome.model},
#'   \code{outcome.samples}, \code{outcome.prior}, \code{outcome.priors} and
#'   \code{link}, since none of them has anything left to describe, and the
#'   propensity score side is supplied as usual. Default is NULL.
#' @param outcome.mu1,outcome.mu0 The same fitted means with the treatment set
#'   to 1 and to 0 for every observation, from which the g-computation (3.7) is
#'   formed. Both have the same dimensions as \code{outcome.mu}, and all three
#'   are supplied together. Default is NULL.
#' @param mc Integer. Number of MCMC iterations for internal sampling. Ignored when using
#'   external samples. Default is 5000.
#' @param bn Integer. Number of burn-in iterations to discard. Applied to both internal
#'   and external samples. Default is 1000.
#' @param thin Integer. Thinning interval for MCMC samples. Applied to both internal
#'   and external samples. Default is 2.
#' @param outcome.priors List of prior distribution parameters to pass to the
#'   outcome model function. For example, when using \code{bayes_lm}, you can specify
#'   \code{list(theta_prior = 1/50, sigma_prior = c(2, 1))} to customize
#'   the prior distributions for regression coefficients and error variance.
#'   Ignored when using external samples. Default is an empty list (uses function defaults).
#' @param ps.priors List of prior distribution parameters to pass to the
#'   propensity score model function. Ignored when using external samples.
#'   Default is an empty list (uses function defaults).
#' @param diagnostics Character. What to do when the posterior draws fail the convergence
#'   criteria: \code{"warn"} (default), \code{"error"}, or \code{"none"} to skip the check.
#'   Both posteriors are assessed before the sequential Monte Carlo sweep, because posterior
#'   coupling is only as trustworthy as the worse of the two.
#' @param na.action Character string specifying how to handle missing values.
#'   Options are "na.omit" (default), "na.fail", or "na.exclude".
#'
#' @return An object of class \code{"DRBayes"}, a list containing:
#' \describe{
#'   \item{g.comp}{Numeric vector of posterior samples for the ATE using g-computation,
#'     that is, the untilted posterior. See \emph{Which average} below for what
#'     the average is taken over.}
#'   \item{pc}{Numeric vector of posterior samples for the ATE using posterior coupling.
#'     This is the doubly robust estimand and the quantity to report.}
#'   \item{family}{Character string indicating the outcome model family used.}
#'   \item{link}{Character string indicating the outcome model link function used.}
#'   \item{ps.link}{Character string indicating the propensity score link function used.}
#'   \item{smc}{List of sequential Monte Carlo diagnostics: the tilting parameter
#'     \code{lambda}, the achieved posterior mean \code{B_mean} of the moment condition,
#'     the tolerance \code{tol} it was compared against, the number of steps taken
#'     \code{n_steps} out of \code{max_steps}, and \code{converged}. When
#'     \code{converged} is FALSE the moment condition was never satisfied and the
#'     \code{pc} draws carry no double robustness guarantee.}
#'   \item{diagnostics}{Data frame of convergence diagnostics, one row per parameter of each
#'     model: rank-normalised split R-hat, bulk and tail effective sample size, and the Monte
#'     Carlo standard error of the mean. NULL when \code{diagnostics = "none"}.}
#'   \item{outcome.prior, ps.prior}{The priors used.}
#'   \item{method}{Which tilting algorithm ran, \code{"smc"} or \code{"is"}.
#'     Worth reading rather than assuming: the nonparametric outcome path can
#'     only run Algorithm 1, so it resolves to \code{"is"} whatever the
#'     default says.}
#'   \item{control}{The resolved \code{\link{drbayes_control}} object.}
#'   \item{particles}{When \code{drbayes_control(keep_particles = TRUE)} asked for
#'     them, a list holding the tilted draws \code{outcome} and \code{ps}, the
#'     untilted draws \code{outcome_untilted} and \code{ps_untilted}, the design
#'     matrices and links in \code{model_data}, and the counterfactual design
#'     matrices \code{Z.lm1} and \code{Z.lm0}. This is what
#'     \code{\link{drbayes_sensitivity}} reweights. On the nonparametric path
#'     there are no outcome coefficients, so it holds \code{ps},
#'     \code{ps_untilted}, \code{model_data} and the resampled draw
#'     \code{index}, and the fit cannot be reweighted. NULL otherwise.}
#'   \item{call}{The matched call.}
#'   \item{data_info}{Summary of the processed data: \code{n_observations},
#'     \code{n_treated}, \code{n_control}, the number of non-intercept columns
#'     \code{n_outcome_terms} and \code{n_ps_terms} in each design matrix, both formulas,
#'     and \code{missing_observations}.}
#' }
#'
#' @details
#' \strong{Formula Interface:}
#'
#' drbayes_pc uses R's standard formula notation, making it intuitive and consistent with
#' other statistical functions like lm() and glm(). The formula interface automatically
#' handles missing values, creates appropriate design matrices, and manages factor variables.
#' \preformatted{
#' # Basic usage
#' result <- drbayes_pc(
#'   outcome.formula = Y ~ A + X1 + X2,
#'   ps.formula = A ~ X1 + X2 + X3,
#'   data = my_data,
#'   outcome.model = bayes_lm,
#'   ps.model = bayes_logit
#' )
#'
#' # With interactions
#' result <- drbayes_pc(
#'   outcome.formula = Y ~ A + X1 + X2 + A:X1,
#'   ps.formula = A ~ X1 + X2 + X3 + I(X1^2),
#'   data = my_data,
#'   outcome.model = bayes_lm,
#'   ps.model = bayes_logit
#' )
#' }
#'
#' \strong{Model Families and Link Functions:}
#'
#' The function supports different outcome model types:
#' \itemize{
#'   \item \strong{Automatic detection} (family=NULL): Determines family based on Y values.
#'   If all Y values are 0 or 1, uses "binomial"; otherwise uses "gaussian".
#'   \item \strong{Gaussian family} (family="gaussian"): For continuous outcomes using linear regression
#'   with identity link. ATE is estimated as the mean difference
#'   \eqn{E[Y | A = 1, X] - E[Y | A = 0, X]}.
#'   \item \strong{Binomial family} (family="binomial"): For binary outcomes using logistic regression
#'   (link="logit") or probit regression (link="probit"). ATE is estimated as the probability
#'   difference P(Y=1|A=1,X) - P(Y=1|A=0,X).
#' }
#'
#' \strong{Sampling Methods:}
#'
#' The function automatically detects the sampling mode based on the provided arguments:
#' \itemize{
#'   \item \strong{Internal sampling}: Provide \code{outcome.model} and \code{ps.model}
#'   \item \strong{External samples}: Provide \code{outcome.samples} and \code{ps.samples}
#' }
#'
#' If both internal and external options are provided, external samples take precedence
#' with a warning message.
#'
#' \strong{Internal Sampling Mode:}
#' The function generates posterior samples from separate outcome and propensity score models
#' using the provided model functions and performs the following steps:
#' \enumerate{
#'   \item Generates posterior samples from separate outcome and propensity score models
#'   \item Applies sequential Monte Carlo (SMC) to enforce the moment condition constraint
#'   \item Computes posterior distributions for the average treatment effect
#' }
#'
#' \strong{External Samples Mode:}
#' When using external posterior samples (e.g., from Stan, JAGS, brms), ensure that:
#' \itemize{
#'   \item Samples are in matrix format with rows = iterations, columns = coefficients
#'   \item The first column contains intercept terms
#'   \item Coefficient ordering matches the column order in the design matrices created from formulas
#'   \item Both outcome.samples and ps.samples have the same number of iterations
#' }
#'
#' The moment condition enforced is:
#' \deqn{E[(A - \pi(X)) \cdot (Y - \mu(X)) / (\pi(X)(1-\pi(X)))] = 0}
#'
#' where \eqn{\pi(X)} is the propensity score and \eqn{\mu(X)} is the outcome model.
#'
#' @section Which average:
#'
#' `g.comp` and `pc` are both posteriors for equation (3.7), the average of the
#' fitted contrast over the covariate vectors that were observed. Section 2.1
#' of the paper defines the estimand as \eqn{\tau = E[Y_1 - Y_0]}, an average
#' over the population those vectors were drawn from. The two are not the same
#' quantity, and the difference is not in their centre but in their spread:
#' averaging over the observed covariates treats their distribution as known,
#' so these intervals do not carry the uncertainty of not knowing it.
#'
#' They coincide when the fitted contrast \eqn{c_i = m_1(X_i) - m_0(X_i)} is the
#' same for every unit, and differ when it is not. Two things make it vary: a
#' treatment by covariate interaction in the outcome formula, and a link that
#' is not the identity. Only the first matters in practice. Writing \eqn{V} for
#' the posterior variance of the average, putting Dirichlet weights on the
#' observations, as a Bayesian bootstrap would, adds \eqn{\mathrm{Var}_i(c_i)/n}
#' to it, so the interval widens by a factor of
#' \eqn{\sqrt{1 + \mathrm{Var}_i(c_i)/(nV)}}. Measured on fits:
#'
#' \tabular{lrr}{
#'   \strong{outcome model} \tab \strong{n} \tab \strong{factor} \cr
#'   `Y ~ A + X1 + X2 + A:X1` \tab 500 \tab 1.11 \cr
#'   `Y ~ A + X1 + X2 + A:X1` \tab 2000 \tab 1.21 \cr
#'   `Y ~ A + X1 + X2`, identity link \tab 500 \tab 1 exactly \cr
#'   `Y ~ A + X1 + X2`, logit link \tab 500 \tab 1.001 \cr
#'   `Y ~ A + X1 + X2`, logit link \tab 2000 \tab 1.001
#' }
#'
#' A logit or probit link does make the contrast vary, but by far too little to
#' matter against the posterior uncertainty; a treatment interaction is the case
#' to watch, and then how large the factor is depends on how much the effect
#' actually varies. The 1.11 and 1.21 above are for a unit level effect of
#' \eqn{2 + 1.5 X_1}; a milder interaction gives a milder factor. For a fit of
#' your own, `drbayes_control(keep_particles = TRUE)` returns what is needed to
#' compute it directly, and the getting started vignette shows how.
#'
#' None of this affects the point estimate or Theorem 1. Dirichlet weights have
#' mean \eqn{1/n}, so the posterior mean is unchanged, and the theorem concerns
#' the posterior mean of an estimator whose two averages differ by
#' \eqn{O_p(n^{-1/2})}. What changes is the credible interval, and which
#' estimand it is a credible interval for.
#'
#' @section Link Functions:
#' \itemize{
#'   \item \strong{Identity} (family="gaussian"): \eqn{\mu = X\beta}
#'   \item \strong{Logit} (family="binomial", link="logit"): \eqn{\mu = \mathrm{logit}^{-1}(X\beta) = \frac{e^{X\beta}}{1 + e^{X\beta}}}
#'   \item \strong{Probit} (family="binomial", link="probit"): \eqn{\mu = \Phi(X\beta)} where \eqn{\Phi} is the standard normal CDF
#' }
#'
#' @section Nonparametric Outcome Models:
#' Section 5.4 of the paper puts BART on the outcome side instead of a linear
#' predictor. Any outcome model can be coupled that way, as long as it reports a
#' posterior draw of its fitted mean: pass those draws as \code{outcome.mu},
#' \code{outcome.mu1} and \code{outcome.mu0} in place of a formula-based
#' sampler. Fitted means from packages such as \code{dbarts} or \code{BART} are
#' the usual source, but nothing here depends on which software produced them.
#' The means are used exactly as supplied, on the scale of the outcome itself,
#' so no link is applied to them.
#'
#' \preformatted{
#' # Draws of the fitted mean from any Bayesian outcome model, S by n
#' fit <- drbayes_pc(
#'   outcome.formula = Y ~ A + X1 + X2,
#'   ps.formula = A ~ X1 + X2,
#'   data = data,
#'   outcome.mu = mu, outcome.mu1 = mu1, outcome.mu0 = mu0,
#'   ps.model = bayes_logit,
#'   chains = 1L, mc = 2 * nrow(mu) + 1000, bn = 1000, thin = 2
#' )
#' }
#'
#' The outcome formula is still needed, since it names the response and the
#' treatment and fixes which rows are complete; its right hand side plays no
#' further part. The moment condition (3.4) reads \eqn{m_{A_i}(X_i)} from
#' \code{outcome.mu}, and the average treatment effect (3.7) is the row mean of
#' \code{outcome.mu1} minus the row mean of \code{outcome.mu0}. All three
#' matrices need one row per propensity score draw that survives \code{bn},
#' \code{thin} and \code{chains} -- all four chains contribute by default, so
#' \code{chains = 1L} is usually what a supplied set of outcome draws wants --
#' and one column per observation used, in the order the rows of \code{data}
#' are in.
#'
#' The draws are resampled systematically, with replacement, which is the
#' standard choice and leaves the resampled cloud an unbiased representation of
#' the tilted posterior. The footnote to Table 2 of the paper reports something
#' different for its own BART results, "importance sampling without
#' replacement", keeping 10,000 of 20,000 draws. Selecting without replacement
#' guarantees that many distinct draws survive, which matters here because
#' nothing replenishes the cloud, but it does so at the cost of inclusion
#' probabilities that are no longer proportional to the weights. The two are
#' not interchangeable, and this package implements the former.
#'
#' Only \code{method = "is"} is available. Algorithm 2 rejuvenates its particles
#' with a Gaussian kernel on the outcome model's parameter vector, and a draw of
#' fitted values has none, so the tilting reweights and resamples the draw index
#' instead, carrying the propensity score draw and the three rows of fitted
#' means together. Nothing then replenishes the cloud, which makes \code{ess} in
#' \code{$smc} the number to watch: a fit resting on fewer than
#' \code{ess_frac} of the draws is reported as not converged.
#'
#' \code{\link{drbayes_sensitivity}} does not apply to a fit made this way,
#' because Algorithm 3 reweights the outcome model's coefficient draws.
#'
#' @section External Software Integration:
#' This function is designed to work with posterior samples from various Bayesian software:
#'
#' \strong{Stan/cmdstanr:}
#' \preformatted{
#' # Extract samples
#' outcome.samples <- fit$draws("beta", format = "matrix")
#' ps.samples <- fit$draws("gamma", format = "matrix")
#'
#' # Use in drbayes_pc
#' result <- drbayes_pc(
#'   outcome.formula = Y ~ A + X1 + X2 + A:X1,
#'   ps.formula = A ~ X1 + X2 + X3,
#'   data = data,
#'   outcome.samples = outcome.samples,
#'   ps.samples = ps.samples
#' )
#' }
#'
#' \strong{brms:}
#' \preformatted{
#' # Extract fixed effects
#' outcome.samples <- as.matrix(brms_fit)[, 1:ncol(model.matrix(outcome.formula, data))]
#' ps.samples <- as.matrix(brms_fit)[, 1:ncol(model.matrix(ps.formula, data))]
#'
#' # Use in drbayes_pc
#' result <- drbayes_pc(
#'   outcome.formula = Y ~ A + X1 + X2,
#'   ps.formula = A ~ X1 + X2 + X3,
#'   data = data,
#'   outcome.samples = outcome.samples,
#'   ps.samples = ps.samples
#' )
#' }
#'
#' \strong{JAGS/R2jags:}
#' \preformatted{
#' # Extract parameter samples
#' outcome.samples <- jags.fit$BUGSoutput$sims.matrix[, grep("beta", colnames(...))]
#' ps.samples <- jags.fit$BUGSoutput$sims.matrix[, grep("gamma", colnames(...))]
#'
#' # Use in drbayes_pc
#' result <- drbayes_pc(
#'   outcome.formula = Y ~ A + X1 + X2,
#'   ps.formula = A ~ X1 + X2 + X3,
#'   data = data,
#'   outcome.samples = outcome.samples,
#'   ps.samples = ps.samples
#' )
#' }
#'
#' @references
#' Orihara, S., Momozaki, T., & Sugasawa, S. (2025). Bayesian Doubly Robust Causal Inference via
#' Posterior Coupling. arXiv preprint arXiv:2506.04868.
#'
#' @examples
#' set.seed(123)
#' n <- 200
#' dat <- data.frame(X1 = rnorm(n), X2 = rnorm(n), X3 = rnorm(n))
#' dat$A <- rbinom(n, 1, plogis(0.2 + 0.5 * dat$X1 - 0.3 * dat$X2))
#' dat$Y <- 1 + 1.5 * dat$A + 0.8 * dat$X1 + 0.6 * dat$X2 +
#'   0.4 * dat$A * dat$X1 + rnorm(n)
#'
#' # The treatment effect varies with X1, so the interaction belongs in the
#' # outcome formula: the g-formula averages the unit level contrast over the
#' # covariates, and it can only average what the model can express.
#' fit <- drbayes_pc(Y ~ A + X1 + X2 + A:X1, A ~ X1 + X2 + X3, dat,
#'                   outcome.model = bayes_lm, ps.model = bayes_logit,
#'                   control = drbayes_control(mc = 1500, bn = 500))
#' fit
#' summary(fit)
#'
#' # Whether the tilting did anything is in `$smc`. `lambda` is zero only when
#' # the untilted posterior already satisfied the moment condition, and
#' # `converged` says whether it was satisfied at all; when it was not, the
#' # draws in `pc` carry no double robustness guarantee.
#' fit$smc[c("lambda", "B_mean", "tol", "converged", "n_steps")]
#'
#' # A binary outcome. The estimand is then a difference of probabilities.
#' dat$Ybin <- rbinom(n, 1, plogis(0.5 + 0.8 * dat$A + 0.4 * dat$X1))
#' drbayes_pc(Ybin ~ A + X1 + X2, A ~ X1 + X2 + X3, dat,
#'            outcome.model = bayes_logit, ps.model = bayes_logit,
#'            control = drbayes_control(mc = 1500, bn = 500))
#'
#' # Transformations are rebuilt through the fitted terms, so the poly() basis
#' # in the two counterfactual designs is the one the model was fitted with
#' # rather than one re-derived from data that has changed.
#' drbayes_pc(Y ~ A + log(abs(X1) + 1) + poly(X3, 2) + A:X1,
#'            A ~ X1 + X2 + X3, dat,
#'            outcome.model = bayes_lm, ps.model = bayes_logit,
#'            control = drbayes_control(mc = 1000, bn = 300))
#'
#' # Draws from elsewhere, supplied instead of fitting. Anything that returns a
#' # draws-by-parameters matrix will do -- Stan, JAGS, brms -- provided its
#' # columns are those of model.matrix() for the same formula, and the
#' # propensity score draws are on the scale `ps.link` names. Here they come
#' # from this package's own samplers, which take the covariates without an
#' # intercept and add one.
#' Zo <- model.matrix(Y ~ A + X1 + X2, dat)[, -1, drop = FALSE]
#' Zp <- model.matrix(A ~ X1 + X2 + X3, dat)[, -1, drop = FALSE]
#' otc <- bayes_lm(dat$Y, Zo, mc = 3000, chains = 2L)
#' ps  <- bayes_logit(dat$A, Zp, mc = 3000, chains = 2L)
#' pool <- function(x) {
#'   matrix(x, ncol = dim(x)[3], dimnames = list(NULL, dimnames(x)[[3]]))
#' }
#' drbayes_pc(Y ~ A + X1 + X2, A ~ X1 + X2 + X3, dat,
#'            outcome.samples = pool(otc), ps.samples = pool(ps),
#'            control = drbayes_control(bn = 1000))
#'
#' @seealso
#' \code{\link{bayes_lm}} for Bayesian linear regression,
#' \code{\link{bayes_logit}} for Bayesian logistic regression,
#' \code{\link{bayes_probit}} for Bayesian probit regression,
#' \code{\link{bayes_lm_hs}} for horseshoe linear regression,
#' \code{\link{bayes_logit_hs}} for horseshoe logistic regression,
#' \code{\link{bayes_probit_hs}} for horseshoe probit regression,
#' \code{\link{drbayes_bb}} for Bayesian bootstrap approach.
#'
#' @importFrom stats cov plogis pnorm model.matrix model.frame terms na.omit na.fail na.exclude
#' @importFrom mvtnorm rmvnorm
#'
#' @export
drbayes_pc <- function(outcome.formula, ps.formula, data,
                       family = NULL, link = NULL, ps.link = "logit",
                       outcome.prior = c("normal", "horseshoe"),
                       ps.prior      = c("normal", "horseshoe"),
                       mc = 5000, bn = 1000, thin = 2, chains = 4L,
                       method = c("smc", "is"),
                       control = drbayes_control(),
                       diagnostics = c("warn", "error", "none"),
                       outcome.model   = NULL, ps.model   = NULL,
                       outcome.samples = NULL, ps.samples = NULL,
                       outcome.mu  = NULL,
                       outcome.mu1 = NULL, outcome.mu0 = NULL,
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

  # The four MCMC settings are also on the control object; naming one directly
  # is the common case and reads better than wrapping it.
  control <- resolve_control(control)
  supplied <- names(as.list(match.call())[-1L])
  for (nm in intersect(c("mc", "bn", "thin", "chains"), supplied)) {
    control[[nm]] <- get(nm)
  }
  control    <- resolve_control(control)
  mc         <- control$mc
  bn         <- control$bn
  thin       <- control$thin
  chains     <- control$chains
  tol        <- control$tol
  rhat_max   <- control$rhat_max
  ess_min    <- control$ess_min

  # Store the call for reference
  call_info <- match.call()

  # ============================================================================
  # 1. Input validation and formula processing
  # ============================================================================

  if (!is.numeric(mc) || length(mc) != 1L || mc != round(mc) || mc < 1L) {
    stop("mc must be a single positive integer")
  }
  if (!is.numeric(bn) || length(bn) != 1L || bn != round(bn) || bn < 0L) {
    stop("bn must be a single non-negative integer")
  }
  if (!is.numeric(thin) || length(thin) != 1L || thin != round(thin) || thin < 1L) {
    stop("thin must be a single positive integer")
  }
  if (mc <= bn + 1L) {
    stop("mc must exceed bn + 1, otherwise no posterior draw survives the ",
         "burn-in (mc = ", mc, ", bn = ", bn, ")")
  }

  md <- prepare_model_data(outcome.formula, ps.formula, data, na.action)

  Y                  <- md$Y
  A                  <- md$A
  X.lm               <- md$X.lm
  X.ps               <- md$X.ps
  Z.lm               <- md$Z.lm
  Z.ps               <- md$Z.ps
  Z.lm1              <- md$Z.lm1
  Z.lm0              <- md$Z.lm0
  treatment_var_name <- md$treatment
  data_info          <- md$data_info

  nn    <- nrow(Z.lm)
  pp.lm <- ncol(X.lm)
  pp.ps <- ncol(X.ps)

  say("Formula processing complete:")
  say("  Observations used: ", nn, " out of ", nrow(data))
  say("  Treatment group: ", sum(A), " | Control group: ", sum(1 - A))
  say("  Outcome model variables: ", pp.lm)
  say("  Propensity score variables: ", pp.ps)

  # ============================================================================
  # 2. Family and link function determination
  # ============================================================================

  # Determine family
  if (is.null(family)) {
    if (all(Y %in% c(0, 1))) {
      family <- "binomial"
      say("Auto-detected binary outcome: using family='binomial'")
    } else {
      family <- "gaussian"
      say("Auto-detected continuous outcome: using family='gaussian'")
    }
  } else {
    # Validate family argument
    family <- match.arg(family, choices = c("gaussian", "binomial"))
  }

  # Determine link function
  if (is.null(link)) {
    link <- switch(family,
                   "gaussian" = "identity",
                   "binomial" = "logit")
  } else {
    # Validate link argument
    valid_links <- switch(family,
                          "gaussian" = "identity",
                          "binomial" = c("logit", "probit"))
    if (!link %in% valid_links) {
      stop(sprintf("Link '%s' not supported for family '%s'. Valid links: %s",
                   link, family, paste(valid_links, collapse = ", ")))
    }
  }

  # Validate outcome values based on family
  if (family == "binomial") {
    if (!all(Y %in% c(0, 1))) {
      stop("For family='binomial', Y must be binary (0 or 1)")
    }
  }

  # Define inverse link function
  inverse_link <- function(eta) {
    switch(link,
           "identity" = eta,
           "logit"    = plogis(eta),
           "probit"   = pnorm(eta))
  }

  # The propensity score has its own link. Without this the moment condition
  # (3.4) would evaluate probit draws with plogis(), so the constraint solved
  # would not be the one the double robustness result is stated for.
  ps.link <- match.arg(ps.link, choices = c("logit", "probit"))
  ps_inverse_link <- function(eta) {
    switch(ps.link,
           "logit"  = plogis(eta),
           "probit" = pnorm(eta))
  }

  # Guard against a sampler whose coefficients are on a different scale from
  # the link they will be evaluated with.
  check_model_link <- function(model, expected, arg, chosen) {
    if (is.null(model)) return(invisible(NULL))
    known <- list(identity = list(bayes_lm, bayes_lm_hs),
                  logit    = list(bayes_logit, bayes_logit_hs),
                  probit   = list(bayes_probit, bayes_probit_hs))
    for (lk in names(known)) {
      if (any(vapply(known[[lk]], identical, logical(1), model)) && lk != expected) {
        stop(arg, " draws coefficients on the ", lk, " scale, but ", chosen,
             " = \"", expected, "\" evaluates them with the ", expected,
             " inverse link.")
      }
    }
    invisible(NULL)
  }
  check_model_link(outcome.model, link, "outcome.model", "link")
  check_model_link(ps.model, ps.link, "ps.model", "ps.link")

  # Given the family and the link, the sampling algorithm is determined, so the
  # user does not have to hand one in. What is genuinely their choice is the
  # prior, and that is what outcome.prior and ps.prior select. An explicitly
  # supplied sampler always wins, as do external draws.
  outcome.prior <- match.arg(outcome.prior)
  ps.prior      <- match.arg(ps.prior)
  method        <- match.arg(method)

  # The nonparametric outcome path of Section 5.4: the outcome model arrives as
  # draws of its fitted means, so there is no linear predictor to evaluate and
  # nothing on the outcome side left for the package to choose or to fit.
  use_nonparametric <- !is.null(outcome.mu) || !is.null(outcome.mu1) ||
    !is.null(outcome.mu0)

  if (use_nonparametric) {
    check_fitted_means(outcome.mu, outcome.mu1, outcome.mu0, nn, nrow(data))

    clash <- intersect(c("outcome.model", "outcome.samples", "outcome.prior",
                         "outcome.priors", "link"), supplied)
    if (length(clash) > 0L) {
      stop("outcome.mu supplies the fitted outcome means themselves, so ",
           toString(clash), if (length(clash) == 1L) " has" else " have",
           " nothing left to describe. Drop ",
           if (length(clash) == 1L) "it" else "them",
           ", or drop outcome.mu and fit the outcome model in the usual way.")
    }

    # The means are used on the scale they arrive on, which is the scale of the
    # outcome, so that is what the fit records having done to them.
    link         <- "identity"
    inverse_link <- function(eta) eta

    if (method != "is" && "method" %in% supplied) {
      stop("method = \"smc\" cannot be used with outcome.mu, because ",
           "Algorithm 2 rejuvenates its particles with a Gaussian kernel on ",
           "the outcome model's parameter vector and a draw of fitted means ",
           "has no parameter vector to smooth. Use method = \"is\", which ",
           "reweights and resamples the draws themselves.")
    }
    method <- "is"
  }

  chosen_outcome <- NULL
  chosen_ps      <- NULL
  if (is.null(outcome.samples) && is.null(ps.samples)) {
    if (is.null(outcome.model) && !use_nonparametric) {
      outcome.model  <- select_sampler(link, outcome.prior)
      chosen_outcome <- sampler_name(link, outcome.prior)
    }
    if (is.null(ps.model)) {
      ps.model  <- select_sampler(ps.link, ps.prior)
      chosen_ps <- sampler_name(ps.link, ps.prior)
    }
  }

  # ============================================================================
  # 3. Sampling method determination
  # ============================================================================

  # Determine sampling method
  use_external <- !is.null(outcome.samples) && !is.null(ps.samples)
  use_internal <- !is.null(outcome.model) && !is.null(ps.model)

  # Check argument consistency
  if (use_external && use_internal) {
    warning("Both external samples and models provided. Using external samples.")
    use_internal <- FALSE
  } else if (use_external && (is.null(outcome.samples) || is.null(ps.samples))) {
    stop("Both outcome.samples and ps.samples must be provided for external sampling")
  } else if (use_internal && (is.null(outcome.model) || is.null(ps.model))) {
    stop("Both outcome.model and ps.model must be provided for internal sampling")
  }

  # ============================================================================
  # 5. Obtain posterior samples
  # ============================================================================

  # Which iterations of externally supplied draws survive the burn-in and the
  # thinning. The nonparametric path shares the rule so that its propensity
  # score draws are treated exactly as they would be alongside an outcome model.
  external_keep <- function(n_iterations) {
    if (n_iterations <= (bn + 100)) {
      warning("External samples may be too few for burn-in. Using all samples.")
      seq_len(n_iterations)
    } else {
      seq(bn + 1L, n_iterations, by = thin)
    }
  }

  if (use_nonparametric) {
    say("Using the supplied draws of the fitted outcome mean...")
    draws.otc <- NULL

    if (!is.null(ps.samples)) {
      if (!is.null(ps.model)) {
        warning("Both external samples and models provided. Using external ",
                "samples.")
      }
      ps.samples <- as_draws_array(ps.samples, "ps.samples")
      if (dim(ps.samples)[3] != ncol(Z.ps)) {
        stop(sprintf(paste("ps.samples must have %d parameters (intercept +",
                           "%d covariates), but has %d"),
                     ncol(Z.ps), pp.ps, dim(ps.samples)[3]))
      }
      keep     <- external_keep(dim(ps.samples)[1])
      draws.ps <- ps.samples[keep, , , drop = FALSE]

    } else {
      say("Generating propensity score samples using ",
          chosen_ps %||% "the supplied propensity score model",
          " (", ps.link, " link, ", ps.prior, " prior).")
      ps.args  <- c(list(Y = A, X = X.ps, mc = mc, chains = chains), ps.priors)
      post.ps  <- do.call(ps.model, ps.args)
      draws.ps <- post.ps[seq(bn + 1L, mc, by = thin), , , drop = FALSE]
      rm(post.ps)   # released as soon as the thinned draws are taken from it
    }

  } else if (use_external) {
    # Use externally provided samples
    say("Using externally provided posterior samples...")

    # Accept the three-dimensional iterations x chains x parameters layout as
    # well as a plain matrix, so draws from Stan, JAGS or brms keep their chain
    # structure and can be diagnosed rather than being silently pooled.
    outcome.samples <- as_draws_array(outcome.samples, "outcome.samples")
    ps.samples      <- as_draws_array(ps.samples, "ps.samples")

    if (dim(outcome.samples)[1] != dim(ps.samples)[1]) {
      stop("outcome.samples and ps.samples must have the same number of iterations")
    }
    if (dim(outcome.samples)[2] != dim(ps.samples)[2]) {
      stop("outcome.samples and ps.samples must have the same number of chains")
    }

    # Check dimensions
    expected_otc_cols <- ncol(Z.lm)   # intercept + covariates
    expected_ps_cols  <- ncol(Z.ps)   # intercept + covariates

    if (dim(outcome.samples)[3] != expected_otc_cols) {
      stop(sprintf("outcome.samples must have %d parameters (intercept + %d covariates), but has %d",
                   expected_otc_cols, pp.lm, dim(outcome.samples)[3]))
    }

    if (dim(ps.samples)[3] != expected_ps_cols) {
      stop(sprintf("ps.samples must have %d parameters (intercept + %d covariates), but has %d",
                   expected_ps_cols, pp.ps, dim(ps.samples)[3]))
    }

    keep      <- external_keep(dim(outcome.samples)[1])
    draws.otc <- outcome.samples[keep, , , drop = FALSE]
    draws.ps  <- ps.samples[keep, , , drop = FALSE]

  } else {
    # Generate samples using internal model functions
    say("Generating posterior samples using provided models...")

    if (!is.null(chosen_outcome) || !is.null(chosen_ps)) {
      say("Sampling with ", chosen_outcome %||% "the supplied outcome model",
          " for the outcome (", link, " link, ", outcome.prior, " prior) and ",
          chosen_ps %||% "the supplied propensity score model",
          " for the propensity score (", ps.link, " link, ", ps.prior,
          " prior).")
    }

    # Prepare arguments for outcome model
    outcome.args <- list(Y = Y, X = X.lm, mc = mc, chains = chains)
    # Starting values, when the caller gave any. drbayes_select() has always
    # forwarded these; this did not, so control$init was validated and then
    # dropped, and a sampler that mixes badly could not be helped.
    if (has_formal(outcome.model, "init")) {
      outcome.args$init <- control_init_for(control, "outcome")
    }
    # A horseshoe outcome model must not shrink the treatment effect. Which
    # columns those are is known here from the formula, and nowhere else.
    if ("unshrunk" %in% names(formals(outcome.model)) &&
        !("unshrunk" %in% names(outcome.priors))) {
      outcome.args$unshrunk <- md$treatment_columns
    }
    outcome.args <- c(outcome.args, outcome.priors)

    # Prepare arguments for propensity score model
    ps.args <- list(Y = A, X = X.ps, mc = mc, chains = chains)
    if (has_formal(ps.model, "init")) {
      ps.args$init <- control_init_for(control, "ps")
    }
    ps.args <- c(ps.args, ps.priors)

    # Call models with additional arguments
    post.otc <- do.call(outcome.model, outcome.args)
    post.ps  <- do.call(ps.model, ps.args)

    # Discard the burn-in and thin. The retained indices are the same ones the
    # external-samples branch keeps, so feeding a sampler's own draws back in
    # through outcome.samples reproduces this branch exactly.
    keep      <- seq(bn + 1L, mc, by = thin)
    draws.otc <- post.otc[keep, , , drop = FALSE]
    draws.ps  <- post.ps[keep, , , drop = FALSE]

    # The full chains are thin times larger than what has just been taken from
    # them and are never read again, so they go now rather than at the end of
    # the call, before the sequential Monte Carlo sweep claims the memory.
    rm(post.otc, post.ps)
  }

  # ============================================================================
  # 5b. Convergence diagnostics
  # ============================================================================

  # Posterior coupling tilts particles from BOTH posteriors, so the doubly
  # robust estimate is only as trustworthy as the worse of the two. Check them
  # before spending the sequential Monte Carlo sweep on draws that cannot
  # support it. Vehtari et al. (2021).
  diagnostics <- match.arg(diagnostics, c("warn", "error", "none"))
  diag_table  <- NULL

  if (diagnostics != "none") {
    # The nonparametric path has no outcome draws to diagnose: whatever fitted
    # the means was run elsewhere and is the caller's to check.
    drawn <- list(outcome = draws.otc, `propensity score` = draws.ps)
    drawn <- drawn[!vapply(drawn, is.null, logical(1))]

    diag_each <- lapply(drawn, convergence_diagnostics)
    for (nm in names(diag_each)) diag_each[[nm]]$model <- nm

    diag_table <- do.call(rbind, unname(diag_each))
    diag_table <- diag_table[, c("model", setdiff(names(diag_table), "model"))]

    failures <- lapply(diag_each, flag_convergence_failures,
                       rhat_max = rhat_max, ess_min = ess_min)

    messages <- character(0)
    for (nm in names(failures)) {
      f <- failures[[nm]]
      if (f$ok) next
      d <- diag_each[[nm]]
      parts <- character(0)
      if (length(f$rhat) > 0) {
        parts <- c(parts, paste0("max R-hat ", format(max(d$rhat, na.rm = TRUE), digits = 4),
                                 " (", toString(f$rhat), ")"))
      }
      if (length(f$ess) > 0) {
        parts <- c(parts, paste0("min ESS ",
                                 format(min(c(d$ess_bulk, d$ess_tail), na.rm = TRUE), digits = 4),
                                 " (", toString(f$ess), ")"))
      }
      messages <- c(messages, paste0(nm, " model: ", paste(parts, collapse = "; ")))
    }

    if (length(messages) > 0) {
      unreliable <- any(vapply(failures, function(f) isTRUE(f$rhat_unreliable), logical(1)))
      advice <- if (unreliable) {
        paste("An effective sample size below", ess_min,
              "also makes R-hat itself unreliable, so a small R-hat here is not",
              "evidence of convergence. Raise mc, or thin less.")
      } else {
        "Raise mc, or use more chains from more dispersed starting values."
      }
      text <- paste0("The posterior draws do not meet the convergence criteria ",
                     "(R-hat < ", rhat_max, ", ESS > ", ess_min, "). ",
                     paste(messages, collapse = " | "), ". ", advice,
                     " Inspect $diagnostics.")
      if (diagnostics == "error") stop(text) else warning(text)
    }
  }

  # The sequential Monte Carlo sweep works on pooled draws: the chain structure
  # has served its purpose once convergence has been assessed. The draws are
  # already stored iteration by chain by parameter in column-major order, so
  # reading the array as a matrix stacks the chains under one another without
  # moving a single value.
  betas.ps  <- matrix(draws.ps, ncol = dim(draws.ps)[3],
                      dimnames = list(NULL, dimnames(draws.ps)[[3]]))
  betas.otc <- if (use_nonparametric) {
    NULL
  } else {
    matrix(draws.otc, ncol = dim(draws.otc)[3],
           dimnames = list(NULL, dimnames(draws.otc)[[3]]))
  }

  num_iterations <- nrow(betas.ps)

  # The tilting resamples a draw index that names a propensity score draw and a
  # row of each matrix of fitted means at once, so the two sides must be draws
  # of the same posterior sample.
  if (use_nonparametric && nrow(outcome.mu) != num_iterations) {
    stop("outcome.mu holds ", nrow(outcome.mu), " draws but the propensity ",
         "score model contributes ", num_iterations, " (", dim(draws.ps)[1],
         " iterations from each of ", dim(draws.ps)[2], " chain(s), after ",
         "discarding bn = ", bn, " and thinning by ", thin, "). Supply one ",
         "row of fitted means per retained propensity score draw.")
  }

  say(sprintf("Using %d posterior samples for analysis...", num_iterations))

  # ============================================================================
  # 6. Posterior coupling with sequential Monte Carlo
  # ============================================================================

  td <- tilting_data(Z.lm = Z.lm, Z.ps = Z.ps, A = A, Y = Y,
                     inverse_link = inverse_link,
                     ps_inverse_link = ps_inverse_link,
                     ps_formula_text = deparse(ps.formula),
                     outcome.mu = outcome.mu)

  tilting <- if (use_nonparametric) {
    run_tilting_np(betas.ps, td, control = control, tol = tol)
  } else {
    run_tilting(betas.ps, betas.otc, td, control = control,
                tol = tol, method = method)
  }

  # ============================================================================
  # 7. Compute posterior of ATE
  # ============================================================================

  # Both estimands average the fitted contrast over the empirical distribution
  # of the covariates, as in equation (3.7), so g.comp and pc differ only
  # through the tilting and remain directly comparable.
  gcomp_ate <- function(betas) {
    fitted_row_means(betas, Z.lm1, inverse_link) -
      fitted_row_means(betas, Z.lm0, inverse_link)
  }

  if (use_nonparametric) {
    # Equation (3.7) over the supplied means. Averaging each matrix over the
    # observations first leaves one value per draw, so the tilted posterior is
    # read off by indexing that vector and neither matrix is ever subset.
    post.ate.g  <- rowMeans(outcome.mu1) - rowMeans(outcome.mu0)
    post.ate.pc <- post.ate.g[tilting$idx]
  } else {
    post.ate.g  <- gcomp_ate(betas.otc)
    post.ate.pc <- gcomp_ate(tilting$betas.otc)
  }

  result <- list(
    g.comp = post.ate.g,
    pc     = post.ate.pc,
    family = family,
    link   = link,
    ps.link = ps.link,
    outcome.prior = outcome.prior,
    ps.prior      = ps.prior,
    method        = method,
    control       = control,
    smc    = tilting$smc,
    diagnostics = diag_table,
    # The per-chain arrays, when asked for. A rank histogram is a statement
    # about how chains interleave, so it cannot be drawn from the tilted
    # particles, which have been resampled and have no chains left.
    draws = if (!isTRUE(control$keep_draws)) {
      NULL
    } else {
      list(outcome = draws.otc, ps = draws.ps)
    },
    # The untilted draws and the design matrices are kept alongside the tilted
    # ones because the sensitivity analysis of Algorithm 3 reweights the
    # original posterior, not the coupled one, and because recovering the model
    # frame later by re-evaluating the call would silently pick up a different
    # data set if the name has since been rebound.
    particles = if (!isTRUE(control$keep_particles)) {
      NULL
    } else if (use_nonparametric) {
      # There are no outcome coefficients to keep, so what identifies the
      # tilted posterior is the resampled draw index. model_data holds a
      # reference to outcome.mu rather than a copy of it.
      list(ps = tilting$betas.ps, ps_untilted = betas.ps,
           index = tilting$idx, model_data = td)
    } else {
      list(outcome = tilting$betas.otc, ps = tilting$betas.ps,
           outcome_untilted = betas.otc, ps_untilted = betas.ps,
           model_data = td, Z.lm1 = Z.lm1, Z.lm0 = Z.lm0)
    },
    call   = call_info,
    data_info = data_info
  )

  class(result) <- c("DRBayes", "list")
  return(result)
}


#' Average the fitted outcome mean over the observations, a block of draws at
#' a time
#'
#' The g-computation of equation (3.7) needs one number per posterior draw, the
#' fitted mean averaged over the observations, and the obvious way to get it
#' holds every fitted value in memory first: at 8,000 draws and 2,000
#' observations that is a 122 MB matrix, and a fit forms four of them, one at
#' each treatment value for each of the untilted and the tilted draws. Only the
#' row means survive, so the fitted values are formed a block of draws at a
#' time and the block is released once its means have been taken.
#'
#' Each row mean is an average of the same `n` fitted values, taken in the same
#' order, however the draws are grouped, so the answer is the one the whole
#' matrix would have given, to the last bit. Splitting by observations instead
#' would not have that property, since it would re-associate the sum inside
#' each mean.
#'
#' The block cannot be replaced by a product against the column means of `Z`.
#' That shortcut holds only for the identity link, where the mean of the fitted
#' values is the fitted value at the mean; a logit or probit link does not
#' commute with the average over observations, so the elementwise transform has
#' to be applied to the fitted values themselves.
#'
#' @param betas Draws by parameters matrix of outcome model coefficients.
#' @param Z Design matrix, observations by parameters, at the treatment value
#'   the fitted means are wanted for.
#' @param inverse_link Inverse link function of the outcome model.
#'
#' @return A numeric vector with one entry per draw.
#'
#' @keywords internal
#' @noRd
fitted_row_means <- function(betas, Z, inverse_link) {
  n_draws <- nrow(betas)

  # Enough draws per block to keep the working matrix near 16 MB, so that the
  # memory a block costs does not grow with the sample size.
  block <- max(1L, min(n_draws, as.integer(2e6 %/% max(nrow(Z), 1L))))

  out  <- numeric(n_draws)
  from <- 1L
  while (from <= n_draws) {
    rows <- seq.int(from, min(from + block - 1L, n_draws))
    out[rows] <-
      rowMeans(inverse_link(tcrossprod(betas[rows, , drop = FALSE], Z)))
    from <- from + block
  }
  names(out) <- rownames(betas)
  out
}


#' Check the fitted means of the nonparametric outcome path
#'
#' These three matrices stand in for the whole outcome model, and nothing
#' downstream can tell a mistake in them from a result: a matrix of the wrong
#' width pairs fitted means with the wrong observations, and one holding a
#' missing value turns the moment condition into `NaN` at that draw. The
#' matching of rows to propensity score draws is checked later, once the
#' burn-in and the thinning have fixed how many of those there are.
#'
#' @param mu,mu1,mu0 The matrices passed to [drbayes_pc()].
#' @param n Observations the models are fitted on, after `na.action`.
#' @param n_data Rows of `data`, named in the error when the two differ.
#'
#' @return `NULL`, invisibly. Called for the errors it raises.
#'
#' @keywords internal
#' @noRd
check_fitted_means <- function(mu, mu1, mu0, n, n_data) {
  given <- list(outcome.mu = mu, outcome.mu1 = mu1, outcome.mu0 = mu0)

  absent <- names(given)[vapply(given, is.null, logical(1))]
  if (length(absent) > 0L) {
    stop("The nonparametric outcome path needs all three of outcome.mu, ",
         "outcome.mu1 and outcome.mu0, and ", toString(absent), " ",
         if (length(absent) == 1L) "is" else "are", " missing. The moment ",
         "condition (3.4) uses the fitted means at the observed treatment and ",
         "the g-computation (3.7) uses those at A = 1 and at A = 0, so each ",
         "is needed on its own.")
  }

  for (nm in names(given)) {
    x <- given[[nm]]
    if (!is.matrix(x) || !is.numeric(x)) {
      stop(nm, " must be a numeric matrix of posterior draws by observations, ",
           "but is ", class(x)[1L], ". Coerce a data frame or a draws object ",
           "with as.matrix().")
    }
    if (length(x) == 0L) {
      stop(nm, " is empty. It must hold one row per posterior draw and one ",
           "column per observation used.")
    }
    # min() and max() see a missing, undefined or infinite entry while walking
    # the matrix in place. is.finite(), or range(), would copy one that can run
    # to tens of megabytes.
    if (!is.finite(min(x)) || !is.finite(max(x))) {
      stop(nm, " holds missing or infinite fitted means, which would leave ",
           "the moment condition (3.4) undefined at those draws. Drop the ",
           "draws or the observations they come from before supplying it.")
    }
  }

  shape <- vapply(given, function(x) paste(dim(x), collapse = " by "),
                  character(1))
  if (length(unique(shape)) > 1L) {
    stop("outcome.mu, outcome.mu1 and outcome.mu0 must have the same ",
         "dimensions, since they are the same draws evaluated at different ",
         "treatments, but are ",
         toString(paste0(names(shape), " ", shape)), ".")
  }

  if (ncol(mu) != n) {
    stop("outcome.mu has ", ncol(mu), " columns but the models are fitted on ",
         n, " observation", if (n == 1L) "" else "s",
         if (n != n_data) {
           paste0(" (", n_data, " rows of data, ", n_data - n,
                  " dropped for missing values by na.action)")
         },
         ". Supply one column per observation used, in the order the rows of ",
         "data are in.")
  }

  invisible(NULL)
}
