#' Bayesian Linear Regression with Horseshoe Prior
#'
#' Fits a Bayesian linear regression model using the horseshoe prior for
#' regularization. This function is designed to work as an outcome model
#' in the \code{\link{drbayes_pc}} framework for causal inference.
#'
#' @param Y Numeric vector of response variables (continuous outcomes).
#' @param X Matrix or data.frame of predictor variables. The intercept is
#'   automatically added, so do not include a column of ones.
#' @param mc Integer. Number of MCMC iterations. Default is 5000.
#' @param chains A positive integer, the number of Markov chains to run
#'   (default: 4). R-hat compares chains, so at least four are recommended.
#' @param init Starting values for the coefficient vector, either NULL
#'   (default) or a list of length \code{chains} whose elements are numeric
#'   vectors of length \code{ncol(X) + 1}, the intercept followed by the
#'   coefficients of the columns of X. When NULL the starting values are
#'   generated from the data and deliberately overdispersed, so that the chains
#'   can reveal sensitivity to where they were started: the ordinary least
#'   squares estimate is computed with \code{\link[stats]{lm.fit}}, chain one
#'   starts there, and each further chain starts at a normal draw centred on it
#'   with a standard deviation of four times that estimate's own standard
#'   error. If the least squares fit is unusable, for instance because the
#'   design is rank deficient, the starting values are drawn from a normal on
#'   the scale of the prior; the horseshoe itself is not used for this, since
#'   its half-Cauchy tails would place chains at values no amount of sampling
#'   would recover from. The scan below draws the coefficients first, from a
#'   full conditional that does not depend on their previous value, so a chain
#'   is started by drawing every variance component from its own full
#'   conditional given these coefficients: a starting value far from the data
#'   begins its chain at a large residual sum of squares, and so at a large
#'   error variance. Chains started in different places therefore stay
#'   different for as long as it takes them to converge, which is what R-hat
#'   needs in order to say anything at all (Vehtari et al., 2021).
#' @param beta0_prior Numeric. Prior precision for the intercept term, on the
#'   absolute scale of the data. NULL (default) asks for a precision scaled to
#'   the data, described under "Details".
#' @param unshrunk Integer column indices into \code{X} whose coefficients are
#'   exempt from the horseshoe prior and instead receive the normal prior set by
#'   \code{unshrunk_prior}. The default \code{integer(0)} shrinks every
#'   coefficient, which is the plain horseshoe and the right choice for a
#'   propensity score model. For an outcome model in a causal analysis, exempt
#'   the treatment and every interaction with it: shrinking those attenuates the
#'   treatment effect. \code{\link{drbayes_pc}} works these columns out from the
#'   formula and passes them for you.
#' @param unshrunk_prior Numeric. Prior precision for the coefficients named in
#'   \code{unshrunk}, on the absolute scale of the data. NULL (default) asks for
#'   a precision scaled to the data, described under "Details".
#' @param p0 Integer. Guess at how many coefficients are non-zero, used only
#'   when \code{tau_prior} is NULL. Default is 5. See \code{tau_prior}.
#' @param sigma_prior Numeric vector of length 2. Inverse gamma prior parameters
#'   (shape, scale) for the error variance. Default is c(1, 1).
#' @param tau_prior Numeric. Scale of the half-Cauchy prior on the global
#'   shrinkage parameter tau, which is itself free of the units of Y because
#'   the horseshoe coefficient scales carry a factor of sigma. If NULL,
#'   defaults to an adaptive value based on problem dimension.
#'
#' @return An array of dimension \code{mc} by \code{chains} by
#'   \code{ncol(X)+1} containing posterior samples for the regression
#'   coefficients, indexed by iteration, then chain, then parameter. This is
#'   the layout the posterior package calls a draws_array. The third dimension
#'   is named with "(Intercept)" followed by the column names of X, or "X1",
#'   "X2", etc. when X has none. The first two dimensions are unnamed.
#'
#' @details
#' The horseshoe prior is a continuous shrinkage prior that provides strong
#' regularization for small coefficients while allowing large coefficients
#' to remain relatively unshrunk. In the context of causal inference, this
#' function applies different priors to different types of coefficients:
#'
#' \deqn{Y | \beta, \sigma^2 \sim N(X\beta, \sigma^2 I)}
#' \deqn{\beta_0 \sim N(0, 1/c_0)}
#' \deqn{\beta_1 \sim N(0, 1/c_1)}
#' \deqn{\beta_j | \lambda_j, \tau, \sigma^2 \sim
#'   N(0, \sigma^2 \tau^2 \lambda_j^2), \quad j = 2, \ldots, p}
#' \deqn{\lambda_j \sim C^{+}(0, 1), \quad j = 2, \ldots, p}
#' \deqn{\tau \sim C^{+}(0, t)}
#' \deqn{\sigma^2 \sim \mathrm{InvGamma}(a, b)}
#'
#' Here \eqn{\beta_1} stands for the coefficients named in \code{unshrunk},
#' which receive no shrinkage, and \eqn{\beta_2, \ldots, \beta_p} for the rest.
#' \eqn{C^{+}(0, s)} denotes the half-Cauchy distribution with scale \eqn{s}.
#' The prior precisions \eqn{c_0} and \eqn{c_1} are set by \code{beta0_prior}
#' and \code{unshrunk_prior}, the global shrinkage scale \eqn{t} is set by
#' \code{tau_prior}, and the shape and scale \eqn{(a, b)} are set by
#' \code{sigma_prior}.
#'
#' The horseshoe coefficient scales carry a factor of \eqn{\sigma}, as in
#' Carvalho, Polson and Scott (2010) and in the sampler of Makalic and Schmidt
#' (2016). That is what makes the amount of shrinkage independent of the units
#' of Y without any reference to the data: the posterior precision of those
#' coefficients is \eqn{(X'X + D)/\sigma^2} with \eqn{D} the prior precision,
#' so \eqn{D} is compared with \eqn{X'X} and not with \eqn{X'X/\sigma^2}. It
#' also leaves \eqn{\tau} free of those units, which is what the default
#' \code{tau_prior} below assumes. Those coefficients, and only those, add
#' degrees of freedom to the full conditional of \eqn{\sigma^2}, which is
#' \eqn{\mathrm{InvGamma}(a + (n + m)/2, b + (\mathrm{RSS} + Q)/2)} for \eqn{m}
#' shrunk coefficients, residual sum of squares \eqn{\mathrm{RSS}} and prior
#' quadratic form \eqn{Q = \sum_j \beta_j^2/(\tau^2\lambda_j^2)}.
#'
#' The intercept and the unshrunk coefficients are treated differently, because
#' they have no shrinkage parameter to adapt to their size. Their priors are
#' fixed variances rather than multiples of \eqn{\sigma^2}, so that they cannot
#' inflate the error variance: on the simulation of Orihara et al. (2025) the
#' intercept is about 100 and the treatment coefficient about 110 while
#' \eqn{\sigma} is 1, and a prior variance of \eqn{100\sigma^2} on each would
#' add 221 to a residual sum of squares of 499. In exchange their precisions
#' have to be set on the scale of the data, and the defaults are the ones
#' \code{\link{bayes_lm}} documents under "Prior Specification": a prior
#' standard deviation of \eqn{100(|\mathrm{mean}(Y)| + \mathrm{sd}(Y))} for the
#' intercept, which is effectively flat, and of
#' \eqn{2.5 \mathrm{sd}(Y)/s_j} for an unshrunk coefficient on a column with
#' standard deviation \eqn{s_j}. Those scales are read for every column and not
#' only for the exempt ones, since they also place the starting values and the
#' reference fit the residual sum of squares is expanded around, so data too far
#' apart in scale for the defaults to be represented stop the fit whichever
#' precisions are named here. \code{\link{bayes_lm}} documents where that is.
#'
#' The default value of \code{tau_prior} follows Piironen and Vehtari (2017):
#' with \eqn{p} coefficients under the horseshoe, \eqn{n} observations and
#' \eqn{p_0} of them guessed to be non-zero, it is
#' \eqn{p_0/(p - p_0) \times 1/\sqrt{n}}, or 1 when \eqn{p \le p_0} leaves
#' nothing to shrink towards.
#'
#' @references
#' Carvalho, C. M., Polson, N. G., & Scott, J. G. (2010). The horseshoe
#' estimator for sparse signals. Biometrika, 97(2), 465-480.
#'
#' Makalic, E., & Schmidt, D. F. (2016). A simple sampler for the horseshoe
#' estimator. IEEE Signal Processing Letters, 23(1), 179-182.
#'
#' Orihara, S., Momozaki, T., & Sugasawa, S. (2025). Bayesian Doubly Robust
#' Causal Inference via Posterior Coupling. arXiv preprint arXiv:2506.04868.
#'
#' Piironen, J., & Vehtari, A. (2017). Sparsity information and
#' regularization in the horseshoe and other shrinkage priors. Electronic
#' Journal of Statistics, 11(2), 5018-5051.
#'
#' Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Burkner, P.-C.
#' (2021). Rank-normalization, folding, and localization: an improved R-hat for
#' assessing convergence of MCMC. Bayesian Analysis, 16(2), 667-718.
#'
#' @examples
#' # Generate synthetic data for causal inference
#' set.seed(1)
#' n <- 100
#' p <- 20
#' W <- matrix(rnorm(n * (p - 1)), n, p - 1,
#'             dimnames = list(NULL, paste0("W", seq_len(p - 1))))
#' A <- rbinom(n, 1, 0.5)  # Treatment
#' X <- cbind(A = A, W)  # Treatment is the first column
#'
#' # Sparse confounding: only the first three confounders matter
#' beta <- c(1.5, 1, -1, 0.5, rep(0, p - 4))  # A, then the confounders
#' Y <- 2 + drop(X %*% beta) + rnorm(n, 0, 0.5)
#'
#' # Fit Bayesian linear regression with horseshoe prior. unshrunk = 1
#' # exempts the treatment, the first column of X, from the horseshoe:
#' # shrinking it would attenuate the effect being estimated.
#' draws <- bayes_lm_hs(Y = Y, X = X, mc = 1000, chains = 4, unshrunk = 1)
#' dim(draws)  # iterations, chains, parameters
#'
#' # Posterior means, pooling the chains: the treatment effect and the three
#' # real confounders are recovered, and the sixteen null coefficients are
#' # shrunk to a few hundredths.
#' beta_hat <- apply(draws, 3, mean)
#' round(beta_hat[1:5], 2)
#' range(beta_hat[6:(p + 1)])
#'
#' # Credible interval for the treatment effect
#' quantile(draws[, , "A"], c(0.025, 0.975))
#'
#' # Pool the chains into an (iterations * chains) by parameters matrix
#' posterior_samples <- apply(draws, 3, as.vector)
#'
#' @seealso \code{\link{drbayes_pc}}, \code{\link{bayes_logit_hs}}
#'
#' @importFrom extraDistr rinvgamma
#' @importFrom stats rnorm
#'
#' @export
bayes_lm_hs <- function(Y, X, mc = 5000, chains = 4L, init = NULL,
                  unshrunk = integer(0),
                  beta0_prior = NULL, unshrunk_prior = NULL,
                  sigma_prior = c(1, 1), tau_prior = NULL, p0 = 5) {

  X      <- validate_sampler_inputs(Y, X, mc)
  chains <- validate_chains(chains)
  Z      <- cbind(1, X)
  nn     <- nrow(X)
  pp     <- ncol(X)

  unshrunk <- validate_unshrunk(unshrunk, pp)
  shrunk   <- setdiff(seq_len(pp), unshrunk)
  n_shrunk <- length(shrunk)
  if (n_shrunk == 0L) {
    stop("Every column of X is listed in unshrunk, so there is nothing for ",
         "the horseshoe prior to shrink. Use bayes_lm instead.")
  }
  # Positions in theta = (intercept, X columns), so column j is theta[j + 1].
  shrunk_theta <- shrunk + 1L

  if (length(sigma_prior) != 2L || !is.numeric(sigma_prior) ||
      any(sigma_prior <= 0)) {
    stop("sigma_prior must be two positive numbers, the shape and scale of ",
         "the inverse gamma prior on the error variance")
  }
  a <- sigma_prior[1]
  b <- sigma_prior[2]

  tau_prior <- resolve_tau_prior(tau_prior, p = n_shrunk, n = nn, p0 = p0,
                                 sigma = NULL)

  # The priors that carry no shrinkage parameter are fixed variances, so their
  # defaults have to come from the data to stay weakly informative.
  scaled <- default_coef_precision(Y, X)
  beta0_precision <- validate_hs_precision(beta0_prior, scaled[1L],
                                           "beta0_prior")
  unshrunk_precision <- validate_hs_precision(unshrunk_prior,
                                              scaled[unshrunk + 1L],
                                              "unshrunk_prior")

  precision <- numeric(pp + 1L)
  precision[1L] <- beta0_precision
  precision[unshrunk + 1L] <- unshrunk_precision

  # Scale of a starting value drawn from the prior. The shrunk coefficients use
  # the global shrinkage scale, on the scale of the data, as a normal standard
  # deviation rather than the horseshoe itself, whose half-Cauchy tails would
  # put a chain somewhere it could not sample its way back from.
  prior_sd <- c(1 / sqrt(beta0_precision), tau_prior / sqrt(scaled[-1L]))
  prior_sd[unshrunk + 1L] <- 1 / sqrt(unshrunk_precision)

  inits <- resolve_inits(init, chains, Y, Z, family = "gaussian",
                         prior_sd = prior_sd)

  post.theta <- draw_storage(mc, chains, X)

  ZtZ <- crossprod(Z)
  ZtY <- drop(crossprod(Z, Y))
  check_cross_products(ZtZ, ZtY, dimnames(post.theta)[[3]])

  # Only the diagonal of the precision moves with the shrinkage parameters, and
  # the diagonal always sits at the same places. diag<-() rebuilds that index
  # vector on every call, which costs several times the addition it is there to
  # perform, so the positions are worked out once.
  diag_index <- seq.int(1L, by = pp + 2L, length.out = pp + 1L)

  # The shrunk coefficients have no fixed precision of their own, so the
  # reference fit borrows the scaled default for them; it only has to be a
  # point near the posterior, and positive definite so that every design has
  # one.
  ridge <- precision
  ridge[shrunk_theta] <- scaled[shrunk_theta]
  ref <- rss_reference(Y, Z, ZtZ, ZtY, diag(ridge, pp + 1L))

  # Only the horseshoe coefficients are scaled by sigma^2 under the prior, so
  # only they add degrees of freedom to the full conditional of sigma^2.
  a_sigma <- a + (nn + n_shrunk) / 2

  # The variance components, drawn from their full conditionals given the
  # coefficients. This is the tail of the Gibbs scan, and running it once
  # before the scan is what starts a chain from its own starting values.
  draw_variances <- function(theta, tau2, lambda2) {
    beta_shrunk <- theta[shrunk_theta]
    rss  <- rss_at(theta, ref)
    quad <- sum(beta_shrunk^2 / lambda2) / tau2
    # Cancellation can leave the residual sum of squares a rounding error
    # below zero when the fit is exact, and a scale has to stay positive.
    b_sigma <- max(b + (rss + quad) / 2, .Machine$double.eps)
    sigma2  <- extraDistr::rinvgamma(1, a_sigma, b_sigma)

    # The squared coefficients enter the scale updates divided by sigma^2,
    # again because that is the unit their prior scales are measured in.
    beta2_shrunk <- beta_shrunk^2 / sigma2

    nu       <- extraDistr::rinvgamma(n_shrunk, 1, 1 + 1 / lambda2)
    b_lambda <- 1 / nu + beta2_shrunk / (2 * tau2)
    lambda2  <- pmax(extraDistr::rinvgamma(n_shrunk, 1, b_lambda), 1e-12)

    xi    <- extraDistr::rinvgamma(1, 1, 1 / tau_prior^2 + 1 / tau2)
    b_tau <- 1 / xi + sum(beta2_shrunk / lambda2) / 2
    tau2  <- max(extraDistr::rinvgamma(1, (n_shrunk + 1) / 2, b_tau), 1e-12)

    list(sigma2 = sigma2, tau2 = tau2, lambda2 = lambda2)
  }

  # The chains are independent and could run in parallel, but they run in a
  # plain loop on purpose: the coupling step of drbayes_pc accounts for some
  # 98% of a call, so worker startup and data copying would cost more than the
  # sampling they replace.
  for (chain in seq_len(chains)) {
    theta <- inits[[chain]]
    state <- draw_variances(theta, tau2 = 1, lambda2 = rep(1, n_shrunk))

    for (iter in seq_len(mc)) {
      sigma2 <- state$sigma2
      # The horseshoe scale is sigma^2 tau^2 lambda_j^2, so its precision picks
      # up the same 1/sigma^2 that Z'Z does; the two fixed priors do not.
      precision[shrunk_theta] <- 1 / (sigma2 * state$tau2 * state$lambda2)
      Q_theta <- ZtZ / sigma2
      Q_theta[diag_index] <- Q_theta[diag_index] + precision

      theta <- draw_normal_precision(Q_theta, ZtY / sigma2, iter)
      state <- draw_variances(theta, state$tau2, state$lambda2)

      post.theta[iter, chain, ] <- theta
    }
  }

  post.theta
}


#' Prior precision of a coefficient the horseshoe leaves alone
#'
#' @param value A single positive number, or NULL for the data-scaled default.
#' @param default Precision used when `value` is NULL, one per coefficient.
#' @param arg Argument name, used in the error message.
#'
#' @return A numeric vector as long as `default`.
#'
#' @keywords internal
#' @noRd
validate_hs_precision <- function(value, default, arg) {
  if (is.null(value)) {
    return(default)
  }
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value <= 0) {
    stop(arg, " is a prior PRECISION, so it must be a single positive number, ",
         "or NULL for a default scaled to the data. Pass 1/variance.")
  }
  rep_len(as.numeric(value), length(default))
}
