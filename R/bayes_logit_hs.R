#' Bayesian Logistic Regression with Horseshoe Prior
#'
#' Fits a Bayesian logistic regression model using the horseshoe prior for
#' regularization. This function is designed to work as an outcome model or
#' propensity score model in the \code{\link{drbayes_pc}} framework for causal inference.
#'
#' @param Y Numeric vector of binary response variables (0 or 1).
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
#'   can reveal sensitivity to where they were started: the maximum likelihood
#'   estimate is computed with \code{\link[stats]{glm.fit}} under the binomial
#'   family with a logit link, chain one starts there, and each further chain
#'   starts at a normal draw centred on it with a standard deviation of four
#'   times that estimate's own standard error. If the fit is unusable, for
#'   instance under separation or a rank deficient design, the starting values
#'   are drawn from a normal on the scale of the prior; the horseshoe itself is
#'   not used for this, since its half-Cauchy tails would place chains at
#'   values no amount of sampling would recover from.
#' @param beta0_prior Numeric. Prior precision for the intercept term. Default is 1/100.
#' @param unshrunk Integer column indices into \code{X} whose coefficients are
#'   exempt from the horseshoe prior and instead receive the normal prior set by
#'   \code{unshrunk_prior}. The default \code{integer(0)} shrinks every
#'   coefficient, which is the plain horseshoe and the right choice for a
#'   propensity score model. For an outcome model in a causal analysis, exempt
#'   the treatment and every interaction with it: shrinking those attenuates the
#'   treatment effect. \code{\link{drbayes_pc}} works these columns out from the
#'   formula and passes them for you.
#' @param unshrunk_prior Numeric. Prior precision for the coefficients named in
#'   \code{unshrunk}. Default is 1/100, i.e. a normal prior with variance 100.
#' @param p0 Integer. Guess at how many coefficients are non-zero, used only when
#'   \code{tau_prior} is NULL. Default is 5. See \code{tau_prior}.
#' @param tau_prior Numeric. Prior parameter for the global shrinkage parameter tau.
#'   If NULL, defaults to an adaptive value based on problem dimension.
#'
#' @return An array of dimension \code{mc} by \code{chains} by
#'   \code{ncol(X)+1} containing posterior samples for the logistic regression
#'   coefficients, indexed by iteration, then chain, then parameter. This is
#'   the layout the posterior package calls a draws_array. The third dimension
#'   is named with "(Intercept)" followed by the column names of X, or "X1",
#'   "X2", etc. when X has none. The first two dimensions are unnamed.
#'
#' @details
#' The horseshoe prior provides strong regularization for small coefficients while
#' allowing large coefficients to remain relatively unshrunk. The model uses the
#' Polya-Gamma data augmentation scheme for efficient sampling:
#'
#' \deqn{Y_i | \theta \sim \mathrm{Bernoulli}(\mathrm{logit}^{-1}(\beta_0 + X_i \beta))}
#' \deqn{\beta_0 \sim N(0, 1/c_0)}
#' \deqn{\beta_j | \lambda_j, \tau \sim N(0, \tau^2 \lambda_j^2)}
#' \deqn{\lambda_j \sim C^{+}(0, 1)}
#' \deqn{\tau \sim C^{+}(0, s_\tau)}
#'
#' where \eqn{c_0} is the prior precision set by \code{beta0_prior},
#' \eqn{s_\tau} is the global shrinkage scale set by \code{tau_prior}, and
#' \eqn{C^{+}(0, s)} denotes the half-Cauchy distribution with scale \eqn{s}.
#'
#' The data augmentation introduces auxiliary variables:
#' \deqn{\omega_i | \theta \sim \mathrm{PG}(1, \beta_0 + X_i \beta)}
#'
#' where PG denotes the Polya-Gamma distribution.
#'
#' The default \code{tau_prior} is set adaptively:
#' - If p <= 5: tau_prior = 1
#' - If p > 5: tau_prior = 5/(p-5) * 2/sqrt(n)
#'
#' where p is the number of coefficients under the horseshoe, that is
#' `ncol(X)` less the length of `unshrunk`, and n is the sample size. The two
#' coincide only when nothing is exempt, and exempting the treatment is the
#' recommended use, so they usually differ.
#'
#' @references
#' Carvalho, C. M., Polson, N. G., & Scott, J. G. (2010). The horseshoe estimator
#' for sparse signals. Biometrika, 97(2), 465-480.
#'
#' Polson, N. G., Scott, J. G., & Windle, J. (2013). Bayesian inference for
#' logistic models using Polya-Gamma latent variables. Journal of the American
#' Statistical Association, 108(504), 1339-1349.
#'
#' Piironen, J., & Vehtari, A. (2017). Sparsity information and regularization
#' in the horseshoe and other shrinkage priors. Electronic Journal of Statistics,
#' 11(2), 5018-5051.
#'
#' @examples
#' # Generate synthetic data: two covariates matter, eighteen do not
#' set.seed(1)
#' n <- 300
#' p <- 20
#' X <- matrix(rnorm(n * p), n, p)
#' theta_true <- c(1, 1.5, -1.5, rep(0, p - 2))  # Sparse coefficients
#' prob <- plogis(theta_true[1] + drop(X %*% theta_true[-1]))
#' Y <- rbinom(n, 1, prob)
#'
#' # Fit Bayesian logistic regression with horseshoe prior
#' draws <- bayes_logit_hs(Y = Y, X = X, mc = 800, chains = 4)
#' dim(draws)  # iterations, chains, parameters
#'
#' # Posterior means, pooling the chains
#' theta_hat <- apply(draws, 3, mean)
#' round(theta_hat[1:3], 2)
#'
#' # Credible intervals for the intercept and the two real signals
#' apply(draws[, , 1:3], 3, quantile, c(0.025, 0.975))
#'
#' # Prediction for new data, from the pooled draws
#' posterior_samples <- apply(draws, 3, as.vector)
#' X_new <- matrix(rnorm(10 * p), 10, p)
#' Z_new <- cbind(1, X_new)
#' prob_samples <- plogis(tcrossprod(posterior_samples, Z_new))
#' prob_mean <- colMeans(prob_samples)
#'
#' # Compare with bayes_logit, which shrinks nothing. Both find the two
#' # signals; the horseshoe holds the eighteen null coefficients several times
#' # closer to zero.
#' draws_standard <- bayes_logit(Y = Y, X = X, mc = 800, chains = 4)
#' round(cbind(
#'   Horseshoe = apply(draws, 3, mean),
#'   Standard  = apply(draws_standard, 3, mean)
#' )[1:3, ], 2)
#' range(theta_hat[4:(p + 1)])  # the null coefficients under the horseshoe
#' range(apply(draws_standard, 3, mean)[4:(p + 1)])  # and without it
#'
#' @seealso \code{\link{drbayes_pc}}, \code{\link{bayes_probit_hs}}, \code{\link{bayes_logit}}
#'
#' @importFrom extraDistr rinvgamma
#' @importFrom pgdraw pgdraw
#' @importFrom stats rnorm
#'
#' @export
bayes_logit_hs <- function(Y, X, mc = 5000, chains = 4L, init = NULL,
                     unshrunk = integer(0),
                     beta0_prior = 1/100, unshrunk_prior = 1/100,
                     tau_prior = NULL, p0 = 5) {

  X      <- validate_sampler_inputs(Y, X, mc, binary = TRUE)
  chains <- validate_chains(chains)
  Z      <- cbind(1, X)
  nn     <- nrow(X)
  pp     <- ncol(X)

  unshrunk <- validate_unshrunk(unshrunk, pp)
  shrunk   <- setdiff(seq_len(pp), unshrunk)
  n_shrunk <- length(shrunk)
  if (n_shrunk == 0L) {
    stop("Every column of X is listed in unshrunk, so there is nothing for the ",
         "horseshoe prior to shrink. Use bayes_logit instead.")
  }
  shrunk_theta <- shrunk + 1L

  tau_prior <- resolve_tau_prior(tau_prior, p = n_shrunk, n = nn, p0 = p0,
                                 sigma = 2)

  # Assembling the precision by index rather than through diag<- means its
  # length check no longer stands behind these, so they are checked here: a
  # vector of the wrong length would otherwise be recycled and silently
  # truncated to its first element.
  beta0_prior    <- validate_hs_precision(beta0_prior, 1 / 100, "beta0_prior")
  unshrunk_prior <- validate_hs_precision(unshrunk_prior, 1 / 100,
                                          "unshrunk_prior")

  precision <- c(beta0_prior, numeric(pp))
  precision[unshrunk + 1L] <- unshrunk_prior

  # Scale of a starting value drawn from the prior. The shrunk coefficients use
  # the global shrinkage scale as a normal standard deviation rather than the
  # horseshoe itself, whose half-Cauchy tails would put a chain somewhere it
  # could not sample its way back from.
  prior_sd <- c(1 / sqrt(beta0_prior), rep(tau_prior, pp))
  prior_sd[unshrunk + 1L] <- 1 / sqrt(unshrunk_prior)

  inits <- resolve_inits(init, chains, Y, Z, family = "binomial",
                         link = "logit", prior_sd = prior_sd)

  post.theta <- draw_storage(mc, chains, X)

  iota <- Y - 0.5
  Ztiota <- drop(crossprod(Z, iota))

  # The weighted cross product Z' diag(w) Z has to be rebuilt at every
  # iteration, because w does, and at any appreciable sample size it is the
  # largest piece of arithmetic in the scan. It is written as t(Z) %*% (w * Z)
  # rather than as crossprod(Z, w * Z): the same product of the same numbers,
  # summed over the observations in the same order, but reaching the general
  # matrix multiply with neither operand transposed, a kernel the reference
  # BLAS that R ships runs roughly twice as fast. The transpose itself is taken
  # once, here, and not mc times inside the loop.
  tZ <- t(Z)

  # Only the diagonal of the precision moves with the shrinkage parameters, and
  # the diagonal always sits at the same places. diag<-() rebuilds that index
  # vector on every call, which costs several times the addition it is there to
  # perform, so the positions are worked out once.
  diag_index <- seq.int(1L, by = pp + 2L, length.out = pp + 1L)

  # The chains are independent and could run in parallel, but they run in a
  # plain loop on purpose: the coupling step of drbayes_pc accounts for some
  # 98% of a call, so worker startup and data copying would cost more than the
  # sampling they replace.
  for (chain in seq_len(chains)) {
    theta   <- inits[[chain]]
    tau2    <- 1
    lambda2 <- rep(1, n_shrunk)

    for (iter in seq_len(mc)) {
      w <- pgdraw::pgdraw(1, drop(Z %*% theta))

      precision[shrunk_theta] <- 1 / (tau2 * lambda2)
      Q_theta <- tZ %*% (w * Z)
      Q_theta[diag_index] <- Q_theta[diag_index] + precision

      theta <- draw_normal_precision(Q_theta, Ztiota, iter)

      beta_shrunk <- theta[shrunk_theta]

      nu       <- extraDistr::rinvgamma(n_shrunk, 1, 1 + 1 / lambda2)
      b_lambda <- 1 / nu + beta_shrunk^2 / (2 * tau2)
      lambda2  <- pmax(extraDistr::rinvgamma(n_shrunk, 1, b_lambda), 1e-12)

      xi    <- extraDistr::rinvgamma(1, 1, 1 / tau_prior^2 + 1 / tau2)
      b_tau <- 1 / xi + sum(beta_shrunk^2 / lambda2) / 2
      tau2  <- max(extraDistr::rinvgamma(1, (n_shrunk + 1) / 2, b_tau), 1e-12)

      post.theta[iter, chain, ] <- theta
    }
  }

  post.theta
}
