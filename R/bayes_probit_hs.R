#' Bayesian Probit Regression with Horseshoe Prior
#'
#' Fits a Bayesian probit regression model using the horseshoe prior for
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
#'   family with a probit link, chain one starts there, and each further chain
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
#'   \code{ncol(X)+1} containing posterior samples for the probit regression
#'   coefficients, indexed by iteration, then chain, then parameter. This is
#'   the layout the posterior package calls a draws_array. The third dimension
#'   is named with "(Intercept)" followed by the column names of X, or "X1",
#'   "X2", etc. when X has none. The first two dimensions are unnamed.
#'
#' @details
#' The horseshoe prior provides strong regularization for small coefficients while
#' allowing large coefficients to remain relatively unshrunk. The model uses the
#' Albert and Chib (1993) latent variable approach for efficient sampling:
#'
#' \deqn{Y_i | \theta \sim \mathrm{Bernoulli}(\Phi(\beta_0 + X_i \beta))}
#' \deqn{W_i | \theta \sim N(\beta_0 + X_i \beta, 1) \mbox{ with } Y_i = \mathbf{1}(W_i > 0)}
#' \deqn{\beta_0 \sim N(0, 1/c_0)}
#' \deqn{\beta_j | \lambda_j, \tau \sim N(0, \tau^2 \lambda_j^2)}
#' \deqn{\lambda_j \sim C^{+}(0, 1)}
#' \deqn{\tau \sim C^{+}(0, s_\tau)}
#'
#' Here \eqn{c_0} is the prior precision for the intercept set by
#' \code{beta0_prior}, \eqn{s_\tau} is the global shrinkage scale set by
#' \code{tau_prior}, and \eqn{C^{+}(0, s)} denotes the half-Cauchy distribution
#' with scale \eqn{s}.
#'
#' The data augmentation introduces latent variables. For \eqn{Y_i = 1} the
#' latent variable is truncated to \eqn{(0, \infty)}:
#' \deqn{W_i | \theta \sim TN(X_i \theta, 1, 0, \infty)}
#' and for \eqn{Y_i = 0} it is truncated to \eqn{(-\infty, 0]}:
#' \deqn{W_i | \theta \sim TN(X_i \theta, 1, -\infty, 0]}
#'
#' where TN denotes the truncated normal distribution and \eqn{\Phi} is the
#' standard normal cumulative distribution function.
#'
#' The default \code{tau_prior} is set adaptively:
#' - If p <= 5: tau_prior = 1
#' - If p > 5: tau_prior = 5/(p-5) * 1/sqrt(n)
#'
#' where p is the number of coefficients under the horseshoe, that is
#' `ncol(X)` less the length of `unshrunk`, and n is the sample size. The
#' scale factor is 1 because the Albert-Chib latent variable has unit
#' variance; the corresponding factor in [bayes_logit_hs()] is 2, which is the
#' logistic one, and the two deliberately differ.
#'
#' @references
#' Albert, J. H., & Chib, S. (1993). Bayesian analysis of binary and polychotomous
#' response data. Journal of the American Statistical Association, 88(422), 669-679.
#'
#' Carvalho, C. M., Polson, N. G., & Scott, J. G. (2010). The horseshoe estimator
#' for sparse signals. Biometrika, 97(2), 465-480.
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
#' theta_true <- c(0.5, 1, -0.8, rep(0, p - 2))  # Sparse coefficients
#' prob <- pnorm(theta_true[1] + drop(X %*% theta_true[-1]))  # Probit link
#' Y <- rbinom(n, 1, prob)
#'
#' # Fit Bayesian probit regression with horseshoe prior
#' draws <- bayes_probit_hs(Y = Y, X = X, mc = 800, chains = 4)
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
#' prob_samples <- pnorm(tcrossprod(posterior_samples, Z_new))
#' prob_mean <- colMeans(prob_samples)
#'
#' # Compare with bayes_probit, which shrinks nothing. Both find the two
#' # signals; the horseshoe holds the eighteen null coefficients several times
#' # closer to zero.
#' draws_standard <- bayes_probit(Y = Y, X = X, mc = 800, chains = 4)
#' round(cbind(
#'   Horseshoe = apply(draws, 3, mean),
#'   Standard  = apply(draws_standard, 3, mean)
#' )[1:3, ], 2)
#' range(theta_hat[4:(p + 1)])  # the null coefficients under the horseshoe
#' range(apply(draws_standard, 3, mean)[4:(p + 1)])  # and without it
#'
#' @seealso \code{\link{drbayes_pc}}, \code{\link{bayes_logit_hs}}, \code{\link{bayes_probit}}
#'
#' @importFrom extraDistr rinvgamma
#' @importFrom RcppTN rtn
#' @importFrom stats rnorm
#'
#' @export
bayes_probit_hs <- function(Y, X, mc = 5000, chains = 4L, init = NULL,
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
         "horseshoe prior to shrink. Use bayes_probit instead.")
  }
  shrunk_theta <- shrunk + 1L

  tau_prior <- resolve_tau_prior(tau_prior, p = n_shrunk, n = nn, p0 = p0,
                                 sigma = 1)

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
                         link = "probit", prior_sd = prior_sd)

  post.theta <- draw_storage(mc, chains, X)

  lower <- ifelse(Y == 1, 0, -Inf)
  upper <- ifelse(Y == 1, Inf, 0)
  sd_one <- rep(1, nn)
  ZtZ <- crossprod(Z)

  # t(Z) %*% W is the same product as crossprod(Z, W), summed over the
  # observations in the same order, but it reaches the BLAS with neither
  # operand transposed, and the reference BLAS that R ships runs that kernel
  # about twice as fast. The transpose is taken once, here, and not mc times
  # inside the loop.
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
      W <- RcppTN::rtn(.mean = drop(Z %*% theta), .sd = sd_one,
                       .low = lower, .high = upper)

      precision[shrunk_theta] <- 1 / (tau2 * lambda2)
      Q_theta <- ZtZ
      Q_theta[diag_index] <- Q_theta[diag_index] + precision

      theta <- draw_normal_precision(Q_theta, drop(tZ %*% W), iter)

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
