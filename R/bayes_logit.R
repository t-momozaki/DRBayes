#' Bayesian Logistic Regression using Polya-Gamma Data Augmentation
#'
#' @title Bayesian Logistic Regression with Polya-Gamma Latent Variables
#'
#' @description
#' Performs Bayesian logistic regression using the Polya-Gamma data augmentation
#' scheme proposed by Polson, Scott, and Windle (2013). This method provides
#' an exact and efficient MCMC algorithm for binary logistic regression by
#' introducing Polya-Gamma auxiliary variables that render the likelihood
#' conditionally Gaussian.
#'
#' @param Y A numeric vector of binary outcomes (0 or 1). Must contain only
#'   values 0 and 1.
#' @param X A matrix or data.frame of covariates. Each row corresponds to an
#'   observation and each column to a predictor variable. An intercept term
#'   will be automatically added.
#' @param mc A positive integer specifying the number of MCMC iterations
#'   (default: 5000).
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
#'   are drawn from the prior instead.
#' @param theta_prior Either a scalar or a square matrix specifying the precision
#'   for the normal prior on regression coefficients. If scalar, assumes
#'   diagonal precision matrix with all elements equal to this value.
#'   If NULL (default), a weakly informative prior with precision
#'   1/100 * I is used, where I is the identity matrix.
#'
#' @return An array of posterior samples with dimensions \code{mc} x
#'   \code{chains} x (\code{ncol(X)} + 1), indexed by iteration, then chain,
#'   then parameter. This is the layout the posterior package calls a
#'   draws_array. The third dimension is named with "(Intercept)" for the
#'   intercept term and with the column names of X, or "X1", "X2", etc. when X
#'   has none. The first two dimensions are unnamed.
#'
#' @details
#' The function implements the Polya-Gamma data augmentation strategy which
#' represents the logistic likelihood as a scale mixture of Gaussians. The
#' algorithm alternates between:
#' \enumerate{
#'   \item Sampling Polya-Gamma auxiliary variables w_i ~ PG(1, X_i^T * theta)
#'   \item Sampling regression coefficients theta from a multivariate normal distribution
#' }
#'
#' The method is exact (no approximation) and often more efficient than
#' alternative data augmentation schemes, especially for hierarchical models.
#'
#' The prior for regression coefficients is theta ~ N(0, V), where V^(-1) is
#' specified by \code{theta_prior}. The default prior is relatively uninformative
#' with precision 0.01 for each coefficient.
#'
#' @note
#' The function requires the \code{pgdraw} package for sampling Polya-Gamma
#' random variables.
#'
#' @examples
#' # Simulate binary logistic regression data
#' set.seed(123)
#' n <- 200
#' X1 <- rnorm(n)
#' X2 <- rnorm(n)
#' X <- cbind(X1, X2)
#'
#' # True coefficients: intercept = 0.5, X1 = 1.2, X2 = -0.8
#' eta <- 0.5 + 1.2 * X1 - 0.8 * X2
#' prob <- plogis(eta)  # logistic transformation
#' Y <- rbinom(n, 1, prob)
#'
#' # Fit Bayesian logistic regression
#' draws <- bayes_logit(Y, X, mc = 1000, chains = 4)
#' dim(draws)  # iterations, chains, parameters
#'
#' # Summarize posterior, pooling the chains
#' apply(draws, 3, mean)  # posterior means
#' apply(draws, 3, sd)  # posterior standard deviations
#'
#' # Pool the chains into an (iterations * chains) by parameters matrix
#' posterior_samples <- apply(draws, 3, as.vector)
#'
#' @references
#' Polson, N. G., Scott, J. G., & Windle, J. (2013). Bayesian inference for
#' logistic models using Polya-Gamma latent variables. Journal of the American
#' Statistical Association, 108(504), 1339-1349.
#'
#' @seealso
#' \code{\link[stats]{glm}} for classical logistic regression,
#' \code{\link[pgdraw]{pgdraw}} for Polya-Gamma random number generation.
#'
#' @export
#' @importFrom pgdraw pgdraw
#' @importFrom stats rnorm

bayes_logit <- function(Y, X, mc = 5000, chains = 4L, init = NULL,
                    theta_prior = NULL) {

  X      <- validate_sampler_inputs(Y, X, mc, binary = TRUE)
  chains <- validate_chains(chains)
  Z      <- cbind(1, X)
  pp     <- ncol(X)

  prior <- validate_prior_precision(theta_prior, pp + 1L)

  inits <- resolve_inits(init, chains, Y, Z, family = "binomial",
                         link = "logit",
                         prior_sd = prior_sd_from_precision(prior))

  iota       <- Y - 0.5
  Ztiota     <- drop(crossprod(Z, iota))
  post.theta <- draw_storage(mc, chains, X)

  # The weighted cross product Z' diag(w) Z has to be rebuilt at every
  # iteration, because w does, and at any appreciable sample size it is the
  # largest piece of arithmetic in the scan. It is written as t(Z) %*% (w * Z)
  # rather than as crossprod(Z, w * Z): the same product of the same numbers,
  # summed over the observations in the same order, but reaching the general
  # matrix multiply with neither operand transposed. The reference BLAS that R
  # ships runs that kernel roughly twice as fast, because its innermost loop
  # then walks down a column and vectorises, while the transposed one carries a
  # running sum it cannot. The transpose itself is taken once, here, and not mc
  # times inside the loop.
  tZ <- t(Z)

  # The chains are independent and could run in parallel, but they run in a
  # plain loop on purpose: the coupling step of drbayes_pc accounts for some
  # 98% of a call, so worker startup and data copying would cost more than the
  # sampling they replace.
  for (chain in seq_len(chains)) {
    theta <- inits[[chain]]

    for (iter in seq_len(mc)) {
      # Polya-Gamma augmentation makes the logistic likelihood conditionally
      # Gaussian (Polson, Scott and Windle 2013).
      w <- pgdraw::pgdraw(1, drop(Z %*% theta))

      Q_theta <- tZ %*% (w * Z) + prior
      theta   <- draw_normal_precision(Q_theta, Ztiota, iter)

      post.theta[iter, chain, ] <- theta
    }
  }

  post.theta
}
