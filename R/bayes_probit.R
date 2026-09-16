#' Bayesian Probit Regression using Latent Variable Approach
#'
#' @title Bayesian Probit Regression with Latent Variables
#'
#' @description
#' Performs Bayesian probit regression using the latent variable data augmentation
#' scheme proposed by Albert and Chib (1993). This method provides an exact and
#' efficient MCMC algorithm for binary probit regression by introducing latent
#' Gaussian variables that render the likelihood conditionally Gaussian.
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
#'   family with a probit link, chain one starts there, and each further chain
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
#' The function implements the Albert and Chib (1993) latent variable approach
#' which introduces continuous latent variables W_i such that:
#' \itemize{
#'   \item Y_i = 1 if W_i > 0, Y_i = 0 if W_i <= 0
#'   \item W_i ~ N(X_i^T * theta, 1)
#' }
#'
#' The algorithm alternates between:
#' \enumerate{
#'   \item Sampling latent variables W_i from truncated normal distributions
#'   \item Sampling regression coefficients theta from a multivariate normal distribution
#' }
#'
#' For observations with Y_i = 1, W_i is sampled from N(X_i^T * theta, 1) truncated
#' to be positive. For Y_i = 0, W_i is sampled from the same distribution truncated
#' to be negative or zero.
#'
#' The prior for regression coefficients is theta ~ N(0, V), where V^(-1) is
#' specified by \code{theta_prior}. The default prior is relatively uninformative
#' with precision 0.01 for each coefficient.
#'
#' @note
#' The function requires the \code{RcppTN} package for sampling truncated
#' normal random variables.
#'
#' @examples
#' # Simulate binary probit regression data
#' set.seed(123)
#' n <- 200
#' X1 <- rnorm(n)
#' X2 <- rnorm(n)
#' X <- cbind(X1, X2)
#'
#' # True coefficients: intercept = 0.5, X1 = 1.2, X2 = -0.8
#' eta <- 0.5 + 1.2 * X1 - 0.8 * X2
#' prob <- pnorm(eta)  # probit transformation
#' Y <- rbinom(n, 1, prob)
#'
#' # Fit Bayesian probit regression
#' draws <- bayes_probit(Y, X, mc = 1000, chains = 4)
#' dim(draws)  # iterations, chains, parameters
#'
#' # Summarize posterior, pooling the chains
#' apply(draws, 3, mean)  # posterior means
#' apply(draws, 3, sd)  # posterior standard deviations
#'
#' # Compare with classical probit regression
#' classical_fit <- glm(Y ~ X, family = binomial(link = "probit"))
#' cbind(
#'   Bayesian = apply(draws, 3, mean),
#'   Classical = coef(classical_fit)
#' )
#'
#' @references
#' Albert, J. H., & Chib, S. (1993). Bayesian analysis of binary and polychotomous
#' response data. Journal of the American Statistical Association, 88(422), 669-679.
#'
#' @seealso
#' \code{\link[stats]{glm}} for classical probit regression,
#' \code{\link{bayes_logit}} for Bayesian logistic regression,
#' \code{\link[RcppTN]{rtn}} for truncated normal random number generation.
#'
#' @export
#' @importFrom RcppTN rtn
#' @importFrom stats rnorm

bayes_probit <- function(Y, X, mc = 5000, chains = 4L, init = NULL,
                     theta_prior = NULL) {

  X      <- validate_sampler_inputs(Y, X, mc, binary = TRUE)
  chains <- validate_chains(chains)
  Z      <- cbind(1, X)
  nn     <- nrow(X)
  pp     <- ncol(X)

  prior <- validate_prior_precision(theta_prior, pp + 1L)

  inits <- resolve_inits(init, chains, Y, Z, family = "binomial",
                         link = "probit",
                         prior_sd = prior_sd_from_precision(prior))

  lower  <- ifelse(Y == 1, 0, -Inf)
  upper  <- ifelse(Y == 1, Inf, 0)
  sd_one <- rep(1, nn)

  post.theta <- draw_storage(mc, chains, X)

  # The latent variance is fixed at 1, so the coefficient precision and its
  # Cholesky factor do not change between iterations, nor between chains.
  Q_theta <- crossprod(Z) + prior
  U_theta <- chol(Q_theta)

  # t(Z) %*% W is the same product as crossprod(Z, W), summed over the
  # observations in the same order, but it reaches the BLAS with neither
  # operand transposed, and the reference BLAS that R ships runs that kernel
  # about twice as fast. The transpose is taken once, here, and not mc times
  # inside the loop.
  tZ <- t(Z)

  # The chains are independent and could run in parallel, but they run in a
  # plain loop on purpose: the coupling step of drbayes_pc accounts for some
  # 98% of a call, so worker startup and data copying would cost more than the
  # sampling they replace.
  for (chain in seq_len(chains)) {
    theta <- inits[[chain]]

    for (iter in seq_len(mc)) {
      # Albert and Chib (1993) latent variables.
      W <- RcppTN::rtn(.mean = drop(Z %*% theta), .sd = sd_one,
                       .low = lower, .high = upper)

      b_theta <- drop(tZ %*% W)
      theta   <- drop(backsolve(U_theta,
                                backsolve(U_theta, b_theta, transpose = TRUE) +
                                  stats::rnorm(pp + 1L)))

      post.theta[iter, chain, ] <- theta
    }
  }

  post.theta
}
