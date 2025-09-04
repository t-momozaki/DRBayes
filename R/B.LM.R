#' Bayesian Linear Regression with Conjugate Priors
#'
#' @title Bayesian Linear Regression using MCMC
#'
#' @description
#' Performs Bayesian linear regression using Gibbs sampling with conjugate
#' normal-inverse-gamma priors. The model assumes Y ~ N(beta_0 1_n + X*theta, sigma^2 I) with
#' theta ~ N(0, V) and sigma^2 ~ InvGamma(a, b), where V^(-1) is the precision matrix
#' specified by the user.
#'
#' @param Y A numeric vector of continuous outcomes (response variable).
#' @param X A matrix or data.frame of covariates. Each row corresponds to an
#'   observation and each column to a predictor variable. An intercept term
#'   will be automatically added.
#' @param mc A positive integer specifying the number of MCMC iterations
#'   (default: 5000).
#' @param theta_prior Either a scalar or a square matrix specifying the precision
#'   for the normal prior on regression coefficients. If scalar, assumes
#'   diagonal precision matrix with all elements equal to this value.
#'   If NULL (default), a weakly informative prior with precision
#'   1/100 * I is used, where I is the identity matrix.
#' @param sigma_prior A numeric vector of length 2 specifying the shape (a) and
#'   rate (b) parameters for the inverse-gamma prior on the error variance sigma^2.
#'   The prior is sigma^2 ~ InvGamma(a, b) with mean b/(a-1) for a > 1.
#'   (default: c(1, 1)).
#'
#' @return A matrix of posterior samples for regression coefficients with
#'   dimensions \code{mc} x (\code{ncol(X)} + 1). Each row represents one MCMC
#'   draw of the regression coefficients. Columns are named with "(Intercept)"
#'   for the intercept term and "X1", "X2", etc. for the predictor coefficients.
#'
#' @details
#' The function implements a Gibbs sampler for Bayesian linear regression with
#' the following hierarchical model:
#' \deqn{Y_i = \beta_0 + X_i^T \beta + \varepsilon_i, \quad \varepsilon_i \sim N(0, \sigma^2)}
#' \deqn{\theta = (\beta_0, \beta)^T \sim N(0, V)}
#' \deqn{\sigma^2 \sim InvGamma(a, b)}
#'
#' The algorithm alternates between:
#' \enumerate{
#'   \item Sampling regression coefficients theta from their full conditional
#'         (multivariate normal)
#'   \item Sampling error variance sigma^2 from its full conditional (inverse-gamma)
#' }
#'
#' The conjugate structure allows for efficient closed-form updates at each step.
#'
#' @section Prior Specification:
#'
#' \strong{Regression Coefficients:} The prior theta ~ N(0, V) can be
#' specified via \code{theta_prior} which sets V^(-1) (the precision matrix):
#' \itemize{
#'   \item Scalar: Creates diagonal precision matrix V^(-1) = theta_prior * I
#'   \item Matrix: Uses the provided matrix as precision V^(-1)
#'   \item NULL: Uses default weakly informative precision 0.01 * I
#' }
#'
#' \strong{Error Variance:} The prior sigma^2 ~ InvGamma(a, b) is specified via
#' \code{sigma_prior = c(a, b)}. Common choices:
#' \itemize{
#'   \item c(1, 1): Improper uniform prior (weakly informative)
#'   \item c(0.01, 0.01): Jeffrey's prior (non-informative)
#'   \item c(2, 1): Proper prior with mean 1
#' }
#'
#' @note
#' This function requires the \code{extraDistr} package for sampling from
#' the inverse-gamma distribution. The current implementation only returns
#' posterior samples for regression coefficients.
#'
#' @examples
#' \dontrun{
#' # Simulate linear regression data
#' set.seed(123)
#' n <- 100
#' X1 <- rnorm(n)
#' X2 <- rnorm(n)
#' X <- cbind(X1, X2)
#'
#' # True coefficients: intercept = 2, X1 = 1.5, X2 = -0.8, sigma^2 = 1
#' Y <- 2 + 1.5 * X1 - 0.8 * X2 + rnorm(n, 0, 1)
#'
#' # Fit Bayesian linear regression with default priors
#' posterior_samples <- B.LM(Y, X, mc = 2000)
#'
#' # Summarize posterior
#' colMeans(posterior_samples)  # posterior means
#' apply(posterior_samples, 2, sd)  # posterior standard deviations
#' apply(posterior_samples, 2, quantile, c(0.025, 0.975))  # 95% credible intervals
#' }
#'
#' @references
#' Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., &
#' Rubin, D. B. (2013). Bayesian data analysis. CRC press.
#'
#' @seealso
#' \code{\link[stats]{lm}} for classical linear regression,
#' \code{\link[extraDistr]{rinvgamma}} for inverse-gamma distribution.
#'
#' @export
#' @importFrom extraDistr rinvgamma
#' @importFrom stats rnorm

B.LM <- function(Y, X, mc = 5000,
                 theta_prior = 1/100, sigma_prior = c(1, 1)) {

  # Input validation
  if (!is.numeric(Y)) {
    stop("Y must be a numeric vector")
  }

  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("X must be a matrix or data.frame")
  }

  # Convert to matrix and add intercept
  X <- as.matrix(X)
  Z <- cbind(1, X)

  nn <- nrow(X)
  pp <- ncol(X)

  if (length(Y) != nn) {
    stop("Number of observations in Y and X must match")
  }

  if (mc <= 0 || !is.numeric(mc)) {
    stop("mc must be a positive integer")
  }

  if (length(sigma_prior) != 2 || !is.numeric(sigma_prior) || any(sigma_prior <= 0)) {
    stop("sigma_prior must be a numeric vector of length 2 with positive values")
  }

  # Prior parameters
  a <- sigma_prior[1]
  b <- sigma_prior[2]

  # Initialize parameters
  sigma2  <- 1

  # Set default prior if not specified
  if (is.null(theta_prior)) {
    theta_prior <- diag(pp+1) * 1/100
  } else if (length(theta_prior) == 1) {
    # If scalar, create diagonal precision matrix
    theta_prior <- diag(pp+1) * theta_prior
  } else if (!is.matrix(theta_prior)) {
    stop("theta_prior must be NULL, a scalar, or a square matrix")
  } else if (nrow(theta_prior) != (pp+1) || ncol(theta_prior) != (pp+1)) {
    stop("theta_prior matrix must have dimensions (p+1) x (p+1) where p is ncol(X)")
  }

  # Storage for posterior samples
  post.theta <- matrix(NA, mc, pp+1)
  colnames(post.theta) <- c("(Intercept)", paste0("X", 1:pp))

  # MCMC sampling
  for (iter in 1:mc) {
    # Sample regression coefficients
    b_theta <- 1/sigma2 * drop( crossprod(Z, Y) )
    Q_theta <- 1/sigma2 * crossprod(Z) + theta_prior
    U_theta <- chol(Q_theta)
    theta   <- drop( base::backsolve(U_theta,
                                     base::backsolve(U_theta, b_theta, transpose=TRUE) + stats::rnorm(pp+1)) )

    # Sample error variance
    a_sigma <- a + nn/2
    b_sigma <- b + crossprod(Y-drop(tcrossprod(theta,Z))) / 2
    sigma2  <- extraDistr::rinvgamma(1, a_sigma, b_sigma)

    # Store samples
    post.theta[iter, ] <- theta
  }

  return(post.theta)
}
