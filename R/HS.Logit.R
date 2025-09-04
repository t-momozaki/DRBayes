#' Bayesian Logistic Regression with Horseshoe Prior
#'
#' Fits a Bayesian logistic regression model using the horseshoe prior for
#' regularization. This function is designed to work as an outcome model or
#' propensity score model in the \code{\link{DRBayes.PC}} framework for causal inference.
#'
#' @param Y Numeric vector of binary response variables (0 or 1).
#' @param X Matrix or data.frame of predictor variables. The intercept is
#'   automatically added, so do not include a column of ones.
#' @param mc Integer. Number of MCMC iterations. Default is 5000.
#' @param beta0_prior Numeric. Prior precision for the intercept term. Default is 1/100.
#' @param tau_prior Numeric. Prior parameter for the global shrinkage parameter tau.
#'   If NULL, defaults to an adaptive value based on problem dimension.
#'
#' @return Matrix of dimension \code{mc} by \code{ncol(X)+1} containing posterior
#'   samples for the logistic regression coefficients. The first column corresponds to
#'   the intercept, and the remaining columns to the coefficients for X.
#'
#' @details
#' The horseshoe prior provides strong regularization for small coefficients while
#' allowing large coefficients to remain relatively unshrunk. The model uses the
#' Pólya-Gamma data augmentation scheme for efficient sampling:
#'
#' \deqn{Y_i | \theta \sim \text{Bernoulli}(\text{logit}^{-1}(\beta_0 + X_i \beta))}
#' \deqn{\beta_0 \sim N(0, 1/\text{beta0_prior})}
#' \deqn{\beta_j | \lambda_j, \tau \sim N(0, \tau^2 \lambda_j^2)}
#' \deqn{\lambda_j \sim \text{Half-Cauchy}(0, 1)}
#' \deqn{\tau \sim \text{Half-Cauchy}(0, \text{tau_prior})}
#'
#' The data augmentation introduces auxiliary variables:
#' \deqn{\omega_i | \theta \sim \text{PG}(1, \beta_0 + X_i \beta)}
#'
#' where PG denotes the Pólya-Gamma distribution.
#'
#' The default tau_prior is set adaptively:
#' - If p ≤ 5: tau_prior = 1
#' - If p > 5: tau_prior = 5/(p-5) * 2/√n
#'
#' where p is the number of predictors and n is the sample size.
#'
#' @references
#' Carvalho, C. M., Polson, N. G., & Scott, J. G. (2010). The horseshoe estimator
#' for sparse signals. Biometrika, 97(2), 465-480.
#'
#' Polson, N. G., Scott, J. G., & Windle, J. (2013). Bayesian inference for
#' logistic models using Pólya–Gamma latent variables. Journal of the American
#' Statistical Association, 108(504), 1339-1349.
#'
#' Piironen, J., & Vehtari, A. (2017). Sparsity information and regularization
#' in the horseshoe and other shrinkage priors. Electronic Journal of Statistics,
#' 11(2), 5018-5051.
#'
#' @examples
#' \dontrun{
#' # Generate synthetic data
#' n <- 100
#' p <- 20
#' X <- matrix(rnorm(n * p), n, p)
#' theta_true <- c(1, 0.5, -0.8, rep(0, p-2))  # Sparse coefficients
#' prob <- plogis(X %*% theta_true[-1] + theta_true[1])
#' Y <- rbinom(n, 1, prob)
#'
#' # Fit Bayesian logistic regression with horseshoe prior
#' posterior_samples <- HS.Logit(Y = Y, X = X, mc = 2000)
#'
#' # Posterior means
#' theta_hat <- colMeans(posterior_samples)
#' print(theta_hat)
#'
#' # Credible intervals
#' apply(posterior_samples, 2, quantile, c(0.025, 0.975))
#'
#' # Prediction for new data
#' X_new <- matrix(rnorm(10 * p), 10, p)
#' Z_new <- cbind(1, X_new)
#' prob_samples <- plogis(tcrossprod(posterior_samples, Z_new))
#' prob_mean <- colMeans(prob_samples)
#' 
#' # Compare with B.Logit (without horseshoe prior)
#' posterior_standard <- B.Logit(Y = Y, X = X, mc = 2000)
#' cbind(
#'   Horseshoe = colMeans(posterior_samples),
#'   Standard = colMeans(posterior_standard)
#' )
#' }
#'
#' @seealso \code{\link{DRBayes.PC}}, \code{\link{HS.Probit}}, \code{\link{B.Logit}}
#'
#' @importFrom extraDistr rinvgamma
#' @importFrom pgdraw pgdraw
#' @importFrom stats rnorm
#'
#' @export
HS.Logit <- function(Y, X, mc = 5000, beta0_prior = 1/100, tau_prior = NULL) {

  # Input validation
  if (!is.numeric(Y)) {
    stop("Y must be a numeric vector")
  }

  if (!all(Y %in% c(0, 1))) {
    stop("Y must be binary (0 or 1)")
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

  # Initialize parameters
  theta   <- rep(0, pp+1)
  tau2    <- 1
  lambda2 <- rep(1, pp)

  # Transform response for Pólya-Gamma augmentation
  iota <- Y - 0.5

  # Set adaptive tau prior if not specified
  if (is.null(tau_prior)) {
    tau_prior <- ifelse(pp <= 5, 1, 5/(pp-5) * 2/sqrt(nn))
  }

  # Storage for posterior samples
  post.theta <- matrix(NA, mc, pp+1)
  colnames(post.theta) <- c("(Intercept)", paste0("X", 1:pp))

  # MCMC sampling
  for (iter in 1:mc) {
    # Step 1: Sample Pólya-Gamma auxiliary variables
    b_w <- drop( tcrossprod(theta, Z) )
    w   <- pgdraw::pgdraw(1, b_w)

    # Step 2: Sample regression coefficients given latent variables and horseshoe prior
    # Precision matrix with horseshoe prior
    invV <- diag( c(beta0_prior, 1/(tau2*lambda2)) )

    b_theta <- drop( crossprod(Z, iota) )
    Q_theta <- crossprod(Z, w*Z) + invV
    U_theta <- chol(Q_theta)
    theta   <- drop( base::backsolve(U_theta,
                                     base::backsolve(U_theta, b_theta, transpose=TRUE) + stats::rnorm(pp+1)) )

    # Step 3: Sample horseshoe prior parameters
    # Sample auxiliary variables nu
    b_nu <- 1 + 1/lambda2
    nu   <- extraDistr::rinvgamma(pp, 1, b_nu)

    # Sample local shrinkage parameters lambda^2
    b_lambda <- 1/nu + theta[-1]^2/(2*tau2)
    lambda2  <- extraDistr::rinvgamma(pp, 1, b_lambda)

    # Sample auxiliary variable xi
    b_xi <- 1/tau_prior^2 + 1/tau2
    xi   <- extraDistr::rinvgamma(1, 1, b_xi)

    # Sample global shrinkage parameter tau^2
    a_tau <- (pp+1)/2
    b_tau <- 1/xi + 1/2*sum(theta[-1]^2/lambda2)
    tau2  <- extraDistr::rinvgamma(1, a_tau, b_tau)

    # Store samples
    post.theta[iter, ] <- theta
  }

  return(post.theta)
}

# Example usage (not run during package build)
if (FALSE) {
  Y <- data$AA
  X <- data[,-(1:2)]

  HS1.Logit <- HS.Logit(Y, X, mc = 100000, tau_prior = 1)
  HS2.Logit <- HS.Logit(Y, X, mc = 100000, tau_prior = NULL)
}
