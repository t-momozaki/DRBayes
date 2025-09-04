#' Bayesian Linear Regression with Horseshoe Prior
#'
#' Fits a Bayesian linear regression model using the horseshoe prior for
#' regularization. This function is designed to work as an outcome model
#' in the \code{\link{DRBayes.PC}} framework for causal inference.
#'
#' @param Y Numeric vector of response variables (continuous outcomes).
#' @param X Matrix or data.frame of predictor variables. The intercept is
#'   automatically added, so do not include a column of ones.
#' @param mc Integer. Number of MCMC iterations. Default is 5000.
#' @param beta0_prior Numeric. Prior precision for the intercept term. Default is 1/100.
#' @param treatment_prior Numeric. Prior precision for the treatment coefficient (first column of X).
#'   This allows the treatment effect to avoid shrinkage. Default is 1/100.
#' @param sigma_prior Numeric vector of length 2. Inverse gamma prior parameters
#'   (shape, rate) for the error variance. Default is c(1, 1).
#' @param tau_prior Numeric. Prior parameter for the global shrinkage parameter tau.
#'   If NULL, defaults to an adaptive value based on problem dimension.
#'
#' @return Matrix of dimension \code{mc} by \code{ncol(X)+1} containing posterior
#'   samples for the regression coefficients. The first column corresponds to
#'   the intercept, and the remaining columns to the coefficients for X.
#'
#' @details
#' The horseshoe prior is a continuous shrinkage prior that provides strong
#' regularization for small coefficients while allowing large coefficients
#' to remain relatively unshrunk. In the context of causal inference, this
#' function applies different priors to different types of coefficients:
#'
#' \deqn{Y | \beta, \sigma^2 \sim N(X\beta, \sigma^2 I)}
#' \deqn{\beta_0 \sim N(0, 1/\text{beta0_prior})}
#' \deqn{\beta_1 \sim N(0, 1/\text{treatment_prior})} (treatment coefficient, no shrinkage)
#' \deqn{\beta_j | \lambda_j, \tau \sim N(0, \tau^2 \lambda_j^2)} for j = 2, ..., p (other covariates)
#' \deqn{\lambda_j \sim \text{Half-Cauchy}(0, 1)} for j = 2, ..., p
#' \deqn{\tau \sim \text{Half-Cauchy}(0, \text{tau_prior})}
#' \deqn{\sigma^2 \sim \text{InvGamma}(\text{sigma_prior})}
#'
#' The default tau_prior is set adaptively:
#' - If p ≤ 5: tau_prior = 1
#' - If p > 5: tau_prior = 5/(p-5) * 2/√n
#'
#' where p is the number of predictors and n is the sample size.
#' 
#' **Important for Causal Inference**: The first column of X is assumed to be 
#' the treatment variable, and it receives a standard normal prior without 
#' horseshoe shrinkage. Only the remaining covariates (confounders) receive 
#' horseshoe regularization.
#'
#' @references
#' Carvalho, C. M., Polson, N. G., & Scott, J. G. (2010). The horseshoe estimator
#' for sparse signals. Biometrika, 97(2), 465-480.
#'
#' Piironen, J., & Vehtari, A. (2017). Sparsity information and regularization
#' in the horseshoe and other shrinkage priors. Electronic Journal of Statistics,
#' 11(2), 5018-5051.
#'
#' @examples
#' \dontrun{
#' # Generate synthetic data for causal inference
#' n <- 100
#' p <- 20
#' X_confounders <- matrix(rnorm(n * (p-1)), n, p-1)  # Confounders
#' A <- rbinom(n, 1, 0.5)  # Treatment
#' X <- cbind(A, X_confounders)  # Treatment is first column
#' 
#' # Sparse confounders: only first few matter
#' beta_true <- c(2, 1.5, 1, -1, 0.5, rep(0, p-4))  # intercept, treatment, confounders
#' Y <- X %*% beta_true[-1] + beta_true[1] + rnorm(n, 0, 0.5)
#'
#' # Fit Bayesian linear regression with horseshoe prior
#' # Treatment effect (column 1) won't be shrunk, confounders will be
#' posterior_samples <- HS.LM(Y = Y, X = X, mc = 2000)
#'
#' # Posterior means
#' beta_hat <- colMeans(posterior_samples)
#' print(beta_hat)
#'
#' # Credible intervals
#' apply(posterior_samples, 2, quantile, c(0.025, 0.975))
#' }
#'
#' @seealso \code{\link{DRBayes.PC}}, \code{\link{HS.Logit}}
#'
#' @importFrom extraDistr rinvgamma
#' @importFrom stats rnorm
#'
#' @export
HS.LM <- function(Y, X, mc = 5000,
                  beta0_prior = 1/100, treatment_prior = 1/100, sigma_prior = c(1, 1), tau_prior = NULL) {

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

  if (pp < 2) {
    stop("X must have at least 2 columns (treatment + at least one confounder)")
  }

  # Prior parameters
  a <- sigma_prior[1]
  b <- sigma_prior[2]

  # Initialize parameters
  sigma2  <- 1
  tau2    <- 1
  # Only confounders (columns 2 onwards in X, which are columns 3 onwards in Z) get horseshoe
  lambda2 <- rep(1, pp - 1)  # pp-1 confounders (excluding treatment)

  # Set adaptive tau prior if not specified
  if (is.null(tau_prior)) {
    # Based on number of confounders, not total predictors
    p_confounders <- pp - 1
    tau_prior <- ifelse(p_confounders <= 5, 1, 5/(p_confounders-5) * 2/sqrt(nn))
  }

  # Storage for posterior samples
  post.theta <- matrix(NA, mc, pp+1)
  colnames(post.theta) <- c("(Intercept)", paste0("X", 1:pp))

  # MCMC sampling
  for (iter in 1:mc) {
    # Sample regression coefficients with different priors:
    # - Intercept: beta0_prior
    # - Treatment (X1): treatment_prior (no horseshoe)
    # - Confounders (X2, X3, ...): horseshoe prior

    invV <- diag(c(beta0_prior,                    # intercept
                   treatment_prior,                # treatment coefficient (no shrinkage)
                   1/(tau2*lambda2)))              # confounder coefficients (with shrinkage)

    b_theta <- 1/sigma2 * drop( crossprod(Z, Y) )
    Q_theta <- 1/sigma2 * crossprod(Z) + invV
    U_theta <- chol(Q_theta)
    theta   <- drop( base::backsolve(U_theta,
                                     base::backsolve(U_theta, b_theta, transpose=TRUE) + stats::rnorm(pp+1)) )

    # Sample error variance
    a_sigma <- a + nn/2
    b_sigma <- b + crossprod(Y-drop(tcrossprod(theta,Z))) / 2
    sigma2  <- extraDistr::rinvgamma(1, a_sigma, b_sigma)

    # Sample auxiliary variables nu for confounders only
    b_nu <- 1 + 1/lambda2
    nu   <- extraDistr::rinvgamma(pp-1, 1, b_nu)  # pp-1 confounders

    # Sample local shrinkage parameters lambda^2 for confounders only
    # theta[3:(pp+1)] are the confounder coefficients
    b_lambda <- 1/nu + theta[3:(pp+1)]^2/(2*tau2)
    lambda2  <- extraDistr::rinvgamma(pp-1, 1, b_lambda)

    # Sample auxiliary variable xi
    b_xi <- 1/tau_prior^2 + 1/tau2
    xi   <- extraDistr::rinvgamma(1, 1, b_xi)

    # Sample global shrinkage parameter tau^2 based on confounders only
    a_tau <- (pp-1+1)/2  # (number of confounders + 1)/2
    b_tau <- 1/xi + 1/2*sum(theta[3:(pp+1)]^2/lambda2)  # sum over confounders
    tau2  <- extraDistr::rinvgamma(1, a_tau, b_tau)

    # Store samples
    post.theta[iter, ] <- theta
  }

  return(post.theta)
}

# Example usage (not run during package build)
if (FALSE) {
  # Example showing treatment effect is not shrunk
  set.seed(123)
  n <- 200
  p_confounders <- 15
  
  # Generate data
  A <- rbinom(n, 1, 0.5)  # Treatment
  X_conf <- matrix(rnorm(n * p_confounders), n, p_confounders)  # Confounders
  X <- cbind(A, X_conf)  # Treatment first, then confounders
  
  # True coefficients: large treatment effect, sparse confounders
  beta_treatment <- 2.0  # Large treatment effect
  beta_conf <- c(1, -0.8, 0.5, rep(0, p_confounders-3))  # Sparse confounders
  
  Y <- 1 + beta_treatment * A + X_conf %*% beta_conf + rnorm(n, 0, 0.8)
  
  # Fit model
  result_hs <- HS.LM(Y, X, mc = 5000, bn = 1000)
  
  # Check that treatment effect is not shrunk towards zero
  cat("True treatment effect:", beta_treatment, "\n")
  cat("Estimated treatment effect:", round(mean(result_hs[,2]), 3), "\n")
  cat("Treatment effect 95% CI:", round(quantile(result_hs[,2], c(0.025, 0.975)), 3), "\n")
  
  # Check that irrelevant confounders are shrunk
  cat("True confounder effects (should be shrunk if zero):\n")
  print(beta_conf)
  cat("Estimated confounder effects:\n")
  print(round(colMeans(result_hs[,4:ncol(result_hs)]), 3))
}