# Bayesian Logistic Regression with Polya-Gamma Latent Variables

Performs Bayesian logistic regression using the Polya-Gamma data
augmentation scheme proposed by Polson, Scott, and Windle (2013). This
method provides an exact and efficient MCMC algorithm for binary
logistic regression by introducing Polya-Gamma auxiliary variables that
render the likelihood conditionally Gaussian.

## Usage

``` r
bayes_logit(Y, X, mc = 5000, chains = 4L, init = NULL, theta_prior = NULL)
```

## Arguments

- Y:

  A numeric vector of binary outcomes (0 or 1). Must contain only values
  0 and 1.

- X:

  A matrix or data.frame of covariates. Each row corresponds to an
  observation and each column to a predictor variable. An intercept term
  will be automatically added.

- mc:

  A positive integer specifying the number of MCMC iterations (default:
  5000).

- chains:

  A positive integer, the number of Markov chains to run (default: 4).
  R-hat compares chains, so at least four are recommended.

- init:

  Starting values for the coefficient vector, either NULL (default) or a
  list of length `chains` whose elements are numeric vectors of length
  `ncol(X) + 1`, the intercept followed by the coefficients of the
  columns of X. When NULL the starting values are generated from the
  data and deliberately overdispersed, so that the chains can reveal
  sensitivity to where they were started: the maximum likelihood
  estimate is computed with
  [`glm.fit`](https://rdrr.io/r/stats/glm.html) under the binomial
  family with a logit link, chain one starts there, and each further
  chain starts at a normal draw centred on it with a standard deviation
  of four times that estimate's own standard error. If the fit is
  unusable, for instance under separation or a rank deficient design,
  the starting values are drawn from the prior instead.

- theta_prior:

  Either a scalar or a square matrix specifying the precision for the
  normal prior on regression coefficients. If scalar, assumes diagonal
  precision matrix with all elements equal to this value. If NULL
  (default), a weakly informative prior with precision 1/100 \* I is
  used, where I is the identity matrix.

## Value

An array of posterior samples with dimensions `mc` x `chains` x
(`ncol(X)` + 1), indexed by iteration, then chain, then parameter. This
is the layout the posterior package calls a draws_array. The third
dimension is named with "(Intercept)" for the intercept term and with
the column names of X, or "X1", "X2", etc. when X has none. The first
two dimensions are unnamed.

## Details

Bayesian Logistic Regression using Polya-Gamma Data Augmentation

The function implements the Polya-Gamma data augmentation strategy which
represents the logistic likelihood as a scale mixture of Gaussians. The
algorithm alternates between:

1.  Sampling Polya-Gamma auxiliary variables w_i ~ PG(1, X_i^T \* theta)

2.  Sampling regression coefficients theta from a multivariate normal
    distribution

The method is exact (no approximation) and often more efficient than
alternative data augmentation schemes, especially for hierarchical
models.

The prior for regression coefficients is theta ~ N(0, V), where V^(-1)
is specified by `theta_prior`. The default prior is relatively
uninformative with precision 0.01 for each coefficient.

## Note

The function requires the `pgdraw` package for sampling Polya-Gamma
random variables.

## References

Polson, N. G., Scott, J. G., & Windle, J. (2013). Bayesian inference for
logistic models using Polya-Gamma latent variables. Journal of the
American Statistical Association, 108(504), 1339-1349.

## See also

[`glm`](https://rdrr.io/r/stats/glm.html) for classical logistic
regression, [`pgdraw`](https://rdrr.io/pkg/pgdraw/man/pgdraw.html) for
Polya-Gamma random number generation.

## Examples

``` r
# Simulate binary logistic regression data
set.seed(123)
n <- 200
X1 <- rnorm(n)
X2 <- rnorm(n)
X <- cbind(X1, X2)

# True coefficients: intercept = 0.5, X1 = 1.2, X2 = -0.8
eta <- 0.5 + 1.2 * X1 - 0.8 * X2
prob <- plogis(eta)  # logistic transformation
Y <- rbinom(n, 1, prob)

# Fit Bayesian logistic regression
draws <- bayes_logit(Y, X, mc = 1000, chains = 4)
dim(draws)  # iterations, chains, parameters
#> [1] 1000    4    3

# Summarize posterior, pooling the chains
apply(draws, 3, mean)  # posterior means
#> (Intercept)          X1          X2 
#>   0.3464994   1.2019891  -1.0427484 
apply(draws, 3, sd)  # posterior standard deviations
#> (Intercept)          X1          X2 
#>   0.1751417   0.2192001   0.2045132 

# Pool the chains into an (iterations * chains) by parameters matrix
posterior_samples <- apply(draws, 3, as.vector)
```
