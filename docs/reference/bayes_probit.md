# Bayesian Probit Regression with Latent Variables

Performs Bayesian probit regression using the latent variable data
augmentation scheme proposed by Albert and Chib (1993). This method
provides an exact and efficient MCMC algorithm for binary probit
regression by introducing latent Gaussian variables that render the
likelihood conditionally Gaussian.

## Usage

``` r
bayes_probit(Y, X, mc = 5000, chains = 4L, init = NULL, theta_prior = NULL)
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
  family with a probit link, chain one starts there, and each further
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

Bayesian Probit Regression using Latent Variable Approach

The function implements the Albert and Chib (1993) latent variable
approach which introduces continuous latent variables W_i such that:

- Y_i = 1 if W_i \> 0, Y_i = 0 if W_i \<= 0

- W_i ~ N(X_i^T \* theta, 1)

The algorithm alternates between:

1.  Sampling latent variables W_i from truncated normal distributions

2.  Sampling regression coefficients theta from a multivariate normal
    distribution

For observations with Y_i = 1, W_i is sampled from N(X_i^T \* theta, 1)
truncated to be positive. For Y_i = 0, W_i is sampled from the same
distribution truncated to be negative or zero.

The prior for regression coefficients is theta ~ N(0, V), where V^(-1)
is specified by `theta_prior`. The default prior is relatively
uninformative with precision 0.01 for each coefficient.

## Note

The function requires the `RcppTN` package for sampling truncated normal
random variables.

## References

Albert, J. H., & Chib, S. (1993). Bayesian analysis of binary and
polychotomous response data. Journal of the American Statistical
Association, 88(422), 669-679.

## See also

[`glm`](https://rdrr.io/r/stats/glm.html) for classical probit
regression,
[`bayes_logit`](https://t-momozaki.github.io/DRBayes/reference/bayes_logit.md)
for Bayesian logistic regression,
[`rtn`](https://rdrr.io/pkg/RcppTN/man/rtn.html) for truncated normal
random number generation.

## Examples

``` r
# Simulate binary probit regression data
set.seed(123)
n <- 200
X1 <- rnorm(n)
X2 <- rnorm(n)
X <- cbind(X1, X2)

# True coefficients: intercept = 0.5, X1 = 1.2, X2 = -0.8
eta <- 0.5 + 1.2 * X1 - 0.8 * X2
prob <- pnorm(eta)  # probit transformation
Y <- rbinom(n, 1, prob)

# Fit Bayesian probit regression
draws <- bayes_probit(Y, X, mc = 1000, chains = 4)
dim(draws)  # iterations, chains, parameters
#> [1] 1000    4    3

# Summarize posterior, pooling the chains
apply(draws, 3, mean)  # posterior means
#> (Intercept)          X1          X2 
#>   0.4066154   1.3067997  -1.0245292 
apply(draws, 3, sd)  # posterior standard deviations
#> (Intercept)          X1          X2 
#>   0.1206760   0.2072669   0.1624851 

# Compare with classical probit regression
classical_fit <- glm(Y ~ X, family = binomial(link = "probit"))
cbind(
  Bayesian = apply(draws, 3, mean),
  Classical = coef(classical_fit)
)
#>               Bayesian  Classical
#> (Intercept)  0.4066154  0.4003846
#> X1           1.3067997  1.2803651
#> X2          -1.0245292 -1.0098641
```
