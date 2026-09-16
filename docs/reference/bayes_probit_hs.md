# Bayesian Probit Regression with Horseshoe Prior

Fits a Bayesian probit regression model using the horseshoe prior for
regularization. This function is designed to work as an outcome model or
propensity score model in the
[`drbayes_pc`](https://t-momozaki.github.io/DRBayes/reference/drbayes_pc.md)
framework for causal inference.

## Usage

``` r
bayes_probit_hs(
  Y,
  X,
  mc = 5000,
  chains = 4L,
  init = NULL,
  unshrunk = integer(0),
  beta0_prior = 1/100,
  unshrunk_prior = 1/100,
  tau_prior = NULL,
  p0 = 5
)
```

## Arguments

- Y:

  Numeric vector of binary response variables (0 or 1).

- X:

  Matrix or data.frame of predictor variables. The intercept is
  automatically added, so do not include a column of ones.

- mc:

  Integer. Number of MCMC iterations. Default is 5000.

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
  the starting values are drawn from a normal on the scale of the prior;
  the horseshoe itself is not used for this, since its half-Cauchy tails
  would place chains at values no amount of sampling would recover from.

- unshrunk:

  Integer column indices into `X` whose coefficients are exempt from the
  horseshoe prior and instead receive the normal prior set by
  `unshrunk_prior`. The default `integer(0)` shrinks every coefficient,
  which is the plain horseshoe and the right choice for a propensity
  score model. For an outcome model in a causal analysis, exempt the
  treatment and every interaction with it: shrinking those attenuates
  the treatment effect.
  [`drbayes_pc`](https://t-momozaki.github.io/DRBayes/reference/drbayes_pc.md)
  works these columns out from the formula and passes them for you.

- beta0_prior:

  Numeric. Prior precision for the intercept term. Default is 1/100.

- unshrunk_prior:

  Numeric. Prior precision for the coefficients named in `unshrunk`.
  Default is 1/100, i.e. a normal prior with variance 100.

- tau_prior:

  Numeric. Prior parameter for the global shrinkage parameter tau. If
  NULL, defaults to an adaptive value based on problem dimension.

- p0:

  Integer. Guess at how many coefficients are non-zero, used only when
  `tau_prior` is NULL. Default is 5. See `tau_prior`.

## Value

An array of dimension `mc` by `chains` by `ncol(X)+1` containing
posterior samples for the probit regression coefficients, indexed by
iteration, then chain, then parameter. This is the layout the posterior
package calls a draws_array. The third dimension is named with
"(Intercept)" followed by the column names of X, or "X1", "X2", etc.
when X has none. The first two dimensions are unnamed.

## Details

The horseshoe prior provides strong regularization for small
coefficients while allowing large coefficients to remain relatively
unshrunk. The model uses the Albert and Chib (1993) latent variable
approach for efficient sampling:

\$\$Y_i \| \theta \sim \mathrm{Bernoulli}(\Phi(\beta_0 + X_i \beta))\$\$
\$\$W_i \| \theta \sim N(\beta_0 + X_i \beta, 1) \mbox{ with } Y_i =
\mathbf{1}(W_i \> 0)\$\$ \$\$\beta_0 \sim N(0, 1/c_0)\$\$ \$\$\beta_j \|
\lambda_j, \tau \sim N(0, \tau^2 \lambda_j^2)\$\$ \$\$\lambda_j \sim
C^{+}(0, 1)\$\$ \$\$\tau \sim C^{+}(0, s\_\tau)\$\$

Here \\c_0\\ is the prior precision for the intercept set by
`beta0_prior`, \\s\_\tau\\ is the global shrinkage scale set by
`tau_prior`, and \\C^{+}(0, s)\\ denotes the half-Cauchy distribution
with scale \\s\\.

The data augmentation introduces latent variables. For \\Y_i = 1\\ the
latent variable is truncated to \\(0, \infty)\\: \$\$W_i \| \theta \sim
TN(X_i \theta, 1, 0, \infty)\$\$ and for \\Y_i = 0\\ it is truncated to
\\(-\infty, 0\]\\: \$\$W_i \| \theta \sim TN(X_i \theta, 1, -\infty,
0\]\$\$

where TN denotes the truncated normal distribution and \\\Phi\\ is the
standard normal cumulative distribution function.

The default `tau_prior` is set adaptively:

- If p \<= 5: tau_prior = 1

- If p \> 5: tau_prior = 5/(p-5) \* 1/sqrt(n)

where p is the number of coefficients under the horseshoe, that is
`ncol(X)` less the length of `unshrunk`, and n is the sample size. The
scale factor is 1 because the Albert-Chib latent variable has unit
variance; the corresponding factor in
[`bayes_logit_hs()`](https://t-momozaki.github.io/DRBayes/reference/bayes_logit_hs.md)
is 2, which is the logistic one, and the two deliberately differ.

## References

Albert, J. H., & Chib, S. (1993). Bayesian analysis of binary and
polychotomous response data. Journal of the American Statistical
Association, 88(422), 669-679.

Carvalho, C. M., Polson, N. G., & Scott, J. G. (2010). The horseshoe
estimator for sparse signals. Biometrika, 97(2), 465-480.

Piironen, J., & Vehtari, A. (2017). Sparsity information and
regularization in the horseshoe and other shrinkage priors. Electronic
Journal of Statistics, 11(2), 5018-5051.

## See also

[`drbayes_pc`](https://t-momozaki.github.io/DRBayes/reference/drbayes_pc.md),
[`bayes_logit_hs`](https://t-momozaki.github.io/DRBayes/reference/bayes_logit_hs.md),
[`bayes_probit`](https://t-momozaki.github.io/DRBayes/reference/bayes_probit.md)

## Examples

``` r
# Generate synthetic data: two covariates matter, eighteen do not
set.seed(1)
n <- 300
p <- 20
X <- matrix(rnorm(n * p), n, p)
theta_true <- c(0.5, 1, -0.8, rep(0, p - 2))  # Sparse coefficients
prob <- pnorm(theta_true[1] + drop(X %*% theta_true[-1]))  # Probit link
Y <- rbinom(n, 1, prob)

# Fit Bayesian probit regression with horseshoe prior
draws <- bayes_probit_hs(Y = Y, X = X, mc = 800, chains = 4)
dim(draws)  # iterations, chains, parameters
#> [1] 800   4  21

# Posterior means, pooling the chains
theta_hat <- apply(draws, 3, mean)
round(theta_hat[1:3], 2)
#> (Intercept)          X1          X2 
#>        0.55        0.85       -0.83 

# Credible intervals for the intercept and the two real signals
apply(draws[, , 1:3], 3, quantile, c(0.025, 0.975))
#>       (Intercept)        X1        X2
#> 2.5%    0.3729065 0.6145877 -1.046227
#> 97.5%   0.7325077 1.0813367 -0.635326

# Prediction for new data, from the pooled draws
posterior_samples <- apply(draws, 3, as.vector)
X_new <- matrix(rnorm(10 * p), 10, p)
Z_new <- cbind(1, X_new)
prob_samples <- pnorm(tcrossprod(posterior_samples, Z_new))
prob_mean <- colMeans(prob_samples)

# Compare with bayes_probit, which shrinks nothing. Both find the two
# signals; the horseshoe holds the eighteen null coefficients several times
# closer to zero.
draws_standard <- bayes_probit(Y = Y, X = X, mc = 800, chains = 4)
round(cbind(
  Horseshoe = apply(draws, 3, mean),
  Standard  = apply(draws_standard, 3, mean)
)[1:3, ], 2)
#>             Horseshoe Standard
#> (Intercept)      0.55     0.58
#> X1               0.85     0.99
#> X2              -0.83    -0.93
range(theta_hat[4:(p + 1)])  # the null coefficients under the horseshoe
#> [1] -0.04090990  0.01614114
range(apply(draws_standard, 3, mean)[4:(p + 1)])  # and without it
#> [1] -0.2210924  0.1280368
```
