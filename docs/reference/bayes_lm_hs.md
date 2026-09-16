# Bayesian Linear Regression with Horseshoe Prior

Fits a Bayesian linear regression model using the horseshoe prior for
regularization. This function is designed to work as an outcome model in
the
[`drbayes_pc`](https://t-momozaki.github.io/DRBayes/reference/drbayes_pc.md)
framework for causal inference.

## Usage

``` r
bayes_lm_hs(
  Y,
  X,
  mc = 5000,
  chains = 4L,
  init = NULL,
  unshrunk = integer(0),
  beta0_prior = NULL,
  unshrunk_prior = NULL,
  sigma_prior = c(1, 1),
  tau_prior = NULL,
  p0 = 5
)
```

## Arguments

- Y:

  Numeric vector of response variables (continuous outcomes).

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
  sensitivity to where they were started: the ordinary least squares
  estimate is computed with
  [`lm.fit`](https://rdrr.io/r/stats/lmfit.html), chain one starts
  there, and each further chain starts at a normal draw centred on it
  with a standard deviation of four times that estimate's own standard
  error. If the least squares fit is unusable, for instance because the
  design is rank deficient, the starting values are drawn from a normal
  on the scale of the prior; the horseshoe itself is not used for this,
  since its half-Cauchy tails would place chains at values no amount of
  sampling would recover from. The scan below draws the coefficients
  first, from a full conditional that does not depend on their previous
  value, so a chain is started by drawing every variance component from
  its own full conditional given these coefficients: a starting value
  far from the data begins its chain at a large residual sum of squares,
  and so at a large error variance. Chains started in different places
  therefore stay different for as long as it takes them to converge,
  which is what R-hat needs in order to say anything at all (Vehtari et
  al., 2021).

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

  Numeric. Prior precision for the intercept term, on the absolute scale
  of the data. NULL (default) asks for a precision scaled to the data,
  described under "Details".

- unshrunk_prior:

  Numeric. Prior precision for the coefficients named in `unshrunk`, on
  the absolute scale of the data. NULL (default) asks for a precision
  scaled to the data, described under "Details".

- sigma_prior:

  Numeric vector of length 2. Inverse gamma prior parameters (shape,
  scale) for the error variance. Default is c(1, 1).

- tau_prior:

  Numeric. Scale of the half-Cauchy prior on the global shrinkage
  parameter tau, which is itself free of the units of Y because the
  horseshoe coefficient scales carry a factor of sigma. If NULL,
  defaults to an adaptive value based on problem dimension.

- p0:

  Integer. Guess at how many coefficients are non-zero, used only when
  `tau_prior` is NULL. Default is 5. See `tau_prior`.

## Value

An array of dimension `mc` by `chains` by `ncol(X)+1` containing
posterior samples for the regression coefficients, indexed by iteration,
then chain, then parameter. This is the layout the posterior package
calls a draws_array. The third dimension is named with "(Intercept)"
followed by the column names of X, or "X1", "X2", etc. when X has none.
The first two dimensions are unnamed.

## Details

The horseshoe prior is a continuous shrinkage prior that provides strong
regularization for small coefficients while allowing large coefficients
to remain relatively unshrunk. In the context of causal inference, this
function applies different priors to different types of coefficients:

\$\$Y \| \beta, \sigma^2 \sim N(X\beta, \sigma^2 I)\$\$ \$\$\beta_0 \sim
N(0, 1/c_0)\$\$ \$\$\beta_1 \sim N(0, 1/c_1)\$\$ \$\$\beta_j \|
\lambda_j, \tau, \sigma^2 \sim N(0, \sigma^2 \tau^2 \lambda_j^2), \quad
j = 2, \ldots, p\$\$ \$\$\lambda_j \sim C^{+}(0, 1), \quad j = 2,
\ldots, p\$\$ \$\$\tau \sim C^{+}(0, t)\$\$ \$\$\sigma^2 \sim
\mathrm{InvGamma}(a, b)\$\$

Here \\\beta_1\\ stands for the coefficients named in `unshrunk`, which
receive no shrinkage, and \\\beta_2, \ldots, \beta_p\\ for the rest.
\\C^{+}(0, s)\\ denotes the half-Cauchy distribution with scale \\s\\.
The prior precisions \\c_0\\ and \\c_1\\ are set by `beta0_prior` and
`unshrunk_prior`, the global shrinkage scale \\t\\ is set by
`tau_prior`, and the shape and scale \\(a, b)\\ are set by
`sigma_prior`.

The horseshoe coefficient scales carry a factor of \\\sigma\\, as in
Carvalho, Polson and Scott (2010) and in the sampler of Makalic and
Schmidt (2016). That is what makes the amount of shrinkage independent
of the units of Y without any reference to the data: the posterior
precision of those coefficients is \\(X'X + D)/\sigma^2\\ with \\D\\ the
prior precision, so \\D\\ is compared with \\X'X\\ and not with
\\X'X/\sigma^2\\. It also leaves \\\tau\\ free of those units, which is
what the default `tau_prior` below assumes. Those coefficients, and only
those, add degrees of freedom to the full conditional of \\\sigma^2\\,
which is \\\mathrm{InvGamma}(a + (n + m)/2, b + (\mathrm{RSS} + Q)/2)\\
for \\m\\ shrunk coefficients, residual sum of squares \\\mathrm{RSS}\\
and prior quadratic form \\Q = \sum_j \beta_j^2/(\tau^2\lambda_j^2)\\.

The intercept and the unshrunk coefficients are treated differently,
because they have no shrinkage parameter to adapt to their size. Their
priors are fixed variances rather than multiples of \\\sigma^2\\, so
that they cannot inflate the error variance: on the simulation of
Orihara et al. (2025) the intercept is about 100 and the treatment
coefficient about 110 while \\\sigma\\ is 1, and a prior variance of
\\100\sigma^2\\ on each would add 221 to a residual sum of squares of
499. In exchange their precisions have to be set on the scale of the
data, and the defaults are the ones
[`bayes_lm`](https://t-momozaki.github.io/DRBayes/reference/bayes_lm.md)
documents under "Prior Specification": a prior standard deviation of
\\100(\|\mathrm{mean}(Y)\| + \mathrm{sd}(Y))\\ for the intercept, which
is effectively flat, and of \\2.5 \mathrm{sd}(Y)/s_j\\ for an unshrunk
coefficient on a column with standard deviation \\s_j\\. Those scales
are read for every column and not only for the exempt ones, since they
also place the starting values and the reference fit the residual sum of
squares is expanded around, so data too far apart in scale for the
defaults to be represented stop the fit whichever precisions are named
here.
[`bayes_lm`](https://t-momozaki.github.io/DRBayes/reference/bayes_lm.md)
documents where that is.

The default value of `tau_prior` follows Piironen and Vehtari (2017):
with \\p\\ coefficients under the horseshoe, \\n\\ observations and
\\p_0\\ of them guessed to be non-zero, it is \\p_0/(p - p_0) \times
1/\sqrt{n}\\, or 1 when \\p \le p_0\\ leaves nothing to shrink towards.

## References

Carvalho, C. M., Polson, N. G., & Scott, J. G. (2010). The horseshoe
estimator for sparse signals. Biometrika, 97(2), 465-480.

Makalic, E., & Schmidt, D. F. (2016). A simple sampler for the horseshoe
estimator. IEEE Signal Processing Letters, 23(1), 179-182.

Orihara, S., Momozaki, T., & Sugasawa, S. (2025). Bayesian Doubly Robust
Causal Inference via Posterior Coupling. arXiv preprint
arXiv:2506.04868.

Piironen, J., & Vehtari, A. (2017). Sparsity information and
regularization in the horseshoe and other shrinkage priors. Electronic
Journal of Statistics, 11(2), 5018-5051.

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Burkner, P.-C.
(2021). Rank-normalization, folding, and localization: an improved R-hat
for assessing convergence of MCMC. Bayesian Analysis, 16(2), 667-718.

## See also

[`drbayes_pc`](https://t-momozaki.github.io/DRBayes/reference/drbayes_pc.md),
[`bayes_logit_hs`](https://t-momozaki.github.io/DRBayes/reference/bayes_logit_hs.md)

## Examples

``` r
# Generate synthetic data for causal inference
set.seed(1)
n <- 100
p <- 20
W <- matrix(rnorm(n * (p - 1)), n, p - 1,
            dimnames = list(NULL, paste0("W", seq_len(p - 1))))
A <- rbinom(n, 1, 0.5)  # Treatment
X <- cbind(A = A, W)  # Treatment is the first column

# Sparse confounding: only the first three confounders matter
beta <- c(1.5, 1, -1, 0.5, rep(0, p - 4))  # A, then the confounders
Y <- 2 + drop(X %*% beta) + rnorm(n, 0, 0.5)

# Fit Bayesian linear regression with horseshoe prior. unshrunk = 1
# exempts the treatment, the first column of X, from the horseshoe:
# shrinking it would attenuate the effect being estimated.
draws <- bayes_lm_hs(Y = Y, X = X, mc = 1000, chains = 4, unshrunk = 1)
dim(draws)  # iterations, chains, parameters
#> [1] 1000    4   21

# Posterior means, pooling the chains: the treatment effect and the three
# real confounders are recovered, and the sixteen null coefficients are
# shrunk to a few hundredths.
beta_hat <- apply(draws, 3, mean)
round(beta_hat[1:5], 2)
#> (Intercept)           A          W1          W2          W3 
#>        1.97        1.57        1.01       -0.99        0.52 
range(beta_hat[6:(p + 1)])
#> [1] -0.03989652  0.02808200

# Credible interval for the treatment effect
quantile(draws[, , "A"], c(0.025, 0.975))
#>     2.5%    97.5% 
#> 1.336625 1.804369 

# Pool the chains into an (iterations * chains) by parameters matrix
posterior_samples <- apply(draws, 3, as.vector)
```
