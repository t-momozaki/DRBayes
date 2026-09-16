# Bayesian Doubly Robust Causal Inference via Posterior Coupling

This function implements Bayesian doubly robust estimation for average
treatment effects (ATE) using posterior coupling. The method
incorporates propensity score information via moment conditions directly
into the posterior distribution, avoiding the feedback problem and
enabling a fully Bayesian interpretation of DR estimation without
requiring two-step estimation.

## Usage

``` r
drbayes_pc(
  outcome.formula,
  ps.formula,
  data,
  family = NULL,
  link = NULL,
  ps.link = "logit",
  outcome.prior = c("normal", "horseshoe"),
  ps.prior = c("normal", "horseshoe"),
  mc = 5000,
  bn = 1000,
  thin = 2,
  chains = 4L,
  method = c("smc", "is"),
  control = drbayes_control(),
  diagnostics = c("warn", "error", "none"),
  outcome.model = NULL,
  ps.model = NULL,
  outcome.samples = NULL,
  ps.samples = NULL,
  outcome.mu = NULL,
  outcome.mu1 = NULL,
  outcome.mu0 = NULL,
  outcome.priors = list(),
  ps.priors = list(),
  seed = NULL,
  verbose = TRUE,
  na.action = "na.omit"
)
```

## Arguments

- outcome.formula:

  Formula specifying the outcome model. For example: Y ~ A + X1 + X2 +
  A:X1. The response variable should be on the left side and predictors
  including treatment variable on the right side. Interactions and
  transformations are fully supported.

- ps.formula:

  Formula specifying the propensity score model. For example: A ~ X1 +
  X2 + X3. The treatment variable should be on the left side and
  confounders on the right side.

- data:

  Data.frame containing all variables specified in the formulas.

- family:

  Character string specifying the outcome model family. If NULL
  (default), automatically detects based on Y values (binary -\>
  "binomial", otherwise -\> "gaussian"). Options are:

  - NULL (default): Automatically detect based on Y values

  - "gaussian": Linear regression for continuous outcomes

  - "binomial": Logistic or probit regression for binary outcomes

- link:

  Character string specifying the link function of the OUTCOME model. If
  NULL (default), automatically determined from family: "identity" for
  gaussian, "logit" for binomial. Options are: "identity", "logit",
  "probit".

- ps.link:

  Character string specifying the link function of the PROPENSITY SCORE
  model, either "logit" (default) or "probit". This is the link the
  moment condition uses to turn the propensity score draws into
  probabilities, so it must match the scale on which `ps.model` or
  `ps.samples` produced them. Using `bayes_probit` or `bayes_probit_hs`
  without setting `ps.link = "probit"` evaluates probit coefficients
  with the logistic inverse link and silently solves a different
  constraint.

- outcome.prior:

  Character. Prior for the outcome model coefficients, either `"normal"`
  (default) or `"horseshoe"`. Given the family and the link, the
  sampling algorithm follows, so the prior is the only choice left to
  make and there is no need to hand in a sampler. Use `"horseshoe"` when
  there are many candidate confounders; the treatment effect and its
  interactions are exempted from shrinkage automatically.

- ps.prior:

  Character. As `outcome.prior`, for the propensity score model. The
  default `"normal"` is the right choice unless the model is high
  dimensional.

- mc:

  Integer. Number of MCMC iterations for internal sampling. Ignored when
  using external samples. Default is 5000.

- bn:

  Integer. Number of burn-in iterations to discard. Applied to both
  internal and external samples. Default is 1000.

- thin:

  Integer. Thinning interval for MCMC samples. Applied to both internal
  and external samples. Default is 2.

- chains:

  Integer. Number of MCMC chains, started from dispersed values. Default
  4, the minimum Vehtari et al. (2021) recommend for the convergence
  diagnostics.

- method:

  Character. How the tilting parameter is found. `"smc"` is Algorithm 2
  of the paper, which walks a grid of tilting parameters and rejuvenates
  the particles at each step. `"is"` is Algorithm 1, a single importance
  sampling step, which is faster but degenerates when the two posteriors
  put little mass where the moment condition holds. Default `"smc"`.

- control:

  List of tuning parameters from
  [`drbayes_control`](https://t-momozaki.github.io/DRBayes/reference/drbayes_control.md),
  covering the sequential Monte Carlo sweep and the diagnostic
  thresholds. `mc`, `bn`, `thin` and `chains` also live there; naming
  one of those directly overrides the control object, which reads better
  for the settings people change most often.

- diagnostics:

  Character. What to do when the posterior draws fail the convergence
  criteria: `"warn"` (default), `"error"`, or `"none"` to skip the
  check. Both posteriors are assessed before the sequential Monte Carlo
  sweep, because posterior coupling is only as trustworthy as the worse
  of the two.

- outcome.model:

  Function for fitting the outcome model, for callers who want a sampler
  the package does not provide. Overrides `outcome.prior`. Should take
  arguments Y, X, and mc. Available functions include
  [`bayes_lm`](https://t-momozaki.github.io/DRBayes/reference/bayes_lm.md)
  for continuous outcomes,
  [`bayes_logit`](https://t-momozaki.github.io/DRBayes/reference/bayes_logit.md)
  and
  [`bayes_probit`](https://t-momozaki.github.io/DRBayes/reference/bayes_probit.md)
  for binary outcomes, and
  [`bayes_lm_hs`](https://t-momozaki.github.io/DRBayes/reference/bayes_lm_hs.md),
  [`bayes_logit_hs`](https://t-momozaki.github.io/DRBayes/reference/bayes_logit_hs.md),
  and
  [`bayes_probit_hs`](https://t-momozaki.github.io/DRBayes/reference/bayes_probit_hs.md)
  for horseshoe priors. Set to NULL when using external samples. Default
  is NULL.

- ps.model:

  Function for fitting the propensity score model. Should take arguments
  Y, X, and mc. Available functions include
  [`bayes_logit`](https://t-momozaki.github.io/DRBayes/reference/bayes_logit.md),
  [`bayes_probit`](https://t-momozaki.github.io/DRBayes/reference/bayes_probit.md),
  [`bayes_logit_hs`](https://t-momozaki.github.io/DRBayes/reference/bayes_logit_hs.md),
  and
  [`bayes_probit_hs`](https://t-momozaki.github.io/DRBayes/reference/bayes_probit_hs.md)
  for binary treatments. Set to NULL when using external samples.
  Default is NULL.

- outcome.samples:

  Matrix of posterior samples for outcome model coefficients. Each row
  represents one MCMC iteration and each column represents a
  coefficient. The first column should be the intercept, followed by
  coefficients for covariates in the outcome model. Set to NULL when
  using internal model functions. Default is NULL.

- ps.samples:

  Matrix of posterior samples for propensity score model coefficients.
  Each row represents one MCMC iteration and each column represents a
  coefficient. The first column should be the intercept, followed by
  coefficients for covariates in the propensity score model. Set to NULL
  when using internal model functions. Default is NULL.

- outcome.mu:

  Matrix of posterior draws of the fitted outcome mean at the observed
  treatment and covariates, with one row per draw and one column per
  observation used. Supplying it, together with `outcome.mu1` and
  `outcome.mu0`, is the nonparametric outcome path of Section 5.4, in
  which the outcome model is whatever produced the draws and no
  coefficients are needed. It is mutually exclusive with
  `outcome.model`, `outcome.samples`, `outcome.prior`, `outcome.priors`
  and `link`, since none of them has anything left to describe, and the
  propensity score side is supplied as usual. Default is NULL.

- outcome.mu1, outcome.mu0:

  The same fitted means with the treatment set to 1 and to 0 for every
  observation, from which the g-computation (3.7) is formed. Both have
  the same dimensions as `outcome.mu`, and all three are supplied
  together. Default is NULL.

- outcome.priors:

  List of prior distribution parameters to pass to the outcome model
  function. For example, when using `bayes_lm`, you can specify
  `list(theta_prior = 1/50, sigma_prior = c(2, 1))` to customize the
  prior distributions for regression coefficients and error variance.
  Ignored when using external samples. Default is an empty list (uses
  function defaults).

- ps.priors:

  List of prior distribution parameters to pass to the propensity score
  model function. Ignored when using external samples. Default is an
  empty list (uses function defaults).

- seed:

  Integer or NULL. When given, the fit is reproducible and the caller's
  random number stream is restored afterwards.

- verbose:

  Logical. Whether to report progress through
  [`message`](https://rdrr.io/r/base/message.html). Default TRUE.

- na.action:

  Character string specifying how to handle missing values. Options are
  "na.omit" (default), "na.fail", or "na.exclude".

## Value

An object of class `"DRBayes"`, a list containing:

- g.comp:

  Numeric vector of posterior samples for the ATE using g-computation,
  that is, the untilted posterior. See *Which average* below for what
  the average is taken over.

- pc:

  Numeric vector of posterior samples for the ATE using posterior
  coupling. This is the doubly robust estimand and the quantity to
  report.

- family:

  Character string indicating the outcome model family used.

- link:

  Character string indicating the outcome model link function used.

- ps.link:

  Character string indicating the propensity score link function used.

- smc:

  List of sequential Monte Carlo diagnostics: the tilting parameter
  `lambda`, the achieved posterior mean `B_mean` of the moment
  condition, the tolerance `tol` it was compared against, the number of
  steps taken `n_steps` out of `max_steps`, and `converged`. When
  `converged` is FALSE the moment condition was never satisfied and the
  `pc` draws carry no double robustness guarantee.

- diagnostics:

  Data frame of convergence diagnostics, one row per parameter of each
  model: rank-normalised split R-hat, bulk and tail effective sample
  size, and the Monte Carlo standard error of the mean. NULL when
  `diagnostics = "none"`.

- outcome.prior, ps.prior:

  The priors used.

- method:

  Which tilting algorithm ran, `"smc"` or `"is"`. Worth reading rather
  than assuming: the nonparametric outcome path can only run Algorithm
  1, so it resolves to `"is"` whatever the default says.

- control:

  The resolved
  [`drbayes_control`](https://t-momozaki.github.io/DRBayes/reference/drbayes_control.md)
  object.

- particles:

  When `drbayes_control(keep_particles = TRUE)` asked for them, a list
  holding the tilted draws `outcome` and `ps`, the untilted draws
  `outcome_untilted` and `ps_untilted`, the design matrices and links in
  `model_data`, and the counterfactual design matrices `Z.lm1` and
  `Z.lm0`. This is what
  [`drbayes_sensitivity`](https://t-momozaki.github.io/DRBayes/reference/drbayes_sensitivity.md)
  reweights. On the nonparametric path there are no outcome
  coefficients, so it holds `ps`, `ps_untilted`, `model_data` and the
  resampled draw `index`, and the fit cannot be reweighted. NULL
  otherwise.

- call:

  The matched call.

- data_info:

  Summary of the processed data: `n_observations`, `n_treated`,
  `n_control`, the number of non-intercept columns `n_outcome_terms` and
  `n_ps_terms` in each design matrix, both formulas, and
  `missing_observations`.

## Details

The function supports both continuous outcomes (linear regression) and
binary outcomes (logistic and probit regression) through the family
argument, using an intuitive formula interface that integrates naturally
with the R ecosystem.

The function supports two modes of operation: internal posterior
sampling using built-in model functions, or external posterior samples
from other Bayesian software packages (e.g., Stan, JAGS, brms).

**Formula Interface:**

drbayes_pc uses R's standard formula notation, making it intuitive and
consistent with other statistical functions like lm() and glm(). The
formula interface automatically handles missing values, creates
appropriate design matrices, and manages factor variables.


    # Basic usage
    result <- drbayes_pc(
      outcome.formula = Y ~ A + X1 + X2,
      ps.formula = A ~ X1 + X2 + X3,
      data = my_data,
      outcome.model = bayes_lm,
      ps.model = bayes_logit
    )

    # With interactions
    result <- drbayes_pc(
      outcome.formula = Y ~ A + X1 + X2 + A:X1,
      ps.formula = A ~ X1 + X2 + X3 + I(X1^2),
      data = my_data,
      outcome.model = bayes_lm,
      ps.model = bayes_logit
    )

**Model Families and Link Functions:**

The function supports different outcome model types:

- **Automatic detection** (family=NULL): Determines family based on Y
  values. If all Y values are 0 or 1, uses "binomial"; otherwise uses
  "gaussian".

- **Gaussian family** (family="gaussian"): For continuous outcomes using
  linear regression with identity link. ATE is estimated as the mean
  difference \\E\[Y \| A = 1, X\] - E\[Y \| A = 0, X\]\\.

- **Binomial family** (family="binomial"): For binary outcomes using
  logistic regression (link="logit") or probit regression
  (link="probit"). ATE is estimated as the probability difference
  P(Y=1\|A=1,X) - P(Y=1\|A=0,X).

**Sampling Methods:**

The function automatically detects the sampling mode based on the
provided arguments:

- **Internal sampling**: Provide `outcome.model` and `ps.model`

- **External samples**: Provide `outcome.samples` and `ps.samples`

If both internal and external options are provided, external samples
take precedence with a warning message.

**Internal Sampling Mode:** The function generates posterior samples
from separate outcome and propensity score models using the provided
model functions and performs the following steps:

1.  Generates posterior samples from separate outcome and propensity
    score models

2.  Applies sequential Monte Carlo (SMC) to enforce the moment condition
    constraint

3.  Computes posterior distributions for the average treatment effect

**External Samples Mode:** When using external posterior samples (e.g.,
from Stan, JAGS, brms), ensure that:

- Samples are in matrix format with rows = iterations, columns =
  coefficients

- The first column contains intercept terms

- Coefficient ordering matches the column order in the design matrices
  created from formulas

- Both outcome.samples and ps.samples have the same number of iterations

The moment condition enforced is: \$\$E\[(A - \pi(X)) \cdot (Y - \mu(X))
/ (\pi(X)(1-\pi(X)))\] = 0\$\$

where \\\pi(X)\\ is the propensity score and \\\mu(X)\\ is the outcome
model.

## Which average

`g.comp` and `pc` are both posteriors for equation (3.7), the average of
the fitted contrast over the covariate vectors that were observed.
Section 2.1 of the paper defines the estimand as \\\tau = E\[Y_1 -
Y_0\]\\, an average over the population those vectors were drawn from.
The two are not the same quantity, and the difference is not in their
centre but in their spread: averaging over the observed covariates
treats their distribution as known, so these intervals do not carry the
uncertainty of not knowing it.

They coincide when the fitted contrast \\c_i = m_1(X_i) - m_0(X_i)\\ is
the same for every unit, and differ when it is not. Two things make it
vary: a treatment by covariate interaction in the outcome formula, and a
link that is not the identity. Only the first matters in practice.
Writing \\V\\ for the posterior variance of the average, putting
Dirichlet weights on the observations, as a Bayesian bootstrap would,
adds \\\mathrm{Var}\_i(c_i)/n\\ to it, so the interval widens by a
factor of \\\sqrt{1 + \mathrm{Var}\_i(c_i)/(nV)}\\. Measured on fits:

|                                  |       |            |
|----------------------------------|-------|------------|
| **outcome model**                | **n** | **factor** |
| `Y ~ A + X1 + X2 + A:X1`         | 500   | 1.11       |
| `Y ~ A + X1 + X2 + A:X1`         | 2000  | 1.21       |
| `Y ~ A + X1 + X2`, identity link | 500   | 1 exactly  |
| `Y ~ A + X1 + X2`, logit link    | 500   | 1.001      |
| `Y ~ A + X1 + X2`, logit link    | 2000  | 1.001      |

A logit or probit link does make the contrast vary, but by far too
little to matter against the posterior uncertainty; a treatment
interaction is the case to watch, and then how large the factor is
depends on how much the effect actually varies. The 1.11 and 1.21 above
are for a unit level effect of \\2 + 1.5 X_1\\; a milder interaction
gives a milder factor. For a fit of your own,
`drbayes_control(keep_particles = TRUE)` returns what is needed to
compute it directly, and the getting started vignette shows how.

None of this affects the point estimate or Theorem 1. Dirichlet weights
have mean \\1/n\\, so the posterior mean is unchanged, and the theorem
concerns the posterior mean of an estimator whose two averages differ by
\\O_p(n^{-1/2})\\. What changes is the credible interval, and which
estimand it is a credible interval for.

## Link Functions

- **Identity** (family="gaussian"): \\\mu = X\beta\\

- **Logit** (family="binomial", link="logit"): \\\mu =
  \mathrm{logit}^{-1}(X\beta) = \frac{e^{X\beta}}{1 + e^{X\beta}}\\

- **Probit** (family="binomial", link="probit"): \\\mu = \Phi(X\beta)\\
  where \\\Phi\\ is the standard normal CDF

## Nonparametric Outcome Models

Section 5.4 of the paper puts BART on the outcome side instead of a
linear predictor. Any outcome model can be coupled that way, as long as
it reports a posterior draw of its fitted mean: pass those draws as
`outcome.mu`, `outcome.mu1` and `outcome.mu0` in place of a
formula-based sampler. Fitted means from packages such as `dbarts` or
`BART` are the usual source, but nothing here depends on which software
produced them. The means are used exactly as supplied, on the scale of
the outcome itself, so no link is applied to them.


    # Draws of the fitted mean from any Bayesian outcome model, S by n
    fit <- drbayes_pc(
      outcome.formula = Y ~ A + X1 + X2,
      ps.formula = A ~ X1 + X2,
      data = data,
      outcome.mu = mu, outcome.mu1 = mu1, outcome.mu0 = mu0,
      ps.model = bayes_logit,
      chains = 1L, mc = 2 * nrow(mu) + 1000, bn = 1000, thin = 2
    )

The outcome formula is still needed, since it names the response and the
treatment and fixes which rows are complete; its right hand side plays
no further part. The moment condition (3.4) reads \\m\_{A_i}(X_i)\\ from
`outcome.mu`, and the average treatment effect (3.7) is the row mean of
`outcome.mu1` minus the row mean of `outcome.mu0`. All three matrices
need one row per propensity score draw that survives `bn`, `thin` and
`chains` – all four chains contribute by default, so `chains = 1L` is
usually what a supplied set of outcome draws wants – and one column per
observation used, in the order the rows of `data` are in.

The draws are resampled systematically, with replacement, which is the
standard choice and leaves the resampled cloud an unbiased
representation of the tilted posterior. The footnote to Table 2 of the
paper reports something different for its own BART results, "importance
sampling without replacement", keeping 10,000 of 20,000 draws. Selecting
without replacement guarantees that many distinct draws survive, which
matters here because nothing replenishes the cloud, but it does so at
the cost of inclusion probabilities that are no longer proportional to
the weights. The two are not interchangeable, and this package
implements the former.

Only `method = "is"` is available. Algorithm 2 rejuvenates its particles
with a Gaussian kernel on the outcome model's parameter vector, and a
draw of fitted values has none, so the tilting reweights and resamples
the draw index instead, carrying the propensity score draw and the three
rows of fitted means together. Nothing then replenishes the cloud, which
makes `ess` in `$smc` the number to watch: a fit resting on fewer than
`ess_frac` of the draws is reported as not converged.

[`drbayes_sensitivity`](https://t-momozaki.github.io/DRBayes/reference/drbayes_sensitivity.md)
does not apply to a fit made this way, because Algorithm 3 reweights the
outcome model's coefficient draws.

## External Software Integration

This function is designed to work with posterior samples from various
Bayesian software:

**Stan/cmdstanr:**


    # Extract samples
    outcome.samples <- fit$draws("beta", format = "matrix")
    ps.samples <- fit$draws("gamma", format = "matrix")

    # Use in drbayes_pc
    result <- drbayes_pc(
      outcome.formula = Y ~ A + X1 + X2 + A:X1,
      ps.formula = A ~ X1 + X2 + X3,
      data = data,
      outcome.samples = outcome.samples,
      ps.samples = ps.samples
    )

**brms:**


    # Extract fixed effects
    outcome.samples <- as.matrix(brms_fit)[, 1:ncol(model.matrix(outcome.formula, data))]
    ps.samples <- as.matrix(brms_fit)[, 1:ncol(model.matrix(ps.formula, data))]

    # Use in drbayes_pc
    result <- drbayes_pc(
      outcome.formula = Y ~ A + X1 + X2,
      ps.formula = A ~ X1 + X2 + X3,
      data = data,
      outcome.samples = outcome.samples,
      ps.samples = ps.samples
    )

**JAGS/R2jags:**


    # Extract parameter samples
    outcome.samples <- jags.fit$BUGSoutput$sims.matrix[, grep("beta", colnames(...))]
    ps.samples <- jags.fit$BUGSoutput$sims.matrix[, grep("gamma", colnames(...))]

    # Use in drbayes_pc
    result <- drbayes_pc(
      outcome.formula = Y ~ A + X1 + X2,
      ps.formula = A ~ X1 + X2 + X3,
      data = data,
      outcome.samples = outcome.samples,
      ps.samples = ps.samples
    )

## References

Orihara, S., Momozaki, T., & Sugasawa, S. (2025). Bayesian Doubly Robust
Causal Inference via Posterior Coupling. arXiv preprint
arXiv:2506.04868.

## See also

[`bayes_lm`](https://t-momozaki.github.io/DRBayes/reference/bayes_lm.md)
for Bayesian linear regression,
[`bayes_logit`](https://t-momozaki.github.io/DRBayes/reference/bayes_logit.md)
for Bayesian logistic regression,
[`bayes_probit`](https://t-momozaki.github.io/DRBayes/reference/bayes_probit.md)
for Bayesian probit regression,
[`bayes_lm_hs`](https://t-momozaki.github.io/DRBayes/reference/bayes_lm_hs.md)
for horseshoe linear regression,
[`bayes_logit_hs`](https://t-momozaki.github.io/DRBayes/reference/bayes_logit_hs.md)
for horseshoe logistic regression,
[`bayes_probit_hs`](https://t-momozaki.github.io/DRBayes/reference/bayes_probit_hs.md)
for horseshoe probit regression,
[`drbayes_bb`](https://t-momozaki.github.io/DRBayes/reference/drbayes_bb.md)
for Bayesian bootstrap approach.

## Examples

``` r
set.seed(123)
n <- 150
dat <- data.frame(X1 = rnorm(n), X2 = rnorm(n), X3 = rnorm(n))
dat$A <- rbinom(n, 1, plogis(0.2 + 0.5 * dat$X1 - 0.3 * dat$X2))
dat$Y <- 1 + 1.5 * dat$A + 0.8 * dat$X1 + 0.6 * dat$X2 +
  0.4 * dat$A * dat$X1 + rnorm(n)

# The treatment effect varies with X1, so the interaction belongs in the
# outcome formula: the g-formula averages the unit level contrast over the
# covariates, and it can only average what the model can express.
fit <- drbayes_pc(Y ~ A + X1 + X2 + A:X1, A ~ X1 + X2 + X3, dat,
                  outcome.model = bayes_lm, ps.model = bayes_logit,
                  control = drbayes_control(mc = 800, bn = 300))
#> Formula processing complete:
#>   Observations used: 150 out of 150
#>   Treatment group: 83 | Control group: 67
#>   Outcome model variables: 4
#>   Propensity score variables: 3
#> Auto-detected continuous outcome: using family='gaussian'
#> Generating posterior samples using provided models...
#> Using 1000 posterior samples for analysis...
fit
#> Bayesian doubly robust causal inference via posterior coupling
#> 
#> Call:
#> drbayes_pc(outcome.formula = Y ~ A + X1 + X2 + A:X1, ps.formula = A ~ 
#>     X1 + X2 + X3, data = dat, control = drbayes_control(mc = 800, 
#>     bn = 300), outcome.model = bayes_lm, ps.model = bayes_logit)
#> 
#> Outcome model:    gaussian family, identity link
#> Propensity score: logit link
#> Observations:     150 (83 treated, 67 control)
#> Posterior draws:  1000
#> 
#> Average treatment effect
#>         mean  2.5% 97.5%
#> pc     1.377 1.054 1.703
#> g.comp 1.366 1.024 1.715
#> 
#> Report pc, the doubly robust estimand. g.comp is the untilted g-formula
#> posterior, which relies on the outcome model alone.
#> 
#> Both are posteriors for the average causal effect over the 150 covariate
#> vectors observed, which is what equation (3.7) averages. See ?drbayes_pc on
#> when that differs from the average over the population they were drawn from.
#> 
#> Posterior coupling: lambda = -0.13, moment condition met, |mean B_n| = 0.000345
#> below the tolerance 0.00573, after 13 of 1000 steps.
#> 
#> Convergence diagnostics (worst over both models)
#>   Largest R-hat: 1.006 (outcome: A:X1)
#>   Smallest ESS:  812 (propensity score: (Intercept))
summary(fit)
#> Bayesian doubly robust causal inference via posterior coupling
#> 
#> Call:
#> drbayes_pc(outcome.formula = Y ~ A + X1 + X2 + A:X1, ps.formula = A ~ 
#>     X1 + X2 + X3, data = dat, control = drbayes_control(mc = 800, 
#>     bn = 300), outcome.model = bayes_lm, ps.model = bayes_logit)
#> 
#> Outcome model:    gaussian family, identity link
#> Propensity score: logit link
#> Observations:     150 (83 treated, 67 control)
#> Posterior draws:  1000
#> Outcome formula:  Y ~ A + X1 + X2 + A:X1
#> PS formula:       A ~ X1 + X2 + X3
#> 
#> Posterior summary of the average treatment effect
#>         mean    sd  2.5%   50% 97.5%
#> pc     1.377 0.171 1.054 1.396 1.703
#> g.comp 1.366 0.172 1.024 1.363 1.715
#> 
#> pc is the doubly robust estimand of equation (3.8) and the one to report.
#> g.comp is the untilted g-formula posterior, lambda = 0, which relies on the
#> outcome model alone. Both average over the same covariates, so their spreads
#> are comparable.
#> 
#> The covariates they average over are the observed ones. Where the fitted
#> contrast varies from unit to unit, which a treatment interaction or a
#> non-identity link both cause, this is not the same quantity as the average over
#> the population, and these intervals do not carry the uncertainty of not knowing
#> the covariate distribution. ?drbayes_pc gives the size of the difference.
#> 
#> Sequential Monte Carlo
#>   Tilting parameter lambda: -0.13
#>   Posterior mean of B_n:    0.000345
#>   Tolerance:                0.00573
#>   Steps taken:              13 of 1000
#>   Moment condition met:     yes
#> 
#> Convergence diagnostics
#>             model   parameter    mean     sd rhat ess_bulk ess_tail mcse_mean
#>           outcome (Intercept)  0.9988 0.1306 1.00      984      950   0.00417
#>           outcome           A  1.3807 0.1726 1.00      949      971   0.00559
#>           outcome          X1  0.6883 0.1212 1.00      934      816   0.00397
#>           outcome          X2  0.6781 0.0912 1.00     1101      848   0.00276
#>           outcome        A:X1  0.6014 0.1791 1.01     1007      830   0.00559
#>  propensity score (Intercept)  0.2867 0.1749 1.00      948      812   0.00568
#>  propensity score          X1  0.4863 0.2010 1.00      877      916   0.00680
#>  propensity score          X2 -0.3760 0.1874 1.00      841      887   0.00645
#>  propensity score          X3  0.0416 0.1757 1.00      931      895   0.00575
#> 
#> All parameters meet R-hat < 1.01 and ESS > 400.

# Whether the tilting did anything is in `$smc`. `lambda` is zero only when
# the untilted posterior already satisfied the moment condition, and
# `converged` says whether it was satisfied at all; when it was not, the
# draws in `pc` carry no double robustness guarantee.
fit$smc[c("lambda", "B_mean", "tol", "converged", "n_steps")]
#> $lambda
#> [1] -0.13
#> 
#> $B_mean
#> [1] 0.0003445655
#> 
#> $tol
#> [1] 0.005730101
#> 
#> $converged
#> [1] TRUE
#> 
#> $n_steps
#> [1] 13
#> 

# A binary outcome. The estimand is then a difference of probabilities.
dat$Ybin <- rbinom(n, 1, plogis(0.5 + 0.8 * dat$A + 0.4 * dat$X1))
drbayes_pc(Ybin ~ A + X1 + X2, A ~ X1 + X2 + X3, dat,
           outcome.model = bayes_logit, ps.model = bayes_logit,
           control = drbayes_control(mc = 800, bn = 300))
#> Formula processing complete:
#>   Observations used: 150 out of 150
#>   Treatment group: 83 | Control group: 67
#>   Outcome model variables: 3
#>   Propensity score variables: 3
#> Auto-detected binary outcome: using family='binomial'
#> Generating posterior samples using provided models...
#> Using 1000 posterior samples for analysis...
#> Bayesian doubly robust causal inference via posterior coupling
#> 
#> Call:
#> drbayes_pc(outcome.formula = Ybin ~ A + X1 + X2, ps.formula = A ~ 
#>     X1 + X2 + X3, data = dat, control = drbayes_control(mc = 800, 
#>     bn = 300), outcome.model = bayes_logit, ps.model = bayes_logit)
#> 
#> Outcome model:    binomial family, logit link
#> Propensity score: logit link
#> Observations:     150 (83 treated, 67 control)
#> Posterior draws:  1000
#> 
#> Average treatment effect
#>          mean    2.5%  97.5%
#> pc     0.0345 -0.0505 0.1241
#> g.comp 0.0275 -0.1202 0.1775
#> 
#> Report pc, the doubly robust estimand. g.comp is the untilted g-formula
#> posterior, which relies on the outcome model alone.
#> 
#> Both are posteriors for the average causal effect over the 150 covariate
#> vectors observed, which is what equation (3.7) averages. See ?drbayes_pc on
#> when that differs from the average over the population they were drawn from.
#> 
#> Posterior coupling: lambda = -1.97, moment condition met, |mean B_n| = 0.00254
#> below the tolerance 0.0026, after 197 of 1000 steps.
#> 
#> Convergence diagnostics (worst over both models)
#>   Largest R-hat: 1.003 (outcome: A)
#>   Smallest ESS:  778 (propensity score: (Intercept))

# Transformations are rebuilt through the fitted terms, so the poly() basis
# in the two counterfactual designs is the one the model was fitted with
# rather than one re-derived from data that has changed.
drbayes_pc(Y ~ A + log(abs(X1) + 1) + poly(X3, 2) + A:X1,
           A ~ X1 + X2 + X3, dat,
           outcome.model = bayes_lm, ps.model = bayes_logit,
           control = drbayes_control(mc = 600, bn = 200))
#> Formula processing complete:
#>   Observations used: 150 out of 150
#>   Treatment group: 83 | Control group: 67
#>   Outcome model variables: 5
#>   Propensity score variables: 3
#> Auto-detected continuous outcome: using family='gaussian'
#> Generating posterior samples using provided models...
#> Using 800 posterior samples for analysis...
#> Bayesian doubly robust causal inference via posterior coupling
#> 
#> Call:
#> drbayes_pc(outcome.formula = Y ~ A + log(abs(X1) + 1) + poly(X3, 
#>     2) + A:X1, ps.formula = A ~ X1 + X2 + X3, data = dat, control = drbayes_control(mc = 600, 
#>     bn = 200), outcome.model = bayes_lm, ps.model = bayes_logit)
#> 
#> Outcome model:    gaussian family, identity link
#> Propensity score: logit link
#> Observations:     150 (83 treated, 67 control)
#> Posterior draws:  800
#> 
#> Average treatment effect
#>         mean  2.5% 97.5%
#> pc     1.300 0.915 1.670
#> g.comp 1.298 0.911 1.708
#> 
#> Report pc, the doubly robust estimand. g.comp is the untilted g-formula
#> posterior, which relies on the outcome model alone.
#> 
#> Both are posteriors for the average causal effect over the 150 covariate
#> vectors observed, which is what equation (3.7) averages. See ?drbayes_pc on
#> when that differs from the average over the population they were drawn from.
#> 
#> Posterior coupling: lambda = -0.04, moment condition met, |mean B_n| = 0.00852
#> below the tolerance 0.00944, after 4 of 1000 steps.
#> 
#> Convergence diagnostics (worst over both models)
#>   Largest R-hat: 1.004 (propensity score: (Intercept))
#>   Smallest ESS:  595 (outcome: (Intercept))

# \donttest{
# Draws from elsewhere, supplied instead of fitting. Anything that returns a
# draws-by-parameters matrix will do -- Stan, JAGS, brms -- provided its
# columns are those of model.matrix() for the same formula, and the
# propensity score draws are on the scale `ps.link` names. Here they come
# from this package's own samplers, which take the covariates without an
# intercept and add one. Two fits rather than one, so it is the slowest
# thing on this page and sits outside the timed examples.
Zo <- model.matrix(Y ~ A + X1 + X2, dat)[, -1, drop = FALSE]
Zp <- model.matrix(A ~ X1 + X2 + X3, dat)[, -1, drop = FALSE]
otc <- bayes_lm(dat$Y, Zo, mc = 1500, chains = 2L)
ps  <- bayes_logit(dat$A, Zp, mc = 1500, chains = 2L)
pool <- function(x) {
  matrix(x, ncol = dim(x)[3], dimnames = list(NULL, dimnames(x)[[3]]))
}
drbayes_pc(Y ~ A + X1 + X2, A ~ X1 + X2 + X3, dat,
           outcome.samples = pool(otc), ps.samples = pool(ps),
           control = drbayes_control(bn = 500))
#> Formula processing complete:
#>   Observations used: 150 out of 150
#>   Treatment group: 83 | Control group: 67
#>   Outcome model variables: 3
#>   Propensity score variables: 3
#> Auto-detected continuous outcome: using family='gaussian'
#> Using externally provided posterior samples...
#> Using 1250 posterior samples for analysis...
#> Bayesian doubly robust causal inference via posterior coupling
#> 
#> Call:
#> drbayes_pc(outcome.formula = Y ~ A + X1 + X2, ps.formula = A ~ 
#>     X1 + X2 + X3, data = dat, control = drbayes_control(bn = 500), 
#>     outcome.samples = pool(otc), ps.samples = pool(ps))
#> 
#> Outcome model:    gaussian family, identity link
#> Propensity score: logit link
#> Observations:     150 (83 treated, 67 control)
#> Posterior draws:  1250
#> 
#> Average treatment effect
#>         mean  2.5% 97.5%
#> pc     1.402 1.172 1.608
#> g.comp 1.344 0.975 1.689
#> 
#> Report pc, the doubly robust estimand. g.comp is the untilted g-formula
#> posterior, which relies on the outcome model alone.
#> 
#> Both are posteriors for the average causal effect over the 150 covariate
#> vectors observed, which is what equation (3.7) averages. See ?drbayes_pc on
#> when that differs from the average over the population they were drawn from.
#> 
#> Posterior coupling: lambda = -4.34, moment condition met, |mean B_n| = 0.00533
#> below the tolerance 0.00553, after 434 of 1000 steps.
#> 
#> Convergence diagnostics (worst over both models)
#>   Largest R-hat: 1.001 (outcome: X2)
#>   Smallest ESS:  1034 (outcome: X2)
# }
```
