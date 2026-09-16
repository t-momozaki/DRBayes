# Confounder Selection with Posterior Coupling

Implements Algorithm 4 of Orihara, Momozaki and Sugasawa (2025). Both
models are fitted under shrinkage priors, the covariates whose
propensity score coefficient survives that shrinkage are selected, and
the two posteriors are then coupled through a moment condition
restricted to the selected set.

## Usage

``` r
drbayes_select(
  outcome.formula,
  ps.formula,
  data,
  threshold = 0.01,
  family = NULL,
  link = NULL,
  ps.link = "logit",
  outcome.prior = c("horseshoe", "normal"),
  ps.prior = c("horseshoe", "normal"),
  mc = 5000,
  bn = 1000,
  thin = 2,
  chains = 4L,
  method = c("smc", "is"),
  rejuvenate = c("all", "selected"),
  control = drbayes_control(),
  diagnostics = c("warn", "error", "none"),
  outcome.model = NULL,
  ps.model = NULL,
  outcome.priors = list(),
  ps.priors = list(),
  seed = NULL,
  verbose = TRUE,
  na.action = "na.omit"
)
```

## Arguments

- outcome.formula:

  Formula for the outcome model, for example `Y ~ A + X1 + X2 + A:X1`,
  with the response on the left and the treatment among the predictors
  on the right.

- ps.formula:

  Formula for the propensity score model, for example
  `A ~ X1 + X2 + X3`, with the treatment on the left and the candidate
  confounders on the right. These candidates are what the selection
  chooses from.

- data:

  Data frame holding every variable used by either formula.

- threshold:

  Numeric. A covariate is selected when the posterior mean of its
  propensity score coefficient is at least `threshold` in absolute
  value, that is \\S = \\j : \|\bar{\alpha}\_j\| \ge\\
  `threshold`\\\\\\. Default 0.01, the value used in the right heart
  catheterization analysis of the paper. Zero selects every candidate.
  The scale is that of the propensity score linear predictor, so it
  depends on how the covariates are scaled.

- family:

  Character string for the outcome model family, either `"gaussian"` or
  `"binomial"`. NULL, the default, detects a binary response as binomial
  and anything else as gaussian.

- link:

  Character string for the link of the OUTCOME model: `"identity"`,
  `"logit"` or `"probit"`. NULL, the default, takes the canonical link
  of `family`.

- ps.link:

  Character string for the link of the PROPENSITY SCORE model, either
  `"logit"` (default) or `"probit"`. The moment condition turns the
  propensity score draws into probabilities with this link.

- outcome.prior, ps.prior:

  Character. Prior for the coefficients of each model in Step 1, either
  `"horseshoe"` (default) or `"normal"`. Step 1 of Algorithm 4 asks for
  a shrinkage prior, "such as the horseshoe prior", which is why the
  horseshoe is the default rather than the only option. The selection
  reads the propensity score coefficients, so `ps.prior = "normal"`
  shrinks nothing and leaves the threshold selecting on unregularised
  posterior means.

- mc, bn, thin, chains:

  Integers. Number of MCMC iterations per chain, the burn-in to discard,
  the thinning interval, and the number of chains. They also live on the
  `control` object; naming one here overrides it.

- method:

  Character. How Step 2 finds the tilting parameter. `"smc"` is
  Algorithm 2, which walks a grid of tilting parameters and rejuvenates
  the particles at each step. `"is"` is Algorithm 1, a single importance
  sampling step, which is faster but degenerates when the two posteriors
  put little mass where the moment condition holds. Step 2 of Algorithm
  4 admits either. Default `"smc"`.

- rejuvenate:

  Character, `"all"` (default) or `"selected"`, choosing which
  coefficients the Liu and West (2001) kernel refreshes at each step of
  the sweep. See "Rejuvenating the carried block" below; the default
  departs deliberately from the literal wording of Algorithm 4, and
  `"selected"` reproduces it. Ignored when `method = "is"`, which moves
  no coefficient at all and records this argument as NA.

- control:

  List of tuning parameters from
  [`drbayes_control`](https://t-momozaki.github.io/DRBayes/reference/drbayes_control.md).
  Its `n_steps` is the one to think about here: it sets how many times
  the sweep resamples, and under `rejuvenate = "selected"` the carried
  block loses ancestry at every one of those resampling steps. The
  default of 1000 is a reasonable choice under the default `rejuvenate`,
  which refreshes the whole particle and so has no ancestry to lose.

- diagnostics:

  Character. What to do when the Step 1 draws fail the convergence
  criteria: `"warn"` (default), `"error"`, or `"none"` to skip the
  check.

- outcome.model, ps.model:

  Functions fitting the Step 1 models, for callers who want a sampler
  the package does not provide, such as
  [`bayes_stan`](https://t-momozaki.github.io/DRBayes/reference/bayes_stan.md)
  or a wrapper around draws obtained elsewhere. Each should take `Y`,
  `X`, `mc` and `chains` and return an iterations by chains by
  parameters array, the layout the built-in samplers return. Given, they
  override `outcome.prior` and `ps.prior`, and the fit records the prior
  it did not use as `"supplied"`. Default NULL, which selects the
  built-in sampler for the link and the prior. There is no
  `outcome.samples` argument here as there is in
  [`drbayes_pc`](https://t-momozaki.github.io/DRBayes/reference/drbayes_pc.md):
  draws already in hand are passed as
  `function(Y, X, mc, chains) draws`.

- outcome.priors, ps.priors:

  Lists of further arguments for the Step 1 samplers, for example
  `list(p0 = 10)` to change the horseshoe's guess at how many
  coefficients are non-zero. The treatment columns of the outcome model
  are exempted from shrinkage automatically; passing `unshrunk` yourself
  replaces that choice.

- seed:

  Integer or NULL. When given, the fit is reproducible and the caller's
  random number stream is restored afterwards.

- verbose:

  Logical. Whether to report progress, including how many covariates
  were selected, through
  [`message`](https://rdrr.io/r/base/message.html). Default TRUE.

- na.action:

  One of `"na.omit"` (default), `"na.fail"` or `"na.exclude"`.

## Value

An object of class `"DRBayes_select"`, which inherits from `"DRBayes"`
and therefore works with
[`summary.DRBayes`](https://t-momozaki.github.io/DRBayes/reference/summary.DRBayes.md)
and
[`plot.DRBayes`](https://t-momozaki.github.io/DRBayes/reference/plot.DRBayes.md).
It holds everything
[`drbayes_pc`](https://t-momozaki.github.io/DRBayes/reference/drbayes_pc.md)
returns, and in addition:

- selected:

  Character vector of the selected covariates, named as the columns of
  the propensity score design matrix, so a factor appears once per dummy
  column.

- threshold:

  The threshold used.

- n_candidates:

  Number of candidate covariates the selection chose from, again
  counting dummy columns separately.

- posterior_means:

  Named numeric vector of the Step 1 posterior means \\\bar{\alpha}\_j\\
  the selection was made on, one per candidate. The intercept is not a
  candidate and is not included.

- method, rejuvenate:

  Which algorithm ran Step 2, and which coefficients its kernel
  refreshed. `rejuvenate` is NA under `method = "is"`, which applies no
  kernel and so refreshes none.

Its `smc` component carries one diagnostic beyond the ones
[`drbayes_pc`](https://t-momozaki.github.io/DRBayes/reference/drbayes_pc.md)
records: `n_ancestors_cumulative`, the number of distinct Step 1 draws
the returned particles still descend from. Unlike `n_ancestors_min`,
which counts the survivors of the worst single resampling step, this
accumulates over the whole sweep, and it is the count that describes a
carried block. It is reported whenever some coordinate of a particle is
a verbatim copy of a Step 1 draw: under `method = "is"`, under
`rejuvenate = "selected"` when the selection left a coefficient outside
the updated block, and when the constraint already held so that no sweep
ran. It is NA otherwise.

`smc$converged` means here what it means in
[`drbayes_pc`](https://t-momozaki.github.io/DRBayes/reference/drbayes_pc.md),
that the moment condition was met on draws numerous enough to describe a
posterior, and it is FALSE when the carried block has collapsed however
well the constraint is satisfied. The `print` method separates the two
cases; [`summary()`](https://rdrr.io/r/base/summary.html) reports only
the flag.

## Details

Applying shrinkage separately to the outcome and the propensity score
model is the obvious thing to do and it is what causes the problem this
function solves: a weak confounder is shrunk away in both models at
once, and the resulting regularization-induced confounding biases the
treatment effect (Hahn et al., 2018). Selecting on the treatment model
and then coupling gives such a covariate a second chance, because the
tilting updates its outcome coefficient jointly with the propensity
score.

**The algorithm:**

*Step 1* generates posterior samples for both models under shrinkage
priors (Carvalho et al., 2010) and forms the posterior means
\\\bar{\alpha}\\ of the propensity score coefficients. The selected set
is \\S = \\j : \|\bar{\alpha}\_j\| \ge \tau\\\\ for the threshold
\\\tau\\.

*Step 2* replaces the moment condition (3.4) with \$\$B_n^S(\alpha_S,
\beta_S) = \frac{1}{n} \sum\_{i=1}^{n} \frac{A_i - e(X\_{Si};
\alpha_S)}{e(X\_{Si}; \alpha_S)(1 - e(X\_{Si}; \alpha_S))} (Y_i -
m\_{A_i}(X_i; \beta_S, \beta\_{S^c}))\$\$ and runs Algorithm 1 or
Algorithm 2 against it. The propensity score is evaluated on the
selected covariates alone, while the outcome mean keeps every covariate.

**Three points the algorithm turns on:**

The selection is made on the coefficients of the PROPENSITY SCORE model,
not of the outcome model. That is the whole point: a covariate that
predicts treatment is a candidate confounder however weakly it predicts
the outcome, and it is exactly such a covariate that separate shrinkage
discards.

The outcome coefficients outside the selected set, \\\beta\_{S^c}\\, are
not dropped, and the g-computation is done with the full \\\beta\\.
Using the selected block alone would silently change the estimand into
the average treatment effect of a different, smaller outcome model.

The treatment is never a candidate for selection and is never shrunk.
Its main effect and every interaction with it are exempted from the
shrinkage in Step 1 and are always part of the block the sweep updates,
since shrinking or freezing them would attenuate the very effect being
estimated.

**Which outcome coefficients a selected covariate reaches:**

A covariate is selected as a column of the propensity score design
matrix, and its outcome coefficients are found through the terms
structure of the two formulas rather than by matching column names. The
propensity score column is traced back to the term it came from, that
term to the variables it is built from, and those variables forward to
every outcome column they contribute to. So selecting `X1` in
`A ~ X1 + X2` reaches `log(X1)`, `poly(X1, 2)` and `A:X1` in the outcome
model, which matching names would not. A selected covariate that reaches
no outcome column at all is reported: the outcome model does not adjust
for it, so the coupling has no coefficient of it to revive.

**Rejuvenating the carried block:**

Algorithm 4 says that in Step 2 "only \\(\alpha_S, \beta_S)\\ are
updated with \\\beta\_{S^c}\\ being unchanged". Taken literally as an
instruction to the sequential Monte Carlo sweep, that means applying the
Liu and West (2001) kernel to the selected block only. The sweep still
resamples whole particles, so the carried block is never refreshed and
its ancestry collapses geometrically. Over the thousand steps
[`drbayes_control`](https://t-momozaki.github.io/DRBayes/reference/drbayes_control.md)
takes by default, a cloud of several thousand particles can come back
tracing to a dozen distinct Step 1 draws, and on some data to one,
leaving \\\beta\_{S^c}\\ a point mass reported as a posterior.

The default `rejuvenate = "all"` therefore applies the kernel to every
coefficient, and departs from the literal wording on purpose. The kernel
is a rejuvenation step: it shrinks each particle towards the cloud mean
and adds back the variance that shrinking removed, so it preserves the
first two moments of the cloud and moves no parameter towards the
constraint. What does move a parameter towards the constraint is the
importance weight, and \\\beta\_{S^c}\\ enters those weights whatever
this argument is set to, because it enters \\B_n^S\\ through
\\m\_{A_i}(X_i; \beta_S, \beta\_{S^c})\\. Rejuvenating the carried block
therefore does not update it in the sense Step 2 excludes, while
freezing it destroys the posterior it is there to represent. The
propensity score coefficients outside \\S\\ are refreshed on the same
grounds, and with more reason: they do not enter \\B_n^S\\ at all, so
the tilting leaves their conditional distribution alone and a frozen
copy of it is the one thing that would not.

`rejuvenate = "selected"` reproduces the literal reading for anyone who
wants it. The kernel then refreshes the selected block alone, so both
\\\beta\_{S^c}\\ and \\\alpha\_{S^c}\\ are carried verbatim and both
collapse together: the propensity score draws in `fit$particles$ps` are
as much a point mass as the outcome ones, and the warning names how many
coefficients of each are affected. It reports
`fit$smc$n_ancestors_cumulative`, and reports the fit as not converged
once the carried block has collapsed onto few enough draws that its
posterior spread is not to be believed. That floor is `ess_frac` of
[`drbayes_control`](https://t-momozaki.github.io/DRBayes/reference/drbayes_control.md),
the same one the tilting holds its own counts to. Fewer steps in the
tilting parameter, through `drbayes_control(n_steps = )`, mean fewer
resampling rounds and so a slower collapse, but they do not remove it.

## References

Orihara, S., Momozaki, T., & Sugasawa, S. (2025). Bayesian Doubly Robust
Causal Inference via Posterior Coupling. arXiv preprint
arXiv:2506.04868.

Carvalho, C. M., Polson, N. G., & Scott, J. G. (2010). The horseshoe
estimator for sparse signals. Biometrika, 97(2), 465-480.

Hahn, P. R., Carvalho, C. M., Puelz, D., & He, J. (2018). Regularization
and confounding in linear regression for treatment effect estimation.
Bayesian Analysis, 13(1), 163-182.

Liu, J., & West, M. (2001). Combined parameter and state estimation in
simulation-based filtering. In A. Doucet, N. de Freitas, & N. Gordon
(Eds.), Sequential Monte Carlo Methods in Practice (pp. 197-223).
Springer.

Ning, Y., Sida, P., & Imai, K. (2020). Robust estimation of causal
effects via a high-dimensional covariate balancing propensity score.
Biometrika, 107(3), 533-554.

## See also

[`drbayes_pc`](https://t-momozaki.github.io/DRBayes/reference/drbayes_pc.md)
for posterior coupling without selection,
[`drbayes_control`](https://t-momozaki.github.io/DRBayes/reference/drbayes_control.md)
for the tuning parameters, and
[`bayes_lm_hs`](https://t-momozaki.github.io/DRBayes/reference/bayes_lm_hs.md)
for the horseshoe samplers Step 1 uses.

## Examples

``` r
set.seed(1)
n <- 150
dat <- data.frame(matrix(rnorm(n * 6), n, 6))
names(dat) <- paste0("X", 1:6)

# X1 predicts treatment strongly and the outcome only weakly, so separate
# shrinkage discards it and the treatment effect absorbs its confounding.
dat$A <- rbinom(n, 1, plogis(1.5 * dat$X1 + 0.8 * dat$X2))
dat$Y <- 2 * dat$A + 0.5 * dat$X1 + 1.2 * dat$X2 + rnorm(n, sd = 3)

covariates <- paste(paste0("X", 1:6), collapse = " + ")
outcome <- as.formula(paste("Y ~ A +", covariates))
ps      <- as.formula(paste("A ~", covariates))

fit <- drbayes_select(
  outcome.formula = outcome,
  ps.formula      = ps,
  data            = dat,
  threshold       = 0.1,
  mc = 800, bn = 200, chains = 2L, seed = 2, verbose = FALSE,
  control = drbayes_control(n_steps = 150)
)
#> Warning: The Step 1 posterior draws do not meet the convergence criteria (R-hat < 1.01, ESS > 400): outcome model: max R-hat 1.012 (X4); min ESS 183.1 (X1, X3, X4) | propensity score model: max R-hat 1.027 (X2, X5); min ESS 100.6 (X1, X2, X5). The selection is made on these draws, so it is not trustworthy either. An effective sample size below 400 also makes R-hat itself unreliable, so a small R-hat here is not evidence of convergence. Raise mc, or thin less. Inspect $diagnostics.
fit
#> Bayesian doubly robust causal inference via posterior coupling
#> 
#> Call:
#> drbayes_select(outcome.formula = outcome, ps.formula = ps, data = dat, 
#>     threshold = 0.1, mc = 800, bn = 200, chains = 2L, control = drbayes_control(n_steps = 150), 
#>     verbose = FALSE)
#> 
#> Outcome model:    gaussian family, identity link
#> Propensity score: logit link
#> Observations:     150 (76 treated, 74 control)
#> Posterior draws:  600
#> 
#> Average treatment effect
#>         mean  2.5% 97.5%
#> pc     1.568 0.585 2.619
#> g.comp 1.764 0.685 2.835
#> 
#> Report pc, the doubly robust estimand. g.comp is the untilted g-formula
#> posterior, which relies on the outcome model alone.
#> 
#> Both are posteriors for the average causal effect over the 150 covariate
#> vectors observed, which is what equation (3.7) averages. See ?drbayes_pc on
#> when that differs from the average over the population they were drawn from.
#> 
#> Posterior coupling: lambda = 0.8, moment condition met, |mean B_n| = 0.0252
#> below the tolerance 0.0258, after 12 of 150 steps.
#> 
#> Convergence diagnostics (worst over both models)
#>   Largest R-hat: 1.027 (propensity score: X2)  <- above 1.01
#>   Smallest ESS:  101 (propensity score: X2)  <- below 400
#>   The draws do not meet the convergence criteria, so the estimates above are
#>   not to be trusted. Raise mc, or use more chains from more dispersed
#>   starting values. See summary() for the full table.
#> 
#> Confounder selection (Algorithm 4)
#>   Selected 3 of 6 candidate covariate(s) at |alpha_bar| >= 0.1
#>     Selected: X1, X2, X5.
#>     Not selected: X3, X4, X6.
#>   The moment condition was restricted to the selected covariates and the
#>   particles were rejuvenated as a whole. Every outcome coefficient is in the
#>   g-computation, so the estimand is unchanged.
fit$selected  # the two covariates that drive treatment
#> [1] "X1" "X2" "X5"

# The literal reading of Step 2, which carries the coefficients outside the
# selected set at their Step 1 draws. Its ancestry count reports how many
# of those draws that block still holds, and fewer steps in the tilting
# parameter leave it more of them. The warning it raises is the point of
# the comparison, not a problem with the data.
literal <- drbayes_select(
  outcome.formula = outcome,
  ps.formula      = ps,
  data            = dat,
  threshold       = 0.1,
  rejuvenate      = "selected",
  mc = 800, bn = 200, chains = 2L, seed = 2, verbose = FALSE,
  control = drbayes_control(n_steps = 150)
)
#> Warning: The Step 1 posterior draws do not meet the convergence criteria (R-hat < 1.01, ESS > 400): outcome model: max R-hat 1.012 (X4); min ESS 183.1 (X1, X3, X4) | propensity score model: max R-hat 1.027 (X2, X5); min ESS 100.6 (X1, X2, X5). The selection is made on these draws, so it is not trustworthy either. An effective sample size below 400 also makes R-hat itself unreliable, so a small R-hat here is not evidence of convergence. Raise mc, or thin less. Inspect $diagnostics.
literal$smc$n_ancestors_cumulative
#> [1] 147
```
