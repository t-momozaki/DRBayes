# Print a confounder selection fit

Everything
[`print.DRBayes`](https://t-momozaki.github.io/DRBayes/reference/print.DRBayes.md)
shows, followed by which covariates Step 1 of Algorithm 4 selected,
which it left out, and what the sweep was allowed to move.

## Usage

``` r
# S3 method for class 'DRBayes_select'
print(x, digits = 3, ...)
```

## Arguments

- x:

  An object of class `"DRBayes_select"` from
  [`drbayes_select`](https://t-momozaki.github.io/DRBayes/reference/drbayes_select.md).

- digits:

  Number of significant digits. Default 3.

- ...:

  Ignored, present for consistency with the generic.

## Value

`x`, invisibly.

## See also

[`drbayes_select`](https://t-momozaki.github.io/DRBayes/reference/drbayes_select.md),
[`print.DRBayes`](https://t-momozaki.github.io/DRBayes/reference/print.DRBayes.md)

## Examples

``` r
set.seed(1)
n <- 250
dat <- data.frame(X1 = rnorm(n), X2 = rnorm(n), X3 = rnorm(n))
dat$A <- rbinom(n, 1, plogis(dat$X1 + 0.5 * dat$X2))
dat$Y <- 2 * dat$A + dat$X1 + dat$X2 + rnorm(n)
fit <- drbayes_select(Y ~ A + X1 + X2 + X3, A ~ X1 + X2 + X3, dat,
                      mc = 1600, bn = 400, seed = 2, verbose = FALSE,
                      control = drbayes_control(n_steps = 150))
print(fit)
#> Bayesian doubly robust causal inference via posterior coupling
#> 
#> Call:
#> drbayes_select(outcome.formula = Y ~ A + X1 + X2 + X3, ps.formula = A ~ 
#>     X1 + X2 + X3, data = dat, mc = 1600, bn = 400, control = drbayes_control(n_steps = 150), 
#>     verbose = FALSE)
#> 
#> Outcome model:    gaussian family, identity link
#> Propensity score: logit link
#> Observations:     250 (126 treated, 124 control)
#> Posterior draws:  2400
#> 
#> Average treatment effect
#>         mean  2.5% 97.5%
#> pc     2.056 1.794 2.320
#> g.comp 2.056 1.794 2.320
#> 
#> Report pc, the doubly robust estimand. g.comp is the untilted g-formula
#> posterior, which relies on the outcome model alone.
#> 
#> Both are posteriors for the average causal effect over the 250 covariate
#> vectors observed, which is what equation (3.7) averages. See ?drbayes_pc on
#> when that differs from the average over the population they were drawn from.
#> 
#> Posterior coupling: lambda = 0, moment condition met, |mean B_n| = 0.00268
#> below the tolerance 0.00277, after 0 of 150 steps.
#> 
#> Convergence diagnostics (worst over both models)
#>   Largest R-hat: 1.003 (outcome: X3)
#>   Smallest ESS:  1037 (propensity score: X2)
#> 
#> Confounder selection (Algorithm 4)
#>   Selected 3 of 3 candidate covariate(s) at |alpha_bar| >= 0.01
#>     Selected: X1, X2, X3.
#>   The moment condition was restricted to the selected covariates and the
#>   particles were rejuvenated as a whole. Every outcome coefficient is in the
#>   g-computation, so the estimand is unchanged.
```
