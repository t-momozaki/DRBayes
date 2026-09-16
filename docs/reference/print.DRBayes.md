# Print a posterior coupling fit

A one screen summary of a
[`drbayes_pc()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_pc.md)
fit: the call, the two models, the data, the posterior mean and a 95
percent credible interval for both estimands, the state of the
sequential Monte Carlo sweep, and the worst convergence diagnostics.

## Usage

``` r
# S3 method for class 'DRBayes'
print(x, digits = 3, ...)
```

## Arguments

- x:

  An object of class `"DRBayes"` from
  [`drbayes_pc()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_pc.md).

- digits:

  Number of significant digits. The table of estimands is shown to the
  number of decimal places that gives the posterior standard deviation
  this many significant digits, so that a wide credible interval around
  a large effect does not round away. Defaults to 3.

- ...:

  Ignored, present for consistency with the generic.

## Value

`x`, invisibly.

## Details

See
[DRBayes-methods](https://t-momozaki.github.io/DRBayes/reference/DRBayes-methods.md)
for what separates `pc` from `g.comp`. Two things in the output are
warnings rather than description. If the sweep did not meet the moment
condition, the `pc` draws are not doubly robust and the message says so;
treat them as no better than `g.comp`. If a convergence threshold is
breached, the offending R-hat or effective sample size is marked, and
the full table is available through
[`summary.DRBayes()`](https://t-momozaki.github.io/DRBayes/reference/summary.DRBayes.md).

## See also

[`summary.DRBayes()`](https://t-momozaki.github.io/DRBayes/reference/summary.DRBayes.md),
[`plot.DRBayes()`](https://t-momozaki.github.io/DRBayes/reference/plot.DRBayes.md)

## Examples

``` r
# generate_dataset() leaves the caller's random number stream alone, so
# seed the coupling itself to make the fit below reproducible.
set.seed(1)
data <- generate_dataset(nn = 200, pp = 0, seed = 1)
fit  <- drbayes_pc(YY ~ AA + WW1 + WW2 + WW3 + WW4,
                   AA ~ WW1 + WW2 + WW3 + WW4,
                   data = data,
                   outcome.model = bayes_lm, ps.model = bayes_logit,
                   mc = 2000, bn = 500, verbose = FALSE,
                   control = drbayes_control(n_steps = 200))
fit
#> Bayesian doubly robust causal inference via posterior coupling
#> 
#> Call:
#> drbayes_pc(outcome.formula = YY ~ AA + WW1 + WW2 + WW3 + WW4, 
#>     ps.formula = AA ~ WW1 + WW2 + WW3 + WW4, data = data, mc = 2000, 
#>     bn = 500, control = drbayes_control(n_steps = 200), outcome.model = bayes_lm, 
#>     ps.model = bayes_logit, verbose = FALSE)
#> 
#> Outcome model:    gaussian family, identity link
#> Propensity score: logit link
#> Observations:     200 (104 treated, 96 control)
#> Posterior draws:  3000
#> 
#> Average treatment effect
#>           mean    2.5%   97.5%
#> pc     109.771 109.434 110.089
#> g.comp 109.792 109.491 110.101
#> 
#> Report pc, the doubly robust estimand. g.comp is the untilted g-formula
#> posterior, which relies on the outcome model alone.
#> 
#> Both are posteriors for the average causal effect over the 200 covariate
#> vectors observed, which is what equation (3.7) averages. See ?drbayes_pc on
#> when that differs from the average over the population they were drawn from.
#> 
#> Posterior coupling: lambda = 1.35, moment condition met, |mean B_n| = 0.0029
#> below the tolerance 0.00314, after 27 of 200 steps.
#> 
#> Convergence diagnostics (worst over both models)
#>   Largest R-hat: 1.001 (outcome: (Intercept))
#>   Smallest ESS:  2447 (propensity score: WW2)
```
