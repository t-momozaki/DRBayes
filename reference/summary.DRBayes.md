# Summarise a posterior coupling fit

Collects the posterior summaries of both estimands into a small data
frame, together with the sequential Monte Carlo diagnostics and the full
convergence diagnostics table, ready to be printed.

## Usage

``` r
# S3 method for class 'DRBayes'
summary(object, ...)

# S3 method for class 'summary.DRBayes'
print(x, digits = 3, ...)
```

## Arguments

- object:

  An object of class `"DRBayes"` from
  [`drbayes_pc()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_pc.md).

- ...:

  Ignored, present for consistency with the generic.

- x:

  An object of class `"summary.DRBayes"`.

- digits:

  Number of significant digits. The table of estimands is shown to the
  number of decimal places that gives the posterior standard deviation
  this many significant digits, so that a wide credible interval around
  a large effect does not round away. Defaults to 3.

## Value

An object of class `"summary.DRBayes"`, a list with components

- estimands:

  Data frame with one row per estimand, `pc` first, and columns `mean`,
  `sd`, `2.5\%`, `50\%` and `97.5\%`.

- call, family, link, ps.link, data_info, smc, diagnostics:

  Taken unchanged from `object`.

- thresholds:

  The R-hat and effective sample size thresholds the diagnostics are
  judged against.

- n_draws:

  Number of posterior draws behind each estimand.

`print.summary.DRBayes` returns `x` invisibly.

## Details

See
[DRBayes-methods](https://t-momozaki.github.io/DRBayes/reference/DRBayes-methods.md)
for what separates `pc` from `g.comp`. The quantiles are the usual
equal-tailed ones, so `2.5\%` and `97.5\%` bound a 95 percent credible
interval and `50\%` is the posterior median.

## See also

[`print.DRBayes()`](https://t-momozaki.github.io/DRBayes/reference/print.DRBayes.md),
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
summary(fit)
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
#> Outcome formula:  YY ~ AA + WW1 + WW2 + WW3 + WW4
#> PS formula:       AA ~ WW1 + WW2 + WW3 + WW4
#> 
#> Posterior summary of the average treatment effect
#>           mean    sd    2.5%     50%   97.5%
#> pc     109.771 0.170 109.434 109.766 110.089
#> g.comp 109.792 0.161 109.491 109.785 110.101
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
#>   Tilting parameter lambda: 1.35
#>   Posterior mean of B_n:    0.0029
#>   Tolerance:                0.00314
#>   Steps taken:              27 of 200
#>   Moment condition met:     yes
#> 
#> Convergence diagnostics
#>             model   parameter    mean     sd  rhat ess_bulk ess_tail mcse_mean
#>           outcome (Intercept) 100.106 0.1112 1.001     3038     2990   0.00203
#>           outcome          AA 109.792 0.1608 1.000     2942     2879   0.00296
#>           outcome         WW1  27.421 0.0844 1.000     3017     3110   0.00154
#>           outcome         WW2  13.577 0.0766 1.000     2999     2745   0.00140
#>           outcome         WW3  13.631 0.0706 1.000     3289     2712   0.00123
#>           outcome         WW4  13.695 0.0690 0.999     2807     2911   0.00130
#>  propensity score (Intercept)   0.106 0.1601 1.000     2955     2947   0.00295
#>  propensity score         WW1   0.765 0.1814 1.000     2519     2798   0.00361
#>  propensity score         WW2  -0.543 0.1653 1.001     2447     2584   0.00335
#>  propensity score         WW3   0.194 0.1475 1.001     2911     2718   0.00273
#>  propensity score         WW4   0.069 0.1472 1.000     2681     2826   0.00284
#> 
#> All parameters meet R-hat < 1.01 and ESS > 400.

# The table on its own, for example to put in a report.
summary(fit)$estimands
#>            mean        sd     2.5%      50%    97.5%
#> pc     109.7711 0.1704165 109.4340 109.7657 110.0889
#> g.comp 109.7924 0.1608122 109.4911 109.7854 110.1005
```
