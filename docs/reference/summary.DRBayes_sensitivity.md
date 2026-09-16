# Summarise a sensitivity analysis

Collects the posterior summaries of the reweighted and the original
average treatment effect into a small data frame, together with the
sensitivity distribution the reweighting assumed, the effective sample
size of the weights it produced, and the coupling that followed it.

## Usage

``` r
# S3 method for class 'DRBayes_sensitivity'
summary(object, ...)

# S3 method for class 'summary.DRBayes_sensitivity'
print(x, digits = 3, ...)
```

## Arguments

- object:

  An object of class `"DRBayes_sensitivity"`.

- ...:

  Ignored, present for consistency with the generic.

- x:

  An object of class `"summary.DRBayes_sensitivity"`.

- digits:

  Number of significant digits. The table of estimands is shown to the
  number of decimal places that gives the posterior standard deviation
  this many significant digits, so that a wide credible interval around
  a large effect does not round away. Defaults to 3.

## Value

An object of class `"summary.DRBayes_sensitivity"`, a list with
components

- estimands:

  Data frame with rows `sensitivity`, `reweighted` and `original` and
  columns `mean`, `sd`, `2.5%`, `50%` and `97.5%`. A quantile the
  effective sample size cannot support is `NA`.

- xi:

  Data frame of the mean, standard deviation and range of the draws from
  the sensitivity distribution.

- call, method, M, ess, xi_ess, n_draws, sigma, smc:

  Taken unchanged from `object`.

`print.summary.DRBayes_sensitivity` returns `x` invisibly.

## See also

[`print.DRBayes_sensitivity()`](https://t-momozaki.github.io/DRBayes/reference/print.DRBayes_sensitivity.md)

## Examples

``` r
# generate_dataset() leaves the caller's random number stream alone, so
# seed the coupling itself to make the fit below reproducible.
set.seed(1)
data <- generate_dataset(nn = 120, pp = 0, seed = 1)
fit  <- drbayes_pc(YY ~ AA + WW1 + WW2 + WW3 + WW4,
                   AA ~ WW1 + WW2 + WW3 + WW4,
                   data = data,
                   outcome.model = bayes_lm, ps.model = bayes_logit,
                   mc = 900, bn = 300, chains = 2L, verbose = FALSE,
                   control = drbayes_control(keep_particles = TRUE,
                                             n_steps = 200))
#> Warning: The posterior draws do not meet the convergence criteria (R-hat < 1.01, ESS > 400). outcome model: min ESS 355.7 (WW1). An effective sample size below 400 also makes R-hat itself unreliable, so a small R-hat here is not evidence of convergence. Raise mc, or thin less. Inspect $diagnostics.
summary(drbayes_sensitivity(fit, xi = xi_triangular(0, 0.1)))
#> Sensitivity analysis for unmeasured confounding
#> 
#> Call:
#> drbayes_sensitivity(object = fit, xi = xi_triangular(0, 0.1))
#> 
#> Sensitivity parameter: one xi shared by every observation (Algorithm A.1)
#> Draws from g():        100
#>                        integral resolved by a median of 98.2 of them, 60 at worst
#> Residual scale sigma:  drawn per posterior draw: mean 1.08, sd 0.0727
#> Posterior draws:       600
#> Effective sample size: 483 of 600 draws (80.5%)
#> Coupled again at:      lambda = 0.2
#> Moment condition:      |mean B_n| = 0.00233 against a tolerance of 0.0114, satisfied
#> Still doubly robust:   yes
#> 
#> Sensitivity parameter xi
#>    mean     sd     min    max
#>  0.0658 0.0245 0.00557 0.0998
#> 
#> Posterior summary of the average treatment effect
#>                mean    sd    2.5%     50%   97.5%
#> sensitivity 110.024 0.235 109.500 110.023 110.453
#> reweighted  110.038 0.216 109.601 110.034 110.449
#> original    110.041 0.227 109.587 110.036 110.462
#> 
#> The sensitivity row is the posterior under the assumed unmeasured confounding;
#> the reweighted row is the same draws before the moment condition was imposed on
#> them again, and holds the whole of the shift; the original row is the fit that
#> was reweighted. See the "Reweight first, then couple" section of
#> ?drbayes_sensitivity for why the first two differ.
```
