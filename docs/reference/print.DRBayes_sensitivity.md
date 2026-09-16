# Print a sensitivity analysis

A one screen summary of a
[`drbayes_sensitivity()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_sensitivity.md)
analysis: the call, the sensitivity distribution, the effective sample
size of the weights, whether the analysis is still doubly robust, and
the average treatment effect before and after.

## Usage

``` r
# S3 method for class 'DRBayes_sensitivity'
print(x, digits = 3, ...)
```

## Arguments

- x:

  An object of class `"DRBayes_sensitivity"`.

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

The diagnostics come before the estimates on purpose. A reweighted
posterior is only as good as the number of draws that carry weight, and
it is doubly robust only if the coupling that follows the reweighting
converged on a reweighting that left it draws to work with. Those are
two lines of the output: the moment condition, and the verdict on the
analysis as a whole.

The credible interval is left out where the effective sample size cannot
support it. A 2.5 percent point with fewer than ten effective draws
beyond it is the most extreme handful of them, not a quantile, which
takes 400 effective draws.

## See also

[`summary.DRBayes_sensitivity()`](https://t-momozaki.github.io/DRBayes/reference/summary.DRBayes_sensitivity.md)

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
drbayes_sensitivity(fit, xi = xi_triangular(0, 0.1))
#> Sensitivity analysis for unmeasured confounding
#> 
#> Call:
#> drbayes_sensitivity(object = fit, xi = xi_triangular(0, 0.1))
#> 
#> Sensitivity parameter: one xi shared by every observation (Algorithm A.1)
#>                        100 draws from g(): mean 0.0658, range 0.00557 to 0.0998
#>                        integral resolved by a median of 98.2 of them, 60 at worst
#> Residual scale sigma:  drawn per posterior draw: mean 1.08, sd 0.0727
#> Posterior draws:       600
#> Effective sample size: 483 of 600 draws (80.5%)
#> Coupled again at:      lambda = 0.2
#> Moment condition:      |mean B_n| = 0.00233 against a tolerance of 0.0114, satisfied
#> Still doubly robust:   yes
#> 
#> Average treatment effect
#>                mean    2.5%   97.5%
#> sensitivity 110.024 109.500 110.453
#> reweighted  110.038 109.601 110.449
#> original    110.041 109.587 110.462
#> 
#> The sensitivity row is the posterior under the assumed unmeasured confounding;
#> the reweighted row is the same draws before the moment condition was imposed on
#> them again, and holds the whole of the shift; the original row is the fit that
#> was reweighted. See the "Reweight first, then couple" section of
#> ?drbayes_sensitivity for why the first two differ.
```
