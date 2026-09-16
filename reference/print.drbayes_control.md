# Print Posterior Coupling Control Settings

Print Posterior Coupling Control Settings

## Usage

``` r
# S3 method for class 'drbayes_control'
print(x, ...)
```

## Arguments

- x:

  An object of class `"drbayes_control"`, as returned by
  [`drbayes_control`](https://t-momozaki.github.io/DRBayes/reference/drbayes_control.md).

- ...:

  Ignored, present for compatibility with the generic.

## Value

`x`, invisibly. Called for the printed output.

## Examples

``` r
print(drbayes_control(n_steps = 200))
#> Control settings for drbayes_pc()
#> 
#>   MCMC         mc = 5000, bn = 1000, thin = 2, chains = 4, init = automatic
#>   Tilting      lambda_max = 10, n_steps = 200 (step 0.05), smoothing = 0.99
#>                tol = Monte Carlo standard error, ridge = 1e-05, ess_frac = 0.1
#>                moment = ipw, pruning = off
#>                newton_steps = 100 (importance sampling only)
#>   Diagnostics  rhat_max = 1.01, ess_min = 400
#>   Output       keep_draws = FALSE, keep_particles = FALSE
```
