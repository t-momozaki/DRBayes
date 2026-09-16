# Triangular sensitivity distributions

Builds the sensitivity distribution \\g()\\ that Section 6.3 of Orihara,
Momozaki and Sugasawa (2025) uses: a triangular density on
`(lower, upper)` peaking at `mode`. The result is a function of `n`
returning `n` draws, ready to hand to the `xi` argument of
[`drbayes_sensitivity()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_sensitivity.md).

## Usage

``` r
xi_triangular(lower, upper, mode = NULL)
```

## Arguments

- lower, upper:

  The support of the sensitivity parameter. Any bias outside
  `(lower, upper)` is being ruled out, so these are the assumption the
  whole analysis rests on.

- mode:

  Where the density peaks. Defaults to whichever of `lower` and `upper`
  is further from zero, the largest bias the support allows.

## Value

A function of `n` returning `n` draws, suitable as the `xi` argument of
[`drbayes_sensitivity()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_sensitivity.md).

## Details

The paper's two scenarios are the triangular distribution on \\(0,
0.5)\\ peaking at \\0.5\\ and the one on \\(-0.5, 0)\\ peaking at
\\-0.5\\, representing unmeasured confounding of a positive and of a
negative sign. In both the peak is at the endpoint further from zero,
which is the pessimistic reading of a bounded bias, and that endpoint is
what `mode` defaults to.

Those endpoints are on the scale of the outcome itself, because the
shift is on the conditional mean: for the binary outcome of Section 6.3
they are risk differences, and 0.5 is not a moderate bias there but a
very large one. That section reports an average treatment effect of
0.010 with 95 percent interval (-0.111, 0.130), a posterior standard
deviation of about 0.06, so \\E\[\xi\] = 1/3\\ is five of them. What a
sensitivity parameter has to be read against is that standard deviation,
and a bias several of them wide is one the original draws cannot
represent, whatever the units of the outcome; see the "Reading the
result" section of
[`drbayes_sensitivity()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_sensitivity.md).

## References

Orihara, S., Momozaki, T., & Sugasawa, S. (2025). Bayesian Doubly Robust
Causal Inference via Posterior Coupling. arXiv preprint
arXiv:2506.04868.

## See also

[`drbayes_sensitivity()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_sensitivity.md)

## Examples

``` r
positive <- xi_triangular(0, 0.5)
negative <- xi_triangular(-0.5, 0)
summary(positive(1000))
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.0194  0.2545  0.3598  0.3379  0.4320  0.4998 

# A bias of unknown sign, symmetric around zero and most likely absent.
either <- xi_triangular(-0.5, 0.5, mode = 0)
```
