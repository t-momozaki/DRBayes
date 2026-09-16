# Bayesian Bootstrap Doubly Robust Estimation of the Average Treatment Effect

Implements the Bayesian bootstrap doubly robust estimator of Saarela et
al. (2016), which the posterior coupling paper uses as its comparator.
Each Dirichlet draw reweights the sample, both models are refitted under
those weights, and the doubly robust contrast is evaluated. The
resulting draws summarise posterior uncertainty about the estimand.

## Usage

``` r
drbayes_bb(
  outcome.formula,
  ps.formula,
  data,
  num_iterations = 3000,
  family = NULL,
  link = NULL,
  ps.link = "logit",
  trim = 0,
  na.action = "na.omit",
  verbose = TRUE
)
```

## Arguments

- outcome.formula:

  Formula for the outcome model, including the treatment and any
  interactions with it, for example `Y ~ A + X1 + A:X1`.

- ps.formula:

  Formula for the propensity score model, with the treatment on the left
  and the confounders on the right, for example `A ~ X1 + X2`. The
  treatment must not appear on the right.

- data:

  Data frame containing every variable used by either formula.

- num_iterations:

  Integer. Number of Dirichlet draws. Default is 3000. These are
  independent draws, not a Markov chain: there is no burn-in, no
  autocorrelation and nothing to thin.

- family:

  Character string, `"gaussian"` or `"binomial"`. If NULL (default) it
  is detected from the outcome.

- link:

  Link function of the outcome model: `"identity"` for gaussian,
  `"logit"` or `"probit"` for binomial. If NULL (default) it is taken
  from the family.

- ps.link:

  Link function of the propensity score model, `"logit"` (default) or
  `"probit"`.

- trim:

  Numeric in \[0, 0.5). Fitted propensity scores are confined to
  `[trim, 1 - trim]`. The default 0 leaves them untouched; a small value
  such as 0.01 bounds the augmentation term when the treatment groups
  overlap poorly. Truncation trades variance for bias, so the count of
  truncated values is reported.

- na.action:

  How to handle missing values: `"na.omit"` (default), `"na.fail"` or
  `"na.exclude"`. Only observations complete in both models are used.

- verbose:

  Logical. Whether to report progress through
  [`message`](https://rdrr.io/r/base/message.html). Default is TRUE.

## Value

A numeric vector of length `num_iterations` holding one draw of the
estimand per Dirichlet weight vector.

## Details

Both regressions are solved on design matrices built once from the
formulas, rather than refitted through the formula interface on every
draw. The Gaussian outcome model has a linear score, so it is solved in
closed form by a pivoted QR; the binomial outcome model and the
propensity score model are solved by iteratively reweighted least
squares. Rank-deficient designs and non-convergent fits raise an error
rather than propagating NA or an arbitrary root into the returned draws.

The weights are Dirichlet(1, ..., 1), drawn as independent unit
exponentials divided by their sum, which is the same distribution.

## Which estimand this targets

This is not the superpopulation average treatment effect. As Orihara,
Momozaki and Sugasawa (2025, section 4.3) put it, Saarela's target is
defined over the population induced by the estimated propensity score,
"a mixed ATE, whereas the latter corresponds to the (super)population
ATE", and "in general, this estimand is not consistent with the ATE".
The two coincide when the propensity score model is correctly specified.
Compare the output against the true ATE only with that caveat in mind,
and prefer
[`drbayes_pc`](https://t-momozaki.github.io/DRBayes/reference/drbayes_pc.md)
when the superpopulation ATE is what is wanted.

## References

Saarela, O., Belzile, L. R., & Stephens, D. A. (2016). A Bayesian view
of doubly robust causal inference. Biometrika, 103(3), 667-681.

Orihara, S., Momozaki, T., & Sugasawa, S. (2025). Bayesian Doubly Robust
Causal Inference via Posterior Coupling. arXiv preprint
arXiv:2506.04868.

## See also

[`drbayes_pc`](https://t-momozaki.github.io/DRBayes/reference/drbayes_pc.md)
for the posterior coupling estimator.

## Examples

``` r
set.seed(123)
n <- 300
d <- data.frame(X1 = rnorm(n), X2 = rnorm(n), X3 = rnorm(n))
d$A <- rbinom(n, 1, plogis(0.3 + 0.5 * d$X1 - 0.4 * d$X2))
d$Y <- 1 + 2 * d$A + d$X1 + d$X2 + 0.5 * d$A * d$X1 + rnorm(n)

# Interactions with the treatment are written in the formula, so the
# counterfactual predictions recompute them.
draws <- drbayes_bb(
  outcome.formula = Y ~ A + X1 + X2 + X3 + A:X1,
  ps.formula      = A ~ X1 + X2 + X3,
  data            = d,
  num_iterations  = 500,
  verbose         = FALSE
)

# The true average treatment effect is 2.
mean(draws)
#> [1] 2.083922
quantile(draws, c(0.025, 0.975))
#>     2.5%    97.5% 
#> 1.835853 2.335726 
```
