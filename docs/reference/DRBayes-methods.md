# Methods for posterior coupling fits

Printing, summarising and plotting for the objects of class `"DRBayes"`
returned by
[`drbayes_pc()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_pc.md).

## Details

Every method reports two posteriors for the average treatment effect,
and it matters which one is which.

`g.comp` is the untilted Bayesian g-formula posterior. It is what the
outcome model alone implies about the average treatment effect, obtained
by averaging the fitted contrast over the observed covariates. It uses
the propensity score model not at all, so it is consistent only if the
outcome model is correctly specified. In the notation of the paper it is
the tilted posterior at `lambda = 0`.

`pc` is the posterior coupling estimand of equation (3.8), and the
quantity to report. It is the same posterior tilted until the doubly
robust moment condition holds in posterior mean, which links the outcome
model draws to the propensity score model draws. It is consistent if
either model is correctly specified, which is the double robustness
property, and it is only doubly robust when the sweep actually met the
moment condition. Both posteriors average over the same empirical
distribution of the covariates, so their locations and their spreads are
directly comparable.

That empirical distribution is the third thing worth being clear about.
Both are posteriors for the average causal effect over the covariate
vectors that were observed, which is the average equation (3.7) takes,
and not for the average over the population those vectors came from. The
two differ whenever the fitted contrast varies from unit to unit;
[`drbayes_pc()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_pc.md)
says how much, and when.

## See also

[`drbayes_pc()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_pc.md)
for the fitting function,
[`rank_plot()`](https://t-momozaki.github.io/DRBayes/reference/rank_plot.md)
and
[`convergence_diagnostics()`](https://t-momozaki.github.io/DRBayes/reference/convergence_diagnostics.md)
for the underlying convergence tools.
