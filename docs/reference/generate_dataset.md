# Generate Synthetic Dataset for Causal Inference Simulation

Generates a synthetic dataset with confounders, treatment assignment,
and outcomes following a specific data generating process designed for
testing causal inference methods. Only the observed outcome is returned;
the potential outcomes are used to build it and then discarded.

## Usage

``` r
generate_dataset(nn = 200, pp = 40, seed = NULL)
```

## Arguments

- nn:

  Integer. Number of observations to generate. Default is 200.

- pp:

  Integer. Number of additional confounders to include beyond the base 4
  confounders. Default is 40. Can be 0 for only the 4 primary
  confounders.

- seed:

  Integer or NULL. Random seed for reproducibility. If NULL, no seed is
  set.

## Value

A data.frame with `nn` rows and `6 + pp` columns:

- AA:

  Binary treatment assignment (0 or 1)

- YY:

  Observed outcome (continuous)

- WW1, WW2, WW3, WW4:

  Primary confounding variables

- WW5, ..., WW(4+pp):

  Irrelevant covariates (only if pp \> 0)

## Details

The data generating process follows these steps:

**Covariates:**

- 4 confounders: \\X_1, X_2, X_3, X_4 \sim N(0, 1)\\

- `pp` further covariates \\X\_{4+j} \sim N(\mu_j, 1)\\ with \\\mu_j
  \sim \mathrm{Unif}(-2, 2)\\. These enter neither the treatment nor the
  outcome model, so they are irrelevant covariates rather than
  confounders, matching the high-dimensional setting of the paper's
  Appendix D.1. Note that the paper draws \\\mu_j \sim \mathrm{Unif}(-1,
  1)\\, so this generator is twice as dispersed as the published one.

**Treatment assignment:** \$\$P(A = 1 \| X) = \mathrm{logit}^{-1}(X_1 -
0.5 X_2 + 0.25 X_3 + 0.1 X_4)\$\$

There is no intercept, so treatment prevalence is one half.

**Potential outcomes:** \$\$Y^{(0)} = 100 + 27.4 X_1 + 13.7(X_2 + X_3 +
X_4) + \epsilon\$\$ \$\$Y^{(1)} = 210 + 27.4 X_1 + 13.7(X_2 + X_3 +
X_4) + \epsilon\$\$

where \\\epsilon \sim N(0, 1)\\. The SAME error appears in both
potential outcomes, as the paper's single \\\epsilon_i\\ implies, so the
unit level treatment effect is exactly 110 for every observation.

**Observed outcome:** \$\$Y = A \cdot Y^{(1)} + (1-A) \cdot Y^{(0)}\$\$

The true average treatment effect (ATE) is \\E\[Y^{(1)} - Y^{(0)}\] =
210 - 100 = 110\\.

## See also

[`run_parallel_simulation`](https://t-momozaki.github.io/DRBayes/reference/run_parallel_simulation.md),
[`drbayes_pc`](https://t-momozaki.github.io/DRBayes/reference/drbayes_pc.md)

## Examples

``` r
# Generate a dataset with only primary confounders
data_minimal <- generate_dataset(nn = 100, pp = 0, seed = 123)
head(data_minimal)
#>   AA       YY         WW1         WW2        WW3        WW4
#> 1  1 206.6657 -0.56047565 -0.71040656  2.1988103 -0.7152422
#> 2  1 215.9273 -0.23017749  0.25688371  1.3124130 -0.7526890
#> 3  1 233.2738  1.55870831 -0.24669188 -0.2651451 -0.9385387
#> 4  1 200.9081  0.07050839 -0.34754260  0.5431941 -1.0525133
#> 5  1 189.7569  0.12928774 -0.95161857 -0.4143399 -0.4371595
#> 6  0 141.7276  1.71506499 -0.04502772 -0.4762469  0.3311792
ncol(data_minimal)  # 6: AA, YY, WW1-WW4
#> [1] 6

# Generate a dataset with additional, unrelated covariates
data_extended <- generate_dataset(nn = 100, pp = 20, seed = 123)
ncol(data_extended)  # 26: AA, YY, WW1-WW24
#> [1] 26

# Check treatment assignment
table(data_minimal$AA)
#> 
#>  0  1 
#> 44 56 

# The unit level treatment effect is exactly 110, so the difference in
# means below is that plus the confounding the covariates carry.
aggregate(YY ~ AA, data = data_minimal, FUN = mean)
#>   AA        YY
#> 1  0  93.74135
#> 2  1 219.04242

# Estimate propensity scores
ps_model <- glm(AA ~ WW1 + WW2 + WW3 + WW4, data = data_minimal,
                family = binomial)
coef(ps_model)
#> (Intercept)         WW1         WW2         WW3         WW4 
#>   0.1297838   1.5863164  -0.6385581   0.1942240  -0.2034123 
```
