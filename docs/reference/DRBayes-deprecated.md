# Functions renamed in DRBayes 0.1.0

The exported functions were renamed before the first CRAN release, to
one naming scheme and away from dots. A dot in a function name is how R
marks an S3 method, so `print.DRBayes` and a function called
`DRBayes.PC` would have sat in the same namespace with no way to tell
them apart.

Each old name still works and forwards to its replacement, with a
warning. They will be removed in the next release.

|  |  |
|----|----|
| **Old** | **New** |
| `DRBayes.PC()` | [`drbayes_pc()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_pc.md) |
| `DRBayes.BB()` | [`drbayes_bb()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_bb.md) |
| `DRBayes.BB.nleqslv()` | [`drbayes_bb()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_bb.md) |
| `B.LM()` | [`bayes_lm()`](https://t-momozaki.github.io/DRBayes/reference/bayes_lm.md) |
| `B.Logit()` | [`bayes_logit()`](https://t-momozaki.github.io/DRBayes/reference/bayes_logit.md) |
| `B.Probit()` | [`bayes_probit()`](https://t-momozaki.github.io/DRBayes/reference/bayes_probit.md) |
| `HS.LM()` | [`bayes_lm_hs()`](https://t-momozaki.github.io/DRBayes/reference/bayes_lm_hs.md) |
| `HS.Logit()` | [`bayes_logit_hs()`](https://t-momozaki.github.io/DRBayes/reference/bayes_logit_hs.md) |
| `HS.Probit()` | [`bayes_probit_hs()`](https://t-momozaki.github.io/DRBayes/reference/bayes_probit_hs.md) |

`DRBayes.BB.nleqslv()` solved the same two weighted regressions with a
general nonlinear root finder. It computed the same estimator as
`DRBayes.BB()` to about 1e-08, the difference being only the root
finder's stopping tolerance, so it forwards there rather than being
reimplemented.

## Usage

``` r
DRBayes.PC(...)

DRBayes.BB(...)

DRBayes.BB.nleqslv(...)

B.LM(...)

B.Logit(...)

B.Probit(...)

HS.LM(...)

HS.Logit(...)

HS.Probit(...)
```

## Arguments

- ...:

  Passed to the replacement function.

## Value

Whatever the replacement function returns.
