// Probit regression, the Stan counterpart of bayes_probit().
//
// The prior is the one that sampler documents, so the two backends target the
// same posterior:
//   theta ~ N(0, prior_precision^-1)
//   Y     ~ Bernoulli(Phi(Z theta))
// Z carries its intercept in column one, matching cbind(1, X) on the R side.

data {
  int<lower=0> N;
  int<lower=1> P;
  matrix[N, P] Z;
  array[N] int<lower=0, upper=1> Y;
  matrix[P, P] prior_precision;
}

transformed data {
  // The prior is data, so its factor is computed once rather than at every
  // gradient evaluation, which is what multi_normal_prec() would cost.
  matrix[P, P] L_prior = cholesky_decompose(prior_precision);

  // Writing the likelihood as Phi(sign * eta) turns both outcomes into one
  // vectorised call and keeps the log density accurate in the far tails, where
  // log(1 - Phi(eta)) computed directly underflows to -Inf.
  vector[N] y_sign;
  for (n in 1 : N) {
    y_sign[n] = Y[n] == 1 ? 1.0 : -1.0;
  }
}

parameters {
  vector[P] theta;
}

model {
  // theta' Q theta = || L' theta ||^2 with Q = L L'.
  target += -0.5 * dot_self(L_prior' * theta);
  target += std_normal_lcdf(y_sign .* (Z * theta) | );
}
