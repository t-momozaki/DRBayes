// Gaussian linear model, the Stan counterpart of bayes_lm().
//
// The priors are the ones that sampler documents, so the two backends target
// the same posterior:
//   theta   ~ N(0, prior_precision^-1)
//   sigma2  ~ InvGamma(sigma_a, sigma_b)
//   Y       ~ N(Z theta, sigma2 I)
// Z carries its intercept in column one, matching cbind(1, X) on the R side.
// The coefficient prior does not involve sigma2: a prior variance measured in
// units of sigma2 shrinks a coefficient by an amount that does not change with
// the units of Y, but it also adds the prior quadratic form to the residual
// sum of squares, which inflates the error variance whenever the coefficients
// are large next to it. bayes_lm() keeps the default weakly informative by
// scaling prior_precision to the data instead.

data {
  int<lower=0> N;
  int<lower=1> P;
  matrix[N, P] Z;
  vector[N] Y;
  matrix[P, P] prior_precision;
  real<lower=0> sigma_a;
  real<lower=0> sigma_b;
}

transformed data {
  // multi_normal_prec() would factorise the precision again at every gradient
  // evaluation. The prior is data, so its factor is computed once here.
  matrix[P, P] L_prior = cholesky_decompose(prior_precision);
}

parameters {
  vector[P] theta;
  real<lower=0> sigma2;
}

model {
  // theta' Q theta = || L' theta ||^2 with Q = L L'. The normalising constant
  // is free of the parameters, so it is dropped.
  target += -0.5 * dot_self(L_prior' * theta);
  sigma2 ~ inv_gamma(sigma_a, sigma_b);
  Y ~ normal_id_glm(Z, 0.0, theta, sqrt(sigma2));
}
