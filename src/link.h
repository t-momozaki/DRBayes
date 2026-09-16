#ifndef DRBAYES_LINK_H
#define DRBAYES_LINK_H

#include <Rmath.h>
#include <cmath>
#include <cstddef>

namespace drbayes {

// The inverse links the kernels recognise. R decides which of these a supplied
// closure is, by the same identical() probe the package has always used, and
// passes the answer down as one of these codes; anything it does not recognise
// never reaches C++ at all.
enum LinkCode {
  LINK_IDENTITY = 0,
  LINK_LOGIT    = 1,
  LINK_PROBIT   = 2
};


// The inverse link, taken from Rmath rather than written out.
//
// 1/(1+exp(-x)) and 0.5*erfc(-x/sqrt(2)) are the same functions mathematically
// and differ from plogis and pnorm in the last bits, worst in the tails, which
// is exactly where the weights are largest. R::plogis and R::pnorm are the
// routines stats::plogis and stats::pnorm call, so using them makes the values
// R's rather than merely close to R's.
inline double inverse_link(double eta, int code) {
  switch (code) {
    case LINK_LOGIT:  return R::plogis(eta, 0.0, 1.0, 1, 0);
    case LINK_PROBIT: return R::pnorm(eta, 0.0, 1.0, 1, 0);
    default:          return eta;
  }
}

inline void apply_inverse_link(double* x, std::size_t n, int code) {
  if (code == LINK_IDENTITY) return;
  for (std::size_t i = 0; i < n; ++i) x[i] = inverse_link(x[i], code);
}


// The inverse probability weight (A - e) / (e (1 - e)) of equation (3.4).
//
// Under a logistic link it collapses exactly. With s = 2A - 1, at A = 1 the
// weight is 1/e = 1 + exp(-eta), and at A = 0 it is -1/(1-e) = -(1 + exp(eta)),
// which is s (1 + exp(-s eta)) in both cases. One exponential replaces a
// plogis, a subtraction, a multiplication and a division, and the denominator
// that could reach zero is gone. The collapse is specific to the logistic link;
// under a probit one the general form is the only form.
// The weight at one observation, under each of the two forms. A caller that
// wants the weights in memory and one that consumes each as it is formed both
// go through these, so the two cannot drift apart. The link is chosen by the
// caller rather than tested here, so that a loop over a column of observations
// carries no branch.
inline double logit_weight_at(double eta, double a) {
  const double s = 2.0 * a - 1.0;
  return s * (1.0 + std::exp(-s * eta));
}

// The rounding order here is R's: the difference, then one minus the score,
// then their product, then the quotient.
inline double general_weight_at(double eta, double a, int ps_link) {
  const double p = inverse_link(eta, ps_link);
  const double numerator = a - p;
  const double complement = 1.0 - p;
  const double denominator = p * complement;
  return numerator / denominator;
}

inline void ipw_weight_logit(double* eta, std::size_t n, const double* a) {
  for (std::size_t i = 0; i < n; ++i) eta[i] = logit_weight_at(eta[i], a[i]);
}

inline void ipw_weight_general(double* eta, std::size_t n, const double* a,
                               int ps_link) {
  for (std::size_t i = 0; i < n; ++i) {
    eta[i] = general_weight_at(eta[i], a[i], ps_link);
  }
}

// Overwrites a block of linear predictors with the weights they imply. The
// weights land in the same buffer because nothing downstream needs the linear
// predictors again, and a second buffer of that size is the allocation this
// design exists to avoid.
inline void ipw_weight(double* eta, std::size_t n, const double* a,
                       int ps_link, std::size_t columns) {
  for (std::size_t j = 0; j < columns; ++j) {
    double* column = eta + j * n;
    if (ps_link == LINK_LOGIT) {
      ipw_weight_logit(column, n, a);
    } else {
      ipw_weight_general(column, n, a, ps_link);
    }
  }
}

}  // namespace drbayes

#endif  // DRBAYES_LINK_H
