// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "link.h"
#include "reduce.h"
#include <algorithm>

// Reports the BLAS-backed matrix product Armadillo will use, so that a build
// which silently failed to link against it fails a test rather than running
// slowly. The kernels that follow all rest on that product being R's own.
// [[Rcpp::export(rng = false)]]
double drb_blas_check(const arma::mat& a, const arma::mat& b) {
  return arma::accu(a * b);
}

// The two reductions and the range tracker, reachable from R so that their
// agreement with colMeans, mean, min and max is a test rather than a claim.
//
// wide says whether R accumulates in long double on this build. It is passed
// rather than tested for here; see the note at the top of reduce.h.
// [[Rcpp::export(rng = false)]]
double drb_mean_one_pass(const Rcpp::NumericVector& x, bool wide) {
  const std::size_t n = static_cast<std::size_t>(x.size());
  return wide ? drbayes::mean_one_pass<long double>(x.begin(), n)
              : drbayes::mean_one_pass<double>(x.begin(), n);
}

// [[Rcpp::export(rng = false)]]
double drb_mean_two_pass(const Rcpp::NumericVector& x, bool wide) {
  const std::size_t n = static_cast<std::size_t>(x.size());
  return wide ? drbayes::mean_two_pass<long double>(x.begin(), n)
              : drbayes::mean_two_pass<double>(x.begin(), n);
}

namespace {

template <typename Acc>
Rcpp::NumericVector row_means_impl(const Rcpp::NumericMatrix& x, int block) {
  const std::size_t rows = static_cast<std::size_t>(x.nrow());
  const std::size_t cols = static_cast<std::size_t>(x.ncol());
  drbayes::RowMeanAccumulator<Acc> acc(rows);
  for (std::size_t lo = 0; lo < cols; lo += static_cast<std::size_t>(block)) {
    const std::size_t take = std::min(static_cast<std::size_t>(block), cols - lo);
    acc.add_block(x.begin() + lo * rows, rows, take);
  }
  Rcpp::NumericVector out(rows);
  acc.write_means(out.begin());
  return out;
}

}  // namespace

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector drb_row_means(const Rcpp::NumericMatrix& x, int block,
                                  bool wide) {
  return wide ? row_means_impl<long double>(x, block)
              : row_means_impl<double>(x, block);
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector drb_range(const Rcpp::NumericVector& x) {
  drbayes::RangeTracker tracker;
  tracker.observe(x.begin(), static_cast<std::size_t>(x.size()));
  if (tracker.empty()) return Rcpp::NumericVector(0);
  return Rcpp::NumericVector::create(tracker.low(), tracker.high());
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector drb_inverse_link(Rcpp::NumericVector x, int code) {
  Rcpp::NumericVector out = Rcpp::clone(x);
  drbayes::apply_inverse_link(out.begin(), static_cast<std::size_t>(out.size()),
                              code);
  return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector drb_ipw_weight(Rcpp::NumericVector eta,
                                   const Rcpp::NumericVector& a, int ps_link) {
  // The loop length and the buffer came from different arguments, so a
  // treatment vector longer than the linear predictors wrote past the end of
  // the buffer. The kernels always pass matching lengths; this entry point
  // exists for the tests, and is where a mismatch would first arrive.
  if (a.size() != eta.size()) {
    Rcpp::stop("the treatment and the linear predictors must be the same "
               "length");
  }
  Rcpp::NumericVector out = Rcpp::clone(eta);
  drbayes::ipw_weight(out.begin(), static_cast<std::size_t>(eta.size()),
                      a.begin(), ps_link, 1);
  return out;
}


namespace {

// The block partition draw_blocks() uses, so that the matrix products here are
// the ones R made block by block.
//
// Matching it is required, not merely tidy. The reduction is indifferent to
// the grouping, since each column is reduced over its own entries in their own
// order however the draws are split; the product is not. A BLAS may reach a
// different kernel for a differently shaped product, and on a reference BLAS
// the same tcrossprod taken whole and taken in halves differs in the last bits
// on about one entry in a hundred. Grouping the draws as R groups them is what
// makes the two agree bit for bit rather than nearly.
inline arma::uword block_size(arma::uword n, double max_cells) {
  const double per_block = max_cells / static_cast<double>(n > 0 ? n : 1);
  const arma::uword size = static_cast<arma::uword>(per_block);
  return size > 0 ? size : 1;
}

// Wraps R-owned memory without copying. The matrices here run to tens of
// megabytes at the scale the paper works at, and Rcpp::as<arma::mat> would
// duplicate every one of them. Nothing writes through these views.
inline arma::mat view(const Rcpp::NumericMatrix& m) {
  return arma::mat(const_cast<double*>(m.begin()), m.nrow(), m.ncol(),
                   false, true);
}

// A block of a working buffer, sized to the draws this pass covers.
//
// The product is assigned to one of these rather than to a submatrix view of a
// longer-lived matrix, because assigning to a submatrix builds the product in
// a fresh matrix first and copies it in. Writing through a header over the
// buffer lets the BLAS write there directly, and at the shapes the sweep uses
// that is about a third off the cost of the product.
inline arma::mat block(double* data, arma::uword rows, arma::uword cols) {
  return arma::mat(data, rows, cols, false, true);
}

}  // namespace


namespace {

template <typename Acc>
Rcpp::List moment_ipw_impl(const Rcpp::NumericMatrix& Z_ps,
                           const Rcpp::NumericMatrix& Z_lm,
                           const Rcpp::NumericVector& A,
                           const Rcpp::NumericVector& Y,
                           const Rcpp::NumericMatrix& betas_ps,
                           const Rcpp::NumericMatrix& betas_otc,
                           int ps_link, int otc_link, double max_cells) {
  const arma::uword n = static_cast<arma::uword>(Y.size());
  const arma::uword S = static_cast<arma::uword>(betas_ps.nrow());

  if (static_cast<arma::uword>(Z_ps.nrow()) != n ||
      static_cast<arma::uword>(Z_lm.nrow()) != n ||
      static_cast<arma::uword>(A.size()) != n) {
    Rcpp::stop("Z.ps, Z.lm, A and Y must agree on the number of observations");
  }
  if (static_cast<arma::uword>(betas_otc.nrow()) != S) {
    Rcpp::stop("betas.ps and betas.otc must have the same number of draws");
  }
  if (static_cast<arma::uword>(betas_ps.ncol()) !=
        static_cast<arma::uword>(Z_ps.ncol()) ||
      static_cast<arma::uword>(betas_otc.ncol()) !=
        static_cast<arma::uword>(Z_lm.ncol())) {
    Rcpp::stop("each set of draws must have one column per design matrix column");
  }

  Rcpp::NumericVector BB(S);
  drbayes::RangeTracker range;

  if (S > 0 && n > 0) {
    const arma::mat Zps = view(Z_ps);
    const arma::mat Zlm = view(Z_lm);
    const arma::mat Bps = view(betas_ps);
    const arma::mat Botc = view(betas_otc);
    const double* a = A.begin();
    const double* y = Y.begin();

    const arma::uword step = block_size(n, max_cells);
    const arma::uword widest = std::min(step, S);
    std::vector<double> eta_buffer(static_cast<std::size_t>(n) * widest);
    std::vector<double> mu_buffer(static_cast<std::size_t>(n) * widest);
    double* eta_data = eta_buffer.data();
    double* mu_data  = mu_buffer.data();

    for (arma::uword lo = 0; lo < S; lo += step) {
      const arma::uword hi = std::min(lo + step, S) - 1;
      const arma::uword width = hi - lo + 1;

      // A * B.t() reaches dgemm with the transpose flag set rather than
      // building the transpose, which is the product tcrossprod() makes.
      arma::mat eta = block(eta_data, n, width);
      arma::mat mu  = block(mu_data, n, width);
      eta = Zps * Bps.rows(lo, hi).t();
      mu  = Zlm * Botc.rows(lo, hi).t();

      range.observe(eta_data, n * width);
      drbayes::apply_inverse_link(mu_data, n * width, otc_link);

      // Each weight is formed where it is used. Writing the whole column of
      // them back over the linear predictors first, as an earlier form did,
      // costs a pass that stores values nothing reads again. The link is
      // settled before the loop over observations, not inside it.
      //
      // The product is named before it is accumulated, and that is not
      // decoration. Left as one expression the compiler contracts it into a
      // fused multiply-add, which keeps the product at full width instead of
      // rounding it to double. R cannot: ww * (Y - mu) is a double matrix
      // before colMeans ever sees it. Written this way the two agreed at every
      // sample size and both links; written as one expression they disagreed
      // in the last bits at every sample size under a probit link, and under a
      // logistic one wherever the sample size was not a multiple of the vector
      // width.
      for (arma::uword j = 0; j < width; ++j) {
        const double* e = eta_data + j * n;
        const double* m = mu_data + j * n;
        Acc sum = 0.0;
        if (ps_link == drbayes::LINK_LOGIT) {
          for (arma::uword i = 0; i < n; ++i) {
            const double term = drbayes::logit_weight_at(e[i], a[i]) *
                                (y[i] - m[i]);
            sum += term;
          }
        } else {
          for (arma::uword i = 0; i < n; ++i) {
            const double term = drbayes::general_weight_at(e[i], a[i], ps_link) *
                                (y[i] - m[i]);
            sum += term;
          }
        }
        BB[lo + j] = static_cast<double>(sum / static_cast<Acc>(n));
      }
    }
  }

  Rcpp::NumericVector eta_range = range.empty()
    ? Rcpp::NumericVector(0)
    : Rcpp::NumericVector::create(range.low(), range.high());

  return Rcpp::List::create(Rcpp::Named("BB") = BB,
                            Rcpp::Named("eta_range") = eta_range);
}

}  // namespace

//' The moment condition of equation (3.4)
//'
//' Two matrix products per block go to the same BLAS R would have used, and
//' everything after them is done in place: the fitted means are overwritten by
//' the inverse link, each weight is formed at the observation it multiplies,
//' and the residual, the product and the column sum never exist as matrices at
//' all. That is where the time goes -- the arithmetic is the same arithmetic,
//' but R forms about seven draws-by-observations matrices per block where this
//' holds two buffers.
//'
//' The positivity guard is not applied here. The extremes of the linear
//' predictor are returned instead, and R raises the error, so that its wording
//' lives in one place and no C++ exception has to travel through it.
//'
//' @return A list with `BB`, one value per draw, and `eta_range`, the smallest
//'   and largest linear predictor seen, or a zero-length vector when there were
//'   no draws to see one in.
//'
//' The accumulator width is settled once, on entry, from what R reports about
//' its own build; see the note at the top of reduce.h.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export(rng = false)]]
Rcpp::List drb_moment_ipw(const Rcpp::NumericMatrix& Z_ps,
                          const Rcpp::NumericMatrix& Z_lm,
                          const Rcpp::NumericVector& A,
                          const Rcpp::NumericVector& Y,
                          const Rcpp::NumericMatrix& betas_ps,
                          const Rcpp::NumericMatrix& betas_otc,
                          int ps_link, int otc_link, double max_cells,
                          bool wide) {
  return wide
    ? moment_ipw_impl<long double>(Z_ps, Z_lm, A, Y, betas_ps, betas_otc,
                                   ps_link, otc_link, max_cells)
    : moment_ipw_impl<double>(Z_ps, Z_lm, A, Y, betas_ps, betas_otc,
                              ps_link, otc_link, max_cells);
}


namespace {

template <typename Acc>
Rcpp::NumericVector moment_subclass_impl(const Rcpp::NumericMatrix& Z_lm,
                                         const Rcpp::NumericVector& Y,
                                         const Rcpp::NumericVector& w,
                                         const Rcpp::NumericMatrix& betas_otc,
                                         int otc_link, double max_cells) {
  const arma::uword n = static_cast<arma::uword>(Y.size());
  const arma::uword S = static_cast<arma::uword>(betas_otc.nrow());

  if (static_cast<arma::uword>(Z_lm.nrow()) != n ||
      static_cast<arma::uword>(w.size()) != n) {
    Rcpp::stop("Z.lm, the stratum weights and Y must agree on the number of "
               "observations");
  }
  if (static_cast<arma::uword>(betas_otc.ncol()) !=
      static_cast<arma::uword>(Z_lm.ncol())) {
    Rcpp::stop("the outcome draws must have one column per design matrix column");
  }

  Rcpp::NumericVector BB(S);
  if (S == 0 || n == 0) return BB;

  const arma::mat Zlm = view(Z_lm);
  const arma::mat Botc = view(betas_otc);
  const double* y = Y.begin();
  const double* weight = w.begin();

  const arma::uword step = block_size(n, max_cells);
  std::vector<double> mu_buffer(static_cast<std::size_t>(n) * std::min(step, S));
  double* mu_data = mu_buffer.data();

  for (arma::uword lo = 0; lo < S; lo += step) {
    const arma::uword hi = std::min(lo + step, S) - 1;
    const arma::uword width = hi - lo + 1;

    arma::mat mu = block(mu_data, n, width);
    mu = Zlm * Botc.rows(lo, hi).t();
    drbayes::apply_inverse_link(mu_data, n * width, otc_link);

    // Four columns at once. A floating point sum is a chain, each addition
    // waiting on the one before it, and this loop has a multiply and a
    // subtraction to do while it waits, which is not enough. Four independent
    // chains give the adder something to interleave, and since each column
    // still accumulates its own terms in its own order the values do not
    // change. Measured against the plain loop, two and a half times faster at
    // n = 500 and seven at n = 8000, identical bit for bit at both.
    //
    // The same rewrite does nothing for equation (3.4), where an exponential
    // per observation already hides the latency: 0.96 to 1.12 times, which is
    // no change. It is not worth the same shape there.
    const Acc divisor = static_cast<Acc>(n);
    arma::uword j = 0;
    for (; j + 4 <= width; j += 4) {
      const double* m0 = mu_data + (j + 0) * n;
      const double* m1 = mu_data + (j + 1) * n;
      const double* m2 = mu_data + (j + 2) * n;
      const double* m3 = mu_data + (j + 3) * n;
      Acc s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
      for (arma::uword i = 0; i < n; ++i) {
        const double residual_weight = weight[i];
        const double outcome = y[i];
        const double t0 = residual_weight * (outcome - m0[i]); s0 += t0;
        const double t1 = residual_weight * (outcome - m1[i]); s1 += t1;
        const double t2 = residual_weight * (outcome - m2[i]); s2 += t2;
        const double t3 = residual_weight * (outcome - m3[i]); s3 += t3;
      }
      BB[lo + j + 0] = static_cast<double>(s0 / divisor);
      BB[lo + j + 1] = static_cast<double>(s1 / divisor);
      BB[lo + j + 2] = static_cast<double>(s2 / divisor);
      BB[lo + j + 3] = static_cast<double>(s3 / divisor);
    }
    for (; j < width; ++j) {
      const double* m = mu_data + j * n;
      Acc sum = 0.0;
      for (arma::uword i = 0; i < n; ++i) {
        const double term = weight[i] * (y[i] - m[i]);
        sum += term;
      }
      BB[lo + j] = static_cast<double>(sum / divisor);
    }
  }
  return BB;
}

}  // namespace

//' The subclassification moment condition of equation (F.2)
//'
//' The stratum treated fraction is empirical and fixed before the sweep starts,
//' so this path never evaluates a propensity score: one matrix product rather
//' than two, and no positivity guard to reach. Neither the propensity score
//' design matrix nor its draws appear in the signature, which makes that a
//' property of the interface rather than a branch someone could later add.
//'
//' @return One value per draw, sized from the outcome draws.
//'
//' The accumulator width is settled once, on entry, from what R reports about
//' its own build; see the note at the top of reduce.h.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector drb_moment_subclass(const Rcpp::NumericMatrix& Z_lm,
                                        const Rcpp::NumericVector& Y,
                                        const Rcpp::NumericVector& w,
                                        const Rcpp::NumericMatrix& betas_otc,
                                        int otc_link, double max_cells,
                                        bool wide) {
  return wide
    ? moment_subclass_impl<long double>(Z_lm, Y, w, betas_otc, otc_link,
                                        max_cells)
    : moment_subclass_impl<double>(Z_lm, Y, w, betas_otc, otc_link, max_cells);
}


namespace {

template <typename Acc>
Rcpp::List ps_rowmeans_impl(const Rcpp::NumericMatrix& Z_ps,
                            const Rcpp::NumericMatrix& betas_ps,
                            int ps_link, double max_cells) {
  const arma::uword n = static_cast<arma::uword>(Z_ps.nrow());
  const arma::uword S = static_cast<arma::uword>(betas_ps.nrow());

  if (static_cast<arma::uword>(betas_ps.ncol()) !=
      static_cast<arma::uword>(Z_ps.ncol())) {
    Rcpp::stop("the propensity score draws must have one column per design "
               "matrix column");
  }

  Rcpp::NumericVector ps_mean(n);
  drbayes::RangeTracker range;

  if (S > 0 && n > 0) {
    const arma::mat Zps = view(Z_ps);
    const arma::mat Bps = view(betas_ps);
    const arma::uword step = block_size(n, max_cells);
    std::vector<double> eta_buffer(static_cast<std::size_t>(n) *
                                   std::min(step, S));
    double* eta_data = eta_buffer.data();
    drbayes::RowMeanAccumulator<Acc> acc(n);

    for (arma::uword lo = 0; lo < S; lo += step) {
      const arma::uword hi = std::min(lo + step, S) - 1;
      const arma::uword width = hi - lo + 1;

      arma::mat eta = block(eta_data, n, width);
      eta = Zps * Bps.rows(lo, hi).t();

      range.observe(eta_data, n * width);
      drbayes::apply_inverse_link(eta_data, n * width, ps_link);
      acc.add_block(eta_data, n, width);
    }
    acc.write_means(ps_mean.begin());
  }

  Rcpp::NumericVector eta_range = range.empty()
    ? Rcpp::NumericVector(0)
    : Rcpp::NumericVector::create(range.low(), range.high());

  return Rcpp::List::create(Rcpp::Named("ps_mean") = ps_mean,
                            Rcpp::Named("eta_range") = eta_range);
}

}  // namespace

//' Mean propensity score per observation, over the draws
//'
//' The strata of Appendix F are cut on this, and it was the largest single
//' allocation in the package: the linear predictors and the scores were both
//' held at full draws-by-observations size at once, 244 MB at two thousand
//' observations and eight thousand draws, to produce a vector of length n.
//' Accumulating across blocks rather than averaging block means keeps the
//' answer the one an unblocked pass would have given, bit for bit.
//'
//' @return A list with `ps_mean`, one value per observation, and `eta_range`.
//'
//' The accumulator width is settled once, on entry, from what R reports about
//' its own build; see the note at the top of reduce.h.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export(rng = false)]]
Rcpp::List drb_ps_rowmeans(const Rcpp::NumericMatrix& Z_ps,
                           const Rcpp::NumericMatrix& betas_ps,
                           int ps_link, double max_cells, bool wide) {
  return wide ? ps_rowmeans_impl<long double>(Z_ps, betas_ps, ps_link, max_cells)
              : ps_rowmeans_impl<double>(Z_ps, betas_ps, ps_link, max_cells);
}


namespace {

template <typename Acc>
Rcpp::List moment_np_ipw_impl(const Rcpp::NumericMatrix& Z_ps,
                              const Rcpp::NumericVector& A,
                              const Rcpp::NumericVector& Y,
                              const Rcpp::NumericMatrix& mu,
                              const Rcpp::NumericMatrix& betas_ps,
                              int ps_link) {
  const arma::uword n = static_cast<arma::uword>(Y.size());
  const arma::uword S = static_cast<arma::uword>(mu.nrow());

  if (static_cast<arma::uword>(Z_ps.nrow()) != n ||
      static_cast<arma::uword>(A.size()) != n ||
      static_cast<arma::uword>(mu.ncol()) != n) {
    Rcpp::stop("Z.ps, A, Y and the fitted means must agree on the number of "
               "observations");
  }
  if (static_cast<arma::uword>(betas_ps.nrow()) < S) {
    Rcpp::stop("there must be at least one propensity score draw per row of "
               "fitted means");
  }

  Rcpp::NumericVector BB(S);
  drbayes::RangeTracker range;

  if (S > 0 && n > 0) {
    const arma::mat Zps = view(Z_ps);
    const arma::mat Bps = view(betas_ps);
    const arma::mat Mu = view(mu);
    const double* a = A.begin();
    const double* y = Y.begin();

    arma::vec eta(n);
    std::vector<double> term(n);

    for (arma::uword s = 0; s < S; ++s) {
      eta = Zps * Bps.row(s).t();
      double* eta_data = eta.memptr();

      range.observe(eta_data, n);
      drbayes::ipw_weight(eta_data, n, a, ps_link, 1);

      for (arma::uword i = 0; i < n; ++i) {
        term[i] = eta_data[i] * (y[i] - Mu(s, i));
      }
      BB[s] = drbayes::mean_two_pass<Acc>(term.data(), n);
    }
  }

  Rcpp::NumericVector eta_range = range.empty()
    ? Rcpp::NumericVector(0)
    : Rcpp::NumericVector::create(range.low(), range.high());

  return Rcpp::List::create(Rcpp::Named("BB") = BB,
                            Rcpp::Named("eta_range") = eta_range);
}

}  // namespace

//' The moment condition from supplied draws of the fitted outcome mean
//'
//' The fitted means arrive as draws by observations, the transpose of every
//' other matrix here, and are read one draw at a time. That is deliberate: R
//' forms each draw's linear predictor with a matrix-vector product, and
//' gathering the draws into a matrix product instead would be a different BLAS
//' kernel accumulating in a different order, moving the last bits for no gain,
//' since this is evaluated once per fit rather than once per sweep step.
//'
//' The reduction is the two-pass one, because R reaches this path through
//' mean() rather than colMeans().
//'
//' @return A list with `BB`, one value per draw, and `eta_range`.
//'
//' The accumulator width is settled once, on entry, from what R reports about
//' its own build; see the note at the top of reduce.h.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export(rng = false)]]
Rcpp::List drb_moment_np_ipw(const Rcpp::NumericMatrix& Z_ps,
                             const Rcpp::NumericVector& A,
                             const Rcpp::NumericVector& Y,
                             const Rcpp::NumericMatrix& mu,
                             const Rcpp::NumericMatrix& betas_ps,
                             int ps_link, bool wide) {
  return wide ? moment_np_ipw_impl<long double>(Z_ps, A, Y, mu, betas_ps, ps_link)
              : moment_np_ipw_impl<double>(Z_ps, A, Y, mu, betas_ps, ps_link);
}


namespace {

template <typename Acc>
Rcpp::NumericVector moment_np_subclass_impl(const Rcpp::NumericVector& Y,
                                            const Rcpp::NumericMatrix& mu,
                                            const Rcpp::NumericVector& w) {
  const arma::uword n = static_cast<arma::uword>(Y.size());
  const arma::uword S = static_cast<arma::uword>(mu.nrow());

  if (static_cast<arma::uword>(mu.ncol()) != n ||
      static_cast<arma::uword>(w.size()) != n) {
    Rcpp::stop("the fitted means and the stratum weights must agree with Y on "
               "the number of observations");
  }

  Rcpp::NumericVector BB(S);
  if (S == 0 || n == 0) return BB;

  const arma::mat Mu = view(mu);
  const double* y = Y.begin();
  const double* weight = w.begin();
  std::vector<double> term(n);

  for (arma::uword s = 0; s < S; ++s) {
    for (arma::uword i = 0; i < n; ++i) {
      term[i] = weight[i] * (y[i] - Mu(s, i));
    }
    BB[s] = drbayes::mean_two_pass<Acc>(term.data(), n);
  }
  return BB;
}

}  // namespace

//' The same, when the stratum weights of Appendix F stand in for the score
//'
//' No propensity score is formed at all, so there is no linear predictor and no
//' guard, and the weights are the same vector at every draw.
//'
//' @return One value per row of fitted means.
//'
//' The accumulator width is settled once, on entry, from what R reports about
//' its own build; see the note at the top of reduce.h.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector drb_moment_np_subclass(const Rcpp::NumericVector& Y,
                                           const Rcpp::NumericMatrix& mu,
                                           const Rcpp::NumericVector& w,
                                           bool wide) {
  return wide ? moment_np_subclass_impl<long double>(Y, mu, w)
              : moment_np_subclass_impl<double>(Y, mu, w);
}
