#ifndef DRBAYES_REDUCE_H
#define DRBAYES_REDUCE_H

#include <vector>
#include <Rmath.h>
#include <cstddef>

namespace drbayes {

// R accumulates sums in LDOUBLE, which is long double on a build that has it
// and double on one that does not, and the kernels here have to follow that
// choice or they stop agreeing with the R they replace.
//
// The choice cannot be read at compile time. It is HAVE_LONG_DOUBLE in R's own
// config.h, which is private: the installed Rconfig.h carries fifteen macros
// and that is not one of them, so an #ifdef on it is never true and would pin
// every platform to double. R is asked instead, through
// capabilities("long.double"), and the answer arrives as an argument. That also
// covers the one case a compile-time test could not have got right in any form:
// CRAN's no-long-double flavour, where the platform has a wide long double and
// R has been built deliberately not to use it.
//
// Each reduction is therefore a template over the accumulator, and each kernel
// picks the instantiation once, on entry, rather than testing per element.


// The mean colMeans computes: one pass, one division at the end.
//
// This and mean_two_pass below are not interchangeable. R uses colMeans for the
// parametric moment conditions and mean() for the nonparametric one, and the
// two disagree in the last bits on ordinary data, so a single reduction here
// would silently move the values on whichever path did not get its own.
template <typename Acc>
inline double mean_one_pass(const double* x, std::size_t n) {
  Acc sum = 0.0;
  for (std::size_t i = 0; i < n; ++i) sum += x[i];
  return static_cast<double>(sum / static_cast<Acc>(n));
}


// The mean mean() computes: a second pass adds back the accumulated rounding
// of the first. R skips the correction when the first pass is not finite,
// because subtracting an infinity from each element would give NaN where the
// answer should be the infinity itself.
template <typename Acc>
inline double mean_two_pass(const double* x, std::size_t n) {
  Acc s = 0.0;
  for (std::size_t i = 0; i < n; ++i) s += x[i];
  s /= static_cast<Acc>(n);

  if (R_FINITE(static_cast<double>(s))) {
    Acc t = 0.0;
    for (std::size_t i = 0; i < n; ++i) t += (x[i] - s);
    s += t / static_cast<Acc>(n);
  }
  return static_cast<double>(s);
}


// Accumulates the row sums of a column-major matrix one block of columns at a
// time, which is how rowMeans walks its input. Carrying the accumulator across
// blocks rather than averaging block means is what makes the blocked result
// equal to the unblocked one rather than merely close to it.
template <typename Acc>
class RowMeanAccumulator {
 public:
  explicit RowMeanAccumulator(std::size_t rows)
      : sums_(rows, Acc(0.0)), columns_(0) {}

  void add_block(const double* block, std::size_t rows, std::size_t cols) {
    for (std::size_t j = 0; j < cols; ++j) {
      const double* column = block + j * rows;
      for (std::size_t i = 0; i < rows; ++i) sums_[i] += column[i];
    }
    columns_ += cols;
  }

  void write_means(double* out) const {
    const Acc divisor = static_cast<Acc>(columns_);
    for (std::size_t i = 0; i < sums_.size(); ++i) {
      out[i] = static_cast<double>(sums_[i] / divisor);
    }
  }

 private:
  std::vector<Acc> sums_;
  std::size_t columns_;
};


// The smallest and largest value seen, with missingness carried rather than
// skipped.
//
// A comparison against NaN is false whichever way it is written, so the obvious
// running minimum steps over a NaN without recording it. R's min() and max()
// return NA when they meet one, and the positivity guard depends on that: a
// propensity score that cannot be evaluated must refuse, not pass.
class RangeTracker {
 public:
  void observe(const double* x, std::size_t n) {
    for (std::size_t i = 0; i < n; ++i) {
      const double v = x[i];
      if (ISNAN(v)) {
        // R keeps the two apart and lets NA win, so min(c(NA, NaN)) is NA
        // while min(c(0.5, NaN)) is NaN. Collapsing them would report a
        // missing propensity score as an unevaluable one, or the reverse.
        if (ISNA(v)) na_ = true; else nan_ = true;
        continue;
      }
      if (!seen_ || v < low_)  low_  = v;
      if (!seen_ || v > high_) high_ = v;
      seen_ = true;
    }
  }

  bool empty() const { return !seen_ && !na_ && !nan_; }
  double low()  const { return missing() ? missing_value() : low_; }
  double high() const { return missing() ? missing_value() : high_; }

 private:
  bool missing() const { return na_ || nan_; }
  double missing_value() const { return na_ ? NA_REAL : R_NaN; }

  bool seen_ = false;
  bool na_ = false;
  bool nan_ = false;
  double low_ = 0.0;
  double high_ = 0.0;
};

}  // namespace drbayes

#endif  // DRBAYES_REDUCE_H
