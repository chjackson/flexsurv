#include <cmath>
#include <Rcpp.h>

namespace {
  inline double restrict_hypot(double x) {
    if (std::isnan(x)) {
      return x;
    }
    x = std::abs(x);
    if (x > 1) {
      const double inv = 1/x;
      return std::sqrt(1 + inv * inv) * x;
    } else {
      return std::sqrt(1 + x * x);
    }
  }

  inline double dexph_work(const double x) {
    if (std::isnan(x)) {
      return x;
    }
    if (x >= 0) {
      return 1 + x/restrict_hypot(x);
    } else {
      const double h = restrict_hypot(x);
      return -1 / h / (x - h);
    }
  }
}

// [[Rcpp::export(rng=false)]]
Rcpp::NumericVector exph(const Rcpp::NumericVector& y) {
  // exp(asinh(x))
  return Rcpp::ifelse(y >= 0,
		      y + Rcpp::sapply(y, restrict_hypot),
		      -1/(y - Rcpp::sapply(y, restrict_hypot)));
}

// [[Rcpp::export(rng=false)]]
Rcpp::NumericVector dexph(const Rcpp::NumericVector& y) {
  return Rcpp::sapply(y, dexph_work);		      
}
