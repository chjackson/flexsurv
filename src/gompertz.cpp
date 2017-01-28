#include <cmath>
#include <exception>

#include <Rcpp.h>

namespace {
  
  inline double exprel(const double x) {
    if (x!=0) {
      return expm1(x) / x;
    } else {
      return 1;
    }
  }
  
  class dgompertz {
  public:
    typedef double result_type;

    inline double operator()(const double x,
			     const double shape,
			     const double rate) const {

      const double scale_x = shape*x;
      const double shift   = x * exprel(scale_x);
      return std::log(rate) + scale_x - rate * shift;
    }
  };

}

// [[Rcpp::export]]
Rcpp::NumericVector
dgompertz_work(const Rcpp::NumericVector& x,
	       const Rcpp::NumericVector& shape,
	       const Rcpp::NumericVector& rate,
	       const bool log) {
  if (log) {
    return Rcpp::mapply(x, shape, rate, dgompertz());
  } else {
    return Rcpp::exp(mapply(x, shape, rate, dgompertz()));
  }
}

