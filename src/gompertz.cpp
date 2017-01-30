#include <cmath>
#include <exception>

#include <Rcpp.h>

#include "distribution.h"

namespace {

  namespace gompertz {
  
    inline double exprel(const double x) {
      if (x!=0) {
	return expm1(x) / x;
      } else {
	return 1;
      }
    }
    
    inline double bad(const double shape, const double rate) {
      if (rate < 0) {
	Rcpp::warning("Negative rate parameter");
	return true;
      } else {
	return false;
      }
    }
    
    class density {
    public:
      typedef double result_type;
      
      inline double operator()(const double x,
			       const double shape,
			       const double rate) const {
	if (bad(shape, rate)) {
	  return NA_REAL;
	}

	if (x < 0) {
	  return R_NegInf;
	}
	
	const double scale_x = shape * x;
	const double shift   = x * exprel(scale_x);
	return std::log(rate) + scale_x - rate * shift;
      }
    };

  }

}

// [[Rcpp::export(rng=false)]]
Rcpp::NumericVector
dgompertz_work(const Rcpp::NumericVector& x,
	       const Rcpp::NumericVector& shape,
	       const Rcpp::NumericVector& rate,
	       const bool log) {
  R_xlen_t size = std::max(x.size(),
			   std::max(shape.size(),
				    rate.size()));
  
  return perhaps_exp(Rcpp::mapply(Rcpp::rep_len(x, size),
				  Rcpp::rep_len(shape, size),
				  Rcpp::rep_len(rate, size),
				  gompertz::density()),
		     log);
}

