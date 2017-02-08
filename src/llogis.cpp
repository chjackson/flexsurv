#include <cmath>

#include <Rcpp.h>

#include "distribution.h"

namespace {

  namespace llogis {
  
    inline bool bad(const double shape,
		    const double scale) {
      const bool bad = (shape<=0) || (scale <= 0);

      if (shape < 0) {
	Rcpp::warning("Non-positive shape parameter");
      }
      if (scale < 0) {
	Rcpp::warning("Non-positive scale parameter");
      }
      
      return bad;
      
    }
    
    class density {
    public:
      typedef double result_type;
      
      inline double operator()(const double x,
			       const double shape,
			       const double scale) const {
	/* 
	   check the parameters 
	*/
	if (bad(shape, scale)) {
	  return NA_REAL;
	};
	
	if (x < 0) {
	  return R_NegInf;
	}

	const double log_shape = std::log(shape);
	const double log_scale = std::log(scale);

	return log_shape - log_scale +
	  (shape - 1) * (std::log(x) - log_scale)
	  - 2 * std::log(1 + std::pow(x/scale, shape));
      }
    };

    class cdf {
    public:
      typedef double result_type;
      
      cdf(bool lower_tail_, bool give_log_) :
	lower_tail(lower_tail_),
	give_log(give_log_) {}
      
      inline double operator()(const double q,
			       const double shape,
			       const double scale) const {
	/* check the arguments */
	if (bad(shape, scale)) {
	  return NA_REAL;
	}
	
	if ( q < 0 ) {
	  return below_distribution(lower_tail, give_log);
	}
	return R::plogis(std::log(q), std::log(scale), 1/shape,
			 lower_tail, give_log);
      }
      
    private:
      bool lower_tail;
      bool give_log;
    };
  }
}



// [[Rcpp::export(rng=false)]]
Rcpp::NumericVector
dllogis_work(const Rcpp::NumericVector& x,
	     const Rcpp::NumericVector& shape,
	     const Rcpp::NumericVector& scale,
	     const bool log) {
  R_xlen_t size
    = std::max(x.size(),
	       std::max(shape.size(),
			scale.size()));
	       
  return perhaps_exp(Rcpp::mapply(Rcpp::rep_len(x, size),
				  Rcpp::rep_len(shape, size),
				  Rcpp::rep_len(scale, size),
				  llogis::density()),
		     log);
}

// [[Rcpp::export(rng=false)]]
Rcpp::NumericVector
pllogis_work(const Rcpp::NumericVector& q,
	     const Rcpp::NumericVector& shape,
	     const Rcpp::NumericVector& scale,
	     const bool lower_tail,
	     const bool give_log) {
  R_xlen_t size
    = std::max(q.size(),
	       std::max(shape.size(),
			scale.size()));

  return mapply(Rcpp::rep_len(q, size),
		Rcpp::rep_len(shape, size),
		Rcpp::rep_len(scale, size),
		llogis::cdf(lower_tail, give_log));
}

// [[Rcpp::export(name="check.llogis", rng=false)]]
Rcpp::LogicalVector check_llogis(const Rcpp::NumericVector& shape,
				 const Rcpp::NumericVector& scale) {
  const R_xlen_t size = shape.size();
  return !Rcpp::mapply(shape,
		       Rcpp::rep_len(scale, size),
		       llogis::bad);
}
