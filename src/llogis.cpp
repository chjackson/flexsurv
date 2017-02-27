#include <cmath>

#include <Rcpp.h>

#include "distribution.h"
#include "rep_len.h"

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

  if (x.size() == 0) { return x; }
  
  const R_xlen_t size
    = std::max(x.size(),
	       std::max(shape.size(),
			scale.size()));
  
  return perhaps_exp(Rcpp::mapply(flexsurv::rep_len(x, size),
				  flexsurv::rep_len(shape, size),
				  flexsurv::rep_len(scale, size),
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
  if (q.size() == 0) { return q; }
  
  const R_xlen_t size
    = std::max(q.size(),
	       std::max(shape.size(),
			scale.size()));

  return mapply(flexsurv::rep_len(q, size),
		flexsurv::rep_len(shape, size),
		flexsurv::rep_len(scale, size),
		llogis::cdf(lower_tail, give_log));
}

// [[Rcpp::export(name="check.llogis", rng=false)]]
Rcpp::LogicalVector check_llogis(const Rcpp::NumericVector& shape,
				 const Rcpp::NumericVector& scale) {
  if ( (0==shape.size()) && (0==scale.size()) ) {
    Rcpp::LogicalVector null_result(0);
    return null_result;
  }
  const R_xlen_t size = shape.size();
  return !Rcpp::mapply(shape,
		       flexsurv::rep_len(scale, size),
		       llogis::bad);
}
