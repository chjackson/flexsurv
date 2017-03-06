#include <cmath>

#include <Rcpp.h>

#include "distribution.h"
#include "rep_len.h"
#include "gengamma.h"
#include "mapply.h"

namespace {

  namespace genf {

    inline bool bad(const double mu,
		    const double sigma,
		    const double Q,
		    const double P) {
      bool res = false;
      if (sigma < 0) {
	Rcpp::warning("Negative scale parameter ""sigma""");
	res = true;
      }
      if (P < 0) {
	Rcpp::warning("Negative shape parameter ""P""");
	res = true;
      }
      return res;
    }
    
    class density {
    public:
      typedef double result_type;
      
      inline double operator()(const double x,
			       const double mu,
			       const double sigma,
			       const double Q,
			       const double P) const {
	/* 
	   check the parameters 
	*/
	if (bad(mu, sigma, Q, P)) {
	  return NA_REAL;
	};
	
	if (x < 0) {
	  return R_NegInf;
	}

	if (P==0) {
	  return gengamma::density()(x, mu, sigma, Q);
	}

	const double tmp = Q * Q + 2 * P;
	const double delta = std::sqrt(tmp);
        const double s1 = 2 / (tmp + Q*delta);
        const double s2 = 2 / (tmp - Q*delta);
	const double expw = std::pow(x, delta/sigma) *
	  std::exp(-mu*delta/sigma);
	return std::log(delta)
	  + s1/sigma*delta*(std::log(x) - mu)
	  + s1*(std::log(s1) - std::log(s2))
	  - std::log(sigma*x)
	  - (s1+s2)*std::log(1 + s1*expw/s2) - R::lbeta(s1, s2);
      }
    };

    class cdf {
    public:
      typedef double result_type;
      
      cdf(bool lower_tail_, bool give_log_) :
	lower_tail(lower_tail_),
	give_log(give_log_) {}
      
      inline double operator()(const double q,
			       const double mu,
			       const double sigma,
			       const double Q,
			       const double P) const {
	if (bad(mu, sigma, Q, P)) {
	  return NA_REAL;
	}

	if ( q < 0 ) {
	  return below_distribution(lower_tail, give_log);
	}

	if (P==0) {
	  return gengamma::cdf(lower_tail, give_log)(q, mu, sigma, Q);
	}

	const double tmp = Q * Q + 2*P;
        const double delta = std::sqrt(tmp);
	const double s1 = 2 / (tmp + Q*delta);
        const double s2 = 2 / (tmp - Q*delta);
        const double expw = std::pow(q, delta/sigma) *
	  std::exp(-mu*delta/sigma);
	return R::pbeta(s2/(s2 + s1*expw), s2, s1, !lower_tail, give_log);
      }
    private:
      bool lower_tail;
      bool give_log;
    };
  }
}

// [[Rcpp::export(rng=false)]]
Rcpp::NumericVector
dgenf_work(const Rcpp::NumericVector& x,
	   const Rcpp::NumericVector& mu,
	   const Rcpp::NumericVector& sigma,
	   const Rcpp::NumericVector& Q,
	   const Rcpp::NumericVector& P,
	   const bool log) {
  if (x.size() == 0) { return x; }

  const R_xlen_t size
    = std::max(P.size(),
	       std::max(std::max(sigma.size(),
				 Q.size()),
			std::max(x.size(),
				 mu.size())));
  return perhaps_exp(mapply(flexsurv::rep_len(x, size),
			    flexsurv::rep_len(mu, size),
			    flexsurv::rep_len(sigma, size),
			    flexsurv::rep_len(Q, size),
			    flexsurv::rep_len(P, size),
			    genf::density()),
		     log);
}

// [[Rcpp::export(rng=false)]]
Rcpp::NumericVector
pgenf_work(const Rcpp::NumericVector& q,
	   const Rcpp::NumericVector& mu,
	   const Rcpp::NumericVector& sigma,
	   const Rcpp::NumericVector& Q,
	   const Rcpp::NumericVector& P,
	   const bool lower_tail,
	   const bool give_log) {
  if (q.size() == 0) { return q; }
  
  const R_xlen_t size
    = std::max(P.size(),
	       std::max(std::max(sigma.size(),
				 Q.size()),
			std::max(q.size(),
				 mu.size())));
  
  return mapply(flexsurv::rep_len(q, size),
		flexsurv::rep_len(mu, size),
		flexsurv::rep_len(sigma, size),
		flexsurv::rep_len(Q, size),
		flexsurv::rep_len(P, size),
		genf::cdf(lower_tail, give_log));
}

namespace Rcpp {
  namespace traits {
    template <typename RESULT_TYPE, typename U1, typename U2, typename U3, typename U4>
    struct result_of< RESULT_TYPE (*)(U1, U2, U3, U4) >{
      typedef RESULT_TYPE type ;
    } ;
  }
}

// [[Rcpp::export(name="check.genf", rng=false)]]
Rcpp::LogicalVector check_genf(const Rcpp::NumericVector& mu,
			       const Rcpp::NumericVector& sigma,
			       const Rcpp::NumericVector& Q,
			       const Rcpp::NumericVector& P) {
  if ( (0==mu.size()) &&
       (0==sigma.size()) &&
       (0==Q.size()) &&
       (0==P.size()) ) {
    Rcpp::LogicalVector null_result(0);
    return null_result;
  }
  
  const R_xlen_t size = mu.size();
  return !mapply(mu,
		 flexsurv::rep_len(sigma, size),
		 flexsurv::rep_len(Q, size),
		 flexsurv::rep_len(P, size),
		 genf::bad);
}
