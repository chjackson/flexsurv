#include <cmath>
#include <exception>

#include <Rcpp.h>

#include "distribution.h"
#include "mapply_4.h"

namespace {

  namespace gengamma {
  
    inline double bad(const double mu,
		      const double sigma,
		      const double Q) {
      if (sigma < 0) {
	Rcpp::warning("Negative scale parameter \"sigma\"");
	return true;
      }
      return false;
    }
    
    class density {
    public:
      typedef double result_type;
      
      inline double operator()(const double x,
			       const double mu,
			       const double sigma,
			       const double Q) const {
	/* 
	   check the parameters 
	*/
	if (bad(mu, sigma, Q)) {
	  return NA_REAL;
	};
	
	if (x < 0) {
	  return R_NegInf;
	}
	
	/*
	  this function relies on interesting levels of cancellation if
	  Q is small or x is small
	*/
	
	if (Q!=0) {
	  const double y = std::log(x);
	  const double w = (y - mu) / sigma;
	  const double abs_q = std::abs(Q);
	  const double qi = 1/(Q * Q);
	  const double qw = Q * w;
	  return -std::log(sigma*x) +
	    std::log(abs_q) * (1 - 2 * qi) +
	    qi * (qw - std::exp(qw)) - R::lgammafn(qi);
	} else {
	  return R::dlnorm(x, mu, sigma, 1);
	}
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
			       const double big_q) const {
	/* check the arguments */
	if (bad(mu, sigma, big_q)) {
	  return NA_REAL;
	}
	
	if ( q < 0 ) {
	  return below_distribution(lower_tail, give_log);
	}
	
	if (big_q!=0) {
	  const double y = std::log(q);
	  const double w = (y - mu) / sigma;
	  const double qq = 1/(big_q * big_q);
	  const double expnu = std::exp(big_q * w) * qq;
	  if (big_q > 0) {
	    return R::pgamma(expnu, qq, 1, lower_tail, give_log);
	  } else {
	    return R::pgamma(expnu, qq, 1, !lower_tail, give_log);
	  }
	} else {
	  return R::plnorm(q, mu, sigma, lower_tail, give_log); 
	}
      }
      
    private:
      bool lower_tail;
      bool give_log;
    };
  }
}



// [[Rcpp::export(rng=false)]]
Rcpp::NumericVector
dgengamma_work(const Rcpp::NumericVector& x,
	       const Rcpp::NumericVector& mu,
	       const Rcpp::NumericVector& sigma,
	       const Rcpp::NumericVector& Q,
	       const bool log) {
  R_xlen_t size
    = std::max(std::max(sigma.size(),
			Q.size()),
	       std::max(x.size(),
			mu.size()));
  return perhaps_exp(mapply(Rcpp::rep_len(x, size),
			    Rcpp::rep_len(mu, size),
			    Rcpp::rep_len(sigma, size),
			    Rcpp::rep_len(Q, size),
			    gengamma::density()),
		     log);
}

// [[Rcpp::export(rng=false)]]
Rcpp::NumericVector
pgengamma_work(const Rcpp::NumericVector& q,
	       const Rcpp::NumericVector& mu,
	       const Rcpp::NumericVector& sigma,
	       const Rcpp::NumericVector& Q,
	       const bool lower_tail,
	       const bool give_log) {
  R_xlen_t size
    = std::max(std::max(sigma.size(),
			Q.size()),
	       std::max(q.size(),
			mu.size()));

  return mapply(Rcpp::rep_len(q, size),
		Rcpp::rep_len(mu, size),
		Rcpp::rep_len(sigma, size),
		Rcpp::rep_len(Q, size),
		gengamma::cdf(lower_tail, give_log));
}

// [[Rcpp::export(name="check.gengamma", rng=false)]]
Rcpp::LogicalVector check_gengamma(const Rcpp::NumericVector& mu,
				   const Rcpp::NumericVector& sigma,
				   const Rcpp::NumericVector& Q) {
  return !(Rcpp::mapply(mu, sigma, Q, gengamma::bad));
}
