#ifndef GENGAMMA_H
#define GENGAMMA_H

#include <cmath>

#include <Rcpp.h>
 
#include "distribution.h"

namespace {

  namespace gengamma {
  
    inline bool bad(const double mu,
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
	
	if (Q!=0) {
	  const double y = std::log(x);
	  const double w = (y - mu) / sigma;
	  const double abs_q = std::abs(Q);
	  const double qw = Q * w;
          if (abs_q < 0.05 && std::isfinite(w)) {
            // The next log-density correction, -Q^10/1260, has magnitude < 8e-17.
            const double q2 = Q * Q;
            double deviance;
            if (std::abs(qw) < 0.5) {
              // Taylor terms through order 16 avoid cancellation for |qw| < 0.5.
              double term = 0.5, sum = term;
              for (int k = 3; k <= 16; ++k) {
                term *= qw / k;
                sum += term;
              }
              deviance = w * w * sum;
            } else {
              deviance = (std::expm1(qw) - qw) / q2;
            }
            // Stirling expansion of lgamma(Q^-2).
            return -std::log(sigma) - std::log(x) - M_LN_SQRT_2PI
              - deviance - q2 / 12 + q2 * q2 * q2 / 360;
          }
          const double qi = 1/(Q * Q);
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

#endif
