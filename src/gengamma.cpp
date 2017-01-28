#include <cmath>
#include <exception>

#include <Rcpp.h>

#include "mapply_4.h"

namespace {
  class dgengamma {
  public:
    typedef double result_type;

    inline double operator()(const double x,
			     const double mu,
			     const double sigma,
			     const double Q) const {

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

  class pgengamma {
  public:
    typedef double result_type;
    
    pgengamma(bool lower_tail_, bool give_log_) :
      lower_tail(lower_tail_),
      give_log(give_log_) {}
    
    inline double operator()(const double q,
			     const double mu,
			     const double sigma,
			     const double big_q) const {
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

// [[Rcpp::export]]
Rcpp::NumericVector
dgengamma_work(const Rcpp::NumericVector& x,
	       const Rcpp::NumericVector& mu,
	       const Rcpp::NumericVector& sigma,
	       const Rcpp::NumericVector& Q,
	       const bool log) {
  if (log) {
    return mapply(x, mu, sigma, Q, dgengamma());
  } else {
    return Rcpp::exp(mapply(x, mu, sigma, Q, dgengamma()));
  }
}

// [[Rcpp::export]]
Rcpp::NumericVector
pgengamma_work(const Rcpp::NumericVector& q,
	       const Rcpp::NumericVector& mu,
	       const Rcpp::NumericVector& sigma,
	       const Rcpp::NumericVector& Q,
	       const bool lower_tail,
	       const bool give_log) {
  return mapply(q, mu, sigma, Q, pgengamma(lower_tail, give_log));
}
