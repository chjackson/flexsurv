#include <cmath>
#include <exception>

#include <Rcpp.h>

namespace {
  inline double dgengamma(const double x,
			  const double mu,
			  const double sigma,
			  const double Q) {
    if (Q!=0) {
      const double y = std::log(x);
      const double w = (y - mu) / sigma;
      const double abs_q = std::abs(Q);
      const double qi = 1/(Q * Q);
      const double qw = Q * w;
      return -std::log(sigma*x) + std::log(abs_q) +
	qi * std::log(qi) + qi * (qw - std::exp(qw)) - R::lgammafn(qi);
    } else {
      return R::dlnorm(x, mu, sigma, 1);
    }
  }

  inline double pgengamma(const double q,
			  const double mu,
			  const double sigma,
			  const double big_q,
			  bool lower_tail,
			  bool give_log) {
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
}

// [[Rcpp::export]]
Rcpp::NumericVector
dgengamma_work(const Rcpp::NumericVector& x,
	       const Rcpp::NumericVector& mu,
	       const Rcpp::NumericVector& sigma,
	       const Rcpp::NumericVector& Q,
	       const bool log) {
  if ( (x.size() != mu.size()) ||
       (x.size() != sigma.size()) ||
       (x.size() != Q.size()) ) {
    throw std::runtime_error("Vector size mismatch in dgengamma_work.");
  }
  Rcpp::NumericVector result(x.size());
  for (R_xlen_t ind=0; ind!=x.size(); ++ind) {
    const double res = dgengamma(x[ind], mu[ind], sigma[ind], Q[ind]);
    if (log) {
      result[ind] = res;
    } else {
      result[ind] = std::exp(res);
    }
  }
  return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector
pgengamma_work(const Rcpp::NumericVector& q,
	       const Rcpp::NumericVector& mu,
	       const Rcpp::NumericVector& sigma,
	       const Rcpp::NumericVector& Q,
	       const bool lower_tail,
	       const bool give_log) {
  if ( (q.size() != mu.size()) ||
       (q.size() != sigma.size()) ||
       (q.size() != Q.size()) ) {
    throw std::runtime_error("Vector size mismatch in pgengamma_work.");
  }
  Rcpp::NumericVector result(q.size());
  for (R_xlen_t ind=0; ind!=q.size(); ++ind) {
    result[ind] = pgengamma(q[ind], mu[ind], sigma[ind], Q[ind],
			    lower_tail, give_log);
  }
  return result;
}
