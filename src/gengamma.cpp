#include <exception>

#include <Rcpp.h>

#include "distribution.h"
#include "gengamma.h"
#include "mapply.h"

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
  const R_xlen_t size = mu.size();
  return !Rcpp::mapply(mu,
		       Rcpp::rep_len(sigma, size),
		       Rcpp::rep_len(Q, size),
		       gengamma::bad);
}
